/*
 * inverse.cc
 * 
 * Jun Korenaga, MIT/WHOI
 * January 1999
 */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <geom.h>
#include <util.h>
#include "inverse.h"
#include "sparse_rect.h"
#include "lsqr.h"

//
// TomographicInversion2d
//
TomographicInversion2d::TomographicInversion2d(SlownessMesh2d& m, const char* datafn,
					       int xorder, int zorder, double crit_len,
					       int nintp, double cg_tol, double br_tol)

    : smesh(m), graph(smesh,xorder,zorder),
      betasp(1,0,nintp), bend(smesh,betasp,cg_tol,br_tol),
      nrefl(0), nnoded(0), do_full_refl(false),
      nnodev(smesh.numNodes()), nx(smesh.Nx()), nz(smesh.Nz()),
      itermax_LSQR(nnodev*100), LSQR_ATOL(1e-3), 
      smooth_velocity(false), logscale_vel(false), do_filter(false),
      wsv_min(0.0), wsv_max(0.0), dwsv(1.0), weight_s_v(0.0),
      smooth_depth(false), logscale_dep(false), wsd_min(0.0),
      wsd_max(0.0), dwsd(1.0), weight_s_d(0.0), corr_dep_p(0),
      damp_velocity(false), damp_depth(false),
      target_dv(0.0), target_dd(0.0),
      damping_is_fixed(false), weight_d_v(0.0), weight_d_d(0.0),
      do_squeezing(false),
      robust(false), crit_chi(1e3), refl_weight(1.0), target_chisq(0.0),
      printLog(false), verbose_level(-1),
      printFinal(false), printTransient(false), out_mask(false), out_level(0),
      gravity(false), out_grav_dws(false), ngravdata(0)
{
    if (crit_len>0) graph.refineIfLong(crit_len);

    // there can be two interfaces at most (seafloor and moho).
    start_i.reserve(2); end_i.reserve(2);
    interp.reserve(2);
    bathyp = new Interface2d(smesh);
    
    read_file(datafn);

    const int max_np = int(2*sqrt(float(nnodev)));
    path.reserve(max_np);
    pp.reserve(max_np+4);
    
    A.resize(ndata); 
    Rv_h.resize(nnodev); Rv_v.resize(nnodev);
    Tv.resize(nnodev);
    data_vec.resize(ndata);
    tmp_data.resize(ndata);
    total_data_vec.reserve(2*ndata+3*nnodev+2*max_np);
    r_dt_vec.resize(ndata);
    int idata=1;
    for (int i=1; i<=obs_dt.size(); i++){
	for (int j=1; j<=obs_dt(i).size(); j++){
	    r_dt_vec(idata) = 1.0/obs_dt(i)(j);
	    idata++;
	}
    }
    modelv.resize(nnodev);
    mvscale.resize(nnodev);
    dmodel_total.resize(nnodev);
    nnode_total = nnodev;
    path_length.resize(ndata);

    Q.resize(nintp);

    nodev_hit.resize(nnodev); 
    tmp_node.resize(nnodev);
    tmp_nodev.resize(nnodev);
}

void TomographicInversion2d::read_file(const char* datafn)
{
    ifstream in(datafn);
    if (!in){
	cerr << "TomographicInversion2d::cannot open " << datafn << "\n";
	exit(1);
    }

    int iline=0;

    int nsrc;
    in >> nsrc; iline++;
    if (nsrc<=0) error("TomographicInversion2d::invalid nsrc");
    src.resize(nsrc); rcv.resize(nsrc);
    raytype.resize(nsrc); obs_ttime.resize(nsrc); obs_dt.resize(nsrc);
    res_ttime.resize(nsrc);
    
    int isrc=0;
    while(in){
	char flag;
	double x, y;
	int nrcv;

	in >> flag >> x >> y >> nrcv; iline++;
	if (flag!='s'){
	    cerr << "TomographicInversion2d::bad input (s) at l."
		 << iline << '\n';
	    exit(1);
	}
	isrc++;
	src(isrc).set(x,y);

	rcv(isrc).resize(nrcv); raytype(isrc).resize(nrcv);
	obs_ttime(isrc).resize(nrcv); obs_dt(isrc).resize(nrcv);
	res_ttime(isrc).resize(nrcv);
	for (int ircv=1; ircv<=nrcv; ircv++){
	    int n;
	    double t_val, dt_val;
	    in >> flag >> x >> y >> n >> t_val >> dt_val; iline++;
	    if (flag!='r'){
		cerr << "TomographicInversion2d::bad input (r) at l."
		     << iline << '\n';
		exit(1);
	    }
	    if (dt_val < 1e-6){
		cerr << "TomographicInversion2d::bad input (zero obs_dt) at l."
		     << iline << '\n';
		exit(1);
	    }
	    rcv(isrc)(ircv).set(x,y);
	    raytype(isrc)(ircv) = n;
	    obs_ttime(isrc)(ircv) = t_val;
	    obs_dt(isrc)(ircv) = dt_val;
	}
	if (isrc==nsrc) break;
    }
    if (isrc != nsrc) error("TomographicInversion2d::mismatch in nsrc");

    ndata = 0;
    for (int isrc=1; isrc<=nsrc; isrc++){
	ndata += rcv(isrc).size();
    }
}

void TomographicInversion2d::doRobust(double a)
{
    if (a>0){
	robust = true;
	crit_chi = a;
    }
}

void TomographicInversion2d::removeOutliers()
{
    Array1d<double> Adm(ndata_valid);
    SparseRectangular sparseA(A,tmp_data,tmp_node,nnode_total);
    sparseA.Ax(dmodel_total,Adm);
    Array1d<int> tmp_data2(tmp_data.size());

    tmp_data2 = 0;
    int ndata_valid2=0;
    for (int i=1; i<=tmp_data.size(); i++){
	int j=tmp_data(i);
	if (j>0){
	    double res=Adm(j)-data_vec(i);
	    if (abs(res)<=crit_chi){
		tmp_data2(i) = ++ndata_valid2;
	    }
	}
    }

    // redefine data node vector and related stats
    ndata_valid = ndata_valid2;
    tmp_data = tmp_data2;
    rms_tres[0]=rms_tres[1]=0;
    init_chi[0]=init_chi[1]=0;
    ndata_in[0]=ndata_in[1]=0;
    int idata=1;
    for (int isrc=1; isrc<=src.size(); isrc++){
	for (int ircv=1; ircv<=rcv(isrc).size(); ircv++){
	    if (tmp_data(idata)>0){
		double res=res_ttime(isrc)(ircv);
		int icode = raytype(isrc)(ircv);
		double res2=res*r_dt_vec(idata);
		double res22=res2*res2;
		rms_tres[icode] += res*res;
		init_chi[icode] += res22;
		++ndata_in[icode];
	    }
	    idata++;
	}
    }
    rms_tres_total = sqrt((rms_tres[0]+rms_tres[1])/ndata_valid);
    init_chi_total = (init_chi[0]+init_chi[1])/ndata_valid;
    for (int i=0; i<=1; i++){
	rms_tres[i] = ndata_in[i]>0 ? sqrt(rms_tres[i]/ndata_in[i]) : 0.0;
	init_chi[i] = ndata_in[i]>0 ? init_chi[i]/ndata_in[i] : 0.0;
    }
}

void TomographicInversion2d::setLSQR_TOL(double a)
{
    if (a<=0){
	cerr << "TomographicInversion2d::setLSQR_TOL - invalid TOL (ignored)\n";
	return;
    }
    LSQR_ATOL = a;
}

void TomographicInversion2d::setLogfile(const char* fn)
{
    printLog = true;
    log_os_p = new ofstream(fn);
    if (!(*log_os_p)){
	cerr << "TomographicInversion2d::setLogfile - can't open "
	     << fn << '\n';
	exit(1);
    }
}

void TomographicInversion2d::setVerbose(int i){ verbose_level=i; }

void TomographicInversion2d::outStepwise(const char* fn, int i)
{
    printTransient = true;
    out_root = fn;
    out_level = i;
}

void TomographicInversion2d::outFinal(const char* fn, int i)
{
    printFinal = true;
    out_root = fn;
    out_level = i;
}

void TomographicInversion2d::targetChisq(double c)
{
    target_chisq = c;
}

void TomographicInversion2d::outMask(const char* fn)
{
    out_mask = true;
    vmesh_os_p = new ofstream(fn);
    if (!(*vmesh_os_p)){
	cerr << "TomographicInversion2d::outMask - can't open "
	     << fn << '\n';
	exit(1);
    }
}

void TomographicInversion2d::addonGravity(double _z, const Array1d<double>& _x,
					  const Array1d<double>& _g,
					  AddonGravityInversion2d* p, double w,
					  const char *dwsfn)
{
    out_grav_dws = true;
    grav_dws_osp = new ofstream(dwsfn);
    if (!(*grav_dws_osp)){
	cerr << "TomographicInversion2d::addonGravity - can't open "
	     << dwsfn << '\n';
	exit(1);
    }
    addonGravity(_z,_x,_g,p,w);
}

void TomographicInversion2d::addonGravity(double _z, const Array1d<double>& _x,
					  const Array1d<double>& _g,
					  AddonGravityInversion2d* p, double w)
{
    gravity = true;
    grav_z0 = _z;
    ngravdata = _x.size();
    grav_x.resize(ngravdata); grav_x = _x;
    obs_grav.resize(ngravdata); obs_grav = _g;
    res_grav.resize(ngravdata);
    B.resize(ngravdata);
    tmp_gravdata.resize(ngravdata);
    for (int i=1; i<=ngravdata; i++) tmp_gravdata(i) = i;
    ginv = p;
    weight_grav = w;
    if (w<0){
	error("TomographicInversion2d::addonGravity - negative weight detected.");
    }
}

void TomographicInversion2d::solve(int niter)
{
    typedef map<int,double>::iterator mapIterator;
    typedef map<int,double>::const_iterator mapBrowser;
    
    const int bend_nfac=1;
    if (printLog){
	*log_os_p << "# strategy " << jumping << " " << robust << " " << crit_chi << '\n';
	*log_os_p << "# ray_trace " << graph.xOrder() << " " << graph.zOrder() << " "
		  << graph.critLength() << " "
		  << betasp.numIntp() << " " << bend.tolerance() << '\n';
	*log_os_p << "# smooth_vel " << smooth_velocity << " "
		  << wsv_min << " " << wsv_max << " " << dwsv << " "
		  << do_filter << '\n';
	*log_os_p << "# smooth_dep " << smooth_depth << " "
		  << wsd_min << " " << wsd_max << " " << dwsd << '\n';
	if (damping_is_fixed){
	    *log_os_p << "# fixed_damping " << weight_d_v << " " << weight_d_d << '\n';
	}else{
	    *log_os_p << "# damp_vel " << damp_velocity << " " << target_dv << '\n';
	    *log_os_p << "# damp_dep " << damp_depth << " " << target_dd << '\n';
	}
	*log_os_p << "# ndata " << ndata << '\n';
	if (gravity){
	    *log_os_p << "# grav_data " << obs_grav.size() << " " << weight_grav << '\n';
	}
	*log_os_p << "# nnodes " << nnodev << " " << nnoded << " " << refl_weight << '\n';
	*log_os_p << "# LSQR " << LSQR_ATOL << endl;
    }

    // construct index mappers
    for (int i=1; i<=nnodev; i++){
	tmp_node(i) = i;
	tmp_nodev(i) = i;
    }
    if (nrefl>0){
	for (int i=1; i<=nnoded; i++){
	    int i_total = i+nnodev;
	    tmp_node(i+nnodev) = i_total;
	    tmp_nodedc(i) = i_total;
	    tmp_nodedr(i) = i;
	}
    }

    Array1d<int> kstart;
    Array1d<double> filtered_model, dx_vec, dxr_vec, dxn_vec;
    if (do_filter){
	filtered_model.resize(dmodel_total.size());
	dx_vec.resize(dmodel_total.size());
	dxr_vec.resize(dmodel_total.size());
	dxn_vec.resize(dmodel_total.size());
	kstart.resize(nx);
	for (int i=1; i<=nx; i++){
	    Point2d p = smesh.nodePos(smesh.nodeIndex(i,1));
	    double boundz = uboundp->z(p.x());
	    kstart(i) = 1;
	    for (int k=1; k<=nz; k++){
		if (smesh.nodePos(smesh.nodeIndex(i,k)).y()>boundz) break;
		kstart(i) = k;
	    }
	}
    }
    
    // pre-allocate temporary arrays
    Array1d<double> tmp_modelv(modelv.size());
    Array1d<double> tmp_modeld(modeld.size());
    if (jumping) dmodel_total_sum.resize(dmodel_total.size());

    bool isFinal=false;
    int iter=1;
    while(iter<=niter){
	if (verbose_level>=0) cerr << "TomographicInversion2d::iter="
				   << iter << "(" << niter << ")\n";
	if (iter==niter) isFinal=true;
	double graph_time=0.0, bend_time=0.0;
    
	smesh.get(modelv);
	if (nrefl==1) reflp->get(modeld);
	if (iter==1 || !jumping){ // set model scaling vectors
	    mvscale = modelv;
	    if (nrefl==1) mdscale = modeld;
	}

	reset_kernel(); nodev_hit=0;
	int idata=1;
	for (int isrc=1; isrc<=src.size(); isrc++){
	    if (verbose_level>=0){
		cerr << "isrc=" << isrc << " nrec=" << rcv(isrc).size()
		     << " ";
	    }
	    ofstream *tres_os_p, *ray_os_p;
	    if (printTransient || (printFinal && isFinal)){
		if (out_level >= 1){
		    char transfn[MaxStr];
		    sprintf(transfn, "%s.tres.%d.%d", out_root, iter, isrc);
		    tres_os_p = new ofstream(transfn);
		}
		if (out_level >= 2){
		    char transfn[MaxStr];
		    sprintf(transfn, "%s.ray.%d.%d", out_root, iter, isrc);
		    ray_os_p = new ofstream(transfn);
		}
	    }

	    // limit range first
	    double xmin=smesh.xmax();
	    double xmax=smesh.xmin();
	    for (int ircv=1; ircv<=rcv(isrc).size(); ircv++){
		double x = rcv(isrc)(ircv).x();
		if (x < xmin) xmin = x;
		if (x > xmax) xmax = x;
	    }
	    double srcx=src(isrc).x();
	    if (srcx < xmin) xmin = srcx;
	    if (srcx > xmax) xmax = srcx;
	    graph.limitRange(xmin,xmax);
	    if (verbose_level>=1){
		cerr << "xrange=(" << xmin << "," << xmax << ") ";
	    }

	    //
	    if (do_full_refl){
		double pwater = 1.0/1.5;
		tmp_modelv = modelv;
		Array1d<int> irefl(nx);
		for (int i=1; i<=nx; i++){
		    Point2d p = smesh.nodePos(smesh.nodeIndex(i,1));
		    double reflz = reflp->z(p.x());
		    irefl(i) = nz;
		    for (int k=1; k<=nz; k++){
			if (smesh.nodePos(smesh.nodeIndex(i,k)).y()>reflz){
			    irefl(i) = k;
			    break;
			}
		    }
		}
		for (int i=1; i<=nx; i++){
		    for (int k=2; k<=nz; k++){
			if (k==irefl(i)){
			    // this is to make the velocity at the node just below the reflector to
			    // be the same sa the velocity just above. If I use pwater for this node,
			    // there would be a unwanted low velocity gradient surrounding the reflector.
			    tmp_modelv(smesh.nodeIndex(i,k)) = tmp_modelv(smesh.nodeIndex(i,k-1));
			}else if (k>irefl(i)){
			    tmp_modelv(smesh.nodeIndex(i,k)) = pwater;
			}
		    }
		}
	    }
	    
	    clock_t start_t=clock();
	    graph.solve(src(isrc));
	    bool is_refl_solved = false;
	    clock_t end_t=clock();
	    graph_time += end_t-start_t;

	    start_t = clock();
	    double graph_refl_time=0;
	    for (int ircv=1; ircv<=rcv(isrc).size(); ircv++){
		Point2d r = rcv(isrc)(ircv);
		double orig_t, final_t;
		int iterbend;
		int icode = raytype(isrc)(ircv);
		if (icode == 0){ // refraction
		    if (smesh.inWater(r)){
			if (verbose_level>=0) cerr << "*";
			int i0, i1;
			graph.pickPathThruWater(r,path,i0,i1);
			start_i.resize(1); end_i.resize(1); interp.resize(1);
			start_i(1) = i0; end_i(1) = i1; interp(1) = bathyp;
			iterbend=bend.refine(path,orig_t,final_t,
					     start_i, end_i, interp);
			if (verbose_level>=1) cerr << "(" << iterbend << ")";
		    }else{
			if (verbose_level>=0) cerr << ".";
			graph.pickPath(rcv(isrc)(ircv),path);
			iterbend=bend.refine(path,orig_t,final_t,bend_nfac);
			if (verbose_level>=1) cerr << "(" << iterbend << ")";
		    }
		    if (iterbend<0){
			cerr << "TomographicInversion2d::iteration = " << iter << '\n';
			cerr << "TomographicInversion2d::too many iterations required\n";
			cerr << "TomographicInversion2d::for bending refinement at (s,r)="
			     << isrc << "," << ircv << '\n';
			exit(1);
		    }
		    add_kernel(idata,path);
		}else if (icode == 1){ // reflection
		    if (nrefl==0){
			error("TomographicInversion2d:: reflector not specified.");
		    }
		    if (do_full_refl){
			// temporarily replace sub-reflector velocity field
			// with water velocity
			smesh.set(tmp_modelv);
		    }
		    if (!is_refl_solved){
			clock_t start_t=clock();
			graph.solve_refl(src(isrc), *reflp);
			clock_t end_t=clock();
			graph_refl_time = end_t-start_t;
			graph_time += graph_refl_time;
			is_refl_solved = true;
		    }
		    int ir0, ir1;
		    if (smesh.inWater(r)){
			if (verbose_level>=0) cerr << "#";
			int i0, i1;
			graph.pickReflPathThruWater(r,path,i0,i1,ir0,ir1);
			start_i.resize(2); end_i.resize(2); interp.resize(2);
			start_i(1) = i0; end_i(1) = i1; interp(1) = bathyp;
			start_i(2) = ir0; end_i(2) = ir1; interp(2) = reflp;
			int iterbend=bend.refine(path,orig_t,final_t,
						 start_i, end_i, interp);
			if (verbose_level>=1) cerr << "(" << iterbend << ")";
		    }else{
			if (verbose_level>=0) cerr << "+";
			graph.pickReflPath(r,path,ir0,ir1);
			start_i.resize(1); end_i.resize(1); interp.resize(1);
			start_i(1) = ir0; end_i(1) = ir1; interp(1) = reflp;
			int iterbend=bend.refine(path,orig_t,final_t,
						 start_i, end_i, interp);
			if (verbose_level>=1) cerr << "(" << iterbend << ")";
		    }
		    if (do_full_refl){
			// revert to the original smesh
			smesh.set(modelv);
		    }
		    add_kernel_refl(idata,path,ir0,ir1);
		}else{
		    error("TomographicInversion2d:: illegal raycode detected.");
		}
		res_ttime(isrc)(ircv)=obs_ttime(isrc)(ircv)-final_t;
		idata++;

		if (printTransient || (printFinal && isFinal)){
		    if (out_level>=1){
			*tres_os_p << rcv(isrc)(ircv).x() << " "
				   << res_ttime(isrc)(ircv) << '\n';
		    }
		    if (out_level>=2){
			*ray_os_p << ">\n";
			printCurve(*ray_os_p, path, betasp);
		    }
		}
	    }
	    if (printTransient || (printFinal && isFinal)){
		if (out_level>=1){ tres_os_p->flush(); delete tres_os_p;}
		if (out_level>=2){ ray_os_p->flush(); delete ray_os_p;}
	    }
		    
	    if (verbose_level>=0) cerr << '\n';
	    end_t = clock();
	    bend_time += end_t-start_t-graph_refl_time;
	}
	graph_time /= CLOCKS_PER_SEC;
	bend_time /= CLOCKS_PER_SEC;

	// construct data vector
	idata=1;
	rms_tres[0]=rms_tres[1]=0;
	init_chi[0]=init_chi[1]=0;
	ndata_in[0]=ndata_in[1]=0;
	ndata_valid=0;
	tmp_data=0;
	for (int isrc=1; isrc<=src.size(); isrc++){
	    for (int ircv=1; ircv<=rcv(isrc).size(); ircv++){
		double res=res_ttime(isrc)(ircv);
		data_vec(idata) = res;
		double res2=res*r_dt_vec(idata);
		double res22=res2*res2;
		int icode = raytype(isrc)(ircv);
		rms_tres[icode] += res*res;
		init_chi[icode] += res22;
		++ndata_in[icode];
		tmp_data(idata) = ++ndata_valid;
		idata++;
	    }
	}
	rms_tres_total = sqrt((rms_tres[0]+rms_tres[1])/ndata_valid);
	init_chi_total = (init_chi[0]+init_chi[1])/ndata_valid;
	for (int i=0; i<=1; i++){
	    rms_tres[i] = ndata_in[i]>0 ? sqrt(rms_tres[i]/ndata_in[i]) : 0.0;
	    init_chi[i] = ndata_in[i]>0 ? init_chi[i]/ndata_in[i] : 0.0;
	}
	if (init_chi_total<target_chisq) isFinal=true;

	// rescale kernel A and data vector
	// note: averaging matrices are scaled upon their construction.
	for (int i=1; i<=ndata; i++){
	    double data_wt=r_dt_vec(i);
	    data_vec(i) *= data_wt;
	    for (mapIterator p=A(i).begin(); p!=A(i).end(); p++){
		int j = p->first;
		double m;
		if (j<=nnodev){
		    m = mvscale(j);
		}else{
		    m = mdscale(j-nnodev)*refl_weight;
		}
		p->second *= m*data_wt;
	    }
	}

	// construct total kernel matrix and data vecter, and solve Ax=b
	if (smooth_velocity) calc_averaging_matrix();
	if (nrefl>0 && smooth_depth) calc_refl_averaging_matrix();
	if (damp_velocity) calc_damping_matrix();
	if (nrefl>0 && damp_depth) calc_refl_damping_matrix();

	// (optional) joint inversion with gravity anomalies
	if (gravity){
	    if (verbose_level>=0) cerr << "calculating residual gravity anomalies...";
	    ginv->calcGravity(grav_z0, grav_x, res_grav);
//	    for (int i=1; i<=res_grav.size(); i++){
//		cerr << grav_x(i) << " " << res_grav(i) << '\n';
//	    }
	    if (verbose_level>=0) cerr << "done.\n";

	    if (verbose_level>=0) cerr << "calculating gravity kernel...";
	    rms_grav = 0.0;
	    for (int i=1; i<=ngravdata; i++){
		if (verbose_level>=1) cerr << i << " ";
		// residual gravity anomalies
		double rg = obs_grav(i)-res_grav(i);
		res_grav(i) = rg*weight_grav;
		rms_grav += rg*rg;

		// gravity kernels
		ginv->calcGravityKernel(grav_z0, grav_x(i), B(i));
		// model normalization
		for (mapIterator p=B(i).begin(); p!=B(i).end(); p++){
		    int j = p->first;
		    double m;
		    if (j<=nnodev){
			m = mvscale(j);
		    }else{
			m = mdscale(j-nnodev)*refl_weight;
		    }
		    p->second *= m*weight_grav;
		}
	    }
	    rms_grav = sqrt(rms_grav/ngravdata);
	    if (verbose_level>=0) cerr << "done.\n";

	    if (printTransient || (printFinal && isFinal)){
		char transfn[MaxStr];
		sprintf(transfn, "%s.rgrav.%d", out_root, iter);
		ofstream os(transfn);
		for (int i=1; i<=ngravdata; i++){
		    os << grav_x(i) << " " << res_grav(i) << " " << res_grav(i)/weight_grav << '\n';
		}
	    }
	}
	
	int iset=0;
	for (double tmp_wsv=wsv_min; tmp_wsv<=wsv_max; tmp_wsv+=dwsv){
	    for (double tmp_wsd=wsd_min; tmp_wsd<=wsd_max; tmp_wsd+=dwsd){
		iset++;

		weight_s_v = logscale_vel ? pow(10.0,tmp_wsv) : tmp_wsv;
		weight_s_d = logscale_dep ? pow(10.0,tmp_wsd) : tmp_wsd;
		
		double wdv,wdd;
		int nlsqr=0, lsqr_iter=0;
		time_t start_t = time(NULL);
		if (damping_is_fixed){
		    fixed_damping(lsqr_iter,nlsqr,wdv,wdd);
		}else{
		    auto_damping(lsqr_iter,nlsqr,wdv,wdd);
		}
		double lsqr_time = difftime(time(NULL),start_t);

		if (do_filter){
		    if (verbose_level>=0) cerr << "filtering velocity perturbation...";
		    for (int i=1; i<=dmodel_total.size(); i++) filtered_model(i) = dmodel_total(i);
		    for (int i=1; i<=nx; i++){
			for (int k=kstart(i); k<=nz; k++){
			    int inode = smesh.nodeIndex(i,k);
			    Point2d p = smesh.nodePos(inode);
			    double Lh, Lv;
			    corr_vel_p->at(p,Lh,Lv);

			    double Lh2 = Lh*Lh;
			    double Lv2 = Lv*Lv;
			    double sum=0.0, beta_sum=0.0;
			    for (int ii=1; ii<=nx; ii++){
				int jnode = smesh.nodeIndex(ii,k);
				double dx = smesh.nodePos(jnode).x()-p.x();
				if (abs(dx)<=Lh){
				    double xexp = exp(-dx*dx/Lh2);
				    for (int kk=kstart(ii); kk<=nz; kk++){
					int knode = smesh.nodeIndex(ii,kk);
					double dz = smesh.nodePos(knode).y()-p.y();
					if (abs(dz)<=Lv){
					    double beta = xexp*exp(-dz*dz/Lv2);
					    sum += beta*dmodel_total(knode);
					    beta_sum += beta;
					}
				    }
				}
			    }
			    filtered_model(inode) = sum/beta_sum;
			}
		    }
		    if (verbose_level>=0) cerr << "done.\n";
		    
		    // conservative filtering (Deal and Nolet, GJI, 1996)
		    if (verbose_level>=0) cerr << "calculating a nullspace shuttle...";
		    dx_vec=0.0;
		    for (int i=1; i<=nnodev; i++) dx_vec(i) = filtered_model(i)-dmodel_total(i);
		    Array1d<double> Adm(ndata_valid);
		    SparseRectangular sparseA(A,tmp_data,tmp_node,nnode_total);
		    sparseA.Ax(dx_vec,Adm);

		    int lsqr_itermax=itermax_LSQR;
		    double chi;
		    dxr_vec=0.0;
		    iterativeSolver_LSQR(sparseA,Adm,dxr_vec,LSQR_ATOL,lsqr_itermax,chi);

		    // note: correction_vec for depth nodes should be zero, so I don't use them.
		    double orig_norm=0, new_norm=0, iprod=0;
		    for (int i=1; i<=nnodev; i++){
			orig_norm += filtered_model(i)*filtered_model(i);
			dmodel_total(i) = filtered_model(i)-dxr_vec(i);
			new_norm += dmodel_total(i)*dmodel_total(i);
			dxn_vec(i) = dx_vec(i)-dxr_vec(i);
			iprod += dxr_vec(i)*dxn_vec(i);
		    }
		    if (verbose_level>=0) cerr << "done.\n";
		    if (printLog){
			*log_os_p << "# a posteriori filter check: "
				  << sqrt(orig_norm/nnodev)
				  << " " << sqrt(new_norm/nnodev) 
				  << " " << iprod << endl;
		    }
		}
		
		// take stats
		double pred_chi = calc_chi();
		double dv_norm = calc_ave_dmv();
		double dd_norm = calc_ave_dmd();
		double Lmvh=-1, Lmvv=-1, Lmd=-1;
		calc_Lm(Lmvh,Lmvv,Lmd);

		if (jumping) dmodel_total_sum += dmodel_total;

		// scale to slowness perturbation
		tmp_modelv = modelv;
		tmp_modeld = modeld;
		for (int i=1; i<=nnodev; i++){
		    tmp_modelv(i) += mvscale(i)*dmodel_total(i);
		}
		smesh.set(tmp_modelv);
		if (nrefl>0){
		    for (int i=1; i<=nnoded; i++){
			int tmp_i = tmp_nodedc(i);
			tmp_modeld(i) += mdscale(i)*dmodel_total(tmp_i)*refl_weight;
		    }
		    reflp->set(tmp_modeld);
		}

		if (printTransient || (printFinal && isFinal)){
		    char transfn[MaxStr];
		    sprintf(transfn, "%s.smesh.%d.%d", out_root, iter, iset);
		    ofstream os(transfn);
		    smesh.outMesh(os);

		    if (nrefl>0){
			char transfn2[MaxStr];
			sprintf(transfn2, "%s.refl.%d.%d", out_root, iter, iset);
			ofstream os2(transfn2);
			os2 << *reflp;
		    }
		}
		if (printLog){
		    *log_os_p << iter << " " << iset << " "
			      << ndata-ndata_valid << " "
			      << rms_tres_total << " " << init_chi_total << " "
			      << ndata_in[0] << " " << rms_tres[0] << " " << init_chi[0] << " "
			      << ndata_in[1] << " " << rms_tres[1] << " " << init_chi[1] << " "
			      << graph_time << " " << bend_time << " " 
			      << weight_s_v << " " << weight_s_d << " "
			      << wdv << " " << wdd << " "
			      << nlsqr << " " << lsqr_iter << " " << lsqr_time << " "
			      << pred_chi << " " << dv_norm << " " << dd_norm << " "
			      << Lmvh << " " << Lmvv << " " << Lmd;
		    if (gravity){
			*log_os_p << " " << rms_grav;
		    }
		    *log_os_p << endl;
		}
	    }
	}

	if (isFinal) break;
	iter++;
    }

    // output DWS for the last iteration
    if (out_mask){
	dws.resize(nnodev);
	dws=0.0;
	typedef map<int,double>::iterator mapIterator;
	for (int i=1; i<=A.size(); i++){
	    for (mapIterator p=A(i).begin(); p!=A(i).end(); p++){
		int inode=p->first;
		if (inode<=dws.size()) dws(inode) += p->second;
	    }
	}
	smesh.printMaskGrid(*vmesh_os_p, dws);
    }
    if (out_grav_dws){
	dws.resize(nnodev);
	dws=0.0;
	typedef map<int,double>::iterator mapIterator;
	for (int i=1; i<=B.size(); i++){
	    for (mapIterator p=B(i).begin(); p!=B(i).end(); p++){
		int inode=p->first;
		if (inode<=dws.size()) dws(inode) += p->second/(weight_grav*mvscale(inode));
	    }
	}
	smesh.printMaskGrid(*grav_dws_osp, dws);
    }
}

void
TomographicInversion2d::addRefl(Interface2d* intfp)
{
    reflp = intfp;
    nrefl = 1;
    nnoded = reflp->numNodes();
    Rd.resize(nnoded);
    Td.resize(nnoded);
    modeld.resize(nnoded);
    mdscale.resize(nnoded);
    dmodel_total.resize(nnodev+nnoded);
    nnode_total = nnodev+nnoded;
    tmp_node.resize(nnodev+nnoded);
    tmp_nodedc.resize(nnoded);
    tmp_nodedr.resize(nnoded);
}

void TomographicInversion2d::doFullRefl()
{
    do_full_refl = true;
    graph.do_refl_downward();
}

void TomographicInversion2d::setReflWeight(double x)
{
    if (x>0){
	refl_weight = x;
    }else{
	cerr << "TomographicInversion2d::setReflWeight - non-positive value ignored.\n";
    }
}

void TomographicInversion2d::SmoothVelocity(const char* fn,
					    double start, double end, double d,
					    bool scale)
{
    smooth_velocity=true;
    wsv_min=start; wsv_max=end; dwsv=d;
    logscale_vel=scale;
    corr_vel_p = new CorrelationLength2d(fn);
}

void TomographicInversion2d::applyFilter(const char *fn)
{
    do_filter=true;
    if (fn[0] != '\0'){
	uboundp = new Interface2d(fn);
    }else{
	uboundp = new Interface2d(smesh); // use bathymetry as upperbound
    }
}

void TomographicInversion2d::SmoothDepth(const char* fn,
					 double start, double end, double d,
					 bool scale)
{
    SmoothDepth(start,end,d,scale);
    corr_dep_p = new CorrelationLength1d(fn);
}

void TomographicInversion2d::SmoothDepth(double start, double end, double d,
					 bool scale)
{
    smooth_depth=true;
    wsd_min=start; wsd_max=end; dwsd=d;
    logscale_dep=scale;
}

void TomographicInversion2d::DampVelocity(double a)
{
    damp_velocity = true;
    target_dv = a/100.0; // % -> fraction
}

void TomographicInversion2d::DampDepth(double a)
{
    damp_depth = true;
    target_dd = a/100.0; // % -> fraction 
}

void TomographicInversion2d::FixDamping(double v, double d)
{
    damping_is_fixed = true;
    damp_velocity = true;
    damp_depth = true;
    weight_d_v = v;
    weight_d_d = d;
}

void TomographicInversion2d::Squeezing(const char* fn)
{
    damping_wt_p = new DampingWeight2d(fn);
    do_squeezing = true;
}

void TomographicInversion2d::doJumping()
{
    jumping = true;
}

void TomographicInversion2d::reset_kernel()
{
    typedef map<int,double>::iterator mapIterator;
    
    for (int i=1; i<=A.size(); i++){
	for (mapIterator p=A(i).begin(); p!=A(i).end(); p++) A(i).erase(p);
    }
    for (int i=1; i<=Rv_h.size(); i++){
	for (mapIterator p=Rv_h(i).begin(); p!=Rv_h(i).end(); p++) Rv_h(i).erase(p);
    }
    for (int i=1; i<=Rv_v.size(); i++){
	for (mapIterator p=Rv_v(i).begin(); p!=Rv_v(i).end(); p++) Rv_v(i).erase(p);
    }
    for (int i=1; i<=Tv.size(); i++){
	for (mapIterator p=Tv(i).begin(); p!=Tv(i).end(); p++) Tv(i).erase(p);
    }
    if (nrefl>0){
	for (int i=1; i<=Rd.size(); i++){
	    for (mapIterator p=Rd(i).begin(); p!=Rd(i).end(); p++) Rd(i).erase(p);
	}
	for (int i=1; i<=Td.size(); i++){
	    for (mapIterator p=Td(i).begin(); p!=Td(i).end(); p++) Td(i).erase(p);
	}
    }
    if (gravity){
	for (int i=1; i<=B.size(); i++){
	    for (mapIterator p=B(i).begin(); p!=B(i).end(); p++) B(i).erase(p);
	}
    }
}

void TomographicInversion2d::add_kernel(int idata, const Array1d<Point2d>& cur_path)
{
    map<int,double>& A_i = A(idata);

    int np = cur_path.size();
    int nintp = betasp.numIntp();
    makeBSpoints(cur_path,pp);

    double path_len=0.0;
    Index2d guess_index = smesh.nodeIndex(smesh.nearest(*pp(1)));
    for (int i=1; i<=np+1; i++){
	int j1=i;
	int j2=i+1;
	int j3=i+2;
	int j4=i+3;
	betasp.interpolate(*pp(j1),*pp(j2),*pp(j3),*pp(j4),Q);
	
	for (int j=2; j<=nintp; j++){
	    Point2d midp = 0.5*(Q(j-1)+Q(j));
	    double dist= Q(j).distance(Q(j-1));
	    path_len+=dist;

	    int jUL, jLL, jLR, jUR;
	    double r, s, rr, ss;
	    int icell=smesh.locateInCell(midp,guess_index,
					 jUL,jLL,jLR,jUR,r,s,rr,ss);
	    if (icell>0){
		A_i[jUL] += rr*ss*dist;
		A_i[jLL] += rr*s*dist;
		A_i[jLR] += r*s*dist;
		A_i[jUR] += r*ss*dist;

		nodev_hit(jUL)++;
		nodev_hit(jLL)++;
		nodev_hit(jLR)++;
		nodev_hit(jUR)++;
	    }
	}
    }
    path_length(idata)=path_len;
}

void TomographicInversion2d::add_kernel_refl(int idata, const Array1d<Point2d>& cur_path,
					     int ir0, int ir1)
{
    map<int,double>& A_i = A(idata);

    int jL, jR;
    reflp->locateInSegment(path(ir0).x(),jL,jR);
    if (jL==jR){
	cerr << "TomographicInversion2d::add_kernel_refl - bottoming point out of bounds\n";
	return; // ignore this ray info
    }
    double x1 = reflp->x(jL);
    double x2 = reflp->x(jR);
    double z1 = reflp->z(x1);
    double z2 = reflp->z(x2);
    double x = path(ir0).x();

    double dx = x2-x1;
    double dz = z2-z1;
    double cos_alpha = dx/sqrt(dx*dx+dz*dz);

    double p1 = smesh.at(path(ir0));

    Point2d a = path(ir0-1)-path(ir0); a /= a.norm();
    Point2d b = path(ir1+1)-path(ir1); b /= b.norm();
    Point2d c = a+b; c /= c.norm();
    double cos_theta = a.inner_product(c);
    
    double com_factor = 2.0*cos_theta*cos_alpha*p1/dx;
    A_i[nnodev+jL] = com_factor*(x2-x);
    A_i[nnodev+jR] = com_factor*(x-x1);

    add_kernel(idata,path);
}    

void TomographicInversion2d::calc_averaging_matrix()
{
    for (int i=1; i<=nx; i++){
	for (int k=1; k<=nz; k++){
	    int inode = smesh.nodeIndex(i,k);
	    Point2d p = smesh.nodePos(inode);
	    double Lh, Lv;
	    corr_vel_p->at(p,Lh,Lv);

	    double Lh2 = Lh*Lh;
	    double Lv2 = Lv*Lv;
	    double beta_sum = 0.0;
	    for (int ii=1; ii<=nx; ii++){
		int jnode = smesh.nodeIndex(ii,k);
		if (jnode!=inode){
		    double dx = smesh.nodePos(jnode).x()-p.x();
		    if (abs(dx)<=Lh){
			double dxL2 = dx*dx/Lh2;
			double beta = exp(-dxL2);
			Rv_h(inode)[jnode] = beta*mvscale(jnode);
			beta_sum += beta;
		    }
		}
	    }
	    Rv_h(inode)[inode] = -beta_sum*mvscale(inode);

	    beta_sum = 0.0;
	    for (int kk=1; kk<=nz; kk++){
		int jnode = smesh.nodeIndex(i,kk);
		if (jnode!=inode){
		    double dz = smesh.nodePos(jnode).y()-p.y();
		    if (abs(dz)<=Lv){
			double dzL2 = dz*dz/Lv2;
			double beta = exp(-dzL2);
			Rv_v(inode)[jnode] = beta*mvscale(jnode);
			beta_sum += beta;
		    }
		}
	    }
	    Rv_v(inode)[inode] = -beta_sum*mvscale(inode);
	}
    }
}

void TomographicInversion2d::calc_refl_averaging_matrix()
{
    for (int i=1; i<=nnoded; i++){
	double x = reflp->x(i);
	double Lh;
	if (corr_dep_p != 0){
	    Lh = corr_dep_p->at(x);
	}else if (corr_vel_p != 0){
	    double Lv;
	    corr_vel_p->at(Point2d(x,reflp->z(x)),Lh,Lv);
	}else{
	    error("TomographicInversion2d::calc_refl_averaging_matrix - no correlation length available");
	}
	double Lh2 = Lh*Lh;
	double beta_sum = 0.0;
	for (int ii=1; ii<=nnoded; ii++){
	    if (ii!=i){
		double dx = reflp->x(ii)-x;
		if (abs(dx)<=Lh){
		    double beta = exp(-dx*dx/Lh2);
		    Rd(i)[ii] = beta*mdscale(ii)*refl_weight;
		    beta_sum += beta;
		}
	    }
	}
	Rd(i)[i] = -beta_sum*mdscale(i)*refl_weight;
    }
}

void TomographicInversion2d::calc_damping_matrix()
{
    Array2d<double> local_T(4,4);
    Array1d<int> j(4);
    
    for (int i=1; i<=smesh.numCells(); i++){
	smesh.cellNodes(i, j(1), j(2), j(3), j(4));
	smesh.cellNormKernel(i, local_T);

	if (do_squeezing){
	     Point2d p = 0.25*(smesh.nodePos(j(1))+smesh.nodePos(j(2))
			       +smesh.nodePos(j(3))+smesh.nodePos(j(4)));
	     double w;
	     damping_wt_p->at(p,w);
	     local_T *= w;
	}
	add_global(local_T, j, Tv);
    }
}

void TomographicInversion2d::calc_refl_damping_matrix()
{
    double dx = reflp->x(2)-reflp->x(1); // assumes uniform interval
    double fac = sqrt(dx);

    for (int i=1; i<=nnoded; i++){
	Td(i)[i] = fac;
    }
}

void TomographicInversion2d::add_global(const Array2d<double>& a,
					const Array1d<int>& j, sparseMat& global)
{
    for (int m=1; m<=4; m++){
	for (int n=1; n<=4; n++){
	    global(j(m))[j(n)] += a(m,n);
	}
    }
}

double TomographicInversion2d::calc_ave_dmv()
{
    double dm_norm=0.0;
    for (int i=1; i<=nnodev; i++){
	dm_norm += dmodel_total(i)*dmodel_total(i);
    }
    return sqrt(dm_norm/nnodev);
}

double TomographicInversion2d::calc_ave_dmd()
{
    if (nrefl==0) return 0.0;
	
    double dm_norm=0.0;
    for (int i=1+nnodev; i<=nnoded+nnodev; i++){
	dm_norm += dmodel_total(i)*dmodel_total(i);
    }
    return sqrt(dm_norm/nnoded)*refl_weight;
}

void TomographicInversion2d::calc_Lm(double& lmh, double& lmv, double& lmd)
{
    if (smooth_velocity){
	Array1d<double> Rvm(nnodev), dmv_vec(nnodev);
	SparseRectangular sparseRv_h(Rv_h,tmp_nodev,tmp_nodev);
	for (int i=1; i<=nnodev; i++) dmv_vec(i) = dmodel_total(i);

	sparseRv_h.Ax(dmv_vec,Rvm);
	double val=0.0;
	for (int i=1; i<=nnodev; i++) val += Rvm(i)*Rvm(i);
	lmh = sqrt(val/nnodev);

	SparseRectangular sparseRv_v(Rv_v,tmp_nodev,tmp_nodev);
	sparseRv_v.Ax(dmv_vec,Rvm);
	val = 0.0;
	for (int i=1; i<=nnodev; i++) val += Rvm(i)*Rvm(i);
	lmv = sqrt(val/nnodev);
    }

    if (nrefl>0 && smooth_depth){
	Array1d<double> Rdm(nnoded), dmd_vec(nnoded);
	SparseRectangular sparseRd(Rd,tmp_nodedr,tmp_nodedr);
	for (int i=1+nnodev, j=1; i<=nnoded+nnodev; i++) dmd_vec(j++) = dmodel_total(i);

	sparseRd.Ax(dmd_vec,Rdm);
	double val = 0.0;
	for (int i=1; i<=nnoded; i++) val += Rdm(i)*Rdm(i);
	lmd = sqrt(val/nnoded);
    }
}
 
double TomographicInversion2d::calc_chi()
{
    Array1d<double> Adm(ndata_valid);
    SparseRectangular sparseA(A,tmp_data,tmp_node,nnode_total);
    sparseA.Ax(dmodel_total,Adm);
    double val=0.0;
    for (int i=1; i<=ndata; i++){
	int j=tmp_data(i);
	if (j>0){
	    double res=Adm(j)-data_vec(i);
	    val += res*res;
	}
    }
    return val/ndata_valid;
}

int TomographicInversion2d::_solve(bool sv, double wsv, bool sd, double wsd,
				   bool dv, double wdv, bool dd, double wdd)
{
    // construct total kernel
    Array1d<const sparseMat*> As;
    Array1d<SparseMatAux> matspec;
    As.push_back(&A);
    matspec.push_back(SparseMatAux(1.0,&tmp_data,&tmp_node));
    int ndata_plus=ndata_valid;
    if (gravity){
	As.push_back(&B);
	matspec.push_back(SparseMatAux(1.0,&tmp_gravdata,&tmp_node));
	ndata_plus += ngravdata;
    }
    if (sv){
	if (wsv<0) error("TomographicInversion2d::_solve - negative wsv");
	As.push_back(&Rv_h);
	matspec.push_back(SparseMatAux(wsv,&tmp_nodev,&tmp_nodev));
	ndata_plus += nnodev;
	As.push_back(&Rv_v);
	matspec.push_back(SparseMatAux(wsv,&tmp_nodev,&tmp_nodev));
	ndata_plus += nnodev;
    }
    if (nrefl>0 && sd){
	if (wsd<0) error("TomographicInversion2d::_solve - negative wsd");
	As.push_back(&Rd);
	matspec.push_back(SparseMatAux(wsd,&tmp_nodedr,&tmp_nodedc));
	ndata_plus += nnoded;
    }
    if (dv){
	if (wdv<0) error("TomographicInversion2d::_solve - negative wdv");
	As.push_back(&Tv);
	matspec.push_back(SparseMatAux(wdv,&tmp_nodev,&tmp_nodev));
	ndata_plus += nnodev;
    }
    if (nrefl>0 && dd){
	if (wdd<0) error("TomographicInversion2d::_solve - negative wdd");
	As.push_back(&Td);
	matspec.push_back(SparseMatAux(wdd,&tmp_nodedr,&tmp_nodedc));
	ndata_plus += nnoded;
    }
    SparseRectangular B(As,matspec,nnode_total);

    // construct total data vector
    total_data_vec.resize(ndata_plus);
    int idata=1;
    for (int i=1; i<=ndata; i++){
	if (tmp_data(i)>0){
	    total_data_vec(idata++) = data_vec(i);
	}
    }
    if (gravity){
	for (int i=1; i<=ngravdata; i++){
	    total_data_vec(idata++) = res_grav(i);
	}
    }
    if (sv){
	if (jumping){
	    Array1d<double> Rvm(nnodev), dmv_vec(nnodev);
	    SparseRectangular sparseRv_h(Rv_h,tmp_nodev,tmp_nodev);
	    for (int i=1; i<=nnodev; i++) dmv_vec(i) = dmodel_total_sum(i);
	    sparseRv_h.Ax(dmv_vec,Rvm);
	    for (int i=1; i<=nnodev; i++) total_data_vec(idata++) = -wsv*Rvm(i); // Lh

	    SparseRectangular sparseRv_v(Rv_v,tmp_nodev,tmp_nodev);
	    sparseRv_v.Ax(dmv_vec,Rvm);
	    for (int i=1; i<=nnodev; i++) total_data_vec(idata++) = -wsv*Rvm(i); // Lv
	}else{
	    for (int i=1; i<=nnodev; i++) total_data_vec(idata++) = 0.0; // Lh
	    for (int i=1; i<=nnodev; i++) total_data_vec(idata++) = 0.0; // Lv
	}
    }
    if (nrefl>0 && sd){
	if (jumping){
	    Array1d<double> Rdm(nnoded), dmd_vec(nnoded);
	    SparseRectangular sparseRd(Rd,tmp_nodedr,tmp_nodedr);
	    for (int i=1+nnodev, j=1; i<=nnoded+nnodev; i++) dmd_vec(j++) = dmodel_total_sum(i);
	    sparseRd.Ax(dmd_vec,Rdm);
	    for (int i=1; i<=nnoded; i++) total_data_vec(idata++) = -wsd*Rdm(i); // Ld
	}else{
	    for (int i=1; i<=nnoded; i++) total_data_vec(idata++) = 0.0;
	}
    }
    if (dv){
	if (jumping){
	    Array1d<double> Tvm(nnodev), dmv_vec(nnodev);
	    SparseRectangular sparseTv(Tv,tmp_nodev,tmp_nodev);
	    for (int i=1; i<=nnodev; i++) dmv_vec(i) = dmodel_total_sum(i);
	    sparseTv.Ax(dmv_vec,Tvm);
	    for (int i=1; i<=nnodev; i++) total_data_vec(idata++) = -wdv*Tvm(i); 
	}else{
	    for (int i=1; i<=nnodev; i++) total_data_vec(idata++) = 0.0;
	}
    }
    if (nrefl>0 && dd){
	if (jumping){
	    Array1d<double> Tdm(nnoded), dmd_vec(nnoded);
	    SparseRectangular sparseTd(Td,tmp_nodedr,tmp_nodedr);
	    for (int i=1+nnodev, j=1; i<=nnoded+nnodev; i++) dmd_vec(j++) = dmodel_total_sum(i);
	    sparseTd.Ax(dmd_vec,Tdm);
	    for (int i=1; i<=nnoded; i++) total_data_vec(idata++) = -wdd*Tdm(i); 
	}else{
	    for (int i=1; i<=nnoded; i++) total_data_vec(idata++) = 0.0;
	}
    }

    int lsqr_itermax=itermax_LSQR;
    double chi;

    int istop=iterativeSolver_LSQR(B,total_data_vec,dmodel_total,LSQR_ATOL,lsqr_itermax,chi);

    return lsqr_itermax;
}

void TomographicInversion2d::fixed_damping(int& iter, int& n, double& wdv, double& wdd)
{
    wdv = weight_d_v;
    wdd = weight_d_d;
    
    iter += _solve(smooth_velocity,weight_s_v,smooth_depth,weight_s_d,
		   true,weight_d_v,true,weight_d_d);
    n++;

    if (robust){
	if (verbose_level>=0){
	    cerr << "TomographicInversion2d:: check outliers... ";
	}
	int ndata_valid_orig=ndata_valid;
	removeOutliers();
	if (verbose_level>=0){
	    cerr << ndata_valid_orig-ndata_valid << " found\n";
	}
	if (ndata_valid < ndata_valid_orig){
	    if (verbose_level>=0){
		cerr << "TomographicInversion2d:: re-inverting without outliers...\n";
	    }
	    iter += _solve(smooth_velocity,weight_s_v,smooth_depth,weight_s_d,
			   true,weight_d_v,true,weight_d_d);
	    n++;
	}
    }
}

void TomographicInversion2d::auto_damping(int& iter, int& n, double& wdv, double& wdd)
{
    // check if damping is necessary
    iter += _solve(smooth_velocity,weight_s_v,smooth_depth,weight_s_d,
		   false,0.0,false,0.0);
    n++;

    if (robust){
	if (verbose_level>=0){
	    cerr << "TomographicInversion2d:: check outliers... ";
	}
	int ndata_valid_orig=ndata_valid;
	removeOutliers();
	if (verbose_level>=0){
	    cerr << ndata_valid_orig-ndata_valid << " found\n";
	}
	if (ndata_valid < ndata_valid_orig){
	    if (verbose_level>=0){
		cerr << "TomographicInversion2d:: re-inverting without outliers...\n";
	    }
	    iter += _solve(smooth_velocity,weight_s_v,smooth_depth,weight_s_d,
			   false,0.0,false,0.0);
	    n++;
	}
    }

    double ave_dmv0 = calc_ave_dmv();
    double ave_dmd0 = calc_ave_dmd();
    if (verbose_level>=0){
	cerr << "\t\tave_dm = " << ave_dmv0*100
	     << "%, " << ave_dmd0*100 << "% with no damping\n";
    }
    wdv = wdd = 1.0;
    if (ave_dmv0<=target_dv || !damp_velocity) wdv = 0.0;
    if (ave_dmd0<=target_dd || !damp_depth) wdd = 0.0; 
    if (wdv==0.0 && wdd==0.0) return;

    if (wdd>0.0) wdd = auto_damping_depth(wdv,iter,n);
    if (wdv>0.0) wdv = auto_damping_vel(wdd,iter,n);
}

double TomographicInversion2d::auto_damping_depth(double wdv, int& iter, int& n)
{
    // secant and bisection search
    if (verbose_level>=0) cerr << "\t\tsearching weight_d_depth...\n";
    const double wd_max = 1e7; // absolute bound
    double wd1 = 1.0; // initial guess
    double wd2 = 1e2;
    const double ddm = target_dd*0.3; // target accuracy = 30 %
    const int secant_itermax=10;
    double rts, xl, swap, dx, fl, f;

    iter += _solve(smooth_velocity,weight_s_v,smooth_depth,weight_s_d,
		   damp_velocity,wdv,true,wd1); n++;
    fl = calc_ave_dmd()-target_dd;
    if (verbose_level>=0){
	cerr << "\t\tave_dm = " << (fl+target_dd)*100 << "% at wd1("
	     << wd1 << ")\n";
    }
    iter += _solve(smooth_velocity,weight_s_v,smooth_depth,weight_s_d,
		   damp_velocity,wdv,true,wd2); n++;
    f = calc_ave_dmd()-target_dd;
    if (verbose_level>=0){
	cerr << "\t\tave_dm = " << (f+target_dd)*100 << "% at wd2("
	     << wd2 << ")\n";
    }

    if (abs(fl) < abs(f)){
	rts = wd1;
	xl = wd2;
	swap=fl; fl=f; f=swap;
    }else{
	xl = wd1; rts = wd2;
    }
    for (int j=1; j<=secant_itermax; j++){
	dx = (xl-rts)*f/(f-fl);
	if (rts+dx<=0){ // switch to bisection
	    xl=rts; fl=f;
	    rts *= 0.5;
	}else if(rts+dx>wd_max){
	    xl=rts; fl=f;
	    rts += (wd_max-rts)*0.5;
	}else{
	    xl=rts; fl=f;
	    rts+=dx;
	}
	iter += _solve(smooth_velocity,weight_s_v,smooth_depth,weight_s_d,
		       damp_velocity,wdv,true,rts); n++;
	f = calc_ave_dmd()-target_dd;
	if (verbose_level>=0){
	    cerr << "\t\tave_dm = " << (f+target_dd)*100
		 << "% at w_d of " << rts << '\n';
	}
	if (rts > wd_max) break;
	if (abs(f) < ddm) break;
    }
    return rts;
}

double TomographicInversion2d::auto_damping_vel(double wdd, int& iter, int& n)
{
    // secant and bisection search
    if (verbose_level>=0) cerr << "\t\tsearching weight_d_vel...\n";
    const double wd_max = 1e5; // absolute bound
    double wd1 = 1.0; // initial guess
    double wd2 = 1e2;
    const double ddm = target_dv*0.3; // target accuracy = 30 %
    const int secant_itermax=10;
    double rts, xl, swap, dx, fl, f;

    iter += _solve(smooth_velocity,weight_s_v,smooth_depth,weight_s_d,
		   true,wd1,damp_depth,wdd); n++;
    fl = calc_ave_dmv()-target_dv;
    if (verbose_level>=0){
	cerr << "\t\tave_dm = " << (fl+target_dv)*100 << "% at wd1("
	     << wd1 << ")\n";
    }
    iter += _solve(smooth_velocity,weight_s_v,smooth_depth,weight_s_d,
		   true,wd2,damp_depth,wdd); n++;
    f = calc_ave_dmv()-target_dv;
    if (verbose_level>=0){
	cerr << "\t\tave_dm = " << (f+target_dv)*100 << "% at wd2("
	     << wd2 << ")\n";
    }

    if (abs(fl) < abs(f)){
	rts = wd1;
	xl = wd2;
	swap=fl; fl=f; f=swap;
    }else{
	xl = wd1; rts = wd2;
    }
    for (int j=1; j<=secant_itermax; j++){
	dx = (xl-rts)*f/(f-fl);
	if (rts+dx<=0){ // switch to bisection
	    xl=rts; fl=f;
	    rts *= 0.5;
	}else{
	    xl=rts; fl=f;
	    rts+=dx;
	}
	iter += _solve(smooth_velocity,weight_s_v,smooth_depth,weight_s_d,
		       true,rts,damp_depth,wdd); n++;
	f = calc_ave_dmv()-target_dv;
	if (verbose_level>=0){
	    cerr << "\t\tave_dm = " << (f+target_dv)*100
		 << "% at w_d of " << rts << '\n';
	}
	if (rts > wd_max) break;
	if (abs(f) < ddm) break;
    }
    return rts;
}

