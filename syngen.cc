/*
 * syngen.cc - forward traveltime calculation 
 * 
 * Jun Korenaga, MIT/WHOI
 * January 1999
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <ctime>
#include "syngen.h"
#include "betaspline.h"
#include "traveltime.h"

SyntheticTraveltimeGenerator2d::SyntheticTraveltimeGenerator2d(
    SlownessMesh2d& m, const char* ifn,
    int xorder, int zorder, double clen, int nintp, double cg_tol, double br_tol)
    : smesh(m), graph(smesh,xorder,zorder),
      betasp(1,0,nintp), bend(smesh,betasp,cg_tol,br_tol),
      nrefl(0), do_full_refl(false),
      outray(false), use_clock(false), verbose_level(-1)
{
    if (clen>0.0) graph.refineIfLong(clen);

    start_i.reserve(2); end_i.reserve(2);
    interp.reserve(2);
    bathyp = new Interface2d(smesh);

    graph_only = false;
    read_file(ifn);
}

void SyntheticTraveltimeGenerator2d::graphOnly()
{
    graph_only = true;
}

void SyntheticTraveltimeGenerator2d::conduct()
{
    ofstream *rout_p;
    if (outray) rout_p = new ofstream(rayfn);
	
    Array1d<Point2d> path;
    const int max_np = int(2*sqrt(float(smesh.numNodes())));
    path.reserve(max_np);

    Array1d<double> modelv, tmp_modelv;
    if (do_full_refl){
	modelv.resize(smesh.numNodes());
	tmp_modelv.resize(smesh.numNodes());
	smesh.get(modelv);
	int nx(smesh.Nx()), nz(smesh.Nz());
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
    
    graph_time=bend_time=0.0;
    
    for (int isrc=1; isrc<=src.size(); isrc++){
	if (verbose_level>=0){
	    cerr << "isrc=" << isrc << " nrec=" << rcv(isrc).size()
		 << " ";
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

	clock_t start_t = clock();
	graph.solve(src(isrc));
	bool is_refl_solved=false;
	clock_t end_t=clock();
	graph_time += end_t-start_t;

	start_t = clock();
	double graph_refl_time=0;
	for (int ircv=1; ircv<=rcv(isrc).size(); ircv++){
	    Point2d r = rcv(isrc)(ircv);
	    double orig_t, final_t;
	    int icode = raytype(isrc)(ircv);
	    if (icode == 0){ // refraction
		if (smesh.inWater(r)){
		    if (verbose_level>=0) cerr << "*";
		    int i0, i1;
		    graph.pickPathThruWater(r,path,i0,i1);
		    start_i.resize(1); end_i.resize(1); interp.resize(1);
		    start_i(1) = i0; end_i(1) = i1; interp(1) = bathyp;
		    if (graph_only){
			final_t = calcTravelTime(smesh,path,betasp.numIntp());
		    }else{
			int iterbend=bend.refine(path,orig_t,final_t,
						 start_i, end_i, interp);
			if (verbose_level>=1) cerr << "(" << iterbend << ")";
		    }
		}else{
		    if (verbose_level>=0) cerr << ".";
		    graph.pickPath(r,path);
		    if (graph_only){
			final_t = calcTravelTime(smesh,path,betasp.numIntp());
		    }else{
			const int nfac=1;
			int iterbend=bend.refine(path,orig_t,final_t,nfac);
			if (verbose_level>=1) cerr << "(" << iterbend << ")";
		    }
		}
	    }else if (icode == 1){ // reflection
		if (nrefl==0){
		    error("SyntheticTraveltimeGenerator2d:: reflector not specified.");
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
		if (smesh.inWater(r)){
		    if (verbose_level>=0) cerr << "#";
		    int i0, i1, ir0, ir1;
		    graph.pickReflPathThruWater(r,path,i0,i1,ir0,ir1);
		    start_i.resize(2); end_i.resize(2); interp.resize(2);
		    start_i(1) = i0; end_i(1) = i1; interp(1) = bathyp;
		    start_i(2) = ir0; end_i(2) = ir1; interp(2) = reflp;
		    if (graph_only){
			final_t = calcTravelTime(smesh,path,betasp.numIntp());
		    }else{
			int iterbend=bend.refine(path,orig_t,final_t,
						 start_i, end_i, interp);
			if (verbose_level>=1) cerr << "(" << iterbend << ")";
		    }
		}else{
		    if (verbose_level>=0) cerr << "+";
		    int ir0, ir1;
		    graph.pickReflPath(r,path,ir0,ir1);
		    start_i.resize(1); end_i.resize(1); interp.resize(1);
		    start_i(1) = ir0; end_i(1) = ir1; interp(1) = reflp;
		    if (graph_only){
			final_t = calcTravelTime(smesh,path,betasp.numIntp());
		    }else{
			int iterbend=bend.refine(path,orig_t,final_t,
						 start_i, end_i, interp);
			if (verbose_level>=1) cerr << "(" << iterbend << ")";
		    }
		}
		if (do_full_refl){
		    // revert to the original smesh
		    smesh.set(modelv);
		}
	    }else{
		error("SyntheticTraveltimeGenerator2d:: illegal raycode detected.");
	    }
	    if (outray){
		*rout_p << ">\n";
		if (graph_only){
		    for (int i=1; i<=path.size(); i++){
			*rout_p << path(i).x() << " " << path(i).y() << '\n';
		    }
		}else{
		    printCurve(*rout_p,path,betasp);
		}
	    }
	    syn_ttime(isrc)(ircv) = final_t;
	}
	if (verbose_level>=0) cerr << '\n';
	end_t = clock();
	bend_time += end_t-start_t-graph_refl_time;
    }
    graph_time /= CLOCKS_PER_SEC;
    bend_time /= CLOCKS_PER_SEC;
    if (use_clock){
	ofstream os(timefn);
	os << "graph_time " << graph_time << " sec " 
	   << "bend_time " << bend_time << " sec \n";
    }
    if (outray) rout_p->flush(); // make sure all rays are printed out.
}

void SyntheticTraveltimeGenerator2d::outputRay(const char* fn)
{ outray=true; rayfn=fn; }

void SyntheticTraveltimeGenerator2d::useClock(const char* fn)
{ use_clock=true; timefn=fn; }

void SyntheticTraveltimeGenerator2d::setVerbose(int i){ verbose_level = i; }

void SyntheticTraveltimeGenerator2d::read_file(const char* ifn)
{
    ifstream in(ifn);
    if (!in){
	cerr << "SyntheticTraveltimeGenerator2d::cannot open " << ifn << "\n";
	exit(1);
    }

    int iline=0;

    string first;
    in >> first; iline++;
    if (!isdigit(*first.c_str()))
	error("SyntheticTraveltimeGenerator2d::first line should be nsrc");
    int nsrc = atoi(first.c_str());
    if (nsrc<=0) error("SyntheticTraveltimeGenerator2d::invalid nsrc");
    src.resize(nsrc); rcv.resize(nsrc);
    raytype.resize(nsrc);
    syn_ttime.resize(nsrc); obs_ttime.resize(nsrc); obs_dt.resize(nsrc);

    int isrc=0;
    while(in){
	char flag;
	double x, y;
	int nrcv;

	in >> flag >> x >> y >> nrcv; iline++;
	if (flag!='s'){
	    cerr << "SyntheticTraveltimeGenerator2d::bad input (s) at l."
		 << iline << '\n';
	    exit(1);
	}
	isrc++;
	src(isrc).set(x,y);

	rcv(isrc).resize(nrcv); raytype(isrc).resize(nrcv);
	syn_ttime(isrc).resize(nrcv);
	obs_ttime(isrc).resize(nrcv); obs_dt(isrc).resize(nrcv);
	for (int ircv=1; ircv<=nrcv; ircv++){
	    int n;
	    double ttime_val, dt_val;
	    in >> flag >> x >> y >> n >> ttime_val >> dt_val; iline++;
	    if (flag!='r'){
		cerr << "SyntheticTraveltimeGenerator2d::bad input (r) at l."
		     << iline << '\n';
		exit(1);
	    }
	    rcv(isrc)(ircv).set(x,y);
	    raytype(isrc)(ircv) = n;
	    obs_ttime(isrc)(ircv) = ttime_val;
	    obs_dt(isrc)(ircv) = dt_val;
	}
	if (isrc==nsrc) break;
    }
    if (isrc != nsrc) error("SyntheticTraveltimeGenerator2d::mismatch in nsrc");
}

void
SyntheticTraveltimeGenerator2d::readRefl(const char* fn)
{
    reflp = new Interface2d(fn);
    nrefl = 1;
}

void
SyntheticTraveltimeGenerator2d::doFullRefl()
{
    do_full_refl = true;
    graph.do_refl_downward();
}
    
void
SyntheticTraveltimeGenerator2d::printSource(ostream& os) const
{
    for (int isrc=1; isrc<=src.size(); isrc++){
	os << src(isrc).x() << " "
	   << src(isrc).y() << '\n';
    }
}

void
SyntheticTraveltimeGenerator2d::printSynTime(ostream& os, double vred) const
{
    for (int isrc=1; isrc<=src.size(); isrc++){
	os << ">\n";
	for (int ircv=1; ircv<=rcv(isrc).size(); ircv++){
	    double redtime = syn_ttime(isrc)(ircv);
	    if (vred!=0){
		redtime -= abs(rcv(isrc)(ircv).x()-src(isrc).x())/vred;
	    }
	    if (ircv>1 &&
		raytype(isrc)(ircv)!=raytype(isrc)(ircv-1)) os << ">\n";
	    os << rcv(isrc)(ircv).x() << " " << redtime << '\n';
	}
    }
}

void
SyntheticTraveltimeGenerator2d::printObsTime(ostream& os, double vred) const
{
    for (int isrc=1; isrc<=src.size(); isrc++){
	os << ">\n";
	for (int ircv=1; ircv<=rcv(isrc).size(); ircv++){
	    double redtime = obs_ttime(isrc)(ircv);
	    if (vred!=0){
		redtime -= abs(rcv(isrc)(ircv).x()-src(isrc).x())/vred;
	    }
	    if (ircv>1 &&
		raytype(isrc)(ircv)!=raytype(isrc)(ircv-1)) os << ">\n";
	    os << rcv(isrc)(ircv).x() << " " << redtime << " "
	       << obs_dt(isrc)(ircv) << '\n';
	}
    }
}

void SyntheticTraveltimeGenerator2d::printDiff(ostream& os, double& misfit, double& chisq) const
{
    misfit = 0.0;
    chisq = 0.0;
    int ndata=0;
    for (int isrc=1; isrc<=src.size(); isrc++){
	for (int ircv=1; ircv<=rcv(isrc).size(); ircv++){
	    double tdiff, tdiff2;
	    tdiff = obs_ttime(isrc)(ircv)-syn_ttime(isrc)(ircv);
	    tdiff2 = tdiff*tdiff;
	    misfit += tdiff2;
	    if (obs_dt(isrc)(ircv)==0.0){
		cerr << "SyntheticTraveltimeGenerator2d::printDiff - zero error found at "
		     << "isrc=" << isrc << ", ircv=" << ircv << '\n';
		exit(1);
	    }
	    double obsdt2 = obs_dt(isrc)(ircv)*obs_dt(isrc)(ircv);
	    chisq += tdiff2/obsdt2;
	    os << tdiff << " " << tdiff/obs_dt(isrc)(ircv) << '\n';
	    ndata++;
	}
    }
    misfit = sqrt(misfit/ndata);
    chisq = chisq/ndata;
    os << "# t_misfit " << misfit << '\n';
    os << "# chisq " << chisq << '\n';
}


ostream&
operator<<(ostream& out, const SyntheticTraveltimeGenerator2d& syn)
{
    const double syn_dt = 0.01; // assume 10 ms error
    
    out << syn.src.size() << '\n'; 
    for (int isrc=1; isrc<=syn.src.size(); isrc++){
	out << 's' << " "
	    << syn.src(isrc).x() << " "
	    << syn.src(isrc).y() << " "
	    << syn.rcv(isrc).size() << '\n';
	for (int ircv=1; ircv<=syn.rcv(isrc).size(); ircv++){
	    out << 'r' << " "
		<< syn.rcv(isrc)(ircv).x() << " "
		<< syn.rcv(isrc)(ircv).y() << " "
		<< syn.raytype(isrc)(ircv) << " "
		<< syn.syn_ttime(isrc)(ircv) << " "
		<< syn_dt << '\n';
	}
    }
}
