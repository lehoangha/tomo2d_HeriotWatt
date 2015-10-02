/*
 * stat_smesh.cc
 *
 * usage: stat_smesh  -Llist_file -C<cmd> [ -Rn ]
 *                 or -Mmesh -D<cmd> [ -Ttopb -Bbotb -m<midb> -P<TPcorr> -U<vrepl> -Xxmin/xmax ]
 *                                   [ -x<cxmin>/<cxmax> -t<ctopb> -b<cbotb> ]
 */

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "array.h"
#include "interface.h"

int main(int argc, char **argv)
{
    bool getList=false, getMesh=false, getListCmd=false, getMeshCmd=false;
    bool average=false, take_rms=false, getTop=false, getBot=false, getMid=false;
    bool haverage=false, hvaverage=false, replV=false;
    bool is_refl=false, applyPTcorr=false, removeCont=false, getCtop=false, getCbot=false;
    bool verbose=false, err=false;
    char *lfn, *ave_fn, *mfn, *tbfn, *bbfn, *mbfn, *ctbfn, *cbbfn;
    double avex, avexmin, avexmax, avedx, wlen;
    double vrepl, Tref, Pref, dVdT, dVdP, aT, bT;
    double abs_xmin=-1e30, abs_xmax=1e30;
    double cxmin, cxmax;
    int np;
    
    for (int i=1; i<argc; i++){
	if (argv[i][0] == '-'){
	    switch(argv[i][1]){
	    case 'M':
		getMesh = true;
		mfn = &argv[i][2];
		break;
	    case 'D':
		getMeshCmd=true;
		getListCmd=false;
		switch(argv[i][2]){
		case 'a': // horizontal average
		    if (sscanf(&argv[i][3],"%lf/%lf",&avex,&wlen)==2){
			haverage=true;
		    }else{
			cerr << "invalid -Da option\n";
			err=true;
		    }
		    break;
		case 'b': // horizontal and vertical average
		    if (sscanf(&argv[i][3],"%lf/%lf/%lf/%lf",&avexmin,&avexmax,&avedx,&wlen)==4){
			hvaverage=true;
		    }else{
			cerr << "invalid -Db option\n";
			err=true;
		    }
		    break;
		default:
		    err=true;
		    break;
		}
		break;
	    case 'P':
		if (sscanf(&argv[i][2], "%lf/%lf/%lf/%lf/%lf/%lf",
			   &Tref, &Pref, &dVdT, &dVdP, &aT, &bT)==6){
		    applyPTcorr=true;
		}else{
		    cerr << "invalid -P option\n";
		    err = true;
		}
		break;
	    case 'U':
		replV = true;
		vrepl = atof(&argv[i][2]);
		break;
	    case 'X':
		if (sscanf(&argv[i][2],"%lf/%lf",&abs_xmin,&abs_xmax)!=2){
		    cerr << "invalid -X option\n";
		    err=true;
		}
		break;
	    case 'T':
		getTop = true;
		tbfn = &argv[i][2];
		break;
	    case 'B':
		getBot = true;
		bbfn = &argv[i][2];
		break;
	    case 'm':
		getMid = true;
		mbfn = &argv[i][2];
		break;
	    case 'x':
		if (sscanf(&argv[i][2],"%lf/%lf",&cxmin,&cxmax)==2){
		    removeCont=true;
		}else{
		    cerr << "invalid -x option\n";
		    err=true;
		}
		break;
	    case 't':
		getCtop = true;
		ctbfn = &argv[i][2];
		break;
	    case 'b':
		getCbot = true;
		cbbfn = &argv[i][2];
		break;
	    case 'L':
		getList = true;
		lfn = &argv[i][2];
		break;
	    case 'C':
	    {
		getListCmd=true;
		getMeshCmd=false;
		switch(argv[i][2]){
		case 'a':
		    average=true;
		    break;
		case 'r':
		    take_rms=true;
		    ave_fn = &argv[i][3];
		    break;
		default:
		    err = true;
		    break;
		}   
		break;
	    }
	    case 'R':
		is_refl = true;
		np = atoi(&argv[i][2]);
		break;
	    case 'V':
		verbose=true;
		break;
	    default:
		err = true;
		break;
	    }
	}else{
	    err = true;
	}
    }

    if (!getList && !getMesh){ err = true; cerr << "1 ";}
    if (getList && !getListCmd){ err=true; cerr << "2 ";}
    if (getMesh && !getMeshCmd){ err=true; cerr << "3 ";}
    if (hvaverage && (!getTop || !getMid || !getBot)){ err = true; cerr << "4 ";}
    if (removeCont && (!getCtop || !getCbot)){ err=true; cerr << "5 ";}
    if (err) error("usage: stat_smesh [ -options ]");

    if (getMesh){
	// read smesh
	ifstream s(mfn);
	if (!s){
	    cerr << "stat_smesh::cannot open " << mfn << "\n";
	    exit(1);
	}
	int nx, nz;
	double v_water, v_air;
	s >> nx >> nz >> v_water >> v_air;
	int nnodes = nx*nz;
    
	Array1d<double> xpos(nx), zpos(nz), topo(nx);
	Array2d<double> vgrid(nx,nz);
	for (int i=1; i<=nx; i++) s >> xpos(i);
	for (int i=1; i<=nx; i++) s >> topo(i);
	for (int k=1; k<=nz; k++) s >> zpos(k);
	for (int i=1; i<=nx; i++){
	    for (int k=1; k<=nz; k++){
		double vel;
		s >> vel;
		vgrid(i,k) = vel;
	    }
	}

	Array1d<int> kstart(nx);
	kstart = 1;
	if (getTop){
	    Interface2d top(tbfn);
	    for (int i=1; i<=nx; i++){
		double boundz = top.z(xpos(i));
		for (int k=1; k<=nz; k++){
		    if (zpos(k)+topo(i)>boundz) break;
		    kstart(i) = k;
		}
		if (kstart(i)>1) kstart(i)++;
	    }
	}

	Array1d<int> kend(nx);
	kend = nz;
	if (getBot){
	    Interface2d bot(bbfn);
	    for (int i=1; i<=nx; i++){
		double boundz = bot.z(xpos(i));
		for (int k=nz; k>=1; k--){
		    kend(i) = k;
		    if (zpos(k)+topo(i)<boundz) break;
		}
	    }
	}

	Array1d<int> kmid(nx);
	kmid = 1;
	if (getMid){
	    Interface2d mid(mbfn);
	    for (int i=1; i<=nx; i++){
		double boundz = mid.z(xpos(i));
		for (int k=1; k<=nz; k++){
		    if (zpos(k)+topo(i)>boundz) break;
		    kmid(i) = k;
		}
		if (kmid(i)>1) kmid(i)++;
	    }
	}
	    
	if (applyPTcorr){
	    for (int i=1; i<=nx; i++){
		double t=topo(i);
		double P=0;
		double g = 9.8;
		if (t>0){
		    P += 1e3*g*t*1e3/1e6; // water column pressure
		}
		for (int k=1; k<=nz; k++){
		    double z = zpos(k);
		    double vel = vgrid(i,k);		    
		    double T = aT*z+bT;  // depth beneath seafloor -> temperature
		    double dVt = dVdT*(T-Tref);
		    double dVp = dVdP*(P-Pref);

//		    if (xpos(i) > 320 && xpos(i) < 325){
//		    cerr << xpos(i) << " " << z+topo(i) << " " << T << " " << P << " "
//			 << vgrid(i,k) << " " << dVt << " " << dVp << " ";
//		    }

		    vgrid(i,k) -= (dVt+dVp);

//		    if (xpos(i) > 320 && xpos(i) < 325){
//		    cerr << vgrid(i,k) << '\n';
//		    }
//    if (i==1){
//	cerr << z << " " << -(dVt+dVp) << '\n';
//    }
			
		    if (k<nz){
			double dz = zpos(k+1)-zpos(k);
			double rho;
			if (vel < 1.6){
			    rho = 1e3;
			}else if (vel >=1.6 && vel < 6.1){
			    rho = (1.32+0.225*vel)*1e3;
			}else if (vel>=6.1){
			    rho = (0.81+0.31*vel)*1e3;
			}
			P += rho*g*dz*1e3/1e6;
		    }

		}
	    }
	}

	if (replV){
	    for (int i=1; i<=nx; i++){
		for (int k=1; k<=nz; k++){
		    if (vgrid(i,k)<vrepl) vgrid(i,k) = vrepl;
		}
	    }
	}

	if (removeCont){
	    Interface2d top(ctbfn); 
	    Interface2d bot(cbbfn); 
	    for (int i=1; i<=nx; i++){
		if (cxmin<=xpos(i) && xpos(i)<=cxmax){
		    double tbz = top.z(xpos(i));
		    double bbz = bot.z(xpos(i));
		    int ckstart=nz, ckend=1;
		    for (int k=1; k<=nz; k++){
			if (zpos(k)+topo(i)>tbz) break;
			ckstart = k;
		    }
		    if (ckstart>1) ckstart++;
		    for (int k=nz; k>=1; k--){
			ckend = k;
			if (zpos(k)+topo(i)<bbz) break;
		    }
		    for (int k=ckstart; k<=ckend; k++){
			vgrid(i,k) = -1;
		    }
		}
	    }
	}

	if (haverage){
	    Array1d<double> avevel(nz), avez(nz);
	    Array1d<int> nentry(nz);
	    double xmin=avex-wlen, xmax=avex+wlen;
	    avevel = avez = 0.0; nentry = 0;
	    for (int i=1; i<=nx; i++){
		if (abs_xmin<=xpos(i) && xpos(i)<=abs_xmax){
		    if (xmin<=xpos(i) && xpos(i)<=xmax){
			int j=1;
			for (int k=kstart(i); k<=nz; k++){
			    double z = topo(i)+zpos(k);
			    avevel(j) += vgrid(i,k);
			    avez(j) += z;
			    nentry(j)++;
			    j++;
			}
		    }
		}
	    }
	    for (int k=1; k<=nz; k++){
		avevel(k) /= nentry(k);
		avez(k) /= nentry(k);
		cout << avez(k) << " " << avevel(k) << '\n';
	    }
	}

	if (hvaverage){
	    for (double ax=avexmin; ax<=avexmax; ax+=avedx){
		double avevel=0, beta=0, avehw=0, aveh=0, alpha=0;
//		double avevel2=0;
		double xmin=ax-wlen, xmax=ax+wlen;
		for (int i=1; i<=nx; i++){
		    if (abs_xmin<=xpos(i) && xpos(i)<=abs_xmax){
			if (xmin<=xpos(i) && xpos(i)<=xmax){
			    double dx = i<nx ? (xpos(i+1)-xpos(i)) : (xpos(i)-xpos(i-1));
			    alpha += dx;
			    for (int k=kstart(i); k<=kend(i); k++){ // whole crust
				if (vgrid(i,k)>0){
				    double dz = k<nz ? (zpos(k+1)-zpos(k)) : (zpos(k)-zpos(k-1));
				    avehw += dx*dz;
				}
			    }
			    for (int k=kmid(i); k<=kend(i); k++){ // lower crust
				if (vgrid(i,k)>0){
				    double dz = k<nz ? (zpos(k+1)-zpos(k)) : (zpos(k)-zpos(k-1));
				    avevel += dx*dz/vgrid(i,k); // harmonic mean
				    aveh += dx*dz;
				    beta += dx*dz;

//				    avevel2 += vgrid(i,k)*dx*dz; // arithmetic mean
				    
//				    if (xpos(i) >320 && xpos(i) < 325){
//					cerr << "hva " << xpos(i) << " " << zpos(k)+topo(i) << " "
//					     << vgrid(i,k) << '\n';
//				    }
				}
			    }
			}
		    }
		}
		cout << ax << " " << beta/avevel << " " << aveh/alpha << " " << avehw/alpha
		     << '\n';
//		     << " " << avevel2/beta << '\n';
	    }
	}
	
	return 0;
    }

    if (getList){
	Array2d<double> avegrid;
	Array1d<double> averefl;
	if (take_rms){
	    ifstream is(ave_fn);
	    if (!is){
		cerr << "stat_smesh::cannot open " << ave_fn << "\n";
		exit(1);
	    }
	    if (is_refl){
		averefl.resize(np);
		double dum;
		for (int i=1; i<=np; i++)
		    is >> dum >> averefl(i);
	    }else{
		int nx, nz;
		double v_water, v_air;
		is >> nx >> nz >> v_water >> v_air;
		Array1d<double> xpos(nx), topo(nx), zpos(nz);
		avegrid.resize(nx,nz);
    
		for (int i=1; i<=nx; i++) is >> xpos(i);
		for (int i=1; i<=nx; i++) is >> topo(i);
		for (int k=1; k<=nz; k++) is >> zpos(k);
		for (int i=1; i<=nx; i++){
		    for (int k=1; k<=nz; k++){
			is >> avegrid(i,k);
		    }
		}
	    }
	}
	
	// read list
	ifstream s(lfn);
	if (!s){
	    cerr << "stat_smesh::cannot open " << lfn << "\n";
	    exit(1);
	}
	string sfn;
	int ifile = 0;
	int nx, nz;
	double v_water, v_air;
	Array1d<double> xpos, zpos, topo;
	Array2d<double> vgrid;
	Array1d<double> refl, range;
	while(s >> sfn){
	    if (verbose) cerr << "reading " << sfn << "...\n";
	    ifile++;
	    ifstream is(sfn.c_str());
	    if (is_refl){
		if (ifile==1){
		    range.resize(np);
		    refl.resize(np);
		}
		double rms_per_file=0;
		for (int i=1; i<=np; i++){
		    double val;
		    is >> range(i) >> val;
		    if (average){
			refl(i) += val;
		    }else if (take_rms){
			double dv = val-averefl(i);
			double dv2 = dv*dv;
			rms_per_file += dv2;
			refl(i) += dv2;
		    }
		}
		if (take_rms && verbose){
		    cerr << "rms=" << sqrt(rms_per_file/np) << '\n';
		}
	    }else{
		is >> nx >> nz >> v_water >> v_air;
		if (ifile==1){
		    xpos.resize(nx); zpos.resize(nz); topo.resize(nx);
		    vgrid.resize(nx,nz);
		    vgrid=0.0;
		}
		for (int i=1; i<=nx; i++) is >> xpos(i);
		for (int i=1; i<=nx; i++) is >> topo(i);
		for (int k=1; k<=nz; k++) is >> zpos(k);
		for (int i=1; i<=nx; i++){
		    for (int k=1; k<=nz; k++){
			double vel;
			is >> vel;
			if (average){
			    vgrid(i,k) += vel;
			}else if (take_rms){
			    double dv = vel-avegrid(i,k);
			    vgrid(i,k) += dv*dv;
			}
		    }
		}
	    }
	}
	if (verbose){
	    cerr << "number of ensembles: " << ifile << '\n';
	}

	if (is_refl){
	    if (average){
		refl /= ifile;
	    }else if (take_rms){
		for (int i=1; i<=np; i++)
		    refl(i) = sqrt(refl(i)/ifile);
	    }
	}else{
	    if (average){
		vgrid /= ifile;
	    }else if (take_rms){
		for (int i=1; i<=nx; i++){
		    for (int k=1; k<=nz; k++){
			vgrid(i,k) = sqrt(vgrid(i,k)/ifile);
		    }
		}
//		v_water = v_air = 0.0;
	    }
	}

	// output smesh
	if (is_refl){
	    for (int i=1; i<=np; i++)
		cout << range(i) << " " << refl(i) << '\n';
	}else{
	    cout << nx << " " << nz << " "
		 << v_water << " " << v_air << '\n';
	    for (int i=1; i<=nx; i++) cout << xpos(i) << " ";
	    cout << '\n';
	    for (int i=1; i<=nx; i++) cout << topo(i) << " ";
	    cout << '\n';
	    for (int k=1; k<=nz; k++) cout << zpos(k) << " ";
	    cout << '\n';
	    for (int i=1; i<=nx; i++){
		for (int k=1; k<=nz; k++){
		    cout << vgrid(i,k) << " ";
		}
		cout << '\n';
	    }
	}
	return 0;
    }
}
