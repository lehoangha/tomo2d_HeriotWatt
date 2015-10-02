/*
 * tt_inverse.cc - traveltime inversion 
 *
 * usage: tt_inverse -Mmesh -Gdata [ -Nxorder/zorder/clen/nintp/bend_cg_tol/bend_br_tol ]
 *                   [ -Frefl -A -Llogfile -Oout [-olevel -l ] -Kdws_file ]
 *                   [ -P -Rcrit_chi -Qlsqr_tol -sbound -Wd_weight -Vlevel ]
 *
 *        [many iterations with a single set of parameters (type1)]
 *                   -Initer [ -Jtarget_chi2 -SVwsv -SDwsd -TVmax_dv -TDmax_dd ]
 *                           [ -DVwdv -DDwdd -DQdamp_v_fn ]
 *        [single iteration with many sets of parameters (type2)
 *         note: with the -X option, smoothing weights will be raised to
 *               the power of 10.]
 *                   -SVwsv_min/wsv_max/dw [ -SDwsd_min/wsd_max/dw
 *                   -TVmax_dv -TDmax_dd -XV -XD ]
 *        [ if smoothing is on, correlation length must be specified. ]
 *                   -CVcorr_v_fn -CDcorr_d_fn
 *
 *        [ joint gravity inversion ]
 *                   -ZGgrav.dat
 *                   -ZXxmin/xmax/zmin/zmax/dx/dz
 *                   -ZCcontup/iconv
 *                   -ZUoceanUup/oceanUlo/iconv
 *                   -ZLoceanLup/iconv
 *                   -ZSsedup/sedlo/iconv
 *                   -ZDdvdp/dvdt/drdp/drdt/dTdz
 *                   -ZRx0/x1
 *                   -ZWweight_g
 *                   -ZZz0
 *                   -ZKgrav_dws
 *                   -ZTcutoff_range/cutoff_val
 *
 * Jun Korenaga, MIT/WHOI
 * January 1999
 * Modified for joint gravity inversion: July 1999
 * Add squeezing option: Oct 1999
 */

#include <iostream>
#include <fstream>
#include <cstdio>
#include <array.h>
#include <util.h>
#include "inverse.h"
#include "jgrav.h"

int main(int argc, char** argv)
{
    bool getMesh=false, getData=false;
    char *meshfn, *datafn;
    int xorder=3, zorder=3, nintp=8;
    double crit_len=-1, bend_cg_tol=1e-4, bend_br_tol=1e-7;
    double crit_chi=-1, lsqr_atol=1e-3;
    double refl_weight=1;
    bool getRefl=false, outLog=false, outMask=false;
    bool doFullRefl=false;
    char *reflfn, *logfn, *maskfn, *outroot="inv_out";
    int outlevel=0;
    bool verbose=false, gotError=false;
    int verbose_level;

    int niter=1;
    bool smooth_vel=false, smooth_dep=false;
    bool apply2Dfilter=false;
    char *boundfn;
    double wsv_min=-1, wsv_max=-1, dwsv=-1;
    double wsd_min=-1, wsd_max=-1, dwsd=-1;
    double max_dv=-1, max_dd=-1;
    double wdv=0, wdd=0;
    double target_chisq=0;
    bool vlogscale=false, dlogscale=false;
    bool getCorrVel=false, getCorrDep=false;
    char *corr_velfn, *corr_depfn;
    bool printFinalOnly=false;
    bool auto_damping=false, fixed_damping=false, jumping=false;
    bool getDamp2d=false;
    char *damp_velfn;

    bool getGrav=false, getGcont=false, getGoceanU=false, getGoceanL=false;
    bool getGsed=false, getGderiv=false, getGgridspec=false, getGrefrange=false;
    bool print_gravDWS=false, getGcutoff=false;
    char *gravfn, contupfn[MaxStr], oceanUupfn[MaxStr], oceanUlofn[MaxStr];
    char oceanLupfn[MaxStr], sedupfn[MaxStr], sedlofn[MaxStr];
    char *grav_dwsfn;
    int icont, ioceanU, ioceanL, ised;
    double dvdp, dvdt, drdp, drdt, dTdz;
    double xmin, xmax, zmin, zmax, dx, dz, x0, x1;
    double weight_grav=1.0, z0=0.0, cutoff_range, cutoff_K;
    
    for (int i=1; i<argc; i++){
	if (argv[i][0] == '-'){
	    switch(argv[i][1]){
	    case 'M':
		meshfn = &argv[i][2];
		getMesh = true;
		break;
	    case 'G':
		datafn = &argv[i][2];
		getData = true;
		break;
	    case 'N':
	    {
		int tmp1a, tmp1b, tmp3;
		double tmp2, tmp4, tmp5;
		if (sscanf(&argv[i][2], "%d/%d/%lf/%d/%lf/%lf",
			   &tmp1a, &tmp1b, &tmp2, &tmp3, &tmp4, &tmp5)==6){
		    if (tmp1a>0 && tmp1b>0 && tmp2>0 && tmp3 > 0 && tmp4>0 && tmp5>0){
			xorder=tmp1a; zorder=tmp1b;
			crit_len=tmp2;
			nintp=tmp3;
			bend_cg_tol=tmp4;
			bend_br_tol=tmp5;
		    }else{
			error("invalid -N option, ignored.\n");
		    }
		}
		break;
	    }
	    case 'P':
		jumping = true; // pure jumping strategy
		break;
	    case 'R':
		crit_chi = atof(&argv[i][2]); // robust inversion
		break;
	    case 'Q':
		lsqr_atol = atof(&argv[i][2]);
		break;
	    case 'A':
		doFullRefl = true;
		break;
	    case 'F':
		getRefl=true;
		reflfn = &argv[i][2];
		break;
	    case 'W':
		refl_weight = atof(&argv[i][2]);
		break;
	    case 'L':
		logfn = &argv[i][2];
		outLog = true;
		break;
	    case 'O':
		outroot = &argv[i][2];
		break;
	    case 'o':
		outlevel = atoi(&argv[i][2]);
		break;
	    case 'l':
		printFinalOnly = true;
		break;
	    case 'K':
		maskfn = &argv[i][2];
		outMask = true;
		break;
	    case 'I':
		niter = atoi(&argv[i][2]);
		break;
	    case 'J':
		target_chisq = atof(&argv[i][2]);
		break;
	    case 'S':
		switch (argv[i][2]){
		case 'V':
		{
		    double a, b, c;
		    int nitem=sscanf(&argv[i][3], "%lf/%lf/%lf", &a, &b, &c);
		    if (nitem==1){ wsv_min=a; wsv_max=a; dwsv=a+1.0; }
		    else if (nitem==3){ wsv_min=a; wsv_max=b; dwsv=c; }
		    else{
			cerr << "invalid -SV option.\n";
			gotError=true;
		    }
		    smooth_vel = true;
		    break;
		}
		case 'D':
		{
		    double a, b, c;
		    int nitem=sscanf(&argv[i][3], "%lf/%lf/%lf", &a, &b, &c);
		    if (nitem==1){ wsd_min=a; wsd_max=a; dwsd=a+1.0; }
		    else if (nitem==3){ wsd_min=a; wsd_max=b; dwsd=c; }
		    else{
			cerr << "invalid -SD option.\n";
			gotError=true;
		    }
		    smooth_dep = true;
		    break;
		}
		default:
		    cerr << "invalid -S option.\n";
		    gotError=true;
		    break;
		}
		break;
	    case 's':
		apply2Dfilter = true;
		boundfn = &argv[i][2];
		break;
	    case 'X':
		switch (argv[i][2]){
		case 'V':
		    vlogscale=true;
		    break;
		case 'D':
		    dlogscale=true;
		    break;
		default:
		    cerr << "invalid -X option.\n";
		    gotError = true;
		    break;
		}
		break;
	    case 'T': // auto damping
		auto_damping = true;
		switch (argv[i][2]){
		case 'V':
		    max_dv = atof(&argv[i][3]);
		    break;
		case 'D':
		    max_dd = atof(&argv[i][3]);
		    break;
		default:
		    cerr << "invalid -T option.\n";
		    gotError = true;
		    break;
		}
		break; 
	    case 'D': // fixed damping
		fixed_damping = true;
		switch (argv[i][2]){
		case 'V':
		    wdv = atof(&argv[i][3]);
		    break;
		case 'D':
		    wdd = atof(&argv[i][3]);
		    break;
		case 'Q':
		    getDamp2d=true;
		    damp_velfn = &argv[i][3];
		    break;
		default:
		    cerr << "invalid -D option.\n";
		    gotError = true;
		    break;
		}
		break; 
	    case 'C':
		switch(argv[i][2]){
		case 'V':
		    getCorrVel=true;
		    corr_velfn = &argv[i][3];
		    break;
		case 'D':
		    getCorrDep=true;
		    corr_depfn = &argv[i][3];
		    break;
		default:
		    cerr << "invalid -C option.\n";
		    gotError = true;
		    break;
		}
		break;
	    case 'Z':
		switch(argv[i][2]){
		case 'G':
		    getGrav=true;
		    gravfn = &argv[i][3];
		    break;
		case 'C':
		    if (sscanf(&argv[i][3], "%[^/]/%d",
			       contupfn, &icont)==2){
			getGcont=true;
		    }else{
			cerr << "invalid -ZC option.\n";
			gotError=true;
		    }
		    break;
		case 'U':
		    if (sscanf(&argv[i][3], "%[^/]/%[^/]/%d",
			       oceanUupfn, oceanUlofn, &ioceanU)==3){
			getGoceanU=true;
		    }else{
			cerr << "invalid -ZU option.\n";
			gotError=true;
		    }
		    break;
		case 'L':
		    if (sscanf(&argv[i][3], "%[^/]/%d",
			       oceanLupfn, &ioceanL)==2){
			getGoceanL=true;
		    }else{
			cerr << "invalid -ZL option.\n";
			gotError=true;
		    }
		    break;
		case 'S':
		    if (sscanf(&argv[i][3], "%[^/]/%[^/]/%d",
			       sedupfn, sedlofn, &ised)==3){
			getGsed=true;
		    }else{
			cerr << "invalid -ZS option.\n";
			gotError=true;
		    }
		    break;
		case 'D':
		    if (sscanf(&argv[i][3], "%lf/%lf/%lf/%lf/%lf",
			       &dvdp, &dvdt, &drdp, &drdt, &dTdz)==5){
			getGderiv=true;
		    }else{
			cerr << "invalid -ZD option.\n";
			gotError=true;
		    }
		    break;
		case 'X':
		    if (sscanf(&argv[i][3], "%lf/%lf/%lf/%lf/%lf/%lf",
			       &xmin, &xmax, &zmin, &zmax, &dx, &dz)==6){
			getGgridspec=true;
		    }else{
			cerr << "invalid -ZX option.\n";
			gotError=true;
		    }
		    break;
		case 'R':
		    if (sscanf(&argv[i][3], "%lf/%lf", &x0, &x1)==2){
			getGrefrange=true;
		    }else{
			cerr << "invalid -ZR option.\n";
			gotError=true;
		    }
		    break;
		case 'W':
		    weight_grav = atof(&argv[i][3]);
		    break;
		case 'Z':
		    z0 = atof(&argv[i][3]);
		    break;
		case 'K':
		    print_gravDWS = true;
		    grav_dwsfn = &argv[i][3];
		    break;
		case 'T':
		    if (sscanf(&argv[i][3], "%lf/%lf",
			       &cutoff_range, &cutoff_K)==2){
			getGcutoff=true;
		    }else{
			cerr << "invalid -ZT option.\n";
			gotError=true;
		    }
		    break;
		default:
		    cerr << "unrecognized -Z sub option.\n";
		    gotError=true;
		    break;
		}
		break;
	    case 'V':
		verbose = true;
		verbose_level = 0;
		if (isdigit(argv[i][2])){
		    verbose_level = atoi(&argv[i][2]);
		}
		break;
	    default:
		gotError = true;
		break;
	    }
	}else{
	    gotError = true;
	}
    }

    if (!getMesh || !getData) gotError=true;
    if (smooth_vel && !getCorrVel){
	cerr << "velocity smoothing needs 2-D correlation length info.\n";
	gotError=true;
    }
    if (smooth_dep && (!getCorrDep && !getCorrVel)){
	cerr << "depth smoothing needs correlation length info.\n";
	gotError=true;
    }
    if (auto_damping && fixed_damping){
	cerr << "-T and -D options are mutually exclusive.\n";
	gotError=true;
    }
    if (getGrav && (!getGcont && !getGoceanU && !getGoceanL && !getGsed)){
	cerr << "at least one domain must be specified for -G option.\n";
	gotError=true;
    }
    if (getGrav && (!getGgridspec || !getGrefrange)){
	cerr << "-GZ & -GR is required for gravity inversion.\n";
	gotError=true;
    }
    if (getGrav && !getRefl){
	cerr << "reflector is needed for joint gravity inversion.\n";
	gotError=true;
    }
    if (gotError) error("usage: tt_inverse ...");

    SlownessMesh2d smesh(meshfn);
    TomographicInversion2d inv(smesh,datafn,xorder,zorder,crit_len,
			       nintp,bend_cg_tol,bend_br_tol);
    if (crit_chi>0) inv.doRobust(crit_chi);
    inv.setLSQR_TOL(lsqr_atol);
    if (printFinalOnly){
	inv.outFinal(outroot,outlevel);
    }else{
	inv.outStepwise(outroot,outlevel);
    }

    if (outLog) inv.setLogfile(logfn);
    if (outMask) inv.outMask(maskfn);
    if (verbose) inv.setVerbose(verbose_level);

    Interface2d *reflp;
    if (getRefl){
	reflp = new Interface2d(reflfn);
	inv.addRefl(reflp);
	inv.setReflWeight(refl_weight);
	if (doFullRefl) inv.doFullRefl();
    }

    if (smooth_vel){
	inv.SmoothVelocity(corr_velfn,wsv_min,wsv_max,dwsv,vlogscale);
	if (apply2Dfilter) inv.applyFilter(boundfn);
    }
    if (smooth_dep){
	if (getCorrDep){
	    inv.SmoothDepth(corr_depfn,wsd_min,wsd_max,dwsd,dlogscale);
	}else{
	    // use velocity's correlation length
	    inv.SmoothDepth(wsd_min,wsd_max,dwsd,dlogscale);
	}
    }
    if (target_chisq>0) inv.targetChisq(target_chisq);
    if (auto_damping){
	if (max_dv>0) inv.DampVelocity(max_dv);
	if (max_dd>0) inv.DampDepth(max_dd);
    }else if (fixed_damping){
	if (wdv>=0 && wdd>=0){
	    inv.FixDamping(wdv,wdd);
	    if (getDamp2d){
		inv.Squeezing(damp_velfn);
	    }
	}
    }
    if (jumping) inv.doJumping();

    Interface2d *contup, *oceanUup, *oceanUlo, *oceanLup;
    Interface2d *sedup, *sedlo;
    AddonGravityInversion2d *ginv;
    Array1d<double> xobs, obsgrav, grav2;
    if (getGrav){
	int ndata = countLines(gravfn);
	xobs.resize(ndata); obsgrav.resize(ndata); 
	ifstream is(gravfn);
	for (int i=1; i<=ndata; i++) is >> xobs(i) >> obsgrav(i);

	ginv = new AddonGravityInversion2d(smesh,*reflp,
					   RegularDomain2d(xmin, xmax, zmin, zmax),
					   dx, dz, x0, x1);
	if (getGcont){
	    contup = new Interface2d(contupfn);
	    ginv->defineContinentRegion(contup, reflp, icont);
	}
	if (getGoceanU){
	    oceanUup = new Interface2d(oceanUupfn);
	    oceanUlo = new Interface2d(oceanUlofn);
	    ginv->defineUpperOceanicRegion(oceanUup, oceanUlo, ioceanU);
	}
	if (getGoceanL){
	    oceanLup = new Interface2d(oceanLupfn);
	    ginv->defineLowerOceanicRegion(oceanLup, reflp, ioceanL);
	}
	if (getGsed){
	    sedup = new Interface2d(sedupfn);
	    sedlo = new Interface2d(sedlofn);
	    ginv->defineSedimentaryRegion(sedup, sedlo, ised);
	}
	if (getGderiv) ginv->setDerivative(dvdp, dvdt, drdp, drdt, dTdz);
	if (getGcutoff) ginv->setThreshold(cutoff_range, cutoff_K);
	    
//	grav2.resize(ndata);
//	ginv->calcGravity(z0, xobs, grav2);
//	ofstream dos("dump.grav");
//	for (int i=1; i<=ndata; i++){
//	    cerr << "tt_inv: " << xobs(i) << " " << grav2(i) << '\n';
//	}
	if (print_gravDWS){ 
	    inv.addonGravity(z0,xobs,obsgrav,ginv,weight_grav,grav_dwsfn);
	}else{
	    inv.addonGravity(z0,xobs,obsgrav,ginv,weight_grav);
	}
    }
    inv.solve(niter);

    return 0;
}
