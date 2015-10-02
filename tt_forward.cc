/*
 * tt_forward.cc - forward traveltime calculation
 *
 * usage: tt_forward -Msmesh [ -Ggeom -Frefl -A ] \
 *                 [ -Nxorder/zorder/clen/nintp/tot1/tot2 -Eelem -g \
 *                   -Tttime -Oobs_ttime -rv0 -Ddiff -Rray -Ssrc -Ivel -iw/e/s/n/dx/dz -n -Cused_time -Vlevel ]
 *
 * if -G is not specified, only operations regarding a slowness mesh will be done.
 *
 * Jun Korenaga, MIT/WHOI
 * January 1999
 */

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cctype>
#include <array.h>
#include "syngen.h"

int main(int argc, char** argv)
{
    bool getMesh=false, getGeom=false, refl=false;
    bool outTime=false, outOTime=false, outSource=false;
    bool outVgrid=false, outRay=false, outElements=false;
    bool outDiff=false, useClock=false, verbose=false, err=false;
    bool graph_only=false, getWESN=false, printAW=true;
    bool doFullRefl=false;
    int verbose_level;
    char *mfn, *gfn, *efn, *ofn, *rfn, *tfn, *sfn, *dfn, *vfn, *reflfn, *cfn;
    double vred=0.0;
    int xorder=4, zorder=4, nintp=8;
    double clen=0.0, bend_cg_tol=1e-4, bend_br_tol=1e-7;
    double west,east,north,south,dx,dz;

    for (int i=1; i<argc; i++){
	if (argv[i][0] == '-'){
	    switch(argv[i][1]){
	    case 'M':
		mfn = &argv[i][2];
		getMesh = true;
		break;
	    case 'G':
		gfn = &argv[i][2];
		getGeom = true;
		break;
	    case 'g':
		graph_only = true;
		break;
	    case 'F':
		reflfn = &argv[i][2];
		refl = true;
		break;
	    case 'A':
		doFullRefl = true;
		break;
	    case 'N':
	    {
		int tmp1a, tmp1b, tmp3;
		double tmp2, tmp4, tmp5;
		if (sscanf(&argv[i][2], "%d/%d/%lf/%d/%lf/%lf",
			   &tmp1a, &tmp1b, &tmp2, &tmp3, &tmp4, &tmp5)==6){
		    if (tmp1a>0 && tmp1b>0 && tmp2>0 && tmp3>0 && tmp4>0 && tmp5>0){
			xorder=tmp1a; zorder=tmp1b;
			clen=tmp2;
			nintp=tmp3;
			bend_cg_tol=tmp4;
			bend_br_tol=tmp5;
		    }else{
			error("invalid -N option, ignored.\n");
		    }
		}
		break;
	    }
	    case 'E':
		efn = &argv[i][2];
		outElements = true;
		break;
	    case 'T':
		tfn = &argv[i][2];
		outTime = true;
		break;
	    case 'O':
		ofn = &argv[i][2];
		outOTime = true;
		break;
	    case 'D':
		dfn =  &argv[i][2];
		outDiff = true;
		break;
	    case 'r':
		vred = atof(&argv[i][2]);
		break;
	    case 'R':
		rfn = &argv[i][2];
		outRay = true;
		break;
	    case 'S':
		sfn = &argv[i][2];
		outSource = true;
		break;
	    case 'C':
		cfn = &argv[i][2];
		useClock = true;
		break;
	    case 'I':
		vfn = &argv[i][2];
		outVgrid = true;
		break;
	    case 'i':
		if (sscanf(&argv[i][2],"%lf/%lf/%lf/%lf/%lf/%lf",
			   &west,&east,&south,&north,&dx,&dz)==6){
		    getWESN=true;
		}else{
		    cerr << "invalid -i option\n";
		}
		break;
	    case 'n':
		printAW = false; // do not print out extra grids for air and water
		break;
	    case 'V':
		verbose = true;
		verbose_level = 0;
		if (isdigit(argv[i][2])){
		    verbose_level = atoi(&argv[i][2]);
		}
		break;
	    default:
		err = true;
		break;
	    }
	}else{
	    err = true;
	}
    }

    if (!getMesh) err=true;
    if (err) error("usage: tt_forward ...");

    SlownessMesh2d smesh(mfn);
    if (outElements){
	ofstream os(efn);
	smesh.printElements(os);
    }
    if (outVgrid){
	ofstream os(vfn);
	if (getWESN){
	    smesh.printVGrid(os,west,east,south,north,dx,dz);
	}else{
	    smesh.printVGrid(os,printAW);
	}
    }
    if (!getGeom) return 0; 

    SyntheticTraveltimeGenerator2d syngen(smesh,gfn,
					  xorder,zorder,clen,
					  nintp,bend_cg_tol,bend_br_tol);
    if (refl){
	syngen.readRefl(reflfn);
	if (doFullRefl) syngen.doFullRefl();
    }
    if (outRay) syngen.outputRay(rfn);
    if (useClock) syngen.useClock(cfn);
    if (graph_only) syngen.graphOnly();
    if (verbose) syngen.setVerbose(verbose_level);

    syngen.conduct();
    cout << syngen;

    if (outSource){
	ofstream os(sfn);
	syngen.printSource(os);
    }
    if (outTime){
	ofstream os(tfn);
	syngen.printSynTime(os,vred);
    }
    if (outOTime){
	ofstream os(ofn);
	syngen.printObsTime(os,vred);
    }
    if (outDiff){
	ofstream os(dfn);
	double ave_misfit, chisq;
	syngen.printDiff(os,ave_misfit,chisq);
	if (verbose){
	    cerr << "average traveltime misfit: " << ave_misfit << '\n';
	    cerr << "chisq: " << chisq << '\n';
	}
    }

    return 0;
}

