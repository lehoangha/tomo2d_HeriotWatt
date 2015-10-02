/*
 * traveltime.cc - travel time integration along a ray path
 *
 * Jun Korenaga, MIT/WHOI
 * January 1999
 */

#include "traveltime.h"
#include "interface.h"
#include <iostream>
#include <fstream>

double calcTravelTime(const SlownessMesh2d& smesh,
		      const Array1d<Point2d>& path,
		      const BetaSpline2d& bs,
		      const Array1d<const Point2d*>& pp,
		      Array1d<Point2d>& Q)
{
    int np = path.size();
    int nintp = bs.numIntp();

    Index2d guess_index = smesh.nodeIndex(smesh.nearest(*pp(1)));
    double ttime=0.0;
    for (int i=1; i<=np+1; i++){
	int j1=i;
	int j2=i+1;
	int j3=i+2;
	int j4=i+3;
	bs.interpolate(*pp(j1),*pp(j2),*pp(j3),*pp(j4),Q);

	if (smesh.inWater(Q(nintp/2))){
	    // rectangular integration
	    for (int j=2; j<=nintp; j++){
		Point2d midp = 0.5*(Q(j-1)+Q(j));
		double u0 = smesh.at(midp,guess_index);
		double dist = Q(j).distance(Q(j-1));
		ttime += u0*dist;
	    }		
	}else{
	    // trapezoidal integration along one segment
	    double u0 = smesh.at(Q(1),guess_index);
	    for (int j=2; j<=nintp; j++){
		double u1 = smesh.at(Q(j),guess_index);
		double dist= Q(j).distance(Q(j-1));
		ttime += 0.5*(u0+u1)*dist;
		u0 = u1;
	    }
	}
    }
    return ttime;
}

double calcTravelTime(const SlownessMesh2d& smesh,
		      const Array1d<Point2d>& path,
		      int ndiv)
{
    int np = path.size();
    Index2d guess_index = smesh.nodeIndex(smesh.nearest(path(1)));
    double dd = 1.0/ndiv;

    double ttime=0.0;
    for (int i=2; i<=np; i++){
	double dist= path(i).distance(path(i-1));

	Point2d midp = 0.5*(path(i-1)+path(i));
	if (smesh.inWater(midp)){
	    // rectangular integration
	    ttime += smesh.atWater()*dist;
	}else{
	    // trapezoidal integration
	    double ddist = dist*dd;
	    double u0 = smesh.at(path(i-1),guess_index);
	    for (int j=1; j<=ndiv; j++){
		double ratio = j*dd;
		Point2d tmpp = (1-ratio)*path(i-1)+ratio*path(i);
		double u1 = smesh.at(tmpp,guess_index);
		ttime += 0.5*(u0+u1)*ddist;
		u0 = u1;
	    }
	}
    }
    return ttime;
}

void calc_dTdV2(const SlownessMesh2d& smesh, const Array1d<Point2d>& path,
 		const BetaSpline2d& bs, Array1d<Point2d>& dTdV,
 		const Array1d<const Point2d*>& pp,
 		Array1d<Point2d>& Q, Array1d<Point2d>& dQdu)
 {
     int np = path.size();
     Array1d<Point2d> path2(np);
     for (int i=1; i<=np; i++) path2(i) = path(i);
     Array1d<const Point2d*> pp2;
     makeBSpoints(path2,pp2);

     double dx=0.0001;
     dTdV(1).set(0.0,0.0);  // end points are fixed
     dTdV(np).set(0.0,0.0);

     for (int i=2; i<np; i++){
 	double orig=path2(i).x();
 	path2(i).x(orig+dx);
 	double t1 = calcTravelTime(smesh,path2,bs,pp2,Q);
 	path2(i).x(orig-dx);
 	double t2 = calcTravelTime(smesh,path2,bs,pp2,Q);
 	dTdV(i).x((t1-t2)/(2*dx));
 	path2(i).x(orig);

 	orig=path2(i).y();
 	path2(i).y(orig+dx);
 	t1 = calcTravelTime(smesh,path2,bs,pp2,Q);
 	path2(i).y(orig-dx);
 	t2 = calcTravelTime(smesh,path2,bs,pp2,Q);
 	dTdV(i).y((t1-t2)/(2*dx));
 	path2(i).y(orig);
     }	
 }

void calc_dTdV3(const SlownessMesh2d& smesh, const Array1d<Point2d>& path,
 		const BetaSpline2d& bs, Array1d<Point2d>& dTdV,
 		const Array1d<const Point2d*>& pp,
 		Array1d<Point2d>& Q, Array1d<Point2d>& dQdu,
		const Array1d<int>& start_i,
		const Array1d<int>& end_i,
		const Array1d<const Interface2d*>& interf)
{
    int np = path.size();
    Array1d<Point2d> path2(np);
    for (int i=1; i<=np; i++) path2(i) = path(i);
    Array1d<const Point2d*> pp2;
    makeBSpoints(path2,pp2);

    double dx=0.0001;
    dTdV(1).set(0.0,0.0);  // end points are fixed
    dTdV(np).set(0.0,0.0);

    int a=1;
    for (int i=2; i<np; i++){
	if (i==start_i(a)){
	    double origx=path2(i).x();
	    double origy=path2(i).y();
	    for (int j=i; j<=end_i(a); j++){
		path2(j).x(origx+dx);
		path2(j).y(interf(a)->z(origx+dx));
	    }
	    double t1 = calcTravelTime(smesh,path2,bs,pp2,Q);
	    for (int j=i; j<=end_i(a); j++){
		path2(j).x(origx-dx);
		path2(j).y(interf(a)->z(origx-dx));
	    }
	    double t2 = calcTravelTime(smesh,path2,bs,pp2,Q);
	    double dtdv = (t1-t2)/(2*dx);
	    dtdv /= 3.0;
 	    cerr << "\n" << "dtdv3: t1=" << t1 << " t2=" << t2
 		 << " dtdv=" << dtdv << '\n';

	    for (int j=i; j<=end_i(a); j++){
		dTdV(j).x(dtdv);
		dTdV(j).y(0.0);
	    }
	    for (int j=i; j<=end_i(a); j++){
		path2(j).x(origx);
		path2(j).y(origy);
	    }
	    i = end_i(a);
	    a++;
	}else{
	    // move x
	    double orig=path2(i).x();
	    path2(i).x(orig+dx);
	    double t1 = calcTravelTime(smesh,path2,bs,pp2,Q);
	    path2(i).x(orig-dx);
	    double t2 = calcTravelTime(smesh,path2,bs,pp2,Q);
	    dTdV(i).x((t1-t2)/(2*dx));
	    path2(i).x(orig);

	    // move y
	    orig=path2(i).y();
	    path2(i).y(orig+dx);
	    t1 = calcTravelTime(smesh,path2,bs,pp2,Q);
	    path2(i).y(orig-dx);
	    t2 = calcTravelTime(smesh,path2,bs,pp2,Q);
	    dTdV(i).y((t1-t2)/(2*dx));
	    path2(i).y(orig);
	}
    }	
}

void calc_dTdV(const SlownessMesh2d& smesh, const Array1d<Point2d>& path,
	       const BetaSpline2d& bs, Array1d<Point2d>& dTdV,
	       const Array1d<const Point2d*>& pp,
	       Array1d<Point2d>& Q, Array1d<Point2d>& dQdu)
{
    int np = path.size();
    int nintp = bs.numIntp();
    double du = 1.0/(nintp-1);
    
    dTdV(1).set(0.0,0.0); // end points are fixed
    dTdV(np).set(0.0,0.0);

    Index2d guess_index = smesh.nodeIndex(smesh.nearest(*pp(1)));
    for (int i=2; i<np; i++){
	double dTdVx=0.0, dTdVz=0.0;

	// loop over four segments affected by the point Vi
	int iseg_start = i-1;
	int iseg_end = i+2; // min(i+2,np-1);
	for (int j=iseg_start; j<=iseg_end; j++){
	    int j1=j;
	    int j2=j+1;
	    int j3=j+2;
	    int j4=j+3;
	    bs.interpolate(*pp(j1),*pp(j2),*pp(j3),*pp(j4),Q,dQdu);
	    int i_j = i-j;

	    // trapezoidal integration along one segment
	    double dudx0, dudz0, dudx1, dudz1;
	    Point2d startp = 0.9*Q(1)+0.1*Q(2); // in order to avoid singularity
	    double u0 = smesh.at(startp,guess_index,dudx0,dudz0);
	    double dQnorm0 = dQdu(1).norm();
	    double integx0, integz0;
	    if (dQnorm0==0.0){
		integx0 = integz0 = 0.0;
	    }else{
		double tmp1 = u0*bs.coeff_dbdu(i_j,1)/dQnorm0;
		double tmp2 = dQnorm0*bs.coeff_b(i_j,1);
		integx0 = tmp1*dQdu(1).x()+tmp2*dudx0;
		integz0 = tmp1*dQdu(1).y()+tmp2*dudz0;
	    }		

	    for (int m=2; m<=nintp; m++){
		double u1;
		if (m==nintp){
		    Point2d endp = 0.1*Q(m-1)+0.9*Q(m); // in order to avoid singularity
		    u1 = smesh.at(endp,guess_index,dudx1,dudz1);
		}else{
		    u1 = smesh.at(Q(m),guess_index,dudx1,dudz1);
		}
		double dQnorm1 = dQdu(m).norm();
		double integx1, integz1;
		if (dQnorm1==0.0){
		    integx1 = integz1 = 0.0;
		}else{
		    double tmp1 = u1*bs.coeff_dbdu(i_j,m)/dQnorm1;
		    double tmp2 = dQnorm1*bs.coeff_b(i_j,m);
		    integx1 = tmp1*dQdu(m).x()+tmp2*dudx1;
		    integz1 = tmp1*dQdu(m).y()+tmp2*dudz1;
		}

		dTdVx += 0.5*(integx0+integx1)*du;
		dTdVz += 0.5*(integz0+integz1)*du;
		
		integx0 = integx1;
		integz0 = integz1;
	    }
	}
	dTdV(i).set(dTdVx,dTdVz);
    }
}
