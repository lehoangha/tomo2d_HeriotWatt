/*
 * bend.cc - ray-bending solver implementation
 *           direct minimization of travel times by conjugate gradient
 *           (rays are represented as beta-splines)
 *
 * Jun Korenaga, MIT/WHOI
 * January 1999
 */

#include "bend.h"
#include "traveltime.h"

BendingSolver2d::BendingSolver2d(const SlownessMesh2d& s,
				 const BetaSpline2d& b,
				 double tol1, double tol2)
    : smesh(s), bs(b), nintp(bs.numIntp()),
      cg_tol(tol1), brent_tol(tol2), eps(1e-10)
{
    Q.resize(nintp);
    dQdu.resize(nintp);

    int max_np = int(2*sqrt(float(smesh.numNodes())));
    new_point.reserve(max_np);
    dTdV.reserve(max_np);
    new_dTdV.reserve(max_np);
    direc.reserve(max_np);

    pp.reserve(max_np+4);
    new_pp.reserve(max_np+4);
}

int BendingSolver2d::refine(Array1d<Point2d>& path,
			    double& orig_time, double& new_time)
{
    int np=check_path_size(path);
    makeBSpoints(path,pp);
    new_point.resize(np); 
    makeBSpoints(new_point,new_pp);

    dTdV.resize(np);
    new_dTdV.resize(np);
    direc.resize(np);
    
    orig_time=calcTravelTime(smesh, path, bs, pp, Q);
    double old_time = orig_time;

    calc_dTdV(smesh, path, bs, dTdV, pp, Q, dQdu);
    direc = -dTdV; // p0
    // conjugate gradient search
    int iter_max = np*10;
    for (int iter=1; iter<=iter_max; iter++){
	new_time=line_min(path, direc);
	if ((old_time-new_time)<=cg_tol){
	    return iter;
	}
	
	calc_dTdV(smesh, path, bs, new_dTdV, pp, Q, dQdu);
	double gg=0.0, dgg=0.0;
	for (int i=1; i<=np; i++){
	    double x0=dTdV(i).x(),y0=dTdV(i).y();
	    double x1=new_dTdV(i).x(),y1=new_dTdV(i).y();
	    gg += x0*x0+y0*y0;
	    dgg += x1*x1+y1*y1;
	}
//	if (gg==0.0){
	if (abs(gg)<1e-10){
	    return iter;
	}else{
	    dgg /= gg;
	}
	for (int i=1; i<=np; i++){
	    direc(i) = -new_dTdV(i)+dgg*direc(i);
	}
	dTdV = new_dTdV;
	old_time = new_time;
    }
    return -iter_max; // too many iterations
}

int BendingSolver2d::refine(Array1d<Point2d>& path,
			    double& orig_time, double& new_time,
			    const Array1d<int>& start_i,
			    const Array1d<int>& end_i,
			    const Array1d<const Interface2d*>& interf)
{
    int np=check_path_size(path,start_i,end_i,interf);

    makeBSpoints(path,pp);
    new_point.resize(np); 
    makeBSpoints(new_point,new_pp);
    
    dTdV.resize(np);
    new_dTdV.resize(np);
    direc.resize(np);

    orig_time=calcTravelTime(smesh, path, bs, pp, Q);
    double old_time = orig_time;

    calc_dTdV(smesh, path, bs, dTdV, pp, Q, dQdu);
    adjust_dTdV(dTdV, path, start_i, end_i, interf);
    direc = -dTdV; // p0
    
    // conjugate gradient search
    int iter_max = np*10;
    for (int iter=1; iter<=iter_max; iter++){
	new_time=line_min(path, direc, start_i, end_i, interf);
	if ((old_time-new_time)<=cg_tol){
	    return iter;
	}
	
	calc_dTdV(smesh, path, bs, new_dTdV, pp, Q, dQdu);
	adjust_dTdV(new_dTdV, path, start_i, end_i, interf);
	double gg=0.0, dgg=0.0;
	for (int i=1; i<=np; i++){
	    double x0=dTdV(i).x(),y0=dTdV(i).y();
	    double x1=new_dTdV(i).x(),y1=new_dTdV(i).y();
	    gg += x0*x0+y0*y0;
	    dgg += x1*x1+y1*y1;
	}
//	if (gg==0.0){ // 2000.5.8 
                      // this may be the cause for "-NaN" reported by
                      // Allegra. (due to division by extremely small gg)
	if (abs(gg)<1e-10){
	    return iter;
	}else{
	    dgg /= gg;
	}
	for (int i=1; i<=np; i++){
	    direc(i) = -new_dTdV(i)+dgg*direc(i);
	}
	dTdV = new_dTdV;
	old_time = new_time;
    }
    return -iter_max; // too many iterations
}

void BendingSolver2d::adjust_dTdV(Array1d<Point2d>& dTdV,
				  const Array1d<Point2d>& path,
				  const Array1d<int>& start_i,
				  const Array1d<int>& end_i,
				  const Array1d<const Interface2d*>& interf)
{
    for (int a=1; a<=interf.size(); a++){
	double slope = interf(a)->dzdx(path(start_i(a)).x());
	double total=0.0;
	for (int i=start_i(a); i<=end_i(a); i++){
	    total += dTdV(i).x()+slope*dTdV(i).y();
//	    total += dTdV(i).x();
	    dTdV(i).y() = 0.0;
	}
	total /= (end_i(a)-start_i(a)+1);
	for (int i=start_i(a); i<=end_i(a); i++){
	    dTdV(i).x() = total; // all points have the same direction now.
	}
    }
}

int BendingSolver2d::check_path_size(Array1d<Point2d>& path)
{
    int np=path.size();
    if (np<2) error("BendingSolver2d::invalid ray path");

    // the following is now taken care of by graph's pickPath()
//    if (np==2){
 	// since end points are fixed, this path wouldn't be
	// refined without an additional point in the middle
// 	Point2d p1=path(1),p2=path(2);
// 	path.resize(3);
// 	path(1)=p1; path(2)=0.5*(p1+p2); path(3)=p2;
//	np=3;
//    }

    return np;
}

int BendingSolver2d::check_path_size(Array1d<Point2d>& path,
				     const Array1d<int>& start_i,
				     const Array1d<int>& end_i,
				     const Array1d<const Interface2d*>& interf)
{
    int np=path.size();
    if (np<2) error("BendingSolver2d::invalid ray path");

    int nintf = interf.size();
    if (start_i.size() != nintf || end_i.size() != nintf)
	error("BendingSolver2d::refine - size mismatch in interface spec");
    for (int i=1; i<=nintf; i++){
	if (start_i(i) > end_i(i) ||
	    start_i(i) < 1 || start_i(i) > np ||
	    end_i(i) < 1 || end_i(i) > np)
	    error("BendingSolver2d::refine - invalid interface points");
    }
    return np;
}
    
int BendingSolver2d::refine(Array1d<Point2d>& path,
			    double& orig_time, double& new_time,
			    int nfac)
{
    if (nfac<1) error("BendingSolver2d::illegal nfac input");
    if (nfac==1){
	return refine(path,orig_time,new_time);
    }
    
    int np=path.size();
    Array1d<Point2d> old_path(np);
    old_path = path;

    int n2=(nfac-1)*(np-1);
    double frac=1.0/nfac;
    path.resize(np+n2);
    for (int i=1; i<np; i++){
	int orig_i = nfac*(i-1)+1;
	path(orig_i) = old_path(i);
	for (int j=1; j<nfac; j++){
	    double ratio=frac*j;
	    path(orig_i+j) = (1-ratio)*old_path(i)+ratio*old_path(i+1);
	}
    }
    path.back()=old_path.back();

    return refine(path,orig_time,new_time);
}

double BendingSolver2d::line_min(Array1d<Point2d>& path, const Array1d<Point2d>& direc)
{
    point_p = &path; // for f1dim()
    direc_p = &direc;

    double ax = 0.0;
    double xx = 1.0;
    double bx, fa, fx, fb, xmin;
    mnbrak(&ax, &xx, &bx, &fa, &fx, &fb, &BendingSolver2d::f1dim);
    double new_val = brent(ax,xx,bx,&xmin, &BendingSolver2d::f1dim);
    for (int i=1; i<=path.size(); i++){
	path(i) += xmin*direc(i);
    }
    return new_val;
}

double BendingSolver2d::f1dim(double x)
{
    for (int i=1; i<=point_p->size(); i++){
	new_point(i) = (*point_p)(i)+x*(*direc_p)(i);
    }
    double d = calcTravelTime(smesh, new_point, bs, new_pp, Q);
    return d;
}

double BendingSolver2d::line_min(Array1d<Point2d>& path, const Array1d<Point2d>& direc,
				 const Array1d<int>& start_i, const Array1d<int>& end_i,
				 const Array1d<const Interface2d*>& interf)
{
    point_p = &path; // for f1dim_interf()
    direc_p = &direc;
    start_i_p = &start_i;
    end_i_p = &end_i;
    interf_p = &interf;

//    double ax = 0.0;
//    double xx = 1.0;
    double ax = -0.1;
    double xx = 0.1;
    double bx, fa, fx, fb, xmin;
    mnbrak(&ax, &xx, &bx, &fa, &fx, &fb, &BendingSolver2d::f1dim_interf);
    double new_val = brent(ax,xx,bx,&xmin, &BendingSolver2d::f1dim_interf);
    for (int i=1; i<=path.size(); i++){
	path(i) += xmin*direc(i);
    }
    for (int a=1; a<=interf.size(); a++){
	double new_z = interf(a)->z(path(start_i(a)).x());
	for (int i=start_i(a); i<=end_i(a); i++){
	    path(i).y() = new_z;
	}
    }
    return new_val;
}

double BendingSolver2d::f1dim_interf(double x)
{
    for (int i=1; i<=point_p->size(); i++){
	new_point(i) = (*point_p)(i)+x*(*direc_p)(i);
    }
    for (int a=1; a<=interf_p->size(); a++){
	double new_z = (*interf_p)(a)->z(new_point((*start_i_p)(a)).x());
	for (int i=(*start_i_p)(a); i<=(*end_i_p)(a); i++){
	    new_point(i).y() = new_z;
	}
    }
    double d = calcTravelTime(smesh, new_point, bs, new_pp, Q);
    return d;
}

