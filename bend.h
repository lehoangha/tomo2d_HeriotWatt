/*
 * bend.h - ray-bending solver interface
 *
 * Jun Korenaga, MIT/WHOI
 * January 1999
 */

#ifndef _TOMO_BEND_H_
#define _TOMO_BEND_H_

#include <array.h> // from mconv
#include <geom.h>
#include "smesh.h"
#include "betaspline.h"
#include "interface.h"

class BendingSolver2d {
public:
    BendingSolver2d(const SlownessMesh2d& s, const BetaSpline2d& bs,
		    double tol1=1e-4, double tol2=1e-7);
    int refine(Array1d<Point2d>& path, double& orig_time, double& new_time);
    int refine(Array1d<Point2d>& path, double& orig_time, double& new_time, int nfac);
    int refine(Array1d<Point2d>& path, double& orig_time, double& new_time,
	       const Array1d<int>& start_i, const Array1d<int>& end_i,
	       const Array1d<const Interface2d*>& interf);
    double tolerance() const { return cg_tol; }
    
private:
    typedef double (BendingSolver2d::*PF1DIM)(double);

    int check_path_size(Array1d<Point2d>& path);
    int check_path_size(Array1d<Point2d>& path,
			const Array1d<int>& start_i,
			const Array1d<int>& end_i,
			const Array1d<const Interface2d*>& interf);
    void adjust_dTdV(Array1d<Point2d>& dTdV,
		     const Array1d<Point2d>& path,
		     const Array1d<int>& start_i,
		     const Array1d<int>& end_i,
		     const Array1d<const Interface2d*>& interf);
    double line_min(Array1d<Point2d>& path, const Array1d<Point2d>& direc);
    double line_min(Array1d<Point2d>& path, const Array1d<Point2d>& direc,
		    const Array1d<int>& start_i, const Array1d<int>& end_i,
		    const Array1d<const Interface2d*>& interf);
    void mnbrak(double *ax, double *bx, double *cx,
		double *fa, double *fb, double *fc, PF1DIM);
    double brent(double ax, double bx, double cx, double *xmin, PF1DIM);
    double f1dim(double x);
    double f1dim_interf(double x);
		    
    const SlownessMesh2d& smesh;
    const BetaSpline2d& bs;
    const int nintp;
    const double cg_tol, brent_tol, eps;

    Array1d<Point2d> Q, dQdu;
    Array1d<const Point2d*> pp, new_pp;
    const Array1d<Point2d> *point_p, *direc_p;
    const Array1d<int> *start_i_p, *end_i_p;
    const Array1d<const Interface2d*> *interf_p;
    Array1d<Point2d> new_point;
    Array1d<Point2d> dTdV, new_dTdV, direc;
};

#endif /* _TOMO_BEND_H_ */
