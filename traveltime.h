/*
 * traveltime.h - traveltime related helper functions
 *
 * Jun Korenaga, MIT/WHOI
 * January 1999
 */

#ifndef _TOMO_TRAVELTIME_H_
#define _TOMO_TRAVELTIME_H_

#include <array.h> // from mconv
#include <geom.h>
#include "smesh.h"
#include "betaspline.h"

double calcTravelTime(const SlownessMesh2d& m, const Array1d<Point2d>& path,
		      const BetaSpline2d& bs,
		      const Array1d<const Point2d*>& pp,
		      Array1d<Point2d>& Q);
double calcTravelTime(const SlownessMesh2d& m,
		      const Array1d<Point2d>& path, int ndiv);
//double calcTravelTime2(const SlownessMesh2d& m, const Array1d<Point2d>& path);
void calc_dTdV(const SlownessMesh2d& m, const Array1d<Point2d>& path,
	       const BetaSpline2d& bs, Array1d<Point2d>& dTdV,
	       const Array1d<const Point2d*>& pp,
	       Array1d<Point2d>& Q, Array1d<Point2d>& dQdu);
void calc_dTdV2(const SlownessMesh2d& m, const Array1d<Point2d>& path,
		const BetaSpline2d& bs, Array1d<Point2d>& dTdV,
		const Array1d<const Point2d*>& pp,
		Array1d<Point2d>& Q, Array1d<Point2d>& dQdu);
void calc_dTdV3(const SlownessMesh2d& smesh, const Array1d<Point2d>& path,
 		const BetaSpline2d& bs, Array1d<Point2d>& dTdV,
 		const Array1d<const Point2d*>& pp,
 		Array1d<Point2d>& Q, Array1d<Point2d>& dQdu,
		const Array1d<int>& start_i,
		const Array1d<int>& end_i,
		const Array1d<const Interface2d*>& interf);
#endif /* _TOMO_TRAVELTIME_H_ */
