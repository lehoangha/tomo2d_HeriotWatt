/*
 * lsqr.h - LSQR definition
 *
 * Jun Korenaga, MIT/WHOI
 * January 1999
 */

#ifndef _TOMO_LSQR_H_
#define _TOMO_LSQR_H_

#include <array.h>
#include "sparse_rect.h"

int iterativeSolver_LSQR(const SparseRectangular& A, const Array1d<double>& b,
			 Array1d<double>& x, double ATOL, int& itermax, double& chi);
double normalize(Array1d<double>& x);

#endif /* _TOMO_LSQR_H_ */
