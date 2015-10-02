/*
 * util.h
 *
 * Jun Korenaga, MIT/WHOI
 * January 1999
 */

#ifndef _JK_UTIL_H_
#define _JK_UTIL_H_

#include "geom.h"
#include "array.h"

const int MaxStr=512;
const char SepChars[]=" ";

void readBC(char* ifn, Array1d<RegularBC2d*>& velBC);
int countLines(const char *fn);
double intp1d(Array1d<double>& x, Array1d<double>& f, double xval);

#endif /* _JK_UTIL_H_ */
