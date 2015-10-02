/*
 * util.cc
 *
 * Jun Korenaga, MIT/WHOI
 * January 1999
 */

#include <iostream>
#include <fstream>
#include "array.h"
#include "util.h"
#include "error.h"

int countLines(const char *fn)
{
    ifstream iin(fn);
    if (!iin){
	cerr << "countLines::can't open " << fn << '\n';
	exit(1);
    }
    int nline=0;
    char line[MaxStr];
    while(iin.getline(line,MaxStr)){
	nline++;
    }
    return nline;
}

double intp1d(Array1d<double>& x, Array1d<double>& f, double xval)
{
    if (x.size() != f.size())
	error("intp1d::size mismatch");

    if (xval < x(1)) return f.front();
    
    for (int i=1; i<x.size(); i++){
	if (xval>=x(i) && xval<x(i+1)){
	    double dx = xval-x(i);
	    double dX = x(i+1)-x(i);
	    double dF = f(i+1)-f(i);
	    return f(i)+dx/dX*dF;
	}
    }

    return f.back();
}
