/*
 * interface.h
 *
 * Jun Korenaga, MIT/WHOI
 * January 1999
 */

#ifndef _TOMO_INTERFACE_H_
#define _TOMO_INTERFACE_H_

#include <array.h>
#include "smesh.h"

class Interface2d {
public:
    Interface2d(){}
    Interface2d(const SlownessMesh2d&);
    Interface2d(const char*);

    double z(double x) const;
    double dzdx(double x) const;
    void locateInSegment(double, int&, int&) const;
    double x(int) const;
    double xmin() const;
    double xmax() const;
    int numNodes() const;
    void set(const Array1d<double>&);
    void get(Array1d<double>&) const;

    friend ostream&
    operator<<(ostream&, const Interface2d&);

private:
    void calc_slope();
    static const double eps;
    Array1d<double> xpos, zpos, slope;
};

#endif

