/*
 * corrlen.h - correlation length functions' interface
 *
 * Jun Korenaga, MIT/WHOI
 * February 1999
 */

#ifndef _TOMO_CORRLEN_H_
#define _TOMO_CORRLEN_H_

#include <array.h>
#include <geom.h>
#include <index.h>

class CorrelationLength1d {
public:
    CorrelationLength1d(const char *fn);
    double at(double) const;

private:
    const double eps;
    Array1d<double> xpos, val;
};

class CorrelationLength2d {
public:
    CorrelationLength2d(const char *fn);
    void at(const Point2d&, double&, double&) const;

private:
    void upperleft(const Point2d& pos, Index2d& index) const;
    void calc_local(const Point2d& pos, int, int,
		    double&, double&, double&, double&) const;

    int nx, nz;
    Array2d<double> hgrid, vgrid;
    Array1d<double> xpos, topo, zpos;
    Array1d<double> rdx_vec, rdz_vec, b_vec;
};

class DampingWeight2d {
public:
    DampingWeight2d(const char *fn);
    void at(const Point2d&, double&) const;

private:
    void upperleft(const Point2d& pos, Index2d& index) const;
    void calc_local(const Point2d& pos, int, int,
		    double&, double&, double&, double&) const;

    int nx, nz;
    Array2d<double> wgrid;
    Array1d<double> xpos, topo, zpos;
    Array1d<double> rdx_vec, rdz_vec, b_vec;
};

#endif /* _TOMO_CORRLEN_H_ */
