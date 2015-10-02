/*
 * geom.cc
 *
 * Jun Korenaga, MIT/WHOI
 * January 1999
 */

#include "error.h"
#include "geom.h"

// Point2d
Point2d::Point2d(double x, double y)
{
    v[0] = x;
    v[1] = y;
}

double Point2d::inner_product(const Point2d& p) const
{
    double ip = v[0]*p.v[0]+v[1]*p.v[1];
    return ip;
}

// RegularDomain2d
double RegularDomain2d::eps=1e-6;

RegularDomain2d::RegularDomain2d(double x1, double x2, double y1, double y2)
{
    if (x1>x2 || y1>y2) error("RegularDomain2d::BadInput"); // i range = 0-1
    
    vmin[0] = x1;
    vmax[0] = x2;
    vmin[1] = y1;
    vmax[1] = y2;
}

// RegularBC2d
RegularBC2d::RegularBC2d(int i, double x1, double x2, double y1, double y2, double d)
    : idof(i), RegularDomain2d(x1,x2,y1,y2), val_(d)
{
    // the following is commented out for VectorMesh2d_plus class
//    if (idof <=0 || idof > 2) error("RegularBC2d::wrong input for degree of freedom");
    if (idof <=0) error("RegularBC2d::wrong input for degree of freedom");
}


