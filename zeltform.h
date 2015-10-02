/*
 * zeltform.h
 *
 * Jun Korenaga, MIT/WHOI
 * January 1999
 */

#ifndef _TOMO_ZELFFORM_H_
#define _TOMO_ZELFFORM_H_

#include <array.h>
#include <geom.h>

struct ZNode2d {
    ZNode2d(){}

    Array1d<double> x;
    Array1d<double> val;
};

class TrapezoidCell2d {
public:
    TrapezoidCell2d(double, double,
		    double, double, double, double,
		    double, double, double, double);
    bool isIn(const Point2d&) const;
    double at(const Point2d&) const;
    void dumpCell(ostream&) const;
    
private:
    double x1, x2, s1, s2, b1, b2;
    double v1, v2, v3, v4;
    double c1, c2, c3, c4, c5, c6, c7;
};

class ZeltVelocityModel2d {
public:
    ZeltVelocityModel2d(char *fn);

    double at(double x, double z) const;
    void getTopo(int i, double dx,
		 Array1d<double>& x, Array1d<double>& topo) const;
    void dumpNodes(const char*) const;
    
private:
    int readLine(char *line, Array1d<double>& tmp);
    double interp(const ZNode2d*, double) const;
    
    Array1d< Array1d<ZNode2d*>* > node_p;
    Array1d<ZNode2d*> depth_node;
    Array1d<ZNode2d*> vupper_node;
    Array1d<ZNode2d*> vlower_node;

    Array1d<TrapezoidCell2d*> cells;
};


#endif /* _TOMO_ZELFFORM_H_ */

