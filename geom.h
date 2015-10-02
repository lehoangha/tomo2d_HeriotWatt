/*
 * geom.h
 *
 * Jun Korenaga, MIT/WHOI
 * January 1999
 */

#ifndef _MCONV_GEOM_H_
#define _MCONV_GEOM_H_

#include <cmath>

// Point2d
class Point2d {
public:
    Point2d(){}
    Point2d(double, double);

    double x() const { return v[0]; }
    double y() const { return v[1]; }
    double& x() { return v[0]; }
    double& y() { return v[1]; }
    void   x(double d) { v[0] = d; }
    void   y(double d) { v[1] = d; }
    void   set(double dx, double dy) { v[0] = dx; v[1] = dy; }
    double distance(const Point2d& p) const;
    double norm() const { return sqrt(v[0]*v[0]+v[1]*v[1]); }
    double inner_product(const Point2d& p) const;

    // unary operators
    Point2d operator-();
    
    // binary operators
    Point2d& operator+=(const Point2d&);
    Point2d& operator-=(const Point2d&);
    Point2d& operator*=(double);
    Point2d& operator/=(double);
    
private:
    double v[2];
};

inline double Point2d::distance(const Point2d& p) const
{
    double dx = p.v[0]-v[0];
    double dy = p.v[1]-v[1];
    return sqrt(dx*dx+dy*dy);
}

inline
Point2d Point2d::operator-()
{
    Point2d neg;

    neg.v[0] = -v[0];
    neg.v[1] = -v[1];

    return neg;
}

inline
Point2d& Point2d::operator+=(const Point2d& a)
{
    v[0] += a.v[0];
    v[1] += a.v[1];
    return *this;
}

inline
Point2d& Point2d::operator-=(const Point2d& a)
{
    v[0] -= a.v[0];
    v[1] -= a.v[1];
    return *this;
}

inline
Point2d& Point2d::operator*=(double a)
{
    v[0] *= a;
    v[1] *= a;
    return *this;
}

inline
Point2d& Point2d::operator/=(double a)
{
    v[0] /= a;
    v[1] /= a;
    return *this;
}

// nonmember functions
inline
Point2d operator+(const Point2d& a, const Point2d& b)
{
    Point2d c=a;
    return c+=b;
}

inline
Point2d operator-(const Point2d& a, const Point2d& b)
{
    Point2d c=a;
    return c-=b;
}

inline
Point2d operator*(const Point2d& a, double val)
{
    Point2d c=a;
    return c*=val;
}

inline
Point2d operator*(double val, const Point2d& a)
{
    return operator*(a,val);
}

inline
Point2d operator/(const Point2d& a, double val)
{
    Point2d c=a;
    return c/=val;
}

// RegularDomain2d
class RegularDomain2d {
  public:
    RegularDomain2d(){ }
    RegularDomain2d(double, double, double, double); // 2-D 

    bool inXRange(double) const;
    bool inYRange(double) const;
    bool inDomain(const Point2d&) const;
    
    double x_min() const { return vmin[0]; }
    double x_max() const { return vmax[0]; }
    double y_min() const { return vmin[1]; }
    double y_max() const { return vmax[1]; }

    void x_min(double d) { vmin[0] = d; }
    void x_max(double d) { vmax[0] = d; }
    void y_min(double d) { vmin[1] = d; }
    void y_max(double d) { vmax[1] = d; }

  private:
    static double eps;

    double vmin[2];
    double vmax[2];
};

inline bool
RegularDomain2d::inXRange(double d) const
{
    if (d >= x_min()-eps && d <= x_max()+eps) return true;
    return false;
}

inline bool
RegularDomain2d::inYRange(double d) const
{
    if (d >= y_min()-eps && d <= y_max()+eps) return true;
    return false;
}

inline bool
RegularDomain2d::inDomain(const Point2d& p) const
{
    if (inXRange(p.x()) && inYRange(p.y())) return true;
    return false;
}

class RegularBC2d : public RegularDomain2d {
  public:
    RegularBC2d(int, double, double, double, double, double); // 2-D
	
    bool   inBoundary(const Point2d& p) const { return inDomain(p); }

    int    iDegOfFreedom() const { return idof; }
    double val() const           { return val_; }
    void   set_iDegOfFreedom(int i) { idof = i; }
    void   set_val(double d)        { val_ = d; }

  private:
    int idof;			// index for deg of freedom
    double val_;
};

#endif /* _MCONV_GEOM_H_ */
