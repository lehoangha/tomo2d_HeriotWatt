/*
 * jgrav.h - for joint inversion of travel times and gravity anomalies
 *
 * Jun Korenaga, MIT/WHOI
 * July 1999
 */

#ifndef _TOMO_JGRAV_H_
#define _TOMO_JGRAV_H_

#include <map>
#include <cmath>
#include "geom.h"
#include "smesh.h"
#include "interface.h"

class AddonGravityInversion2d {
public:
    AddonGravityInversion2d(const SlownessMesh2d&, const Interface2d&,
			    const RegularDomain2d&, double, double, double, double);
    
    void calcGravity(double, const Array1d<double>&, Array1d<double>&);
    void calcGravityKernel(double, double, map<int,double>&);
    void setThreshold(double, double);

    void defineContinentRegion(const Interface2d*, const Interface2d*, int);
    void defineUpperOceanicRegion(const Interface2d*, const Interface2d*, int);
    void defineLowerOceanicRegion(const Interface2d*, const Interface2d*, int);
    void defineSedimentaryRegion(const Interface2d*, const Interface2d*, int);
    void setDerivative(double, double, double, double, double);
    
private:
    void setVel();
    void setFlag();
    void setKernelFlag();
    void vel2rho();
    double cont_vel2rho(double, int) const;
    double oceanU_vel2rho(double, int) const;
    double oceanL_vel2rho(double, int) const;
    double sed_vel2rho(double, int) const;
    double cont_drhods(const Point2d&, int) const;
    double oceanU_drhods(const Point2d&, int) const;
    double oceanL_drhods(const Point2d&, int) const;
    bool isIn(double, double, const Interface2d*, const Interface2d*) const;
    void padAndMirror();

    static const double G, PI, v_mantle, rho_mantle; 
    
    const SlownessMesh2d& smesh;
    const Interface2d& moho;

    const double xmin, xmax, zmin, zmax;
    const double ref_x0, ref_x1;
    const double dx, dz;
    const int nx, nz;
    Array1d<double> x, z;
    Array1d< Array1d<double> > rho;
    enum { Water, Cont, Sed, OceanU, OceanL, Mantle};
    bool is_flagged;
    Array1d< Array1d<int> > flag;

    bool getContBound, getOceanUBound, getOceanLBound, getSedBound;
    const Interface2d *contup, *contlo, *oceanUup, *oceanUlo;
    const Interface2d *oceanLup, *oceanLlo, *sedup, *sedlo;
    int icontconv, ioceanUconv, ioceanLconv, isedconv;
    double dvdp, dvdt, drdp, drdt, dTdz;

    int nfft, hnfft, npad;
    Array1d<double> kx, coeff, grav; 
    Array1d< Array1d<double> > rhofft;

    bool initKernel;
    int ndepth;
    double twoG, PI05, cutoff_range, cutoff_K;
    Array1d<int> kernel_flag;
    Array1d<double> cell_dx, cell_dz, depth_dx;
    Array1d<Point2d> cell_center;
};

#endif /* _TOMO_JGRAV_H_ */
