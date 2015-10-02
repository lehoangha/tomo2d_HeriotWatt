/*
 * jgrav.cc
 *
 * Jun Korenaga, MIT/WHOI
 * July 1999
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include "jgrav.h"
#include "d_nr.h"

const double AddonGravityInversion2d::G = 6.67e-8;
const double AddonGravityInversion2d::PI = 3.141592653;
const double AddonGravityInversion2d::v_mantle = 8.0;
const double AddonGravityInversion2d::rho_mantle = 3.3;

AddonGravityInversion2d::AddonGravityInversion2d(const SlownessMesh2d& m,
						 const Interface2d& f,
						 const RegularDomain2d& d,
						 double _dx, double _dz,
						 double x0, double x1)
    : smesh(m), moho(f),
      xmin(d.x_min()), xmax(d.x_max()), zmin(d.y_min()), zmax(d.y_max()),
      dx(_dx), dz(_dz), nx(int((xmax-xmin)/dx)), nz(int((zmax-zmin)/dz)),
      is_flagged(false), getContBound(false), getOceanUBound(false),
      getOceanLBound(false), getSedBound(false),
      dvdp(0.0), dvdt(0.0), drdp(0.0), drdt(0.0), dTdz(0.0),
      ref_x0(x0), ref_x1(x1),
      initKernel(true), cutoff_range(50.0), cutoff_K(1e-3)
{
    x.resize(nx);
    z.resize(nz);
    for (int i=1; i<=nx; i++) x(i) = xmin+0.5*dx+(i-1)*dx; // pixel registration
    for (int i=1; i<=nz; i++) z(i) = zmin+0.5*dz+(i-1)*dz;
    rho.resize(nz); flag.resize(nz); 
    for (int k=1; k<=nz; k++){
	rho(k).resize(nx); flag(k).resize(nx); 
    }

    int pow2 = int(log(double(nx*4))/log(2.0));
    if (4*nx==int(pow(2.0,pow2))){
	nfft = nx*4;
    }else{
	nfft = int(pow(2.0,pow2+1));
    }
    hnfft = nfft/2;
    npad = (hnfft-nx)/2;
    kx.resize(hnfft);
    coeff.resize(hnfft);
    grav.resize(nfft);
    rhofft.resize(nz);
    for (int k=1; k<=nz; k++) rhofft(k).resize(nfft);

    ndepth = moho.numNodes();
    twoG = 2.0*G;
    PI05 = 0.5*PI;
}

void AddonGravityInversion2d::setDerivative(double vp, double vt,
					    double rp, double rt, double tz)
{
    dvdp = vp; dvdt = vt; drdp = rp; drdt = rt; dTdz = tz;
}

void AddonGravityInversion2d::defineContinentRegion(const Interface2d* top,
						    const Interface2d* bot,
						    int conv)
{ getContBound = true; contup = top; contlo = bot; icontconv = conv; }

void AddonGravityInversion2d::defineUpperOceanicRegion(const Interface2d* top,
						       const Interface2d* bot,
						       int conv)
{ getOceanUBound = true; oceanUup = top; oceanUlo = bot; ioceanUconv = conv; }

void AddonGravityInversion2d::defineLowerOceanicRegion(const Interface2d* top,
						       const Interface2d* bot,
						       int conv)
{ getOceanLBound = true; oceanLup = top; oceanLlo = bot; ioceanLconv = conv; }

void AddonGravityInversion2d::defineSedimentaryRegion(const Interface2d* top,
						      const Interface2d* bot,
						      int conv)
{ getSedBound = true; sedup = top; sedlo = bot; isedconv = conv; }

void AddonGravityInversion2d::setVel()
{
    Index2d guess=smesh.nodeIndex(smesh.nearest(Point2d(x(1),z(1))));
    for (int i=1; i<=nx; i++){
	double tmpx = x(i);
	double mohoz = moho.z(tmpx);
	for (int k=1; k<=nz; k++){
	    double tmpz = z(k);
	    double vel;
	    if (tmpz<=mohoz){
		vel = 1.0/smesh.at(Point2d(tmpx,tmpz),guess);
	    }else{
		vel = v_mantle;
	    }
	    rho(k)(i) = vel;
	}
    }
}


void AddonGravityInversion2d::setFlag()
{
	for (int k=1; k<=nz; k++){
	    double tmpz = z(k);
	    for (int i=1; i<=nx; i++){
		double tmpx = x(i);
		if (rho(k)(i) <=1.5){
		    flag(k)(i) = Water;
		}else if (getContBound && isIn(tmpx,tmpz,contup,contlo)){
		    flag(k)(i) = Cont;
		}else if (getOceanUBound && isIn(tmpx,tmpz,oceanUup,oceanUlo)){
		    flag(k)(i) = OceanU;
		}else if (getOceanLBound && isIn(tmpx,tmpz,oceanLup,oceanLlo)){
		    flag(k)(i) = OceanL;
		}else if (getSedBound && isIn(tmpx,tmpz,sedup,sedlo)){
		    flag(k)(i) = Sed;
		}else{
		    flag(k)(i) = Mantle;
		}
	    }
	}
}

void AddonGravityInversion2d::vel2rho()
{
    for (int k=1; k<=nz; k++){
	double tmpz = z(k);
	for (int i=1; i<=nx; i++){
	    switch(flag(k)(i)){
	    case Water:
		rho(k)(i) = 1.0;
		break;
	    case Cont:
		rho(k)(i) = cont_vel2rho(rho(k)(i), icontconv);
		break;
	    case OceanU:
	    {
		double origvel, dp, dT;
		origvel = rho(k)(i);
		dT = 25-tmpz*dTdz;
		dp = 1000-tmpz*33.3333; /* km -> MPa */
		origvel += (dT*dvdt+dp*dvdp); /* now at room temperature and 10 kbar */
		rho(k)(i) = oceanU_vel2rho(origvel, ioceanUconv);
		rho(k)(i) -= (dT*drdt+dp*drdp); /* convert to in situ value */
		break;		
	    }
	    case OceanL:
	    {
		double origvel, dp, dT;
		origvel = rho(k)(i);
		dT = 25-tmpz*dTdz;
		dp = 1000-tmpz*33.3333; /* km -> MPa */
		origvel += (dT*dvdt+dp*dvdp); /* now at room temperature and 10 kbar */
		rho(k)(i) = oceanL_vel2rho(origvel, ioceanLconv);
		rho(k)(i) -= (dT*drdt+dp*drdp); /* convert to in situ value */
		break;
	    }
	    case Sed:
		rho(k)(i) = sed_vel2rho(rho(k)(i), isedconv);
		break;
	    case Mantle:
		rho(k)(i) = rho_mantle;
		break;
	    default:
		cerr << "AddonGravityInversion2d::vel2rho unrecognized flag ("
		     << flag(k)(i) << ") detected\n";
		exit(1);
		break;
	    }
	}
    }
}

double AddonGravityInversion2d::cont_vel2rho(double vel, int mode) const
{
    double tmp;
    
    switch(mode){
    case 1:
	// Christensen+Mooney, 1995 (z=20km), nonlinear regression 
	tmp = 5.055-14.094/vel; 
	break;
    case 2:
	// Christensen+Mooney, 1995 (z=20km) 
	// linear regression for all rocks except volcanic and monomineralic rocks 
	tmp = 0.444+0.375*vel; 	
	break;
    case 100:
	tmp = 2.6;
	break;
    default:
	cerr <<  "[cont_vel2rho]: illegal mode (", mode, ") detected.\n";
	exit(1);
	break;
    }

    if (mode>0){
	if (tmp <= 0) tmp = 2.5; // minimum value (safety) 
	if (tmp >=4.0) tmp = 4.0; // maximum value 
    }
    return tmp;
}

double AddonGravityInversion2d::sed_vel2rho(double vel, int mode) const
{
    double tmp;
    switch(mode){
    case 1:
	if (vel<=1.5){
	    tmp = 1.0; // water
	}else{
	    tmp = 1+1.18*pow((vel-1.5), 0.22); // sediments 
	}
	break;
    case 100:
	tmp = 1.5;
	break;
    default:
	cerr << "[sed_vel2rho]: illegal mode (", mode, ") detected.\n";
	exit(1);
	break;
    }
    if (tmp >= 2.6 && mode>0) tmp = 2.6;  // for safety 
    return tmp;
}

double AddonGravityInversion2d::oceanU_vel2rho(double vel, int mode) const
{
    double tmp;
    switch(mode){
    case 1:
	tmp = 3.81-6.0/vel; // layer 2 (Carlson+Herrick,1990) 
	break;
    case 100:
	tmp = 2.0; 
	break;
    default:
	cerr << "[oceanU_vel2rho]: illegal mode (", mode, ") detected.\n";
	exit(1);
	break;
    }
    if (tmp <= 1.0) tmp = 1.0; // for safety 
    return tmp;
}

double AddonGravityInversion2d::oceanL_vel2rho(double vel, int mode) const
{
    double tmp;
    switch(mode){
    case 1:
	// Birch, 1961, solution 7 (diabase,gabbro,eclogite)
	// at 10 kbar & room temperature 
	tmp = (vel+1.0)/2.67; 
	break;
    case 100:
	tmp = 3.0;
	break;
    default:
	cerr << "[oceanL_vel2rho]: illegal mode (", mode, ") detected.\n";
	exit(1);
	break;
    }
    return tmp;
}

bool AddonGravityInversion2d::isIn(double px, double pz,
				   const Interface2d* top, const Interface2d* bot) const
{
    if (px>=top->xmin() && px<=top->xmax()
	&& px>=bot->xmin() && px<=bot->xmax()){
	double topz = top->z(px);
	double botz = bot->z(px);
	if (pz>=topz && pz<=botz) return true;
    }
    return false;
}

void AddonGravityInversion2d::padAndMirror()
{
    for (int k=1; k<=nz; k++){
	 for (int i=1; i<=npad; i++){
	      rhofft(k)(i) = rho(k)(1);
	 }
	 int j=1;
	 for (int i=npad+1; i<=npad+nx; i++,j++){
	      rhofft(k)(i) = rho(k)(j);
	 }
	 for (int i=npad+nx+1; i<=hnfft; i++){
	      rhofft(k)(i) = rho(k)(nx);
	 }
	 j=hnfft;
	 for (int i=hnfft+1; i<=nfft; i++,j--){
	      rhofft(k)(i) = rhofft(k)(j);
	 }
    }
}

void AddonGravityInversion2d::calcGravity(double z0,
					  const Array1d<double>& xobs,
					  Array1d<double>& gravobs)
{
    for (int i=1; i<=xobs.size(); i++){
	double tmpx = xobs(i);
	if (tmpx<xmin && tmpx>xmax){
	    cerr << "AddonGravityInversion2d::calcGravity::"
		 << "out of bound is detected for xobs at i=" << i << '\n';
	    exit(1);
	}
    }

    setVel();
    setFlag();
    vel2rho();
    padAndMirror();
    
    double dkx = 2*PI/(nfft*dx);	
    double pi2g = 2*PI*G;
    for (int i=1; i<=hnfft; i++){
	kx(i) = i*dkx;
	coeff(i) = pi2g*(1-exp(-kx(i)*dz))/kx(i);
    }

    grav = 0.0; 
    for (int k=1; k<=nz; k++){
	double tmpz = z(k)-z0;
	if (tmpz<0) continue;	// ignore mass above the observation level 
				// i.e., calculate Bouguer gravity anomaly for on-land
	d_realft(rhofft(k).toRecipe(), nfft, 1);
	int i=1;
	for (int j=3; j<nfft; j+=2,i++){
	    double fac = exp(-kx(i)*tmpz)*coeff(i);
	    grav(j) += fac*rhofft(k)(j);
	    grav(j+1) += fac*rhofft(k)(j+1);
	}
	grav(1) += pi2g*rhofft(k)(1)*dz; // infinite slab formula for k=0 
	grav(2) += exp(-kx(hnfft)*tmpz)*coeff(hnfft)*rhofft(k)(2);
    }
    d_realft(grav.toRecipe(), nfft, -1);
    grav /= hnfft;

    int nin=0;
    double ave=0;
    for (int i=1; i<=nx; i++){
	if (x(i)>=ref_x0 && x(i)<=ref_x1){
	    ave += grav(i+npad);
	    nin++;
	}
    }
    ave /= nin;
    
    // linear interpolation for given xobs
    for (int i=1; i<=xobs.size(); i++){
	double tmpi = (xobs(i)-x(1))/dx+1;
	int i0 = int(tmpi); // round off
	double r = tmpi-i0;
	gravobs(i) = (1-r)*grav(i0+npad)+r*grav(i0+1+npad);
	gravobs(i) -= ave;
	gravobs(i) *= 1e8; // in mGal
    }
}

void AddonGravityInversion2d::calcGravityKernel(double z0, double x0,
						map<int,double>& B_i)
{
    int smesh_nx = smesh.Nx();
    int smesh_nz = smesh.Nz();
    int nvelnodes = smesh.numNodes();
    double smesh_xmin = smesh.xmin();
    double smesh_xmax = smesh.xmax();

    if (initKernel){
	smesh.nodalCellVolume(cell_dx, cell_dz, cell_center);
	    
	depth_dx.resize(ndepth);
	depth_dx(1) = moho.x(2)-moho.x(1);
	for (int i=2; i<ndepth; i++) depth_dx(i) = 0.5*(moho.x(i+1)-moho.x(i-1));
	depth_dx(ndepth) = moho.x(ndepth)-moho.x(ndepth-1);

	kernel_flag.resize(nvelnodes);
	setKernelFlag();
	initKernel = false;
    }else{
        setKernelFlag();
    }

    // velocity node's kernel
    for (int i=1; i<=smesh_nx; i++){
	double tmpx = cell_center(smesh.nodeIndex(i,1)).x();
	double range = abs(tmpx-x0);
	if (range > cutoff_range) continue;
	for (int k=1; k<=smesh_nz; k++){
	    int inode = smesh.nodeIndex(i,k);
	    double drhods;
	    switch(kernel_flag(inode)){
	    case Water:
	    case Sed:
	    case Mantle:
		drhods = 0.0;
		break;
	    case Cont:
		drhods = cont_drhods(cell_center(inode), icontconv);
		break;
	    case OceanU:
		drhods = oceanU_drhods(cell_center(inode), ioceanUconv);
		break;
	    case OceanL:
		drhods = oceanL_drhods(cell_center(inode), ioceanLconv);
		break;
	    default:
		cerr << "AddonGravityInversion2d::calcGravityKernel-unrecognized flag ("
		     << kernel_flag(inode) << ") detected\n";
		exit(1);
		break;
	    }
	    if (drhods==0.0) continue;
	    
	    double xx = cell_center(inode).x()-x0;
	    double zz = cell_center(inode).y()-z0;
	    double K1 = zz*cell_dx(inode)/(xx*xx+zz*zz);
	    if (i==1){ // left edge
		K1 += (PI05+atan((smesh_xmin-x0)/zz));
	    }else if (i==smesh_nx){ // right edge
		K1 += (PI05+atan((x0-smesh_xmax)/zz));
	    }
	    K1 *= twoG*drhods*cell_dz(inode)*1e8;
	    if (abs(K1)>=cutoff_K) B_i[inode] = K1;
//	    B_i[inode] = zz/(xx*xx+zz*zz);
	}
    }

    // depth node's kernel
    for (int i=1; i<=ndepth; i++){
	int inode = i+nvelnodes;
	double xd = moho.x(i);
	double zd = moho.z(xd);
	int ixd = int((xd-xmin)/dx+1);
	int izd = int((zd-zmin)/dz+1);
	if (ixd>=1 && ixd<=nx && izd>=1 && izd<=nz){
	    double rho_c = rho(izd)(ixd); // assumes rho()() is already set by calcGravity()
	    double xx = xd-x0;
	    double zz = zd-z0;
	    double K2 = twoG*(rho_c-rho_mantle)*zz*depth_dx(i)/(xx*xx+zz*zz)*1e8;
//	    if (abs(K2)>=cutoff_K) B_i[inode] = K2;
	    B_i[inode] = K2;
	}
    }
}

void AddonGravityInversion2d::setThreshold(double range, double val)
{
    cutoff_range = range;
    cutoff_K = val;
}

void AddonGravityInversion2d::setKernelFlag()
{
    for (int i=1; i<=kernel_flag.size(); i++){
	Point2d p=smesh.nodePos(i);
	double tmpx = p.x();
	double tmpz = p.y();
	if (getContBound && isIn(tmpx,tmpz,contup,contlo)){
	    kernel_flag(i) = Cont;
	}else if (getOceanUBound && isIn(tmpx,tmpz,oceanUup,oceanUlo)){
	    kernel_flag(i) = OceanU;
	}else if (getOceanLBound && isIn(tmpx,tmpz,oceanLup,oceanLlo)){
	    kernel_flag(i) = OceanL;
	}else if (getSedBound && isIn(tmpx,tmpz,sedup,sedlo)){
	    kernel_flag(i) = Sed;
	}else{
	    kernel_flag(i) = Mantle;
	}
    }
}

double AddonGravityInversion2d::cont_drhods(const Point2d& p, int mode) const
{
    double tmp;
    switch(mode){
    case 1:
	// Christensen+Mooney, 1995 (z=20km), nonlinear regression 
	tmp = -14.094;
	break;
    case 2:
    {
	// Christensen+Mooney, 1995 (z=20km) 
	// linear regression for all rocks except volcanic and monomineralic rocks
	double v = 1.0/smesh.at(p);
	tmp = -0.375*v*v;
	break;
    }
    case 100:
	tmp = 0.0;
	break;
    default:
	cerr <<  "[cont_drhods]: illegal mode (", mode, ") detected.\n";
	exit(1);
	break;
    }
    return tmp;
}

double AddonGravityInversion2d::oceanU_drhods(const Point2d& p, int mode) const
{
    double tmp;
    switch(mode){
    case 1:
	tmp = -6.0; // layer 2 (Carlson+Herrick,1990) 
	break;
    case 100:
	tmp = 0.0; 
	break;
    default:
	cerr << "[oceanU_drhods]: illegal mode (", mode, ") detected.\n";
	exit(1);
	break;
    }
    return tmp;
}

double AddonGravityInversion2d::oceanL_drhods(const Point2d& p, int mode) const
{
    double tmp;
    switch(mode){
    case 1:
    {
	// Birch, 1961, solution 7 (diabase,gabbro,eclogite)
	// at 10 kbar & room temperature 
	double v = 1.0/smesh.at(p);
	tmp = -v*v/2.67;
	break;
    }
    case 100:
	tmp = 0;
	break;
    default:
	cerr << "[oceanL_drhods]: illegal mode (", mode, ") detected.\n";
	exit(1);
	break;
    }
    return tmp;
}
