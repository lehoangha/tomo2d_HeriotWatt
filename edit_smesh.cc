/*
 * edit_smesh.cc
 *
 * usage: edit_smesh smesh_file -C<cmd>
 *
 *        -C<cmd>
 *          cmd = 'a' - set 1-D average
 *                'p<smesh>' - paste <smesh> on the original smesh
 *                'P<1dprof>' - paste 1d-prof 
 *                'sh/v' - apply Gaussian smoothing operator with an window of h(horizontal)
 *                         and v(vertical)
 *                'rmx/mz' - refine mesh by mx for x-direction and by mz for z-direction
 *                'cA/h/v' - add checkerboard pattern
 *                           (A:amplitude [%], h:horizontal, v:vertical)
 *                'dA/xmin/xmax/zmin/zmax' - add rectangular anomaly
 *                           (A:amplitude [%])
 *                'gA/x0/z0/Lh/Lv' - add gaussian anomaly
 *                'l' - remove low velocity zone
 *                'R<seed>/A/nrand' - randomize the velocity field
 *                'S<seed>/A/xmin/xmax/dx/zmin/zmax/dz' 
 *                'G<seed>/A/N/xmin/xmax/zmin/zmax'
 *                'm<v>/mohofile' - set sub-Moho velocity as <v>
 *
 *        -Lvcorr.file
 *        -Uupperbound.file
 *
 * Jun Korenaga, MIT/WHOI
 * February 1999
 */

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <array.h>
#include <util.h>
#include "smesh.h"
#include "corrlen.h"
#include "interface.h"

int main(int argc, char **argv)
{
    bool getMesh=false, getCmd=false, setAverage=false;
    bool pasteSmesh=false, smooth=false, checkerboard=false;
    bool paste1dprof=false;
    bool refineMesh=false, rectangular=false, removeLVZ=false;
    bool gaussian=false, randomize=false, randomize2=false, randomize3=false;
    bool verbose=false, err=false;
    bool getCorr=false, getBound=false;
    char *sfn, *pfn, *corrfn, *prof1dfn;
    double Lh, Lv, A, ch, cv;
    double dxmin, dxmax, dzmin, dzmax;
    double x0, z0;
    int mx, mz, seed, nrand, N;
    double rand_xmax, rand_xmin, rand_dx, rand_zmax, rand_zmin, rand_dz;
    bool subMohoVel=false;
    double vsubmoho;
    char uboundfn[MaxStr];
    
    for (int i=1; i<argc; i++){
	if (argv[i][0] == '-'){
	    switch(argv[i][1]){
	    case 'C':
		getCmd = true;
		switch(argv[i][2]){
		case 'a':
		    setAverage = true;
		    break;
		case 'p':
		    pasteSmesh = true;
		    pfn = &argv[i][3];
		    break;
		  case 'P':
		    paste1dprof = true;
		    prof1dfn = &argv[i][3];
		    break;
		case 's':
		    if (sscanf(&argv[i][3], "%lf/%lf", &Lh, &Lv)==2){
			smooth = true;
		    }else{
			cerr << "invalid -Cs option\n";
			err = true;
		    }
		    break;
		case 'c':
		    if (sscanf(&argv[i][3], "%lf/%lf/%lf", &A, &ch, &cv)==3){
			checkerboard = true;
		    }else{
			cerr << "invalid -Cc option\n";
			err = true;
		    }
		    break;
		case 'd':
		    if (sscanf(&argv[i][3], "%lf/%lf/%lf/%lf/%lf",
			       &A, &dxmin, &dxmax, &dzmin, &dzmax)==5){
			rectangular=true;
		    }else{
			cerr << "invalid -Cd option\n";
			err = true;
		    }
		    break;
		case 'g':
		    if (sscanf(&argv[i][3], "%lf/%lf/%lf/%lf/%lf",
			       &A, &x0, &z0, &Lh, &Lv)==5){
			gaussian=true;
		    }else{
			cerr << "invalid -Cg option\n";
			err = true;
		    }
		    break;
		case 'r':
		    if (sscanf(&argv[i][3], "%d/%d", &mx, &mz)==2){
			refineMesh = true;
		    }else{
			cerr << "invalid -Cr option\n";
			err = true;
		    }
		    break;
		case 'l':
		    removeLVZ = true;
		    break;
		case 'R':
		    if (sscanf(&argv[i][3], "%d/%lf/%d",
			       &seed, &A, &nrand)==3){
			randomize = true;
		    }else{
			cerr << "invalid -CR option\n";
			err = true;
		    }
		    break;
		  case 'S':
		    if (sscanf(&argv[i][3], "%d/%lf/%lf/%lf/%lf/%lf/%lf/%lf",
			       &seed, &A, &rand_xmin, &rand_xmax, &rand_dx,
			       &rand_zmin, &rand_zmax, &rand_dz)==8){
		      randomize2=true;
		    }else{
		      cerr << "invalid -CS option\n";
		      err = true;
		    }
		    break;
		  case 'G':
		    if (sscanf(&argv[i][3], "%d/%lf/%d/%lf/%lf/%lf/%lf",
			       &seed, &A, &N, &rand_xmin, &rand_xmax,
			       &rand_zmin, &rand_zmax)==7){
		      randomize3=true;
		    }else{
		      cerr << "invalid -CG option\n";
		      err = true;
		    }
		    break;
		case 'm':
		    if (sscanf(&argv[i][3], "%lf/%[^/]", &vsubmoho, uboundfn)==2){
			subMohoVel = true;
		    }else{
			cerr << "invalid -Cm option\n";
			err = true;
		    }
		    break;
		default:
		    cerr << "invalid -C option\n";
		    err = true;
		    break;
		}
		break;
	    case 'L':
		getCorr = true;
		corrfn = &argv[i][2];
		break;
	    case 'U':
		getBound = true;
		sscanf(&argv[i][2], "%s", uboundfn);
		break;
	    default:
		err = true;
		break;
	    }
	}else{
	    getMesh = true;
	    sfn = argv[i];
	}
    }

    if (!getMesh || !getCmd) err=true;
    if (err) error("usage: edit_smesh [ -options ]");

    // read smesh
    ifstream s(sfn);
    if (!s){
	cerr << "edit_smesh::cannot open " << sfn << "\n";
	exit(1);
    }

    int nx, nz;
    double v_water, v_air;
    s >> nx >> nz >> v_water >> v_air;
    int nnodes = nx*nz;
    
    Array1d<double> xpos(nx), zpos(nz), topo(nx);
    Array2d<double> vgrid(nx,nz);
    for (int i=1; i<=nx; i++) s >> xpos(i);
    for (int i=1; i<=nx; i++) s >> topo(i);
    for (int k=1; k<=nz; k++) s >> zpos(k);
    for (int i=1; i<=nx; i++){
	for (int k=1; k<=nz; k++){
	    double vel;
	    s >> vel;
	    vgrid(i,k) = vel;
	}
    }

    // 2-D correlation length
    CorrelationLength2d *corr_p;
    if (getCorr){
	corr_p = new CorrelationLength2d(corrfn);
    }

    // upper bound
    Interface2d *ubound_p;
    Array1d<int> kstart(nx);
    kstart = 1;
    if (getBound || subMohoVel){
	ubound_p = new Interface2d(uboundfn);
	for (int i=1; i<=nx; i++){
	    double boundz = ubound_p->z(xpos(i));
	    for (int k=1; k<=nz; k++){
		if (zpos(k)+topo(i)>boundz) break;
		kstart(i) = k;
	    }
	    if (kstart(i)>1) kstart(i)++;
	}
    }
    
    // do something
    if (setAverage){
	for (int k=1; k<=nz; k++){
	    double ave=0;
	    for (int i=1; i<=nx; i++){
		ave += vgrid(i,k);
	    }
	    ave /= nx;
	    for (int i=1; i<=nx; i++){
		vgrid(i,k) = ave;
	    }
	}
    }else if (pasteSmesh){
	SlownessMesh2d smesh1(pfn);
	for (int i=1; i<=nx; i++){
	    double x = xpos(i);
	    double t = topo(i);
	    for (int k=1; k<=nz; k++){
		double z = t+zpos(k);
		Point2d p(x,z);
		if (!smesh1.inWater(p) && !smesh1.inAir(p)){
		    vgrid(i,k) = 1.0/smesh1.at(p); // override with smesh1
		}
	    }
	}
    }else if (paste1dprof){
        Interface2d prof1d(prof1dfn); // this usage of interface1d class
	                              // is quite ad hoc (caution)
	for (int i=1; i<=nx; i++){
	    double z0 = zpos(kstart(i));
	    for (int k=kstart(i); k<=nz; k++){
	        double z = zpos(k)-z0;
		double vval = prof1d.z(z);
		vgrid(i,k) = vval; 
            }
        }
    }else if (smooth){
	double Lh2 = Lh*Lh;
	double Lv2 = Lv*Lv;
	Array2d<double> new_vgrid(nx,nz);
	new_vgrid = vgrid;
	for (int i=1; i<=nx; i++){
	    double x = xpos(i);
	    double t = topo(i);
	    for (int k=kstart(i); k<=nz; k++){
		double z = t+zpos(k);

		if (getCorr){
		    Point2d p(x,z);
		    corr_p->at(p,Lh,Lv);
		    Lh2 = Lh*Lh;
		    Lv2 = Lv*Lv;
		}
		
		double sum = 0.0, beta_sum = 0.0;
		for (int ii=1; ii<=nx; ii++){
		    double xx = xpos(ii);
		    double tt = topo(ii);
		    double dx = x-xx;
		    if (abs(dx)<=Lh){
			double xexp = exp(-dx*dx/Lh2);
			for (int kk=kstart(ii); kk<=nz; kk++){
			    double zz = tt+zpos(kk);
			    double dz = z-zz;
			    if (abs(dz)<=Lv){
				double beta = xexp*exp(-dz*dz/Lv2);
				sum += beta*vgrid(ii,kk);
				beta_sum += beta;
			    }
			}
		    }
		}
		new_vgrid(i,k) = sum/beta_sum;
	    }
	}
	vgrid = new_vgrid;	
    }else if (checkerboard){
	double pi = 3.1412;
	double coeffx = 2.0*pi/ch;
	double coeffz = 2.0*pi/cv;
	for (int i=1; i<=nx; i++){
	    double x = xpos(i);
	    double t = topo(i);
	    for (int k=1; k<=nz; k++){
		double z = t+zpos(k);

		double origvel = vgrid(i,k);
		double newvel = origvel*(1.0+0.01*A*sin(coeffx*x)*sin(coeffz*z));
		vgrid(i,k) = newvel;
	    }
	}
    }else if (rectangular){
	for (int i=1; i<=nx; i++){
	    double x = xpos(i);
	    double t = topo(i);
	    for (int k=1; k<=nz; k++){
		double z = t+zpos(k);
		if (x>=dxmin && x<=dxmax && z>=dzmin && z<=dzmax){
		    double origvel = vgrid(i,k);
		    double newvel = origvel*(1.0+0.01*A);
		    vgrid(i,k) = newvel;
		}
	    }
	}
    }else if (gaussian){
	double rLh2 = 1.0/(Lh*Lh);
	double rLv2 = 1.0/(Lv*Lv);
	A *= 0.01; // percent to fraction
	for (int i=1; i<=nx; i++){
	    double x = xpos(i);
	    double dx = x-x0;
	    double dx2 = dx*dx;
	    double t = topo(i);
	    for (int k=1; k<=nz; k++){
		double z = t+zpos(k);
		double dz = z-z0;
		double dz2 = dz*dz;
		double coeff = exp(-dx2*rLh2-dz2*rLv2);
		vgrid(i,k) *= (1.0+A*coeff);
	    }
	}
    }else if (removeLVZ){
	for (int i=1; i<=nx; i++){
	    double vmax = 0.0;
	    for (int k=1; k<=nz; k++){
		double origvel = vgrid(i,k);
		if (origvel > vmax) vmax = origvel;
		if (origvel < vmax) vgrid(i,k) = vmax;
	    }
	}
    }else if (randomize){
	double coeff = 1.0/RAND_MAX;
	A *= 0.01; // percent to fraction
	Array2d<double> pur(nx,nz), new_pur(nx,nz);
	
	srand(seed);
	int ik=0;
	for (int i=1; i<=nx; i++){
	    double x = xpos(i);
	    double t = topo(i);
	    for (int k=kstart(i); k<=nz; k++){
		double z = t+zpos(k);
		double amp=0;
		if (ik%nrand==0)
		    amp = A*2.0*(coeff*rand()-0.5);;
//		if (getCorr){
//		    Point2d p(x,z);
//		    double h,v;
//		    corr_p->at(p,h,v);
//		    amp *= h*v;
//		}
		pur(i,k) = amp;
		ik++;
	    }
	}
	if (getCorr){
	    for (int i=1; i<=nx; i++){
		double x = xpos(i);
		double t = topo(i);
		for (int k=kstart(i); k<=nz; k++){
		    double z = t+zpos(k);
		    Point2d p(x,z);
		    corr_p->at(p,Lh,Lv);
		    double Lh2 = Lh*Lh;
		    double Lv2 = Lv*Lv;
		
		    double sum = 0.0, beta_sum = 0.0;
		    for (int ii=1; ii<=nx; ii++){
			double xx = xpos(ii);
			double tt = topo(ii);
			double dx = x-xx;
			if (abs(dx)<=Lh){
			    double xexp = exp(-dx*dx/Lh2);
			    for (int kk=kstart(ii); kk<=nz; kk++){
				double zz = tt+zpos(kk);
				double dz = z-zz;
				if (abs(dz)<=Lv){
				    double beta = xexp*exp(-dz*dz/Lv2);
				    sum += beta*pur(ii,kk);
				    beta_sum += beta;
				}
			    }
			}
		    }
		    new_pur(i,k) = sum/beta_sum;
		}
	    }
	    pur = new_pur;
	}
	for (int i=1; i<=nx; i++){
	    for (int k=kstart(i); k<=nz; k++){
		vgrid(i,k) *= (1.0+pur(i,k));
	    }
	}
    }else if (randomize2){
	double coeff = 1.0/RAND_MAX;
	A *= 0.01*2; // percent to fraction
	A = sqrt(A); 
	srand(seed);

	double pi2 = 3.1415926*2.0;
	Array1d<double> xamp, zamp, xph, zph, kx, kz;
	for (double x=rand_xmin; x<=rand_xmax; x+=rand_dx){
	    xamp.push_back(A*(coeff*rand()-0.5));
	    kx.push_back(pi2/x);
	    xph.push_back(pi2*(coeff*rand()-0.5));
//cerr << xamp.back() << ' ' << xph.back() << '\n';
        }
	for (double z=rand_zmin; z<=rand_zmax; z+=rand_dz){
	    zamp.push_back(A*(coeff*rand()-0.5));
	    kz.push_back(pi2/z);
	    zph.push_back(pi2*(coeff*rand()-0.5));
//cerr << zamp.back() << ' ' << zph.back() << '\n';

        }

	for (int i=1; i<=nx; i++){
	    double x = xpos(i);
            double t = topo(i);
	    for (int k=kstart(i); k<=nz; k++){
	        double z = t+zpos(k);
		double purx=0.0, purz=0.0;
		for (int ii=1; ii<=xamp.size(); ii++){
		    purx += xamp(ii)*sin(kx(ii)*x+xph(ii));
	        }
		for (int kk=1; kk<=zamp.size(); kk++){
                    purz += zamp(kk)*sin(kz(kk)*z+zph(kk));
	        }
		
		vgrid(i,k) *= (1.0+purx*purz);
	    }
         }
      }else if (randomize3){
	double coeff = 1.0/RAND_MAX;
	A *= 0.01*2; // percent to fraction
	srand(seed);

	for (int nn=1; nn<=N; nn++){
	  double Lh = rand_xmin+coeff*rand()*(rand_xmax-rand_xmin);
	  double Lv = rand_zmin+coeff*rand()*(rand_zmax-rand_zmin);
	  double rLh2 = 1.0/(Lh*Lh);
	  double rLv2 = 1.0/(Lv*Lv);
	  double amp = A*(coeff*rand()-0.5);
	  int ci = int(nx*coeff*rand()); if (ci==0) ci=1;
	  int ck = int(nz*coeff*rand()); if (ck==0) ck=1;

	  for (int i=1; i<=nx; i++){
	    double x = xpos(i);
	    double dx = x-xpos(ci);
	    double dx2 = dx*dx;
            double t = topo(i);
	    for (int k=kstart(i); k<=nz; k++){
	        double z = t+zpos(k);
		double dz = z-zpos(ck);
		double dz2 = dz*dz;
		double gauval = amp*exp(-dx2*rLh2-dz2*rLv2);
		vgrid(i,k) *= (1.0+gauval);
	    }
	  }
	}
    }else if (refineMesh){
	int rnx = (nx-1)*mx+1;
	int rnz = (nz-1)*mz+1;
	double dr = 1.0/mx;
	double ds = 1.0/mz;
	Array1d<double> rxpos(rnx), rzpos(rnz), rtopo(rnx);
	Array2d<double> rvgrid(rnx,rnz);

	for (int i=1; i<nx; i++){
	    int ri=i+(i-1)*(mx-1);
	    double x1=xpos(i), x2=xpos(i+1);
	    double t1=topo(i), t2=topo(i+1);
	    for (int ii=0; ii<=mx; ii++){
		double r=ii*dr;
		rxpos(ri+ii) = (1-r)*x1+r*x2;
		rtopo(ri+ii) = (1-r)*t1+r*t2;
	    }
	}
	for (int k=1; k<nz; k++){
	    int rk=k+(k-1)*(mz-1);
	    double z1=zpos(k), z2=zpos(k+1);
	    for (int kk=0; kk<=mz; kk++){
		double s=kk*ds;
		rzpos(rk+kk) = (1-s)*z1+s*z2;
	    }
	}
	rvgrid = -1;
	for (int i=1; i<nx; i++){
	    for (int k=1; k<nz; k++){
		int ri=i+(i-1)*(mx-1);
		int rk=k+(k-1)*(mz-1);
		double v1=vgrid(i,k), v2=vgrid(i,k+1);
		double v3=vgrid(i+1,k+1), v4=vgrid(i+1,k);
		for (int ii=0; ii<=mx; ii++){
		    double r=ii*dr;
		    for (int kk=0; kk<=mz; kk++){
			double s=kk*ds;
			if (rvgrid(ri+ii,rk+kk)<0){
			    rvgrid(ri+ii,rk+kk)
				= v1*(1-r)*(1-s)+v2*(1-r)*s
				+v3*r*s+v4*r*(1-s);
			}
		    }
		}
	    }
	}
	// output smesh
	cout << rnx << " " << rnz << " "
	     << v_water << " " << v_air << '\n';
	for (int i=1; i<=rnx; i++) cout << rxpos(i) << " ";
	cout << '\n';
	for (int i=1; i<=rnx; i++) cout << rtopo(i) << " ";
	cout << '\n';
	for (int k=1; k<=rnz; k++) cout << rzpos(k) << " ";
	cout << '\n';
	for (int i=1; i<=rnx; i++){
	    for (int k=1; k<=rnz; k++){
		cout << rvgrid(i,k) << " ";
	    }
	    cout << '\n';
	}
	return 0; // end the program here.
    }else if (subMohoVel){
	for (int i=1; i<=nx; i++){
	    for (int k=kstart(i); k<=nz; k++){
		vgrid(i,k) = vsubmoho;
	    }
	}
    }

    // output smesh
    cout << nx << " " << nz << " "
	 << v_water << " " << v_air << '\n';
    for (int i=1; i<=nx; i++) cout << xpos(i) << " ";
    cout << '\n';
    for (int i=1; i<=nx; i++) cout << topo(i) << " ";
    cout << '\n';
    for (int k=1; k<=nz; k++) cout << zpos(k) << " ";
    cout << '\n';
    for (int i=1; i<=nx; i++){
	for (int k=1; k<=nz; k++){
	    cout << vgrid(i,k) << " ";
	}
	cout << '\n';
    }
}
