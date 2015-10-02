/*
 * corrlen.cc - correlation length functions' implementation
 * 
 * Jun Korenaga, MIT/WHOI
 * February 1999
 */

#include <iostream>
#include <fstream>
#include <util.h>
#include "corrlen.h"

CorrelationLength1d::CorrelationLength1d(const char* fn)
    : eps(1e-6)
{
    int n = countLines(fn);
    ifstream in(fn);
    xpos.resize(n); val.resize(n);
    for (int i=1; i<=n; i++){
	in >> xpos(i) >> val(i);
    }
}

double CorrelationLength1d::at(double x) const
{
    if (x<=xpos.front()) return val.front();
    if (x>=xpos.back()) return val.back();
    
    int nx=xpos.size();
    for (int i=1; i<nx; i++){
	double xnode=xpos(i);
	double xnode2=xpos(i+1);

	if (abs(x-xnode)<eps){ // on a grid
	    return val(i);
	}else if (abs(x-xnode2)<eps){ // on a grid
	    return val(i+1); 
	}else if (x > xnode && x < xnode2){
	    double r = (x-xnode)/(xnode2-xnode);
	    double dv = val(i+1)-val(i);
	    return val(i)+r*dv;
	}
    }

    error("CorrelationLengh1d::at() - impossible!");
}

CorrelationLength2d::CorrelationLength2d(const char* fn)
{
    ifstream s(fn);
    if (!s){
	cerr << "CorrelationLength2d::cannot open " << fn << "\n";
	exit(1);
    }

    s >> nx >> nz;
    xpos.resize(nx);
    topo.resize(nx);
    zpos.resize(nz);
    hgrid.resize(nx,nz);
    vgrid.resize(nx,nz);

    for (int i=1; i<=nx; i++) s >> xpos(i);
    for (int i=1; i<=nx; i++) s >> topo(i);
    for (int k=1; k<=nz; k++) s >> zpos(k);
    for (int i=1; i<=nx; i++){
	for (int k=1; k<=nz; k++){
	    s >> hgrid(i,k);
	}
    }
    for (int i=1; i<=nx; i++){
	for (int k=1; k<=nz; k++){
	    s >> vgrid(i,k);
	}
    }

    rdx_vec.resize(nx-1);
    rdz_vec.resize(nz-1);
    b_vec.resize(nx-1);
    for (int i=1; i<nx; i++){
	rdx_vec(i) = 1.0/(xpos(i+1)-xpos(i));
	b_vec(i) = topo(i+1)-topo(i);
    }
    for (int k=1; k<nz; k++){
	rdz_vec(k) = 1.0/(zpos(k+1)-zpos(k));
    }
}

void CorrelationLength2d::upperleft(const Point2d& p, Index2d& index) const
{
    // note: this code guarantees that the final guess index is bounded by
    //       valid range (1...nx-1)(1...nz-1).
    int i_guess=1, k_guess=1;

    if (xpos(i_guess)<=p.x()){
	if (i_guess < nx) i_guess++;
	while (xpos(i_guess)<=p.x() && i_guess<nx) i_guess++;
	i_guess--;
    }else{
	while (xpos(i_guess)>p.x() && i_guess>1) i_guess--;
    }
    double r = (p.x()-xpos(i_guess))*rdx_vec(i_guess);
    double dz = r*b_vec(i_guess)+topo(i_guess);

    if (zpos(k_guess)+dz<=p.y()){
	if (k_guess < nz) k_guess++;
	while (zpos(k_guess)+dz<=p.y() && k_guess<nz) k_guess++;
	k_guess--;
    }else{
	while (zpos(k_guess)+dz>p.y() && k_guess>1) k_guess--;
    }

    index.set(i_guess,k_guess);
}

void CorrelationLength2d::calc_local(const Point2d& pos, int i, int k,
				     double& r, double& s, double& rr, double& ss) const
{
    double rDX = rdx_vec(i);
    double rDZ = rdz_vec(k);
    double B = b_vec(i);

    // global coordinates with an origin at the upperleft node
    double x = pos.x()-xpos(i);
    double z = pos.y()-(zpos(k)+topo(i));
	
    // local coordinates
    r = x*rDX;
    s = (z-B*r)*rDZ;
    rr=1-r;
    ss=1-s;
}

void CorrelationLength2d::at(const Point2d& p, double& h, double& v) const
{
    Index2d index;
    upperleft(p,index);
    int i=index.i(), k=index.k();
    double r,s,rr,ss;
    calc_local(p,i,k,r,s,rr,ss);

    h=rr*ss*hgrid(i,k)+rr*s*hgrid(i,k+1)
	+r*s*hgrid(i+1,k+1)+r*ss*hgrid(i+1,k);
    v=rr*ss*vgrid(i,k)+rr*s*vgrid(i,k+1)
	+r*s*vgrid(i+1,k+1)+r*ss*vgrid(i+1,k);
}

DampingWeight2d::DampingWeight2d(const char* fn)
{
    ifstream s(fn);
    if (!s){
	cerr << "DampingWeight2d::cannot open " << fn << "\n";
	exit(1);
    }

    s >> nx >> nz;
    xpos.resize(nx);
    topo.resize(nx);
    zpos.resize(nz);
    wgrid.resize(nx,nz);

    for (int i=1; i<=nx; i++) s >> xpos(i);
    for (int i=1; i<=nx; i++) s >> topo(i);
    for (int k=1; k<=nz; k++) s >> zpos(k);
    for (int i=1; i<=nx; i++){
	for (int k=1; k<=nz; k++){
	    s >> wgrid(i,k);
	}
    }

    rdx_vec.resize(nx-1);
    rdz_vec.resize(nz-1);
    b_vec.resize(nx-1);
    for (int i=1; i<nx; i++){
	rdx_vec(i) = 1.0/(xpos(i+1)-xpos(i));
	b_vec(i) = topo(i+1)-topo(i);
    }
    for (int k=1; k<nz; k++){
	rdz_vec(k) = 1.0/(zpos(k+1)-zpos(k));
    }
}

void DampingWeight2d::upperleft(const Point2d& p, Index2d& index) const
{
    // note: this code guarantees that the final guess index is bounded by
    //       valid range (1...nx-1)(1...nz-1).
    int i_guess=1, k_guess=1;

    if (xpos(i_guess)<=p.x()){
	if (i_guess < nx) i_guess++;
	while (xpos(i_guess)<=p.x() && i_guess<nx) i_guess++;
	i_guess--;
    }else{
	while (xpos(i_guess)>p.x() && i_guess>1) i_guess--;
    }
    double r = (p.x()-xpos(i_guess))*rdx_vec(i_guess);
    double dz = r*b_vec(i_guess)+topo(i_guess);

    if (zpos(k_guess)+dz<=p.y()){
	if (k_guess < nz) k_guess++;
	while (zpos(k_guess)+dz<=p.y() && k_guess<nz) k_guess++;
	k_guess--;
    }else{
	while (zpos(k_guess)+dz>p.y() && k_guess>1) k_guess--;
    }

    index.set(i_guess,k_guess);
}

void DampingWeight2d::calc_local(const Point2d& pos, int i, int k,
				 double& r, double& s, double& rr, double& ss) const
{
    double rDX = rdx_vec(i);
    double rDZ = rdz_vec(k);
    double B = b_vec(i);

    // global coordinates with an origin at the upperleft node
    double x = pos.x()-xpos(i);
    double z = pos.y()-(zpos(k)+topo(i));
	
    // local coordinates
    r = x*rDX;
    s = (z-B*r)*rDZ;
    rr=1-r;
    ss=1-s;
}

void DampingWeight2d::at(const Point2d& p, double& dw) const
{
    Index2d index;
    upperleft(p,index);
    int i=index.i(), k=index.k();
    double r,s,rr,ss;
    calc_local(p,i,k,r,s,rr,ss);

    dw=rr*ss*wgrid(i,k)+rr*s*wgrid(i,k+1)
	+r*s*wgrid(i+1,k+1)+r*ss*wgrid(i+1,k);
}

