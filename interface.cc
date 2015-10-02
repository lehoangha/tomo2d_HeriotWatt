/*
 * interface.cc
 *
 * Jun Korenaga, MIT/WHOI
 * January 1999
 */

#include "interface.h"
#include <iostream>
#include <fstream>
#include <error.h>
#include <util.h>

const double Interface2d::eps = 1e-6;

Interface2d::Interface2d(const SlownessMesh2d& smesh)
{
    int n=smesh.xpos.size();
    xpos.resize(n); zpos.resize(n);
    xpos = smesh.xpos;
    zpos = smesh.topo;
    calc_slope();
}

Interface2d::Interface2d(const char* fn)
{
    int n = countLines(fn);
    ifstream in(fn);
    xpos.resize(n); zpos.resize(n);
    for (int i=1; i<=n; i++){
	in >> xpos(i) >> zpos(i);
    }
    // sanity check
    for (int i=2; i<=n; i++){
	if (xpos(i)<=xpos(i-1)){
	    cerr << "Interface2d: illegal ordering of x nodes at i=" << i << '\n';
	    exit(1);
	}
    }
    
    calc_slope();
}

void Interface2d::calc_slope()
{
    slope.resize(xpos.size()-1);
    for (int i=1; i<=slope.size(); i++){
	slope(i) = (zpos(i+1)-zpos(i))/(xpos(i+1)-xpos(i));
    }
}

double Interface2d::xmin() const { return xpos.front(); }
double Interface2d::xmax() const { return xpos.back(); }

double Interface2d::z(double x) const
{
    if (x<=xpos.front()) return zpos.front();
    if (x>=xpos.back()) return zpos.back();
    
    int nx=xpos.size();
    for (int i=1; i<nx; i++){
	double xnode=xpos(i);
	double xnode2=xpos(i+1);

	if (abs(x-xnode)<eps){ // on a grid
	    return zpos(i);
	}else if (abs(x-xnode2)<eps){ // on a grid
	    return zpos(i+1); 
	}else if (x > xnode && x < xnode2){
	    double r = (x-xnode)/(xnode2-xnode);
	    double dz = zpos(i+1)-zpos(i);
	    return zpos(i)+r*dz;
	}
    }

    error("Interface2d::z - impossible!");
}

double Interface2d::dzdx(double x) const
{
    if (x<=xpos.front() || x>=xpos.back()) return 0.0; // out of bounds

    int nx=xpos.size();
    for (int i=1; i<nx; i++){
	double xnode=xpos(i);
	double xnode2=xpos(i+1);

	if (abs(x-xnode)<eps){ // on a grid
	    return 0.5*(slope(i)+slope(i+1));
	}else if (x > xnode && x < xnode2){
	    return slope(i);
	}
    }

    cerr << "at x=" << x << " nx=" << nx << '\n';
    error("Interface2d::dzdx - impossible!");
}

void Interface2d::locateInSegment(double x, int& j1, int& j2) const
{
    int nx=xpos.size();

    if (x<xpos.front()){ j1=j2=1; return; }
    if (x>xpos.back()){ j1=j2=nx; return; } // out of bounds

    for (int i=1; i<nx; i++){
	double xnode=xpos(i);
	double xnode2=xpos(i+1);
	if (x >= xnode && x < xnode2){
	    j1 = i;
	    j2 = i+1;
	    return;
	}
    }
    j1 = nx-1;
    j2 = nx;
    return;
}

double Interface2d::x(int i) const
{
    if (i>0 && i<=xpos.size()) return xpos(i);
    error("Interface2d::x - subscript out of range");
}

int Interface2d::numNodes() const { return xpos.size(); }

void Interface2d::set(const Array1d<double>& a)
{
    if (a.size() != zpos.size())
	error("Interface2d::set - size mismatch");
    for (int i=1; i<=a.size(); i++){
	zpos(i) = a(i);
    }
    calc_slope();
}

void Interface2d::get(Array1d<double>& a) const
{
    if (a.size() != zpos.size())
	error("Interface2d::get - size mismatch");
    for (int i=1; i<=a.size(); i++){
	a(i) = zpos(i);
    }
}

ostream&
operator<<(ostream& os, const Interface2d& itf)
{
    for (int i=1; i<=itf.xpos.size(); i++){
	os << itf.xpos(i) << " " <<  itf.zpos(i) << '\n';
    }
    return os;
}
    
