/*
 * betaspline.cc - beta spline implementation
 *
 * Jun Korenaga, MIT/WHOI
 * January 1999
 */

#include "betaspline.h"
#include <error.h> // from mconv

BetaSpline2d::BetaSpline2d(double beta1, double beta2, int n)
    : nintp(n)
{
    if (nintp<=2) error("BetaSpline2d::invalid nintp");
    calc_c(beta1,beta2);
    calc_u();
    calc_B();
    calc_dBdu();
}

void BetaSpline2d::interpolate(const Point2d& p1, const Point2d& p2,
			       const Point2d& p3, const Point2d& p4,
			       Array1d<Point2d>& Q) const
{
    if (Q.size() != nintp) error("BetaSpline2d::size mismatch");

    for (int i=1; i<=nintp; i++){
	Q(i) = b[0](i)*p1+b[1](i)*p2+b[2](i)*p3+b[3](i)*p4;
    }
}

void BetaSpline2d::interpolate(const Point2d& p1, const Point2d& p2,
			       const Point2d& p3, const Point2d& p4,
			       Array1d<Point2d>& Q, Array1d<Point2d>& dQdu) const
{
    if (Q.size() != nintp) error("BetaSpline2d::size mismatch");

    for (int i=1; i<=nintp; i++){
	Q(i) = b[0](i)*p1+b[1](i)*p2+b[2](i)*p3+b[3](i)*p4;
	dQdu(i) = dbdu[0](i)*p1+dbdu[1](i)*p2+dbdu[2](i)*p3+dbdu[3](i)*p4;
    }
}

void BetaSpline2d::resetNIntp(int n)
{
    if (n<=2) error("BetaSpline2d::invalid nintp");
    nintp = n;
    calc_u();
    calc_B();
    calc_dBdu();
}

void BetaSpline2d::calc_c(double beta1, double beta2)
{
    double beta1_2 = beta1*beta1;
    double beta1_3 = beta1_2*beta1;
    double delta = 2.0*beta1_3+4.0*(beta1_2+beta1)+beta2+2.0;
    double rdelta = 1.0/delta;

    c[0][0] = 2.0*beta1_3*rdelta;
    c[0][1] = (4.0*(beta1_2+beta1)+beta2)*rdelta;
    c[0][2] = 2.0*rdelta;
    c[0][3] = 0.0;

    c[1][0] = -6.0*beta1_3*rdelta;
    c[1][1] = 6.0*beta1*(beta1_2-1.0)*rdelta;
    c[1][2] = 6.0*beta1*rdelta;
    c[1][3] = 0.0;

    c[2][0] = -1.0*c[1][0];
    c[2][1] = 3.0*(-2.0*beta1_3-2.0*beta1_2-beta2)*rdelta;
    c[2][2] = 3.0*(2.0*beta1_2+beta2)*rdelta;
    c[2][3] = 0.0;

    c[3][0] = c[1][0]/3.0;
    c[3][1] = 2.0*(beta1_3+beta1_2+beta1+beta2)*rdelta;
    c[3][2] = -2.0*(beta1_2+beta1+beta2+1.0)*rdelta;
    c[3][3] = c[0][2];
}

void BetaSpline2d::calc_u()
{
    u.resize(nintp);

    // set local coordinate
    double du = 1.0/(nintp-1);
    for (int i=1; i<=nintp; i++){
	u(i) = (i-1)*du;
    }
}

void BetaSpline2d::calc_B()
{
    for (int k=0; k<4; k++) b[k].resize(nintp,0.0);

    for (int i=1; i<=nintp; i++){
	double ug=1.0;
	for (int j=0; j<4; j++){
	    for (int k=0; k<4; k++){
		b[k](i) += c[j][k]*ug;
	    }
	    ug *= u(i);
	}
    }
}    

void BetaSpline2d::calc_dBdu()
{
    for (int k=0; k<4; k++) dbdu[k].resize(nintp,0.0);

    for (int i=1; i<=nintp; i++){
	double ug=1.0;
	for (int j=1; j<=3; j++){
	    for (int k=0; k<4; k++){
		dbdu[k](i) += c[j][k]*j*ug;
	    }
	    ug *= u(i);
	}
    }
}    

// helper functions
void makeBSpoints(const Array1d<Point2d>& orig, Array1d<const Point2d*>& pp)
{
    int np=orig.size();
    pp.resize(np+4);

    pp(1) = pp(2) = &orig(1);
    for (int i=1; i<=np; i++){
	pp(i+2) = &orig(i);
    }
    pp(np+3) = pp(np+4) = &orig(np);
}

void makeBSpoints(const list<Point2d>& orig, Array1d<const Point2d*>& pp)
{
    int np=orig.size();
    pp.resize(np+4);

    list<Point2d>::const_iterator pt=orig.begin();
    pp(1) = pp(2) = &(*pt);
    int i=3;
    while(pt != orig.end()){
	pp(i) = &(*pt);
	i++;
	pt++; 
    }
    pp(np+3) = pp(np+4) = &(*pt);
}

void printCurve(ostream& os,
		const Array1d<Point2d>& orig, const BetaSpline2d& bs)
{
    Array1d<const Point2d*> pp;
    makeBSpoints(orig,pp);

    int np=orig.size();
    int nintp=bs.numIntp();
    Array1d<Point2d> Q(nintp);
    for (int i=1; i<=np+1; i++){
	int j1=i;
	int j2=i+1;
	int j3=i+2;
	int j4=i+3;
	bs.interpolate(*pp(j1),*pp(j2),*pp(j3),*pp(j4),Q);
	for (int j=1; j<=nintp; j++){
	    os << Q(j).x() << " " << Q(j).y() << '\n';
	}
    }
}

void printCurve(ostream& os,
		const list<Point2d>& orig, const BetaSpline2d& bs)
{
    Array1d<const Point2d*> pp;
    makeBSpoints(orig,pp);

    int np=orig.size();
    int nintp=bs.numIntp();
    Array1d<Point2d> Q(nintp);
    for (int i=1; i<=np+1; i++){
	int j1=i;
	int j2=i+1;
	int j3=i+2;
	int j4=i+3;
	bs.interpolate(*pp(j1),*pp(j2),*pp(j3),*pp(j4),Q);
	for (int j=1; j<=nintp; j++){
	    os << Q(j).x() << " " << Q(j).y() << '\n';
	}
    }
}


