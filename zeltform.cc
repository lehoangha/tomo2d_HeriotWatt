/*
 * zeltform.cc
 * 
 * Jun Korenaga, MIT/WHOI
 * January 1999
 */

#include "zeltform.h"
#include <iostream>
#include <fstream>
#include <list>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <util.h>
#include <error.h>

TrapezoidCell2d::TrapezoidCell2d(double _x1, double _x2,
				 double _z1, double _z2, double _z3, double _z4,
				 double _v1, double _v2, double _v3, double _v4)
    : x1(_x1), x2(_x2), v1(_v1), v2(_v2), v3(_v3), v4(_v4)
{
    double rdx = 1.0/(x1-x2);
    s1 = (_z1-_z2)*rdx;
    b1 = _z1-s1*x1;
    s2 = (_z3-_z4)*rdx;
    b2 = _z3-s2*x1;

    c1 = s2*(x2*v1-x1*v2)+b2*(v2-v1)-s1*(x2*v3-x1*v4)-b1*(v4-v3);
    c2 = s2*(v2-v1)-s1*(v4-v3);
    c3 = x1*v2-x2*v1+x2*v3-x1*v4;
    c4 = v1-v2+v4-v3;
    c5 = b2*(x2*v1-x1*v2)-b1*(x2*v3-x1*v4);
    c6 = (s2-s1)*(x2-x1);
    c7 = (b2-b1)*(x2-x1);

// check
//    double vv1 = at(Point2d(x1,_z1));
//    double vv2 = at(Point2d(x2,_z2));
//    double vv3 = at(Point2d(x1,_z3));
//    double vv4 = at(Point2d(x2,_z4));
//    cerr << "input: " << v1 << " " << v2 << " " << v3 << " " << v4 << '\n';
//    cerr << "out: " << vv1 << " " << vv2 << " " << vv3 << " " << vv4 << '\n';
//    cerr << '\n';
    
}

bool TrapezoidCell2d::isIn(const Point2d& p) const
{
    double x=p.x(), z=p.y();
    if (x<x1 || x>x2) return false;

    double zup = s1*x+b1;
    double zdown = s2*x+b2;
    if (z<zup || z>zdown) return false;

    return true;
}

double TrapezoidCell2d::at(const Point2d& p) const
{
    double x=p.x(), z=p.y();

    double val1 = (c1+c2*x)*x + (c3+c4*x)*z + c5;
    double val2 = c6*x+c7;

    return val1/val2;
}

void TrapezoidCell2d::dumpCell(ostream& os) const
{
    os << x1 << " " << x1*s1+b1 << '\n';
    os << x2 << " " << x2*s1+b1 << '\n';
    os << x2 << " " << x2*s2+b2 << '\n';
    os << x1 << " " << x1*s2+b2 << '\n';
}

ZeltVelocityModel2d::ZeltVelocityModel2d(char *fn)
{
    ifstream s(fn);
    if (!s){
	cerr << "ZeltVelocityModel2d::can't open " << fn << '\n';
	exit(1);
    }

    node_p.resize(3);
    node_p(1) = &depth_node;
    node_p(2) = &vupper_node;
    node_p(3) = &vlower_node;
    
    char line[MaxStr];
    Array1d<double> tmp_val;

    int itype = 1;
    int ilayer = 1;
    node_p(itype)->push_back(new ZNode2d);
    while (s.getline(line,MaxStr)){
	// read x coordinate line
	int nx=readLine(line,tmp_val);
	for (int i=2; i<=tmp_val.size(); i++){
	    (*node_p(itype))(ilayer)->x.push_back(tmp_val(i));
	}

	// read value line
	s.getline(line,MaxStr);
	int nv=readLine(line,tmp_val);
	if (nx!=nv) error("ZeltVelocityModel2d::invalid input file");
	for (int i=2; i<=tmp_val.size(); i++){
	    (*node_p(itype))(ilayer)->val.push_back(tmp_val(i));
	}
	int cont_flag = int(tmp_val(1));

	// read inversion flag line (do nothing for this line)
	if (!s.getline(line,MaxStr)) break; // end of v.in

	if (cont_flag == 0){
	    if (itype == 3){
		itype = 1; ilayer++;
	    }else{
		itype++;
	    }
	    node_p(itype)->push_back(new ZNode2d);
	}
    }
    int nlayer=ilayer;

    for (ilayer=1; ilayer<nlayer; ilayer++){
	list<double> xs;
	for (int itype=1; itype<=3; itype++){
	    const Array1d<double>* xp = &((*node_p(itype))(ilayer)->x);
	    for (int i=1; i<=xp->size(); i++){
		xs.push_back((*xp)(i));
	    }
	}
	const Array1d<double>* xp = &((*node_p(1))(ilayer+1)->x);
	for (int i=1; i<=xp->size(); i++){
	    xs.push_back((*xp)(i));
	}
	xs.sort();
	xs.unique();

	list<double>::const_iterator p1 = xs.begin();
	list<double>::const_iterator p2 = xs.begin(); p2++;
	while(p2!=xs.end()){
	    double x1 = *p1++; 
	    double x2 = *p2++; 
	    double z1 = interp((*node_p(1))(ilayer), x1);
	    double z2 = interp((*node_p(1))(ilayer), x2);
	    double z3 = interp((*node_p(1))(ilayer+1), x1);
	    double z4 = interp((*node_p(1))(ilayer+1), x2);
	    double v1 = interp((*node_p(2))(ilayer), x1);
	    double v2 = interp((*node_p(2))(ilayer), x2);
	    double v3 = interp((*node_p(3))(ilayer), x1);
	    double v4 = interp((*node_p(3))(ilayer), x2);
	    cells.push_back(new TrapezoidCell2d(x1,x2,z1,z2,z3,z4,
						v1,v2,v3,v4));
	}
    }
    
    // check the above (by printing out)
//     for (int i=1; i<=node_p(1)->size(); i++){
// 	cout << "ilayer=" << i << '\n';
// 	int itype_max;
// 	if (i<node_p(1)->size()){
// 	    itype_max = 3;
// 	}else{
// 	    itype_max = 1;
// 	}
// 	for (int itype=1; itype<=itype_max; itype++){
// 	    cout << "itype=" << itype << '\n';
// 	    Array1d<double>* xp = &((*node_p(itype))(i)->x);
// 	    Array1d<double>* vp = &((*node_p(itype))(i)->val);
// 	    for (int j=1; j<=xp->size(); j++){
// 		cout << "(" << (*xp)(j) << "," << (*vp)(j) << ") ";
// 	    }
// 	    cout << '\n';
// 	}
//     }
	
}
    
double ZeltVelocityModel2d::at(double x, double z) const
{
    Point2d p(x,z);
    for (int i=1; i<=cells.size(); i++){
	if (cells(i)->isIn(p)) return cells(i)->at(p);
    }

    cerr << "ZeltVelocityModel2d::at - failed at ("
	 << x << "," << z << ")\n";
    exit(1);
}

void ZeltVelocityModel2d::getTopo(int ilayer, double dx,
				  Array1d<double>& xs,
				  Array1d<double>& ts) const
{
    if (ilayer<1 || ilayer>node_p(1)->size())
	error("ZeltVelocityModel2d::getTopo - invalid layer number");
    
    const Array1d<double>* xp = &((*node_p(1))(ilayer)->x);
    const Array1d<double>* vp = &((*node_p(1))(ilayer)->val);

    xs.resize(0); ts.resize(0);
    for (int i=1; i<xp->size(); i++){
	double DX = (*xp)(i+1)-(*xp)(i);
	int ndiv = int(DX/dx+0.5);
	double ddx = DX/ndiv;
	for (int j=0; j<ndiv; j++){
	    double x = (*xp)(i)+j*ddx;
	    double t = interp((*node_p(1))(ilayer), x);
	    xs.push_back(x);
	    ts.push_back(t);
	}
    }
    xs.push_back(xp->back());
    ts.push_back(vp->back());
}

void ZeltVelocityModel2d::dumpNodes(const char* fn) const
{
    char dfn[MaxStr], vfn[MaxStr], cfn[MaxStr];
    sprintf(dfn, "%s.dnodes", fn);
    sprintf(vfn, "%s.vnodes", fn);
    sprintf(cfn, "%s.cells", fn);

    ofstream ds(dfn), vs(vfn), cs(cfn);

    for (int ilayer=1; ilayer<=node_p(1)->size(); ilayer++){
	const Array1d<double>* xp = &((*node_p(1))(ilayer)->x);
	const Array1d<double>* vp = &((*node_p(1))(ilayer)->val);

	for (int i=1; i<=xp->size(); i++){
	    double x = (*xp)(i);
	    double z = (*vp)(i);
	    ds << x << " " << z << '\n';
	}
    }
    for (int ilayer=1; ilayer<=node_p(2)->size(); ilayer++){
	const Array1d<double>* xp = &((*node_p(2))(ilayer)->x);

	for (int i=1; i<=xp->size(); i++){
	    double x = (*xp)(i);
	    double z = interp((*node_p(1))(ilayer),x);
	    vs << x << " " << z << '\n';
	}
    }
    for (int ilayer=1; ilayer<=node_p(3)->size(); ilayer++){
	const Array1d<double>* xp = &((*node_p(3))(ilayer)->x);

	for (int i=1; i<=xp->size(); i++){
	    double x = (*xp)(i);
	    double z = interp((*node_p(1))(ilayer),x);
	    vs << x << " " << z << '\n';
	}
    }
    for (int i=1; i<=cells.size(); i++){
	cs << ">\n";
	cells(i)->dumpCell(cs);
    }
}

int ZeltVelocityModel2d::readLine(char *line, Array1d<double>& tmp)
{
    tmp.resize(0);
    char *word = strtok(line, " ");
    while (word != NULL){
	tmp.push_back(atof(word));
	word = strtok(NULL, " ");
    }

    return tmp.size();
}

double ZeltVelocityModel2d::interp(const ZNode2d* p, double x) const
{
    const double eps = 1e-6;
    const int nx=p->x.size();

    if (nx==1){
	return p->val(1);
    }
    for (int i=1; i<nx; i++){
	double xnode=p->x(i);
	double xnode2=p->x(i+1);
	
	if (i==1 && (x-xnode)<=eps){ // to the left
	    return p->val.front();
	}else if (i==(nx-1) && (x-xnode2)>=eps){ // to the right
	    return p->val.back();
	}else if (abs(x-xnode)<eps){ // on a grid
	    return p->val(i);
	}else if (abs(x-xnode2)<eps){ // on a grid
	    return p->val(i+1); 
	}else if (x > xnode && x < xnode2){ // interpolate
	    double r = (x-xnode)/(xnode2-xnode);
	    double dv = (p->val(i+1))-(p->val(i));
	    return p->val(i)+r*dv;
	}
    }

    cerr << "x=" << x << " nx=" << nx << '\n';
    for (int i=1; i<=nx; i++){
	cerr << p->x(i) << ", ";
    }
    cerr << '\n';
    error("ZeltVelocityModel2d::interp - impossible!");
}
