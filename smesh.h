/*
 * smesh.h - slowness mesh interface
 *
 * Jun Korenaga, MIT/WHOI
 * December 1998
 */

#ifndef _TOMO_SMESH_H_
#define _TOMO_SMESH_H_

#include <list>
#include <array.h> // from mconv source
#include <geom.h>
#include "heap_deque.h"
#include "index.h"

class Interface2d; // forward declaration

class SlownessMesh2d {
public:
    SlownessMesh2d(const char*); // file input

    void set(const Array1d<double>&);
    void get(Array1d<double>&) const;
    void vget(Array1d<double>&) const;
    
    // general bookkeeping functions
    int numNodes() const { return nnodes; }
    int Nx() const { return nx; }
    int Nz() const { return nz; }
    double xmin() const;
    double xmax() const;
    double zmin() const;
    double zmax() const;
    const Index2d& nodeIndex(int i) const;
    int nodeIndex(int i, int k) const;
    Point2d nodePos(int i) const;

    // slowness interpolation
    double at(const Point2d& pos) const;
    double at(const Point2d& pos, Index2d& guess) const;
    double at(const Point2d& pos, Index2d& guess,
	      double& dudx, double& dudz) const;
    double atWater() const { return p_water; }

    // for graph traveltime calculation
    double calc_ttime(int node_src, int node_rcv) const; 
									   
    // cell-oriented functions
    int numCells() const { return ncells; }
    void cellNodes(int, int&, int&, int&, int&) const;
    int cellIndex(int i, int k) const { return index2cell(i,k); }
    void cellGradientKernel(int, Array2d<double>&, double, double) const;
    void cellNormKernel(int, Array2d<double>&) const;
    int locateInCell(const Point2d&, Index2d&,
		     int&, int&, int&, int&,
		     double&, double&, double&, double&) const;
    void nodalCellVolume(Array1d<double>&,
			 Array1d<double>&, Array1d<Point2d>&) const; // for gravity
    
    // miscellaneous
    int nearest(const Point2d& src) const;
    void nearest(const Interface2d& itf, Array1d<int>& inodes) const;
    bool inWater(const Point2d& pos) const;
    bool inAir(const Point2d& pos) const;

    // output functions
    void outMesh(ostream&) const;
    void printElements(ostream&) const;
    void printVGrid(ostream&,bool) const;
    void printVGrid(ostream&,
		    double,double,double,double,double,double) const;
    void printMaskGrid(ostream&, const Array1d<int>&) const;
    void printMaskGrid(ostream&, const Array1d<double>&) const;

    friend class Interface2d;
    
private:
    void upperleft(const Point2d& pos, Index2d& guess) const;
    void calc_local(const Point2d& pos, int, int,
		    double&, double&, double&, double&) const;
    void commonGradientKernel();
    void commonNormKernel();
    bool in_water(const Point2d& pos, const Index2d& guess) const;
    bool in_air(const Point2d& pos, const Index2d& guess) const;
//    double almost_exact_ttime(double v1, double v2,
//			      double dpath) const;

    int nx, nz, nnodes, ncells;
    double p_water; /* slowness of water column */
    double p_air; /* slowness of air column */
    Array2d<double> pgrid, vgrid;
    Array2d<int> ser_index, index2cell;
    Array1d<Index2d> node_index;
    Array1d<int> cell_index;
    Array1d<double> xpos, topo, zpos;
    Array1d<double> rdx_vec, rdz_vec, b_vec;
    Array1d<double> dx_vec, dz_vec;
    Array2d<double> Sm_H1, Sm_H2, Sm_V, T_common;
    const double eps;
};

#endif /* _TOMO_SMESH_H_ */
