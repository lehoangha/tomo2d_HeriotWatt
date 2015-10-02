/*
 * graph.h - graph method related classes
 *
 * Jun Korenaga, MIT/WHOI
 * January 1999
 */

#ifndef _TOMO_GRAPH_H_
#define _TOMO_GRAPH_H_

#include <list>
#include <array.h> // from mconv source
#include <geom.h>
#include "heap_deque.h"
#include "index.h"
#include "interface.h"

class TTCmp {
public:
    TTCmp(Array1d<double>& a) : data(a) {}
    bool operator()(int i, int j) const { return data(i)>data(j); }
private:
    Array1d<double>& data;
};

class GraphSolver2d {
public:
    GraphSolver2d(const SlownessMesh2d&, int xorder, int zorder);

    int xOrder() const { return fs_xorder; }
    int zOrder() const { return fs_zorder; }
    void solve(const Point2d& src);
    void solve_refl(const Point2d& src, const Interface2d& itf);
    void do_refl_downward(); 
    void limitRange(double xmin, double xmax);
    void delimitRange();
    void refineIfLong(double);
    double critLength() const { return crit_len; }
    
    void pickPath(const Point2d& rcv, Array1d<Point2d>& path) const;
    void pickPathThruWater(const Point2d& rcv, Array1d<Point2d>& path,
			   int&, int&) const;
    void pickReflPath(const Point2d& rcv, Array1d<Point2d>& path,
		      int&, int&) const;
    void pickReflPathThruWater(const Point2d& rcv, Array1d<Point2d>& path,
			       int&, int&, int&, int&) const;
    
    void printPath(ostream&) const;

private:
    void findRayEntryPoint(const Point2d& rcv, int& cur, 
			   const Array1d<double>& nodetime) const;
    
    const SlownessMesh2d& smesh;
    int nnodes;
    bool is_solved, is_refl_solved, solve_refl_downward;
    Point2d src;
    int fs_xorder, fs_zorder;
    Array1d<int> prev_node, prev_node_down, prev_node_up;
    Array1d<double> ttime, ttime_down, ttime_up;
    list<int> C;
    TTCmp cmp, cmp_down, cmp_up;
    heap_deque<int,TTCmp> B, B_down, B_up;
    const Interface2d* itfp;
    typedef list<int>::const_iterator nodeBrowser;
    typedef list<int>::iterator nodeIterator;
    typedef heap_deque<int,TTCmp>::const_iterator heapBrowser;

    const double local_search_radius;
    double xmin, xmax;
    bool refine_if_long;
    double crit_len;
};

class ForwardStar2d {
public:
    ForwardStar2d(const Index2d& id, int ix, int iz);
    bool isIn(const Index2d&) const;

private:
    const Index2d& orig;
    int xorder, zorder;
};

#endif /* _TOMO_GRAPH_H_ */
