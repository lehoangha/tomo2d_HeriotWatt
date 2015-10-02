/*
 * graph.cc - graph solver implementation
 *
 * Jun Korenaga, MIT/WHOI
 * January 1999
 */

#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include "smesh.h"
#include "graph.h"

GraphSolver2d::GraphSolver2d(const SlownessMesh2d& m, int xorder, int zorder)
    : smesh(m), nnodes(m.numNodes()),
      is_solved(false), is_refl_solved(false), solve_refl_downward(false),
      fs_xorder(xorder), fs_zorder(zorder),
      cmp(ttime), B(cmp),
      cmp_down(ttime_down), B_down(cmp_down),
      cmp_up(ttime_up), B_up(cmp_up),
      local_search_radius(10.0), xmin(smesh.xmin()), xmax(smesh.xmax()),
      refine_if_long(false)
{
    prev_node.resize(nnodes);
    ttime.resize(nnodes);
}

void GraphSolver2d::limitRange(double x1, double x2)
{
    if (x1>x2) error("GraphSolver2d::limitRange - invalid input");
    xmin=x1-local_search_radius;
    xmax=x2+local_search_radius;
}

void GraphSolver2d::delimitRange()
{
    xmin = smesh.xmin();
    xmax = smesh.xmax();
}

void GraphSolver2d::refineIfLong(double x)
{
    refine_if_long = true;
    crit_len = x;
}

void GraphSolver2d::solve(const Point2d& s)
{
    double ttime_inf = 1e30;

    src = s;
    if (smesh.inWater(src) || smesh.inAir(src)){
	error("GraphSolver2d::solve - currently unsupported source configuration detected.\n");
    }
    int Nsrc = smesh.nearest(src);
    
    // initialization
    B.resize(0);
    prev_node(Nsrc) = Nsrc;
    ttime(Nsrc) = 0.0;
    B.push(Nsrc); 
	
    C.resize(0);
    for (int i=1; i<=nnodes; i++){
	if (i!=Nsrc){
	    double x = smesh.nodePos(i).x();
	    if (x>=xmin && x<=xmax){
		C.push_back(i);
		ttime(i) = ttime_inf;
		prev_node(i) = Nsrc;
	    }
	}
    }

    while(B.size() != 0){
	// find a node with minimum traveltime from B,
	// (and transfer it to A)
	int N0=B.top();
	B.pop();

	// form a forward star
	ForwardStar2d fstar(smesh.nodeIndex(N0), fs_xorder, fs_zorder);

	// loop over [FS and B] to update traveltime
	heapBrowser pB=B.begin();
	while(pB!=B.end()){
	    if (fstar.isIn(smesh.nodeIndex(*pB))){
		double tmp_ttime=ttime(N0)+smesh.calc_ttime(N0,*pB);
		double orig_ttime=ttime(*pB);
 		if (orig_ttime>tmp_ttime){
		    ttime(*pB) = tmp_ttime;
		    prev_node(*pB) = N0;
		    B.promote(pB);
		}
	    }
	    pB++;  
	} 
	
	// loop over [FS and C] to calculate traveltime
	// and transfer them to B
	nodeIterator pC=C.begin();
	nodeIterator epC;
	while(pC!=C.end()){
	    if (fstar.isIn(smesh.nodeIndex(*pC))){
		ttime(*pC) = ttime(N0)+smesh.calc_ttime(N0,*pC);
		prev_node(*pC) = N0;
		B.push(*pC);

		epC = pC;
		pC++;
		C.erase(epC);
	    }else{ 
		pC++;
	    }
	}
    }
    is_solved = true;
}

void GraphSolver2d::do_refl_downward(){ solve_refl_downward = true; }

void GraphSolver2d::solve_refl(const Point2d& s, const Interface2d& itf)
{
    if (!is_solved && !solve_refl_downward)
	error("GraphSolver2d::not ready for solve_refl()");
    
    double ttime_inf = 1e30;

    if (prev_node_down.size() != nnodes) prev_node_down.resize(nnodes);
    if (ttime_down.size() != nnodes) ttime_down.resize(nnodes);
    if (prev_node_up.size() != nnodes) prev_node_up.resize(nnodes);
    if (ttime_up.size() != nnodes) ttime_up.resize(nnodes);
    itfp = &itf;
    Array1d<int> itfnodes;
    smesh.nearest(itf,itfnodes);

    if (solve_refl_downward){
	// initialization (downgoing)
	src = s;
	if (smesh.inWater(src) || smesh.inAir(src)){
	    error("GraphSolver2d::solve_refl - currently unsupported source configuration detected.\n");
	}
	int Nsrc = smesh.nearest(src);
    
	B_down.resize(0);
	prev_node_down(Nsrc) = Nsrc;
	ttime_down(Nsrc) = 0.0;
	B_down.push(Nsrc); 
	
	C.resize(0);
	for (int i=1; i<=nnodes; i++){
	    if (i!=Nsrc){
		double x = smesh.nodePos(i).x();
		double z = smesh.nodePos(i).y();
		double itfz = itf.z(x);
		if (x>=xmin && x<=xmax){
		    bool is_ok = false;
		    if (z<=itfz){
			is_ok = true;
		    }else{
			for (int n=1; n<=itfnodes.size(); n++){
			    if (itfnodes(n)==i){
				is_ok = true; break;
			    }
			}
		    }
		    if (is_ok){
			C.push_back(i);
			ttime_down(i) = ttime_inf;
			prev_node_down(i) = Nsrc;
		    }
		}
	    }
	}

	while(B_down.size() != 0){
	    // find a node with minimum traveltime from B,
	    // (and transfer it to A)
	    int N0=B_down.top();
	    B_down.pop();

	    // form a forward star
	    ForwardStar2d fstar(smesh.nodeIndex(N0), fs_xorder, fs_zorder);

	    // loop over [FS and B] to update traveltime
	    heapBrowser pB=B_down.begin();
	    while(pB!=B_down.end()){
		if (fstar.isIn(smesh.nodeIndex(*pB))){
		    double tmp_ttime=ttime_down(N0)+smesh.calc_ttime(N0,*pB);
		    double orig_ttime=ttime_down(*pB);
		    if (orig_ttime>tmp_ttime){
			ttime_down(*pB) = tmp_ttime;
			prev_node_down(*pB) = N0;
			B_down.promote(pB);
		    }
		}
		pB++;  
	    } 
	
	    // loop over [FS and C] to calculate traveltime
	    // and transfer them to B
	    nodeIterator pC=C.begin();
	    nodeIterator epC;
	    while(pC!=C.end()){
		if (fstar.isIn(smesh.nodeIndex(*pC))){
		    ttime_down(*pC) = ttime_down(N0)+smesh.calc_ttime(N0,*pC);
		    prev_node_down(*pC) = N0;
		    B_down.push(*pC);

		    epC = pC;
		    pC++;
		    C.erase(epC);
		}else{ 
		    pC++;
		}
	    }
	}
    }else{
	// use refraction's graph solution
	prev_node_down = prev_node;
	ttime_down = ttime;
    }
    
    // initialization (upgoing)
    B_up.resize(0);
    C.resize(0);
    for (int i=1; i<=nnodes; i++){
	bool is_itf = false;
	double x=smesh.nodePos(i).x();
	double z=smesh.nodePos(i).y();
	double itfz=itf.z(x);
	if (x>=xmin && x<=xmax){
	    for (int n=1; n<=itfnodes.size(); n++){
		if (itfnodes(n)==i){
		    is_itf = true;
		    prev_node_up(i) = -1;
		    ttime_up(i) = ttime_down(i);
		    B_up.push(i);
		    break;
		}
	    }
	    if (!is_itf && z<=itfz){ // take nodes above the reflector
		ttime_up(i) = ttime_inf;
		prev_node_up(i) = 0;
		C.push_back(i);
	    }
	}
    }

    while(B_up.size() != 0){
	// find a node with minimum traveltime from B,
	// (and transfer it to A)
	int N0=B_up.top();
	B_up.pop();

	// form a forward star
	ForwardStar2d fstar(smesh.nodeIndex(N0), fs_xorder, fs_zorder);

	// loop over [FS and B] to update traveltime
	heapBrowser pB=B_up.begin();
	while(pB!=B_up.end()){
	    if (fstar.isIn(smesh.nodeIndex(*pB))){
		double tmp_ttime=ttime_up(N0)+smesh.calc_ttime(N0,*pB);
		double orig_ttime=ttime_up(*pB);
 		if (orig_ttime>tmp_ttime){
		    ttime_up(*pB) = tmp_ttime;
		    prev_node_up(*pB) = N0;
		    B_up.promote(pB);
		}
	    }
	    pB++;  
	} 
	
	// loop over [FS and C] to calculate traveltime
	// and transfer them to B
	nodeIterator pC=C.begin();
	nodeIterator epC;
	while(pC!=C.end()){
	    if (fstar.isIn(smesh.nodeIndex(*pC))){
		ttime_up(*pC) = ttime_up(N0)+smesh.calc_ttime(N0,*pC);
		prev_node_up(*pC) = N0;
		B_up.push(*pC);
		epC = pC;
		pC++;
		C.erase(epC);
	    }else{ 
		pC++;
	    }
	}
    }
    
    is_refl_solved = true;
}

void GraphSolver2d::pickPath(const Point2d& rcv, Array1d<Point2d>& path) const
{
    if (!is_solved) error("GraphSolver2d::not yet solved");

    list<Point2d> tmp_path;
    tmp_path.push_back(rcv);

    int cur = smesh.nearest(rcv);
    int prev;
    int nadd=0;
    while ((prev=prev_node(cur)) != cur){
	Point2d p = smesh.nodePos(prev);
	tmp_path.push_back(p); nadd++;
	cur = prev;
    }

    if (nadd>0){
	tmp_path.pop_back();
    }
    if (nadd==1){
	// insert mid point for later bending refinement
	tmp_path.push_back(0.5*(src+rcv));
    }
    tmp_path.push_back(src);

    if (refine_if_long){
	list<Point2d>::iterator pt=tmp_path.begin();
	Point2d prev_p = *pt;
	pt++; // skip the first point
	while(pt!=tmp_path.end()){
	    double seg_len = prev_p.distance(*pt);
	    if (seg_len>crit_len){
		int ndiv = int(seg_len/crit_len+1.0);
		double frac = 1.0/ndiv;
		for (int j=1; j<ndiv; j++){
		    double ratio = frac*j;
		    Point2d new_p = (1.0-ratio)*prev_p+ratio*(*pt);
		    tmp_path.insert(pt,new_p);
		}
	    }
	    prev_p = *pt++;
	}
    }

    path.resize(tmp_path.size());
    list<Point2d>::const_iterator pt=tmp_path.begin();
    int i=1;
    while(pt!=tmp_path.end()){
	path(i++) = *pt++;
    }
}

void GraphSolver2d::pickPathThruWater(const Point2d& rcv, Array1d<Point2d>& path,
				      int& i0, int& i1) const
{
    if (!is_solved) error("GraphSolver2d::not yet solved");

    list<Point2d> tmp_path;
    tmp_path.push_back(rcv);

    int cur = smesh.nearest(rcv);
    findRayEntryPoint(rcv,cur,ttime);

    Point2d entryp=smesh.nodePos(cur);
    tmp_path.push_back(entryp);
    tmp_path.push_back(entryp);
    tmp_path.push_back(entryp);
    i0 = 2; i1 = 4;
    
    int prev;
    int nadd=0;
    while ((prev=prev_node(cur)) != cur){
	Point2d p = smesh.nodePos(prev);
	tmp_path.push_back(p); nadd++;
	cur = prev;
    }
    
    if (nadd>0){
	tmp_path.pop_back();
    }
    if (nadd==1){
	// insert mid point for later bending refinement
	tmp_path.push_back(0.5*(src+entryp));
    }
    tmp_path.push_back(src);

    if (refine_if_long){
	list<Point2d>::iterator pt=tmp_path.begin();
	pt++; pt++; pt++; // now at the entry point
	Point2d prev_p = *pt;
	pt++; // skip the entry point
	while(pt!=tmp_path.end()){
	    double seg_len = prev_p.distance(*pt);
	    if (seg_len>crit_len){
		int ndiv = int(seg_len/crit_len+1.0);
		double frac = 1.0/ndiv;
		for (int j=1; j<ndiv; j++){
		    double ratio = frac*j;
		    Point2d new_p = (1.0-ratio)*prev_p+ratio*(*pt);
		    tmp_path.insert(pt,new_p);
		}
	    }
	    prev_p = *pt++;
	}
    }
    
    path.resize(tmp_path.size());
    list<Point2d>::const_iterator pt=tmp_path.begin();
    int i=1;
    while(pt!=tmp_path.end()){
	path(i++) = *pt++;
    }
}

void GraphSolver2d::pickReflPath(const Point2d& rcv, Array1d<Point2d>& path,
				 int& ir0, int& ir1) const
{
    if (!is_refl_solved) error("GraphSolver2d::not yet solved");

    list<Point2d> tmp_path;
    tmp_path.push_back(rcv);

    int cur = smesh.nearest(rcv);
    int prev;
    int nadd=0;
    while ((prev=prev_node_up(cur)) != -1){
	Point2d p = smesh.nodePos(prev);
	tmp_path.push_back(p); nadd++; 
	cur = prev;
    }
    double reflx = smesh.nodePos(cur).x();
    double reflz = itfp->z(reflx);
    Point2d onrefl(reflx,reflz);

    if (nadd>0){
	tmp_path.pop_back();
    }
    if (nadd==1){
	// insert mid point for later bending refinement
	tmp_path.push_back(0.5*(onrefl+rcv));
    }

    tmp_path.push_back(onrefl); ir0 = tmp_path.size();
    tmp_path.push_back(onrefl);
    tmp_path.push_back(onrefl); ir1 = tmp_path.size();
    
    nadd=0;
    while ((prev=prev_node_down(cur)) != cur){
	Point2d p = smesh.nodePos(prev);
	tmp_path.push_back(p); nadd++;
	cur = prev;
    }
    if (nadd>0){
	tmp_path.pop_back();
    }
    if (nadd==1){
	// insert mid point for later bending refinement
	tmp_path.push_back(0.5*(src+rcv));
    }
    tmp_path.push_back(src);

    if (refine_if_long){
	int old_i=1, nadd_pre=0;
	list<Point2d>::iterator pt=tmp_path.begin();
	Point2d prev_p = *pt;
	pt++; old_i++; // skip the first point
	while(pt!=tmp_path.end()){
	    double seg_len = prev_p.distance(*pt);
	    if (seg_len>crit_len){
		int ndiv = int(seg_len/crit_len+1.0);
		double frac = 1.0/ndiv;
		for (int j=1; j<ndiv; j++){
		    double ratio = frac*j;
		    Point2d new_p = (1.0-ratio)*prev_p+ratio*(*pt);
		    tmp_path.insert(pt,new_p);
		    if (old_i<=ir0) nadd_pre++;
		}
	    }
	    prev_p = *pt++; old_i++;
	}
	ir0 += nadd_pre;
	ir1 += nadd_pre;
    }
    
    path.resize(tmp_path.size());
    list<Point2d>::const_iterator pt=tmp_path.begin();
    int i=1;
    while(pt!=tmp_path.end()){
	path(i++) = *pt++;
    }
}

void GraphSolver2d::pickReflPathThruWater(const Point2d& rcv, Array1d<Point2d>& path,
					  int& i0, int& i1, int& ir0, int& ir1) const
{
    if (!is_refl_solved) error("GraphSolver2d::not yet solved");

    list<Point2d> tmp_path;
    tmp_path.push_back(rcv);

    int cur = smesh.nearest(rcv);
    findRayEntryPoint(rcv,cur,ttime_up);

    Point2d entryp=smesh.nodePos(cur);
    tmp_path.push_back(entryp);
    tmp_path.push_back(entryp);
    tmp_path.push_back(entryp);
    i0 = 2; i1 = 4;
    
    int prev;
    int nadd=0;
    while ((prev=prev_node_up(cur)) != -1){
	Point2d p = smesh.nodePos(prev);
	tmp_path.push_back(p); nadd++; 
	cur = prev;
    }
    double reflx = smesh.nodePos(cur).x();
    double reflz = itfp->z(reflx);
    Point2d onrefl(reflx,reflz);

    if (nadd>0){
	tmp_path.pop_back();
    }
    if (nadd==1){
	// insert mid point for later bending refinement
	tmp_path.push_back(0.5*(onrefl+rcv));
    }

    tmp_path.push_back(onrefl); ir0 = tmp_path.size();
    tmp_path.push_back(onrefl);
    tmp_path.push_back(onrefl); ir1 = tmp_path.size();
    
    nadd=0;
    while ((prev=prev_node_down(cur)) != cur){
	Point2d p = smesh.nodePos(prev);
	tmp_path.push_back(p); nadd++;
	cur = prev;
    }
    if (nadd>0){
	tmp_path.pop_back();
    }
    if (nadd==1){
	// insert mid point for later bending refinement
	tmp_path.push_back(0.5*(src+rcv));
    }
    tmp_path.push_back(src);

    if (refine_if_long){
	int nadd_pre=0;
	list<Point2d>::iterator pt=tmp_path.begin();
	pt++; pt++; pt++; // now at the entry point
	Point2d prev_p = *pt;
	pt++; // skip the entry point
	int old_i = 5;
	while(pt!=tmp_path.end()){
	    double seg_len = prev_p.distance(*pt);
	    if (seg_len>crit_len){
		int ndiv = int(seg_len/crit_len+1.0);
		double frac = 1.0/ndiv;
		for (int j=1; j<ndiv; j++){
		    double ratio = frac*j;
		    Point2d new_p = (1.0-ratio)*prev_p+ratio*(*pt);
		    tmp_path.insert(pt,new_p);
		    if (old_i<=ir0) nadd_pre++;
		}
	    }
	    prev_p = *pt++; old_i++;
	}
	ir0 += nadd_pre;
	ir1 += nadd_pre;
    }

    path.resize(tmp_path.size());
    list<Point2d>::const_iterator pt=tmp_path.begin();
    int i=1;
    while(pt!=tmp_path.end()){
	path(i++) = *pt++;
    }
}

void GraphSolver2d::findRayEntryPoint(const Point2d& rcv, int& cur, const Array1d<double>& nodetime) const
{
    Point2d p0 = smesh.nodePos(cur);
    double pwater = smesh.atWater();
    double ttime0 = pwater*rcv.distance(p0)+nodetime(cur);

    // local grid search for best ray entry point (+/- 10km radius)
    Point2d pmin(rcv.x()-local_search_radius,rcv.y());
    Point2d pmax(rcv.x()+local_search_radius,rcv.y()); 
    Index2d indmin = smesh.nodeIndex(smesh.nearest(pmin)); 
    Index2d indmax = smesh.nodeIndex(smesh.nearest(pmax)); 
    for (int i=indmin.i()+1; i<indmax.i(); i++){
	int inode = smesh.nodeIndex(i,1);
	Point2d p = smesh.nodePos(inode);
	double tmp_ttime = pwater*rcv.distance(p)+nodetime(inode);
	if (tmp_ttime < ttime0){
	    ttime0 = tmp_ttime;
	    cur = inode; 
	}
    }
}

void GraphSolver2d::printPath(ostream& os) const
{
    for (int cur=1; cur<=nnodes; cur++){
	int prev = prev_node(cur);
	Point2d cur_p = smesh.nodePos(cur);
	Point2d prev_p = smesh.nodePos(prev);
	    os << ">\n"
	       << cur_p.x() << " "  << cur_p.y() << "\n"
	       << prev_p.x() << " " << prev_p.y() << '\n';
    }
}

//
// forward star
//
ForwardStar2d::ForwardStar2d(const Index2d& id, int ix, int iz)
    : orig(id)
{
    if (ix>0 && iz>0){ xorder=ix; zorder=iz; }
    else{ error("ForwardStar2d::non-positive order detected."); }
}

bool ForwardStar2d::isIn(const Index2d& node) const
{
    int diff;
    int di = (diff=orig.i()-node.i()) > 0 ? diff : -diff;
    int dk = (diff=orig.k()-node.k()) > 0 ? diff : -diff;

    if (di>xorder || dk>zorder) return false;
    if (di==1 || dk==1) return true;
    if (di==0 || dk==0) return false; // NB: (0,1),(1,0) are already returned true.

    // now di && dk must be >=2 
    int mod = dk>di ? dk%di : di%dk;
    return mod==0 ? false : true;
}
