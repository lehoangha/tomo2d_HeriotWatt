/*
 * syngen.h - forward traveltime calculation
 *
 * Jun Korenaga, MIT/WHOI
 * January 1999
 */

#ifndef _TOMO_SYNGEN_H_
#define _TOMO_SYNGEN_H_

#include <array.h>
#include <geom.h>
#include "smesh.h"
#include "graph.h"
#include "bend.h"
#include "interface.h"

class SyntheticTraveltimeGenerator2d {
public:
    SyntheticTraveltimeGenerator2d(SlownessMesh2d& m, const char* ifn,
				   int xorder=3, int zorder=3, double clen=0.0,
				   int nintp=8, double cg_tol=1e-4, double br_tol=1e-7);

    void conduct();
    void readRefl(const char *);
    void doFullRefl();
    void outputRay(const char* fn);
    void useClock(const char* fn);
    void graphOnly();
    void setVerbose(int i=0);
    
    void printSource(ostream&) const;
    void printSynTime(ostream&, double) const;
    void printObsTime(ostream&, double) const;
    void printDiff(ostream&, double& misfit, double& chisq) const;

    friend ostream& operator<<(ostream&,
			       const SyntheticTraveltimeGenerator2d&);
    
private:
    void read_file(const char* ifn);
    
    SlownessMesh2d& smesh;
    GraphSolver2d graph;
    BetaSpline2d betasp;
    BendingSolver2d bend;
    
    Array1d<Point2d> src;
    Array1d< Array1d<Point2d> > rcv;
    Array1d< Array1d<int> > raytype;
    Array1d< Array1d<double> > obs_ttime, syn_ttime;
    Array1d< Array1d<double> > obs_dt;

    int nrefl;
    Array1d<int> start_i, end_i;
    Array1d<const Interface2d*> interp;
    const Interface2d *bathyp, *reflp;
    bool do_full_refl;
    
    bool outray;
    const char* rayfn;

    bool use_clock, graph_only;
    const char* timefn;
    double graph_time, bend_time;
    
    int verbose_level;
};

ostream&
operator<<(ostream&, const SyntheticTraveltimeGenerator2d&);
    
#endif /* _TOMO_SYNGEN_H_ */
