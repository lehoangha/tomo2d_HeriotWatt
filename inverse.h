/*
 * inverse.h
 *
 * Jun Korenaga, MIT/WHOI
 * January 1999
 */

#ifndef _TOMO_INVERSE_H_
#define _TOMO_INVERSE_H_

#include <iostream>
#include <fstream>
#include <map>
#include <array.h>
#include "smesh.h"
#include "graph.h"
#include "betaspline.h"
#include "bend.h"
#include "sparse_rect.h"
#include "interface.h"
#include "corrlen.h"
#include "jgrav.h"

class TomographicInversion2d {
public:
    TomographicInversion2d(SlownessMesh2d& m, const char* datafn,
			   int xorder, int zorder, double crit_len,
			   int nintp, double cg_tol, double br_tol);
    void solve(int niter);
    void doRobust(double);
    void removeOutliers();
    void setLSQR_TOL(double);
    void addRefl(Interface2d* intfp);
    void doFullRefl();
    void setReflWeight(double);
    
    void SmoothVelocity(const char*, double, double, double,
			bool logscale=false);
    void SmoothDepth(double, double, double,
		     bool logscale=false);
    void SmoothDepth(const char*, double, double, double,
		     bool logscale=false);
    void applyFilter(const char*);
    void DampVelocity(double);
    void DampDepth(double);
    void FixDamping(double, double);
    void Squeezing(const char*);
    void targetChisq(double);
    void doJumping();
    void addonGravity(double, const Array1d<double>&, const Array1d<double>&,
		      AddonGravityInversion2d*, double);
    void addonGravity(double, const Array1d<double>&, const Array1d<double>&,
		      AddonGravityInversion2d*, double, const char*);

    void outStepwise(const char*, int level);
    void outFinal(const char*, int level);
    void setLogfile(const char*);
    void outMask(const char*);
    void setVerbose(int);
    
private:
    typedef Array1d< map<int,double> > sparseMat;

    void read_file(const char* ifn);

    void reset_kernel();
    void add_kernel(int, const Array1d<Point2d>&);
    void add_kernel_refl(int, const Array1d<Point2d>&, int, int);
    void calc_averaging_matrix();
    void calc_damping_matrix();
    void calc_refl_averaging_matrix();
    void calc_refl_damping_matrix();
    void add_global(const Array2d<double>&, const Array1d<int>&,
		    sparseMat&);
    int _solve(bool,double,bool,double,bool,double,bool,double);
    void auto_damping(int&, int&, double&, double&);
    void fixed_damping(int&, int&, double&, double&);
    double auto_damping_depth(double, int&, int&);
    double auto_damping_vel(double, int&, int&);
    double calc_ave_dmv();
    double calc_ave_dmd();
    void calc_Lm(double&, double&, double&);
    double calc_chi();

    SlownessMesh2d& smesh;
    GraphSolver2d graph;
    BetaSpline2d betasp;
    BendingSolver2d bend;
    
    Array1d<Point2d> src;
    Array1d< Array1d<Point2d> > rcv;
    Array1d< Array1d<int> > raytype;
    Array1d< Array1d<double> > obs_ttime;
    Array1d< Array1d<double> > obs_dt;
    Array1d< Array1d<double> > res_ttime;
    Array1d<double> r_dt_vec, path_wt;
    Array1d<double> path_length;

    int nrefl;
    Array1d<int> start_i, end_i;
    Array1d<const Interface2d*> interp;
    const Interface2d *bathyp;
    Interface2d *reflp;
    double refl_weight;
    bool do_full_refl;

    int nnodev, nnoded, ndata, ndata_valid, nx, nz;
    double rms_tres[2], init_chi[2], rms_tres_total, init_chi_total;
    int ndata_in[2];
    bool robust;
    double crit_chi;
    sparseMat A, Rv_h, Rv_v, Rd, Tv, Td;
    Array1d<double> data_vec, total_data_vec;
    Array1d<double> modelv, modeld, dmodel_total;
    int itermax_LSQR;
    double LSQR_ATOL;

    bool jumping;
    Array1d<double> dmodel_total_sum, mvscale, mdscale;

    bool smooth_velocity, logscale_vel;
    double wsv_min, wsv_max, dwsv;
    double weight_s_v;
    CorrelationLength2d *corr_vel_p;
    bool do_filter;
    Interface2d *uboundp;

    bool smooth_depth, logscale_dep;
    double wsd_min, wsd_max, dwsd;
    double weight_s_d;

    CorrelationLength1d *corr_dep_p;
    bool damp_velocity, damp_depth;
    double target_dv, target_dd;
    bool damping_is_fixed; 
    double weight_d_v, weight_d_d;
    bool do_squeezing;
    DampingWeight2d *damping_wt_p;

    int nnode_total;
    Array1d<int> nodev_hit;
    Array1d<int> tmp_node, tmp_nodev, tmp_data;
    Array1d<int> tmp_nodedc, tmp_nodedr;
    bool out_mask;
    Array1d<double> dws;
    ofstream* vmesh_os_p;

    double target_chisq;
    
    bool printLog;
    ofstream* log_os_p;
    int verbose_level;
    bool printTransient, printFinal;
    const char* out_root;
    int out_level;

    bool gravity, out_grav_dws;
    AddonGravityInversion2d *ginv;
    int ngravdata;
    sparseMat B;
    double weight_grav, grav_z0, rms_grav;
    Array1d<int> tmp_gravdata;
    Array1d<double> grav_x, obs_grav, res_grav;
    ofstream *grav_dws_osp;
    
    // variables of temporary use
    Array1d<Point2d> path, Q;
    Array1d<const Point2d*> pp;

    ofstream* dump_os_p;
};

#endif /* _TOMO_INVERSE_H_ */
