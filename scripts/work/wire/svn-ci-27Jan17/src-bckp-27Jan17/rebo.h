/*
  rebo.h
  by Wei Cai  caiwei@stanford.edu
  Last Modified : Thu Jul 19 18:51:42 2007

  FUNCTION  :  MD simulation package of C using REBO potential
  Reference:   

  Note:        Implemented based on LAMMPS (7/24/09)

*/

#ifndef _REBO_H
#define _REBO_H

#include "mdparallel.h"

class REBOFrame : public MDPARALLELFrame /* Carbon, H with Brenner potential */
{
  /* general parameters */

  double acut, acutsq;
    
  /* REBO potential parameters (from LAMMPS) */

  double cut3rebo;                 // maximum distance for 3rd REBO neigh

  int *REBO_numneigh;              // # of pair neighbors for each atom
  int **REBO_firstneigh;           // ptr to 1st neighbor of each atom
  char *REBO_firstneigh_mem;       // memory for REBO_firstneigh
  double *closestdistsq;           // closest owned atom dist to each ghost
  double *nC,*nH;                  // sum of weighting fns with REBO neighs

  double smin,Nmin,Nmax,NCmin,NCmax,thmin,thmax;
  double rcmin[2][2],rcmax[2][2],rcmaxsq[2][2],rcmaxp[2][2];
  double Q[2][2],alpha[2][2],A[2][2],rho[2][2],BIJc[2][2][3],Beta[2][2][3];
  double rcLJmin[2][2],rcLJmax[2][2],rcLJmaxsq[2][2],bLJmin[2][2],bLJmax[2][2];
  double epsilon[2][2],sigma[2][2],epsilonT[2][2];

  // spline coefficients

  double gCdom[5],gC1[4][6],gC2[4][6],gHdom[4],gH[3][6];
  double pCCdom[2][2],pCHdom[2][2],pCC[4][4][16],pCH[4][4][16];
  double piCCdom[3][2],piCHdom[3][2],piHHdom[3][2];
  double piCC[4][4][9][64],piCH[4][4][9][64],piHH[4][4][9][64];
  double Tijdom[3][2],Tijc[4][4][9][64];

  // spline knot values

  double PCCf[5][5],PCCdfdx[5][5],PCCdfdy[5][5],PCHf[5][5];
  double PCHdfdx[5][5],PCHdfdy[5][5];
  double piCCf[5][5][11],piCCdfdx[5][5][11];
  double piCCdfdy[5][5][11],piCCdfdz[5][5][11];
  double piCHf[5][5][11],piCHdfdx[5][5][11];
  double piCHdfdy[5][5][11],piCHdfdz[5][5][11];
  double piHHf[5][5][11],piHHdfdx[5][5][11];
  double piHHdfdy[5][5][11],piHHdfdz[5][5][11];
  double Tf[5][5][10],Tdfdx[5][5][10],Tdfdy[5][5][10],Tdfdz[5][5][10];

  /* internal functions */

  double bondorder(int, int, double *, double, double, Vector3 *, int);

  double Sp(double, double, double, double &);
  double Sp2(double, double, double, double &);

  double gSpline(double, double, int, double *, double *);
  double PijSpline(double, double, int, int, double *);
  double piRCSpline(double, double, double, int, int, double *);
  double TijSpline(double, double, double, double *);

  double kronecker(int, int);

  double Sp5th(double, double *, double *);
  double Spbicubic(double, double, double *, double *);
  double Sptricubic(double, double, double, double *, double *);
  void spline_init();
  void read_file(char *);

private:
    char rebofile[1000]; /* potential file */

protected:

  int eflag;          // these variables determine whether or not to
  int evflag;         // accumulate energy and virial (default 1)

  void ev_tally(int, int, double, double, double,
		double, double, double);
  void v_tally2(int, int, double, double *);
  void v_tally3(int, int, int, double *, double *, double *, double *);
  void v_tally4(int, int, int, int, double *, double *, double *,
		double *, double *, double *);

public:
    REBOFrame():acut(3.0),eflag(1),evflag(1) {};

    double eng_vdwl,eng_coul;
    void rebo();
    void REBO_neigh();
    int readAIREBO();
    virtual void Alloc();
    virtual void potential();
    virtual void initvars();
    virtual void initparser();
    
    virtual int exec(const char *nam);

};

#endif // _REBO_H

