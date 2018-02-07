/*
  airebo.h
  by Keonwook Kang  kwkang75@yonsei.ac.kr
  Last Modified : Thurs Feb 04 11:28:54 2014

  FUNCTION  :  MD simulation package of C using AI-REBO potential
  Reference:   
  (1) D. W. Brenner et al, J. Phys.: Condens. Matter 14 (2002) 783

  Note:        Implemented based on LAMMPS (02/04/14)

*/

#ifndef _AIREBO_H
#define _AIREBO_H

#include "mdparallel.h"

class AIREBOFrame : public MDPARALLELFrame /* Carbon, H with Brenner potential */
{
  /* general parameters */

  double acut, acutsq;              // cut-off for REBO, in unit of Angstrom
    
  /* REBO potential parameters (from LAMMPS) */

  double cut3rebo;                 // maximum distance for 3rd REBO neigh

  int *REBO_numneigh;              // # of pair neighbors for each atom
  int **REBO_firstneigh;           // ptr to 1st neighbor of each atom
  char *REBO_firstneigh_mem;       // memory for REBO_firstneigh
  double *closestdistsq;           // closest owned atom dist to each ghost
  double *nC,*nH;                  // sum of weighting fns with REBO neighs
  class Vector3 *REBO_R;
  
  double smin,Nmin,Nmax,NCmin,NCmax,thmin,thmax;
  double rcmin[2][2],rcmax[2][2],rcmaxsq[2][2],rcmaxp[2][2];
  double Q[2][2],alpha[2][2],A[2][2],rho[2][2],BIJc[2][2][3],Beta[2][2][3];

  // LJ parameters
  int    ljflag;                    // If ljflag=1, lj contribution considered
  double cutlj;                     // cut-off for LJ, in unit of sigma
  double rcLJmin[2][2],rcLJmax[2][2],rcLJmaxsq[2][2],bLJmin[2][2],bLJmax[2][2];
  double epsilon[2][2],sigma[2][2],epsilonT[2][2];
  double cutljsq[2][2],lj1[2][2],lj2[2][2],lj3[2][2],lj4[2][2];

  // TORSION parameters
  // If torsionflag=1, torsion contribution considered
  int    torsionflag;                    
  
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
  double bondorderLJ(int, int, double *, double, double, double *, double, Vector3 *, int);

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
  void ev_tally4(int, int, int, int, double, double *, double *,
		double *, double *, double *, double *);
  void v_tally2(int, int, double, double *);
  void v_tally3(int, int, int, double *, double *, double *, double *);
  void v_tally4(int, int, int, int, double *, double *, double *,
		double *, double *, double *);
  void virial_fdotr_compute();

public:
    AIREBOFrame():acut(3.0),cutlj(3.0),eflag(1),evflag(1),
                  ljflag(1),torsionflag(1) {};

    double eng_vdwl,eng_coul;
    void airebo(),rebo(),flj(),torsion();
    void REBO_neigh();
    int readAIREBO();
    virtual void Alloc();
    virtual void potential();
    virtual void initvars();
    virtual void initparser();
    
    virtual int exec(const char *nam);

    static inline double powint(const double &x, const int n) 
    {
        double yy,ww;

        if (x == 0.0) return 0.0;
        int nn = (n > 0) ? n : -n;
        ww = x;

        for (yy = 1.0; nn != 0; nn >>= 1, ww *=ww)
            if (nn & 1) yy *= ww;

        return (n > 0) ? yy : 1.0/yy;
    }

};

#endif // _REBO_H

