/*
  meam-baskes.h
  by Wei Cai  caiwei@stanford.edu
  Last Modified : Mon Dec  3 16:34:05 2007

  FUNCTION  : MD++ with MEAM potential (using MDCASK Fortran Codes of Mike Baskes baskes@lanl.gov)
*/

#ifndef _MEAM_BASKES_H 
#define _MEAM_BASKES_H 

#include "mdparallel.h"
/* This is no longer necessary because we renamed cal_phiid to calphiid
 * In this case only one _ will be added in the end, which is consistent
 * with the convention.
 */
//#ifdef _CYGWIN
//#define cal_phiid_ cal_phiid__
//#endif

extern "C" int force_();
extern "C" int inter_();             /* read MEAM parameters */
extern "C" int phif_(double *,int *,int *);
extern "C" int phiid_(double *,int *);
extern "C" int calphiid_(double *,int *,double *);

/* need to match utils/param.inc */
#define natmax  10000
#define neimax  150
#define ngrid   5010
#define ngridar 5010
#define ntable  10000
#define nelmax  3
#define natcmax 500
#define icmax   30
#define icmay   30
#define icmaz   30
#define neimx   3000

/* utils/utils.inc */
extern struct INTERACT { double rcutsq, rcut; } interact_;
extern struct LATTICE  { double perub[3], perlb[3], perlen[3], alat, xbound[2], ybound[2], zbound[2];
                    int natoms; double latty; } lattice_;
extern struct PARTICLE { double rv[natmax][6]; int itype[natmax]; } particle_;
extern struct SCALED   { double y[natmax][6]; } scaled_;
extern struct FORCES   { double f[natmax][3], e[natmax], stresst[3][3], slocal[natmax][3][3],
                    avslocal[natmax][3][3]; int noutput; } forces_;
extern struct TYPES    { double amass[nelmax]; int ielement[nelmax], netype[nelmax], ntypes ; } types_;
extern struct CNEIGH   { double rctsqn, dradn; int nneigh, nneips, nmeth, nnindx[natmax+1],
                    newlst, nlstmx, nnlst[natmax*neimax]; } cneigh_;

/* meam/meam.inc */
extern struct CMEAM { double omegas[nelmax], res[nelmax][nelmax],
                    esubs[nelmax][nelmax], alats[nelmax], zs[nelmax],
                    alphas[nelmax][nelmax],deltas[nelmax][nelmax],
                    betas[nelmax][4],ts[nelmax][4],asubs[nelmax],
                    rozros[nelmax],z0s[nelmax],
                    repuls[nelmax][nelmax],attrac[nelmax][nelmax];
                float legend; int ibarr[nelmax],irhos[nelmax];
                char enames[nelmax][4],lattces[nelmax][4]; } cmeam_;

extern struct SCR { double rscrmax, cmin[nelmax][nelmax][nelmax],
                    cmax[nelmax][nelmax][nelmax]; int nscree,ipscree; } scr_;

/* meam/dyn.inc */
extern struct TWONN { double a2nn[nelmax], b2nn[nelmax]; 
                      bool nn; } twonn_;  // 2NN

/* meam/dyn88.f */
extern struct MEA {char meamf[80],meafile[80]; } mea_;
extern struct KODES {char kodes[nelmax][7]; } kodes_;
extern struct MYOUTPUT {int fid; } myoutput_;

class MEAMFrame : public MDPARALLELFrame
{ /* MEAM potential */
public:
    double cmin0[nelmax];
    
    MEAMFrame()  {};

    virtual void potential ();

    void readMEAM();
    void MEAM();

    virtual void Alloc();
    virtual void initvars();
    virtual void initparser();
    
    virtual int exec(const char *nam);

    void printpairpot();
};


#endif // _MEAM-BASKES_H

