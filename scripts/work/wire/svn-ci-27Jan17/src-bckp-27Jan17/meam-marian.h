/*
  meam-marian.h
  by Wei Cai  caiwei@stanford.edu
  Last Modified : Sat Dec 30 21:19:15 2006

  FUNCTION  : MD++ with MEAM potential (using MDCASK F90 Codes)
*/

#ifndef _MEAM_MARIAN_H
#define _MEAM_MARIAN_H

#include "mdparallel.h"

extern "C" int forces_();
extern "C" int meamdefine_();
extern "C" int input_();

#define nmmax  100000
#define nlcmax 100000
#define ndim  3
#define nnbrs 400
#define nnbrs2 500

extern struct {char meamfile[100]; } meamdatafile_;
extern struct {int nlcx, nlcy, nlcz, nlc; } linklist_;
extern struct {double Lx, Ly, Lz; } box_;
extern struct {int nm; } np_;
extern struct {double x0[nmmax], y0[nmmax], z0[nmmax];} adim_;
extern struct {double atpe[nmmax], atpe3b[nmmax], pe, petrip; } ppote_;
extern struct {double fx[nmmax], fy[nmmax], fz[nmmax]; } area6_;
extern struct {double zsmeam, alphas, betas[4], esubs, asubs, ts[4],
               rozros, rcutmeam, cmin, cmax, replus, attrac, legend;
               double equidist;  } meamread_;

//extern struct {double xi[nmmax][3], x[nmmax][ndim]; } pos_;
//extern struct {double rhotot[nmmax], embf[nmmax], embfp[nmmax];} embed_;
//extern struct {double fplus[nmmax][ndim], fmin[nmmax][ndim], f[nmmax][ndim]; } force_;
//extern struct {double xnbr[nnbrs], ynbr[nnbrs], znbr[nnbrs]; } area18_;
//extern struct {int ltop[nlcmax], jaddr[nnbrs], linkmp[nmmax]; } lnmap_;
//extern struct {int ibarr, noscr, ncrys; } meamread_int_;
//extern struct {double tav[nmmax][3], c8a[nmmax], dang2[nmmax], dang1[nmmax], rhsq[nmmax][3], cg8c[nmmax], a8b[nmmax][3], ag[nmmax][3], b8c[nmmax][3][3], d8d[nmmax][3][3][3], res, rcutmeam2; } meamforces_;
//extern struct {double scrnab[nnbrs2][nmmax]; } screendata_;
//extern struct {double dscrnab[nnbrs2][nmmax][3]; } screendife_;
//extern struct {int numneigh[nmmax], neighlabel[nnbrs2][nmmax]; } neighborlist_;
//extern struct {int id; } identifier_;

class MEAMFrame : public MDPARALLELFrame
{ /* MEAM potential */
public:
    
    MEAMFrame()  {};

    virtual void potential ();

    void readMEAM();
    void MEAM();

    virtual void Alloc();
    virtual void initvars();
    virtual void initparser();
    
    virtual int exec(const char *nam);

};


#endif // _MEAM-MARIAN_H

