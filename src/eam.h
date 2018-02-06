/*
  eam.h
  by Wei Cai  caiwei@mit.edu
  Last Modified : Mon Sep  8 18:41:56 2008

  FUNCTION  : MD++ implementation of EAM/MEAM potential
*/

#ifndef _EAM_H
#define _EAM_H

#include "mdparallel.h"

class EAMFrame : public MDPARALLELFrame
{ /* EAM/MEAM potential */
public:
    int pottype; /* 1: EAM, 2: MEAM */
    int eamfiletype; /* 1: single element, 2: binary element */
    int eamgrid;
    /* EAM potential parameters */
#define NGRID 10000
    char title3[500];
    int ntyp;
    double rmass, rlatt, drar, drhoar, actual, actual2, rmin;
    double rho[4][NGRID], rhop[4][NGRID], phi[2][NGRID], phip[2][NGRID];
    double phix[NGRID], phipx[NGRID];
    double frho[2][NGRID], frhop[2][NGRID];
    double rho_spline[4][NGRID][4], phi_spline[2][NGRID][4];
    double phix_spline[NGRID][4], frho_spline[2][NGRID][4];
    double rval[NGRID], rhoval[NGRID];
    double petrip, rhocon, rhomax;
    double *atpe3b, *rhotot, *embf, *embfp;
    double *atpe3b1, *rhotot1, *embf1;
    double *atpe3bUMB, *rhototUMB, *embfUMB;
    double MC_maxd_trial;
    int MC_need, *nbst, *nbst1;
   
    int NumOfChanged,*ListOfChanged;
    
    EAMFrame() : pottype(1), eamfiletype(0), eamgrid(NGRID),
                 ntyp(0), rmass(0), rlatt(0), drar(0), drhoar(0),
                 actual(0), actual2(0), rmin(0),
                 petrip(0), rhocon(0), rhomax(0),
                 atpe3b(0), rhotot(0), embf(0), embfp(0),
                 atpe3b1(0), rhotot1(0), embf1(0),
                 atpe3bUMB(0), rhototUMB(0), embfUMB(0),
                 MC_maxd_trial(0), MC_need(0),
                 nbst(0), nbst1(0) {};

    virtual double potential_energyonly_before(int iatom);
    virtual double potential_energyonly_after(int iatom); 
    virtual void potential ();
    virtual void potential_energyonly();
    virtual double potential_energyonly(int iatom);
    virtual double potential_energyonly_change(int iatom); 
    void kraeam();
    void rhoeam();
    void kraMEAM();
    void rhoMEAM();
    int readeam();
    int readMEAM();

    virtual void Alloc();
    virtual void MC_Alloc();
    virtual void MC_Recover();
    virtual void MC_Update();

    virtual void MC_AllocUMB();
    virtual void MC_RecoverUMB();
    virtual void MC_UpdateUMB();

    virtual int MC_FFSprocess();

    virtual void initvars();
    virtual void initparser();
    
    virtual int exec(const char *nam);

#ifdef _PARALLEL
    void Broadcast_EAM_Param();
#endif
};


#endif // _EAM_H

