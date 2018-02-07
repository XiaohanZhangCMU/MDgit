/*
  alglue.h
  by Wei Cai  caiwei@mit.edu, Maurice de Koning maurice@mit.edu
  Last Modified : Sat Dec 30 21:07:50 2006

  FUNCTION  : MD++ implementation of Al Glue potential
              (Ercolessi-Adams Glue potential)
*/

#ifndef _ALGLUE_H
#define _ALGLUE_H

#include "mdparallel.h"

extern void v2(double, double *, double *, double *);
extern void rh(double, double *, double *, double *);
extern void uu(double, double *, double *, double *);

class ALGLUEFrame : public MDPARALLELFrame
{ /* ALGLUE potential */
public:
    int pottype; /* 1: EAM, 2: MEAM */

    /* ALGLE potential parameters */
#define NGRID 20001
    double invdeltarsq,invdeltacoord,deltacoord;
    double phitab[NGRID],dphitab[NGRID],rhotab[NGRID],drhotab[NGRID];
    double utab[NGRID],dutab[NGRID];
    double rcutoff,rmin,rsqmin,deltarsq,coordmax;

    double *deru, *coord;
    Vector3 *rijstore;
    
    ALGLUEFrame():deru(0),coord(0),rijstore(0) {};

    virtual void potential ();
    int initglue();
    int al_glue();

    virtual void Alloc();
    virtual void initvars();
    virtual void initparser();
    
    virtual int exec(const char *nam);

};


#endif // _EAM_H

