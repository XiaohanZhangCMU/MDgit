/*
  siedip.h
  by Wei Cai  caiwei@mit.edu
  Last Modified : Thu Jul 19 18:54:22 2007

  FUNCTION  :  MD simulation package of Si using EDIP potential
*/

#ifndef _SIEDIP_H
#define _SIEDIP_H

#include "mdparallel.h"

class SIEDIPFrame : public MDPARALLELFrame /* Si with EDIP potential */
{
    /* DEBUGGING FLAGS - systematically turn on various pieces of the potential. */
#define V2_on	1	
#define V2Z_on	1	
#define V3_on 	1	
#define V3g_on	1	
#define V3h_on	1		
#define V3Z_on	1 
#define Zfast	1
    
    /* EDIP potential parameters */
    double A,B,rh,sig,a;
    double lam,gam,b,c;
    double mu,Qo,bet,alp;
    double bg;   /* cutoff for g(r) */
    double palp;    /* justo prefactor for bond order - delete later */
    double delta,eta,zet;

    /* tau(Z) (Ismail & Kaxiras, 1993) */
    double u1, u2, u3, u4;

public:
    SIEDIPFrame(){};
    void edip();
    virtual void potential();
    virtual void initvars();
};

#endif // _SIEDIP_H

