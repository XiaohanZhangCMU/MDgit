/*
  silica-bksmod.cpp
  by Wei Cai  caiwei@mit.edu
  Last Modified : Fri Feb 16 09:59:52 2007

  FUNCTION  :  Modified BKS model for Silica (SiO2)
               Si: species = 0
               O : species = 1

  The original BKS potential is described in
   B. W. H. van Beest, G. J. Kramer and R. A. van Santen
   Phys. Rev. Lett. 16, 1955-1958 (1990). 
   Force fields for silicas and aluminophosphates based on 
    ab initio calculations

  The modified version is described in
   S. S. Mahajan, G. Subbarayan and B. G. Sammakia,
   Phys. Rev. E 76, 056701 (2007).
   Estimating thermal conductivity of amorphous silica nanoparticles 
    and nanowires using molecular dynamics simulations

  Short range repulsion is added.  Coulomb interaction is truncated by erfc().
*/

#ifndef _SILICA_BKSMOD_H
#define _SILICA_BKMMOD_H

#include "mdparallel.h"

#define P_E 1.6021892e-19  // (C) electron charge
#define N_A 6.0221415e+23  // Avogadro's number  

class BKSFrame : public MDPARALLELFrame 
{
    /* BKS parameters */
    double _A_00,_B_00,_C_00; /* Si-Si */
    double _A_01,_B_01,_C_01; /* Si-O  */
    double _A_11,_B_11,_C_11; /*  O-O  */
    double _Q_0,_Q_1;         /* charge for Si and O */ 
    double _ALPHA_BKS;        /* cut-off parameter for Coulomb interaction */

    double BKS_RC;            /* cut-off radius */

    double _SIG_00,_EPS_00;
    double _SIG_01,_EPS_01;
    double _SIG_11,_EPS_11;

public:
    BKSFrame(){};
    void bks_mod();
    virtual void potential();
    virtual void initvars();

    virtual void initparser();
    virtual int exec(const char *name);
    virtual void calcprop();
    
    void initBKS();
    
};



#endif // _SILICA_BKSMOD_H

