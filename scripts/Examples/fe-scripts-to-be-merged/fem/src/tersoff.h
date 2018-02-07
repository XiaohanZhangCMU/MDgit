/*
  tersoff.h
  by Wei Cai  caiwei@stanford.edu
  Last Modified : Thu Jul 19 18:51:42 2007

  FUNCTION  :  MD simulation package of Si/Ge/C using Tersoff potential
  Reference:   J. Tersoff, Modeling solid-state chemistry: Interatomic
               potentials for multicomponent systems,
               Phys. Rev. B, 39, 5566 (1989).

  Note:        Implemented with the help of the code from KJ Cho's group
               created: 08/26/2001 Ram (panchram@stanford.edu)
               revised: Byeongchan Lee (bclee@stanford.edu)
*/

#ifndef _TERSOFF_H
#define _TERSOFF_H

#include "mdparallel.h"

class TersoffFrame : public MDPARALLELFrame /* Si with Tersoff potential */
{
    /* general parameters */
    double acut, acutsq;
    
    /* Tersoff potential parameters */
    double A,A_i,A_j;
    double B,B_i,B_j;
    double Lam,Lam_i,Lam_j;
    double Mu,Mu_i,Mu_j;   /* Same as Lam2 in the ref #7 paper. */
    double R,R_i,R_j,R_k;
    double S,S_i,S_j,S_k;
    double Beta;
    double nTf;
    double cTf;
    double dTf;
    double hTf;
    

    /* Fitting parameters for Si */
    double A_Si;
    double B_Si;
    double Lam_Si;
    double Mu_Si;
    double R_Si;
    double S_Si;
    double Beta_Si;
    double nTf_Si;
    double cTf_Si;
    double dTf_Si;
    double hTf_Si;

    /* Fitting parameters for C */
    double A_C;
    double B_C;
    double Lam_C;
    double Mu_C;
    double R_C;
    double S_C;
    double Beta_C;
    double nTf_C;
    double cTf_C;
    double dTf_C;
    double hTf_C;
                                                                                                                     
    /* Fitting parameters for Ge */                                                                              
    double A_Ge;
    double B_Ge;
    double Lam_Ge;
    double Mu_Ge;
    double R_Ge;
    double S_Ge;
    double Beta_Ge;
    double nTf_Ge;
    double cTf_Ge;
    double dTf_Ge;
    double hTf_Ge;

    /* Inter species parameter */                                                                                
    double Chi;                                                                                                    
    double Chi_C_Si;
    double Chi_Si_Ge;
    double Chi_C_Ge;

public:
    TersoffFrame(){};
    void tersoff();
    virtual void potential();
    virtual void initvars();
};

#endif // _TERSOFF_H

