/*
  tersoff88.h
  by Keonwook Kang kwkang75@yonsei.ac.kr
  Last Modified : 

  FUNCTION  :  
  Reference:   J. Tersoff, New empirical approach for the structure and 
               energy of covalent systems
               Phys. Rev. B, 37, 6991 (1988), Si(B)
               J. Tersoff, Empirical interatomic potential for silicon 
               with improved elastic properties
               Phys. Rev. B, 38, 9902 (1988), Si(C)

  Note:        Implemented by modifying tersoff.h
*/

#ifndef _TERSOFF_H
#define _TERSOFF_H

#include "mdparallel.h"

class Tersoff88Frame : public MDPARALLELFrame /* Si with Tersoff potential */
{
    /* general parameters */
    double acut, acutsq;
    
    /* Tersoff potential parameters */
    double A,A_i,A_j;
    double B,B_i,B_j;
    double Lam,Lam_i,Lam_j;
    double Mu,Mu_i,Mu_j;   /* Same as Lam2 in the ref #7 paper. */
    double Lam3; 
    double R,R_i,R_j,R_k;
    double S,S_i,S_j,S_k;
    double Alpha;
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
    double Lam3_Si;
    double R_Si;
    double S_Si;
    double Alpha_Si;
    double Beta_Si;
    double nTf_Si;
    double cTf_Si;
    double dTf_Si;
    double hTf_Si;

    /* Inter species parameter */                                                                                
    double Chi, Chi2;    // Chi2 is the mixing parameter for a_ij                                                                                                
    double Chi_Si_Si, Chi2_Si_Si;
    
public:
    Tersoff88Frame(){};
    void tersoff();
    virtual void potential();
    virtual void initvars();
};

#endif // _TERSOFF_H

