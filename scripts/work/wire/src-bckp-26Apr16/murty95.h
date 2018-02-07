/*
  murty95.h
  by Keonwook Kang kwkang75@yonsei.ac.kr
  Last Modified : 

  FUNCTION  :  MD simulation package of Si-H interaction
  Reference:   M. V. Ramana Murty and Harry A. Atwater
               "Empirical interatomic potential for Si-H interactions"
			   Phys. Rev. B, 51, 4889 (1995)

  Note:        
*/

#ifndef _MURTY95_H
#define _MURTY95_H

#ifndef NNM /* need to be taken care of later */
#define NNM 60
#endif

#include "mdparallel.h"

#define NGRD 5 /* the number of grid points in F1, F2, H */
class Murty95Frame : public MDPARALLELFrame /* Si with Tersoff potential */
{
    /* general parameters */
    double acut, acutsq;
    double F1_csp[NGRD], F2_csp[NGRD], H_csp[NGRD]; // cubic spline coefficients
	
    /* Muryy potential parameters */
    double A,B,Lam,Mu;  /* Mu is same as Lam2 in the ref paper. */
    double R0,R,D;
    double Alpha,Beta;
    double cTf,dTf,hTf;
	double eta,delta;
	double F1,F2,H;
    
    /* Fitting parameters for Si-Si interaction */
    double A_SiSi,B_SiSi;
    double Lam_SiSi,Mu_SiSi;
	double Alpha_SiSi,Beta_SiSi;
    double R0_SiSi,R_SiSi,D_SiSi;
    double eta_SiSi,delta_SiSi;
    double cTf_SiSi,dTf_SiSi,hTf_SiSi;
    double F1_SiSi,F2_SiSi;

    /* Fitting parameters for H-H interaction */
    double A_HH,B_HH;
    double Lam_HH,Mu_HH;
	double Alpha_HH,Beta_HH;
    double R0_HH,R_HH,D_HH;
    double eta_HH,delta_HH;
    double cTf_HH,dTf_HH,hTf_HH;
    double F1_HH,F2_HH;
	
    /* Fitting parameters for Si-H interaction */
    double A_SiH,B_SiH;
    double Lam_SiH,Mu_SiH;
	double Alpha_SiH,Beta_SiH,Alpha_SiH2,Beta_SiH2;
    double R0_SiH,R_SiH,D_SiH;
    double eta_SiH,delta_SiH;
    double cTf_SiH,dTf_SiH,hTf_SiH[NGRD];
    double cTf_SiH2,dTf_SiH2,hTf_SiH2;
    double F1_SiH[NGRD],F2_SiH[NGRD];
	
public:
    Murty95Frame(){};
    void murty();
    virtual void potential();
    virtual void initvars();
};

#endif // _MURTY95_H

