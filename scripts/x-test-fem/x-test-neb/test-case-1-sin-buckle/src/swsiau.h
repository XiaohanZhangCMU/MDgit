/*
  swsiau.h
  by Wei Cai  caiwei@mit.edu
  Last Modified : Sat Dec 30 21:20:07 2006

  FUNCTION  :  MD simulation package of Ge using Stillinger-Weber potential
*/

#ifndef _SWSIAU_H
#define _SWSIAU_H

#include "mdparallel.h"

class SWFrame : public MDPARALLELFrame /* Ge with Stillinger-Weber potential */
{
    /* Stillinger-Weber potential parameters */
    double psig_si, pepsi_si, aa_si, bb_si, plam_si, pgam_si, acut_si, pss_si;
    double psig_au, pepsi_au, aa_au, bb_au, plam_au, pgam_au, acut_au, pss_au;
    double rho, rho1, acut, acutsq;

public:
    SWFrame(){};
    void stillinger_weber();
    //    void stillinger_weber_energyonly();
    //    double stillinger_weber_energyonly(int iatom);
    virtual void potential();
    //virtual void potential_energyonly();
    //virtual double potential_energyonly(int iatom);
    virtual void initvars();
};

#endif // _SWSIAU_H

