/*
  swlj.h
  by Wei Cai  caiwei@mit.edu
  Last Modified : Sat Jul 28 15:17:09 2007

  FUNCTION  :  MD simulation package of Si using Stillinger-Weber potential
               together with Lennard-Jones potential, depending on group ID
               group 0 and group 0 interact with SW potential
               group 1 and group 1 interact with SW potential
               group 0 and group 1 interact with LJ potential

               This model is developed for Nano-indentation simulation.
*/

#ifndef _SWLJ_H
#define _SWLJ_H

#include "mdparallel.h"

#define LJ_ENERGY 2.314999998 //silicon energy
#define LJ_LENGTH 2.0951      //silicon length sigma
#define LJ_RC     3.77118     //cuttoff for LJ

class SWFrame : public MDPARALLELFrame /* Si with Stillinger-Weber potential */
{
    /* Stillinger-Weber potential parameters */
    double psig, pepsi, aa, bb, plam, pgam, acut, pss, rho;
    double rho1, acutsq;

    /* Lennard-Jones parameters */
    double _ALJ,_BLJ; /* relative depth of the energy well of surface interaction  */
    double  ALJ,BLJ;  /* relative depth of the energy well of surface interaction  */
    double  Uc,DUDRc;

    Vector3 tipforce;
    
public:
    SWFrame(){};
    void stillinger_weber();
    void stillinger_weber_energyonly();
    double stillinger_weber_energyonly(int iatom);
    virtual void potential();
    virtual void potential_energyonly();
    virtual double potential_energyonly(int iatom);
    virtual void initvars();

    virtual void initparser();
    virtual int exec(const char *name);
    virtual void calcprop();
    
    void initLJ();
    
};



#endif // _SWLJ_H

