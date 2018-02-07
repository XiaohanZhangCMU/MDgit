/*
  ljdimer.h
  by Wei Cai  caiwei@mit.edu
  Last Modified : Fri Feb 23 18:31:19 2007

  FUNCTION  :  Lennard-Jones potential and dimer
*/

#ifndef _LJDIMER_H
#define _LJDIMER_H

#include "mdparallel.h"

//All physical constants starts with P_
#define P_KB 1.380662e-23  // (kg*m^2/s^2/K) Boltzmann constant
#define P_E 1.6021892e-19  // (C) electron charge

#define LJ_ENERGY (119.8*P_KB/P_E) // argon energy (epsilon)
#define LJ_LENGTH 3.405            // argon length (sigma)
#define LJ_RC     (1.12246204830937*LJ_LENGTH) //2^(1/6) *(sigma)

class LJDIMERFrame : public MDPARALLELFrame /* LJ and dimer */
{
    /* Lennard-Jones parameters */
    double _ALJ_00,_BLJ_00, ALJ_00,BLJ_00, Uc_00,DUDRc_00;

    /* Dimer parameter (JCP, 110, 6617, 1999) */
    double _H_DIMER, _W_DIMER, H_DIMER, W_DIMER, R_DIMER, F_DIMER,
           F_DIMER_INTERNAL, F_DIMER_EXTERNAL;
    
public:
    LJDIMERFrame():_ALJ_00(4),_BLJ_00(4),
                   _H_DIMER(6),_W_DIMER(0.25), R_DIMER(0), 
                   F_DIMER(0), F_DIMER_INTERNAL(0), F_DIMER_EXTERNAL(0){};
    void lennard_jones_dimer();
    void find_dimer_indices(int *ind0,int *ind1);
    virtual void potential();
    virtual void SWITCHpotential_user(double lambda);
    void lennard_jones_dimer_constrained(double R);
    virtual void initvars();

    virtual void initparser();
    virtual int exec(const char *name);
    virtual void calcprop();
    
    void initLJ();
    
};



#endif // _LJDIMER_H

