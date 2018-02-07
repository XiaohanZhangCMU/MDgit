/*
  lj2.h
  by Wei Cai  caiwei@mit.edu
  Last Modified : Fri Feb 16 09:59:52 2007

  FUNCTION  :  Two component Lennard-Jones potential, depending on group ID
*/

#ifndef _LJ2_H
#define _LJ2_H

#include "mdparallel.h"

//All physical constants starts with P_
#define P_KB 1.380662e-23  // (kg*m^2/s^2/K) Boltzmann constant
#define P_E 1.6021892e-19  // (C) electron charge

#define LJ_ENERGY (119.8*P_KB/P_E) // argon energy (sigma)
#define LJ_LENGTH 3.405            // argon length (epsilon)
#define LJ_RC     (2.37343077641*LJ_LENGTH) //4 neighbor

class LJ2Frame : public MDPARALLELFrame /* Si with Stillinger-Weber potential */
{
    /* Lennard-Jones parameters */
    double _ALJ_00,_BLJ_00, ALJ_00,BLJ_00, Uc_00,DUDRc_00;
    double _ALJ_11,_BLJ_11, ALJ_11,BLJ_11, Uc_11,DUDRc_11;
    double _ALJ_01,_BLJ_01, ALJ_01,BLJ_01, Uc_01,DUDRc_01;
    
public:
    LJ2Frame(){};
    void lennard_jones_2();
    virtual void potential();
    virtual void initvars();

    virtual void initparser();
    virtual int exec(const char *name);
    virtual void calcprop();
    
    void initLJ();
    
};



#endif // _LJ2_H

