/*
  lj.h
  by Wei Cai  caiwei@mit.edu
  Last Modified : Fri Feb 16 09:59:47 2007

  FUNCTION  :  Argon MD simulation using Lennard-Jones Potential

  To Do     :  test potential integrety
*/

#ifndef _LJ_H
#define _LJ_H

#include "mdparallel.h"

//All physical constants starts with P_
#define P_C 2.99792458e8  // (m/s) the speed of light (EXACT)
#define P_HBAR 1.0545919e-34 // (kg*m^2/s) hbar plank constant
#define P_KB 1.380662e-23  // (kg*m^2/s^2/K) Boltzmann constant
#define P_E 1.6021892e-19  // (C) electron charge
#define P_ME 9.109558e-31  // (kg) electron mass
#define P_G 6.6732e-11      // (m^3/s^2/kg) Gravitational constant
#define P_NA 6.022169e23   // (1/mol) Avogadro constant
#define P_U 1.660531e-27    // (kg) atomic mass unit
#define P_MU0 (4*M_PI*1e-7) // (C^2*s^4/kg/m^5)
#define P_EPSILON0 (1/P_MU0/P_C/P_C) // (kg*m^3/C^2/s^2)
//Equality sqrt(P_MU0*P_EPSILON0)*P_C=1

#define LJ_ENERGY (119.8*P_KB/P_E) //argon energy

#define LJ_LENGTH 3.405 //argon length
#define LJ_MASS 39.948  //argon mass

#define LJ_DT (2e-3/UTIME_IN_PS)
//#define LJ_RC (3.5*LJ_LENGTH) //a0 = 5.2688, NNM = 360
//#define LJ_RC (2.5*LJ_LENGTH) //5 neighbor, a0 = 5.3050, NNM = 120
#define LJ_RC (2.37343077641*LJ_LENGTH) //4 neighbor, a0 = 5.3142, NNM = 120

const double A=4*LJ_ENERGY*POW12(LJ_LENGTH);
const double B=4*LJ_ENERGY*POW6(LJ_LENGTH);
const double Uc=(A*POW6(1/LJ_RC)-B)*POW6(1/LJ_RC);
const double DUDRc=-(12*A/POW6(LJ_RC)-6*B)/POW6(LJ_RC)/LJ_RC;
//const int ArgonFlag=F_NEEDKAPPA;

class LJFrame: public MDPARALLELFrame
{
public:
    virtual void initvars();
    virtual void potential();
    void lennard_jones();
};        



#endif // _LJ_H

