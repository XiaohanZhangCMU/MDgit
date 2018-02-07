/*
  rebolj.h
  by Wei Cai  caiwei@stanford.edu
  Last Modified : Aug 9 2011

  FUNCTION  :  Add LJ potential on top of REBO potential

*/

#ifndef _REBOLJ_H
#define _REBOLJ_H

#include "rebo.h"

//All physical constants starts with P_
#define P_KB 1.380662e-23  // (kg*m^2/s^2/K) Boltzmann constant
#define P_E 1.6021892e-19  // (C) electron charge

#define LJ_ENERGY (240.554*P_KB/P_E) // argon energy (sigma)
#define LJ_LENGTH 3.3333   // argon length (epsilon)
//#define LJ_RC     (2.37343077641*LJ_LENGTH) //4 neighbor

/* let MAXSP be the same as MAXSPECIES defined in md.h */
#define MAXSP MAXSPECIES

class REBOLJFrame : public REBOFrame /* Carbon, H with Brenner potential */
{
    /* Lennard-Jones parameters (see ljbond.h) */
    double   _ALJ[MAXSP][MAXSP];
    double   _BLJ[MAXSP][MAXSP];
    double    ALJ[MAXSP][MAXSP];
    double    BLJ[MAXSP][MAXSP];
    double     Uc[MAXSP][MAXSP];
    double  DUDRc[MAXSP][MAXSP];
    double    LJ_RC;

    int _NUM_REBO_GROUPS;

    /* substrate parameters */
    double _SUBSTRATE_Z, _SUBSTRATE_REP, _SUBSTRATE_ATR; /* parameter for smooth substrate potential */

public:
    REBOLJFrame():LJ_RC(0),_NUM_REBO_GROUPS(0),_SUBSTRATE_Z(0),_SUBSTRATE_REP(0),_SUBSTRATE_ATR(0) {};

    void lennard_jones();
    void substrate_potential();
    void initLJ();

    virtual void potential();
    virtual void initvars();
    virtual void initparser();
    
    virtual int exec(const char *nam);

};

#endif // _REBOLJ_H

