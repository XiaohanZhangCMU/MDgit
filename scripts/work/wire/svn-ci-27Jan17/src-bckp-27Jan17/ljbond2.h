/*
  ljbond.h
  by Wei Cai  caiwei@mit.edu
  Last Modified : Fri Feb 16 09:59:52 2007

  FUNCTION  :  Two component Lennard-Jones potential, depending on group ID
*/

#ifndef _LJBOND_H
#define _LJBOND_H

#include "mdparallel.h"

//All physical constants starts with P_
#define P_KB 1.380662e-23  // (kg*m^2/s^2/K) Boltzmann constant
#define P_E 1.6021892e-19  // (C) electron charge

#define LJ_ENERGY (240.554*P_KB/P_E) //(1.347*P_KB/P_E) // argon energy (sigma)
#define LJ_LENGTH 3.3333   // argon length (epsilon)
//#define LJ_RC     (2.37343077641*LJ_LENGTH) //4 neighbor

/* In principle, MAXSP can be the same as MAXSPECIES defined in md.h
 * here we make MAXSP smaller than MAXSPECIES to save memory 
 */
#define MAXSP 4
#define MAXNUMBOND 3
class LJBONDFrame : public MDPARALLELFrame /* Si with Stillinger-Weber potential */
{
public:
    /* Lennard-Jones parameters */
    double   _ALJ[MAXSP][MAXSP];
    double   _BLJ[MAXSP][MAXSP];
    double    ALJ[MAXSP][MAXSP];
    double    BLJ[MAXSP][MAXSP];
    double     Uc[MAXSP][MAXSP];
    double  DUDRc[MAXSP][MAXSP];
    double    BOND_R0,BOND_K,LJ_RC;
    double  ATTRAC2_Y0,ATTRAC2_K;
    double  WALL3_X0MAX, WALL3_X0, WALL3_K, EWALL3, TNPT_WALL3;
    double  LMAX,XMAX,YMAX,ZMAX,LAMBDA;
    double  dUdLAM_ATT2,dUdLAM_WALL3,dUdLAM_L,dUdLAM_ALKANE,dUdLAM,dUdLAM_XYZ,dUdLAM_X,dUdLAM_Y,dUdLAM_Z;
    
    int *num_bonds;
    int *bond_index;
     
    LJBONDFrame(): LJ_RC(0),BOND_R0(0),BOND_K(0),WALL3_X0MAX(0),WALL3_X0(0),WALL3_K(0),ATTRAC2_Y0(0),ATTRAC2_K(0),LMAX(0),XMAX(0),
                   YMAX(0),ZMAX(0),LAMBDA(0),TNPT_WALL3(0),dUdLAM_ATT2(0),dUdLAM_WALL3(0),dUdLAM_X(0),dUdLAM_Y(0),
                   dUdLAM_Z(0),dUdLAM_XYZ(0),EWALL3(0),dUdLAM_ALKANE(0),dUdLAM_L(0),dUdLAM(0){};
    void lennard_jones_bond();
    virtual void potential();
    virtual void initvars();

    virtual void initparser();
    virtual int exec(const char *name);
    virtual void calcprop();
    
    virtual void Alloc();
    virtual void SWITCHpotential_user(double lambda); /* user defined function for switching */

    void initLJ();
    void Alloc1();
    void makelipids();
    void linklipids();
    void makealkanes();
    void linkalkanes();
    
};



#endif // _LJBOND_H

