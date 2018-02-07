/*
  ljbond3.h
  by Wei Cai  caiwei@stanford.edu
  Last Modified : Thu Dec 16 10:05:46 PST 2010

  FUNCTION  :  MD simulation package of Lennard-Jones (LJ) potential
               LJ parameters depend on species
               Atoms can also form covalent bond with neighbors
               covalent bonds have stretching energy and bending energy

               Disabled free energy calculation in this model

  Test cases: scripts/work/fibers/fibers.tcl
*/

#ifndef _LJBOND_H
#define _LJBOND_H

#include "mdparallel.h"

/* All physical constants starts with P_ */
#define P_KB 1.380662e-23            // (kg*m^2/s^2/K) Boltzmann constant
#define P_E 1.6021892e-19            // (C) electron charge

#define LJ_ENERGY (240.554*P_KB/P_E) //(1.347*P_KB/P_E) // argon energy (sigma)
#define LJ_LENGTH 3.3333             // argon length (epsilon)

/* In principle, MAXSP can be the same as MAXSPECIES defined in md.h
 * here we make MAXSP smaller than MAXSPECIES to save memory 
 */
#define MAXSP 4
#define MAXNUMBOND 3
class LJBONDFrame : public virtual MDPARALLELFrame /* Si with Stillinger-Weber potential */
{
public:
    /* Lennard-Jones parameters */
    double   _ALJ[MAXSP][MAXSP];
    double   _BLJ[MAXSP][MAXSP];
    double    ALJ[MAXSP][MAXSP];
    double    BLJ[MAXSP][MAXSP];
    double     Uc[MAXSP][MAXSP];
    double  DUDRc[MAXSP][MAXSP];
    double    LJ_RC;             /* cut-off radius */
    double    BOND_R0,BOND_K;    /* bond stretching stiffness */
    double    BOND_B;            /* bond bending stiffness */
    
    int *num_bonds;
    int *bond_index;
     
    char usrfile[200];           /* atomeye usr file name */

    /* substrate parameters */
    double _SUBSTRATE_Y, _SUBSTRATE_REP, _SUBSTRATE_ATR, _SUBSTRATE_CON; /* parameter for smooth substrate potential */

    LJBONDFrame(): LJ_RC(0),BOND_R0(0),BOND_K(0),BOND_B(0),
                   num_bonds(0), bond_index(0),
                   _SUBSTRATE_Y(0),_SUBSTRATE_REP(0),_SUBSTRATE_ATR(0),_SUBSTRATE_CON(0) {};

    void lennard_jones_bond();
    void substrate_potential();
    virtual void potential();
    virtual void initvars();

    virtual void initparser();
    virtual int exec(const char *name);
    
    virtual void Alloc();
    virtual void plot();

    void initLJ();
    void Alloc1();
    virtual void makefibers();
    virtual void linkfibers();
    
    void writeatomeyeusrfile(char *fname);
};



#endif // _LJBOND_H

