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
#define BONDINFOSIZE 4
//<<<<<<< .mine
class LJBONDFrame : public virtual MDPARALLELFrame /* Si with Stillinger-Weber potential */
//=======
//class LJBONDFrame : public MDPARALLELFrame /* Si with Stillinger-Weber potential */
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
    double    BOND_CGSD[5];      /* effective interaction through polymer chains */
    double    DENSITY_PARAM[5];  /* parameter for local density potential */
    double    RING_AREA0, RING_AREA_B; /* ring area potential parameters */
    int       RING_CHAINLEN; 
    int       _ENABLE_2D;
    
    int *num_bonds;              /* number of bonds per atom */
    int *bond_index;             /* id of atoms bonded to each atom */

    int total_num_bonds;         /* total number of bonds */
    int *bond_info;              /* i, j, bond_type, n_ij */
     
    char usrfile[200];           /* atomeye usr file name */
    char bndfile[200];           /* file containing bond information */

    /* substrate parameters */
    double _SUBSTRATE_Y, _SUBSTRATE_REP, _SUBSTRATE_ATR, _SUBSTRATE_CON; /* parameter for smooth substrate potential */

    /* wall (channel) parameters */
    double _WALL_K, _WALL_SLOPE, _WALL_CHANNEL_IN, _WALL_CHANNEL_OUT;

    LJBONDFrame(): LJ_RC(0),BOND_R0(0),BOND_K(0),BOND_B(0),RING_AREA0(0),RING_AREA_B(0),RING_CHAINLEN(0),_ENABLE_2D(0),
                   num_bonds(0), bond_index(0),
                   _SUBSTRATE_Y(0),_SUBSTRATE_REP(0),_SUBSTRATE_ATR(0),_SUBSTRATE_CON(0),
                   _WALL_K(0),_WALL_SLOPE(0),_WALL_CHANNEL_IN(0),_WALL_CHANNEL_OUT(0) {};

    void lennard_jones_bond();
    void substrate_potential();
    void wall_potential();
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

    void Alloc_Bonds();
    int read_bonds(const char *name);
    int write_bndfile_from_pairdist(const char *name);

    void calc_local_density();
    
    void writeatomeyeusrfile(char *fname);
};



#endif // _LJBOND_H

