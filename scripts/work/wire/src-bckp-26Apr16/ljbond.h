/*
  ljbond.h
  by Wei Cai  caiwei@stanford.edu
  Last Modified : Fri Feb 16 09:59:52 2007

  FUNCTION  :  Two component Lennard-Jones potential, depending on group ID
*/

#ifndef _LJBOND_H
#define _LJBOND_H

#include "mdparallel.h"

//All physical constants starts with P_
#define P_KB 1.380662e-23  // (kg*m^2/s^2/K) Boltzmann constant
#define P_E 1.6021892e-19  // (C) electron charge

#define LJ_ENERGY (240.554*P_KB/P_E) // argon energy (sigma)
#define LJ_LENGTH 3.3333   // argon length (epsilon)
//#define LJ_RC     (2.37343077641*LJ_LENGTH) //4 neighbor

/* let MAXSP be the same as MAXSPECIES defined in md.h */
#define MAXSP MAXSPECIES
#define MAXNUMBOND 3
class LJBONDFrame : public MDPARALLELFrame /* Si with Stillinger-Weber potential */
{
public:
    /* Lennard-Jones parameters */
    double   _C12[MAXSP][MAXSP];
    double   _C9 [MAXSP][MAXSP];
    double   _C6 [MAXSP][MAXSP];
    double    C12[MAXSP][MAXSP];
    double    C9 [MAXSP][MAXSP];
    double    C6 [MAXSP][MAXSP];
    double     Uc[MAXSP][MAXSP];
    double  DUDRc[MAXSP][MAXSP];
    double    LJ_RC;
    double  BOND_R0, BOND_K;

    double LX0, LX1, LY0, ALKANE_SY0;
    double ATTRAC_EPOT, ATTRAC_W, ATTRAC_K; 
    double PROBE1_EPOT, PROBE1_WX, PROBE1_WY, PROBE1_DY, PROBE1_SY0, PROBE1_K;
    double PROBE2_EPOT, PROBE2_WX, PROBE2_WY, PROBE2_DY, PROBE2_SY0, PROBE2_K, PROBE2_KP;
    double PROBE3_EPOT, PROBE3_WX0,PROBE3_WX1,PROBE3_WY, PROBE3_DY,  PROBE3_SY0,PROBE3_K;
    double PROBE4_EPOT, PROBE4_WX0,PROBE4_WX1,PROBE4_WY, PROBE4_DY,  PROBE4_SY0,PROBE4_K;
    double PROBEC_EPOT, PROBEC_WX, PROBEC_WY, PROBEC_DY, PROBEC_SY0, PROBEC_K;

    double LJ_EPOT, BOND_EPOT;
    double LJ_FACTOR, BOND_FACTOR, ALKANE_FACTOR, LX_FACTOR;
    double ATTRAC_K_FACTOR, PROBE1_WX_FACTOR, PROBE2_WX_FACTOR, PROBE3_WX_FACTOR, PROBE4_WX_FACTOR, PROBEC_WX_FACTOR;
    double dEdLJ_FACTOR, dEdBOND_FACTOR, dEdALKANE_FACTOR, dEdLX_FACTOR, dEdATTRAC_K_FACTOR;
    double dEdPROBE1_WX_FACTOR, dEdPROBE2_WX_FACTOR, dEdPROBE3_WX_FACTOR, dEdPROBE4_WX_FACTOR, dEdPROBEC_WX_FACTOR;
    char   LJ_FACTOR_SPEC[20], BOND_FACTOR_SPEC[20], ALKANE_FACTOR_SPEC[20];
    char   LX_FACTOR_SPEC[20], ATTRAC_K_FACTOR_SPEC[20];
    char   PROBE1_WX_FACTOR_SPEC[20], PROBE2_WX_FACTOR_SPEC[20];
    char   PROBE3_WX_FACTOR_SPEC[20], PROBE4_WX_FACTOR_SPEC[20], PROBEC_WX_FACTOR_SPEC[20];

    int *num_bonds;
    int *bond_index;
    double lipid_shape_spec[10];
    int LIPID_CHAINLEN;
 
    char usrfile[200];           /* atomeye usr file name */

    LJBONDFrame(): LJ_RC(0),
                   BOND_R0(0),BOND_K(0), LX0(10.0), LX1(10.0), LY0(10.0), ALKANE_SY0(0.0),

                   ATTRAC_EPOT(0), ATTRAC_W(0),  ATTRAC_K(0), 
                   PROBE1_EPOT(0), PROBE1_WX(0), PROBE1_WY(0), PROBE1_DY(1), PROBE1_SY0(0), PROBE1_K(0),
                   PROBE2_EPOT(0), PROBE2_WX(0), PROBE2_WY(0), PROBE2_DY(1), PROBE2_SY0(0), PROBE2_K(0), PROBE2_KP(100),
                   PROBE3_EPOT(0), PROBE3_WX0(0),PROBE3_WX1(0),PROBE3_WY(0), PROBE3_DY(1),  PROBE3_SY0(0),PROBE3_K(0),
                   PROBE4_EPOT(0), PROBE4_WX0(0),PROBE4_WX1(0),PROBE4_WY(0), PROBE4_DY(1),  PROBE4_SY0(0),PROBE4_K(0),
                   PROBEC_EPOT(0), PROBEC_WX(0), PROBEC_WY(0), PROBEC_DY(1), PROBEC_SY0(0), PROBEC_K(0), 

                   LJ_EPOT(0), BOND_EPOT(0),
                   LJ_FACTOR(1.0),  BOND_FACTOR(1.0),  ALKANE_FACTOR(1.0),
                   LX_FACTOR(0.0),  ATTRAC_K_FACTOR(1.0),  
                   PROBE1_WX_FACTOR(1.0),  PROBE2_WX_FACTOR(1.0),  
                   PROBE3_WX_FACTOR(0.0),  PROBE4_WX_FACTOR(0.0), PROBEC_WX_FACTOR(1.0),
                   dEdLJ_FACTOR(0), dEdBOND_FACTOR(0), dEdALKANE_FACTOR(0),
                   dEdLX_FACTOR(0), dEdATTRAC_K_FACTOR(0), 
                   dEdPROBE1_WX_FACTOR(0), dEdPROBE2_WX_FACTOR(0), 
                   dEdPROBE3_WX_FACTOR(0), dEdPROBE4_WX_FACTOR(0), dEdPROBEC_WX_FACTOR(0),
                   LIPID_CHAINLEN(0)

                   {};

    void lennard_jones_bond();
    virtual void potential();
    virtual void initvars();

    virtual void initparser();
    virtual int exec(const char *name);
    virtual void calcprop();
    
    virtual void Alloc();
    virtual void plot();

    void assign_from_lambda(const char *spec, double *value);
    void add_to_dEdlambda  (const char *spec, double value);

    void initLJ();
    void makelipids();
    void linklipids();

    void update_groupID_by_position();
    void update_groupID_by_energy();
    
    void writeatomeyeusrfile(char *fname);

    /* from ljbond2.cpp */
    virtual void SWITCHpotential_user(double lambda);
    void Alloc1();
    void makealkanes();
    void linkalkanes();

    /* addition from Jichul Kim */
    void makeproteins();
    void linkproteins();
};



#endif // _LJBOND_H

