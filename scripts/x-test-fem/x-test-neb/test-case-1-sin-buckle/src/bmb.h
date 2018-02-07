/*
  bmb.h
  by Wei Cai  caiwei@mit.edu
  Last Modified : Mon Aug 27 11:43:03 2007

  FUNCTION  :  BMB potential (for YSZ)
  
  Test cases: scripts/work/ysz/zro2.tcl
*/

#ifndef _BMB_H
#define _BMB_H

#include "mdparallel.h"

//All physical constants starts with P_
#define P_KB 1.380662e-23  // (kg*m^2/s^2/K) Boltzmann constant
#define P_E 1.6021892e-19  // (C) electron charge
#define P_BMB 1.0*1000.0/(EV*AVO) // (eV) conversion factor for BMB parameters (From wikipedia)

/* In principle, MAXSP can be the same as MAXSPECIES defined in md.h
 * here we make MAXSP smaller than MAXSPECIES to save memory 
 */
#define MAXSP 4
#define MAXNUMBOND 3

//-------------------------------------------------------------------
// Interatomic Parameters
//-------------------------------------------------------------------
// BMB Potential Parameters
//-------------------------------------------------------------------
// [ Zr-Zr | Zr-O | O-O | Zr-Y | O-Y | None | Y-Y ] [ A | B | C ]
//-------------------------------------------------------------------
const double BMB_POT[7][3] = {{    0.0*P_BMB,        0.0*P_BMB,        1.0 },  // Zr-Zr
			      {    0.0*P_BMB,    95120.0*P_BMB, 2.65957450 },  // Zr-O
			      { 2691.0*P_BMB,  2196000.0*P_BMB, 6.71140940 },  // O-O
			      {    0.0*P_BMB,        0.0*P_BMB,        1.0 },  // Zr-Y
			      {    0.0*P_BMB,   129780.0*P_BMB, 2.86450870 },  // O-Y
			      {    0.0*P_BMB,        0.0*P_BMB,        1.0 },  // dummy variables
			      {    0.0*P_BMB,        0.0*P_BMB,        1.0 }}; // Y-Y

class BMBFrame : public MDPARALLELFrame /* Born Meyer Buckingham potential (requires Ewald) */
{
public:
    double BMB_Rc;
    
    BMBFrame(){};
    void born_meyer_buckingham();
    void born_meyer_buckingham_Ewald();
    virtual void potential();
    virtual void initvars();

    virtual void initparser();
    virtual int exec(const char *name);
    virtual void calcprop();
    
    virtual void Alloc();

    void dope_Yttria();
};



#endif // _BMB_H

