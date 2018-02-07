/*
  rods.h
  by Wei Cai  caiwei@mit.edu
  Last Modified : Wed May  2 18:37:42 2007

  FUNCTION  :  nodes connected by a set of rods
*/

#ifndef _RODS_H
#define _RODS_H

#include "mdparallel.h"

#define MAXSP 4
#define MAXNUMBOND 3
class RODSFrame : public MDPARALLELFrame /* Si with Stillinger-Weber potential */
{
public:
    double *BOND_R0, BOND_K, REPUL;
    int *BOND_INDEX;
    int NUM_BONDS;
 
    RODSFrame():BOND_R0(0),BOND_K(1),REPUL(0),BOND_INDEX(0),NUM_BONDS(0){};

    void rods_potential();
    virtual void potential();
    virtual void initvars();

    virtual void initparser();
    virtual int exec(const char *name);
    virtual void calcprop();
    
    virtual void Alloc();
    virtual void plot();
    //virtual void writeatomeyecfgfile(char *fname);

    void init_structure();
    void final_structure();
    void connect_nodes();
    void Alloc_Bonds();
    void zero_com_rotation();
};


#endif // _RODS_H

