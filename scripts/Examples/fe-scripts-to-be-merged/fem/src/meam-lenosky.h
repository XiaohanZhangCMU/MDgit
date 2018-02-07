/*
  meam-lenosky.h
  by Wei Cai  caiwei@stanford.edu
  Last Modified : Sat Dec 30 21:18:55 2006

  FUNCTION  : MD++ with Lenosky potential for Si

  Note: The Lenosky potential is for Si only
        Only rectangular supercell is allowed
        Virial stress is not computed
*/

#ifndef _MEAM_LENOSKY_H 
#define _MEAM_LENOSKY_H 

#include "mdparallel.h"

extern "C" int lenosky_(int *nat, double *alat, double *rxyz0,
                        double *fxyz, double *ener_ind, double *ener,
                        double *coord, double *ener_var, double *coord_var, int *count);

class MEAMFrame : public MDPARALLELFrame
{ /* MEAM potential */
public:
    int count;
    MEAMFrame():count(0)  {};

    virtual void potential ();

    void MEAM();

    virtual void Alloc();
    virtual void initvars();
    virtual void initparser();
    
    virtual int exec(const char *nam);

};


#endif // _MEAM_LENOSKY_H

