/*
  fs.h
  by Wei Cai  caiwei@mit.edu
  Last Modified : Sat Mar  1 06:36:44 2008

  Modified for the new version of MD++: Thu Sep 29 2005
  Keonwook Kang kwkang@stanford.edu

  FUNCTION  :  MD simulation package using Finnis-Sinclair potential
*/

#ifndef _FS_H
#define _FS_H

#include "mdparallel.h"

class FSParam : public SCParser 
{ /* Finnis-Sinclair Potential class (Mo, Ta, W) */
public:
    double const_d, const_A, const_beta, const_c,
        const_c0, const_c1, const_c2, const_d_sq, const_c_sq,
        rcut, const_B, const_alpha, const_b0;
    double rlist,skin;
public:
    FSParam(){initparser();};
    void initparser();
    int readfile(char *ppmfile);
    virtual int exec(const char *name);
    virtual int assignvar(int offset=0);
    void print();
};


class FSFrame : public MDPARALLELFrame 
{ /* Finnis-Sinclair (Mo,Ta,W,Fe) */
public:
    class FSParam _FSPM;

    double *rhoh, *af;

    FSFrame(){};

    virtual void potential();
    virtual void potential_energyonly();

    void finnis_sinclair();
    void finnis_sinclair_energyonly();

    virtual void Alloc();

    int readpot();
    virtual void initvars();
    virtual void initparser();
    virtual int exec(const char *name);
#ifdef _PARALLEL
    void Broadcast_FS_Param();
#endif

};    


#endif // _FS_H

