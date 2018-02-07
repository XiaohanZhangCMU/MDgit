/*
  meam-lammps.h
  by Wei Cai  caiwei@stanford.edu
  Last Modified : Fri Apr 13 14:25:50 2007

  FUNCTION  : MD++ with MEAM potential (using part of LAMMPS Fortran codes)
*/

#ifndef _MEAM_LAMMPS_H
#define _MEAM_LAMMPS_H 

#include "mdparallel.h"

class Vector6
{
public:
    double x[6];
    Vector6() {for(int i=0;i<6;i++) x[i]=0;}
    double operator [] (int i) const { return x[i]; }
    double & operator [] (int i) { return x[i]; }
    void clear() {for(int i=0;i<6;i++) x[i]=0;}
};

class Vector10
{
public:
    double x[10];
    Vector10() {for(int i=0;i<10;i++) x[i]=0;}
    double operator [] (int i) const { return x[i]; }
    double & operator [] (int i) { return x[i]; }
    void clear() {for(int i=0;i<10;i++) x[i]=0;}
};

/*  LAMMPS pair_meam.h */

/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef _CYGWIN
#define meam_setup_global_ meam_setup_global__
#define meam_setup_param_ meam_setup_param__
#define meam_setup_done_ meam_setup_done__
#define meam_dens_init_ meam_dens_init__
#define meam_dens_final_ meam_dens_final__
#define meam_force_ meam_force__
//#define phif_ phif__
#endif

extern "C" {
  void meam_setup_global_(int *, int *, double *, int *, double *, double *,
			 double *, double *, double *, double *, double *,
			 double *, double *, double *, double *, double *,
			 double *, double *, int *);
  void meam_setup_param_(int *, double *, int *, int *, int *);
  void meam_setup_done_(double *);

  void meam_dens_init_(int *, int *, int *, int *, int *, double *,
		       int *, int *, int *, int *,
		       double *, double *, double *, double *, double *, double *,
		       double *, double *, double *, double *, double *, int *);
  
  void meam_dens_final_(int *, int *, int *, int *, int *, 
                        double *, double *, int *, int *, int *,
			double *, double *, double *, double *, double *, double *,
			double *, double *, double *, double *, double *, double *, 
			double *, double *, double *, double *, double *,
                        int *, int *);

#ifdef _TORSION_OR_BENDING 
  void meam_force_(int *, int *, int *, int *, double *, int *, int *, int *, double *,
		   int *, int *, int *, int *, double *, double *, double *,
		   double *, double *, double *, double *, double *, double *, double *, double *,
		   double *, double *, double *, double *, double *, double *, double *, double *,
                   int *, double *, 
                   int *, double *,
                   int *, double *, double *,
                   int *);
#else
  void meam_force_(int *, int *, int *, int *, int *,
                   int *, double *, double *, int *, int *, int *, double *,
		   int *, int *, int *, int *,
		   double *, double *, double *, double *, double *, double *,
		   double *, double *, double *, double *, double *, double *,
		   double *, double *, double *, double *, double *, double *,
                   double *, double *, int *);

#endif

  void phif_(int *, int *, double *, double *, double *);
}

enum{FCC,BCC,HCP,DIM,DIAMOND,B1,C11,L12};
int nkeywords = 21;
char *keywords[] = {"Ec","alpha","rho0","delta","lattce",
                    "attrac","repuls","nn2","Cmin","Cmax","rc","delr",
                    "augt1","gsmooth_factor","re","ialloy","mix_ref_t","erose_form",
                    "zbl","emb_lin_neg","bkgd_dyn"};

class MEAMFrame : public MDPARALLELFrame
{ /* MEAM potential */
 private:

    char meamfile[1000], meafile[1000]; /* potential file */
    
//  double cutmax;                // max cutoff for all elements
//  int nelements;                // # of unique elements
//  char **elements;              // names of unique elements
//  double *mass;                 // mass of each element
//  int reverse_flag;             // which pass of reverse comm is being done
//
//  int *map;                     // mapping from atom types to elements
//  int *fmap;                    // Fortran version of map array for MEAM lib
    
    double *scrfcn,*dscrfcn,*fcpair;
    double *rho,*rho0,*rho1,*rho2,*rho3,*frhop;
    double *gamma,*dgamma1,*dgamma2,*dgamma3,*arho2b;
    Vector6 *vatom;
    Vector3 *arho1; Vector6 *arho2; Vector10 *arho3;
    Vector3 *arho3b, *t_ave, *tsq_ave, *rtmp;
    int *type, *fmap;
    double rcut;

    /* neighbor list */
    int maxneigh;
    int *num_neigh_full, **ind_neigh_full; char *ind_neigh_full_mem;
    int *num_neigh_half, **ind_neigh_half; char *ind_neigh_half_mem;
    
//  void allocate();
//  void read_files(char *, char *);

  void neigh_f2c(int *numn, int **firstn);
  void neigh_c2f(int *numn, int **firstn);

public:
    
    MEAMFrame():
                scrfcn(0),dscrfcn(0),fcpair(0),
                rho(0),rho0(0),rho1(0),rho2(0),rho3(0),frhop(0),
                gamma(0),dgamma1(0),dgamma2(0),dgamma3(0),arho2b(0),
                arho1(0),arho2(0),arho3(0),arho3b(0),t_ave(0),tsq_ave(0),
                rtmp(0),fmap(0),rcut(0), maxneigh(0),
                num_neigh_full(0),ind_neigh_full(0),ind_neigh_full_mem(0),
                num_neigh_half(0),ind_neigh_half(0),ind_neigh_half_mem(0)
    {};

    virtual void potential ();
    virtual void NbrList_reconstruct(int iatom=-1);
    void NbrList_translate();
    
    int readMEAM();
    void read_files(char *, char *);
    
    void MEAM();

    virtual void Alloc();
    virtual void initvars();
    virtual void initparser();
    
    virtual int exec(const char *nam);

    void compute(int, int);
    void settings(int, char **);
    void coeff(int, char **);
    double init_one(int, int);
    void init_style();

    void printpairpot();
    
#ifdef _PARALLEL
    void Broadcast_MEAM_Param();
#endif

//  int pack_comm(int, int *, double *, int *);
//  void unpack_comm(int, int, double *);
//  int pack_reverse_comm(int, int, double *);
//  void unpack_reverse_comm(int, int *, double *);
//  int memory_usage();


};

#endif // _MEAM-LAMMPS_H

