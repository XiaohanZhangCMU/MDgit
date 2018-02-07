/*
  fem.h
  by Wei Cai  caiwei@stanford.edu
  Last Modified : Fri Mar 11 21:07:00 2016

  FUNCTION  :  FEM model of hyperelastic material with brick elements
  NOTE :       All elements must be the same shape and size for now.
*/

#ifndef _FEM_H
#define _FEM_H

#include <cassert>
#include "mdparallel.h"

class FEMFrame : public virtual MDPARALLELFrame /* FEM with brick elements */
{
 protected:
    /* FEM parameters */
    int _NELE;              /* number of elements */
    int _NNODE_PER_ELEMENT; /* number of nodes per element */
    int _NINT_PER_ELEMENT;  /* number of Gauss integration points per element */
    int _NDIM;              /* dimension of elements */
    int *elements;
    int *map23;
    char ELEMENT_TYPE[100], ELEMENT_INT_TYPE[100];
    char elements_file[MAXFILENAMELEN];
    char fem_coeff_file[MAXFILENAMELEN];

    double *gauss_weight;   /* weight of each Gauss integration point */
    double *dFdu;           /* derivative of F (deformation gradient) wrt nodal displacement */
    Vector3 *_SRref;        /* reference nodal position (scaled coordinates) */
    Vector3 *_Rref;         /* reference nodal position (real coordinates) */

public:
    FEMFrame():_NELE(0),_NNODE_PER_ELEMENT(0),_NINT_PER_ELEMENT(0),
      elements(0),gauss_weight(0),dFdu(0),_SRref(0),_Rref(0),
      n_bdy_nodes(0),n_bdy_tags(0), n_existed_atoms(_NP), nfemNodes(0) {};

    void RtoRref();
    int read_elements(const char*);
    int read_fem_coeff(const char*);

    void fem_energy_force();
    virtual void potential();

    //void fem_energy_only();
    //virtual void potential_energyonly();

    void Alloc_Elements();
    void Alloc_Element_Coeff();
    
    virtual void plot();
    virtual void Alloc();

    virtual void initvars();
    virtual void initparser();
    virtual int exec(const char *nam);

    // The following is from xfem
public:
  char fem_bdy_nodes_file[MAXFILENAMELEN];

  int *bNds, *bTags;
  double *bXs, *bYs, *bZs;
  int n_bdy_nodes;
  int n_bdy_tags;
  int n_existed_atoms;
  int nfemNodes;

  virtual int read_bdy_nodes(const char*);
  //int read_elements(const char*);
#ifdef _PARALLEL
  void Broadcast_FEM_Param();
#endif

  int read2cn();
  void shift_fem_node_id(int);
};

#endif // _FEM_H

