/*
  wire.h
  by Xiaohan Zhang, Wei Cai, caiwei@stanford.edu
  Last Modified : Sun Mar 29 12:07 2016

  FUNCTION  : Couple fiber model with hyperelastic FEM 
  NOTE      : xfem and ljbond are both virtually inherited from mdparallel
              which solves the ``diamond'' problem in multiple-inheritance 
 */

#ifndef _WIRE_H
#define _WIRE_H

#include <cassert>

#include "xfem.h"
#include "ljbond3.h"
#include "linalg3.h"

class WIREFrame : public XFEMFrame, public LJBONDFrame
{
 protected:
  int *atom_fem_links;
  double *atom_fem_weights;
  int _NANCHOR;
  int beads_are_docked;
  double _K_ANCHOR;
 public:
 WIREFrame() : _NANCHOR(0), _K_ANCHOR(0), beads_are_docked(0) {}

  ~WIREFrame() 
    {
      if (cached_EPOT_IND != NULL)   free(cached_EPOT_IND);
      if (cached_EPOT_RMV != NULL)   free(cached_EPOT_RMV);
      if (cached_F != NULL)          free(cached_F);
      if (cached_VIRIAL_IND != NULL) free(cached_VIRIAL_IND);
    }

  void fiber_fem_interaction();
  virtual void potential();

  virtual void plot();
  virtual void Alloc();
  virtual void initvars();
  virtual void initparser();
  virtual int  exec(const char*nam);

 private:
  double cached_EPOT, *cached_EPOT_IND, *cached_EPOT_RMV;
  class Vector3 *cached_F;
  class Matrix33 *cached_VIRIAL_IND;
  class Matrix33 cached_VIRIAL;

  // private helper functions

  void clcvars();
  void clvars();
  void swapvars();
  void accumvars(int );

  void dockbeads();
  bool select_bead(int, int);
  int belongs_to(int);
  int get_weight(int,double*);

  virtual void makefibers();
  virtual void linkfibers();

};

#endif // _WIRE_H
