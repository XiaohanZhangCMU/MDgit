#ifndef _XFEM_H
#define _XFEM_H
#include <cassert>

#include "fem.h"

class XFEMFrame : public virtual FEMFrame /* FEM with brick elements */
{
public:
  char fem_bdy_nodes_file[MAXFILENAMELEN];

  int *bNds, *bTags;
  double *bXs, *bYs, *bZs;
  int n_bdy_nodes;
  int n_bdy_tags;
  int n_existed_atoms;
  int nfemNodes;

  //assert this class is created only after fiber is created
 XFEMFrame():n_bdy_nodes(0),n_bdy_tags(0), n_existed_atoms(_NP), nfemNodes(0) { };

  int read_bdy_nodes(const char*);
  //int read_elements(const char*);
  virtual void initparser();
  virtual int exec(const char*);
  int read2cn();
  void shift_fem_node_id(int);
};

#endif // _XFEM_H

