/*
  fem.h
  by Wei Cai  caiwei@stanford.edu
  Last Modified : Fri Mar 11 21:07:00 2016

  FUNCTION  :  FEM model of hyperelastic material with brick elements
  NOTE :       All elements must be the same shape and size for now.
*/

#ifndef _FEM_H
#define _FEM_H
#include <fstream>      // std::ofstream

#include <cassert>
#include "mdparallel.h"

#define _USECUDA
#include "linalg3_cu.h"

#define _MAX_NELEM_SHARE_NODE 10

class FEMFrame : public virtual MDPARALLELFrame /* FEM with brick elements */
{
 protected:
    /* FEM parameters */
    int _NELE;              /* number of elements */
    int _NNODE_PER_ELEMENT; /* number of nodes per element */
    int _NINT_PER_ELEMENT;  /* number of Gauss integration points per element */
    int _NFACE_PER_ELEMENT;
    int _NDIM;              /* dimension of elements */
    int _ReadColorID;
    int _EquationType;
    int *elements;
    int *inv_elements;
    int *colorids;
    int *map23;

#ifdef _USECUDA
    int *_d_elements;
    int *_d_inv_elements;
    double *_d_gauss_weight; 
    int *_d_map23;
    G_Vector3 *_d_SR;
    G_Vector3 *_d_SRref;       
    G_Vector3 *_d_Rref;        
    G_Vector3 *_d_F;
    G_Vector3 *_d_F_padding; //for padding forces of elements
    double *_d_dFdu;        
    G_Matrix33 _d_H;

    //The following device pointers are not yet memcpied. 
    int *_d_colorids;
    int *_d_fixed;
    double *_d_EPOT;
    double *_d_EPOT_IND;
    double* _d_EPOT_RMV;
    class G_Matrix33 *_d_VIRIAL_IND;
#endif

    char ELEMENT_TYPE[100], ELEMENT_INT_TYPE[100];
    char elements_file[MAXFILENAMELEN];
    char fem_coeff_file[MAXFILENAMELEN];

    double y_eigen_strain;
    double x_eigen_strain;
    double x_eigen_zbound_min;
    double x_eigen_zbound_max;
    double y_eigen_zbound_min;
    double y_eigen_zbound_max;

    double *gauss_weight;   /* weight of each Gauss integration point */
    double *dFdu;           /* derivative of F (deformation gradient) wrt nodal displacement */
    Vector3 *_SRref;        /* reference nodal position (scaled coordinates) */
    Vector3 *_Rref;         /* reference nodal position (real coordinates) */

public:
    FEMFrame():_NELE(0),_NNODE_PER_ELEMENT(0),_NINT_PER_ELEMENT(0),
      elements(0),gauss_weight(0),dFdu(0),_SRref(0),_Rref(0), _ReadColorID(0),
      n_bdy_nodes(0),n_bdy_tags(0), n_existed_atoms(_NP), nfemNodes(0) {};

    void RtoRref();
    int read_elements(const char*);
    int read_fem_coeff(const char*);
    int read_uniform_fem_coeff(const char*);
    int read_element_wise_fem_coeff(const char*, int);
    int read_1stfem_to_Alloc(const char*);
    Matrix33 host_getEigenF(Vector3 p, Matrix33 Fdef);

    void WriteStressCoord();

    void test_saxpy();
    void beam_fem_energy_force();
    void create_inverse_connectivities_matrix();
    void islam_fem_energy_force();
    void snap_fem_energy_force();
    void snap_fem_energy_force_1();

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
  char contact_file[MAXFILENAMELEN];

  int *bNds, *bTags;
  int *cNds;
  double *bXs, *bYs, *bZs;
  int n_bdy_nodes;
  int n_bdy_tags;  
  int n_existed_atoms;
  int nfemNodes;

  int n_contact_nodes;

  int read_bdy_nodes(const char*);
  int read_contact(const char*);
  //int read_elements(const char*);
#ifdef _PARALLEL
  void Broadcast_FEM_Param();
#endif

  int read2cn();
  void shift_fem_node_id(int);

#ifdef _USECUDA
  void cuda_memcpy_init_all(void);
  void cuda_memory_alloc(int);
  void cuda_memory_alloc_elements(int);
  void cuda_memory_alloc_element_coeff(int, int);
#endif
};

#endif // _FEM_H

