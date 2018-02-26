/*
  eam.h
  by Wei Cai  caiwei@mit.edu, Xiaohan zhang (GPU)
  Last Modified : Mon Sep  8 18:41:56 2008

  FUNCTION  : MD++ implementation of EAM/MEAM potential
*/

#ifndef _EAM_H
#define _EAM_H

#include <assert.h>
#include "mdparallel.h"

#define _USECUDA
#ifdef _USECUDA
#include <cuda_runtime.h>
#include "linalg3_cu.h"
//#define DEBUG_USECUDA
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code,const char *file,int line,bool abort = true) {
  if (code != cudaSuccess) {
    fprintf(stderr, "GPUassert: %s %s %dn",cudaGetErrorString(code),file,line);
    if (abort) { getchar(); exit(code); } } } 
#define _MAX_NELEM_SHARE_NODE 10
#endif


class EAMFrame : public MDPARALLELFrame
{ /* EAM/MEAM potential */
public:
    int pottype; /* 1: EAM, 2: MEAM */
    int eamfiletype; /* 1: single element, 2: binary element */
    int eamgrid;
    /* EAM potential parameters */
#define NGRID 6000
    char title3[500];
    int ntyp;
    double rmass, rlatt, drar, drhoar, actual, actual2, rmin;
    double rho[2][NGRID], rhop[2][NGRID], phi[2][NGRID], phip[2][NGRID];
    double phix[NGRID], phipx[NGRID];
    double frho[2][NGRID], frhop[2][NGRID];
    double rho_spline[2][NGRID][4], phi_spline[2][NGRID][4];
    double phix_spline[NGRID][4], frho_spline[2][NGRID][4];
    double rval[NGRID], rhoval[NGRID];
    double petrip, rhocon, rhomax;
    double *atpe3b, *rhotot, *embf, *embfp;
    double *atpe3b1, *rhotot1, *embf1;
    double *atpe3bUMB, *rhototUMB, *embfUMB;
    double MC_maxd_trial;
    int MC_need, *nbst, *nbst1;

#ifdef _USECUDA
    int *_d_nbst;
    int *_d_nn;
    int *_d_nindex;
    int *_d_fixed;
    int *_d_species;

    double *_d_rho, *_d_rhop, *_d_phi, *_d_phip;
    double *_d_phix, *_d_phipx;
    double *_d_frho, *_d_frhop;
    double *_d_rho_spline, *_d_phi_spline;
    double *_d_phix_spline, *_d_frho_spline;
    double *_d_rval, *_d_rhoval;
    double *_d_atpe3b, *_d_rhotot, *_d_embf, *_d_embfp;

    double *_d_H_element;
    double *_d_VIRIAL_IND_element;
    double *_d_EPOT_IND;
    G_Vector3 *_d_SR;
    G_Vector3 *_d_F;
    /* order:
     *  0     1     2     3      4      5      6     7      8      9 
     * rmass,rlatt,drar,drhoar,actual,actual2,rmin,petrip,rhocon,rhomax
     */
    double *fscalars, *_d_fscalars;  
#endif

#ifdef DEBUG_USECUDA
    double _h_d_rho[2][NGRID], _h_d_rhop[2][NGRID], _h_d_phi[2][NGRID], _h_d_phip[2][NGRID];
    double _h_d_phix[NGRID], _h_d_phipx[NGRID];
    double _h_d_frho[2][NGRID], _h_d_frhop[2][NGRID];
    double _h_d_rho_spline[2][NGRID][4], _h_d_phi_spline[2][NGRID][4];
    double _h_d_phix_spline[NGRID][4], _h_d_frho_spline[2][NGRID][4];
    double _h_d_rval[NGRID], _h_d_rhoval[NGRID];

    double *_h_d_atpe3b, *_h_d_rhotot, *_h_d_embf, *_h_d_embfp;
    int *_h_d_nbst;
    Matrix33 _h_d_H;
    Matrix33 _h_d_VIRIAL;
    Matrix33 *_h_d_VIRIAL_IND;
    double *_h_d_EPOT_IND;
    int *_h_d_species;
    int **_h_d_nindex;
    int *_h_d_nn;
    Vector3 *_h_d_SR;
    Vector3 *_h_d_F;
    double *_h_d_fscalars;
    int check_host_device_memory_transfer(); 
#endif
   
    int NumOfChanged,*ListOfChanged;
    
    EAMFrame() : pottype(1), eamfiletype(0), eamgrid(NGRID),
                 ntyp(0), rmass(0), rlatt(0), drar(0), drhoar(0),
                 actual(0), actual2(0), rmin(0),
                 petrip(0), rhocon(0), rhomax(0),
                 atpe3b(0), rhotot(0), embf(0), embfp(0),
                 atpe3b1(0), rhotot1(0), embf1(0),
                 atpe3bUMB(0), rhototUMB(0), embfUMB(0),
                 MC_maxd_trial(0), MC_need(0),
                 nbst(0), nbst1(0) {};
    ~EAMFrame() {
#ifdef _USECUDA 
      free_device_ptr();
#endif
    }

    virtual double potential_energyonly_before(int iatom);
    virtual double potential_energyonly_after(int iatom); 
    virtual void potential ();
    virtual void potential_energyonly();
    virtual double potential_energyonly(int iatom);
    virtual double potential_energyonly_change(int iatom); 
    void kraeam();
    void rhoeam();
    void kraMEAM();
    void rhoMEAM();
    int readeam();
    int readMEAM();

    virtual void Alloc();
    virtual void MC_Alloc();
    virtual void MC_Recover();
    virtual void MC_Update();

    virtual void MC_AllocUMB();
    virtual void MC_RecoverUMB();
    virtual void MC_UpdateUMB();

    virtual int MC_FFSprocess();

    virtual void initvars();
    virtual void initparser();
    
    virtual int exec(const char *nam);
#ifdef _USECUDA
    void kraeam_cuda(void);
    void rhoeam_cuda(void);
    void cuda_memory_alloc(void);
    void cuda_memcpy_all(void);
    int test_saxpy(void);
    void free_device_ptr(void);

#endif

#ifdef _PARALLEL
    void Broadcast_EAM_Param();
#endif
};

#endif // _EAM_H

