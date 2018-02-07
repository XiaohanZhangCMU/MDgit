/*

To Do:
 1. Make thread block 3D rectangle and use shared memory
 2. Profiling to tune performance
 3. Multi-phase field implementation
*/

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include<thrust/reduce.h>
#include<thrust/extrema.h>
#include <thrust/inner_product.h>
#include <thrust/transform.h>
#include <thrust/functional.h>
#include <thrust/device_malloc.h>
#include <thrust/device_free.h>

#ifdef _USECUDA

#include "phasefield.h"

#ifdef _USEFLOAT
#define double float
#endif


__global__ void kernel_FreeEnergy_single_step1 (int NX, int NY, int NZ, double gridsize,
                    double U, double EPS, double EPS1, double EPS2, double EPS3, double *R,
                    double *PHI,
                    double *fden, double *dFdPHI, 
                    double *EPS_loc, double *tmp_x, double *tmp_y, double *tmp_z)
/* PHI   : phase field
   fden  : free energy density
   dFdPHI: variational derivative
   dPHIdt: phase field time derivative
*/
{
    int i, j, k, ip, im, jp, jm, kp, km, ind;
    double h2, f, g, phi, dphidx, dphidy, dphidz, d2phi, dphi_SQR;
    double absdphi, EPS_prefactor, triplesum;
    double nx_lab, ny_lab, nz_lab, nx_cryst, ny_cryst, nz_cryst;
    double dnx_lab_dphidx, dny_lab_dphidx, dnz_lab_dphidx; 
    double dnx_lab_dphidy, dny_lab_dphidy, dnz_lab_dphidy; 
    double dnx_lab_dphidz, dny_lab_dphidz, dnz_lab_dphidz; 
    double dEPS_dnx_lab, dEPS_dny_lab, dEPS_dnz_lab;
    double dEPS_dnx_cryst, dEPS_dny_cryst, dEPS_dnz_cryst;
    double dEPS_dphidx, dEPS_dphidy, dEPS_dphidz;

    double my_fden, my_dFdPHI;

    my_fden = 0;
    my_dFdPHI = 0;

    i = blockIdx.x; j = blockIdx.y; k = threadIdx.x;
    
    /* compute spatial gradients */
    h2 = gridsize*gridsize;

    ip = (i+1)%NX; im = (i-1+NX)%NX;
    jp = (j+1)%NY; jm = (j-1+NY)%NY;
    kp = (k+1)%NZ; km = (k-1+NZ)%NZ;

    ind = i*NY*NZ+j*NZ+k;
  
    phi = PHI[ind];

    d2phi = (PHI[ip*NY*NZ+j*NZ+k] + PHI[im*NY*NZ+j*NZ+k]
            +PHI[i*NY*NZ+jp*NZ+k] + PHI[i*NY*NZ+jm*NZ+k]
            +PHI[i*NY*NZ+j*NZ+kp] + PHI[i*NY*NZ+j*NZ+km]
            -PHI[i*NY*NZ+j*NZ+k] * 6.0 ) / h2;

    dphidx = (PHI[ip*NY*NZ+j*NZ+k] - PHI[im*NY*NZ+j*NZ+k]) / (2.0*gridsize);
    dphidy = (PHI[i*NY*NZ+jp*NZ+k] - PHI[i*NY*NZ+jm*NZ+k]) / (2.0*gridsize);
    dphidz = (PHI[i*NY*NZ+j*NZ+kp] - PHI[i*NY*NZ+j*NZ+km]) / (2.0*gridsize);

    f = phi*phi*(1-phi)*(1-phi);
    g = 2.0*phi*(1-phi)*(1-2.0*phi);

    my_fden += f * U; 
    my_dFdPHI = g * U;


    if ((EPS1 == 0)&&(EPS2 == 0)&&(EPS3 == 0))
    { /* isotropic interface energy */
         EPS_loc[ind] = EPS;
         EPS_prefactor = 0.5*SQR(EPS_loc[ind]);

         my_fden += (dphidx*dphidx + dphidy*dphidy + dphidz*dphidz) * EPS_prefactor;
         my_dFdPHI += - 2.0 * d2phi * EPS_prefactor;

    }
    else
    { /* cubic anisotropic interface energy */
         dphi_SQR   = SQR(dphidx) + SQR(dphidy) + SQR(dphidz);
         //absdphi = sqrt(dphi_SQR + 1e-24);
         absdphi = sqrt(dphi_SQR + 1e-24); /* single precision */
         nx_lab = dphidx / absdphi;
         ny_lab = dphidy / absdphi;
         nz_lab = dphidz / absdphi;

         nx_cryst = R[0*3+0]*nx_lab + R[0*3+1]*ny_lab + R[0*3+2]*nz_lab;
         ny_cryst = R[1*3+0]*nx_lab + R[1*3+1]*ny_lab + R[1*3+2]*nz_lab;
         nz_cryst = R[2*3+0]*nx_lab + R[2*3+1]*ny_lab + R[2*3+2]*nz_lab;

         dnx_lab_dphidx = 1./absdphi - dphidx * dphidx / CUBE(absdphi);
         dny_lab_dphidx =            - dphidy * dphidx / CUBE(absdphi);
         dnz_lab_dphidx =            - dphidz * dphidx / CUBE(absdphi);
   
         dnx_lab_dphidy =            - dphidx * dphidy / CUBE(absdphi);
         dny_lab_dphidy = 1./absdphi - dphidy * dphidy / CUBE(absdphi);
         dnz_lab_dphidy =            - dphidz * dphidy / CUBE(absdphi);
   
         dnx_lab_dphidz =            - dphidx * dphidz / CUBE(absdphi);
         dny_lab_dphidz =            - dphidy * dphidz / CUBE(absdphi);
         dnz_lab_dphidz = 1./absdphi - dphidz * dphidz / CUBE(absdphi);
   
         triplesum = SQR(nx_cryst)*SQR(ny_cryst) + SQR(ny_cryst)*SQR(nz_cryst) + SQR(nz_cryst)*SQR(nx_cryst);
         EPS_loc[ind] = EPS + EPS1*triplesum + EPS2*SQR(nx_cryst*ny_cryst*nz_cryst) + EPS3*SQR(triplesum);

         my_fden += 0.5 * SQR(EPS_loc[ind]) * dphi_SQR;

         /* dFdPHI += -d/dx ( SQR(EPS_loc) * dphidx + dphi_SQR * EPS_loc * d(EPS_loc)/d(dphidx) )
                      -d/dy ( SQR(EPS_loc) * dphidy + dphi_SQR * EPS_loc * d(EPS_loc)/d(dphidy) )
                      -d/dz ( SQR(EPS_loc) * dphidz + dphi_SQR * EPS_loc * d(EPS_loc)/d(dphidz) ) ;
         */

         dEPS_dnx_cryst = 2.*EPS1*nx_cryst*(SQR(ny_cryst)+SQR(nz_cryst))
                        + 2.*EPS2*nx_cryst*SQR(ny_cryst)*SQR(nz_cryst)
                        + 4.*EPS3*nx_cryst*(SQR(ny_cryst)+SQR(nz_cryst))*triplesum;
         dEPS_dny_cryst = 2.*EPS1*ny_cryst*(SQR(nz_cryst)+SQR(nx_cryst))
                        + 2.*EPS2*ny_cryst*SQR(nz_cryst)*SQR(nx_cryst)
                        + 4.*EPS3*ny_cryst*(SQR(nz_cryst)+SQR(nx_cryst))*triplesum;
         dEPS_dnz_cryst = 2.*EPS1*nz_cryst*(SQR(nx_cryst)+SQR(ny_cryst))
                        + 2.*EPS2*nz_cryst*SQR(nx_cryst)*SQR(ny_cryst)
                        + 4.*EPS3*nz_cryst*(SQR(nx_cryst)+SQR(ny_cryst))*triplesum;
  
         dEPS_dnx_lab = R[0*3+0]*dEPS_dnx_cryst + R[1*3+0]*dEPS_dny_cryst + R[2*3+0]*dEPS_dnz_cryst;
         dEPS_dny_lab = R[0*3+1]*dEPS_dnx_cryst + R[1*3+1]*dEPS_dny_cryst + R[2*3+1]*dEPS_dnz_cryst;
         dEPS_dnz_lab = R[0*3+2]*dEPS_dnx_cryst + R[1*3+2]*dEPS_dny_cryst + R[2*3+2]*dEPS_dnz_cryst;
 
         dEPS_dphidx = dEPS_dnx_lab*dnx_lab_dphidx + dEPS_dny_lab*dny_lab_dphidx + dEPS_dnz_lab*dnz_lab_dphidx;
         dEPS_dphidy = dEPS_dnx_lab*dnx_lab_dphidy + dEPS_dny_lab*dny_lab_dphidy + dEPS_dnz_lab*dnz_lab_dphidy;
         dEPS_dphidz = dEPS_dnx_lab*dnx_lab_dphidz + dEPS_dny_lab*dny_lab_dphidz + dEPS_dnz_lab*dnz_lab_dphidz;

         //tmp_x[ind] = SQR(EPS_loc)*dphidx + dphi_SQR*EPS_loc*dEPS_dphidx;
         //tmp_y[ind] = SQR(EPS_loc)*dphidy + dphi_SQR*EPS_loc*dEPS_dphidy;
         //tmp_z[ind] = SQR(EPS_loc)*dphidz + dphi_SQR*EPS_loc*dEPS_dphidz;

         /* new scheme for computing variational derivative, 2012/04/27  */
         tmp_x[ind] = dphi_SQR*EPS_loc[ind]*dEPS_dphidx;
         tmp_y[ind] = dphi_SQR*EPS_loc[ind]*dEPS_dphidy;
         tmp_z[ind] = dphi_SQR*EPS_loc[ind]*dEPS_dphidz;

         my_dFdPHI += - 1.0 * SQR(EPS_loc[ind]) * d2phi;
     }

     fden[ind] = my_fden;
     dFdPHI[ind] = my_dFdPHI;
}

__global__ void kernel_FreeEnergy_single_step2 (int NX, int NY, int NZ, double gridsize,
                    double EPS1, double EPS2, double EPS3, double dynamics_type, double K,
                    double *PHI,double *EPS_loc, double *tmp_x, double *tmp_y, double *tmp_z,
                    double *dFdPHI, double *dPHIdt0, double *abs_dPHIdt0) 
                    
{
    int i, j, k, ip, im, jp, jm, kp, km, ind;
    double dphidx, dphidy, dphidz;
    //double d2dFdPHI;

    double my_dFdPHI;

    i = blockIdx.x; j = blockIdx.y; k = threadIdx.x;
    
    /* compute spatial gradients */
    //h2 = gridsize*gridsize;

    ip = (i+1)%NX; im = (i-1+NX)%NX;
    jp = (j+1)%NY; jm = (j-1+NY)%NY;
    kp = (k+1)%NZ; km = (k-1+NZ)%NZ;

    ind = i*NY*NZ+j*NZ+k;
    my_dFdPHI = dFdPHI[ind];
   
    dphidx = (PHI[ip*NY*NZ+j*NZ+k] - PHI[im*NY*NZ+j*NZ+k]) / (2.0*gridsize);
    dphidy = (PHI[i*NY*NZ+jp*NZ+k] - PHI[i*NY*NZ+jm*NZ+k]) / (2.0*gridsize);
    dphidz = (PHI[i*NY*NZ+j*NZ+kp] - PHI[i*NY*NZ+j*NZ+km]) / (2.0*gridsize);

    /* this section must be moved to another kernel call */
    if ((EPS1 != 0)||(EPS2 != 0)||(EPS3 != 0))
    { /* cubic anisotropic interface energy */
           my_dFdPHI +=
                 - (tmp_x[ip*NY*NZ+j*NZ+k]-tmp_x[im*NY*NZ+j*NZ+k]) / (2.0*gridsize)
                 - (tmp_y[i*NY*NZ+jp*NZ+k]-tmp_y[i*NY*NZ+jm*NZ+k]) / (2.0*gridsize)
                 - (tmp_z[i*NY*NZ+j*NZ+kp]-tmp_z[i*NY*NZ+j*NZ+km]) / (2.0*gridsize) ;


           /* new scheme for computing variational derivative, 2012/04/27  */
           my_dFdPHI +=
               - (SQR(EPS_loc[ip*NY*NZ+j*NZ+k])-SQR(EPS_loc[im*NY*NZ+j*NZ+k])) / (2.0*gridsize) * dphidx
               - (SQR(EPS_loc[i*NY*NZ+jp*NZ+k])-SQR(EPS_loc[i*NY*NZ+jm*NZ+k])) / (2.0*gridsize) * dphidy
               - (SQR(EPS_loc[i*NY*NZ+j*NZ+kp])-SQR(EPS_loc[i*NY*NZ+j*NZ+km])) / (2.0*gridsize) * dphidz;
    }

    dFdPHI[ind] = my_dFdPHI;


    if ((dynamics_type == 0) || (dynamics_type == 2)) /* Cahn-Hillard equation */
    {
       dPHIdt0[ind] = -1.0 * K * dFdPHI[i*NY*NZ+j*NZ+k];
       abs_dPHIdt0[ind] = fabs((double) dPHIdt0[ind]);
    }
}

__global__ void kernel_FreeEnergy_single_step3 (int NX, int NY, int NZ, int dynamics_type, double gridsize, double K, double *dFdPHI, double *dPHIdt0, double *dPHIdt, double *abs_dPHIdt0)
{
    int i, j, k, ip, im, jp, jm, kp, km, ind;
    double h2, d2dfdphi;
    i = blockIdx.x; j = blockIdx.y; k = threadIdx.x;
    ind = i*NY*NZ+j*NZ+k;
    h2 = gridsize*gridsize;
    
    
    ip = (i+1)%NX; im = (i-1+NX)%NX;
    jp = (j+1)%NY; jm = (j-1+NY)%NY;
    kp = (k+1)%NZ; km = (k-1+NZ)%NZ; 
     
 if (dynamics_type == 1)
    {      
        /* Cahn-Hillard equation */ 
         d2dfdphi = 
                   (dFdPHI[ip*NY*NZ+j*NZ+k] + dFdPHI[im*NY*NZ+j*NZ+k]
                   +dFdPHI[i*NY*NZ+jp*NZ+k] + dFdPHI[i*NY*NZ+jm*NZ+k]
                   +dFdPHI[i*NY*NZ+j*NZ+kp] + dFdPHI[i*NY*NZ+j*NZ+km]
                   -dFdPHI[ind] * 6.0 ) / h2;

           dPHIdt0[ind] = -1.0 * K * d2dfdphi;
           abs_dPHIdt0[ind] = fabs((double) dPHIdt0[ind]);
    }

     dPHIdt[ind] = dPHIdt0[ind];
}

__global__ void kernel_FreeEnergy_single_step4 (int NX, int NY, int NZ, double avg_dPHIdt0, double *dPHIdt0, double *dPHIdt)
{
     int i, j, k, ind;
     i = blockIdx.x; j = blockIdx.y; k = threadIdx.x;
     ind = i*NY*NZ+j*NZ+k;
     
     dPHIdt[ind] = dPHIdt0[ind] - avg_dPHIdt0;
}


#ifdef _USEFLOAT
float PhaseFieldFrame::cuda_FreeEnergy_single() /* traditional (single) phase field model */
#else
double PhaseFieldFrame::cuda_FreeEnergy_single() /* traditional (single) phase field model */
#endif
{
    double *d_PHI, *d_fden, *d_dFdPHI, *d_dPHIdt, *d_dPHIdt0, *d_EPS_loc, *d_tmp_x, *d_tmp_y, *d_tmp_z;
    double *d_R, *d_abs_dPHIdt0;
    double my_mob;
    
    int size_box = _NX*_NY*_NZ;
    
    if (num_fields != 1)
      FATAL("cuda_FreeEnergy_single: num_fields = "<<num_fields);

    if (gridsize<=0) 
      FATAL("cuda_FreeEnergy_single: gridesize = "<<gridsize);

    cudaMalloc( &d_PHI,    _NX*_NY*_NZ*sizeof(double) );
    cudaMalloc( &d_fden,   _NX*_NY*_NZ*sizeof(double) );
    cudaMalloc( &d_dFdPHI, _NX*_NY*_NZ*sizeof(double) );
    cudaMalloc( &d_dPHIdt, _NX*_NY*_NZ*sizeof(double) );
    cudaMalloc( &d_abs_dPHIdt0, _NX*_NY*_NZ*sizeof(double) );
    cudaMalloc( &d_dPHIdt0,_NX*_NY*_NZ*sizeof(double) );
    cudaMalloc( &d_EPS_loc,_NX*_NY*_NZ*sizeof(double) );
    cudaMalloc( &d_tmp_x,  _NX*_NY*_NZ*sizeof(double) );
    cudaMalloc( &d_tmp_y,  _NX*_NY*_NZ*sizeof(double) );
    cudaMalloc( &d_tmp_z,  _NX*_NY*_NZ*sizeof(double) );
    cudaMalloc( &d_R,      3*3*sizeof(double) );
    
  
    cudaMemcpy( d_PHI,  _PHI[0], _NX*_NY*_NZ*sizeof(double), cudaMemcpyHostToDevice);
#ifndef _USEFLOAT
    cudaMemcpy( d_R,_ROT_MATRIX.element, 3*3*sizeof(double), cudaMemcpyHostToDevice);
#else
    double floatR[3][3];
    floatR[0][0] = _ROT_MATRIX[0][0]; floatR[0][1] = _ROT_MATRIX[0][1]; floatR[0][2] = _ROT_MATRIX[0][2];
    floatR[1][0] = _ROT_MATRIX[1][0]; floatR[1][1] = _ROT_MATRIX[1][1]; floatR[1][2] = _ROT_MATRIX[1][2];
    floatR[2][0] = _ROT_MATRIX[2][0]; floatR[2][1] = _ROT_MATRIX[2][1]; floatR[2][2] = _ROT_MATRIX[2][2];
    cudaMemcpy( d_R,floatR, 3*3*sizeof(double), cudaMemcpyHostToDevice);
#endif

    dim3 blocks(_NX,_NY);
    dim3 threads(_NZ);
    kernel_FreeEnergy_single_step1<<<blocks, threads>>>
                   (_NX, _NY, _NZ, gridsize, _U[0][0], _EPS[0][0], _EPS1[0][0], _EPS2[0][0],
                    _EPS3[0][0], d_R,
                    d_PHI,
                    d_fden, d_dFdPHI, 
                    d_EPS_loc, d_tmp_x, d_tmp_y, d_tmp_z);
    
    _F = 0;
    _Fraw = 0;
    
    /* Use Thrust for free energy reduction */    
    thrust::device_ptr<double> t_fden = thrust::device_pointer_cast(d_fden);
    _F = thrust::reduce(t_fden,t_fden+size_box);
    
    _F *= CUBE(gridsize);
    _Fraw = _F;

    if ((dynamics_type == 0) || (dynamics_type == 2)) /* Cahn-Hillard equation */
    {
          my_mob = Mob_GL;
    } 
    else if (dynamics_type == 1) /* Ginzburg-Landau equation */
    {
          my_mob = Mob_D;
    }
    kernel_FreeEnergy_single_step2<<<blocks, threads>>>
                   (_NX, _NY, _NZ, gridsize,  _EPS1[0][0], _EPS2[0][0], _EPS3[0][0], 
                    dynamics_type, my_mob,
                    d_PHI,d_EPS_loc, d_tmp_x, d_tmp_y, d_tmp_z,
                    d_dFdPHI, d_dPHIdt0, d_abs_dPHIdt0);
    
    if (dynamics_type == 0 || dynamics_type == 1)
    {
    kernel_FreeEnergy_single_step3<<<blocks, threads>>>
                   (_NX, _NY, _NZ, dynamics_type, gridsize, my_mob, d_dFdPHI, d_dPHIdt0, d_dPHIdt, d_abs_dPHIdt0);              
    }         
    else if (dynamics_type == 2) 
    {   

         double avg_dPHIdt0 = 0;
         thrust::device_ptr<double> t_dphidt0 = thrust::device_pointer_cast(d_dPHIdt0);
         avg_dPHIdt0 = thrust::reduce(t_dphidt0,t_dphidt0+size_box);
         avg_dPHIdt0 /= (_NX*_NY*_NZ);
         kernel_FreeEnergy_single_step4<<<blocks, threads>>>
                   (_NX, _NY, _NZ, avg_dPHIdt0, d_dPHIdt0, d_dPHIdt);
     }
     else
     {
        ERROR("unknown dynamics_type = "<<dynamics_type);
     }
     
     thrust::device_ptr<double> t_abs_dPHIdt0 = thrust::device_pointer_cast(d_abs_dPHIdt0);
     thrust::device_ptr<double> max_ptr = thrust::max_element(t_abs_dPHIdt0,t_abs_dPHIdt0+size_box);
     _G = max_ptr[0];
    
     timestep = dtmax/fmax(_G*dtmax/dphimax,1);
    
     cudaMemcpy( _dPHIdt[0], d_dPHIdt, _NX*_NY*_NZ*sizeof(double), cudaMemcpyDeviceToHost);

    cudaFree(d_PHI);
    cudaFree(d_fden);
    cudaFree(d_dFdPHI);
    cudaFree(d_dPHIdt0);
    cudaFree(d_abs_dPHIdt0);
    cudaFree(d_dPHIdt);
    cudaFree(d_EPS_loc);
    cudaFree(d_tmp_x);
    cudaFree(d_tmp_y);
    cudaFree(d_tmp_z);
    cudaFree(d_R);

    cudaError_t error = cudaGetLastError();
    if (error!=cudaSuccess)
    {
       printf("CUDA error: %s\n", cudaGetErrorString(error));
       exit(-1);
    }
    
    return _F;
}




///////////////////////////////////////////////////////////////////////////////////////////////
/* Implementation of the multi-phase field model */



__global__ void kernel_FreeEnergy_multi_interface (int NX, int NY, int NZ, int Ntot, double Mob_D, double *PHI, double *K_S, double *K_L, double *K_LS, double *K_LV, double *K_SV, double *K_other, double *K_D, double *CPHI)
{
int i, j, k, ind;
int ip, im, jp, jm, kp, km;
int ip2, im2, jp2, jm2, kp2, km2;
int ip3, im3, jp3, jm3, kp3, km3;

i = blockIdx.x; j = blockIdx.y; k = threadIdx.x;

ip = (i+1)%NX; im = (i-1+NX)%NX;
jp = (j+1)%NY; jm = (j-1+NY)%NY;
kp = (k+1)%NZ; km = (k-1+NZ)%NZ;

ip2 = (i+2)%NX; im2 = (i-2+NX)%NX;
jp2 = (j+2)%NY; jm2 = (j-2+NY)%NY;
kp2 = (k+2)%NZ; km2 = (k-2+NZ)%NZ;

ip3 = (i+3)%NX; im3 = (i-3+NX)%NX;
jp3 = (j+3)%NY; jm3 = (j-3+NY)%NY;
kp3 = (k+3)%NZ; km3 = (k-3+NZ)%NZ;

ind = i*NY*NZ+j*NZ+k;


if ((((PHI[Ntot + ind]-0.5)*(PHI[Ntot + ip*NY*NZ+j*NZ+k]-0.5)<=0)
|| ((PHI[Ntot + ind]-0.5)*(PHI[Ntot + im*NY*NZ+j*NZ+k]-0.5)<=0)
|| ((PHI[Ntot + ind]-0.5)*(PHI[Ntot + ip2*NY*NZ+j*NZ+k]-0.5)<=0)
|| ((PHI[Ntot + ind]-0.5)*(PHI[Ntot + im2*NY*NZ+j*NZ+k]-0.5)<=0)
|| ((PHI[Ntot + ind]-0.5)*(PHI[Ntot + ip3*NY*NZ+j*NZ+k]-0.5)<=0)
|| ((PHI[Ntot + ind]-0.5)*(PHI[Ntot + im3*NY*NZ+j*NZ+k]-0.5)<=0)
|| ((PHI[Ntot + ind]-0.5)*(PHI[Ntot + i*NY*NZ+jp*NZ+k]-0.5)<=0)
|| ((PHI[Ntot + ind]-0.5)*(PHI[Ntot + i*NY*NZ+jm*NZ+k]-0.5)<=0)
|| ((PHI[Ntot + ind]-0.5)*(PHI[Ntot + i*NY*NZ+jp2*NZ+k]-0.5)<=0)
|| ((PHI[Ntot + ind]-0.5)*(PHI[Ntot + i*NY*NZ+jm2*NZ+k]-0.5)<=0)
|| ((PHI[Ntot + ind]-0.5)*(PHI[Ntot + i*NY*NZ+jp3*NZ+k]-0.5)<=0)
|| ((PHI[Ntot + ind]-0.5)*(PHI[Ntot + i*NY*NZ+jm3*NZ+k]-0.5)<=0)
|| ((PHI[Ntot + ind]-0.5)*(PHI[Ntot + i*NY*NZ+j*NZ+kp]-0.5)<=0)
|| ((PHI[Ntot + ind]-0.5)*(PHI[Ntot + i*NY*NZ+j*NZ+km]-0.5)<=0)
|| ((PHI[Ntot + ind]-0.5)*(PHI[Ntot + i*NY*NZ+j*NZ+kp2]-0.5)<=0)
|| ((PHI[Ntot + ind]-0.5)*(PHI[Ntot + i*NY*NZ+j*NZ+km2]-0.5)<=0)
|| ((PHI[Ntot + ind]-0.5)*(PHI[Ntot + i*NY*NZ+j*NZ+kp3]-0.5)<=0)
|| ((PHI[Ntot + ind]-0.5)*(PHI[Ntot + i*NY*NZ+j*NZ+km3]-0.5)<=0)) )
{
K_S[ind] = 1;
}
else
{
K_S[ind] = 0;
}

if ((((PHI[ind]-0.5)*(PHI[ip*NY*NZ+j*NZ+k]-0.5)<=0)
|| ((PHI[ind]-0.5)*(PHI[im*NY*NZ+j*NZ+k]-0.5)<=0)
|| ((PHI[ind]-0.5)*(PHI[ip2*NY*NZ+j*NZ+k]-0.5)<=0)
|| ((PHI[ind]-0.5)*(PHI[im2*NY*NZ+j*NZ+k]-0.5)<=0)
|| ((PHI[ind]-0.5)*(PHI[ip3*NY*NZ+j*NZ+k]-0.5)<=0)
|| ((PHI[ind]-0.5)*(PHI[im3*NY*NZ+j*NZ+k]-0.5)<=0)
|| ((PHI[ind]-0.5)*(PHI[i*NY*NZ+jp*NZ+k]-0.5)<=0)
|| ((PHI[ind]-0.5)*(PHI[i*NY*NZ+jm*NZ+k]-0.5)<=0)
|| ((PHI[ind]-0.5)*(PHI[i*NY*NZ+jp2*NZ+k]-0.5)<=0)
|| ((PHI[ind]-0.5)*(PHI[i*NY*NZ+jm2*NZ+k]-0.5)<=0)
|| ((PHI[ind]-0.5)*(PHI[i*NY*NZ+jp3*NZ+k]-0.5)<=0)
|| ((PHI[ind]-0.5)*(PHI[i*NY*NZ+jm3*NZ+k]-0.5)<=0)
|| ((PHI[ind]-0.5)*(PHI[i*NY*NZ+j*NZ+kp]-0.5)<=0)
|| ((PHI[ind]-0.5)*(PHI[i*NY*NZ+j*NZ+km]-0.5)<=0)
|| ((PHI[ind]-0.5)*(PHI[i*NY*NZ+j*NZ+kp2]-0.5)<=0)
|| ((PHI[ind]-0.5)*(PHI[i*NY*NZ+j*NZ+km2]-0.5)<=0)
|| ((PHI[ind]-0.5)*(PHI[i*NY*NZ+j*NZ+kp3]-0.5)<=0)
|| ((PHI[ind]-0.5)*(PHI[i*NY*NZ+j*NZ+km3]-0.5)<=0)) )
{
K_L[ind] = 1;
}
else
{
K_L[ind] = 0;
}

if (( K_L[ind]==1 ) && ( K_S[ind]==1 ))
{
K_LS[ind] = 1;
}
else
{
K_LS[ind] = 0;
}

if (( K_L[ind]==0 ) && ( K_S[ind]==1 ))
{
K_SV[ind] = 1;
}
else
{
K_SV[ind] = 0;
}

if (( K_L[ind]==1 ) && ( K_S[ind]==0 ))
{
K_LV[ind] = 1;
}
else
{
K_LV[ind] = 0;
}

K_other[ind] = 1 - K_LV[ind] - K_SV[ind] - K_LS[ind];
    
K_D[ind] = Mob_D;
    
CPHI[ind] = tanh((PHI[ind]-0.5)/0.10)/2 + 0.5;

}


__global__ void kernel_FreeEnergy_multi_freeenergy (int NX, int NY, int NZ, int Ntot, double gridsize, int num_fields, double Fext_x, double Fext_y, double Penalty_3phase, double Penalty_NW, double *PHI, double *fden, double *fden_raw, double *dFdPHIi, double *d2PHI, double *EPS_loc, double *dPHIdx, double *dPHIdy, double *dPHIdz, double *tmp_x, double *tmp_y, double *tmp_z, double *matrix_x, double *matrix_y, double *U, double *EPS, double *EPS1, double *EPS2, double *EPS3, double *MU, double *R, double *NW_orig)
{
    int n, i, j, k, ind;
    int ip, im, jp, jm, kp, km;
    double h2, my_fden, my_fden_raw, my_dfdphii, fp;
    double d2phi, dphidx, dphidy, dphidz, grad_term, absdphi;
    double nx_lab, ny_lab, nz_lab, nx_cryst, ny_cryst, nz_cryst;
    double eps1, eps2, eps3, triplesum, eps_loc;
    double EPS_prefactor;
    double dnx_lab_dphidx, dny_lab_dphidx, dnz_lab_dphidx;
    double dnx_lab_dphidy, dny_lab_dphidy, dnz_lab_dphidy;
    double dnx_lab_dphidz, dny_lab_dphidz, dnz_lab_dphidz;
    double dEPS_dnx_lab, dEPS_dny_lab, dEPS_dnz_lab;
    double dEPS_dnx_cryst, dEPS_dny_cryst, dEPS_dnz_cryst;
    double dEPS_dphidx, dEPS_dphidy, dEPS_dphidz;
    double dEPS_dphidx_loc, dEPS_dphidy_loc, dEPS_dphidz_loc, grad_loc_SQR;

    i = blockIdx.x; j = blockIdx.y; k = threadIdx.x;

    ip = (i+1)%NX; im = (i-1+NX)%NX;
    jp = (j+1)%NY; jm = (j-1+NY)%NY;
    kp = (k+1)%NZ; km = (k-1+NZ)%NZ;

    h2 = gridsize*gridsize;
    ind = i*NY*NZ+j*NZ+k;
    my_fden = 0; my_fden_raw = 0;

    matrix_x[ind] = (i - NX/2)*gridsize;
    matrix_y[ind] = (j - NY/2)*gridsize;


    for(n=0;n<num_fields;n++)
    {
        my_dfdphii = 0;
        
        d2phi = (PHI[n*Ntot+ip*NY*NZ+j*NZ+k] + PHI[n*Ntot+im*NY*NZ+j*NZ+k]
        +PHI[n*Ntot+i*NY*NZ+jp*NZ+k] + PHI[n*Ntot+i*NY*NZ+jm*NZ+k]
        +PHI[n*Ntot+i*NY*NZ+j*NZ+kp] + PHI[n*Ntot+i*NY*NZ+j*NZ+km]
        -PHI[n*Ntot+ind] * 6.0 ) / h2;

        d2PHI[n*Ntot+ind] = d2phi;

        dphidx = (PHI[n*Ntot+ip*NY*NZ+j*NZ+k] - PHI[n*Ntot+im*NY*NZ+j*NZ+k]) / (2.0*gridsize);
        dphidy = (PHI[n*Ntot+i*NY*NZ+jp*NZ+k] - PHI[n*Ntot+i*NY*NZ+jm*NZ+k]) / (2.0*gridsize);
        dphidz = (PHI[n*Ntot+i*NY*NZ+j*NZ+kp] - PHI[n*Ntot+i*NY*NZ+j*NZ+km]) / (2.0*gridsize);

        dPHIdx[n*Ntot+ind] = dphidx;
        dPHIdy[n*Ntot+ind] = dphidy;
        dPHIdz[n*Ntot+ind] = dphidz;

        grad_term = SQR(dphidx) + SQR(dphidy) + SQR(dphidz) ;

        fp = tanh((PHI[n*Ntot+ind]-0.5)/0.10)/2 + 0.5;
        my_fden += SQR(PHI[n*Ntot+ind]) * SQR(1-PHI[n*Ntot+ind])*U[n*num_fields+n] + fp * MU[n];
        my_fden_raw += SQR(PHI[n*Ntot+ind]) * SQR(1-PHI[n*Ntot+ind])*U[n*num_fields+n];
        
        if ( n == 0 )
        {
            my_fden += Fext_x * fp * matrix_x[ind] + Fext_y * fp * matrix_y[ind];
        }

        // Penalty term for three phase co-existing
        my_fden += Penalty_3phase/num_fields*SQR(PHI[ind])*SQR(PHI[Ntot+ind])*SQR(PHI[2*Ntot+ind]);
        
        my_fden += Penalty_NW/num_fields*SQR(PHI[Ntot+ind] - NW_orig[ind]);

        if ( (EPS1[n*num_fields+n] == 0)&&(EPS2[n*num_fields+n] == 0)&&(EPS3[n*num_fields+n] == 0) )
        {
            /* isotropic interface energy */
            eps_loc = 1.0;
            EPS_prefactor = SQR(eps_loc);
            my_fden += EPS[n*num_fields+n] * EPS_prefactor * grad_term;
            my_fden_raw += EPS[n*num_fields+n] * EPS_prefactor * grad_term;

            /* the equation of motion should follow the Steinbach 1999 formulation */
            my_dfdphii += U[n*num_fields+n]*(4.0*PHI[n*Ntot+ind]*PHI[n*Ntot+ind]*PHI[n*Ntot+ind]
                       + 2.0*PHI[n*Ntot+ind] - 6.0*PHI[n*Ntot+ind]*PHI[n*Ntot+ind])
                       + MU[n]*(5.0-5.0*SQR(tanh(10*PHI[n*Ntot+ind]-5.0)))
                       - 2.0*EPS[n*num_fields+n]*EPS_prefactor*d2PHI[n*Ntot+ind];

            // with the external force
            if ( n==0 )
            {
                my_dfdphii += Fext_x * (5.0-5.0*SQR(tanh(10*PHI[n*Ntot+ind]-5.0)))*matrix_x[ind]
                            + Fext_y * (5.0-5.0*SQR(tanh(10*PHI[n*Ntot+ind]-5.0)))*matrix_y[ind];
            }

        }
        else
        {
            /* anisotropic interface energy */

            /* surface normal defined by solid phase gradient */

            absdphi = sqrt(grad_term + 1e-24);

            nx_lab = dphidx / absdphi;
            ny_lab = dphidy / absdphi;
            nz_lab = dphidz / absdphi;
                
            nx_cryst = R[0*3+0]*nx_lab + R[0*3+1]*ny_lab + R[0*3+2]*nz_lab;
            ny_cryst = R[1*3+0]*nx_lab + R[1*3+1]*ny_lab + R[1*3+2]*nz_lab;
            nz_cryst = R[2*3+0]*nx_lab + R[2*3+1]*ny_lab + R[2*3+2]*nz_lab;

            eps1 = EPS1[n*num_fields+n];
            eps2 = EPS2[n*num_fields+n];
            eps3 = EPS3[n*num_fields+n];

            triplesum = SQR(nx_cryst)*SQR(ny_cryst) + SQR(ny_cryst)*SQR(nz_cryst) + SQR(nz_cryst)*SQR(nx_cryst);

            /* Notice the difference from the Matlab implementation */
            eps_loc   = 1.0 + eps1*triplesum + eps2*SQR(nx_cryst*ny_cryst*nz_cryst)
                      + eps3*SQR(triplesum);
            
            EPS_loc[n*Ntot+ind] = eps_loc;

            my_fden += EPS[n*num_fields+n] * SQR(eps_loc) * grad_term;
            my_fden_raw += EPS[n*num_fields+n] * SQR(eps_loc) * grad_term;

            dnx_lab_dphidx = 1./absdphi - dphidx * dphidx / CUBE(absdphi);
            dny_lab_dphidx =            - dphidy * dphidx / CUBE(absdphi);
            dnz_lab_dphidx =            - dphidz * dphidx / CUBE(absdphi);

            dnx_lab_dphidy =            - dphidx * dphidy / CUBE(absdphi);
            dny_lab_dphidy = 1./absdphi - dphidy * dphidy / CUBE(absdphi);
            dnz_lab_dphidy =            - dphidz * dphidy / CUBE(absdphi);

            dnx_lab_dphidz =            - dphidx * dphidz / CUBE(absdphi);
            dny_lab_dphidz =            - dphidy * dphidz / CUBE(absdphi);
            dnz_lab_dphidz = 1./absdphi - dphidz * dphidz / CUBE(absdphi);

            dEPS_dnx_cryst = 2.*eps1*nx_cryst*(SQR(ny_cryst)+SQR(nz_cryst))
            + 2.*eps2*nx_cryst*SQR(ny_cryst)*SQR(nz_cryst)
            + 4.*eps3*nx_cryst*(SQR(ny_cryst)+SQR(nz_cryst))*triplesum;
            dEPS_dny_cryst = 2.*eps1*ny_cryst*(SQR(nz_cryst)+SQR(nx_cryst))
            + 2.*eps2*ny_cryst*SQR(nz_cryst)*SQR(nx_cryst)
            + 4.*eps3*ny_cryst*(SQR(nz_cryst)+SQR(nx_cryst))*triplesum;
            dEPS_dnz_cryst = 2.*eps1*nz_cryst*(SQR(nx_cryst)+SQR(ny_cryst))
            + 2.*eps2*nz_cryst*SQR(nx_cryst)*SQR(ny_cryst)
            + 4.*eps3*nz_cryst*(SQR(nx_cryst)+SQR(ny_cryst))*triplesum;

            dEPS_dnx_lab = R[0*3+0]*dEPS_dnx_cryst + R[1*3+0]*dEPS_dny_cryst + R[2*3+0]*dEPS_dnz_cryst;
            dEPS_dny_lab = R[0*3+1]*dEPS_dnx_cryst + R[1*3+1]*dEPS_dny_cryst + R[2*3+1]*dEPS_dnz_cryst;
            dEPS_dnz_lab = R[0*3+2]*dEPS_dnx_cryst + R[1*3+2]*dEPS_dny_cryst + R[2*3+2]*dEPS_dnz_cryst;

            dEPS_dphidx = dEPS_dnx_lab*dnx_lab_dphidx + dEPS_dny_lab*dny_lab_dphidx + dEPS_dnz_lab*dnz_lab_dphidx;
            dEPS_dphidy = dEPS_dnx_lab*dnx_lab_dphidy + dEPS_dny_lab*dny_lab_dphidy + dEPS_dnz_lab*dnz_lab_dphidy;
            dEPS_dphidz = dEPS_dnx_lab*dnx_lab_dphidz + dEPS_dny_lab*dny_lab_dphidz + dEPS_dnz_lab*dnz_lab_dphidz;

            dEPS_dphidx_loc = dEPS_dphidx;
            dEPS_dphidy_loc = dEPS_dphidy;
            dEPS_dphidz_loc = dEPS_dphidz;
     
            grad_loc_SQR = SQR(dphidx) + SQR(dphidy) + SQR(dphidz);

            /* new scheme for computing variational derivative, 2012/04/27  */
            tmp_x[n*Ntot+ind] = grad_loc_SQR*eps_loc*dEPS_dphidx_loc;
            tmp_y[n*Ntot+ind] = grad_loc_SQR*eps_loc*dEPS_dphidy_loc;
            tmp_z[n*Ntot+ind] = grad_loc_SQR*eps_loc*dEPS_dphidz_loc;

            my_dfdphii += U[n*num_fields+n]*(4*PHI[n*Ntot+ind]*PHI[n*Ntot+ind]*PHI[n*Ntot+ind]
                        + 2*PHI[n*Ntot+ind] - 6*PHI[n*Ntot+ind]*PHI[n*Ntot+ind])
                        + MU[n]*(5.0-5.0*SQR(tanh(10*PHI[n*Ntot+ind]-5.0)));

            my_dfdphii += (-2.0) * EPS[n*num_fields+n]*SQR(eps_loc) * d2phi;

            // with the external force term
            if ( n==0 )
            {
                my_dfdphii += Fext_x * (5.0-5.0*SQR(tanh(10*PHI[n*Ntot+ind]-5.0)))*matrix_x[ind]
                            + Fext_y * (5.0-5.0*SQR(tanh(10*PHI[n*Ntot+ind]-5.0)))*matrix_y[ind];
            }
        }
        dFdPHIi[n*Ntot+ind] = my_dfdphii;
    }
    
    fden[ind] = my_fden;
    fden_raw[ind] = my_fden_raw;
}

__global__ void kernel_FreeEnergy_multi_freeenergy2 (int NX, int NY, int NZ, int Ntot, int num_fields, double gridsize, double Penalty_3phase, double Penalty_NW, double *PHI, double *dFdPHIi, double *tmp_x, double *tmp_y, double *tmp_z, double *EPS_loc, double *dPHIdx, double *dPHIdy, double *dPHIdz, double *EPS, double *EPS1, double *EPS2, double *EPS3, double *NW_orig)
{
    int n, i, j, k, ind;
    int ip, im, jp, jm, kp, km;
    double grad_loc_x, grad_loc_y, grad_loc_z;
    double my_dfdphii;
    double my_dfdphi1, my_dfdphi2, my_dfdphi3;
    
    i = blockIdx.x; j = blockIdx.y; k = threadIdx.x;
    
    ip = (i+1)%NX; im = (i-1+NX)%NX;
    jp = (j+1)%NY; jm = (j-1+NY)%NY;
    kp = (k+1)%NZ; km = (k-1+NZ)%NZ;
    
    ind = i*NY*NZ+j*NZ+k;
    
    for(n=0;n<num_fields;n++)
    {
        my_dfdphii = dFdPHIi[n*Ntot+ind];
        
        if ( (EPS1[n*num_fields+n] != 0)||(EPS2[n*num_fields+n] != 0)||(EPS3[n*num_fields+n] != 0) )
        {
            my_dfdphii +=
            + EPS[n*num_fields+n]*(-2.0)*(tmp_x[n*Ntot+ip*NY*NZ+j*NZ+k]-tmp_x[n*Ntot+im*NY*NZ+j*NZ+k]) / (2.0*gridsize)
            + EPS[n*num_fields+n]*(-2.0)*(tmp_y[n*Ntot+i*NY*NZ+jp*NZ+k]-tmp_y[n*Ntot+i*NY*NZ+jm*NZ+k]) / (2.0*gridsize)
            + EPS[n*num_fields+n]*(-2.0)*(tmp_z[n*Ntot+i*NY*NZ+j*NZ+kp]-tmp_z[n*Ntot+i*NY*NZ+j*NZ+km]) / (2.0*gridsize);
            
            /* new scheme for computing variational derivative, 2012/04/27  */
            grad_loc_x = dPHIdx[n*Ntot+ind];
            grad_loc_y = dPHIdy[n*Ntot+ind];
            grad_loc_z = dPHIdz[n*Ntot+ind];
            
            my_dfdphii +=
            +EPS[n*num_fields+n]*(-4.0)*EPS_loc[n*Ntot+ind]*(EPS_loc[n*Ntot+ip*NY*NZ+j*NZ+k]-EPS_loc[n*Ntot+im*NY*NZ+j*NZ+k])/ (2.0*gridsize) * grad_loc_x
            +EPS[n*num_fields+n]*(-4.0)*EPS_loc[n*Ntot+ind]*(EPS_loc[n*Ntot+i*NY*NZ+jp*NZ+k]-EPS_loc[n*Ntot+i*NY*NZ+jm*NZ+k])/ (2.0*gridsize) * grad_loc_y
            +EPS[n*num_fields+n]*(-4.0)*EPS_loc[n*Ntot+ind]*(EPS_loc[n*Ntot+i*NY*NZ+j*NZ+kp]-EPS_loc[n*Ntot+i*NY*NZ+j*NZ+km])/ (2.0*gridsize) * grad_loc_z;
        }
        
        dFdPHIi[n*Ntot+ind] = my_dfdphii;
    }
    
    my_dfdphi1 = dFdPHIi[ind];
    my_dfdphi2 = dFdPHIi[Ntot+ind];
    my_dfdphi3 = dFdPHIi[2*Ntot+ind];
    
    my_dfdphi1 += 2*Penalty_3phase*PHI[ind]*SQR(PHI[Ntot+ind])*SQR(PHI[2*Ntot+ind]);
    my_dfdphi2 += 2*Penalty_3phase*PHI[Ntot+ind]*SQR(PHI[ind])*SQR(PHI[2*Ntot+ind]);
    my_dfdphi3 += 2*Penalty_3phase*PHI[2*Ntot+ind]*SQR(PHI[ind])*SQR(PHI[Ntot+ind]);
    
    my_dfdphi2 += 2*Penalty_NW*(PHI[Ntot+ind] - NW_orig[ind]);
    
    dFdPHIi[ind] = my_dfdphi1;
    dFdPHIi[Ntot+ind] = my_dfdphi2;
    dFdPHIi[2*Ntot+ind] = my_dfdphi3;
}

__global__ void kernel_FreeEnergy_multi_EOM (int NX, int NY, int NZ, int Ntot, int dynamics_type, double Mob_GL, double Mob_LV, double Mob_LS, double Mob_SV, double *dPHIdt0, double *dPHIdtCH, double *dFdPHIi, double *K_LS, double *K_LV, double *K_SV, double *K_other, double *K_D, double *abs_dPHIdt)
{
    int i, j, k, ind;
    double my_dfdphi1, my_dfdphi2, my_dfdphi3;

    
    i = blockIdx.x; j = blockIdx.y; k = threadIdx.x;
    
    ind = i*NY*NZ+j*NZ+k;
    
    my_dfdphi1 = dFdPHIi[ind];
    my_dfdphi2 = dFdPHIi[Ntot+ind];
    my_dfdphi3 = dFdPHIi[2*Ntot+ind];
    
    /* Steinback 1999 formulation */
    
    dPHIdt0[ind] = - K_SV[ind]*(Mob_GL*(my_dfdphi1 - my_dfdphi2)
                                    + Mob_GL*(my_dfdphi1 - my_dfdphi3))
    - K_LV[ind]*(Mob_GL*((my_dfdphi1+my_dfdphi3)/2 - my_dfdphi2)
                  + Mob_LV*(my_dfdphi1 - my_dfdphi3))
    - K_LS[ind]*(Mob_GL*((my_dfdphi1+my_dfdphi2)/2 - my_dfdphi3)
                  + Mob_LS*(my_dfdphi1 - my_dfdphi2))
    - K_other[ind]*(Mob_GL*(my_dfdphi1 - my_dfdphi2)
                     + Mob_GL*(my_dfdphi1 - my_dfdphi3));
    
    dPHIdt0[Ntot+ind] = - K_LV[ind]*(Mob_GL*(my_dfdphi2 - my_dfdphi3)
                                     + Mob_GL*(my_dfdphi2 - my_dfdphi1))
    - K_SV[ind]*(Mob_GL*((my_dfdphi2+my_dfdphi3)/2 - my_dfdphi1)
                  + Mob_SV*(my_dfdphi2 - my_dfdphi3))
    - K_LS[ind]*(Mob_GL*((my_dfdphi1+my_dfdphi2)/2 - my_dfdphi3)
                  + Mob_LS*(my_dfdphi2 - my_dfdphi1))
    - K_other[ind]*(Mob_GL*(my_dfdphi2 - my_dfdphi1)
                     + Mob_GL*(my_dfdphi2 - my_dfdphi3));
    
    dPHIdt0[2*Ntot+ind] = - K_LS[ind]*(Mob_GL*(my_dfdphi3 - my_dfdphi2)
                                     + Mob_GL*(my_dfdphi3 - my_dfdphi1))
    - K_LV[ind]*(Mob_GL*((my_dfdphi1+my_dfdphi3)/2 - my_dfdphi2)
                  + Mob_LV*(my_dfdphi3 - my_dfdphi1))
    - K_SV[ind]*(Mob_GL*((my_dfdphi2+my_dfdphi3)/2 - my_dfdphi1)
                  + Mob_SV*(my_dfdphi3 - my_dfdphi2))
    - K_other[ind]*(Mob_GL*(my_dfdphi3 - my_dfdphi1)
                     + Mob_GL*(my_dfdphi3 - my_dfdphi2));
    
    /* C-H equation of motion */
    dPHIdtCH[ind] = - K_D[ind]*(my_dfdphi1 - my_dfdphi2)
    - K_D[ind]*(my_dfdphi1 - my_dfdphi3);
    
    dPHIdtCH[Ntot+ind] = - K_D[ind]*(my_dfdphi2 - my_dfdphi1)
    - K_D[ind]*(my_dfdphi2 - my_dfdphi3);
    
    dPHIdtCH[2*Ntot+ind] = - K_D[ind]*(my_dfdphi3 - my_dfdphi1)
    - K_D[ind]*(my_dfdphi3 - my_dfdphi2);

    
    if (dynamics_type != 1)
    {
        abs_dPHIdt[ind] = fabs((double) dPHIdt0[ind]);
        abs_dPHIdt[Ntot+ind] = fabs((double) dPHIdt0[Ntot+ind]);
        abs_dPHIdt[2*Ntot+ind] = fabs((double) dPHIdt0[2*Ntot+ind]);
    }
    else
    {
        abs_dPHIdt[ind] = fabs((double) dPHIdtCH[ind]);
        abs_dPHIdt[Ntot+ind] = fabs((double) dPHIdtCH[Ntot+ind]);
        abs_dPHIdt[2*Ntot+ind] = fabs((double) dPHIdtCH[2*Ntot+ind]);
    }
}

__global__ void kernel_FreeEnergy_multi_dFdmu (int NX, int NY, int NZ, int Ntot, double Mob_GL, double Mob_LS, double Mob_SV, double Mob_LV, double *PHI, double *Cprime0, double *Cprime2, double *K_LS, double *K_LV, double *K_SV, double *K_other, double *dGLdmuL, double *dGLdmuV, double *dGVdmuL, double *dGVdmuV, double *dGSdmuL, double *dGSdmuV)
{
    int i, j, k, ind;
    double cprime0, cprime2;
    i = blockIdx.x; j = blockIdx.y; k = threadIdx.x;
    ind = i*NY*NZ+j*NZ+k;
    
    cprime0 = (5.0-5.0*SQR(tanh(10*PHI[ind]-5.0)));
    cprime2 = (5.0-5.0*SQR(tanh(10*PHI[2*Ntot+ind]-5.0)));
    Cprime0[ind] = cprime0;
    Cprime2[ind] = cprime2;
    
    dGLdmuL[ind] = - K_LS[ind]*(Mob_GL/2+Mob_LS)*cprime0
                 - K_SV[ind]*(Mob_GL+Mob_GL)*cprime0
                 - K_LV[ind]*(Mob_GL/2+Mob_LV)*cprime0
                 - K_other[ind]*(Mob_GL+Mob_GL)*cprime0;
    
    dGLdmuV[ind] = - K_LS[ind]*(-Mob_GL)*cprime2
                 - K_SV[ind]*(-Mob_GL)*cprime2
                 - K_LV[ind]*(Mob_GL/2-Mob_LV)*cprime2
                 - K_other[ind]*(-Mob_GL)*cprime2;
    
    dGVdmuL[ind] = - K_LS[ind]*(-Mob_GL)*cprime0
                 - K_SV[ind]*(-Mob_GL)*cprime0
                 - K_LV[ind]*(Mob_GL/2-Mob_LV)*cprime0
                 - K_other[ind]*(-Mob_GL)*cprime0;
    
    dGVdmuV[ind] = - K_LS[ind]*(Mob_GL+Mob_GL)*cprime2
                 - K_SV[ind]*(Mob_GL/2+Mob_SV)*cprime2
                 - K_LV[ind]*(Mob_GL/2+Mob_LV)*cprime2
                 - K_other[ind]*(Mob_GL+Mob_GL)*cprime2;
    
    dGSdmuL[ind] = - K_LS[ind]*(Mob_GL/2-Mob_LS)*cprime0
                 - K_SV[ind]*(-Mob_GL)*cprime0
                 - K_LV[ind]*(-Mob_GL)*cprime0
                 - K_other[ind]*(-Mob_GL)*cprime0;
    
    dGSdmuV[ind] = - K_LS[ind]*(-Mob_GL)*cprime2
                 - K_SV[ind]*(Mob_GL/2-Mob_SV)*cprime2
                 - K_LV[ind]*(-Mob_GL)*cprime2
                 - K_other[ind]*(-Mob_GL)*cprime2;
}

__global__ void kernel_FreeEnergy_multi_dynamics0 (int NX, int NY, int NZ, int Ntot, int num_fields, double *dPHIdt0, double *dPHIdt)
{
    int n, i, j, k, ind;
    double dphidt0;
    i = blockIdx.x; j = blockIdx.y; k = threadIdx.x;
    ind = i*NY*NZ+j*NZ+k;
    for(n=0;n<num_fields;n++)
    {
        dphidt0 = dPHIdt0[n*Ntot+ind];
        dPHIdt[n*Ntot+ind] = dphidt0;
    }
    
}

__global__ void kernel_FreeEnergy_multi_dynamics1 (int NX, int NY, int NZ, int Ntot, int num_fields, double gridsize, double *dPHIdtCH, double *dPHIdt)
{
    int n, i, j, k, ind;
    int ip, im, jp, jm, kp, km;
    double h2, d2dphidt0;
    
    i = blockIdx.x; j = blockIdx.y; k = threadIdx.x;
    
    ip = (i+1)%NX; im = (i-1+NX)%NX;
    jp = (j+1)%NY; jm = (j-1+NY)%NY;
    kp = (k+1)%NZ; km = (k-1+NZ)%NZ;
    
    ind = i*NY*NZ+j*NZ+k;
    h2 = gridsize*gridsize;
    
    for(n=0;n<num_fields;n++)
    {
        d2dphidt0 = 0;
        
        d2dphidt0 = (dPHIdtCH[n*Ntot+ip*NY*NZ+j*NZ+k] + dPHIdtCH[n*Ntot+im*NY*NZ+j*NZ+k]
                  + dPHIdtCH[n*Ntot+i*NY*NZ+jp*NZ+k] + dPHIdtCH[n*Ntot+i*NY*NZ+jm*NZ+k]
                  + dPHIdtCH[n*Ntot+i*NY*NZ+j*NZ+kp] + dPHIdtCH[n*Ntot+i*NY*NZ+j*NZ+km]
                  - dPHIdtCH[n*Ntot+ind] * 6.0 ) / h2;
        
        dPHIdt[n*Ntot+ind] = - d2dphidt0;
    }
}

__global__ void kernel_FreeEnergy_multi_mprime(int NX, int NY, int NZ, int Ntot, double Mob_GL, double Mob_LS, double Mob_SV, double Mob_LV, double *PHI, double *mprime, double *K_LS, double *K_SV, double *K_LV, double *K_other, double *Cprime0)
{
    int i, j, k, ind;
    double cprime0;

    i = blockIdx.x; j = blockIdx.y; k = threadIdx.x;
    ind = i*NY*NZ+j*NZ+k;

    cprime0 = (5.0-5.0*SQR(tanh(10*PHI[ind]-5.0)));

    Cprime0[ind] = cprime0;


    mprime[ind] = - K_LS[ind]*(Mob_GL/2+Mob_LS)*(5.0-5.0*SQR(tanh(10*PHI[ind]-5.0)))
                  - K_SV[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*PHI[ind]-5.0)))
                  - K_LV[ind]*(Mob_GL/2+Mob_LV)*(5.0-5.0*SQR(tanh(10*PHI[ind]-5.0)))
                  - K_other[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*PHI[ind]-5.0)));
}

__global__ void kernel_FreeEnergy_multi_dynamics2(int NX, int NY, int NZ, int Ntot, double Mob_GL, double Mob_LS, double Mob_SV, double Mob_LV, double dmu1, double *PHI, double *dPHIdt0, double *dPHIdt, double *mprime, double *K_LS, double *K_SV, double *K_LV, double *K_other)
{
    int i, j, k, ind;
    i = blockIdx.x; j = blockIdx.y; k = threadIdx.x;
    ind = i*NY*NZ+j*NZ+k;

    dPHIdt[ind] = dPHIdt0[ind] + dmu1 * mprime[ind];

    dPHIdt[Ntot+ind] = dPHIdt0[Ntot+ind] + dmu1 * (K_LS[ind]*(-Mob_GL/2+Mob_LS)+K_LV[ind]*Mob_GL+K_SV[ind]*Mob_GL+K_other[ind]*Mob_GL)*(5.0-5.0*SQR(tanh(10.0*PHI[ind]-5.0)));

    dPHIdt[2*Ntot+ind] = dPHIdt0[2*Ntot+ind] + dmu1 * (K_LS[ind]*Mob_GL+K_LV[ind]*(-Mob_GL/2+Mob_LV)+K_SV[ind]*Mob_GL+K_other[ind]*Mob_GL)*(5.0-5.0*SQR(tanh(10.0*PHI[ind]-5.0)));
}

__global__ void kernel_FreeEnergy_multi_dynamics3(int NX, int NY, int NZ, int Ntot, double gridsize, double Mob_GL, double Mob_LS, double Mob_SV, double Mob_LV, double dmu1, double *PHI, double *dPHIdt0, double *dPHIdt, double *dPHIdtCH, double *mprime, double *K_LS, double *K_SV, double *K_LV, double *K_other)
{
    int i, j, k, ind;
    int ip, im, jp, jm, kp, km;
    double h2, d2dphidt0_1, d2dphidt0_2, d2dphidt0_3;

    i = blockIdx.x; j = blockIdx.y; k = threadIdx.x;
    
    ip = (i+1)%NX; im = (i-1+NX)%NX;
    jp = (j+1)%NY; jm = (j-1+NY)%NY;
    kp = (k+1)%NZ; km = (k-1+NZ)%NZ;
    
    ind = i*NY*NZ+j*NZ+k;
    h2 = gridsize*gridsize;
    
    d2dphidt0_1 = 0; d2dphidt0_2 = 0; d2dphidt0_3 = 0;
    
        
    d2dphidt0_1 = (dPHIdtCH[ip*NY*NZ+j*NZ+k] + dPHIdtCH[im*NY*NZ+j*NZ+k]
                + dPHIdtCH[i*NY*NZ+jp*NZ+k] + dPHIdtCH[i*NY*NZ+jm*NZ+k]
                + dPHIdtCH[i*NY*NZ+j*NZ+kp] + dPHIdtCH[i*NY*NZ+j*NZ+km]
                - dPHIdtCH[ind] * 6.0 ) / h2;
        
    d2dphidt0_2 = (dPHIdtCH[Ntot+ip*NY*NZ+j*NZ+k] + dPHIdtCH[Ntot+im*NY*NZ+j*NZ+k]
                + dPHIdtCH[Ntot+i*NY*NZ+jp*NZ+k] + dPHIdtCH[Ntot+i*NY*NZ+jm*NZ+k]
                + dPHIdtCH[Ntot+i*NY*NZ+j*NZ+kp] + dPHIdtCH[Ntot+i*NY*NZ+j*NZ+km]
                - dPHIdtCH[Ntot+ind] * 6.0 ) / h2;
    
    d2dphidt0_3 = (dPHIdtCH[2*Ntot+ip*NY*NZ+j*NZ+k] + dPHIdtCH[2*Ntot+im*NY*NZ+j*NZ+k]
                + dPHIdtCH[2*Ntot+i*NY*NZ+jp*NZ+k] + dPHIdtCH[2*Ntot+i*NY*NZ+jm*NZ+k]
                + dPHIdtCH[2*Ntot+i*NY*NZ+j*NZ+kp] + dPHIdtCH[2*Ntot+i*NY*NZ+j*NZ+km]
                - dPHIdtCH[2*Ntot+ind] * 6.0 ) / h2;

    
    dPHIdt[ind] = dPHIdt0[ind] + dmu1* mprime[ind] - d2dphidt0_1;
    
    dPHIdt[Ntot+ind] = dPHIdt0[Ntot+ind] + dmu1 * (K_LS[ind]*(-Mob_GL/2+Mob_LS)+K_LV[ind]*Mob_GL+K_SV[ind]*Mob_GL+K_other[ind]*Mob_GL)*(5.0-5.0*SQR(tanh(10.0*PHI[ind]-5.0))) - d2dphidt0_2;
    
    dPHIdt[2*Ntot+ind] = dPHIdt0[2*Ntot+ind] + dmu1 * (K_LS[ind]*Mob_GL+K_LV[ind]*(-Mob_GL/2+Mob_LV)+K_SV[ind]*Mob_GL+K_other[ind]*Mob_GL)*(5.0-5.0*SQR(tanh(10.0*PHI[ind]-5.0))) - d2dphidt0_3;
}

__global__ void kernel_FreeEnergy_multi_dynamics4(int NX, int NY, int NZ, int Ntot, double dmuL, double dmuV, double *dPHIdt0, double *dPHIdt, double *dGLdmuL, double *dGSdmuL, double *dGVdmuL, double *dGLdmuV, double *dGSdmuV, double *dGVdmuV)
{
    int i, j, k, ind;
    i = blockIdx.x; j = blockIdx.y; k = threadIdx.x;
    ind = i*NY*NZ+j*NZ+k;
    
    dPHIdt[ind] = dPHIdt0[ind] + dmuL*dGLdmuL[ind] + dmuV*dGLdmuV[ind];
    
    dPHIdt[Ntot+ind] = dPHIdt0[Ntot+ind] + dmuL*dGSdmuL[ind] + dmuV*dGSdmuV[ind];
    
    dPHIdt[2*Ntot+ind] = dPHIdt0[2*Ntot+ind] + dmuL*dGVdmuL[ind] + dmuV*dGVdmuV[ind];
}

__global__ void kernel_FreeEnergy_multi_dynamics5(int NX, int NY, int NZ, int Ntot, double dmuL, double dmuV, double dfY, double *dPHIdt0, double *dPHIdt, double *dGLdmuL, double *dGSdmuL, double *dGVdmuL, double *dGLdmuV, double *dGSdmuV, double *dGVdmuV, double *matrix_y)
{
    int i, j, k, ind;
    i = blockIdx.x; j = blockIdx.y; k = threadIdx.x;
    ind = i*NY*NZ+j*NZ+k;
    
    dPHIdt[ind] = dPHIdt0[ind] + dmuL*dGLdmuL[ind] + dmuV*dGLdmuV[ind] + dfY*dGLdmuL[ind]*matrix_y[ind];
    
    dPHIdt[Ntot+ind] = dPHIdt0[Ntot+ind] + dmuL*dGSdmuL[ind] + dmuV*dGSdmuV[ind] + dfY*dGSdmuL[ind]*matrix_y[ind];
    
    dPHIdt[2*Ntot+ind] = dPHIdt0[2*Ntot+ind] + dmuL*dGVdmuL[ind] + dmuV*dGVdmuV[ind] + dfY*dGVdmuL[ind]*matrix_y[ind];
}

__global__ void kernel_FreeEnergy_multi_dynamics6(int NX, int NY, int NZ, int Ntot, double dmuL, double dmuV, double dfX, double *dPHIdt0, double *dPHIdt, double *dGLdmuL, double *dGSdmuL, double *dGVdmuL, double *dGLdmuV, double *dGSdmuV, double *dGVdmuV, double *matrix_x)
{
    int i, j, k, ind;
    i = blockIdx.x; j = blockIdx.y; k = threadIdx.x;
    ind = i*NY*NZ+j*NZ+k;
    
    dPHIdt[ind] = dPHIdt0[ind] + dmuL*dGLdmuL[ind] + dmuV*dGLdmuV[ind] + dfX*dGLdmuL[ind]*matrix_x[ind];
    
    dPHIdt[Ntot+ind] = dPHIdt0[Ntot+ind] + dmuL*dGSdmuL[ind] + dmuV*dGSdmuV[ind] + dfX*dGSdmuL[ind]*matrix_x[ind];
    
    dPHIdt[2*Ntot+ind] = dPHIdt0[2*Ntot+ind] + dmuL*dGVdmuL[ind] + dmuV*dGVdmuV[ind] + dfX*dGVdmuL[ind]*matrix_x[ind];
}

__global__ void kernel_FreeEnergy_multi_dynamics7(int NX, int NY, int NZ, int Ntot, double dmuL, double dmuV, double dfX, double dfY, double *dPHIdt0, double *dPHIdt, double *dGLdmuL, double *dGSdmuL, double *dGVdmuL, double *dGLdmuV, double *dGSdmuV, double *dGVdmuV, double *matrix_x, double *matrix_y)
{
    int i, j, k, ind;
    i = blockIdx.x; j = blockIdx.y; k = threadIdx.x;
    ind = i*NY*NZ+j*NZ+k;
    
    dPHIdt[ind] = dPHIdt0[ind] + dmuL*dGLdmuL[ind] + dmuV*dGLdmuV[ind]
                + dfX*dGLdmuL[ind]*matrix_x[ind] + dfY*dGLdmuL[ind]*matrix_y[ind];
    
    dPHIdt[Ntot+ind] = dPHIdt0[Ntot+ind] + dmuL*dGSdmuL[ind] + dmuV*dGSdmuV[ind]
                     + dfX*dGSdmuL[ind]*matrix_x[ind] + dfY*dGSdmuL[ind]*matrix_y[ind];
    
    dPHIdt[2*Ntot+ind] = dPHIdt0[2*Ntot+ind] + dmuL*dGVdmuL[ind] + dmuV*dGVdmuV[ind]
                       + dfX*dGVdmuL[ind]*matrix_x[ind] + dfY*dGVdmuL[ind]*matrix_y[ind];
}

__global__ void kernel_FreeEnergy_multi_dynamics8(int NX, int NY, int NZ, int Ntot, double dmuL, double dfX, double dfY, double *dPHIdt0, double *dPHIdt, double *dGLdmuL, double *dGSdmuL, double *dGVdmuL, double *matrix_x, double *matrix_y)
{
    int i, j, k, ind;
    i = blockIdx.x; j = blockIdx.y; k = threadIdx.x;
    ind = i*NY*NZ+j*NZ+k;
    
    dPHIdt[ind] = dPHIdt0[ind] + dmuL*dGLdmuL[ind]
    + dfX*dGLdmuL[ind]*matrix_x[ind] + dfY*dGLdmuL[ind]*matrix_y[ind];
    
    dPHIdt[Ntot+ind] = dPHIdt0[Ntot+ind] + dmuL*dGSdmuL[ind]
    + dfX*dGSdmuL[ind]*matrix_x[ind] + dfY*dGSdmuL[ind]*matrix_y[ind];
    
    dPHIdt[2*Ntot+ind] = dPHIdt0[2*Ntot+ind] + dmuL*dGVdmuL[ind]
    + dfX*dGVdmuL[ind]*matrix_x[ind] + dfY*dGVdmuL[ind]*matrix_y[ind];
}

__global__ void kernel_FreeEnergy_multi_Euler(int NX, int NY, int NZ, int Ntot, int num_fields, double timestep,    double *PHI, double *dPHIdt)
{
    int i, j, k, n, ind;
    double phi_tmp;
    i = blockIdx.x; j = blockIdx.y; k = threadIdx.x;
    ind = i*NY*NZ+j*NZ+k;
    
    for(n = 0; n < num_fields; n++)
    {
        phi_tmp = PHI[n*Ntot + ind];
        phi_tmp += dPHIdt[n*Ntot + ind]*timestep;
        PHI[n*Ntot + ind] = phi_tmp;
    }
}

#ifdef _USEFLOAT
float PhaseFieldFrame::cuda_FreeEnergy_multi()
#else
double PhaseFieldFrame::cuda_FreeEnergy_multi()
#endif
{

int  i, j, n;
double dmm, dmu1;
double vol_incr, x_incr, y_incr;
double A_11, A_12, A_13, A_14, A_21, A_22, A_23, A_24, A_31, A_32, A_33, A_34, A_41, A_42, A_43, A_44;
double B_1, B_2, B_3, B_4;
double dmuL, dmuV, dfX, dfY, divider;

if (gridsize<=0)
        FATAL("FreeEnergy: gridesize = "<<gridsize);
    
if (num_fields!=3)
        FATAL("FreeEnergy: num_field = "<<num_fields);
    
    
    
/* Reset liquid chemical potential, will be updated by constraint */
if (dynamics_type == 2 || dynamics_type == 3 || dynamics_type == 8)
{
    _M[0][1] = 0; _M[0][2] = _M[1][2];
    _M[1][0] = 0; _M[2][0] = _M[2][1];
    _MU[0] = _MU[1];
}

if (dynamics_type == 4 || dynamics_type == 5 || dynamics_type == 6 || dynamics_type == 7)
{
    // _MU[1] = 0;
    _MU[0] = _MU[1];
    _MU[2] = _MU[1];
}


/* copy the PHI value from the host to device */
cudaMemcpy( d_PHI,  _PHI[0], _NX*_NY*_NZ*sizeof(double), cudaMemcpyHostToDevice);
cudaMemcpy( d_PHI + Ntot,  _PHI[1], _NX*_NY*_NZ*sizeof(double), cudaMemcpyHostToDevice);
cudaMemcpy( d_PHI + Ntot*2,  _PHI[2], _NX*_NY*_NZ*sizeof(double), cudaMemcpyHostToDevice);


dim3 blocks(_NX,_NY);
dim3 threads(_NZ);
    
kernel_FreeEnergy_multi_interface<<<blocks, threads>>>
(_NX, _NY, _NZ, Ntot, Mob_D, d_PHI, d_K_S, d_K_L, d_K_LS, d_K_LV, d_K_SV, d_K_other, d_K_D, d_CPHI);
    
thrust::device_ptr<double> t_CPHI = thrust::device_pointer_cast(d_CPHI);
liquid_volume = thrust::reduce(t_CPHI,t_CPHI+_NX*_NY*_NZ);

_F = 0; _Fraw = 0;

kernel_FreeEnergy_multi_freeenergy<<<blocks, threads>>>
(_NX, _NY, _NZ, Ntot, gridsize, num_fields, Fext_x, Fext_y, Penalty_3phase, Penalty_NW, d_PHI, d_fden, d_fden_raw, d_dFdPHIi, d_d2PHI, d_EPS_loc, d_dPHIdx, d_dPHIdy, d_dPHIdz, d_tmp_x, d_tmp_y, d_tmp_z, d_matrix_x, d_matrix_y, d_U, d_EPS, d_EPS1, d_EPS2, d_EPS3, d_MU, d_R, d_NW_orig);
    
/* Use Thrust for free energy reduction */
thrust::device_ptr<double> t_fden = thrust::device_pointer_cast(d_fden);
_F = thrust::reduce(t_fden,t_fden+_NX*_NY*_NZ);
_F *= CUBE(gridsize);
thrust::device_ptr<double> t_fden_raw = thrust::device_pointer_cast(d_fden_raw);
_Fraw = thrust::reduce(t_fden_raw,t_fden_raw+_NX*_NY*_NZ);
_Fraw *= CUBE(gridsize);
    
    
kernel_FreeEnergy_multi_freeenergy2<<<blocks, threads>>>
(_NX, _NY, _NZ, Ntot, num_fields, gridsize, Penalty_3phase, Penalty_NW, d_PHI, d_dFdPHIi, d_tmp_x, d_tmp_y, d_tmp_z, d_EPS_loc, d_dPHIdx, d_dPHIdy, d_dPHIdz, d_EPS, d_EPS1, d_EPS2, d_EPS3, d_NW_orig);

kernel_FreeEnergy_multi_EOM<<<blocks, threads>>>
(_NX, _NY, _NZ, Ntot, dynamics_type, Mob_GL, Mob_LV, Mob_LS, Mob_SV, d_dPHIdt0, d_dPHIdtCH, d_dFdPHIi, d_K_LS, d_K_LV, d_K_SV, d_K_other, d_K_D, d_abs_dPHIdt);

_G = 0;
    
thrust::device_ptr<double> t_abs_dPHIdt = thrust::device_pointer_cast(d_abs_dPHIdt);
thrust::device_ptr<double> max_ptr = thrust::max_element(t_abs_dPHIdt,t_abs_dPHIdt+num_fields*_NX*_NY*_NZ);
_G = max_ptr[0];
    
timestep = dtmax/fmax(_G*dtmax/dphimax,1);

vol_incr = vol_incr0/timestep;
x_incr = x_incr0*liquid_volume/timestep;
y_incr = y_incr0*liquid_volume/timestep;
    
/* Equation of motions part */
    
if (dynamics_type == 0) /* no additional constraints */
{
/* 0: Ginzburg-Landau equation: not conserved */
    kernel_FreeEnergy_multi_dynamics0<<<blocks, threads>>>
    (_NX, _NY, _NZ, Ntot, num_fields, d_dPHIdt0, d_dPHIdt);
}
else if (dynamics_type == 1) /* no additional constraints */
{
    kernel_FreeEnergy_multi_dynamics1<<<blocks, threads>>>
    (_NX, _NY, _NZ, Ntot, num_fields, gridsize, d_dPHIdtCH, d_dPHIdt);
}
else if (dynamics_type == 2)
{

    kernel_FreeEnergy_multi_mprime<<<blocks, threads>>>
    (_NX, _NY, _NZ, Ntot, Mob_GL, Mob_LS, Mob_SV, Mob_LV, d_PHI, d_mprime, d_K_LS, d_K_SV, d_K_LV, d_K_other, d_Cprime0);


    /* Update DM */
    dmm = 0; dmu1 = 0;
        
    thrust::device_ptr<double> t_mprime = thrust::device_pointer_cast(d_mprime);
    thrust::device_ptr<double> t_Cprime0 = thrust::device_pointer_cast(d_Cprime0);
    thrust::device_ptr<double> t_dPHIdt0 = thrust::device_pointer_cast(d_dPHIdt0);
    dmm = thrust::inner_product(t_mprime, t_mprime+_NX*_NY*_NZ, t_Cprime0, 0.0);
    dmu1 = thrust::inner_product(t_Cprime0, t_Cprime0+_NX*_NY*_NZ, t_dPHIdt0, 0.0);
    dmu1 = -dmu1/dmm;
    
    kernel_FreeEnergy_multi_dynamics2<<<blocks, threads>>>
    (_NX, _NY, _NZ, Ntot, Mob_GL, Mob_LS, Mob_SV, Mob_LV, dmu1, d_PHI, d_dPHIdt0, d_dPHIdt, d_mprime, d_K_LS, d_K_SV, d_K_LV, d_K_other);

    /* update chemical potential */
    _MU[0] = _MU[1] + dmu1;

    for (i=0;i<3;i++)
        for (j=0;j<3;j++)
            _M[i][j] = _MU[i]-_MU[j];

}
else if (dynamics_type == 3)
{

    kernel_FreeEnergy_multi_mprime<<<blocks, threads>>>
    (_NX, _NY, _NZ, Ntot, Mob_GL, Mob_LS, Mob_SV, Mob_LV, d_PHI, d_mprime, d_K_LS, d_K_SV, d_K_LV, d_K_other, d_Cprime0);
    
    
    /* Update DM */
    dmm = 0; dmu1 = 0;
    
    thrust::device_ptr<double> t_mprime = thrust::device_pointer_cast(d_mprime);
    thrust::device_ptr<double> t_Cprime0 = thrust::device_pointer_cast(d_Cprime0);
    thrust::device_ptr<double> t_dPHIdt0 = thrust::device_pointer_cast(d_dPHIdt0);
    dmm = thrust::inner_product(t_mprime, t_mprime+_NX*_NY*_NZ, t_Cprime0, 0.0);
    dmu1 = thrust::inner_product(t_Cprime0, t_Cprime0+_NX*_NY*_NZ, t_dPHIdt0, 0.0);
    dmu1 = -dmu1/dmm;
    
    kernel_FreeEnergy_multi_dynamics3<<<blocks, threads>>>
    (_NX, _NY, _NZ, Ntot, gridsize, Mob_GL, Mob_LS, Mob_SV, Mob_LV, dmu1, d_PHI, d_dPHIdt0, d_dPHIdt, d_dPHIdtCH, d_mprime, d_K_LS, d_K_SV, d_K_LV, d_K_other);
    /* update chemical potential */
    _MU[0] = _MU[1] + dmu1;

    for (i=0;i<3;i++)
        for (j=0;j<3;j++)
            _M[i][j] = _MU[i]-_MU[j];
}
else if (dynamics_type == 4)
{
    _MU[1] = 0;
    
    kernel_FreeEnergy_multi_dFdmu<<<blocks, threads>>>
    (_NX, _NY, _NZ, Ntot, Mob_GL, Mob_LS, Mob_SV, Mob_LV, d_PHI, d_Cprime0, d_Cprime2, d_K_LS, d_K_LV, d_K_SV, d_K_other, d_dGLdmuL, d_dGLdmuV, d_dGVdmuL, d_dGVdmuV, d_dGSdmuL, d_dGSdmuV);
   
    /* initialize the parameters */
    A_11=0; A_12=0; A_21=0; A_22=0; B_1=0; B_2=0; dmuL=0; dmuV=0;
    
    
    thrust::device_ptr<double> t_dGLdmuL = thrust::device_pointer_cast(d_dGLdmuL);
    thrust::device_ptr<double> t_dGLdmuV = thrust::device_pointer_cast(d_dGLdmuV);
    thrust::device_ptr<double> t_dGVdmuL = thrust::device_pointer_cast(d_dGVdmuL);
    thrust::device_ptr<double> t_dGVdmuV = thrust::device_pointer_cast(d_dGVdmuV);
    thrust::device_ptr<double> t_Cprime0 = thrust::device_pointer_cast(d_Cprime0);
    thrust::device_ptr<double> t_Cprime2 = thrust::device_pointer_cast(d_Cprime2);
    thrust::device_ptr<double> t_dPHIdt0 = thrust::device_pointer_cast(d_dPHIdt0);
    
    A_11 = thrust::inner_product(t_dGLdmuL, t_dGLdmuL+_NX*_NY*_NZ, t_Cprime0, 0.0);
    A_12 = thrust::inner_product(t_dGLdmuV, t_dGLdmuV+_NX*_NY*_NZ, t_Cprime0, 0.0);
    A_21 = thrust::inner_product(t_dGVdmuL, t_dGVdmuL+_NX*_NY*_NZ, t_Cprime2, 0.0);
    A_22 = thrust::inner_product(t_dGVdmuV, t_dGVdmuV+_NX*_NY*_NZ, t_Cprime2, 0.0);
    B_1 =  -thrust::inner_product(t_dPHIdt0, t_dPHIdt0+Ntot, t_Cprime0, 0.0);
    B_2 =  -thrust::inner_product(t_dPHIdt0+2*Ntot, t_dPHIdt0+3*Ntot, t_Cprime2, 0.0);
    B_2 += vol_incr;

    dmuL = (A_12*B_2 - A_22*B_1)/(A_12*A_21-A_11*A_22);
    dmuV = -(A_11*B_2 - A_21*B_1)/(A_12*A_21-A_11*A_22);
    
    kernel_FreeEnergy_multi_dynamics4<<<blocks, threads>>>
    (_NX, _NY, _NZ, Ntot, dmuL, dmuV, d_dPHIdt0, d_dPHIdt, d_dGLdmuL, d_dGSdmuL, d_dGVdmuL, d_dGLdmuV, d_dGSdmuV, d_dGVdmuV);
    
    /* update chemical potential */
    _MU[0] = _MU[1] + dmuL;
    _MU[2] = _MU[1] + dmuV;
    
    for (i=0;i<3;i++)
        for (j=0;j<3;j++)
            _M[i][j] = _MU[i]-_MU[j];
}
else if (dynamics_type == 5)
{
    _MU[1] = 0;
    
    kernel_FreeEnergy_multi_dFdmu<<<blocks, threads>>>
    (_NX, _NY, _NZ, Ntot, Mob_GL, Mob_LS, Mob_SV, Mob_LV, d_PHI, d_Cprime0, d_Cprime2, d_K_LS, d_K_LV, d_K_SV, d_K_other, d_dGLdmuL, d_dGLdmuV, d_dGVdmuL, d_dGVdmuV, d_dGSdmuL, d_dGSdmuV);
    
    /* initialize the parameters */
    A_11=0; A_12=0; A_13=0; A_21=0; A_22=0; A_23=0; A_31=0; A_32=0; A_33=0;
    B_1=0; B_2=0; B_3=0; dmuL=0; dmuV=0; dfY=0; divider=0;
    
    thrust::device_ptr<double> t_dGLdmuL = thrust::device_pointer_cast(d_dGLdmuL);
    thrust::device_ptr<double> t_dGLdmuV = thrust::device_pointer_cast(d_dGLdmuV);
    thrust::device_ptr<double> t_dGVdmuL = thrust::device_pointer_cast(d_dGVdmuL);
    thrust::device_ptr<double> t_dGVdmuV = thrust::device_pointer_cast(d_dGVdmuV);
    thrust::device_ptr<double> t_Cprime0 = thrust::device_pointer_cast(d_Cprime0);
    thrust::device_ptr<double> t_Cprime2 = thrust::device_pointer_cast(d_Cprime2);
    thrust::device_ptr<double> t_dPHIdt0 = thrust::device_pointer_cast(d_dPHIdt0);
    thrust::device_ptr<double> t_matrix_x = thrust::device_pointer_cast(d_matrix_x);
    thrust::device_ptr<double> t_matrix_y = thrust::device_pointer_cast(d_matrix_y);
    thrust::device_ptr<double> inter_array = thrust::device_malloc<double>(_NX*_NY*_NZ);
    thrust::device_ptr<double> inter_array2 = thrust::device_malloc<double>(_NX*_NY*_NZ);
    thrust::multiplies<double> multi_op;
    
    A_11 = thrust::inner_product(t_dGLdmuL, t_dGLdmuL+_NX*_NY*_NZ, t_Cprime0, 0.0);
    A_12 = thrust::inner_product(t_dGLdmuV, t_dGLdmuV+_NX*_NY*_NZ, t_Cprime0, 0.0);
    thrust::transform(t_matrix_y, t_matrix_y + _NX*_NY*_NZ, t_dGLdmuL, inter_array, multi_op);
    A_13 = thrust::inner_product(inter_array, inter_array+_NX*_NY*_NZ, t_Cprime0, 0.0);
    A_21 = thrust::inner_product(t_dGVdmuL, t_dGVdmuL+_NX*_NY*_NZ, t_Cprime2, 0.0);
    A_22 = thrust::inner_product(t_dGVdmuV, t_dGVdmuV+_NX*_NY*_NZ, t_Cprime2, 0.0);
    thrust::transform(t_matrix_y, t_matrix_y + _NX*_NY*_NZ, t_dGVdmuL, inter_array, multi_op);
    A_23 = thrust::inner_product(inter_array, inter_array+_NX*_NY*_NZ, t_Cprime2, 0.0);
    A_31 = A_13;
    thrust::transform(t_matrix_y, t_matrix_y + _NX*_NY*_NZ, t_dGLdmuV, inter_array, multi_op);
    A_32 = thrust::inner_product(inter_array, inter_array+_NX*_NY*_NZ, t_Cprime0, 0.0);
    thrust::transform(t_matrix_y, t_matrix_y + _NX*_NY*_NZ, t_matrix_y, inter_array, multi_op);
    thrust::transform(inter_array, inter_array + _NX*_NY*_NZ, t_dGLdmuL, inter_array2, multi_op);
    A_33 = thrust::inner_product(inter_array2, inter_array2+_NX*_NY*_NZ, t_Cprime0, 0.0);
    B_1 = -thrust::inner_product(t_Cprime0, t_Cprime0+_NX*_NY*_NZ, t_dPHIdt0, 0.0);
    B_2 = -thrust::inner_product(t_Cprime2, t_Cprime2+_NX*_NY*_NZ, t_dPHIdt0+2*Ntot, 0.0);
    thrust::transform(t_matrix_y, t_matrix_y + _NX*_NY*_NZ, t_Cprime0, inter_array, multi_op);
    B_3 = -thrust::inner_product(inter_array, inter_array+_NX*_NY*_NZ, t_dPHIdt0, 0.0);
    
    B_2 += vol_incr;
    B_3 += y_incr;
    
    divider = A_11*A_22*A_33-A_11*A_23*A_32-A_12*A_21*A_33+A_12*A_23*A_31+A_13*A_21*A_32-A_13*A_22*A_31;
    dmuL =  (A_12*A_23*B_3-A_12*A_33*B_2-A_13*A_22*B_3+A_13*A_32*B_2+A_22*A_33*B_1-A_23*A_32*B_1)/divider;
    dmuV = -(A_11*A_23*B_3-A_11*A_33*B_2-A_13*A_21*B_3+A_13*A_31*B_2+A_21*A_33*B_1-A_23*A_31*B_1)/divider;
    dfY  =  (A_11*A_22*B_3-A_11*A_32*B_2-A_12*A_21*B_3+A_12*A_31*B_2+A_21*A_32*B_1-A_22*A_31*B_1)/divider;
    
    thrust::device_free(inter_array);
    thrust::device_free(inter_array2);
    
    kernel_FreeEnergy_multi_dynamics5<<<blocks, threads>>>
    (_NX, _NY, _NZ, Ntot, dmuL, dmuV, dfY, d_dPHIdt0, d_dPHIdt, d_dGLdmuL, d_dGSdmuL, d_dGVdmuL, d_dGLdmuV, d_dGSdmuV, d_dGVdmuV, d_matrix_y);
    
    /* update chemical potential */
    _MU[0] = _MU[1] + dmuL;
    _MU[2] = _MU[1] + dmuV;
    Fext_y = 0.0 + dfY;
    
    for (i=0;i<3;i++)
        for (j=0;j<3;j++)
            _M[i][j] = _MU[i]-_MU[j];
    
}
else if (dynamics_type == 6)
{
    _MU[1]=0;
    
    kernel_FreeEnergy_multi_dFdmu<<<blocks, threads>>>
    (_NX, _NY, _NZ, Ntot, Mob_GL, Mob_LS, Mob_SV, Mob_LV, d_PHI, d_Cprime0, d_Cprime2, d_K_LS, d_K_LV, d_K_SV, d_K_other, d_dGLdmuL, d_dGLdmuV, d_dGVdmuL, d_dGVdmuV, d_dGSdmuL, d_dGSdmuV);
    
    /* initialize the parameters */
    A_11=0; A_12=0; A_13=0; A_21=0; A_22=0; A_23=0; A_31=0; A_32=0; A_33=0;
    B_1=0; B_2=0; B_3=0; dmuL=0; dmuV=0; dfX=0; divider=0; _MU[1]=0;
    
    thrust::device_ptr<double> t_dGLdmuL = thrust::device_pointer_cast(d_dGLdmuL);
    thrust::device_ptr<double> t_dGLdmuV = thrust::device_pointer_cast(d_dGLdmuV);
    thrust::device_ptr<double> t_dGVdmuL = thrust::device_pointer_cast(d_dGVdmuL);
    thrust::device_ptr<double> t_dGVdmuV = thrust::device_pointer_cast(d_dGVdmuV);
    thrust::device_ptr<double> t_Cprime0 = thrust::device_pointer_cast(d_Cprime0);
    thrust::device_ptr<double> t_Cprime2 = thrust::device_pointer_cast(d_Cprime2);
    thrust::device_ptr<double> t_dPHIdt0 = thrust::device_pointer_cast(d_dPHIdt0);
    thrust::device_ptr<double> t_matrix_x = thrust::device_pointer_cast(d_matrix_x);
    thrust::device_ptr<double> t_matrix_y = thrust::device_pointer_cast(d_matrix_y);
    thrust::device_ptr<double> inter_array = thrust::device_malloc<double>(_NX*_NY*_NZ);
    thrust::device_ptr<double> inter_array2 = thrust::device_malloc<double>(_NX*_NY*_NZ);
    thrust::multiplies<double> multi_op;
    
    A_11 = thrust::inner_product(t_dGLdmuL, t_dGLdmuL+_NX*_NY*_NZ, t_Cprime0, 0.0);
    A_12 = thrust::inner_product(t_dGLdmuV, t_dGLdmuV+_NX*_NY*_NZ, t_Cprime0, 0.0);
    thrust::transform(t_matrix_x, t_matrix_x + _NX*_NY*_NZ, t_dGLdmuL, inter_array, multi_op);
    A_13 = thrust::inner_product(inter_array, inter_array+_NX*_NY*_NZ, t_Cprime0, 0.0);
    A_21 = thrust::inner_product(t_dGVdmuL, t_dGVdmuL+_NX*_NY*_NZ, t_Cprime2, 0.0);
    A_22 = thrust::inner_product(t_dGVdmuV, t_dGVdmuV+_NX*_NY*_NZ, t_Cprime2, 0.0);
    thrust::transform(t_matrix_x, t_matrix_x + _NX*_NY*_NZ, t_dGVdmuL, inter_array, multi_op);
    A_23 = thrust::inner_product(inter_array, inter_array+_NX*_NY*_NZ, t_Cprime2, 0.0);
    A_31 = A_13;
    thrust::transform(t_matrix_x, t_matrix_x + _NX*_NY*_NZ, t_dGLdmuV, inter_array, multi_op);
    A_32 = thrust::inner_product(inter_array, inter_array+_NX*_NY*_NZ, t_Cprime0, 0.0);
    thrust::transform(t_matrix_x, t_matrix_x + _NX*_NY*_NZ, t_matrix_x, inter_array, multi_op);
    thrust::transform(inter_array, inter_array + _NX*_NY*_NZ, t_dGLdmuL, inter_array2, multi_op);
    A_33 = thrust::inner_product(inter_array2, inter_array2+_NX*_NY*_NZ, t_Cprime0, 0.0);
    B_1 = -thrust::inner_product(t_Cprime0, t_Cprime0+_NX*_NY*_NZ, t_dPHIdt0, 0.0);
    B_2 = -thrust::inner_product(t_Cprime2, t_Cprime2+_NX*_NY*_NZ, t_dPHIdt0+2*Ntot, 0.0);
    thrust::transform(t_matrix_x, t_matrix_x + _NX*_NY*_NZ, t_Cprime0, inter_array, multi_op);
    B_3 = -thrust::inner_product(inter_array, inter_array+_NX*_NY*_NZ, t_dPHIdt0, 0.0);
    
    B_2 += vol_incr;
    B_3 += x_incr;
    
    divider = A_11*A_22*A_33-A_11*A_23*A_32-A_12*A_21*A_33+A_12*A_23*A_31+A_13*A_21*A_32-A_13*A_22*A_31;
    dmuL =  (A_12*A_23*B_3-A_12*A_33*B_2-A_13*A_22*B_3+A_13*A_32*B_2+A_22*A_33*B_1-A_23*A_32*B_1)/divider;
    dmuV = -(A_11*A_23*B_3-A_11*A_33*B_2-A_13*A_21*B_3+A_13*A_31*B_2+A_21*A_33*B_1-A_23*A_31*B_1)/divider;
    dfX  =  (A_11*A_22*B_3-A_11*A_32*B_2-A_12*A_21*B_3+A_12*A_31*B_2+A_21*A_32*B_1-A_22*A_31*B_1)/divider;
    
    thrust::device_free(inter_array);
    thrust::device_free(inter_array2);
    
    kernel_FreeEnergy_multi_dynamics6<<<blocks, threads>>>
    (_NX, _NY, _NZ, Ntot, dmuL, dmuV, dfX, d_dPHIdt0, d_dPHIdt, d_dGLdmuL, d_dGSdmuL, d_dGVdmuL, d_dGLdmuV, d_dGSdmuV, d_dGVdmuV, d_matrix_x);
    
    /* update chemical potential */
    _MU[0] = _MU[1] + dmuL;
    _MU[2] = _MU[1] + dmuV;
    Fext_x = 0.0 + dfX;
    
    for (i=0;i<3;i++)
        for (j=0;j<3;j++)
            _M[i][j] = _MU[i]-_MU[j];
    
}
else if (dynamics_type == 7)
{
    kernel_FreeEnergy_multi_dFdmu<<<blocks, threads>>>
    (_NX, _NY, _NZ, Ntot, Mob_GL, Mob_LS, Mob_SV, Mob_LV, d_PHI, d_Cprime0, d_Cprime2, d_K_LS, d_K_LV, d_K_SV, d_K_other, d_dGLdmuL, d_dGLdmuV, d_dGVdmuL, d_dGVdmuV, d_dGSdmuL, d_dGSdmuV);
    
    /* initialize the parameters */
    A_11=0; A_12=0; A_13=0; A_14=0; A_21=0; A_22=0; A_23=0; A_24=0; A_31=0; A_32=0; A_33=0; A_34=0; A_41=0; A_42=0; A_43=0; A_44=0;
    B_1=0; B_2=0; B_3=0; B_4=0; dmuL=0; dmuV=0; dfX=0; dfY=0; divider=0; _MU[1]=0;
    
    thrust::device_ptr<double> t_dGLdmuL = thrust::device_pointer_cast(d_dGLdmuL);
    thrust::device_ptr<double> t_dGLdmuV = thrust::device_pointer_cast(d_dGLdmuV);
    thrust::device_ptr<double> t_dGVdmuL = thrust::device_pointer_cast(d_dGVdmuL);
    thrust::device_ptr<double> t_dGVdmuV = thrust::device_pointer_cast(d_dGVdmuV);
    thrust::device_ptr<double> t_Cprime0 = thrust::device_pointer_cast(d_Cprime0);
    thrust::device_ptr<double> t_Cprime2 = thrust::device_pointer_cast(d_Cprime2);
    thrust::device_ptr<double> t_dPHIdt0 = thrust::device_pointer_cast(d_dPHIdt0);
    thrust::device_ptr<double> t_matrix_x = thrust::device_pointer_cast(d_matrix_x);
    thrust::device_ptr<double> t_matrix_y = thrust::device_pointer_cast(d_matrix_y);
    thrust::device_ptr<double> inter_array = thrust::device_malloc<double>(_NX*_NY*_NZ);
    thrust::device_ptr<double> inter_array2 = thrust::device_malloc<double>(_NX*_NY*_NZ);
    thrust::multiplies<double> multi_op;
    
    A_11 = thrust::inner_product(t_dGLdmuL, t_dGLdmuL+_NX*_NY*_NZ, t_Cprime0, 0.0);
    A_12 = thrust::inner_product(t_dGLdmuV, t_dGLdmuV+_NX*_NY*_NZ, t_Cprime0, 0.0);
    thrust::transform(t_matrix_x, t_matrix_x + _NX*_NY*_NZ, t_dGLdmuL, inter_array, multi_op);
    A_13 = thrust::inner_product(inter_array, inter_array+_NX*_NY*_NZ, t_Cprime0, 0.0);
    thrust::transform(t_matrix_y, t_matrix_y + _NX*_NY*_NZ, t_dGLdmuL, inter_array, multi_op);
    A_14 = thrust::inner_product(inter_array, inter_array+_NX*_NY*_NZ, t_Cprime0, 0.0);
    A_21 = thrust::inner_product(t_dGVdmuL, t_dGVdmuL+_NX*_NY*_NZ, t_Cprime2, 0.0);
    A_22 = thrust::inner_product(t_dGVdmuV, t_dGVdmuV+_NX*_NY*_NZ, t_Cprime2, 0.0);
    thrust::transform(t_matrix_x, t_matrix_x + _NX*_NY*_NZ, t_dGVdmuL, inter_array, multi_op);
    A_23 = thrust::inner_product(inter_array, inter_array+_NX*_NY*_NZ, t_Cprime2, 0.0);
    thrust::transform(t_matrix_y, t_matrix_y + _NX*_NY*_NZ, t_dGVdmuL, inter_array, multi_op);
    A_24 = thrust::inner_product(inter_array, inter_array+_NX*_NY*_NZ, t_Cprime2, 0.0);
    A_31 = A_13;
    thrust::transform(t_matrix_x, t_matrix_x + _NX*_NY*_NZ, t_dGLdmuV, inter_array, multi_op);
    A_32 = thrust::inner_product(inter_array, inter_array+_NX*_NY*_NZ, t_Cprime0, 0.0);
    thrust::transform(t_matrix_x, t_matrix_x + _NX*_NY*_NZ, t_matrix_x, inter_array, multi_op);
    thrust::transform(inter_array, inter_array + _NX*_NY*_NZ, t_dGLdmuL, inter_array2, multi_op);
    A_33 = thrust::inner_product(inter_array2, inter_array2+_NX*_NY*_NZ, t_Cprime0, 0.0);
    thrust::transform(t_matrix_x, t_matrix_x + _NX*_NY*_NZ, t_matrix_y, inter_array, multi_op);
    thrust::transform(inter_array, inter_array + _NX*_NY*_NZ, t_dGLdmuL, inter_array2, multi_op);
    A_34 = thrust::inner_product(inter_array2, inter_array2+_NX*_NY*_NZ, t_Cprime0, 0.0);
    A_41 = A_14;
    thrust::transform(t_matrix_y, t_matrix_y + _NX*_NY*_NZ, t_dGLdmuV, inter_array, multi_op);
    A_42 = thrust::inner_product(inter_array, inter_array+_NX*_NY*_NZ, t_Cprime0, 0.0);
    A_43 = A_34;
    thrust::transform(t_matrix_y, t_matrix_y + _NX*_NY*_NZ, t_matrix_y, inter_array, multi_op);
    thrust::transform(inter_array, inter_array + _NX*_NY*_NZ, t_dGLdmuL, inter_array2, multi_op);
    A_44 = thrust::inner_product(inter_array2, inter_array2+_NX*_NY*_NZ, t_Cprime0, 0.0);
    
    B_1 = -thrust::inner_product(t_Cprime0, t_Cprime0+_NX*_NY*_NZ, t_dPHIdt0, 0.0);
    B_2 = -thrust::inner_product(t_Cprime2, t_Cprime2+_NX*_NY*_NZ, t_dPHIdt0+2*Ntot, 0.0);
    thrust::transform(t_matrix_x, t_matrix_x + _NX*_NY*_NZ, t_Cprime0, inter_array, multi_op);
    B_3 = -thrust::inner_product(inter_array, inter_array+_NX*_NY*_NZ, t_dPHIdt0, 0.0);
    thrust::transform(t_matrix_y, t_matrix_y + _NX*_NY*_NZ, t_Cprime0, inter_array, multi_op);
    B_4 = -thrust::inner_product(inter_array, inter_array+_NX*_NY*_NZ, t_dPHIdt0, 0.0);
    
    B_2 += vol_incr;
    B_3 += x_incr;
    B_4 += y_incr;
    
    divider = A_11*A_22*A_33*A_44-A_11*A_22*A_34*A_43-A_11*A_23*A_32*A_44+A_11*A_23*A_34*A_42+A_11*A_24*A_32*A_43
    - A_11*A_24*A_33*A_42-A_12*A_21*A_33*A_44+A_12*A_21*A_34*A_43+A_12*A_23*A_31*A_44-A_12*A_23*A_34*A_41
    - A_12*A_24*A_31*A_43+A_12*A_24*A_33*A_41+A_13*A_21*A_32*A_44-A_13*A_21*A_34*A_42-A_13*A_22*A_31*A_44
    + A_13*A_22*A_34*A_41+A_13*A_24*A_31*A_42-A_13*A_24*A_32*A_41-A_14*A_21*A_32*A_43+A_14*A_21*A_33*A_42
    + A_14*A_22*A_31*A_43-A_14*A_22*A_33*A_41-A_14*A_23*A_31*A_42+A_14*A_23*A_32*A_41;
    dmuL = -(A_12*A_23*A_34*B_4-A_12*A_23*A_44*B_3-A_12*A_24*A_33*B_4+A_12*A_24*A_43*B_3
             +A_12*A_33*A_44*B_2-A_12*A_34*A_43*B_2-A_13*A_22*A_34*B_4+A_13*A_22*A_44*B_3
             +A_13*A_24*A_32*B_4-A_13*A_24*A_42*B_3-A_13*A_32*A_44*B_2+A_13*A_34*A_42*B_2
             +A_14*A_22*A_33*B_4-A_14*A_22*A_43*B_3-A_14*A_23*A_32*B_4+A_14*A_23*A_42*B_3
             +A_14*A_32*A_43*B_2-A_14*A_33*A_42*B_2-A_22*A_33*A_44*B_1+A_22*A_34*A_43*B_1
             +A_23*A_32*A_44*B_1-A_23*A_34*A_42*B_1-A_24*A_32*A_43*B_1+A_24*A_33*A_42*B_1)/divider;
    dmuV =  (A_11*A_23*A_34*B_4-A_11*A_23*A_44*B_3-A_11*A_24*A_33*B_4+A_11*A_24*A_43*B_3
             +A_11*A_33*A_44*B_2-A_11*A_34*A_43*B_2-A_13*A_21*A_34*B_4+A_13*A_21*A_44*B_3
             +A_13*A_24*A_31*B_4-A_13*A_24*A_41*B_3-A_13*A_31*A_44*B_2+A_13*A_34*A_41*B_2
             +A_14*A_21*A_33*B_4-A_14*A_21*A_43*B_3-A_14*A_23*A_31*B_4+A_14*A_23*A_41*B_3
             +A_14*A_31*A_43*B_2-A_14*A_33*A_41*B_2-A_21*A_33*A_44*B_1+A_21*A_34*A_43*B_1
             +A_23*A_31*A_44*B_1-A_23*A_34*A_41*B_1-A_24*A_31*A_43*B_1+A_24*A_33*A_41*B_1)/divider;
    dfX  = -(A_11*A_22*A_34*B_4-A_11*A_22*A_44*B_3-A_11*A_24*A_32*B_4+A_11*A_24*A_42*B_3
             +A_11*A_32*A_44*B_2-A_11*A_34*A_42*B_2-A_12*A_21*A_34*B_4+A_12*A_21*A_44*B_3
             +A_12*A_24*A_31*B_4-A_12*A_24*A_41*B_3-A_12*A_31*A_44*B_2+A_12*A_34*A_41*B_2
             +A_14*A_21*A_32*B_4-A_14*A_21*A_42*B_3-A_14*A_22*A_31*B_4+A_14*A_22*A_41*B_3
             +A_14*A_31*A_42*B_2-A_14*A_32*A_41*B_2-A_21*A_32*A_44*B_1+A_21*A_34*A_42*B_1
             +A_22*A_31*A_44*B_1-A_22*A_34*A_41*B_1-A_24*A_31*A_42*B_1+A_24*A_32*A_41*B_1)/divider;
    dfY  =  (A_11*A_22*A_33*B_4-A_11*A_22*A_43*B_3-A_11*A_23*A_32*B_4+A_11*A_23*A_42*B_3
             +A_11*A_32*A_43*B_2-A_11*A_33*A_42*B_2-A_12*A_21*A_33*B_4+A_12*A_21*A_43*B_3
             +A_12*A_23*A_31*B_4-A_12*A_23*A_41*B_3-A_12*A_31*A_43*B_2+A_12*A_33*A_41*B_2
             +A_13*A_21*A_32*B_4-A_13*A_21*A_42*B_3-A_13*A_22*A_31*B_4+A_13*A_22*A_41*B_3
             +A_13*A_31*A_42*B_2-A_13*A_32*A_41*B_2-A_21*A_32*A_43*B_1+A_21*A_33*A_42*B_1
             +A_22*A_31*A_43*B_1-A_22*A_33*A_41*B_1-A_23*A_31*A_42*B_1+A_23*A_32*A_41*B_1)/divider;
    
    thrust::device_free(inter_array);
    thrust::device_free(inter_array2);
    
    kernel_FreeEnergy_multi_dynamics7<<<blocks, threads>>>
    (_NX, _NY, _NZ, Ntot, dmuL, dmuV, dfX, dfY, d_dPHIdt0, d_dPHIdt, d_dGLdmuL, d_dGSdmuL, d_dGVdmuL, d_dGLdmuV, d_dGSdmuV, d_dGVdmuV, d_matrix_x, d_matrix_y);
    
    /* update chemical potential */
    _MU[0] = _MU[1] + dmuL;
    _MU[2] = _MU[1] + dmuV;
    Fext_x = 0.0 + dfX;
    Fext_y = 0.0 + dfY;
    
    for (i=0;i<3;i++)
        for (j=0;j<3;j++)
            _M[i][j] = _MU[i]-_MU[j];
    
}
else if (dynamics_type == 8)
{
    kernel_FreeEnergy_multi_dFdmu<<<blocks, threads>>>
    (_NX, _NY, _NZ, Ntot, Mob_GL, Mob_LS, Mob_SV, Mob_LV, d_PHI, d_Cprime0, d_Cprime2, d_K_LS, d_K_LV, d_K_SV, d_K_other, d_dGLdmuL, d_dGLdmuV, d_dGVdmuL, d_dGVdmuV, d_dGSdmuL, d_dGSdmuV);
    
    /* initialize the parameters */
    A_11=0; A_12=0; A_13=0; A_21=0; A_22=0; A_23=0; A_31=0; A_32=0; A_33=0;
    B_1=0; B_2=0; B_3=0; dmuL=0; dfX=0; dfY=0; divider=0; _MU[1]=0;
    
    thrust::device_ptr<double> t_dGLdmuL = thrust::device_pointer_cast(d_dGLdmuL);
    thrust::device_ptr<double> t_dGLdmuV = thrust::device_pointer_cast(d_dGLdmuV);
    thrust::device_ptr<double> t_dGVdmuL = thrust::device_pointer_cast(d_dGVdmuL);
    thrust::device_ptr<double> t_dGVdmuV = thrust::device_pointer_cast(d_dGVdmuV);
    thrust::device_ptr<double> t_Cprime0 = thrust::device_pointer_cast(d_Cprime0);
    thrust::device_ptr<double> t_Cprime2 = thrust::device_pointer_cast(d_Cprime2);
    thrust::device_ptr<double> t_dPHIdt0 = thrust::device_pointer_cast(d_dPHIdt0);
    thrust::device_ptr<double> t_matrix_x = thrust::device_pointer_cast(d_matrix_x);
    thrust::device_ptr<double> t_matrix_y = thrust::device_pointer_cast(d_matrix_y);
    thrust::device_ptr<double> inter_array = thrust::device_malloc<double>(_NX*_NY*_NZ);
    thrust::device_ptr<double> inter_array2 = thrust::device_malloc<double>(_NX*_NY*_NZ);
    thrust::multiplies<double> multi_op;
    
    A_11 = thrust::inner_product(t_dGLdmuL, t_dGLdmuL+_NX*_NY*_NZ, t_Cprime0, 0.0);
    thrust::transform(t_matrix_x, t_matrix_x + _NX*_NY*_NZ, t_dGLdmuL, inter_array, multi_op);
    A_12 = thrust::inner_product(inter_array, inter_array+_NX*_NY*_NZ, t_Cprime0, 0.0);
    thrust::transform(t_matrix_y, t_matrix_y + _NX*_NY*_NZ, t_dGLdmuL, inter_array, multi_op);
    A_13 = thrust::inner_product(inter_array, inter_array+_NX*_NY*_NZ, t_Cprime0, 0.0);
    A_21 = A_12;
    thrust::transform(t_matrix_x, t_matrix_x + _NX*_NY*_NZ, t_matrix_x, inter_array, multi_op);
    thrust::transform(inter_array, inter_array + _NX*_NY*_NZ, t_dGLdmuL, inter_array2, multi_op);
    A_22 = thrust::inner_product(inter_array2, inter_array2+_NX*_NY*_NZ, t_Cprime0, 0.0);
    thrust::transform(t_matrix_x, t_matrix_x + _NX*_NY*_NZ, t_matrix_y, inter_array, multi_op);
    thrust::transform(inter_array, inter_array + _NX*_NY*_NZ, t_dGLdmuL, inter_array2, multi_op);
    A_23 = thrust::inner_product(inter_array2, inter_array2+_NX*_NY*_NZ, t_Cprime0, 0.0);
    A_31 = A_13;
    A_32 = A_23;
    thrust::transform(t_matrix_y, t_matrix_y + _NX*_NY*_NZ, t_matrix_y, inter_array, multi_op);
    thrust::transform(inter_array, inter_array + _NX*_NY*_NZ, t_dGLdmuL, inter_array2, multi_op);
    A_33 = thrust::inner_product(inter_array2, inter_array2+_NX*_NY*_NZ, t_Cprime0, 0.0);

    B_1 = -thrust::inner_product(t_Cprime0, t_Cprime0+_NX*_NY*_NZ, t_dPHIdt0, 0.0);
    thrust::transform(t_matrix_x, t_matrix_x + _NX*_NY*_NZ, t_Cprime0, inter_array, multi_op);
    B_2 = -thrust::inner_product(inter_array, inter_array+_NX*_NY*_NZ, t_dPHIdt0, 0.0);
    thrust::transform(t_matrix_y, t_matrix_y + _NX*_NY*_NZ, t_Cprime0, inter_array, multi_op);
    B_3 = -thrust::inner_product(inter_array, inter_array+_NX*_NY*_NZ, t_dPHIdt0, 0.0);


    divider = A_11*A_22*A_33-A_11*A_23*A_32-A_12*A_21*A_33+A_12*A_23*A_31+A_13*A_21*A_32-A_13*A_22*A_31;
    dmuL =  (A_12*A_23*B_3-A_12*A_33*B_2-A_13*A_22*B_3+A_13*A_32*B_2+A_22*A_33*B_1-A_23*A_32*B_1)/divider;
    dfX  = -(A_11*A_23*B_3-A_11*A_33*B_2-A_13*A_21*B_3+A_13*A_31*B_2+A_21*A_33*B_1-A_23*A_31*B_1)/divider;
    dfY  =  (A_11*A_22*B_3-A_11*A_32*B_2-A_12*A_21*B_3+A_12*A_31*B_2+A_21*A_32*B_1-A_22*A_31*B_1)/divider;
    
    thrust::device_free(inter_array);
    thrust::device_free(inter_array2);
    
    kernel_FreeEnergy_multi_dynamics8<<<blocks, threads>>>
    (_NX, _NY, _NZ, Ntot, dmuL, dfX, dfY, d_dPHIdt0, d_dPHIdt, d_dGLdmuL, d_dGSdmuL, d_dGVdmuL, d_matrix_x, d_matrix_y);
    
    /* update chemical potential */
    _MU[0] = _MU[1] + dmuL;
    Fext_x = 0.0 + dfX;
    Fext_y = 0.0 + dfY;
    
    for (i=0;i<3;i++)
        for (j=0;j<3;j++)
            _M[i][j] = _MU[i]-_MU[j];

}
else
{
    ERROR("unknown dynamics_type = "<<dynamics_type);
}

thrust::device_ptr<double> t_PHI = thrust::device_pointer_cast(d_PHI);
thrust::device_ptr<double> t_matrix_x = thrust::device_pointer_cast(d_matrix_x);
thrust::device_ptr<double> t_matrix_y = thrust::device_pointer_cast(d_matrix_y);
for (n =0 ; n< num_fields; n++)
{
    phase_volume[n] = thrust::reduce(t_PHI+n*Ntot,t_PHI+(n+1)*Ntot);
    phase_volume[n] /= (_NX * _NY * _NZ);
}

     com_x = thrust::inner_product(t_matrix_x, t_matrix_x+_NX*_NY*_NZ, t_CPHI, 0.0);
     com_x /= liquid_volume;
     com_y = thrust::inner_product(t_matrix_y, t_matrix_y+_NX*_NY*_NZ, t_CPHI, 0.0);
     com_y /= liquid_volume;

/*
    cudaFree(d_PHI);
    cudaFree(d_d2PHI);
    cudaFree(d_dFdPHIi);
    cudaFree(d_dPHIdt);
    cudaFree(d_abs_dPHIdt);
    cudaFree(d_dPHIdt0);
    cudaFree(d_dPHIdtCH);
    cudaFree(d_tmp_x);
    cudaFree(d_tmp_y);
    cudaFree(d_tmp_z);
    cudaFree(d_EPS_loc);
    cudaFree(d_dPHIdx);
    cudaFree(d_dPHIdy);
    cudaFree(d_dPHIdz);
    cudaFree(d_U);
    cudaFree(d_EPS);
    cudaFree(d_EPS1);
    cudaFree(d_EPS2);
    cudaFree(d_EPS3);
    cudaFree(d_R);
    cudaFree(d_MU);
    cudaFree(d_matrix_x);
    cudaFree(d_matrix_y);
    cudaFree(d_fden);
    cudaFree(d_fden_raw);
    cudaFree(d_mprime);
    cudaFree(d_K_S);
    cudaFree(d_K_L);
    cudaFree(d_K_LS);
    cudaFree(d_K_LV);
    cudaFree(d_K_SV);
    cudaFree(d_K_D);
    cudaFree(d_K_other);
    cudaFree(d_NW_orig);
    cudaFree(d_dGLdmuL);
    cudaFree(d_dGLdmuV);
    cudaFree(d_dGVdmuL);
    cudaFree(d_dGVdmuV);
    cudaFree(d_dGSdmuL);
    cudaFree(d_dGSdmuV);
    cudaFree(d_Cprime0);
    cudaFree(d_Cprime2);
*/

return _F;


}

void PhaseFieldFrame::cuda_memory_alloc()
{
cudaMalloc( &d_PHI, num_fields*_NX*_NY*_NZ*sizeof(double) );
cudaMalloc( &d_dFdPHIi, num_fields*_NX*_NY*_NZ*sizeof(double) );
cudaMalloc( &d_dPHIdt0, num_fields*_NX*_NY*_NZ*sizeof(double) );
cudaMalloc( &d_dPHIdt, num_fields*_NX*_NY*_NZ*sizeof(double) );
cudaMalloc( &d_abs_dPHIdt, num_fields*_NX*_NY*_NZ*sizeof(double) );
cudaMalloc( &d_dPHIdtCH, num_fields*_NX*_NY*_NZ*sizeof(double) );
cudaMalloc( &d_d2PHI, num_fields*_NX*_NY*_NZ*sizeof(double) );
cudaMalloc( &d_tmp_x, num_fields*_NX*_NY*_NZ*sizeof(double) );
cudaMalloc( &d_tmp_y, num_fields*_NX*_NY*_NZ*sizeof(double) );
cudaMalloc( &d_tmp_z, num_fields*_NX*_NY*_NZ*sizeof(double) );
cudaMalloc( &d_EPS_loc, num_fields*_NX*_NY*_NZ*sizeof(double) );
cudaMalloc( &d_dPHIdx, num_fields*_NX*_NY*_NZ*sizeof(double) );
cudaMalloc( &d_dPHIdy, num_fields*_NX*_NY*_NZ*sizeof(double) );
cudaMalloc( &d_dPHIdz, num_fields*_NX*_NY*_NZ*sizeof(double) );

cudaMalloc( &d_CPHI,  _NX*_NY*_NZ*sizeof(double) );
cudaMalloc( &d_NW_orig,  _NX*_NY*_NZ*sizeof(double) );
cudaMalloc( &d_fden,  _NX*_NY*_NZ*sizeof(double) );
cudaMalloc( &d_fden_raw,  _NX*_NY*_NZ*sizeof(double) );
cudaMalloc( &d_mprime,  _NX*_NY*_NZ*sizeof(double) );
cudaMalloc( &d_matrix_x,  _NX*_NY*_NZ*sizeof(double) );
cudaMalloc( &d_matrix_y,  _NX*_NY*_NZ*sizeof(double) );
cudaMalloc( &d_K_S,  _NX*_NY*_NZ*sizeof(double) );
cudaMalloc( &d_K_L,  _NX*_NY*_NZ*sizeof(double) );
cudaMalloc( &d_K_LS,  _NX*_NY*_NZ*sizeof(double) );
cudaMalloc( &d_K_SV,  _NX*_NY*_NZ*sizeof(double) );
cudaMalloc( &d_K_LV,  _NX*_NY*_NZ*sizeof(double) );
cudaMalloc( &d_K_other,  _NX*_NY*_NZ*sizeof(double) );
cudaMalloc( &d_K_D,  _NX*_NY*_NZ*sizeof(double) );

cudaMalloc( &d_dGLdmuL,  _NX*_NY*_NZ*sizeof(double) );
cudaMalloc( &d_dGLdmuV,  _NX*_NY*_NZ*sizeof(double) );
cudaMalloc( &d_dGVdmuL,  _NX*_NY*_NZ*sizeof(double) );
cudaMalloc( &d_dGVdmuV,  _NX*_NY*_NZ*sizeof(double) );
cudaMalloc( &d_dGSdmuL,  _NX*_NY*_NZ*sizeof(double) );
cudaMalloc( &d_dGSdmuV,  _NX*_NY*_NZ*sizeof(double) );
cudaMalloc( &d_Cprime0,  _NX*_NY*_NZ*sizeof(double) );
cudaMalloc( &d_Cprime2,  _NX*_NY*_NZ*sizeof(double) );

cudaMalloc( &d_U, num_fields*num_fields*sizeof(double) );
cudaMalloc( &d_EPS, num_fields*num_fields*sizeof(double) );
cudaMalloc( &d_EPS1, num_fields*num_fields*sizeof(double) );
cudaMalloc( &d_EPS2, num_fields*num_fields*sizeof(double) );
cudaMalloc( &d_EPS3, num_fields*num_fields*sizeof(double) );
cudaMalloc( &d_MU, num_fields*sizeof(double) );
cudaMalloc( &d_R, 3*3*sizeof(double) );

}

void PhaseFieldFrame::cuda_init()
{
/* copy all the related energetic parameters from the host to device */
cudaMemcpy( d_NW_orig, _NW_orig, _NX*_NY*_NZ*sizeof(double), cudaMemcpyHostToDevice);
cudaMemcpy( d_U, _U, num_fields*num_fields*sizeof(double), cudaMemcpyHostToDevice);
cudaMemcpy( d_EPS, _EPS, num_fields*num_fields*sizeof(double), cudaMemcpyHostToDevice);
cudaMemcpy( d_EPS1, _EPS1, num_fields*num_fields*sizeof(double), cudaMemcpyHostToDevice);
cudaMemcpy( d_EPS2, _EPS2, num_fields*num_fields*sizeof(double), cudaMemcpyHostToDevice);
cudaMemcpy( d_EPS3, _EPS3, num_fields*num_fields*sizeof(double), cudaMemcpyHostToDevice);
cudaMemcpy( d_MU, _MU, num_fields*sizeof(double), cudaMemcpyHostToDevice);
cudaMemcpy( d_R,_ROT_MATRIX.element, 3*3*sizeof(double), cudaMemcpyHostToDevice);

}

void PhaseFieldFrame::cuda_integrator()
{
    dim3 blocks(_NX,_NY);
    dim3 threads(_NZ);
    kernel_FreeEnergy_multi_Euler<<<blocks, threads>>>
    (_NX, _NY, _NZ, Ntot, num_fields, timestep, d_PHI, d_dPHIdt);
    
    cudaMemcpy( _PHI[0], d_PHI, _NX*_NY*_NZ*sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy( _PHI[1], d_PHI + Ntot, _NX*_NY*_NZ*sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy( _PHI[2], d_PHI + 2*Ntot, _NX*_NY*_NZ*sizeof(double), cudaMemcpyDeviceToHost);
}


#endif //_USECUDA

