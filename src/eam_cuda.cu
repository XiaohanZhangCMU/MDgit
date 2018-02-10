#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <thrust/extrema.h>
#include <thrust/inner_product.h>
#include <thrust/transform.h>
#include <thrust/functional.h>
#include <thrust/device_malloc.h>
#include <thrust/device_free.h>

/*
 Todo:
 1) The surface tension is not calcualted! 
 -> if (SURFTEN==1 && (curstep%savepropfreq)==1) AddnvvtoPtPn(_SR[jpt],rij,rij,-fp);
 2) It only supports MD and relaxation calculations. MC can be impelemented in a similar way.
 */

#include "eam.h"

void EAMFrame::cuda_memory_alloc() {
  int size = _NP*allocmultiple;
  gpuErrchk(cudaMalloc(&_d_atpe3b,sizeof(double)*size));
  gpuErrchk(cudaMalloc(&_d_rhotot,sizeof(double)*size));
  gpuErrchk(cudaMalloc(&_d_embf,  sizeof(double)*size));
  gpuErrchk(cudaMalloc(&_d_embfp, sizeof(double)*size));
  gpuErrchk(cudaMalloc(&_d_nbst,  sizeof(int)*size));
  /* eam grid allocation */
  gpuErrchk(cudaMalloc(&_d_rho,   sizeof(double)*4*NGRID));
  gpuErrchk(cudaMalloc(&_d_rhop,  sizeof(double)*4*NGRID));
  gpuErrchk(cudaMalloc(&_d_phi,   sizeof(double)*2*NGRID));
  gpuErrchk(cudaMalloc(&_d_phip,  sizeof(double)*2*NGRID));
  gpuErrchk(cudaMalloc(&_d_phix,   sizeof(double)*NGRID));
  gpuErrchk(cudaMalloc(&_d_phipx,  sizeof(double)*NGRID));
  gpuErrchk(cudaMalloc(&_d_frho,   sizeof(double)*2*NGRID));
  gpuErrchk(cudaMalloc(&_d_frhop,  sizeof(double)*2*NGRID));
  gpuErrchk(cudaMalloc(&_d_rho_spline,   sizeof(double)*4*NGRID*4));
  gpuErrchk(cudaMalloc(&_d_phi_spline,   sizeof(double)*2*NGRID*4));
  gpuErrchk(cudaMalloc(&_d_phix_spline,   sizeof(double)*4*NGRID));
  gpuErrchk(cudaMalloc(&_d_frho_spline,   sizeof(double)*2*NGRID*4));
  gpuErrchk(cudaMalloc(&_d_rval,          sizeof(double)*NGRID));
  gpuErrchk(cudaMalloc(&_d_rhoval,        sizeof(double)*NGRID));
  /* md data allocation */
  gpuErrchk(cudaMalloc(&_d_H_element,      sizeof(double)*3*3));
  gpuErrchk(cudaMalloc(&_d_VIRIAL_element, sizeof(double)*3*3));
  gpuErrchk(cudaMalloc(&_d_EPOT_IND,   sizeof(double)*size));
  gpuErrchk(cudaMalloc(&_d_species,    sizeof(int)*size));
  /* size of nindex is obtained from md.cpp:NbrList_init. mx = size, mz = NNM,
     but we don`t need memory for pointers (shft2) since cuda wants a 1d array. 
   */
  gpuErrchk(cudaMalloc(&_d_nindex,     sizeof(int)*size*_NNM));
  gpuErrchk(cudaMalloc(&_d_nn,         sizeof(int)*size));
  gpuErrchk(cudaMalloc(&_d_SR,         sizeof(G_Vector3)*size));
  gpuErrchk(cudaMalloc(&_d_F,          sizeof(G_Vector3)*size));
  gpuErrchk(cudaMalloc(&_d_F_padding,  sizeof(G_Vector3)*size*_NNM));
  gpuErrchk(cudaMalloc(&_d_VIRIAL_IND_element, sizeof(double)*9*size));
  gpuErrchk(cudaMalloc(&_d_VIRIAL_IND_element_padding, sizeof(double)*_NNM*9*size));
}

void EAMFrame::cuda_memcpy_all() {
  int size = _NP*allocmultiple;
  gpuErrchk(cudaMemcpy(_d_atpe3b,      atpe3b,      sizeof(double)*size,      cudaMemcpyHostToDevice) );
  gpuErrchk(cudaMemcpy(_d_rhotot,      rhotot,      sizeof(double)*size,      cudaMemcpyHostToDevice) );
  gpuErrchk(cudaMemcpy(_d_embf,        embf,        sizeof(double)*size,      cudaMemcpyHostToDevice) );
  gpuErrchk(cudaMemcpy(_d_embfp,       embfp,       sizeof(double)*size,      cudaMemcpyHostToDevice) );
  gpuErrchk(cudaMemcpy(_d_nbst,        nbst,        sizeof(int)*size,         cudaMemcpyHostToDevice) );

  /* copy eam data on grid to device */
  gpuErrchk(cudaMemcpy(_d_rho,         rho,         sizeof(double)*4*NGRID,   cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(_d_rhop,        rhop,        sizeof(double)*4*NGRID,   cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(_d_phi,         phi,         sizeof(double)*2*NGRID,   cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(_d_phip,        phip,        sizeof(double)*2*NGRID,   cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(_d_phix,        phix,        sizeof(double)*NGRID,     cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(_d_phipx,       phipx,       sizeof(double)*NGRID,     cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(_d_frho,        frho,        sizeof(double)*2*NGRID,   cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(_d_frhop,       frhop,       sizeof(double)*2*NGRID,   cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(_d_rho_spline,  rho_spline,  sizeof(double)*4*NGRID*4, cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(_d_phi_spline,  phi_spline,  sizeof(double)*2*NGRID*4, cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(_d_phix_spline, phix_spline, sizeof(double)*4*NGRID,   cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(_d_frho_spline, frho_spline, sizeof(double)*2*NGRID*4, cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(_d_rval,        rval,        sizeof(double)*NGRID,     cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(_d_rhoval,      rhoval,      sizeof(double)*NGRID,     cudaMemcpyHostToDevice));

  /* initial copy of runtime cpu data to device */
  gpuErrchk(cudaMemcpy(_d_nindex,      nindex,     sizeof(int)*size*_NNM,    cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(_d_nn,          nn,          sizeof(int)*size,         cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(_d_SR,          _SR,         sizeof(G_Vector3)*size,   cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(_d_F,           _F,          sizeof(G_Vector3)*size,   cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(_d_VIRIAL_IND_element, _VIRIAL_IND[0].element,sizeof(double)*9*size, cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(_d_fscalars,    fscalars,    sizeof(double)*10,        cudaMemcpyHostToDevice));
}

__device__ double spline(double* _d_spline_coeff,int ind, double qq)
{
  double a, b, c, d, qq2, qq3, f;
  a = _d_spline_coeff[ind*4+0];
  b = _d_spline_coeff[ind*4+1];
  c = _d_spline_coeff[ind*4+2];
  d = _d_spline_coeff[ind*4+3];
  qq2=qq*qq; qq3=qq2*qq;
  f=a+b*qq+c*qq2+d*qq3;
  return f;
}  
__device__ double spline1(double* _d_spline_coeff,int ind, double qq)
{
  double b, c, d, qq2, fp;//, qq3;
  //double a = _d_spline_coeff[ind*4+0];
  b = _d_spline_coeff[ind*4+1];
  c = _d_spline_coeff[ind*4+2];
  d = _d_spline_coeff[ind*4+3];
  qq2=qq*qq; //qq3=qq2*qq;
  fp=b+2*c*qq+3*d*qq2;
  return fp;
}


/* 
   Given ipt, return a list of jpt, each of jpt`s neighbor list nindex[jpt] contain ipt 
   Every time neighborhood is refreshed, this matrix must be redone.
*/
//void EAMFrame::create_inverse_nn() { 
//	for(int i = 0; i< _NP; )
//
//}

__global__ void kernel_rhoeam(int _NP, int _NNM, int eamfiletype, int eamgrid,
                              double *_d_rho,double *_d_rhop,double *_d_phi,double *_d_phip,
                              double *_d_phix,double *_d_phipx,
                              double *_d_frho,double *_d_frhop,
                              double *_d_rho_spline, double *_d_phi_spline,
                              double *_d_phix_spline,double *_d_frho_spline,
                              double *_d_rval,double *_d_rhoval,
                              double *_d_atpe3b,double *_d_rhotot, double *_d_embf, double *_d_embfp,
			      int    *_d_nbst,
                              double *_d_EPOT,
                              double *_d_H_element,
                              double *_d_VIRIAL_element,
                              double *_d_EPOT_IND,
                              int *_d_species,
                              int *_d_nindex,
                              int *_d_nn,
                              G_Vector3 *_d_SR,
                              G_Vector3 *_d_F,
                              G_Vector3 *_d_F_padding,
                              double * _d_VIRIAL_IND_element,
                              double * _d_VIRIAL_IND_element_padding,
                              double *_d_fscalars)
{
    int i, j, l, jpt, idx, jdx, ind;
    G_Vector3 sij,rij;
    double r2ij, rmagg, qq, qr;
    double rhoi, rhoj;
    G_Matrix33 _d_VIRIAL(_d_VIRIAL_element);
    G_Matrix33 _d_H(_d_H_element);

    //double _d_rmass   = _d_fscalars[0];
    //double _d_rlatt   = _d_fscalars[1];
    double _d_drar    = _d_fscalars[2];
    double _d_drhoar  = _d_fscalars[3];
    //double _d_actual  = _d_fscalars[4];
    double _d_actual2 = _d_fscalars[5];
    double _d_rmin    = _d_fscalars[6];
    double _d_petrip  = _d_fscalars[7];
    double _d_rhocon  = _d_fscalars[8];
    double _d_rhomax  = _d_fscalars[9];

    for(i=blockDim.x*blockIdx.x+threadIdx.x;i<_NP;i+=blockDim.x*gridDim.x)
    {
        /*INFO("rhoeam");*/
        _d_petrip = 0.0;
        //rhocon = 1e-10;
        _d_rhocon = 0.0;            /* replaced by Keonwook Kang, Apr 29, 2011 */
        _d_rhomax = eamgrid*_d_drhoar; 

        _d_atpe3b[i]=0;
        _d_rhotot[i]=0;
        _d_nbst[i]=0;

        _d_F[i].clear(); _d_EPOT_IND[i]=0;
	for(l = 0; l<9; l++) {
	  _d_VIRIAL_IND_element[i*9+l] = 0; 
	  _d_VIRIAL_IND_element_padding[i*9+l] = 0;
	 }
        _d_EPOT=0;
        _d_VIRIAL.clear();
    
        /* modify here for binary systems (0/1) */
        idx = _d_species[i]; /* type of atom (i) */
        /* do on j-particles */
        for(j=0;j<_d_nn[i];j++)
        {
            /* modify here for binary systems (0/1) */
            jpt=_d_nindex[i*_NNM+j]; /* nindex[i][j]; */
            jdx = _d_species[jpt]; /* type of atom (j) */
            if(i>=jpt) continue;
            sij=_d_SR[jpt]-_d_SR[i];
            sij.subint();
            rij=_d_H*sij;
            r2ij=rij.norm2();
            if(r2ij>_d_actual2) continue;
            rmagg=sqrt(r2ij)-_d_rmin;
            
            _d_nbst[i]++;
            ind = (int)(rmagg/_d_drar);

            if(ind>=eamgrid)
            {
                ind=eamgrid-1;
                printf("ind = %d r=%f in RHOEAM\n",ind, rmagg+_d_rmin);
            }
            else if(ind<0)
            {
                ind=0;
                printf("ind = %d r=%f in RHOEAM\n",ind, rmagg+_d_rmin);
            }
            qq=rmagg-_d_rval[ind];

            if(idx==jdx)
            {
#ifndef _CUBICSPLINE
            rhoi=_d_rho[jdx*NGRID+ind]+ qq*_d_rhop[jdx*NGRID+ind];
#else
            //rhoi = interp(rho[jdx],rhop[jdx],drar,ind,qq);
            rhoi = spline(_d_rho_spline[jdx],ind,qq);
#endif
            _d_rhotot[i]+=rhoi;
            _d_rhotot[jpt]+=rhoi;
            }
            else
            {
              if (eamfiletype == 2)
              {
#ifndef _CUBICSPLINE
              rhoi=_d_rho[jdx*NGRID+ind]+ qq*_d_rhop[jdx*NGRID+ind];
              rhoj=_d_rho[idx*NGRID+ind]+ qq*_d_rhop[idx*NGRID+ind];
#else
              //rhoi = interp(rho[jdx],rhop[jdx],drar,ind,qq);
              //rhoj = interp(rho[idx],rhop[idx],drar,ind,qq);
              rhoi = spline(_d_rho_spline[jdx*NGRID*4],ind,qq);
              rhoj = spline(_d_rho_spline[idx*NGRID*4],ind,qq);
#endif
              } else if (eamfiletype == 4)
              {
#ifndef _CUBICSPLINE
              rhoi=_d_rho[(idx+2)*NGRID+ind]+ qq*_d_rhop[(idx+2)*NGRID+ind];
              rhoj=_d_rho[(jdx+2)*NGRID+ind]+ qq*_d_rhop[(jdx+2)*NGRID+ind];
#else
              rhoi = spline(_d_rho_spline[(idx+2)*NGRID*4],ind,qq);
              rhoj = spline(_d_rho_spline[(jdx+2)*NGRID*4],ind,qq);
#endif
              }
              _d_rhotot[i]+=rhoi;
              _d_rhotot[jpt]+=rhoj;
            }
        }
    }

    
    for(i=blockDim.x*blockIdx.x+threadIdx.x;i<_NP;i+=blockDim.x*gridDim.x)
    {
        /* modify here for binary systems (0/1) */
        idx = _d_species[i]; /* type of atom (i) */
        if(_d_rhotot[i]<_d_rhocon)
        {
            _d_rhotot[i]=_d_rhocon;
        }
        if(_d_rhotot[i]>_d_rhomax)
        {
            _d_rhotot[i]=_d_rhomax;
        }
        ind = (int)(_d_rhotot[i]/_d_drhoar);
        if(ind>=eamgrid-1) ind=eamgrid-1;
        qr = _d_rhotot[i] - _d_rhoval[ind];

#ifndef _CUBICSPLINE
        _d_embf[i] = _d_frho[idx*NGRID+ind] + qr*_d_frhop[idx*NGRID+ind];
        _d_embfp[i] = _d_frhop[idx*NGRID+ind] +
          qr*(_d_frhop[idx*NGRID+ind+1]-_d_frhop[idx*NGRID+ind])/_d_drhoar;
#else
        //embf[i] = interp(frho[idx],frhop[idx],drhoar,ind,qr);
        //embfp[i] = interp1(frho[idx],frhop[idx],drhoar,ind,qr);
        _d_embf[i] = spline(_d_frho_spline[idx*4*NGRID],ind,qr);
        _d_embfp[i] = spline1(_d_frho_spline[idx*4*NGRID],ind,qr);
#endif
        _d_atpe3b[i] = _d_embf[i];
        _d_EPOT_IND[i]+=_d_atpe3b[i];
        _d_EPOT[i]+=_d_atpe3b[i];
    }    

    /* debug */
#if 0
    for(i=0;i<_NP;i++)
    {
        INFO_Printf("atom[%d] embf=%e\n",i,embf[i]);
    }
#endif
}

void EAMFrame::rhoeam_cuda() {
  gpuErrchk(cudaMemcpy(_d_SR,_SR,_NP*sizeof(G_Vector3), cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(_d_H_element,_H.element,9*sizeof(double),cudaMemcpyHostToDevice));
  kernel_rhoeam<<< 1,1 >>>(_NP, _NNM, eamfiletype, eamgrid,
                           _d_rho, _d_rhop, _d_phi, _d_phip,
                           _d_phix,_d_phipx,
                           _d_frho,_d_frhop,
                           _d_rho_spline, _d_phi_spline,
                           _d_phix_spline,_d_frho_spline,
                           _d_rval,_d_rhoval,
                           _d_atpe3b,_d_rhotot, _d_embf,_d_embfp,
			   _d_nbst,
                           _d_EPOT,
                           _d_H_element,
                           _d_VIRIAL_element,
                           _d_EPOT_IND,
                           _d_species,
                           _d_nindex,
                           _d_nn,
                           _d_SR,
                           _d_F,
                           _d_F_padding,
                           _d_VIRIAL_IND_element,
                           _d_VIRIAL_IND_element_padding,
                           _d_fscalars);
  thrust::device_ptr<double> t_EPOT = thrust::device_pointer_cast(_d_EPOT);
  _EPOT = thrust::reduce(t_EPOT,t_EPOT+_NP); 
}



__global__ void kernel_assemble_back_force(int _NP, int _NNM, int *_d_nn,int *_d_nindex,G_Vector3 *_d_F,G_Vector3 *_d_F_padding, double* _d_VIRIAL_IND_element, double* _d_VIRIAL_IND_element_padding) { 
  int i, j, jpt, l, m;
  for(i = blockDim.x * blockIdx.x + threadIdx.x; i<_NP;i+=blockDim.x*gridDim.x) {
    _d_F[i].clear();  
    for (j = 0; j<_d_nn[i]; j++) { 
      jpt = _d_nindex[i*_NNM+j];
      _d_F[i] += _d_F_padding[ jpt ];
      for(l = 0;l<3;l++) for(m= 0;m<3;m++)
        _d_VIRIAL_IND_element[i*9+l*3+m] += _d_VIRIAL_IND_element_padding[jpt*9+l*3+m];
    }
  }
}


__global__ void kernel_kraeam(int _NP, int _NNM, int eamfiletype, int eamgrid,
                              double *_d_rho, double *_d_rhop, double *_d_phi, double *_d_phip,
                              double *_d_phix, double *_d_phipx,
                              double *_d_frho, double *_d_frhop,
                              double *_d_rho_spline, double *_d_phi_spline,
                              double *_d_phix_spline, double *_d_frho_spline,
                              double *_d_rval, double *_d_rhoval,
                              double *_d_atpe3b,double *_d_rhotot, double *_d_embf, double *_d_embfp,
			      int    *_d_nbst,
                              double *_d_EPOT,
                              double *_d_H_element,
                              double *_d_VIRIAL_element,
                              double *_d_EPOT_IND,
                              int *_d_species,
                              int *_d_nindex,
                              int *_d_nn,
                              G_Vector3 *_d_SR,
                              G_Vector3 *_d_F,
                              G_Vector3 *_d_F_padding,
                              double *_d_VIRIAL_IND_element,
                              double *_d_VIRIAL_IND_element_padding,
                              double *_d_fscalars)
{
    int i, j, l, m, jpt, idx, jdx, ind;
    G_Vector3 sij,rij,fij;
    G_Matrix33 _d_H(_d_H_element);
    G_Matrix33 _d_VIRIAL(_d_VIRIAL_element);
    double r2ij, rmagg, pp, qq, fcp, fpp, fp, denspi, denspj;

    //double _d_rmass   = _d_fscalars[0];
    //double _d_rlatt   = _d_fscalars[1];
    double _d_drar    = _d_fscalars[2];
    //double _d_drhoar  = _d_fscalars[3];
    //double _d_actual  = _d_fscalars[4];
    double _d_actual2 = _d_fscalars[5];
    double _d_rmin    = _d_fscalars[6];
    //double _d_petrip  = _d_fscalars[7];
    //double _d_rhocon  = _d_fscalars[8];
    //double _d_rhomax  = _d_fscalars[9];
    
    /*	start calculation
     *  calculate pair contribution to total energy and forces
     */
    
    /* do on i-particles */
    for(i=blockDim.x*blockIdx.x+threadIdx.x;i<_NP;i+=blockDim.x*gridDim.x)
    {
	    G_Matrix33 _d_VIRIAL_IND_i(_d_VIRIAL_IND_element+i*9);
	    G_Matrix33 _d_VIRIAL_IND_jpt(_d_VIRIAL_IND_element+jpt*9);

        /* modify here for binary systems (0/1) */
        idx = _d_species[i]; /* type of atom (i) */
        /* do on j-particles */
        for(j=0;j<_d_nn[i];j++)
        {
            /* modify here for binary systems (0/1) */
            jpt=_d_nindex[i*_NNM+j];
            jdx = _d_species[jpt]; /* type of atom (j) */
            if(i>=jpt) continue;
            sij=_d_SR[jpt]-_d_SR[i];
            sij.subint();
            rij=_d_H*sij;
            r2ij=rij.norm2();
            if(r2ij>_d_actual2) continue;
            rmagg=sqrt(r2ij)-_d_rmin;

            ind = int(rmagg/_d_drar);
            
            if(ind>=eamgrid)
            {
                ind=eamgrid-1;
                printf("ind = %d in RHOEAM\n",ind);
            }
            else if(ind<0)
            {
                ind=0;
                printf("ind = %d in RHOEAM\n",ind);
            }
            qq=rmagg-_d_rval[ind];

            if(idx==jdx)
            {
#ifndef _CUBICSPLINE
              pp = _d_phi[idx*NGRID+ind] + qq*_d_phip[idx*NGRID+ind];
              fpp = _d_phip[idx*NGRID+ind] +
                qq*(_d_phip[idx*NGRID+ind+1]-_d_phip[idx*NGRID+ind])/_d_drar;
#else
              //pp = interp(phi[idx],phip[idx],drar,ind,qq);
              //fpp = interp1(phi[idx],phip[idx],drar,ind,qq);
              pp = spline(_d_phi_spline[idx*NGRID*4],ind,qq);
              fpp = spline1(_d_phi_spline[idx*NGRID*4],ind,qq);
#endif
            }
            else
            {
#ifndef _CUBICSPLINE
              pp =_d_phix[ind] + qq*_d_phipx[ind];
              fpp = _d_phipx[ind] + qq*(_d_phipx[ind+1]-_d_phipx[ind])/_d_drar;
#else
              //pp = interp(phix,phipx,drar,ind,qq);
              //fpp = interp1(phix,phipx,drar,ind,qq);
              pp = spline(_d_phix_spline,ind,qq);
              fpp = spline1(_d_phix_spline,ind,qq);
#endif
            }
            //INFO_Printf("phi = %20.18e\n",pp);
            if ( (idx==jdx) || (eamfiletype==2) )
            {
#ifndef _CUBICSPLINE
            denspi = _d_rhop[idx*NGRID+ind] +
                qq*(_d_rhop[idx*NGRID+ind+1]-_d_rhop[idx*NGRID+ind])/_d_drar ;
            denspj = _d_rhop[jdx*NGRID+ind] +
                qq*(_d_rhop[jdx*NGRID+ind+1]-_d_rhop[jdx*NGRID+ind])/_d_drar ; /* typo idx fixed to jdx ! */
#else
            //denspi = interp1(rho[idx],rhop[idx],drar,ind,qq);
            //denspj = interp1(rho[jdx],rhop[jdx],drar,ind,qq);
            denspi = spline1(_d_rho_spline[idx],ind,qq);
            denspj = spline1(_d_rho_spline[jdx],ind,qq);
#endif
            } else if ( (idx!=jdx) && (eamfiletype==4) ) {
#ifndef _CUBICSPLINE
            denspi = _d_rhop[(jdx+2)*NGRID+ind] +
                qq*(_d_rhop[(jdx+2)*NGRID+ind+1]-_d_rhop[(jdx+2)*NGRID+ind])/_d_drar ;
            denspj = _d_rhop[(idx+2)*NGRID+ind] +
                qq*(_d_rhop[(idx+2)*NGRID+ind+1]-_d_rhop[(idx+2)*NGRID+ind])/_d_drar ;
#else
            denspi = spline1(_d_rho_spline[(jdx+2)*NGRID*4],ind,qq);
            denspj = spline1(_d_rho_spline[(idx+2)*NGRID*4],ind,qq);
#endif
            }

            fcp = denspj * _d_embfp[i] + denspi * _d_embfp[jpt];
            fp = (fpp + fcp) / (rmagg+_d_rmin);
            _d_EPOT_IND[i]+= 0.5*pp;
            _d_EPOT_IND[jpt]+= 0.5*pp;
            _d_EPOT[i]+=pp;
            
            fij=rij*fp;
            _d_F[i]+=fij;
            _d_F_padding[jpt]-=fij;
#if 0
            if(!(fixed[i]||fixed[jpt]))
                _VIRIAL.addnvv(-fp,rij,rij);
            else if(!(fixed[i]&&fixed[jpt]))
                _VIRIAL.addnvv(-0.5*fp,rij,rij);
#else
	    _d_VIRIAL_IND_i.addnvv(-.5*fp,rij,rij);
	    _d_VIRIAL_IND_jpt.addnvv(-.5*fp,rij,rij);
            //_VIRIAL.addnvv(-fp,rij,rij);
            //assert(SURFTEN==0);
#if 0
            if (SURFTEN==1 && (curstep%savepropfreq)==1) AddnvvtoPtPn(_SR[jpt],rij,rij,-fp);
#endif
#endif
        }
	for(l = 0; l<3; l++)  for(m = 0; m<3; m++) { 
	  _d_VIRIAL_IND_element[i*9+l*3+m] = _d_VIRIAL_IND_i[l][m]; 
	  _d_VIRIAL_IND_element_padding[jpt*9+l*3+m]=_d_VIRIAL_IND_jpt[l][m]; 
	}
        _d_VIRIAL+=_d_VIRIAL_IND_i;
    }
    for(l = 0; l<3; l++) for (m = 0;m<3; m++)
      _d_VIRIAL_element[l*3+m] = _d_VIRIAL[l][m];
}

void EAMFrame::kraeam_cuda() {
  kernel_kraeam<<< 1,1 >>>(_NP,  _NNM, eamfiletype,  eamgrid,
                           _d_rho, _d_rhop, _d_phi, _d_phip,
                           _d_phix, _d_phipx,
                           _d_frho, _d_frhop,
                           _d_rho_spline, _d_phi_spline,
                           _d_phix_spline, _d_frho_spline,
                           _d_rval, _d_rhoval,
                           _d_atpe3b, _d_rhotot, _d_embf, _d_embfp,
			   _d_nbst,
                           _d_EPOT,
                           _d_H_element,
                           _d_VIRIAL_element,
                           _d_EPOT_IND,
                           _d_species,
                           _d_nindex,
                           _d_nn,
                           _d_SR,
                           _d_F,
                           _d_F_padding,
                           _d_VIRIAL_IND_element,
                           _d_VIRIAL_IND_element_padding,
                           _d_fscalars);
  kernel_assemble_back_force<<<1,1>>>(_NP, _NNM, _d_nn, _d_nindex, _d_F, _d_F_padding, _d_VIRIAL_IND_element, _d_VIRIAL_IND_element_padding);
  /* Copy force (Vector3 *) back to CPU for relax function to call */ 
  cudaMemcpy(_F, _d_F, _NP*sizeof(G_Vector3), cudaMemcpyDeviceToHost);
}


void EAMFrame::free_device_ptr() {
  cudaFree(_d_rho);
  cudaFree(_d_rhop);
  cudaFree(_d_phi);
  cudaFree(_d_phip);

  cudaFree(_d_phix);
  cudaFree(_d_phipx);
  cudaFree(_d_frho);
  cudaFree(_d_frhop);
  cudaFree(_d_rho_spline);
  cudaFree(_d_phi_spline);
  cudaFree(_d_phix_spline);
  cudaFree(_d_frho_spline);
  cudaFree(_d_rval);
  cudaFree(_d_rhoval);

  cudaFree(_d_atpe3b);
  cudaFree(_d_rhotot);
  cudaFree(_d_embf);
  cudaFree(_d_embfp);
  cudaFree(_d_H_element);
  cudaFree(_d_EPOT_IND);
  cudaFree(_d_species);        
  cudaFree(_d_nindex);
  cudaFree(_d_nn);
  cudaFree(_d_SR);
  cudaFree(_d_F);
  cudaFree(_d_F_padding);
  cudaFree(_d_VIRIAL_IND_element);
  cudaFree(_d_VIRIAL_element);
}

/* This is a simple test for GPU. run the function to see if maxErro == 0. If not, GPU device is not set correctly */
__global__ void saxpy(int n, float a, const float *x, float *y)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  if (i < n) y[i] = a*x[i] + y[i];
}
int EAMFrame::test_saxpy(void)
{
  int N = 1<<20;
  float *x, *y, *d_x, *d_y;
  x = (float*)malloc(N*sizeof(float));
  y = (float*)malloc(N*sizeof(float));
  cudaMalloc(&d_x, N*sizeof(float));
  cudaMalloc(&d_y, N*sizeof(float));
  for (int i = 0; i < N; i++) { x[i] = 1.0f; y[i] = 2.0f; }
  cudaMemcpy(d_x, x, N*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_y, y, N*sizeof(float), cudaMemcpyHostToDevice);

  // Perform SAXPY on 1M elements
  saxpy<<<(N+255)/256, 256>>>(N, 2.0f, d_x, d_y);
  cudaMemcpy(y, d_y, N*sizeof(float), cudaMemcpyDeviceToHost);
  float maxError = 0.0f;
  for (int i = 0; i < N; i++) maxError = max(maxError, abs(y[i]-4.0f));
  INFO_Printf("Max error: %f\n", maxError);
  cudaFree(d_x);
  cudaFree(d_y);
  free(x);
  free(y);
  return 0;
}

