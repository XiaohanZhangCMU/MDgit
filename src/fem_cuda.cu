#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <thrust/extrema.h>
#include <thrust/inner_product.h>
#include <thrust/transform.h>
#include <thrust/functional.h>
#include <thrust/device_malloc.h>
#include <thrust/device_free.h>

#include "fem.h"

void FEMFrame::cuda_memory_alloc() {
  assert(allocmultiple == 1);
  assert(_NP > 0);
  cudaMalloc( &_d_map23,         2 * sizeof(int));
  cudaMalloc( &_d_SR,           _NP*sizeof(G_Vector3));
  cudaMalloc( &_d_SRref,        _NP*sizeof(G_Vector3));
  cudaMalloc( &_d_Rref,         _NP*sizeof(G_Vector3));
  cudaMalloc( &_d_F,            _NP*sizeof(G_Vector3));
  cudaMalloc( &_d_fixed,        _NP*sizeof(int));
  cudaMalloc( &_d_EPOT_IND,     _NP*sizeof(double));
  cudaMalloc( &_d_EPOT_RMV,     _NP*sizeof(double));

  cudaMemset(_d_F_padding, 0,   _NP*_MAX_NELEM_SHARE_NODE*sizeof(double));
  cudaMemset(_d_EPOT,      0,   _NP*sizeof(double));

  gpuErrchk(cudaMalloc(&_d_H_element, 3*3*sizeof(double))); 
}

void FEMFrame::cuda_memory_alloc_elements() {
  assert(_NELE > 0);
  cudaMalloc( &_d_EPOT,             _NELE*sizeof(double));
  cudaMalloc( &_d_elements,         _NELE*_NNODE_PER_ELEMENT*sizeof(int));
  cudaMalloc( &_d_inv_elements,     _NP*_MAX_NELEM_SHARE_NODE*sizeof(int));
}

void FEMFrame::cuda_memory_alloc_element_coeff() { 
  assert(_NINT_PER_ELEMENT > 0);
  assert(_NNODE_PER_ELEMENT > 0);
  assert(_NDIM > 0);
  assert(_NELE > 0);
  int size =  _NINT_PER_ELEMENT*_NDIM*_NDIM*_NDIM*_NNODE_PER_ELEMENT*_NELE;
  cudaMalloc( &_d_gauss_weight, _NINT_PER_ELEMENT*sizeof(double));
  cudaMalloc( &_d_dFdu,         size*sizeof(double));
}

void FEMFrame::create_inverse_connectivities_matrix() { 
  assert(_NP > 0);
  int i, j, k, jpt;
  Realloc(inv_elements,int,_MAX_NELEM_SHARE_NODE*_NP*sizeof(int));
  memset(inv_elements,-1,  _MAX_NELEM_SHARE_NODE*_NP*sizeof(int));

  INFO_Printf("NP = %d", _NP);
 // for(int i = 0;i<_MAX_NELEM_SHARE_NODE*_NP;i++) { 
 //   INFO_Printf("inv_elements[%d] = %d\n",i, inv_elements[i]) ;
 // }
   
  //for (i = 0; i< _NP; i++) for (k = 0;k<_MAX_NELEM_SHARE_NODE;k++) inv_elements[i*_MAX_NELEM_SHARE_NODE+k] = -1;

  for (i = 0; i< _NELE; i++) {
    for (j = 0; j< _NNODE_PER_ELEMENT; j++) {
      jpt = elements[_NNODE_PER_ELEMENT*i + j];
      k = 0;
      for (k = 0; k < _MAX_NELEM_SHARE_NODE; k++) {
          //INFO_Printf("here 1 : k = %d, i = %d, j=  %d, inv_ind=%d\n", k, i, j, inv_elements[jpt*_MAX_NELEM_SHARE_NODE+k]);
       
        if (inv_elements[jpt*_MAX_NELEM_SHARE_NODE+k] == -1) { 
          //INFO_Printf("here 2 : k = %d, i = %d, j=  %d\n", k, i, j);
          inv_elements[jpt*_MAX_NELEM_SHARE_NODE+k] = i*_NNODE_PER_ELEMENT +j;
          break;
        }
      }
      if (k == _MAX_NELEM_SHARE_NODE) FATAL("Too many elements share a node!");
    }
  } 
}

void FEMFrame::cuda_memcpy_all() {
  assert(sizeof(G_Vector3) == sizeof(Vector3));
  assert(sizeof(G_Matrix33) == sizeof(Matrix33));
  assert(_NELE > 0); assert(_NNODE_PER_ELEMENT>0);assert(_NINT_PER_ELEMENT>0);
  assert(_H[0][0]>0 && _H[1][1]>0 && _H[2][2]>0);
  int size1 = _NELE*_NNODE_PER_ELEMENT;
  int size2 = _NP*_MAX_NELEM_SHARE_NODE;

  cudaMemcpy( _d_map23,   map23,      2*sizeof(int),         cudaMemcpyHostToDevice);
  cudaMemcpy( _d_SR,      _SR,        _NP*sizeof(G_Vector3), cudaMemcpyHostToDevice);
  cudaMemcpy( _d_SRref,   _SRref,     _NP*sizeof(G_Vector3), cudaMemcpyHostToDevice);
  cudaMemcpy( _d_Rref,    _Rref,      _NP*sizeof(G_Vector3), cudaMemcpyHostToDevice);
  cudaMemcpy( _d_F,       _F,         _NP*sizeof(G_Vector3), cudaMemcpyHostToDevice);
  cudaMemcpy( _d_fixed,   fixed,      _NP*sizeof(int),       cudaMemcpyHostToDevice);
  cudaMemcpy( _d_EPOT_IND, _EPOT_IND, _NP*sizeof(double),    cudaMemcpyHostToDevice);
  cudaMemcpy( _d_EPOT_RMV, _EPOT_RMV, _NP*sizeof(double),    cudaMemcpyHostToDevice);

  cudaMemcpy( _d_H_element, _H.element, 3*3*sizeof(double),  cudaMemcpyHostToDevice); //???????????
  cudaMemcpy(_d_elements,    elements,    size1*sizeof(int),cudaMemcpyHostToDevice);
  cudaMemcpy(_d_inv_elements,inv_elements,size2*sizeof(int),cudaMemcpyHostToDevice);
  cudaMemcpy(_d_gauss_weight,gauss_weight,_NINT_PER_ELEMENT*sizeof(double),cudaMemcpyHostToDevice);

  create_inverse_connectivities_matrix();
}

__device__ G_Matrix33 getEigenF(G_Vector3 p, G_Matrix33 Fdef, 
                                double y_eigen_zbound_max, double y_eigen_zbound_min,
                                double x_eigen_zbound_max, double x_eigen_zbound_min,
                                double y_eigen_strain, double x_eigen_strain ) {
  G_Matrix33 I;    
  I[0][0]=I[1][1]=I[2][2]=1;
   if (p[2] <= y_eigen_zbound_max && p[2] >= y_eigen_zbound_min ){ 
    I[1][1] = y_eigen_strain;
  }
  if (p[2] <= x_eigen_zbound_max && p[2] >= x_eigen_zbound_min ){ 
    I[0][0] = x_eigen_strain;
  }
  return I;
}

__global__ void kernel_snap_fem_energy_force(int _NDIM, int _NELE, int _NINT_PER_ELEMENT, int _NNODE_PER_ELEMENT, int *_d_elements, int *_d_inv_elements, int *_d_map23, int *_d_fixed, double *_d_gauss_weight, double *_d_dFdu, double *_d_EPOT, G_Vector3 *_d_SR, G_Vector3 *_d_SRref, G_Vector3 *_d_Rref, G_Vector3 *_d_F, G_Vector3 *_d_F_padding, double* _d_H_element, 
    double y_eigen_zbound_max, double y_eigen_zbound_min,
    double x_eigen_zbound_max, double x_eigen_zbound_min,
    double y_eigen_strain, double x_eigen_strain ) {

  int iele = blockDim.x * blockIdx.x + threadIdx.x; 
  if (iele < _NELE) { 
//  printf("I am here 1 : iele = %d\n", iele);
    int i,iele,j,jpt,iA, in, ip, iq, ind, p, q, r;
    G_Vector3 dsj, drj, elem_center;
    G_Matrix33 Fe, Fdef, B, C, Cinv, FCinvtran, dEdF, hinv;
    G_Matrix33 eigF, invEigF;
    double Eele, Eint, Jet, trace, detB, J2;
    double lambda, mu, temp;
    G_Matrix33 E, E2, I, pk2, temp1, temp2;    
    mu = 1; lambda = 1.95; 

    G_Matrix33 _d_H(_d_H_element);
    I.eye(); //I[0][0]=I[1][1]=I[2][2]=1;
    /* energy of this element */
    hinv = _d_H.inv();
    Eele = 0;
  //  printf("_d_H= %g, %g, %g\n", _d_H[0][0], _d_H[1][1], _d_H[2][2]);
  //  printf("hinv= %g, %g, %g\n", hinv[0][0], hinv[1][1], hinv[2][2]);

    /* center of the element */
    elem_center.clear();elem_center[0] = elem_center[1]=elem_center[2]= 0;
    for(j=0;j<_NNODE_PER_ELEMENT;j++) {
	    jpt=_d_elements[iele*_NNODE_PER_ELEMENT+j];	
      for (i = 0;i<_NDIM;i++)  {
	      elem_center[i] += 1.0/_NNODE_PER_ELEMENT *_d_Rref[jpt][i];
	    }
    }

    for(iA=0;iA<_NINT_PER_ELEMENT;iA++) {
	    /* energy of this Gauss integration point */
	    Eint = 0;
	    /* deformation gradient */
	    Fdef.clear(); Fdef[0][0] = Fdef[1][1] = Fdef[2][2] = 1.0;
	    for(j=0;j<_NNODE_PER_ELEMENT;j++) {
	      jpt=_d_elements[iele*_NNODE_PER_ELEMENT+j];
//    printf("_d_Rref= %g, %g, %g\n", _d_Rref[jpt][0], _d_Rref[jpt][1], _d_Rref[jpt][2]);
	      _d_SRref[jpt ] = hinv*_d_Rref[jpt];
//    printf("_d_SRref= %g, %g, %g\n", _d_SRref[jpt][0], _d_SRref[jpt][1], _d_SRref[jpt][2]);
	      dsj = _d_SR[jpt] - _d_SRref[jpt];
	      dsj.subint();
	      _d_H.multiply(dsj, drj);
	  
	      /* Add contribution from drj to F using dFdu */
	      for(ip=0;ip<_NDIM;ip++) {
	        for(iq=0;iq<_NDIM;iq++) {
	          for(in=0;in<_NDIM;in++) {
	   	        ind=(((iA*_NDIM+ip)*_NDIM+iq)*_NDIM+in)*_NNODE_PER_ELEMENT+j;
//	                printf("_d_dFdu[%d] = %g",ind, _d_dFdu[ind]);
//	                printf("drj[%d] = %g", in, drj[in]);
		          if (_NDIM == 2)
           		  Fdef[ip][iq] += _d_dFdu[ind]*drj[_d_map23[in]]; 
              else //_NDIM == 3
           		  Fdef[ip][iq] += _d_dFdu[ind]*drj[in];
	      } } }
      }

	    eigF = getEigenF(elem_center, Fdef, y_eigen_zbound_max, y_eigen_zbound_min, x_eigen_zbound_max, x_eigen_zbound_min, y_eigen_strain,     x_eigen_strain );

//	    printf("eigF = %g, %g, %g, %g, %g, %g", eigF[0][0], eigF[1][1], eigF[2][2], eigF[0][1], eigF[0][2], eigF[1][2]);
	    invEigF = eigF.inv();
	    Fe = Fdef*invEigF;     
//	    printf("Fe = %g, %g, %g, %g, %g, %g", Fe[0][0], Fe[1][1], Fe[2][2], Fe[0][1], Fe[0][2], Fe[1][2]);
	    E = Fe.tran()*Fe-I;
	    B = Fe*Fe.tran();
	    C = Fe.tran()*Fe;
	    Cinv = C.inv();
	    Jet = Fdef.det();  J2 = Jet*Jet;
	    detB = B.det(); 
	    Eint = 0;

    	/* Add contribution from F to Eint */ 
    	if (_NDIM == 2) {
    	  trace = B[0][0] + B[1][1];
    	  Eint = 0.5*(trace + 1.0/detB - 3); /* multiply MU and V0 later */
    	  dEdF = FCinvtran*(-1.0/J2) + Fdef; /* multiply MU later */
    	  /* Add contribution from drj to F using dFdu */
    	  for(j=0;j<_NNODE_PER_ELEMENT;j++) {
    	    //jpt=elements[iele*_NNODE_PER_ELEMENT+j];
    	    jpt=iele*_NNODE_PER_ELEMENT+j;
    	    for(ip=0;ip<_NDIM;ip++) {
    	      for(iq=0;iq<_NDIM;iq++) {
    		      for(in=0;in<_NDIM;in++) {
    		        ind=(((iA*_NDIM+ip)*_NDIM+iq)*_NDIM+in)*_NNODE_PER_ELEMENT+j;
    		        if (_d_fixed[jpt] == 0)				  
    		          _d_F_padding[jpt][_d_map23[in] ] -= dEdF[ip][iq]*_d_dFdu[ind]; 
    		} } } }
    	}
      
    	/* Add contribution from F to Eint */ 
	    if (_NDIM == 3) {
		
#if defined NeoHooken
        Eint = 1.0/8.0 * (0.5*lambda*log(Jet)*log(Jet) - mu*log(Jet) + 0.5*mu*(C.trace()-3));
	      temp1 = I;
	      temp1*= lambda * log(Jet);
	      temp2 = I;
	      temp2 = temp2-Cinv;
	      temp2*= mu;
	      pk2  = temp1 + temp2;
	      dEdF = Fdef * pk2;
#else
        E2 = E*E;
 	      dEdF.clear();
	      for (i = 0;i<_NDIM;i++)
	        for (j = 0;j<_NDIM;j++)
	          for (p = 0;p<_NDIM;p++)
		          for (r = 0;r<_NDIM;r++) {
		            temp = 0;
		            for (q = 0;q<_NDIM;q++)
		              temp += 2*mu*invEigF[j][p]*Fdef[i][r]*invEigF[r][q]*E[p][q] + 
		                2*mu*invEigF[r][p]*Fdef[i][r]*invEigF[j][q]*E[p][q];
		            dEdF[i][j] += 0.5*(2*lambda*E.trace()*invEigF[j][p]*Fdef[i][r] * invEigF[r][p] + temp);
		          }
        Eint = 0.5*lambda*(E.trace())*(E.trace()) + mu*(E2.trace());
#endif
	      /* Add contribution from drj to F using dFdu */

    	  for(j=0;j<_NNODE_PER_ELEMENT;j++) {
    	    //jpt=_d_elements[iele*_NNODE_PER_ELEMENT+j];
    	    jpt=iele*_NNODE_PER_ELEMENT+j;
    	    for(ip=0;ip<_NDIM;ip++) {
    	      for(iq=0;iq<_NDIM;iq++) {
    		      for(in=0;in<_NDIM;in++) {
    		        ind=(((iA*_NDIM+ip)*_NDIM+iq)*_NDIM+in)*_NNODE_PER_ELEMENT+j;
    		        if (_d_fixed[jpt] == 0)				  
    		          _d_F_padding[jpt][in ] -= dEdF[ip][iq]*_d_dFdu[ind]; 				
    		} } } }
      }
      /* Add contribution from Eint to Eele */
	    printf("gs_weight[%d]=%g, Eint = %g\n", iA, _d_gauss_weight[iA], Eint);
	    Eele += _d_gauss_weight[iA]*Eint;
    }
    printf("Eele = %g, iele = %d\n", Eele, iele);
    _d_EPOT[iele] = Eele;
  }
}

__global__ void kernel_assemble_back_force(int _NP, int *_d_inv_elements, G_Vector3 *_d_F, G_Vector3 *_d_F_padding,  double *_d_EPOT_IND, double *_d_EPOT_RMV, G_Matrix33* _d_VIRIAL_IND) { 
  int ind = blockDim.x * blockIdx.x + threadIdx.x; 
  if (ind < _NP) { 
    _d_F[ind].clear();  _d_EPOT_IND[ind]=0; 
    _d_EPOT_RMV[ind]=0; _d_VIRIAL_IND[ind].clear();

    for (int k = 0; k<_MAX_NELEM_SHARE_NODE; k++) { 
      int indice = ind *_MAX_NELEM_SHARE_NODE+k;
      if (_d_inv_elements[ indice ] > 0)
        _d_F[ind] += _d_F_padding[ _d_inv_elements[indice] ];
    }
  }
}

void FEMFrame::cuda_snap_fem_energy_force() {
  DUMP("FEM");

  _EPOT=0; //put this back

  _VIRIAL.clear();
  INFO_Printf("I am here, _EPOT = %g", _EPOT);
  kernel_snap_fem_energy_force<<<(_NELE/255),256>>>(_NDIM, _NELE, _NINT_PER_ELEMENT, _NNODE_PER_ELEMENT, _d_elements, _d_inv_elements, _d_map23, _d_fixed, _d_gauss_weight, _d_dFdu, _d_EPOT, _d_SR, _d_SRref, _d_Rref, _d_F, _d_F_padding, _d_H_element,
 y_eigen_zbound_max, y_eigen_zbound_min, 
 x_eigen_zbound_max, x_eigen_zbound_min,
 y_eigen_strain,     x_eigen_strain );

  kernel_assemble_back_force<<<(_NP/255),256>>>(_NP, _d_inv_elements, _d_F, _d_F_padding, _d_EPOT_IND, _d_EPOT_RMV, _d_VIRIAL_IND);

  Realloc( _h_EPOT, double, _NELE);
  cudaMemcpy( _h_EPOT,  _d_EPOT, _NELE*sizeof(double), cudaMemcpyDeviceToHost);
  for (int i = 0; i<_NELE; i ++) _EPOT += _h_EPOT[i];

  //thrust::device_ptr<double> t_EPOT = thrust::device_pointer_cast(_d_EPOT);
  //INFO_Printf("I am here 1, _EPOT = %g, _NELE = %d", _EPOT, _NELE);

  //_EPOT = thrust::reduce(t_EPOT,t_EPOT+_NELE);
  INFO_Printf("I am here 2, _EPOT = %g", _EPOT);
  
  cudaMemcpy( _F,       _d_F,         _NP*sizeof(G_Vector3), cudaMemcpyDeviceToHost);

}

void FEMFrame::free_device_ptr() {
  cudaFree(_d_elements);
  cudaFree(_d_inv_elements);
  cudaFree(_d_map23);
  cudaFree(_d_fixed);
  cudaFree(_d_colorids);
  cudaFree(_d_VIRIAL_IND);
  cudaFree(_d_gauss_weight); 
  cudaFree( _d_dFdu);
  cudaFree( _d_EPOT);
  cudaFree( _d_EPOT_IND);
  cudaFree(  _d_EPOT_RMV);
  cudaFree( _d_SR);
  cudaFree( _d_SRref);
  cudaFree( _d_Rref);        
  cudaFree( _d_F);
  cudaFree( _d_F_padding); //for padding forces of elements
}

__global__ void saxpy(int n, float a, const float *x, float *y)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  if (i < n) y[i] = a*x[i] + y[i];
  //printf("In GPU: n= %d: [%d,%d,%d]\n", n, blockIdx.x, blockDim.x, threadIdx.x);
}

int FEMFrame::test_saxpy(void)
{
  int N = 1<<20;
  float *x, *y, *d_x, *d_y;
  x = (float*)malloc(N*sizeof(float));
  y = (float*)malloc(N*sizeof(float));

  cudaMalloc(&d_x, N*sizeof(float));
  cudaMalloc(&d_y, N*sizeof(float));

  for (int i = 0; i < N; i++) {
    x[i] = 1.0f;
    y[i] = 2.0f;
  }
//  std::cout<<"This is a test"<<std::endl;
  cudaMemcpy(d_x, x, N*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_y, y, N*sizeof(float), cudaMemcpyHostToDevice);

  // Perform SAXPY on 1M elements
  saxpy<<<(N+255)/256, 256>>>(N, 2.0f, d_x, d_y);

  cudaMemcpy(y, d_y, N*sizeof(float), cudaMemcpyDeviceToHost);

  float maxError = 0.0f;
  for (int i = 0; i < N; i++)
    maxError = max(maxError, abs(y[i]-4.0f));
  printf("Max error: %f\n", maxError);

  cudaFree(d_x);
  cudaFree(d_y);
  free(x);
  free(y);

  return 0;
}



























