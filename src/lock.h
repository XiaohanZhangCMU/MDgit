#ifndef __LOCK__H
#define __LOCK_H__

#include <cuda_runtime.h>
struct Lock {
  int *mutex;
  Lock() { 
    int state = 0;
    cudaMalloc ((void**) &mutex, sizeof(int));
    cudaMemcpy (mutex, &state, sizeof(int), cudaMemcpyHostToDevice); 
  }
  ~Lock() {
    cudaFree (mutex); 
  }
  __device__ void lock(){
    while(atomicCAS(mutex, 0, 1)!=0); 
  }
  __device__ void unlock(){
   atomicExch(mutex, 0); 
  }
};

#endif
