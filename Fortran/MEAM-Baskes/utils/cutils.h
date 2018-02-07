#include <stddef.h>
#include <stdlib.h>
#define FREE_ARG char*

int *ivector(int size);
int **imatrix(int size1, int size2);
int ***itensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
double *dvector(int size);
double **dmatrix(int size1, int size2);
double ***dtensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
float *fvector(int size);
float **fmatrix(int size1, int size2);
float ***ftensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_ivector(int *v);
void free_fvector(float *v);
void free_dvector(double *v);
void free_imatrix(int **m);
void free_fmatrix(float **m);
void free_dmatrix(double **m);
void free_itensor(int ***m);
void free_dtensor(double ***m);
void free_ftensor(float ***m);
double Integrate (double *f, int nf);
 

