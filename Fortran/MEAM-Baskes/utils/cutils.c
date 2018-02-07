#include "cutils.h"
#include <time.h>
#include <stdio.h>

/* trying to enable ctime in fortran */
char* ctime_(const time_t *timep) {
  return ctime(timep);
}

double Integrate (double *f, int nf){
    double s;
    int i;
    /* printf("We are integrating.\n");
       for (i=0;i<nf;i++)printf("%lf\t%d\n",f[i],i);*/
    s = 0.5 * (f[0] + f[nf-1]);
    for (i = 1; i < nf - 1; i ++){/*printf("%lf\n",s);*/ s = s + f[i];}
    return (s);
}

int *ivector(int size){
   int *v;
   v= (int *) malloc(size*sizeof(int));
   return v;
}

int **imatrix(int size1, int size2){
   int **m;
   int k;
   m=(int **) malloc (size1*sizeof(int *));
   m[0]=(int *) malloc (size1*size2*sizeof (int));
   for (k=1;k<size1; k++) m[k]=m[k-1]+size2;
   return m;
}

int ***itensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)

{	long i, j, nrow=nrh-nrl+1, ncol=nch-ncl+1,ndep=ndh-ndl+1;
	int ***t;

	t=(int ***)malloc((size_t)((nrow+1)*sizeof(int**)));
	t +=1;
	t -=nrl;

        t[nrl]=(int **)malloc((size_t)((nrow*ncol+1)*sizeof(int*)));
	t[nrl] += 1;
	t[nrl] -= ncl;

	t[nrl][ncl]=(int *)malloc((size_t)((nrow*ncol*ndep+1)*sizeof(int)));
	t[nrl][ncl] += 1;
	t[nrl][ncl] -= ndl;

	for(j=ncl;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++)
	{	t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
        }

	return t;
}

double ***dtensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)

{	long i, j, nrow=nrh-nrl+1, ncol=nch-ncl+1,ndep=ndh-ndl+1;
	double ***t;

	t=(double ***)calloc((size_t)(nrow+1),sizeof(double**));
if (!t) printf("allocation failure 1 in f3tensor()\n");
	t +=1;
	t -=nrl;

        t[nrl]=(double **)calloc((size_t)(nrow*ncol+1),sizeof(double*));
       	if (!t[nrl])printf("allocation failure 2 in f3tensor()\n");
        t[nrl] += 1;
	t[nrl] -= ncl;

	t[nrl][ncl]=(double *)calloc((size_t)(nrow*ncol*ndep+1),sizeof(double));
	if (!t[nrl][ncl]) printf("allocation failure 3 in f3tensor()\n");
	t[nrl][ncl] += 1;
	t[nrl][ncl] -= ndl;

	for(j=ncl;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++)
	{	t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
        }

	return t;
}

float ***ftensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)

{	long i, j, nrow=nrh-nrl+1, ncol=nch-ncl+1,ndep=ndh-ndl+1;
	float ***t;

	t=(float ***)malloc((size_t)((nrow+1)*sizeof(float**)));
	t +=1;
	t -=nrl;

        t[nrl]=(float **)malloc((size_t)((nrow*ncol+1)*sizeof(float*)));
	t[nrl] += 1;
	t[nrl] -= ncl;

	t[nrl][ncl]=(float *)malloc((size_t)((nrow*ncol*ndep+1)*sizeof(float)));
	t[nrl][ncl] += 1;
	t[nrl][ncl] -= ndl;

	for(j=ncl;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++)
	{	t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
        }

	return t;
}



double *dvector(int size){
   double *v;
   v= (double *) malloc(size*sizeof(double));
   return v;
}

double **dmatrix(int size1, int size2){
   double **m;
   int k;
   m=(double **) malloc (size1*sizeof(double *));
   m[0]=(double *) malloc (size1*size2*sizeof (double));
   for (k=1;k<size1; k++) m[k]=m[k-1]+size2;
   return m;
}


float *fvector(int size){
   float *v;
   v= (float *) malloc(size*sizeof(float));
   return v;
}

float **fmatrix(int size1, int size2){
   float **m;
   int k;
   m=(float **) malloc (size1*sizeof(float *));
   m[0]=(float *) malloc (size1*size2*sizeof (float));
   for (k=1;k<size1; k++) m[k]=m[k-1]+size2;
   return m;
}



void free_ivector(int *v){
   free((FREE_ARG)(v));
}

void free_fvector(float *v){
   free((FREE_ARG)(v));
}

void free_dvector(double *v){
    free((FREE_ARG)(v));
}

void free_imatrix(int **m){
    free((FREE_ARG)(m[0]));
    free((FREE_ARG)(m));
}

void free_fmatrix(float **m){
    free((FREE_ARG)(m[0]));
    free((FREE_ARG)(m));
}

void free_dmatrix(double **m){
    free((FREE_ARG)(m[0]));
    free((FREE_ARG)(m));
}

void free_itensor(int ***m){
    free((FREE_ARG)(m[0][0]));
    free((FREE_ARG)(m[0]));
    free((FREE_ARG)(m));
}

void free_dtensor(double ***m){
    free((FREE_ARG)(m[0][0]));
    free((FREE_ARG)(m[0]));
    free((FREE_ARG)(m));
}

void free_ftensor(float ***m){
    free((FREE_ARG)(m[0][0]));
    free((FREE_ARG)(m[0]));
    free((FREE_ARG)(m));
}

