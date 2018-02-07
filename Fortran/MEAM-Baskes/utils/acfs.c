#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stddef.h>
#include "cutils.h"

double Sqr (double number);
double CalcTemp (double rv[],int itype[],double amass[], int natoms);

void evalvacf_ (int *nAtom, double rv[], double ra[], double rfAtom[],
         double *deltaT, double aMass[], int iType[],int *nBuffAcf, 
         int *nValAcf, int *nType, int totType[], int *limitAcfAv, 
         double *dens, double e[]) {
    double thermVec[3], v[3], viscVec[3 + 1];
    double fac;
    static double **thermAcfOrg, **thermAcf, *thermAcfAv;
    static double **rvAcfOrg, **rvAcf, **rvAcfAv,*tmpAcfAv;
    static double **viscAcfOrg, **viscAcf, *viscAcfAv;
    static int *indexAcf;
    static double *enAtom;
    static int nflag=1, cnave,nave=0;
    static int countAcfAv=0;
    static double CTemp;
    static double *AcfInt;
    int i, j, nb, ni;
    int k, n;       
    double tVal, temp;
    static char *fle="AcfData.txt";
    FILE *fp;
    
    /* printf("I am within evalvacf. %d\n",nave);*/
    nave++;
    if(nflag){
      /* Allocating Memory for static variables */
       indexAcf=ivector(*nBuffAcf);
       thermAcfAv=dvector(*nValAcf); 
       thermAcfOrg=dmatrix(*nBuffAcf,3);
       thermAcf=dmatrix(*nBuffAcf,*nValAcf);
       tmpAcfAv=dvector(*nValAcf);
       rvAcfAv=dmatrix(*nType,*nValAcf);
       rvAcf=dmatrix(*nType*(*nBuffAcf),*nValAcf);
       rvAcfOrg=dmatrix(*nType*(*nBuffAcf),3 *(*nAtom));
       viscAcfAv=dvector(*nValAcf);
       viscAcf=dmatrix(*nBuffAcf,*nValAcf);   
       viscAcfOrg=dmatrix(*nBuffAcf,3);
       AcfInt=dvector(*nType+2);
       enAtom=dvector(*nAtom);
       for(nb=0;nb<*nBuffAcf;nb++){
	 indexAcf[nb]=-nb*(*nValAcf)/(*nBuffAcf)-1;
         }
       for(j=0;j<*nValAcf;j++){
         for(k=0;k<*nType;k++){
	   rvAcfAv[k][j]=0.;}
         viscAcfAv[j]=thermAcfAv[j]=0.;
       }
       nflag=0;
       CTemp=0;
       cnave=nave;
       }
       CTemp= CTemp+CalcTemp(rv,iType,aMass,*nAtom);
    for (k = 0; k < 3; k ++) viscVec[k] = thermVec[k] = 0.;
    for (n = 0; n < *nAtom; n ++) {
      enAtom[n]=e[n];
      for (k = 0; k < 3; k ++){
         v[k] = rv[k+3+n*6] - 0.5 * ra[k+n*3] * (*deltaT);
                                enAtom[n]=enAtom[n]+ 0.5 * aMass[iType[n]-1] *Sqr(v[k]);
      }
#if MEAM
      viscVec[0] = viscVec[0] + aMass[iType[n]-1] * v[1] * v[2] + 0.5 * rfAtom[n*9+5];
      viscVec[1] = viscVec[1] + aMass[iType[n]-1] * v[2] * v[0] + 0.5 * rfAtom[n*9+2];
      viscVec[2] = viscVec[2] + aMass[iType[n]-1] * v[0] * v[1] + 0.5 * rfAtom[n*9+1];
      viscVec[0] = viscVec[0] + 0.5 * rfAtom[n*9+7];
      viscVec[1] = viscVec[1] + 0.5 * rfAtom[n*9+6];
      viscVec[2] = viscVec[2] + 0.5 * rfAtom[n*9+3];
#else
      viscVec[0] = viscVec[0] + aMass[iType[n]-1] * v[1] * v[2] + rfAtom[n*9+5];
      viscVec[1] = viscVec[1] + aMass[iType[n]-1] * v[2] * v[0] + rfAtom[n*9+2];
      viscVec[2] = viscVec[2] + aMass[iType[n]-1] * v[0] * v[1] + rfAtom[n*9+1];
#endif     
      for (k = 0; k < 3; k ++) {
        thermVec[k] = thermVec[k] + v[k] * enAtom[n];
        for (j = 0; j < 3; j ++)
           thermVec[k] = thermVec[k] +  v[j] *(
              rfAtom[k+j*3+9*n]);
    } }
    for (nb = 0; nb < *nBuffAcf; nb ++) {
      indexAcf[nb] = indexAcf[nb] + 1;
      if (indexAcf[nb] < 0){ continue;}
      if (indexAcf[nb] == 0) {
        for (k = 0; k < 3; k ++) {
          viscAcfOrg[nb][k] = viscVec[k];
          thermAcfOrg[nb][k] = thermVec[k];
        }/*printf("Initializing rvAcfOrg.\n");*/
        for (n = 0; n < *nAtom; n ++) {
	  /*	  printf("%d\t%d\t%lf\t%lf\t%lf\t%lf\n",n,iType[n]-1,rv[5+n*6],ra[2+n*3],*deltaT,rvAcfOrg[0][0]); */
	for (k = 0; k < 3; k ++){ 
             rvAcfOrg[(iType[n]-1)*(*nBuffAcf)+nb][3*n+k]= rv[k+3+n*6] -
	       0.5 * ra[k+n*3] * (*deltaT);}
	}/*printf("Done initializing rvAcfOrg.\n");*/
      }
      
      ni = indexAcf[nb]; 
         
      for(k=0;k<*nType;k++)
         rvAcf[k*(*nBuffAcf)+nb][ni]= 0.;
      for (n = 0; n < *nAtom; n ++) {
	/*printf("%d\t%d\t%lf\t%lf\t%lf\t%lf\n",n,iType[n]-1,rv[5+n*6],ra[2+n*3],*deltaT,rvAcfOrg[0][n]);*/
        for (k = 0; k < 3; k ++) 
           rvAcf[(iType[n]-1)*(*nBuffAcf)+nb][ni]= rvAcf[(iType[n]-1)*(*nBuffAcf)+nb][ni] + rvAcfOrg[(iType[n]-1)*(*nBuffAcf)+nb][3*n + k] *
              (rv[n*6+k+3] - 0.5 * ra[k+n*3] * (*deltaT));
      }
      viscAcf[nb][ni]= thermAcf[nb][ni] = 0.;
      for (k = 0; k < 3; k ++) {
        viscAcf[nb][ni] = viscAcf[nb][ni] + viscAcfOrg[nb][k] * viscVec[k];
        thermAcf[nb][ni] = thermAcf[nb][ni] + thermAcfOrg[nb][k] * thermVec[k];
      }}
    for (nb = 0; nb < *nBuffAcf; nb ++) {
      if (indexAcf[nb] == (*nValAcf-1)) {
        for (j = 0; j < *nValAcf; j ++) {
          for(k=0;k<*nType;k++)
             rvAcfAv[k][j] = rvAcfAv[k][j] + rvAcf[(k*(*nBuffAcf)+nb)][j];
          viscAcfAv[j] = viscAcfAv[j] + viscAcf[nb][j];
          thermAcfAv[j] = thermAcfAv[j] + thermAcf[nb][j];
        }
        indexAcf[nb] =-1;    
        countAcfAv = countAcfAv + 1;
        if (countAcfAv == *limitAcfAv) {
          countAcfAv=0;
          for(k=0;k<*nType;k++){
          fac = 1. / (3 * totType[k] * (*limitAcfAv));
          for(i=0;i<*nValAcf;i++)
	    tmpAcfAv[i]=rvAcfAv[k][i];
          AcfInt[k] = 0.0001 * fac * (*deltaT) * Integrate (tmpAcfAv, *nValAcf);
          for (j = 1; j < *nValAcf; j ++) rvAcfAv[k][j] = rvAcfAv[k][j] / rvAcfAv[k][0];
          rvAcfAv[k][0] = 1.;
          }
          temp=CTemp/(nave-cnave);
          fac = *dens / (3. * (temp)*(*nAtom) * (*limitAcfAv));
          AcfInt[*nType] = 0.19257971014 * fac * (*deltaT) * Integrate (viscAcfAv, (*nValAcf));
          for (j = 1; j < *nValAcf; j ++)
             viscAcfAv[j] = viscAcfAv[j] / viscAcfAv[0];
          viscAcfAv[0] = 1.;
          fac = *dens / (3. * Sqr(temp)*(*nAtom) * (*limitAcfAv));
          AcfInt[*nType+1] = 0.1380658 * fac * (*deltaT) * Integrate (thermAcfAv, *nValAcf);
          for (j = 1; j < *nValAcf; j ++)
             thermAcfAv[j] = thermAcfAv[j] / thermAcfAv[0];
          thermAcfAv[0] = 1.;
          fp=fopen(fle,"a");
          fprintf (fp, "acf\n");
          CTemp=0.;
          cnave=nave;
          for (j = 0; j < *nValAcf; j ++) {
             tVal = (j) * (*deltaT);
             fprintf (fp, "%8.4f", tVal);
             for(k=0;k<*nType;k++)
               fprintf (fp,"  %8.4f", rvAcfAv[k][j]);
             fprintf(fp,"%8.4f  %8.4f\n", viscAcfAv[j], thermAcfAv[j]); 
          }
          printf("Acf Integrals:\t");
          fprintf (fp, "acf integrals:\t");
          for(k=0;k<(*nType+2);k++){
	    fprintf (fp, "%.8e\t",AcfInt[k]);
            printf ("%.8e\t",AcfInt[k]);
	  }
          printf("%lf\n",temp);
          fprintf(fp,"  %lf\n",temp); fclose(fp);
          for(j=0;j<*nValAcf;j++){
            for(k=0;k<*nType;k++)
	       rvAcfAv[k][j]=0.;
            viscAcfAv[j]=thermAcfAv[j]=0.;
          }
       }
    }
  }
}

double Sqr (double number){
  double sq;
  sq=number*number;
  return sq;
}

double CalcTemp (double rv[],int itype[],double amass[], int natoms)
{
  double ke, tmp, vvsum, vv;
  int k, n;
  vvsum=0.0;
  for(n=0;n<natoms;n++){
    vv=0.0;
    for (k=0;k<3;k++){
      vv=vv+amass[itype[n]-1]*Sqr(rv[k+3+n*6]);
    }
    vvsum=vvsum+vv;
  }
  ke=0.5*vvsum/natoms;
  tmp=2.*ke/3.;
  return tmp;
}








