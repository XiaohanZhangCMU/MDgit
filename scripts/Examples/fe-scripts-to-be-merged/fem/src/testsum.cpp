#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#define NMAX 1000000

int main(void)
{
     int *a, sp;
     double b[10], sum1, sum2;
     int i;

     a = (int *) malloc(sizeof(int)*NMAX);
     for(i=0;i<NMAX;i++)
     {
         a[i] = (int) floor(drand48()*4);
     }

     b[0] = 1.0;  b[1] = 2.0; b[3] = -3.0; b[4] = 4.0;

     sum1 = 0;
     for(i=0;i<NMAX;i++)
     {
        sum1 += b[a[i]] * b[a[i]];
     }
     printf("sum1 = %20.12e\n",sum1);

     sum2 = 0;
     for(i=0;i<NMAX;i++)
     {
        sp = a[i]; 
        sum2 += b[sp] * b[sp];
     }
     printf("sum2 = %20.12e\n",sum2);

     if(sum1==sum2) return 0;
     else return 1;
} 
          
