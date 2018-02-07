/**************************************************************************
 *
 *      Module:  This module contains the functions needed for
 *               calculating the solid angle of an ellipse.
 *
 *      To test:
 *
 *        x1=rand; y1=rand; z1=rand; x2=rand; y2=rand; z2=rand; x3=rand; y3=rand; z3=rand; xp=rand; yp=rand; zp=rand;
 *        omega1=SolidAngleTriangle  (x1,y1,z1,x2,y2,z2,x3,y3,z3,xp,yp,zp)
 *      Compare with Matlab implementation: 
 *        omega2=solid_angle_triangle([x1 y1 z1],[x2 y2 z2],[x3 y3 z3],xp,yp,zp)
 *
 *************************************************************************/
#include <math.h>
#include <mex.h>

static void SolidAngleTriangle(double x1, double y1, double z1,
                               double x2, double y2, double z2,
                               double x3, double y3, double z3,
                               double xp, double yp, double zp,
                               double *Omega);

void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
    {
    
    double x1,y1,z1,x2,y2,z2,x3,y3,z3,xp,yp,zp; 
    double Omega;
    int i;
    
    /* Check for proper number of input arguments. */
    if (nrhs!=12) mexErrMsgTxt("12 inputs required!");
    if (nlhs!=1 ) mexErrMsgTxt("1 outputs required!");
        
    x1=mxGetScalar(prhs[0]);
    y1=mxGetScalar(prhs[1]);
    z1=mxGetScalar(prhs[2]);
    x2=mxGetScalar(prhs[3]);
    y2=mxGetScalar(prhs[4]);
    z2=mxGetScalar(prhs[5]);
    x3=mxGetScalar(prhs[6]);
    y3=mxGetScalar(prhs[7]);
    z3=mxGetScalar(prhs[8]);
    xp=mxGetScalar(prhs[9]);
    yp=mxGetScalar(prhs[10]);
    zp=mxGetScalar(prhs[11]);
    
    SolidAngleTriangle(x1,y1,z1,x2,y2,z2,x3,y3,z3,xp,yp,zp,&Omega);
                
     for (i=0; i<nlhs; i++)  {
	     plhs[i]=mxCreateDoubleMatrix(1,1,mxREAL);}
	 *mxGetPr(plhs[0])=Omega;
   
    }


/*-------------------------------------------------------------------------
 *
 *      Function:       SolidAngleTriangle
 *      Description:    Calculate the solid angle of a triangle.
 *
 *      Arguments:
 *              x1,y1,z1     First  vertex of the triangle
 *              x2,y2,z2     Second vertex of the triangle
 *              x3,y3,z3     Third  vertex of the triangle
 *              xp,yp,zp     Point at which solid angle is computed
 *              Omega*       pointers to locations in which to return solid angle
 *-----------------------------------------------------------------------*/
void SolidAngleTriangle(double x1, double y1, double z1,
                        double x2, double y2, double z2,
                        double x3, double y3, double z3, 
                        double xp, double yp, double zp,
                        double *Omega)
{
    double xs1,ys1,zs1,xs2,ys2,zs2,xs3,ys3,zs3;
    double RR1,RR2,RR3,dotR1R2,dotR2R3,dotR3R1;
    double numer,denom;
    
    xs1 = x1 - xp; ys1 = y1 - yp; zs1 = z1 - zp;
    xs2 = x2 - xp; ys2 = y2 - yp; zs2 = z2 - zp;
    xs3 = x3 - xp; ys3 = y3 - yp; zs3 = z3 - zp;

    RR1 = sqrt(xs1*xs1 + ys1*ys1 + zs1*zs1);
    RR2 = sqrt(xs2*xs2 + ys2*ys2 + zs2*zs2);
    RR3 = sqrt(xs3*xs3 + ys3*ys3 + zs3*zs3);
    
    numer = xs1*ys2*zs3 + ys1*zs2*xs3 + zs1*xs2*ys3
           -xs1*zs2*ys3 - ys1*xs2*zs3 - zs1*ys2*xs3;

    dotR1R2 = xs1*xs2 + ys1*ys2 + zs1*zs2;
    dotR2R3 = xs2*xs3 + ys2*ys3 + zs2*zs3;
    dotR3R1 = xs3*xs1 + ys3*ys1 + zs3*zs1;
    
    denom = RR1*RR2*RR3 + dotR1R2*RR3 + dotR2R3*RR1 + dotR3R1*RR2;
  
    *Omega = -2 * atan2(numer, denom);
}


