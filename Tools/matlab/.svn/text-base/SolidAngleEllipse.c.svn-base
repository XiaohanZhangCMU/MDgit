/**************************************************************************
 *
 *      Module:  This module contains the functions needed for
 *               calculating the solid angle of an ellipse.
 *
 *      To test:
 *
 *        a=rand; b=rand; x=rand; y=rand; z=rand;
 *        omega1=SolidAngleEllipse  (a,b,x,y,z)
 *      Compare with Matlab implementation: 
 *        omega2=solid_angle_ellipse(a,b,x,y,z)
 *
 *************************************************************************/
#include <math.h>
#include <mex.h>

static void SolidAngleEllipse(double a, double b, double x, double y, double z,
                              double *Omega);

void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
    {
    
    double a,b,x,y,z; 
    double Omega;
    int i;
    
    /* Check for proper number of input arguments. */
    if (nrhs!=5) mexErrMsgTxt("5 inputs required!");
    if (nlhs!=1) mexErrMsgTxt("1 outputs required!");
    
    
    a=mxGetScalar(prhs[0]);
    b=mxGetScalar(prhs[1]);
    x=mxGetScalar(prhs[2]);
    y=mxGetScalar(prhs[3]);
    z=mxGetScalar(prhs[4]);
    
    SolidAngleEllipse(a,b,x,y,z,&Omega);
                
     for (i=0; i<nlhs; i++)  {
	     plhs[i]=mxCreateDoubleMatrix(1,1,mxREAL);}
	 *mxGetPr(plhs[0])=Omega;
   
    }


/*-------------------------------------------------------------------------
 *
 *      Function:       SolidAngleEllipse
 *      Description:    Calculate the solid angle of an ellipse.
 *
 *      Arguments:
 *              a,b          Axis of the ellipse (x/a)^2+(y/b)^2=1
 *              x,y,z        Point at which solid angle is computed
 *              Omega*       pointers to locations in which to return solid angle
 *-----------------------------------------------------------------------*/
void SolidAngleEllipse(double a, double b, double x, double y, double z,
                       double *Omega)
{
        int     Nint, maxiter, iter, j, converged;
        double  eps;
        double  omega, oldomega, theta, dtheta, weight, integrand, integral;
        double  dx, dy1, dy2, r1, r2;
        
        eps = 1e-4;            

        if (fabs(z)<eps)
        {
            if ( (x*x/a/a+y*y/b/b)<1 )
                *Omega = ((z>0)?(1):(-1))*(M_PI*2);
            else
                *Omega = 0;
            return;
        }
        
        omega = 0.0;
        Nint = 4;
        maxiter = 10;
        for (iter=0;iter<maxiter;iter++)
        {
            if (iter>0) oldomega = omega;
            else oldomega = 0.0;
            
            integral = 0.0;
            for (j=0;j<=Nint;j++)
            {
                if ((j==0)||(j==Nint)) weight = 0.5;
                else weight = 1.0;
                
                theta = ((j*1.0)/Nint - 0.5)*M_PI;
                dx=a*sin(theta)-x; 
                dy1=b*cos(theta)-y; dy2=-b*cos(theta)-y;
                r1=sqrt(dx*dx+dy1*dy1+z*z); r2=sqrt(dx*dx+dy2*dy2+z*z); 
                integrand = cos(theta)/(dx*dx+z*z)*( dy1/r1 - dy2/r2 );
                integral += integrand * weight;
            }
            dtheta = M_PI/Nint;
            integral *= dtheta;
            omega = a*z*integral;
            
            Nint = Nint*4;
    
            /* test convergence */
            converged = 0;
            if (iter>0)
              if (fabs(oldomega-omega) < fabs(eps*omega))
                  converged = 1;
            if (converged)
                break;
    
            if (iter==maxiter)
                fprintf(stderr,"error: solid_angle_ellipse: maxiter exceeded\n");
        }/* end of for(iter) */
       *Omega = omega;
}


