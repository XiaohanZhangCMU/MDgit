/*
  meam-lenosky.cpp
  by Wei Cai  caiwei@mit.edu
  Last Modified : Mon Jan  1 21:27:48 2007

  FUNCTION  : MD++ with Lenosky potential for Si

  Note: The Lenosky potential is for Si only
        Only rectangular supercell is allowed
        Virial stress is not computed
*/

#include "meam-lenosky.h"

void MEAMFrame::potential()
{
    MEAM();
}

void MEAMFrame::MEAM()
{
    double alat[3], coord, ener_var, coord_var, eps;

    alat[0]=_H[0][0];
    alat[1]=_H[1][1];
    alat[2]=_H[2][2];

    SHtoR();

    eps=1e-8;
    if ( (fabs(_H[0][1])>eps)||(fabs(_H[0][2])>eps)||
         (fabs(_H[1][0])>eps)||(fabs(_H[1][2])>eps)||
         (fabs(_H[2][0])>eps)||(fabs(_H[2][1])>eps) )
    {
        WARNING(" Only rectangular supercell is allowed!");
    }
    lenosky_(&_NP, alat, (double *)_R, (double *)_F, _EPOT_IND, &_EPOT,
             &coord, &ener_var, &coord_var, &count);
}

void MEAMFrame::Alloc()
{
    MDPARALLELFrame::Alloc();
}    

void MEAMFrame::initvars()
{
    MDPARALLELFrame::initvars();
}

void MEAMFrame::initparser()
{
    MDPARALLELFrame::initparser();
}

int MEAMFrame::exec(const char *name)
{
    if(MDPARALLELFrame::exec(name)==0) return 0;
    return -1;
}
   
#ifdef _TEST

/* Main Program Begins */
class MEAMFrame sim;

/* The main program is defined here */
#include "main.cpp"

#endif//_TEST






























































/* old code */
#if 0

inline double interp(double func[],double deriv[],double dr,int ind,double qq)
{
    double f, a, b, c, d, f1, p1, A1, A2, dr2, dr3, qq2, qq3;
//    f = func[ind] + qq*deriv[ind];
    dr2=dr*dr; dr3=dr2*dr;
    qq2=qq*qq; qq3=qq2*qq;
    a = func[ind];
    b = deriv[ind];
    f1 = func[ind+1];
    p1 = deriv[ind+1];
    A1 = f1-a-b*dr;
    A2 = (p1-b)*dr;
    d = (A2-2*A1)/dr3;
    c = (3*A1-A2)/dr2;
    f=a+b*qq+c*qq2+d*qq3;
    return f;
}
    
    
inline double interp1(double func[],double deriv[],double dr,int ind,double qq)
{
    double fp, a, b, c, d, f1, p1, A1, A2, dr2, dr3, qq2, qq3;
//    fp = deriv[ind] + qq1/dr*(deriv[ind+1]-deriv[ind]);
    dr2=dr*dr; dr3=dr2*dr;
    qq2=qq*qq; qq3=qq2*qq;
    a = func[ind];
    b = deriv[ind];
    f1 = func[ind+1];
    p1 = deriv[ind+1];
    A1 = f1-a-b*dr;
    A2 = (p1-b)*dr;
    d = (A2-2*A1)/dr3;
    c = (3*A1-A2)/dr2;
    fp=b+2*c*qq+3*d*qq2;
    return fp;
}
#endif
