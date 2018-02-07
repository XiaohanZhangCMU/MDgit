/*
  meam-marian.cpp
  by Wei Cai  caiwei@mit.edu
  Last Modified : Mon Jan  1 21:28:21 2007

  FUNCTION  : MD++ with MEAM potential (using MDCASK F90 Codes from Jaime Marian)
*/

#include "meam-marian.h"


void MEAMFrame::Alloc()
{
    MDPARALLELFrame::Alloc();
}    

void MEAMFrame::potential()
{
    MEAM();
}

void MEAMFrame::readMEAM()
{
    meamdefine_();
    _RLIST=meamread_.rcutmeam;
    _SKIN=0;
}

void MEAMFrame::MEAM()
{
    Vector3 h;
    int i;
    double eps;
    
    /*refreshneighborlist();*/
    h=_H.height();
    linklist_.nlcx=(int)floor(h[0]/((_RLIST)*1.05));if(linklist_.nlcx==0)linklist_.nlcx=1;
    linklist_.nlcy=(int)floor(h[1]/((_RLIST)*1.05));if(linklist_.nlcy==0)linklist_.nlcy=1;
    linklist_.nlcz=(int)floor(h[2]/((_RLIST)*1.05));if(linklist_.nlcz==0)linklist_.nlcz=1;
    linklist_.nlc=linklist_.nlcx*linklist_.nlcy*linklist_.nlcz;
    
    //INFO_Printf("nlcx=%d nlcy=%d nlcz=%d\n",linklist_.nlcx,linklist_.nlcy,linklist_.nlcz);
    
    np_.nm=_NP;
    for(i=0;i<_NP;i++)
    {
        adim_.x0[i]=_SR[i].x;
        adim_.y0[i]=_SR[i].y;
        adim_.z0[i]=_SR[i].z;
    }
    box_.Lx=_H[0][0];
    box_.Ly=_H[1][1];
    box_.Lz=_H[2][2];

    eps=1e-8;
    if ( (fabs(_H[0][1])>eps)||(fabs(_H[0][2])>eps)||
         (fabs(_H[1][0])>eps)||(fabs(_H[1][2])>eps)||
         (fabs(_H[2][0])>eps)||(fabs(_H[2][1])>eps) )
    {
        WARNING(" Only rectangular supercell is allowed!");
    }
    
    //INFO_Printf("np_.nm=%d\n",np_.nm);
    forces_();

    for(i=0;i<_NP;i++)
    {
        _F[i].set(area6_.fx[i],area6_.fy[i],area6_.fz[i]);
        _EPOT_IND[i] = ppote_.atpe[i]+ppote_.atpe3b[i];        
    }

    _EPOT=ppote_.pe+ppote_.petrip;
    //INFO_Printf("p2=%20.12e p3=%20.12e\n",ppote_.pe,ppote_.petrip);
}

void MEAMFrame::initvars()
{
    MDPARALLELFrame::initvars();
}

void MEAMFrame::initparser()
{
    MDPARALLELFrame::initparser();
    bindvar("meamfile",meamdatafile_.meamfile,STRING);
}

int MEAMFrame::exec(const char *name)
{
    if(MDPARALLELFrame::exec(name)==0) return 0;
    bindcommand(name,"readMEAM",readMEAM());
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
