/*
  heaziz.cpp
  by Wei Cai  caiwei@stanford.edu
  Last Modified : Aug 22 2011 
  Notes: with Maurice de Koning

  FUNCTION  :  He MD simulation using Aziz Potential
               Phys. Rev. Lett. 74, 1586 (1995) 

  To Do:  1. Add proper cut-off scheme (constant + linear shift)
*/

#include "heaziz.h"

#define _EFFECTIVE_POTENTIAL

void HEAZIZFrame::initvars()
{
    _RLIST=AZIZ_RC+1.1;
    _SKIN=_RLIST-AZIZ_RC;
    INFO("_RLIST="<<_RLIST<<"  AZIZ_RC="<<AZIZ_RC);
    MDPARALLELFrame::initvars();
}

void HEAZIZFrame::he_aziz()
{/*
   U= potential
   Vir= -1/r(dU/dr)
   F= 1/r^2(d^2U/dr^2-1/r(dU/dr))

   x = r / rm
   U(x) = epsilon*(A exp(-alpha*x+beta*x^2) - (C6/x^6 + C8/x^8 + C10/x^10)*H(x))
   
   H(x) = exp( -(D/x-1)^2 )  x <= D
          1                  x >= D
 */

    int i,j,ipt,jpt;
    double U,r,r2,ri6,x,H,dHdx,Dx,xi,xi2,xi6,xi8,xi10,dUdx,expterm,polyterm,dUdr;
    Vector3 sij, rij, fij;
    DUMP(HIG"Lennard Jones"NOR);
        
    refreshneighborlist();
    
    _EPOT=0;

    for(i=0;i<_NP;i++)
    {_F[i].clear(); _EPOT_IND[i]=0;}
    _VIRIAL.clear();

    for(ipt=0;ipt<_NP;ipt++)
    {
        for(j=0;j<nn[ipt];j++)
        {
            jpt=nindex[ipt][j];
            if(ipt>jpt) continue;
            sij=_SR[jpt]-_SR[ipt];
            sij.subint();
            rij=_H*sij;
            r2=rij.norm2();
            r=sqrt(r2);
            if(r<=AZIZ_RC)
            {
#ifndef _EFFECTIVE_POTENTIAL
                x = r / _RM;
                xi = 1.0/x;
                xi2 = xi*xi;
                xi6 = xi2*xi2*xi2;
                xi8 = xi6*xi2;
                xi10= xi8*xi2;
                if (x >= _D) {
                   H = 1.0;
                   dHdx = 0.0;
                }
                else {
                   Dx = _D/x - 1;
                   H = exp( - Dx*Dx );
                   dHdx = 2.0*H*Dx*_D*xi2;
                } 
                expterm  = exp(-_ALPHA*x+_BETA*x*x);
                polyterm = (_C6*xi6 + _C8*xi8 + _C10*xi10);
                U = _A*expterm - polyterm*H*1;

                dUdx = _A*expterm*(-_ALPHA+2.0*_BETA*x) 
                     + (_C6*6.0*xi6+_C8*8.0*xi8+_C10*10.0*xi10)*xi*H - polyterm*dHdx;
                U *= _ESCALE;
                dUdx *= _ESCALE;
                dUdr = dUdx / _RM;

#else /* use _EFFECTIVE_POTENTIAL */
                double _EA, _EALPHA, ri, riA;
                /*_EA = 4.0;  _EALPHA = 9; */ /* this gives P = 3.86 MPa at a = 3.577 A and T = 0 K */
                _EA = 2.59; _EALPHA = 9;  /* this gives P = 2.50 MPa at a = 3.577 A and T = 0 K */
                /*_EA = 8.66; _EALPHA = 10;*/ /* this gives P = 2.50 MPa at a = 3.577 A and T = 0 K */
                /*_EA = 92.32; _EALPHA = 12;*/ /* this gives P = 2.50 MPa at a = 3.577 A and T = 0 K */
                /*_EA = 1012.5; _EALPHA = 14;*/ /* this gives P = 2.50 MPa at a = 3.577 A and T = 0 K */
                ri = 1/r;
                riA = pow(ri,_EALPHA);
                U = _EA*riA;
                dUdr = -_EA*_EALPHA*riA*ri;
#endif

                fij = rij * (-dUdr/r);

                _F[ipt]-=fij;
                _F[jpt]+=fij;
                _EPOT_IND[ipt]+=U*0.5;
                _EPOT_IND[jpt]+=U*0.5;
                _EPOT+=U;
                _VIRIAL.addnvv(1.,fij,rij);
            }
        }
    }
}


void HEAZIZFrame::potential()
{
    he_aziz();
}


/* Main Program Begins */
#ifdef _TEST
class HEAZIZFrame sim;

/* The main program is defined here */
#include "main.cpp"

#endif//_TEST


    


