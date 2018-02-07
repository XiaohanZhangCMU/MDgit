/*
  lj.cpp
  by Wei Cai  caiwei@mit.edu
  Last Modified : Mon Jan  1 21:27:02 2007

  FUNCTION  :  Argon MD simulation using Lennard-Jones Potential
*/

#include "lj.h"


void LJFrame::initvars()
{
    _RLIST=LJ_RC+1.1;
    _SKIN=_RLIST-LJ_RC;
//    _RLIST = 4.2;
//    _SKIN = 0.4;
    INFO("_RLIST="<<_RLIST<<"  LJ_RC="<<LJ_RC);
    MDPARALLELFrame::initvars();
}

void LJFrame::lennard_jones()
{/*
   U= potential
   Vir= -1/r(dU/dr)
   F= 1/r^2(d^2U/dr^2-1/r(dU/dr))
 */

    /* !!! Parallel version not implemented !!! */
    int i,j,ipt,jpt;
    double U,r,r2,ri6;
    Vector3 sij, rij, fij;
    DUMP(HIG"Lennard Jones"NOR);
        
    refreshneighborlist();
    
    _EPOT=0;

    for(i=0;i<_NP;i++)
    {_F[i].clear(); _EPOT_IND[i]=0;}
    _VIRIAL.clear();

    for(ipt=0;ipt<_NP;ipt++)
    {
//        INFO("ipt="<<ipt<<"   nn="<<nn[ipt]);
        for(j=0;j<nn[ipt];j++)
        {
            jpt=nindex[ipt][j];
            if(ipt>jpt) continue;
//            INFO("ipt="<<ipt<<"  jpt="<<jpt);
            sij=_SR[jpt]-_SR[ipt];
            sij.subint();
            rij=_H*sij;
            r2=rij.norm2();
            r=sqrt(r2);
            if(r<=LJ_RC)
            {
                ri6=1./(r2*r2*r2);
                U=(A*ri6-B)*ri6-r*DUDRc+(LJ_RC*DUDRc-Uc);
//                INFO("U  "<<U<<"  r="<<r);
                fij=rij*((12*A*ri6-6*B)*ri6/r2+DUDRc/r);
                
                _F[ipt]-=fij;
                _F[jpt]+=fij;
                _EPOT_IND[ipt]+=U*0.5;
                _EPOT_IND[jpt]+=U*0.5;
                _EPOT+=U;
                _VIRIAL.addnvv(1.,fij,rij);
            }
        }
    }
//    INFO("F[0]="<<_F[0]);
}


void LJFrame::potential()
{
    lennard_jones();
}


/* Main Program Begins */
#ifdef _TEST
class LJFrame sim;

/* The main program is defined here */
#include "main.cpp"

#endif//_TEST


    


