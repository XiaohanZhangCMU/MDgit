/*
  lj2.cpp
  by Wei Cai  caiwei@mit.edu
  Last Modified : Mon Jan  1 21:26:52 2007

  FUNCTION  :  MD simulation package of Si using Stillinger-Weber potential
               and Lennard-Jones potential (depending on group ID)
*/

#include "lj2.h"

void LJ2Frame::initparser()
{
    MDPARALLELFrame::initparser();

    /* input */
    bindvar("C12_00",&_ALJ_00,DOUBLE);
    bindvar("C6_00",&_BLJ_00,DOUBLE);

    bindvar("C12_11",&_ALJ_11,DOUBLE);
    bindvar("C6_11",&_BLJ_11,DOUBLE);

    bindvar("C12_01",&_ALJ_01,DOUBLE);
    bindvar("C6_01",&_BLJ_01,DOUBLE);
    
}

int LJ2Frame::exec(const char *name)
{
    if(MDPARALLELFrame::exec(name)==0) return 0;
    bindcommand(name,"initLJ",initLJ());
    return -1;
}

void LJ2Frame::initvars()
{
    _RLIST=LJ_RC+1.1;
    _SKIN=_RLIST-LJ_RC;
    MDPARALLELFrame::initvars();
}

void LJ2Frame::initLJ()
{
    ALJ_00=_ALJ_00*LJ_ENERGY*POW12(LJ_LENGTH);
    BLJ_00=_BLJ_00*LJ_ENERGY*POW6(LJ_LENGTH);
    Uc_00=(ALJ_00*POW6(1/LJ_RC)-BLJ_00)*POW6(1/LJ_RC);
    DUDRc_00=-(12*ALJ_00/POW6(LJ_RC)-6*BLJ_00)/POW6(LJ_RC)/LJ_RC;

    ALJ_11=_ALJ_11*LJ_ENERGY*POW12(LJ_LENGTH);
    BLJ_11=_BLJ_11*LJ_ENERGY*POW6(LJ_LENGTH);
    Uc_11=(ALJ_11*POW6(1/LJ_RC)-BLJ_11)*POW6(1/LJ_RC);
    DUDRc_11=-(12*ALJ_11/POW6(LJ_RC)-6*BLJ_11)/POW6(LJ_RC)/LJ_RC;

    ALJ_01=_ALJ_01*LJ_ENERGY*POW12(LJ_LENGTH);
    BLJ_01=_BLJ_01*LJ_ENERGY*POW6(LJ_LENGTH);
    Uc_01=(ALJ_01*POW6(1/LJ_RC)-BLJ_01)*POW6(1/LJ_RC);
    DUDRc_01=-(12*ALJ_01/POW6(LJ_RC)-6*BLJ_01)/POW6(LJ_RC)/LJ_RC;

    INFO_Printf("ALJ_00 = %e BLJ_00 = %e Uc_00 = %e DUDRc_00 = %e\n",ALJ_00,BLJ_00,Uc_00,DUDRc_00);
    INFO_Printf("ALJ_11 = %e BLJ_11 = %e Uc_11 = %e DUDRc_11 = %e\n",ALJ_11,BLJ_11,Uc_11,DUDRc_11);
    INFO_Printf("ALJ_01 = %e BLJ_01 = %e Uc_01 = %e DUDRc_01 = %e\n",ALJ_01,BLJ_01,Uc_01,DUDRc_01);
}

void LJ2Frame::calcprop()
{
    MDPARALLELFrame::calcprop();
}

void LJ2Frame::lennard_jones_2()
{
/*
   U= potential
   Vir= -1/r(dU/dr)
   F= 1/r^2(d^2U/dr^2-1/r(dU/dr))
 */

    /* !!! Parallel version not implemented !!! */
    int i,j,ipt,jpt;
    double U,r,r2,ri6;
    double ALJ,BLJ,Uc,DUDRc;   
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
            if(r<=LJ_RC)
            {
                if((species[ipt]==0)&&(species[jpt]==0))
                {
                    ALJ=ALJ_00; BLJ=BLJ_00; Uc=Uc_00; DUDRc=DUDRc_00;
                }
                else if((species[ipt]==1)&&(species[jpt]==1))
                {
                    ALJ=ALJ_11; BLJ=BLJ_11; Uc=Uc_11; DUDRc=DUDRc_11;
                }
                else
                {
                    ALJ=ALJ_01; BLJ=BLJ_01; Uc=Uc_01; DUDRc=DUDRc_01;
                }

                //INFO_Printf("ALJ = %e BLJ = %e Uc = %e DUDRc = %e\n",ALJ,BLJ,Uc,DUDRc);
                
                ri6=1./(r2*r2*r2);
                
                U=(ALJ*ri6-BLJ)*ri6-r*DUDRc+(LJ_RC*DUDRc-Uc);
                fij=rij*((12*ALJ*ri6-6*BLJ)*ri6/r2+DUDRc/r);

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


void LJ2Frame::potential()
{
    lennard_jones_2();
}

#ifdef _TEST

/* Main Program Begins */
class LJ2Frame sim;

/* The main program is defined here */
#include "main.cpp"

#endif//_TEST

