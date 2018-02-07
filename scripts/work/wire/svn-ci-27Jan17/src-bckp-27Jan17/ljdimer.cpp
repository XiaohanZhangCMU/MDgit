/*
  ljdimer.cpp
  by Wei Cai  caiwei@mit.edu
  Last Modified : Fri Feb 23 18:32:33 2007

  FUNCTION  :  Lennard-Jones potential and dimer
*/

#include "ljdimer.h"

void LJDIMERFrame::initparser()
{
    MDPARALLELFrame::initparser();

    /* input */
    bindvar("C12_00",&_ALJ_00,DOUBLE);
    bindvar("C6_00",&_BLJ_00,DOUBLE);
    bindvar("R_DIMER",&R_DIMER,DOUBLE);
    bindvar("F_DIMER",&F_DIMER,DOUBLE);
    bindvar("F_DIMER_INTERNAL",&F_DIMER_INTERNAL,DOUBLE);
    bindvar("F_DIMER_EXTERNAL",&F_DIMER_EXTERNAL,DOUBLE);
}

int LJDIMERFrame::exec(const char *name)
{
    if(MDPARALLELFrame::exec(name)==0) return 0;
    bindcommand(name,"initLJ",initLJ());
    return -1;
}

void LJDIMERFrame::initvars()
{
    _RLIST=LJ_RC+1.1;
    _SKIN=_RLIST-LJ_RC;
    MDPARALLELFrame::initvars();
}

void LJDIMERFrame::initLJ()
{
    ALJ_00=_ALJ_00*LJ_ENERGY*POW12(LJ_LENGTH);
    BLJ_00=_BLJ_00*LJ_ENERGY*POW6(LJ_LENGTH);
    Uc_00=(ALJ_00*POW6(1/LJ_RC)-BLJ_00)*POW6(1/LJ_RC);
    DUDRc_00=-(12*ALJ_00/POW6(LJ_RC)-6*BLJ_00)/POW6(LJ_RC)/LJ_RC;

    H_DIMER=_H_DIMER*LJ_ENERGY;
    W_DIMER=_W_DIMER*LJ_LENGTH;

    INFO_Printf("ALJ_00 = %e BLJ_00 = %e Uc_00 = %e DUDRc_00 = %e\n",ALJ_00,BLJ_00,Uc_00,DUDRc_00);
    INFO_Printf("H_DIMER = %e W_DIMER = %e\n",H_DIMER,W_DIMER);
}

void LJDIMERFrame::calcprop()
{
    MDPARALLELFrame::calcprop();
}

void LJDIMERFrame::find_dimer_indices(int *ind0,int *ind1)
{
    int ipt, i0, i1;
    /* find the two atoms (species[i]==1) in the dimer */
    i0 = -1; i1 = -1;
    for(ipt=0;ipt<_NP;ipt++)
    {
        if(species[ipt]==1)
        {
            if(i0<0)
            { /* first atom */
                i0=ipt;
            }
            else if(i1<0)
            { /* second atom */
                i1=ipt;
            }
            else
            {
                FATAL("more than two atoms with species 1");
            }
        }
    }
    *ind0 = i0;
    *ind1 = i1;
}

void LJDIMERFrame::lennard_jones_dimer()
{
/*
   U= potential
   Vir= -1/r(dU/dr)
   F= 1/r^2(d^2U/dr^2-1/r(dU/dr))
 */

    int i,j,ipt,jpt, i0, i1;
    double U,r,r2,ri6,tmp;
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
                if((species[ipt]==0)||(species[jpt]==0))
                {/* any interaction involving solvent is LJ */
                    ALJ=ALJ_00; BLJ=BLJ_00; Uc=Uc_00; DUDRc=DUDRc_00;
                    
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
    /* find the two atoms (species[i]==1) in the dimer */
    find_dimer_indices(&i0,&i1);
//    INFO_Printf("i0=%d i1=%d\n",i0,i1);

    if((i0>=0)&&(i1>=0))
    {
        /* internal interactions within dimer */
        sij=_SR[i1]-_SR[i0];
        sij.subint();
        rij=_H*sij;
        r2=rij.norm2();
        r=sqrt(r2);
        R_DIMER = r;
        
        tmp = 1 - SQR((r-LJ_RC-W_DIMER)/W_DIMER);    
        U = H_DIMER * SQR( tmp );
        F_DIMER_INTERNAL = 2*H_DIMER* tmp * ( 2*(r-LJ_RC-W_DIMER)/ SQR(W_DIMER) );
        fij = rij *( F_DIMER_INTERNAL / r); 
        
        _F[i0]-=fij;
        _F[i1]+=fij;
        _EPOT_IND[i0]+=U*0.5;
        _EPOT_IND[i1]+=U*0.5;
        _EPOT+=U;
        _VIRIAL.addnvv(1.,fij,rij);
    }
}

void LJDIMERFrame::lennard_jones_dimer_constrained(double R)
{
/* the distance between the two atoms in the dimer is constrained at R
 * before evaluating force: make sure the length of the dimer is exactly at R
 * after  evaluating force: project out any force component that will change R
 */

    int i0, i1;
    double r, vn, fn;
    Vector3 sij, rij, snij, nij, svij, vij, fij;
    
    /* constrain length of dimer to r */
    find_dimer_indices(&i0,&i1);
    sij = _SR[i1]-_SR[i0];
    sij.subint(); rij = _H*sij; r = rij.norm();
    //INFO_Printf("before: R = %g r = %g i0 = %d i1 = %d\n",R,r,i0,i1);

    if((r!=R)&&(r>0))
    { /* find center of mass position */
      _SR[i1] += sij * ((R/r-1)/2);
      _SR[i0] -= sij * ((R/r-1)/2);
    }

    if(r==0) FATAL("lennard_jones_dimer_constrained: r=0!");

    sij = _SR[i1]-_SR[i0];
    sij.subint(); rij = _H*sij; r = rij.norm(); nij = rij/r; snij = sij/r;
    //INFO_Printf("after:  R = %g r = %g\n",R,r);

    /* project out velocity component that will change r */
    svij = _VSR[i1] - _VSR[i0]; vij = _H*svij;
    vn = dot(vij,nij); 
    //INFO_Printf("before: vn = %g\n",vn);

    _VSR[i1] -= snij*(vn/2);
    _VSR[i0] += snij*(vn/2);
    svij = _VSR[i1] - _VSR[i0]; vij = _H*svij;
    vn = dot(vij,nij); 
    //INFO_Printf("after:  vn = %g\n",vn);
 
    lennard_jones_dimer();

    /* project out force component that will change r */
    fij = _F[i1] - _F[i0]; fn = dot(fij,nij); 
    F_DIMER = fn/2;  F_DIMER_EXTERNAL = F_DIMER - F_DIMER_INTERNAL;
    dEdlambda = F_DIMER_EXTERNAL;
    //INFO_Printf("before: fn = %g\n",fn);

    _F[i1] -= nij*(fn/2);
    _F[i0] += nij*(fn/2);
    //fij = _F[i1] - _F[i0]; fn = dot(fij,nij); 
    //INFO_Printf("after:  fn = %g\n",fn);
}

void LJDIMERFrame::potential()
{
    lennard_jones_dimer();
}

void LJDIMERFrame::SWITCHpotential_user(double lambda)
{
    lennard_jones_dimer_constrained(lambda);
}


#ifdef _TEST

/* Main Program Begins */
class LJDIMERFrame sim;

/* The main program is defined here */
#include "main.cpp"

#endif//_TEST

