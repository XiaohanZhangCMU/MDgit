/*
  swlj.cpp
  by Wei Cai  caiwei@mit.edu
  Last Modified : Thu Apr 22 16:41:30 2010

  FUNCTION  :  MD simulation package of Si using Stillinger-Weber potential
               and Lennard-Jones potential (depending on group ID)
               group 0 and group 0 interact with SW potential
               group 1 and group 1 interact with SW potential
               group 0 and group 1 interact with LJ potential

               This model is developed for Nano-indentation simulation.
*/

#include "swlj.h"

void SWFrame::initparser()
{
    MDPARALLELFrame::initparser();

    /* input */
    bindvar("C12",&_ALJ,DOUBLE);
    bindvar("C6",&_BLJ,DOUBLE);

    /* output */
    bindvar("tipforce",&(tipforce[0]),DOUBLE);
}

int SWFrame::exec(const char *name)
{
    if(MDPARALLELFrame::exec(name)==0) return 0;
    bindcommand(name,"initLJ",initLJ());
    return -1;
}

void SWFrame::initvars()
{
    /* modified Si parameters PRB 46, 2250 (1992) */
    aa=16.31972277; bb=11.60319228; plam=48.61499998;
    pgam=2.51412;  acut=3.77118; pss=2.0951; rho=4.;

    
    _RLIST=acut*1.1;
    _SKIN=_RLIST-acut;
    rho1=rho+1; acutsq=acut*acut;

    strcpy(incnfile,"../si.cn");
    DUMP(HIG"SWFrame initvars"NOR);
    MDPARALLELFrame::initvars();
}

void SWFrame::initLJ()
{
    ALJ=_ALJ*LJ_ENERGY*POW12(LJ_LENGTH);
    BLJ=_BLJ*LJ_ENERGY*POW6(LJ_LENGTH);
    Uc=(ALJ*POW6(1/LJ_RC)-BLJ)*POW6(1/LJ_RC);
    DUDRc=-(12*ALJ/POW6(LJ_RC)-6*BLJ)/POW6(LJ_RC)/LJ_RC;
}

void SWFrame::calcprop()
{
    MDPARALLELFrame::calcprop();
}

void SWFrame::stillinger_weber()
{
#ifndef NNM /* need to be taken care of later */
#define NNM 60
#endif
    /* Multi process function */
    int i,j,k,ipt,jpt,kpt;
    int neg;
    Vector3 sij,rij,f3i,f3j,f3k,f2;
    Vector3 rg[NNM],rt[NNM];
    double rr[NNM],ri[NNM],rdi[NNM],rdi2[NNM],xpon[NNM],xpon2[NNM];
    double rri6[NNM];
    int list2[NNM];
    double r2ij,rrij;
    double ang,tm1,tm2,tm3,dhij,dhik,dhmu,tmij,tmik,tm2ij,tm2ik,
        tm3ik,tm3ij,eterm,tote3,au,pplm,ff2,eterm2,tote2;
    int n0, n1;
    int group_01;
    
    DUMP("SW");

    _SKIN=_RLIST-acut;
    refreshneighborlist();

    /* _EPOT, _F, _EPOT_IND, _VIRIAL all local */
    _EPOT=0;tote2=tote3=0;
    for(i=0;i<_NP;i++)
    { _F[i].clear(); _EPOT_IND[i]=0;}
    _VIRIAL.clear();

    n0=0;
    n1=_NP;

    for(ipt=n0;ipt<n1;ipt++)
    {        
        /* if fixed == -1, simply ignore this atom */
        if(fixed[ipt]==-1) continue;
        neg=0;
        for(j=0;j<nn[ipt];j++)
        {
            jpt=nindex[ipt][j];
            /* if fixed == -1, simply ignore this atom */
            if(fixed[jpt]==-1) continue;
            
            /*sij=_SR[jpt]-_SR[ipt];*/
            sij.subtract(_SR[jpt],_SR[ipt]);
            sij.subint();
            /*rij=_H*sij;*/
            _H.multiply(sij,rij);
            
            r2ij=rij.norm2();
            rrij=sqrt(r2ij);
            if(r2ij>acutsq) continue;
            rg[neg]=rij;
            /*rt[neg]=rij/rrij;*/
            rt[neg]=rij; rt[neg]/=rrij;
            rr[neg]=rrij;
            ri[neg]=1./rrij;
            rdi[neg]=1./(rrij-acut);
            rdi2[neg]=rdi[neg]*rdi[neg];

            rri6[neg]=1./(r2ij*r2ij*r2ij);
            
            if(fabs(pgam*rdi[neg])>30) xpon[neg]=0;/* avoid denormalization*/
            else xpon[neg]=exp(pgam*rdi[neg]);
            if(fabs(pss*rdi[neg])>30) xpon2[neg]=0;/* avoid denormalization*/
            else xpon2[neg]=exp(pss*rdi[neg]);
            list2[neg]=jpt;
            neg++;
        }
        /* second inner loop */
        for(j=0;j<neg;j++)
        {
            /* three body loop */
            jpt=list2[j];
            /* if fixed == -1, simply ignore this atom */
            if(fixed[jpt]==-1) continue;

            /* between group 0 and 1, no three-body interaction*/
            if(((group[ipt]==0) && (group[jpt]==1))||
               ((group[ipt]==1) && (group[jpt]==0)))
                group_01 = 1;
            else
                group_01 = 0;


            if(!group_01)
            for(k=j+1;k<neg;k++)
            {
                kpt=list2[k];
                /* if fixed == -1, simply ignore this atom */
                if(fixed[kpt]==-1) continue;

                /* between group 0 and 1, no three-body interaction*/
                if(((group[ipt]==0) && (group[kpt]==1))||
                   ((group[ipt]==1) && (group[kpt]==0)))
                    continue;
                
                ang=rg[j]*rg[k]/(rr[j]*rr[k]);
                tm1=ang+1./3.;
                tm2=plam*pgam*tm1*tm1;
                tm3=xpon[j]*xpon[k];

                dhij=-1.*tm2*rdi[j]*rdi[j]*tm3;
                dhik=-1.*tm2*rdi[k]*rdi[k]*tm3;
                dhmu=2.*plam*tm3*tm1;
                    
                tmij=dhij+dhmu*(ri[k]-ang*ri[j]);
                tmik=dhik+dhmu*(ri[j]-ang*ri[k]);
                tm2ij=-dhij+dhmu*ang*ri[j];
                tm2ik=-dhmu*ri[j];
                tm3ik=-dhik+dhmu*ang*ri[k]; 
                tm3ij=-dhmu*ri[k];

                /*f3i=rt[j]*tmij+rt[k]*tmik;*/
                f3i.clear();f3i.addnv(tmij,rt[j]);f3i.addnv(tmik,rt[k]);
                _F[ipt]+=f3i;
                    
                /*f3j=rt[j]*tm2ij+rt[k]*tm2ik;*/
                f3j.clear();f3j.addnv(tm2ij,rt[j]);f3j.addnv(tm2ik,rt[k]);
                _F[jpt]+=f3j;

                /*f3k=rt[k]*tm3ik+rt[j]*tm3ij;*/
                f3k.clear();f3k.addnv(tm3ik,rt[k]);f3k.addnv(tm3ij,rt[j]);
                _F[kpt]+=f3k;
                if(fixedatomenergypartition==0)
                {
                    _VIRIAL.addnvv(1.,f3j,rg[j]);
                    _VIRIAL.addnvv(1.,f3k,rg[k]);
                }
                else
                {
                    if(!(fixed[ipt]||fixed[jpt]))
                        _VIRIAL.addnvv(1.,f3j,rg[j]);
                    else if(!(fixed[ipt]&&fixed[jpt]))
                        _VIRIAL.addnvv(0.5,f3j,rg[j]);
                    if(!(fixed[ipt]||fixed[kpt]))
                        _VIRIAL.addnvv(1.,f3k,rg[k]);
                    else if(!(fixed[ipt]&&fixed[kpt]))
                        _VIRIAL.addnvv(0.5,f3k,rg[k]);
                }
                /*_VIRIAL.addnvv(1.,f3k,rg[k]);*/

                eterm=tm2*tm3/pgam;

                if(fixedatomenergypartition==0)
                {
                    tote3+=eterm;
                }
                else
                {
                    if(!(fixed[ipt]&&fixed[jpt]&&fixed[kpt]))
                        tote3+=eterm;
//                    tote3+=eterm*(1.0-(fixed[ipt]+fixed[jpt]+fixed[kpt])/3.0);
                    /*tote3+=eterm;*/
                }
                _EPOT_IND[ipt]+=eterm/3;
                _EPOT_IND[jpt]+=eterm/3;
                _EPOT_IND[kpt]+=eterm/3;
            }

            /* two body terms */            
            if(!Bond(ipt,jpt)) continue;

            if(group_01)
            {
                eterm2=(ALJ*rri6[j]-BLJ)*rri6[j]-rr[j]*DUDRc+(LJ_RC*DUDRc-Uc);
                f2.clear();   
                f2=rt[j]*(((12*ALJ*rri6[j]-6*BLJ)*rri6[j])*ri[j]+DUDRc);
            }
            else
            {
                au=aa*xpon2[j];
                pplm=pss*(bb*pow(ri[j],rho)-1.);
                eterm2=aa*(bb*pow(ri[j],rho)-1.)*xpon2[j];
                ff2=au*(rho*bb*pow(ri[j],rho1)+pplm*rdi2[j]);
                f2.clear();f2.addnv(ff2,rt[j]);
            }
            
            _F[ipt]-=f2;
            _F[jpt]+=f2;


            
            if(fixedatomenergypartition==0)
            {
                tote2+=eterm2;
            }
            else
            {
                if(!(fixed[ipt]&&fixed[jpt]))
                    tote2+=eterm2;
            }
//            tote2+=eterm2*(1.0-(fixed[ipt]+fixed[jpt])/2.0);
            /*tote2+=eterm2;*/
            _EPOT_IND[ipt]+=eterm2/2;
            _EPOT_IND[jpt]+=eterm2/2;
            
            if(fixedatomenergypartition==0)
            {
                _VIRIAL.addnvv(1.,f2,rg[j]);
            }
            else
            {
                if(!(fixed[ipt]||fixed[jpt]))
                    _VIRIAL.addnvv(1.,f2,rg[j]);
                else if(!(fixed[ipt]&&fixed[jpt]))
                    _VIRIAL.addnvv(0.5,f2,rg[j]);
                /*_VIRIAL.addnvv(1.,f2,rg[j]);*/
            }
        }
    }          
    _EPOT=tote2+tote3;

//    _EPOT_TIP=tipe2+tipe3;
//    _EPOT_SUB=sube2+sube3;
    
    tipforce.clear();    
    for(i=0;i<_NP;i++)
    {
        if((group[i]==0)||(group[i]==4))
            tipforce+=_F[i];
    }

    
    /* zero out the forces on fixed atoms */
//    for(i=0;i<_NP;i++)
//    {
//        if(fixed[i])
//            _F[i].clear();
//    }

    DUMP("SW complete");
}

void SWFrame::stillinger_weber_energyonly()
{
    ERROR("SWFrame: SWLJ energy only is not implemented!");
}

double SWFrame::stillinger_weber_energyonly(int iatom)
{             
    ERROR("SWFrame: SWLJ energy only is not implemented!");
    return 0;
}

void SWFrame::potential()
{
    stillinger_weber();
}

void SWFrame::potential_energyonly()
{
    stillinger_weber_energyonly();
}

double SWFrame::potential_energyonly(int iatom)
{
    return stillinger_weber_energyonly(iatom);
}

#ifdef _TEST

/* Main Program Begins */
class SWFrame sim;

/* The main program is defined here */
#include "main.cpp"

#endif//_TEST

