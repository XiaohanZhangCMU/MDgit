/*
  swsige.cpp
  by Wei Cai  caiwei@mit.edu
  Last Modified : Thu Apr 22 16:41:58 2010

  FUNCTION  :  MD simulation package of Si using Stillinger-Weber potential
*/

#include "swsige.h"

void SWFrame::initvars()
{
    /* Si */
#ifdef _SW_ORIG    
    /* original Si version PRB 31, 5262 (1985) */
       aa_si=15.27991323; bb_si=11.60319228; plam_si=45.51575;
       pgam_si=2.51412; acut_si=3.77118; pss_si=2.0951; 
#else   
    /* modified Si parameters PRB 46, 2250 (1992) */
       aa_si=16.31972277; bb_si=11.60319228; plam_si=48.61499998;
       pgam_si=2.51412;  acut_si=3.77118; pss_si=2.0951;
#endif    

    /* Ge version Ding and Anderson, PRB 34, 6987 (1986)
       (in Balamanne notation) */
       aa_ge=13.60564361; bb_ge=13.62639971; plam_ge=59.83;
       pgam_ge=2.6172; acut_ge=3.9258; pss_ge=2.181; 

       rho=4.0; 

       acut=acut_ge;
       _RLIST=acut*1.2;

       _SKIN=_RLIST-acut;
       rho1=rho+1; acutsq=acut*acut;

    strcpy(incnfile,"../si.cn");
    DUMP(HIG"SWFrame initvars"NOR);
    MDPARALLELFrame::initvars();
}

void SWFrame::stillinger_weber()
{
#ifndef NNM /* need to be taken care of later */
#define NNM 60
#endif
    /* Multi process function */
    double aa, bb, pgam, pss, acut_sp;
    double plam_i, pgam_ij, pgam_ik;
    int i,j,k,ipt,jpt,kpt;
    int neg;
    Vector3 sij,rij,f3i,f3j,f3k,f2;
    Vector3 rg[NNM],rt[NNM];
    double rr[NNM],ri[NNM],rdi[NNM],rdi2[NNM],xpon[NNM],xpon2[NNM];
    int list2[NNM];
    double r2ij,rrij;
    double ang,tm1,tm2,tm3,dhij,dhik,dhmu,tmij,tmik,tm2ij,tm2ik,
        tm3ik,tm3ij,eterm,tote3,au,pplm,ff2,eterm2,tote2;
    int n0, n1;

    DUMP("SW");
    
    refreshneighborlist();

    /* _EPOT, _F, _EPOT_IND, _VIRIAL all local */
    _EPOT=0;tote2=tote3=0;
    for(i=0;i<_NP;i++)
    { _F[i].clear(); _EPOT_IND[i]=0;}
    _VIRIAL.clear();

    /* divide jobs among CPU's
       TODO: put n0, n1 into MD.cpp
     */
    n0=0;
    n1=_NP;

    for(ipt=n0;ipt<n1;ipt++)
    {
        /* if fixed == -1, simply ignore this atom */
        if(fixed[ipt]==-1) continue;
        neg=0;
	
	if(species[ipt]==0)
	  plam_i=plam_si;
	else
	  plam_i=plam_ge;
	
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
            
	    if((species[ipt]==0)&&(species[jpt]==0))
            {
	      pgam=pgam_si; pss=pss_si; acut_sp = acut_si;
	    }
	    else if((species[ipt]==1)&&(species[jpt]==1))
	    {
	      pgam=pgam_ge; pss=pss_ge; acut_sp = acut_ge;
	    }
	    else 
	    {
	      pgam=(pgam_ge+pgam_si)/2; pss=(pss_ge+pss_si)/2;
	      acut_sp = (acut_si+acut_ge)/2;
	    }

            r2ij=rij.norm2();
            rrij=sqrt(r2ij);

            if(r2ij>acut_sp*acut_sp) continue;

            rg[neg]=rij;
            /*rt[neg]=rij/rrij;*/
            rt[neg]=rij; rt[neg]/=rrij;
            rr[neg]=rrij;
            ri[neg]=1./rrij;

            rdi[neg]=1./(rrij-acut_sp);
            rdi2[neg]=rdi[neg]*rdi[neg];

            if(fabs(pgam*rdi[neg])>30) xpon[neg]=0;/* avoid denormalization*/
            else xpon[neg]=exp(pgam*rdi[neg]);
            if(fabs(pss*rdi[neg])>30) xpon2[neg]=0;/* avoid denormalization*/
            else xpon2[neg]=exp(pss*rdi[neg]);
            list2[neg]=jpt;
            neg++;
        }
        /*
        for(j=0;j<neg;j++)
            INFO_Printf("xpon[%d]=%e rdi=%e pgam=%e\n",
                        j,xpon[j],rdi[j],pgam);
        */
        /* second inner loop */
        for(j=0;j<neg;j++)
        {
            /* three body loop */
            jpt=list2[j];

	    if((species[ipt]==0)&&(species[jpt]==0))
            {
	      pgam_ij=pgam_si;
	    }
	    else if((species[ipt]==1)&&(species[jpt]==1))
	    {
	      pgam_ij=pgam_ge;
	    }
	    else 
	    {
	      pgam_ij=(pgam_ge+pgam_si)/2;
	    }

            /* if fixed == -1, simply ignore this atom */
            if(fixed[jpt]==-1) continue;
            
            for(k=j+1;k<neg;k++)
            {
                kpt=list2[k];

		if((species[ipt]==0)&&(species[kpt]==0))
		{
		    pgam_ik=pgam_si;
		}
		else if((species[ipt]==1)&&(species[kpt]==1))
	        {
		    pgam_ik=pgam_ge;
	        }
   	        else 
	        {
	            pgam_ik=(pgam_ge+pgam_si)/2;
  	        }

                /* if fixed == -1, simply ignore this atom */
                if(fixed[kpt]==-1) continue;
                ang=rg[j]*rg[k]/(rr[j]*rr[k]);
                tm1=ang+1./3.;

                tm2=plam_i*tm1*tm1;
                tm3=xpon[j]*xpon[k];

                eterm=tm2*tm3;

                dhij=-1.*tm2*pgam_ij*rdi[j]*rdi[j]*tm3;


                dhik=-1.*tm2*pgam_ik*rdi[k]*rdi[k]*tm3;
                dhmu=2.*plam_i*tm3*tm1;
                    
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
//            if(ipt>jpt) continue;
            if(!Bond(ipt,jpt)) continue;
            if((species[ipt]==0)&&(species[jpt]==0))
            {
	      aa=aa_si; bb=bb_si; pss=pss_si; 
	    }
	    else if((species[ipt]==1)&&(species[jpt]==1))
	    {
	      aa=aa_ge; bb=bb_ge; pss=pss_ge;
	    }
	    else 
	    {
	      aa=(aa_ge+aa_si)/2; bb=(bb_ge+bb_si)/2; pss=(pss_ge+pss_si)/2; 
	    }
            au=aa*xpon2[j];
            pplm=pss*(bb*pow(ri[j],rho)-1.);
            ff2=au*(rho*bb*pow(ri[j],rho1)+pplm*rdi2[j]);
            /*f2=rt[j]*ff2;*/
            f2.clear();f2.addnv(ff2,rt[j]);

            _F[ipt]-=f2;
            _F[jpt]+=f2;

            eterm2=aa*(bb*pow(ri[j],rho)-1.)*xpon2[j];

            
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
    //_EPOT=tote2;
    /* zero out the forces on fixed atoms */
//    for(i=0;i<_NP;i++)
//    {
//        if(fixed[i])
//            _F[i].clear();
//    }
    
    DUMP("SW complete");
}



void SWFrame::potential()
{
    stillinger_weber();
}

//void SWFrame::potential_energyonly()
//{
//    stillinger_weber_energyonly();
//}

//double SWFrame::potential_energyonly(int iatom)
//{
//    return stillinger_weber_energyonly(iatom);
//}

#ifdef _TEST

/* Main Program Begins */
class SWFrame sim;

/* The main program is defined here */
#include "main.cpp"

#endif//_TEST

