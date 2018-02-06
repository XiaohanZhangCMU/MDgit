/*
  tersoff.h
  by Wei Cai  caiwei@stanford.edu
  Last Modified : Thu Apr 22 16:42:07 2010

  FUNCTION  :  MD simulation package of Si/Ge/C using Tersoff potential
  Reference:   J. Tersoff, Modeling solid-state chemistry: Interatomic
               potentials for multicomponent systems,
               Phys. Rev. B, 39, 5566 (1989).

  Note:        Implemented with the help of the code from KJ Cho's group
               created: 08/26/2001 Ram (panchram@stanford.edu)
               revised: Byeongchan Lee (bclee@stanford.edu)

  Atom species:
               0 - Si
               1 - Ge
               2 - C

  Fixed box CGrelax working
  Free box (Parrinello-Raman) not working yet, need to debug Virial
*/

#include "tersoff.h"

void TersoffFrame::initvars()
{
    /* Fitting parameters for Si */
    A_Si               = 1.8308e3;
    B_Si               = 4.7118e2;
    Lam_Si             = 2.4799;
    Mu_Si              = 1.7322;
    R_Si               = 2.7;
    S_Si               = 3.0;
    Beta_Si            = 1.1e-6;
    nTf_Si             = 0.78734;
    cTf_Si             = 1.0039e5;
    dTf_Si             = 16.217;
    hTf_Si             = -0.59825;
    
    /* Fitting parameters for C */
    A_C                = 1.3936e3;
    B_C                = 3.467e2;
    Lam_C              = 3.4879;
    Mu_C               = 2.2119;
    R_C                = 1.8;    
    S_C                = 2.1;
    Beta_C             = 1.5724e-7; 
    nTf_C              = 0.72751; 
    cTf_C              = 3.8049e4;
    dTf_C              = 4.3484;  
    hTf_C              = -0.57058;
    
    /* Fitting parameters for Ge */ 
    A_Ge               = 1.769e3;   
    B_Ge               = 4.1923e2;
    Lam_Ge             = 2.4451;  
    Mu_Ge              = 1.7047;  
    R_Ge               = 2.8;     
    S_Ge               = 3.1;     
    Beta_Ge            = 9.0166e-7;
    nTf_Ge             = 0.75627;  
    cTf_Ge             = 1.0643e5; 
    dTf_Ge             = 15.652;   
    hTf_Ge             = -0.43884; 
    
    /* Inter species parameter */   
    Chi_C_Si           = 0.9776;    
    Chi_Si_Ge          = 1.00061;   
    Chi_C_Ge           = 1.0;       

    acut=S_Ge*1.0;
    _RLIST=acut*1.1;
    _SKIN=_RLIST-acut;
    acutsq=acut*acut;

    strcpy(incnfile,"../si.cn");
    DUMP(HIG"TersoffFrame initvars"NOR);
    MDPARALLELFrame::initvars();
}

void TersoffFrame::tersoff()
{
#ifndef NNM /* need to be taken care of later */
#define NNM 60
#endif
    /* Multi process function */
    int i,j,k,ipt,jpt,kpt,neg;
    Vector3 sij,rij,f3i,fZ,fTheta;
    double cSqr, dSqr, cdSqr, invn2, fC, fA, fR, dfCdr, dfAdr, dfRdr;
    double Arg, temp, Z, CosTheta, dCosTheta, gTheta, BetaZN, BetaZN1, bij;
    double r2ij, rrij, r;
    double V, dVdr, dVdZ;
    double rr[NNM], Asave[NNM], Bsave[NNM], Musave[NNM], Lamsave[NNM];
    double Rsave[NNM], Ssave[NNM];
    double fCsave[NNM], dfCdrsave[NNM], fAsave[NNM], dfAdrsave[NNM];
    double fRsave[NNM], dfRdrsave[NNM], dZdr[NNM], dZdCosTheta[NNM];
    Vector3 rg[NNM], rt[NNM];
    int list2[NNM];
    int n0, n1;

    DUMP("Tersoff");
    
    refreshneighborlist();

    /* _EPOT, _F, _EPOT_IND, _VIRIAL all local */
    _EPOT=0;
    for(i=0;i<_NP;i++)
    { _F[i].clear(); _EPOT_IND[i]=0;}
    _VIRIAL.clear();

#ifdef _TORSION_OR_BENDING
    if (_TORSIONSIM) _TORQUE = 0;
    if (_BENDSIM)    _BENDMOMENT = 0;
#endif

    n0=0;
    n1=_NP;

    for(ipt=n0;ipt<n1;ipt++)
    {
        if(fixed[ipt]==-1) continue; /* ignore atom if -1 */

        /* set potential parameters */
        switch (species[ipt]) {
        case 0: A_i=A_Si; B_i=B_Si; Lam_i=Lam_Si;
            Mu_i=Mu_Si; R_i=R_Si; S_i=S_Si; Beta=Beta_Si;
            nTf=nTf_Si; cTf=cTf_Si; dTf=dTf_Si; hTf=hTf_Si; break;
        case 1: A_i=A_Ge; B_i=B_Ge; Lam_i=Lam_Ge;
            Mu_i=Mu_Ge; R_i=R_Ge; S_i=S_Ge; Beta=Beta_Ge;
            nTf=nTf_Ge; cTf=cTf_Ge; dTf=dTf_Ge; hTf=hTf_Ge; break;
        case 2: A_i=A_C; B_i=B_C; Lam_i=Lam_C;
            Mu_i=Mu_C; R_i=R_C; S_i=S_C; Beta=Beta_C;
            nTf=nTf_C; cTf=cTf_C; dTf=dTf_C; hTf=hTf_C; break;
        }
        
        cSqr  = cTf*cTf;
        dSqr  = dTf*dTf;
        cdSqr = cSqr/dSqr;
        invn2 = -1.0/2.0/nTf;

        neg = 0;
        for(j=0;j<nn[ipt];j++)
        {
            jpt=nindex[ipt][j];
            if(fixed[jpt]==-1) continue; /* ignore atom if -1 */
            
            sij.subtract(_SR[jpt],_SR[ipt]);
            sij.subint();
            _H.multiply(sij,rij);
            
            r2ij=rij.norm2();
            if(r2ij>acutsq) continue;
            rrij=sqrt(r2ij);

            /* set potential parameters */
            switch (species[jpt]) {
            case 0: A_j=A_Si; B_j=B_Si; Lam_j=Lam_Si;
                Mu_j=Mu_Si; R_j=R_Si; S_j=S_Si; break;
            case 1: A_j=A_Ge; B_j=B_Ge; Lam_j=Lam_Ge;
                Mu_j=Mu_Ge; R_j=R_Ge; S_j=S_Ge; break;
            case 2: A_j=A_C; B_j=B_C; Lam_j=Lam_C;
                Mu_j=Mu_C; R_j=R_C; S_j=S_C; break;
            }
            
            /* bond ij potential parameters */
            A       = sqrt(A_i*A_j);
            B       = sqrt(B_i*B_j);
            Lam     = (Lam_i + Lam_j)*0.5;
            Mu      = (Mu_i + Mu_j)*0.5;
            R       = sqrt(R_i*R_j);
            S       = sqrt(S_i*S_j);

            r = rrij;

            if(r>=S) continue;

            if(r<R)
            {
                fC = 1;
                dfCdr = 0;
            }
            else if (r<S)
            {
                Arg = M_PI*(r-R)/(S-R);
                fC = 0.5+0.5*cos(Arg);
                dfCdr = -0.5*M_PI*sin(Arg)/(S-R);
            }
            else
            {
                fC = 0;
                dfCdr = 0;
            }

            fR = A*exp(-Lam*r);
            dfRdr = -Lam*fR;
            fA = -B*exp(-Mu*r);
            dfAdr = -Mu*fA;

            /* save all neighbors of i */
            rg[neg]=rij;                /* displacement vector */
            rt[neg]=rij; rt[neg]/=rrij; /* unit vector */
            rr[neg]=rrij;               /* distance */
            list2[neg]=jpt;             /* neighbor index */

            Asave[neg]=A_j;
            Bsave[neg]=B_j;
            Lamsave[neg]=Lam_j;
            Musave[neg]=Mu_j;
            Rsave[neg]=R_j;
            Ssave[neg]=S_j;
            
            fCsave[neg]=fC;
            dfCdrsave[neg]=dfCdr;
            fAsave[neg]=fA;
            dfAdrsave[neg]=dfAdr;
            fRsave[neg]=fR;
            dfRdrsave[neg]=dfRdr;
            
            neg++;
        }
//        if(ipt==0|ipt!=0)
//        {
//        INFO_Printf("atom %d neg=%d\n",ipt,neg);
//        for(j=0;j<neg;j++)
//            INFO_Printf(" neighbor j=%d fC=%f fR=%e fA=%e \n",j,fCsave[j],fRsave[j], fAsave[j]);
//        }
        
        /* going through the neighbors of atom i again */
        for(j=0;j<neg;j++)
        {
            jpt=list2[j];
            if(fixed[jpt]==-1) continue; /* ignore atom if -1 */

            Z = 0; /* bond order */
            /* three body loop */
            for(k=0;k<neg;k++)
            {
                if(k==j) continue; /* k and j are different neighbors */
                kpt=list2[k];
                /* if fixed == -1, simply ignore this atom */
                if(fixed[kpt]==-1) continue; /* ignore atom if -1 */

                /* cos(theta_ijk): angle between bond ij and bond ik */
                CosTheta=rt[j]*rt[k];

                dCosTheta = dSqr+(hTf-CosTheta)*(hTf-CosTheta);
                gTheta = 1+cdSqr-cSqr/dCosTheta;
                Z += fCsave[k]*gTheta;
                dZdr[k] = dfCdrsave[k]*gTheta;
                dZdCosTheta[k] = fCsave[k]*(-2.0*cSqr*(hTf-CosTheta)/(dCosTheta*dCosTheta));
            }
//            if(ipt==0)
//                INFO_Printf("atom %d-%d Z=%f\n",ipt,jpt,Z);
            
            fC = fCsave[j];
            dfCdr = dfCdrsave[j];
            fA = fAsave[j];
            dfAdr = dfAdrsave[j];
            fR = fRsave[j];
            dfRdr = dfRdrsave[j];
            
            BetaZN = pow((Beta*Z),nTf);
//            if(Z==0)
//            {
//                ERROR("Z=0 in Tersoff()");
//            }
            BetaZN1 = pow(Beta,nTf)*pow(Z,(nTf-1.0));

            Chi = 1.0;
            if ((species[ipt]==0)&&(species[jpt]==1))
                Chi = Chi_Si_Ge; /* Si Ge */
            else if ((species[ipt]==1)&&(species[jpt]==0))
                Chi = Chi_Si_Ge; /* Si Ge */
            else if ((species[ipt]==0)&&(species[jpt]==2))
                Chi = Chi_C_Si;  /* Si C */
            else if ((species[ipt]==2)&&(species[jpt]==0))
                Chi = Chi_C_Si;  /* Si C */
            
            bij = Chi*pow((1.0+BetaZN),invn2);
            V = 0.5*fC*(fR+bij*fA);
            _EPOT += V;
            _EPOT_IND[ipt] += V;

            dVdr = 0.5*dfCdr*(fR+bij*fA) + 0.5*fC*(dfRdr+bij*dfAdr);
            f3i=rt[j]*dVdr;
            _F[ipt]+=f3i;
            _F[jpt]-=f3i;
            _VIRIAL.addnvv(-1.,f3i,rg[j]);
            if (SURFTEN==1 && (curstep%savepropfreq)==1) AddnvvtoPtPn(_SR[jpt],f3i,rg[j],-1.0);

//            INFO_Printf(" bij = %e V = %e dVdr = %e f3i=(%e,%e,%e)\n",bij,V,dVdr,f3i.x,f3i.y,f3i.z);            
//            INFO_Printf(" bij = %e rt[j]=(%e,%e,%e)\n",bij,rt[j].x,rt[j].y,rt[j].z);            

            /* debug: changed from 0.5 to 0.25 */
            dVdZ = 0.25*fC*fA*(-bij*BetaZN1/(1+BetaZN));

            /* three body loop */
            for(k=0;k<neg;k++)
            {
                if(k==j) continue; /* k and j are different neighbors */
                kpt=list2[k];
                if(fixed[kpt]==-1) continue; /* ignore atom if -1 */

                /* cos(theta_ijk): angle between bond ij and bond ik */
                CosTheta=rt[j]*rt[k];

                /* Forces due to derivatives with respect to Z, then rik */
                temp = dVdZ*dZdr[k];

                fZ = rt[k]*temp;
                _F[ipt]+=fZ;
                _F[kpt]-=fZ;
                _VIRIAL.addnvv(-1.,fZ,rg[k]);
		if (SURFTEN==1 && (curstep%savepropfreq)==1) AddnvvtoPtPn(_SR[kpt],fZ,rg[k],-1.0);

                
                /* Forces due to derivatives with respect to cosTheta */
                temp = dVdZ*dZdCosTheta[k];

                fTheta = (rt[j]*CosTheta-rt[k])*(-temp/rr[j]);
                _F[ipt]+=fTheta;
                _F[jpt]-=fTheta;
                _VIRIAL.addnvv(-1.,fTheta,rg[j]);
		if (SURFTEN==1 && (curstep%savepropfreq)==1) AddnvvtoPtPn(_SR[jpt],fTheta,rg[j],-1.0);


                fTheta = (rt[k]*CosTheta-rt[j])*(-temp/rr[k]);
                _F[ipt]+=fTheta;
                _F[kpt]-=fTheta;
                _VIRIAL.addnvv(-1.,fTheta,rg[k]);
		if (SURFTEN==1 && (curstep%savepropfreq)==1) AddnvvtoPtPn(_SR[kpt],fTheta,rg[k],-1.0);

            }
        }          
    }
    /* zero out the forces on fixed atoms */
//    for(i=0;i<_NP;i++)
//    {
//        if(fixed[i])
//            _F[i].clear();
//    }

#if 0    
    DUMP("Tersoff Collect");
    /* _EPOT, _F, _EPOT_IND, _VIRIAL will be summed */
    MP_AccumulateTo(_shmEPOTptr,&(_EPOT),1);
    MP_AccumulateTo((double *)_shmF,(double *)_F,3*_NP);
    MP_AccumulateTo(_shmEPOT_IND,_EPOT_IND,_NP);
    MP_AccumulateTo((double *)((*_shmVIRIALptr)[0]),
                    (double *)(_VIRIAL[0]),9);
    
    StrMemCpy(&(_EPOT),_shmEPOTptr,sizeof(double));
    StrMemCpy((double *)_F,(double *)_shmF,3*_NP*sizeof(double));
    StrMemCpy(_EPOT_IND,_shmEPOT_IND,_NP*sizeof(double));
    StrMemCpy((double *)(_VIRIAL[0]),
              (double *)((*_shmVIRIALptr)[0]),9*sizeof(double));
    
    DUMP("Tersoff complete");
#endif
}

void TersoffFrame::potential()
{
    tersoff();
}

#ifdef _TEST

/* Main Program Begins */
class TersoffFrame sim;

/* The main program is defined here */
#include "main.cpp"

#endif//_TEST

