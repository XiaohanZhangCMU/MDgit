/*
  tersoff88.h
  by Keonwook Kang kwkang75@yonsei.ac.kr
  Last Modified : 

  FUNCTION  :  
  Reference:   J. Tersoff, New empirical approach for the structure and 
               energy of covalent systems
               Phys. Rev. B, 37, 6991 (1988), Si(B)
               J. Tersoff, Empirical interatomic potential for silicon 
               with improved elastic properties
               Phys. Rev. B, 38, 9902 (1988), Si(C)

  Note:        Implemented by modifying tersoff.cpp
               In cut-off function fC, 
                   R_{ij} = R-D and S_{ij} = R+D or
                   R = (R_{ij} + S_{ij})/2, D = - (R_{ij} - S_{ij})/2

  Fixed box CGrelax working
  Free box (Parrinello-Raman) not working yet, need to debug Virial
*/

#include "tersoff88.h"

void Tersoff88Frame::initvars()
{
#if 1
    /* Fitting parameters for Si(B), PRB 38 9902 (1988) 
       Fitting parameters for Si,    PRB 37 6991 (1988)
	   C11 = 121, C12 = 86, C44 = 10 (GPa) poor reproduction */
    A_Si               = 3.2647e3;          // A (eV)
    B_Si               = 9.5373e1;          // B (eV)
    Lam_Si             = 3.2394;            // lambda_1 (1/\AA)
    Mu_Si              = 1.3258;            // lambda_2 (1/\AA)
    Lam3_Si            = 1.3258;            // lambda_3 (1/\AA) = lambda_2
    R_Si               = 2.8;               // R-D (\AA)
    S_Si               = 3.2;               // R+D (\AA)
    Alpha_Si           = 0.0;
    Beta_Si            = 0.33675;           // beta
    nTf_Si             = 22.956;            // n
    cTf_Si             = 4.8381;            // c
    dTf_Si             = 2.0417;            // d
    hTf_Si             = 0.0;               // h
    m_Si               = 3.0;               // default
#else 
    /* Fitting parameters for Si(C), PRB 38 9902 (1988) 
       Si(C) has better reproduction of elastic constants.
       C11 = 150, C12 = 80, C44 = 70 (GPa) */
    A_Si               = 1.8308e3;          // A (eV)
    B_Si               = 4.7118e2;          // B (eV)
    Lam_Si             = 2.4799;            // lambda_1 (1/\AA)
    Mu_Si              = 1.7322;            // lambda_2 (1/\AA)
    Lam3_Si            = 1.7322;            // lambda_3 (1/\AA) = lambda_2
    R_Si               = 2.7;               // R-D (\AA)
    S_Si               = 3.0;               // R+D (\AA)
    Alpha_Si           = 0.0;
    Beta_Si            = 1.0999e-6;         // beta
    nTf_Si             = 7.8734e-1;         // n
    cTf_Si             = 1.0039e5;          // c
    dTf_Si             = 1.6218e1;          // d
    hTf_Si             =-5.9826e-1;         // h
    m_Si               = 3.0;               // default
#endif
    
    /* Inter species parameter */   
    Chi_Si_Si          = 1.0;
    Chi2_Si_Si         = 1.0;    

    acut=S_Si*1.0;
    _RLIST=acut*1.1;
    _SKIN=_RLIST-acut;
    acutsq=acut*acut;

    strcpy(incnfile,"../si.cn");
    DUMP(HIG"TersoffFrame initvars"NOR);
    MDPARALLELFrame::initvars();
}

void Tersoff88Frame::tersoff()
{
#ifndef NNM /* need to be taken care of later */
#define NNM 60
#endif
    /* Multi process function */
    int i,j,k,ipt,jpt,kpt,neg;
    Vector3 sij,rij,f3i,fZ,fTheta;
    double cSqr, dSqr, cdSqr, invn2, fC, fA, fR, dfCdr, dfAdr, dfRdr;
    double Arg, temp, Z, CosTheta, dCosTheta, gTheta, BetaZN, BetaZN1, bij;
    double explr3, dexplr3dr, E, AlphaEN, AlphaEN1, aij;
    double r2ij, rrij, r;
    double V, dVdr, dVdE, dVdZ;
    double rr[NNM];
    //double Asave[NNM], Bsave[NNM], Musave[NNM], Lamsave[NNM];
    //double Rsave[NNM], Ssave[NNM];
    double fCsave[NNM], dfCdrsave[NNM], fAsave[NNM], dfAdrsave[NNM];
    double fRsave[NNM], dfRdrsave[NNM], dZdCosTheta[NNM];
    double dEdrij, dZdrij, dEdrik[NNM], dZdrik[NNM]; 
    Vector3 rg[NNM], rt[NNM];
    int list2[NNM];
    int n0, n1;

    DUMP("Tersoff");
    
    refreshneighborlist();

    /* _EPOT, _F, _EPOT_IND, _VIRIAL all local */
    _EPOT=0;
    for(i=0;i<_NP;i++)
    { _F[i].clear(); _EPOT_IND[i]=0.0;}
    _VIRIAL.clear();

    n0=0;
    n1=_NP;

    /* set potential parameters */
    A = A_Si; B = B_Si; Lam = Lam_Si; Mu = Mu_Si; Lam3 = Lam3_Si;
    R = R_Si; S = S_Si; Alpha = Alpha_Si; Beta = Beta_Si;
    nTf=nTf_Si; cTf=cTf_Si; dTf=dTf_Si; hTf=hTf_Si; m=m_Si;
    cSqr  = cTf*cTf; dSqr  = dTf*dTf; cdSqr = cSqr/dSqr;
    invn2 = -1.0/2.0/nTf;
    Chi = Chi_Si_Si; Chi2 = Chi2_Si_Si;

    for(ipt=n0;ipt<n1;ipt++)
    {
        if(fixed[ipt]==-1) continue; /* ignore atom if -1 */
      
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
            r = rrij;

            if(r>=S) continue;

            if(r<R)
            {
                fC = 1;
                dfCdr = 0;
            }
            else if (r<S)
            {
                Arg = M_PI*(r-.5*(S+R))/(S-R);
                fC = 0.5-0.5*sin(Arg);
                dfCdr = -0.5*M_PI*cos(Arg)/(S-R);
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

            dEdrij=0.0; dZdrij=0.0;
            Z = 0; /* bond order for attractive term */
            E = 0; /* bond order for repulsive term */
            /* three body loop */
            for(k=0;k<neg;k++)
            {
                if(k==j) continue; /* k and j are different neighbors */
                kpt=list2[k];
                /* if fixed == -1, simply ignore this atom */
                if(fixed[kpt]==-1) continue; /* ignore atom if -1 */

                /* cos(theta_ijk): angle between bond ij and bond ik */
                CosTheta=rt[j]*rt[k];
                /* d^2 + (h-cos(theta_ijk))^2 */
                dCosTheta = dSqr+(hTf-CosTheta)*(hTf-CosTheta);
                /* g(theta) */
                gTheta = 1+cdSqr-cSqr/dCosTheta;
                //explr3 = exp(pow(Lam3,3)*pow((rr[j]-rr[k]),3));
                explr3 = exp(pow(Lam3,m)*pow((rr[j]-rr[k]),m));
                E += fCsave[k]*explr3;
                Z += fCsave[k]*gTheta*explr3;
                /* d(explr3)/d(r_ik) */
                //dexplr3dr = explr3*pow(Lam3,3)*(-3.0)*pow((rr[j] - rr[k]),2);
                /* d(explr3)/d(r_ij) */
                //dexplr3dr = explr3*pow(Lam3,3)*(3.0)*pow((rr[j]-rr[k]),2);
                dexplr3dr = explr3*pow(Lam3,m)*m*pow((rr[j]-rr[k]),(m-1.0));
                dEdrij+=(fCsave[k]*dexplr3dr);
                dZdrij+=(fCsave[k]*gTheta*dexplr3dr);
                dEdrik[k]=(dfCdrsave[k]*explr3-fCsave[k]*dexplr3dr);
                dZdrik[k]=dEdrik[k]*gTheta;
                dZdCosTheta[k] = (fCsave[k]*explr3)*(-2.0*cSqr*(hTf-CosTheta)/(dCosTheta*dCosTheta));
            }
//            if(ipt==0)
//                INFO_Printf("atom %d-%d Z=%f\n",ipt,jpt,Z);
            
            fC = fCsave[j];
            dfCdr = dfCdrsave[j];
            fA = fAsave[j];
            dfAdr = dfAdrsave[j];
            fR = fRsave[j];
            dfRdr = dfRdrsave[j];
            
            AlphaEN = pow((Alpha*E),nTf);    // = (alpha * eta)^n
            BetaZN = pow((Beta*Z),nTf);     // = (beta * zeta)^n
            AlphaEN1 = pow(Alpha,nTf)*pow(E,(nTf-1.0)); // = alpha^n * eta^(n-1)
            BetaZN1 = pow(Beta,nTf)*pow(Z,(nTf-1.0));   // = beta^n * zeta^(n-1)
            if(neg==1) {AlphaEN1=0.0; BetaZN1=0.0;}
            //Chi = 1.0; Chi2 = 1.0;
            // The term aij can have species dependence e.g. Chi2
            aij = pow((1.0+AlphaEN),invn2);     
            bij = pow((1.0+BetaZN),invn2);
            V = 0.5*fC*(aij*fR+bij*fA);
            _EPOT += V;
            _EPOT_IND[ipt] += V;

            dVdE = 0.25*fC*fR*(-aij*AlphaEN1/(1+AlphaEN));
            dVdZ = 0.25*fC*fA*(-bij*BetaZN1/(1+BetaZN));
            dVdr = 0.5*dfCdr*(aij*fR+bij*fA) + 0.5*fC*(aij*dfRdr+bij*dfAdr) + dVdE*dEdrij + dVdZ*dZdrij;

            f3i=rt[j]*dVdr;
            _F[ipt]+=f3i;
            _F[jpt]-=f3i;
            _VIRIAL.addnvv(-1.,f3i,rg[j]);
            if (SURFTEN==1 && (curstep%savepropfreq)==1) AddnvvtoPtPn(_SR[jpt],f3i,rg[j],-1.0);
           
            /* three body loop */
            for(k=0;k<neg;k++)
            {
                if(k==j) continue; /* k and j are different neighbors */
                kpt=list2[k];
                if(fixed[kpt]==-1) continue; /* ignore atom if -1 */

                /* Forces due to derivatives with respect to Z, then rik, 
                   dV/d(zeta)*d(zeta)/dr_ik */
                temp = (dVdE*dEdrik[k] + dVdZ*dZdrik[k]);

                fZ = rt[k]*temp;
                _F[ipt]+=fZ;
                _F[kpt]-=fZ;
                _VIRIAL.addnvv(-1.,fZ,rg[k]);
		if (SURFTEN==1 && (curstep%savepropfreq)==1) AddnvvtoPtPn(_SR[kpt],fZ,rg[k],-1.0);
                /* cos(theta_ijk): angle between bond ij and bond ik */
                CosTheta=rt[j]*rt[k];
                
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
    

void Tersoff88Frame::potential()
{
    tersoff();
}

#ifdef _TEST

/* Main Program Begins */
class Tersoff88Frame sim;

/* The main program is defined here */
#include "main.cpp"

#endif//_TEST

