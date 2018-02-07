/*
  silica-bksmod.cpp
  by Wei Cai  caiwei@mit.edu
  Last Modified : Mon Jan  1 21:26:52 2007

  FUNCTION  :  Modified BKS model for Silica (SiO2)
               Si: species = 0
               O : species = 1
*/

#include "silica-bksmod.h"

void BKSFrame::initparser()
{
    MDPARALLELFrame::initparser();

    /* input */
    bindvar("ALPHA_BKS",&_ALPHA_BKS,DOUBLE);
}

int BKSFrame::exec(const char *name)
{
    if(MDPARALLELFrame::exec(name)==0) return 0;
    bindcommand(name,"initBKS",initBKS());
    return -1;
}

void BKSFrame::initvars()
{
    initBKS();

    _RLIST=BKS_RC+1.1;
    _SKIN=_RLIST-BKS_RC;
    MDPARALLELFrame::initvars();
}

void BKSFrame::initBKS()
{
    _A_00 = 0; _B_00 = 0; _C_00 = 0; /* Si-Si */
    _A_01 = 18003.7572; /* eV */     /* Si-O  */
    _B_01 = 4.87318;    /* A^-1 */
    _C_01 = 133.5381;   /* eV A^6 */
    _A_11 = 1388.7730;  /* eV */     /*  O-O  */
    _B_11 = 2.760;      /* A^-1 */
    _C_11 = 175.00;     /* eV A^6 */
    _Q_0  = 2.4;                     /* Si */
    _Q_1  =-1.2;                     /* O  */  

    BKS_RC = 8.0;       /* A */

    /* Modified BKS parameters */
    /* Coulomb truncation */
    _ALPHA_BKS = 0.192;     /* A^-1 */  

    /* short-range */
    _ALPHA_BKS = 0.5;     /* A^-1 */  

    /* short range repulsion */
    _EPS_00 = 1219.45 *1e3/P_E/N_A; /* kJ / mol -> eV */  /* Si-Si */
    _SIG_00 = 0.42;                 /* A */ 
    _EPS_01 = 1.083   *1e3/P_E/N_A; /* kJ / mol -> eV */  /* Si-O  */
    _SIG_01 = 1.31;                 /* A */ 
    _EPS_11 = 0.0344  *1e3/P_E/N_A; /* kJ / mol -> eV */  /*  O-O  */
    _SIG_11 = 2.20;                 /* A */ 

#if 0 /* debugging, turn off short-range correction */
    _EPS_00 = _EPS_01 = _EPS_11 = 0.0;
#endif 
}

void BKSFrame::calcprop()
{
    MDPARALLELFrame::calcprop();
}

void BKSFrame::bks_mod()
{
/* 
   U = q_i q_j erfc(alpha*r_ij) / r_ij +
       A_ij exp(-B_ij*r_ij) - C_ij / r_ij^6 +
       4 EPS_ij [ (SIG_ij/r_ij)^24 - (SIG_ij/r_ij)^6 ] 
*/
    /* !!! Parallel version not implemented !!! */
    int i,j,ipt,jpt;
    double U,r,r2,ri6;
    double Aij,Bij,Cij,SIGij,EPSij,QiQj;   
    double SIG2, SIG6, SIG_ri6, SIG_ri12, SIG_ri24, U_coulomb, tmp;
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
            if(r<=BKS_RC)
            {
                if((species[ipt]==0)&&(species[jpt]==0))
                {
                    Aij = _A_00; Bij = _B_00; Cij = _C_00; 
                    SIGij = _SIG_00;  EPSij = _EPS_00;
                    QiQj = _Q_0*_Q_0; 
                }
                else if((species[ipt]==1)&&(species[jpt]==1))
                {
                    Aij = _A_11; Bij = _B_11; Cij = _C_11; 
                    SIGij = _SIG_11;  EPSij = _EPS_11;
                    QiQj = _Q_1*_Q_1; 
                }
                else
                {
                    Aij = _A_01; Bij = _B_01; Cij = _C_01; 
                    SIGij = _SIG_01;  EPSij = _EPS_01;
                    QiQj = _Q_0*_Q_1; 
                }
               
                SIG2=SIGij*SIGij;
                SIG6=SIG2*SIG2*SIG2;                

                ri6=1./(r2*r2*r2);
                SIG_ri6 = SIG6 * ri6;
                SIG_ri12 = SIG_ri6 * SIG_ri6;
                SIG_ri24 = SIG_ri12 * SIG_ri12;

                U_coulomb = QiQj * erfc(_ALPHA_BKS * r) / r ;

                U = U_coulomb +
                    Aij * exp(-Bij * r) - Cij * ri6 +
                    4*EPSij * ( SIG_ri24 - SIG_ri6 );

                tmp = (U_coulomb + QiQj*2*_ALPHA_BKS/sqrt(M_PI) * 
                            exp(-_ALPHA_BKS*_ALPHA_BKS*r2)) / r2 +
                      Aij*Bij*exp(-Bij * r) / r - 6*Cij*ri6 / r2 +
                      4*EPSij * ( 24*SIG_ri24 - 6*SIG_ri6 ) / r2;

                fij = rij*tmp;  

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


void BKSFrame::potential()
{
    int i;
    bks_mod();

    /* Classical Ewald */
    if(Ewald_CE_or_PME==0){CE();}
    /* Paricle-Mesh Ewald */
    else{PME();}

    /* Sum up energy, force and stress */
    _EPOT += _EPOT_Ewald;
    _VIRIAL += _VIRIAL_Ewald;
    for(i=0;i<_NP;i++)
    {
        _F[i]+=_F_Ewald[i];
        _EPOT_IND[i]+=_EPOT_IND_Ewald[i];
    }

}

#ifdef _TEST

/* Main Program Begins */
class BKSFrame sim;

/* The main program is defined here */
#include "main.cpp"

#endif//_TEST

