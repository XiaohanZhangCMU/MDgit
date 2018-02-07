/*
  rebolj.cpp
  by Wei Cai  caiwei@stanford.edu
  Last Modified : Aug 9 2011

  FUNCTION  :  MD simulation package of C using REBO + LJ potential

  To Do:
   1. Possible error in gradient when relaxing with substrate (LJ interaction)

*/

#include "rebolj.h"

void REBOLJFrame::initvars()
{
    REBOFrame::initvars();

    int i, j;
    for(i=0;i<MAXSP;i++)
       for(j=i;j<MAXSP;j++)
       {
           _ALJ[i][j] = 0;
           _BLJ[i][j] = 0;
       }
}

void REBOLJFrame::initparser()
{
    REBOFrame::initparser();

    int i, j; char s[100];
    for(i=0;i<MAXSP;i++)
       for(j=i;j<MAXSP;j++)
       {
           sprintf(s,"C12_%d%d",i,j);
           bindvar(s,&(_ALJ[i][j]),DOUBLE);
           sprintf(s,"C6_%d%d",i,j);
           bindvar(s,&(_BLJ[i][j]),DOUBLE);
       }
    bindvar("LJ_Rcut",&LJ_RC,DOUBLE);

    bindvar("NUM_REBO_GROUPS",&_NUM_REBO_GROUPS,INT);

    bindvar("substrate_Z",  &_SUBSTRATE_Z,DOUBLE);
    bindvar("substrate_REP",&_SUBSTRATE_REP,DOUBLE);
    bindvar("substrate_ATR",&_SUBSTRATE_ATR,DOUBLE);
}

int REBOLJFrame::exec(const char *name)
{
    if(REBOFrame::exec(name)==0) return 0;
    bindcommand(name,"initLJ",initLJ());
    return -1;
}

void REBOLJFrame::initLJ()
{
    /* should be called after read_AIREBO, which sets acut */
    if (LJ_RC > (_RLIST-_SKIN) ) {
      _RLIST = LJ_RC * 1.1;
      _SKIN = _RLIST - LJ_RC;
    }

    int i, j;
    for(i=0;i<MAXSP;i++)
       for(j=i;j<MAXSP;j++)
       {
         /* always use LJ interactions */
         ALJ[i][j]=_ALJ[i][j]*LJ_ENERGY*POW12(LJ_LENGTH);
         BLJ[i][j]=_BLJ[i][j]*LJ_ENERGY*POW6(LJ_LENGTH);
         Uc[i][j]=(ALJ[i][j]*POW6(1/LJ_RC)-BLJ[i][j])*POW6(1/LJ_RC);
         DUDRc[i][j]=-(12*ALJ[i][j]/POW6(LJ_RC)-6*BLJ[i][j])/POW6(LJ_RC)/LJ_RC;
       }
    /* symmetrize LJ parameter matrix */
    for(i=1;i<MAXSP;i++)
       for(j=0;j<i;j++)
       {
	   _ALJ[i][j] = _ALJ[j][i];
           _BLJ[i][j] = _BLJ[j][i];
            ALJ[i][j] =  ALJ[j][i];
            BLJ[i][j] =  BLJ[j][i];
             Uc[i][j] =   Uc[j][i];
          DUDRc[i][j] =DUDRc[j][i];
       }

}

void REBOLJFrame::lennard_jones()
{
   int ipt, j, jpt;
   Vector3 si, ri, fi;
   Vector3 sij, rij, fij;
   double ALJ_local,BLJ_local,Uc_local,DUDRc_local;   
   double U,r,r2,ri6, R0,U0,K0,ri02,ri06;
  
   /* called immediately after rebo, no need to set _EPOT, _F to zero */

   /* LJ interaction */
    for(ipt=0;ipt<_NP;ipt++)
    {
        for(j=0;j<nn[ipt];j++)
        {
            jpt=nindex[ipt][j];
            if(ipt>jpt) continue;

            if (_NUM_REBO_GROUPS >= 1) 
            { /* only include LJ interactions between different groups */
              if (group[ipt] == group[jpt]) continue;
            }

            sij=_SR[jpt]-_SR[ipt];
            sij.subint();
            rij=_H*sij;
            r2=rij.norm2();
            r=sqrt(r2);

            if(r<=LJ_RC)
            {
                ALJ_local =   ALJ[species[ipt]][species[jpt]];
                BLJ_local =   BLJ[species[ipt]][species[jpt]];
                Uc_local =    Uc[species[ipt]][species[jpt]];
		DUDRc_local = DUDRc[species[ipt]][species[jpt]];

                R0 = 0.8;  /* short range cut-off */
                if(r<R0)
                {
                  ri02  = 1./(R0*R0);
                  ri06 = ri02*ri02*ri02;
                  U0 = (7*ALJ_local*ri06 - 4*BLJ_local)*ri06; 
                  K0 = (6*ALJ_local*ri06 - 3*BLJ_local)*ri06*ri02;
                  U = U0 - K0*r2 - r*DUDRc_local + (LJ_RC*DUDRc_local-Uc_local);
                  fij = rij*(2*K0 + DUDRc_local/r);
                }
                else
		{
                    /* always use LJ interactions */
		    ri6=1./(r2*r2*r2);
		    U=(ALJ_local*ri6-BLJ_local)*ri6
                      -r*DUDRc_local+(LJ_RC*DUDRc_local-Uc_local);
                    fij=rij*((12*ALJ_local*ri6-6*BLJ_local)*ri6/r2+DUDRc_local/r);
		}

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

void REBOLJFrame::substrate_potential()
{
   int ipt;
   Vector3 si, ri, fi;
   double z, esub, fsub;

   /* interaction with smooth substrate */
   if ((_SUBSTRATE_REP != 0) || (_SUBSTRATE_ATR != 0))
   {
    SHtoR(); 
    for (ipt=0;ipt<_NP;ipt++)
    {
        si=_SR[ipt]; 
        ri=_H*si;
        z=ri.z;
 
        if (z < _SUBSTRATE_Z) {
           esub = 0.5*_SUBSTRATE_REP*(z-_SUBSTRATE_Z)*(z-_SUBSTRATE_Z);
           fsub = _SUBSTRATE_REP*(_SUBSTRATE_Z-z); 
        } else {
           esub = _SUBSTRATE_ATR*(1.0-exp(-0.5*(z-_SUBSTRATE_Z)*(z-_SUBSTRATE_Z)));
           fsub = _SUBSTRATE_ATR*((_SUBSTRATE_Z-z)*exp(-0.5*(z-_SUBSTRATE_Z)*(z-_SUBSTRATE_Z))); 
        }
        _F[ipt].z += fsub;
        _EPOT_IND[ipt] += esub;
        _EPOT += esub;
   }
  }
}

void REBOLJFrame::potential()
{
    int igroup, ipt;

    _EPOT=0; _VIRIAL.clear();
    for(ipt=0;ipt<_NP;ipt++) { _F[ipt].clear(); _EPOT_IND[ipt]=0; }

    if (_NUM_REBO_GROUPS <= 1) 
    {
      REBOFrame::rebo();
      lennard_jones();
    }
    else
    {
      refreshneighborlist();

      for(igroup=0;igroup<_NUM_REBO_GROUPS;igroup++)
      {
         /* assuming all atoms have fixed[ipt] = 0 initially */
         for(ipt=0;ipt<_NP;ipt++)
         {
            if(group[ipt]==igroup) fixed[ipt] = 0;  /* mark atom as free */
            else fixed[ipt] = -1;                   /* mark atom as removed */
         }
         REBOFrame::rebo();
      }
      for(ipt=0;ipt<_NP;ipt++) fixed[ipt] = 0;  /* mark atom as free */
      lennard_jones();
    }

    substrate_potential();
}


#ifdef _TEST

/* Main Program Begins */
class REBOLJFrame sim;

/* The main program is defined here */
#include "main.cpp"

#endif//_TEST

