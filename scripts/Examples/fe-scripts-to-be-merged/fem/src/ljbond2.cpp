/*
  ljbond.cpp
  by Wei Cai  caiwei@mit.edu
  Last Modified : Fri Mar  2 18:51:39 2007

  FUNCTION  :  MD simulation package of Lennard-Jones (LJ) potential
               LJ parameters depend on species
               Atoms can also form covalent bond with neighbors

  Test cases: scripts/work/ar/lj2-md.script, ljbond.tcl
*/

#include "ljbond2.h"

void LJBONDFrame::initparser()
{
    int i, j;
    char s[100];
    MDPARALLELFrame::initparser();

    /* input */
    for(i=0;i<MAXSP;i++)
       for(j=i;j<MAXSP;j++)
       {
           sprintf(s,"C12_%d%d",i,j);
           bindvar(s,&(_ALJ[i][j]),DOUBLE);
           sprintf(s,"C6_%d%d",i,j);
           bindvar(s,&(_BLJ[i][j]),DOUBLE);
       }

    bindvar("BOND_R0",&BOND_R0,DOUBLE);
    bindvar("BOND_K", &BOND_K, DOUBLE);

    bindvar("ATTRAC2_Y0",&ATTRAC2_Y0,DOUBLE);
    bindvar("ATTRAC2_K",&ATTRAC2_K,DOUBLE);
    bindvar("WALL3_X0MAX",&WALL3_X0MAX,DOUBLE);
    bindvar("WALL3_X0",&WALL3_X0,DOUBLE);
    bindvar("WALL3_K",&WALL3_K,DOUBLE);
    bindvar("LMAX",&LMAX,DOUBLE);
    bindvar("XMAX",&XMAX,DOUBLE);
    bindvar("YMAX",&YMAX,DOUBLE);
    bindvar("ZMAX",&ZMAX,DOUBLE);
    bindvar("LAMBDA",&LAMBDA,DOUBLE);
    bindvar("TNPT_WALL3",&TNPT_WALL3,DOUBLE);

    bindvar("dUdLAM_ATT2",&dUdLAM_ATT2,DOUBLE);
    bindvar("dUdLAM_WALL3",&dUdLAM_WALL3,DOUBLE);
    bindvar("dUdLAM_L",&dUdLAM_L,DOUBLE);
    bindvar("dUdLAM_ALKANE", &dUdLAM_ALKANE, DOUBLE);
    bindvar("dUdLAM",&dUdLAM,DOUBLE);
    bindvar("dUdLAM_X",&dUdLAM_X,DOUBLE);
    bindvar("dUdLAM_Y",&dUdLAM_Y,DOUBLE);
    bindvar("dUdLAM_Z",&dUdLAM_Z,DOUBLE);
    bindvar("dUdLAM_XYZ",&dUdLAM_XYZ,DOUBLE);
    bindvar("EWALL3",&EWALL3,DOUBLE);

    bindvar("Rcut",&LJ_RC,DOUBLE);
}

int LJBONDFrame::exec(const char *name)
{
    if(MDPARALLELFrame::exec(name)==0) return 0;
    bindcommand(name,"initLJ",initLJ());
    bindcommand(name,"makelipids",makelipids());
    bindcommand(name,"linklipids",linklipids());
    bindcommand(name,"makealkanes",makealkanes());
    bindcommand(name,"linkalkanes",linkalkanes());
    return -1;
}

void LJBONDFrame::Alloc()
{
    int size;
    MDPARALLELFrame::Alloc();

    size = _NP*allocmultiple;
    Realloc(num_bonds, int,size);
    Realloc(bond_index,int,size*MAXNUMBOND);

    memset(num_bonds, 0,sizeof(int)*size);
    memset(bond_index,0,sizeof(int)*size*MAXNUMBOND);
}

void LJBONDFrame::Alloc1()
{
    int size;

    size = _NP*allocmultiple;
    Realloc(num_bonds, int,size);
    Realloc(bond_index,int,size*MAXNUMBOND);

    memset(num_bonds, 0,sizeof(int)*size);
    memset(bond_index,0,sizeof(int)*size*MAXNUMBOND);
}

void LJBONDFrame::initvars()
{
    int i, j;
    for(i=0;i<MAXSP;i++)
       for(j=i;j<MAXSP;j++)
       {
           _ALJ[i][j] = 4;
           _BLJ[i][j] = 4;
       }

    BOND_R0 = LJ_LENGTH;
    BOND_K  = 2;
    LJ_RC = 2.37343077641 * LJ_LENGTH; /* initial value, 4th neighbor */
    _RLIST=LJ_RC+1.1;
    _SKIN=_RLIST-LJ_RC;
    MDPARALLELFrame::initvars();
}

void LJBONDFrame::initLJ()
{
    int i, j;

    for(i=0;i<MAXSP;i++)
       for(j=i;j<MAXSP;j++)
       {
	 if ((i==0 && j==1) || (i==0 && j==3) || (i==1 && j==2) || (i==2 && j==3)) /*check for w-t and t-h interaction, use Usc and not ULJ*/
	   {
	     ALJ[i][j]=_ALJ[i][j]*LJ_ENERGY*CUBE(CUBE((1.05*LJ_LENGTH)));
	     BLJ[i][j]=_BLJ[i][j]*LJ_ENERGY*POW6(1.05*LJ_LENGTH);
	     Uc[i][j]=ALJ[i][j]*CUBE(CUBE(1/LJ_RC))-BLJ[i][j]*POW6(1/LJ_RC); /*Uc[i][j]=(ALJ[i][j]*POW12(1/LJ_RC)-BLJ[i][j])*POW6(1/LJ_RC)*/
	     DUDRc[i][j]=-9*ALJ[i][j]/CUBE(CUBE(LJ_RC))/LJ_RC-6*BLJ[i][j]/POW6(LJ_RC)/LJ_RC;
	   }
	 else
	   {
	     ALJ[i][j]=_ALJ[i][j]*LJ_ENERGY*POW12(LJ_LENGTH);
	     BLJ[i][j]=_BLJ[i][j]*LJ_ENERGY*POW6(LJ_LENGTH);
	     Uc[i][j]=(ALJ[i][j]*POW6(1/LJ_RC)-BLJ[i][j])*POW6(1/LJ_RC);
	     DUDRc[i][j]=-(12*ALJ[i][j]/POW6(LJ_RC)-6*BLJ[i][j])/POW6(LJ_RC)/LJ_RC;
	   }
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

    _RLIST=LJ_RC+1.1;
    _SKIN=_RLIST-LJ_RC;

}

void LJBONDFrame::lennard_jones_bond()
{
    int i,j,ipt,jpt,k, bonded;
    double U,r,r2,ri6, R0,U0,K0,ri02,ri06;
    double ALJ_local,BLJ_local,Uc_local,DUDRc_local;   
    Vector3 sij, rij, fij;
    
    DUMP(HIG"Lennard Jones"NOR);
        
    refreshneighborlist();
    dUdLAM_ALKANE=0;
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

            /* test if atoms ipt and jpt have a covalent bond */
            bonded = 0;
            for(k=0;k<num_bonds[ipt];k++)
            {
                if(bond_index[ipt*MAXNUMBOND+k] == jpt)
                {
                    bonded = 1;
                    break;
                }
            }
            /* exclude covalently bonded pairs from LJ interaction */
            if(bonded) continue;

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
		  if ((species[ipt]==0 && species[jpt]==1) ||(species[ipt]==1 && species[jpt]==0) ||(species[ipt]==1 && species[jpt]==2) ||(species[ipt]==2 && species[jpt]==1)||(species[ipt]==0 && species[jpt]==3)||(species[ipt]==3 && species[jpt]==0)||(species[ipt]==2 && species[jpt]==3)||(species[ipt]==3 && species[jpt]==2))
		  {
		    ri6=1./(r2*r2*r2);
		    U=ALJ_local*ri6/CUBE(r)-BLJ_local*ri6-r*DUDRc_local+(LJ_RC*DUDRc_local-Uc_local);
                    fij=rij*(9*ALJ_local*ri6*ri6*r-6*BLJ_local*ri6/r2+DUDRc_local/r);
		  }
		  else
		  {
		    ri6=1./(r2*r2*r2);
		    U=(ALJ_local*ri6-BLJ_local)*ri6-r*DUDRc_local+(LJ_RC*DUDRc_local-Uc_local);
                    fij=rij*((12*ALJ_local*ri6-6*BLJ_local)*ri6/r2+DUDRc_local/r);
		  }
                }

                _F[ipt]-=fij;
                _F[jpt]+=fij;
                if (((species[ipt]==3)&&(species[jpt]!=3)) || ((species[jpt]==3)&&(species[ipt]!=3)))
                 {
                    dUdLAM_ALKANE+=U;
                 } 
                _EPOT_IND[ipt]+=U*0.5;
                _EPOT_IND[jpt]+=U*0.5;
                _EPOT+=U;
                _VIRIAL.addnvv(1.,fij,rij);
            }
        }
    }

    /* covalently bonded interaction */
  if(BOND_K != 0)
  {
    for(ipt=0;ipt<_NP;ipt++)
    {
       //INFO_Printf("ljbond: ipt=%d num_bonds = %d\n",ipt,num_bonds[ipt]);
       for(j=0;j<num_bonds[ipt];j++)
       {
          jpt = bond_index[ipt*MAXNUMBOND + j];
          //INFO_Printf("ljbond: ipt=%d jpt=%d\n",ipt,jpt);

          if(ipt>jpt) continue;

          sij=_SR[jpt]-_SR[ipt];
          sij.subint();
          rij=_H*sij;
          r2=rij.norm2();
          r=sqrt(r2);
                
          U = BOND_K * SQR(r-BOND_R0);
          fij = rij * (-2*BOND_K*(1-BOND_R0/r)); /* bug fix from old ljbond.cpp */

          _F[ipt]-=fij;
          _F[jpt]+=fij;
          _EPOT_IND[ipt]+=U*0.5;
          _EPOT_IND[jpt]+=U*0.5;
          _EPOT+=U;
          _VIRIAL.addnvv(1.,fij,rij);
       }
       //INFO_Printf("ljbond: F[%d] = (%e %e %e)\n",ipt,_F[ipt].x,_F[ipt].y,_F[ipt].z);       
    }
  }

  /* add external potential */
  Vector3 si, ri, fi;
  double x, y, fattrac2, eattrac2, fwall3, ewall3;
  double du_att2,du_wall3;
  dUdLAM_ATT2=0;
  dUdLAM_WALL3=0;
  EWALL3=0;
  dUdLAM_L=0;
  dUdLAM=0;
  TNPT_WALL3=0;
   
  if ((ATTRAC2_Y0 > 0 ) && (ATTRAC2_K > 0 )) {
   for (ipt=0;ipt<_NP;ipt++)
   {
       if (species[ipt]!=1) continue; 
       si=_SR[ipt]; si.subint();
       ri=_H*si;
       y=ri.y;

    if (fabs(y) < ATTRAC2_Y0)
    {
       eattrac2=-1.0*ATTRAC2_K*pow( (1.0-pow(1.0*(y/ATTRAC2_Y0),2.0)),2.0);
       fattrac2=-1.0*ATTRAC2_K*(1.0-pow(1.0* (y/ATTRAC2_Y0),2.0))*4.0*y/pow(1.0* (ATTRAC2_Y0),2.0);
       _F[ipt].y+=fattrac2;
       _EPOT_IND[ipt]+=eattrac2;
       _EPOT+=eattrac2;

       _VIRIAL[1][1] += fattrac2*y;

       if (LAMBDA>0)
       du_att2=(eattrac2)/LAMBDA;
       else
       du_att2=0;
       dUdLAM_ATT2=dUdLAM_ATT2+du_att2;
    }

   }
  }

  if ((WALL3_X0 > 0) && (WALL3_K >0)) {
   for (ipt=0;ipt<_NP;ipt++)
   {
       /* multiply a factor depending on group ID */
       double fact;
       if (group[ipt]!=0) fact = pow(0.3, group[ipt]);
       else               fact = 1.0;
       si=_SR[ipt]; si.subint();
       ri=_H*si;
       x=ri.x;
  
    if (fabs(x) < WALL3_X0) 
    {
       ewall3=WALL3_K*pow((pow(WALL3_X0,2)-pow(x,2)),2)/pow(WALL3_X0,2);
       fwall3=WALL3_K*4*(pow(WALL3_X0,2)-pow(x,2))*x/pow(WALL3_X0,2);

       /* multiply a factor depending on group ID */
       ewall3*=fact;
       fwall3*=fact;

       _F[ipt].x+=fwall3;
       _EPOT_IND[ipt]+=ewall3;
       _EPOT+=ewall3;
       
       if (species[ipt]!=3)  
       {
	    _VIRIAL[0][0] += fwall3*x;
       }
       ++TNPT_WALL3;

       if ((LAMBDA>0) && (species[ipt]!=3)) 
       du_wall3=fwall3*pow(WALL3_X0,2)/(x*LAMBDA)-2*ewall3/LAMBDA;
       else
       du_wall3=0;

       dUdLAM_WALL3=dUdLAM_WALL3+du_wall3;
       EWALL3=EWALL3+ewall3;
    }

   }
  }

  dUdLAM=dUdLAM_ATT2+dUdLAM_WALL3;

}

void LJBONDFrame::linklipids()
{
    int nsolvent, nlipid, nalkane, chainlen, alkanelen, headid;
    int ipt, i, j;
    double bond_len;
    Vector3 dr, dr1, ds;
    Matrix33 hinv;

    nsolvent = (int) input[0];
    nlipid = (int) input[1];
    chainlen = (int) input[2];
    headid = (int) input[3];
    bond_len = input[4];
    nalkane = (int) input[5];
    alkanelen = (int) input[6];

    for(ipt=0;ipt<_NP;ipt++) group[ipt] = 0;

    ipt = 0;
    for(i=0;i<nlipid;i++)
    {
        for(j=0;j<chainlen;j++)
        {
           if(j==0)
           {
               /* only one bond */
               num_bonds[ipt] = 1;
               bond_index[ipt*MAXNUMBOND + 0] = ipt+1;
               
           } else if (j==(chainlen-1))
           {   /* only one bond */
               num_bonds[ipt]=1;
               bond_index[ipt*MAXNUMBOND + 0] = ipt-1;
           } else
           {   /* two bonds */
               num_bonds[ipt]=2;
               bond_index[ipt*MAXNUMBOND + 0] = ipt-1;
               bond_index[ipt*MAXNUMBOND + 1] = ipt+1;
           }
           group[ipt] = j;
           ipt ++;
        }
    }
}

void LJBONDFrame::linkalkanes()
{
    int nsolvent, nlipid, nalkane, chainlen, alkanelen, headid;
    int ipt, i, j;
    double bond_len;
    Vector3 dr, dr1, ds;
    Matrix33 hinv;

    nsolvent = (int) input[0];
    nlipid = (int) input[1];
    chainlen = (int) input[2];
    headid = (int) input[3];
    bond_len = input[4];
    nalkane = (int) input[5];
    alkanelen = (int) input[6];
    ipt = nlipid*chainlen+nsolvent;
   
    for(i=0;i<nalkane;i++)
    {
        for(j=0;j<alkanelen;j++)
        {
           if(j==0)
           {
               /* only one bond */
               num_bonds[ipt] = 1;
               bond_index[ipt*MAXNUMBOND + 0] = ipt+1;
               
           } else if (j==(alkanelen-1))
           {   /* only one bond */
               num_bonds[ipt]=1;
               bond_index[ipt*MAXNUMBOND + 0] = ipt-1;
           } else
           {   /* two bonds */
               num_bonds[ipt]=2;
               bond_index[ipt*MAXNUMBOND + 0] = ipt-1;
               bond_index[ipt*MAXNUMBOND + 1] = ipt+1;
           }
           ipt ++;
        }
    }
}


void LJBONDFrame::makelipids()
{
    int nsolvent, nlipid, nalkane, chainlen, alkanelen, headid;
    int ipt, i, j;
    double bond_len,d_by_two,offset;
    Vector3 dr, dr1, ds;
    Matrix33 hinv;

    nsolvent = (int) input[0];
    nlipid = (int) input[1];
    chainlen = (int) input[2];
    headid = (int) input[3];
    bond_len = input[4];
    nalkane = (int) input[5];
    alkanelen = (int) input[6];
    d_by_two=0.04;
    offset=0.2;
    
    _H.clear();
    _H[0][0] = input[7];
    _H[1][1] = input[8];
    _H[2][2] = input[9];

    _NP = nlipid*chainlen + nsolvent;
    Alloc();

    hinv = _H.inv();
    ipt = 0;
    for(i=0;i<nlipid;i++)
    {
        for(j=0;j<chainlen;j++)
        {
           if(j==0)
           {
               /* randomly position the first atom */
               _SR[ipt].set(drand48()-0.5,drand48()-0.5,drand48()-0.5);

               /* randomly choose a growth direction */
               dr.set(drand48()-0.5,drand48()-0.5,drand48()-0.5);
               dr/=dr.norm();  dr*=bond_len; ds = hinv*dr;

           } else {
               /* add a new atom in the randomly chosen direction */
               _SR[ipt] = _SR[ipt-1] + ds;
           }

           if(j==headid)
           {
               species[ipt] = 2; /* head particle (hydrophilic) */

               /* choose a slightly different growth direction */
               dr1.set(drand48()-0.5,drand48()-0.5,drand48()-0.5);
               dr1/=dr1.norm(); dr1*=(bond_len*0.5); dr*=-1; dr+=dr1;
               dr/=dr.norm();  dr*=bond_len; ds = hinv*dr;

           } else {
               species[ipt] = 1; /* tail particle (hydrophobic) */
           }

           ipt ++;
        }
    }

    linklipids();

    for(i=0;i<nsolvent;i++)
    {
        /* randomly position the solvent particle */
        _SR[ipt].set(drand48()-0.5,drand48()-0.5,drand48()-0.5);
        species[ipt] = 0;     /* solvent particle (water) */
        ipt ++;
    }    }
    
void LJBONDFrame::makealkanes()
{
    int nsolvent, nlipid, nalkane, chainlen, alkanelen, headid;
    int ipt, i, j;
    double bond_len,d_by_two,offset;
    Vector3 dr, dr1, ds;
    Matrix33 hinv;

    nsolvent = (int) input[0];
    nlipid = (int) input[1];
    chainlen = (int) input[2];
    headid = (int) input[3];
    bond_len = input[4];
    nalkane = (int) input[5];
    alkanelen = (int) input[6];
    d_by_two=0.032;

    offset=0.2;

    INFO_Printf("old NP=%d",_NP);
    _NP = nlipid*chainlen + nalkane*alkanelen + nsolvent;
    INFO_Printf("old NP=%d",_NP);
    Alloc1();


    hinv = _H.inv();
    ipt = nlipid*chainlen+nsolvent;
  
    for(i=0;i<(nalkane/4);i++)
    {
        for(j=0;j<alkanelen;j++)
        {
           if(j==0)
           {
               /* randomly position the first atom */
               _SR[ipt].set(0.158,d_by_two, 0.5-0.0973*i);
/*(WALL3_X0/_H[1][1])+3.33/_H[1][1]*(alkanelen-1),d_by_two, 0.5-0.1*i);

               /* randomly choose a growth direction */
               dr.set(1,0,0);
               dr/=dr.norm();  dr*=bond_len; ds = hinv*dr;

           } else {
               /* add a new atom in the randomly chosen direction */
               _SR[ipt] = _SR[ipt-1] + ds;
           }

           if(j==0)
           {
               species[ipt] = 3; /* head particle (hydrophilic) */
	       fixed[ipt]=1; /*fix head particle*/

           } else {
               species[ipt] = 3; /* tail particle (hydrophobic) */
           }

           ipt ++;
        }
    }

    for(i=0;i<(nalkane/4);i++)
    {
        for(j=0;j<alkanelen;j++)
        {
           if(j==0)
           {
               /* randomly position the first atom */
               _SR[ipt].set(0.158,-d_by_two, 0.5-0.0973*i);
/*(WALL3_X0/_H[1][1])+3.33/_H[1][1]*(alkanelen-1),-d_by_two, 0.5-0.1*i);

               /* randomly choose a growth direction */
               dr.set(1,0,0);
               dr/=dr.norm();  dr*=bond_len; ds = hinv*dr;

           } else {
               /* add a new atom in the randomly chosen direction */
               _SR[ipt] = _SR[ipt-1] + ds;
           }

           if(j==0)
           {
               species[ipt] = 3; /* head particle (hydrophilic) */
	       fixed[ipt]=1; /*fix head particle*/

           } else {
               species[ipt] = 3; /* tail particle (hydrophobic) */
           }

           ipt ++;
        }
    }

    for(i=0;i<(nalkane/4);i++)
    {
        for(j=0;j<alkanelen;j++)
        {
           if(j==0)
           {
               /* randomly position the first atom */
               _SR[ipt].set(-0.158,d_by_two, 0.45-0.0973*i);
/*-(WALL3_X0/_H[1][1])+3.33/_H[1][1]*(alkanelen-1),d_by_two, 0.5-0.1*i);

               /* randomly choose a growth direction */
               dr.set(-1,0,0);
               dr/=dr.norm();  dr*=bond_len; ds = hinv*dr;

           } else {
               /* add a new atom in the randomly chosen direction */
               _SR[ipt] = _SR[ipt-1] + ds;
           }

           if(j==0)
           {
               species[ipt] = 3; /* head particle (hydrophilic) */
	       fixed[ipt]=1; /*fix head particle*/

           } else {
               species[ipt] = 3; /* tail particle (hydrophobic) */
           }

           ipt ++;
        }
    }

    for(i=0;i<(nalkane/4);i++)
    {
        for(j=0;j<alkanelen;j++)
        {
           if(j==0)
           {
               /* randomly position the first atom */
               _SR[ipt].set(-0.158,-d_by_two, 0.45-0.0973*i);
/*-(WALL3_X0/_H[1][1])+3.33/_H[1][1]*(alkanelen-1) ,-d_by_two, 0.5-0.1*i);

               /* randomly choose a growth direction */
               dr.set(-1,0,0);
               dr/=dr.norm();  dr*=bond_len; ds = hinv*dr;

           } else {
               /* add a new atom in the randomly chosen direction */
               _SR[ipt] = _SR[ipt-1] + ds;
           }

           if(j==0)
           {
               species[ipt] = 3; /* head particle (hydrophilic) */
	       fixed[ipt]=1; /*fix head particle*/

           } else {
               species[ipt] = 3; /* tail particle (hydrophobic) */
           }

           ipt ++;
        }
    }
    linkalkanes();
              
}

void LJBONDFrame::potential()
{
    int i; double tmp;
    lennard_jones_bond();

    /* compute stress, copied from md.cpp calcprop() */
    tmp=0.5*(_ATOMMASS[0]*MASSCONVERT)/(_TIMESTEP*_TIMESTEP);
    _KATOM=0; _PSTRESS.clear();
    _EPOT_BOX=0; _KBOX=0;
    _OMEGA=_H.det();
    
    if(_VR!=NULL)
    {
      for(i=0;i<_NP;i++)
      {
          if(fixed[i]) continue;
          _NPfree ++;
          _VR[i]=_H*_VSR[i];
          if(species[i]==0)
          {
              _KATOM+=tmp*_VR[i].norm2();
              _PSTRESS.addnvv(2.0*tmp,_VR[i],_VR[i]);
          }
          else
          {
              if(_ATOMMASS[species[i]]<=0)
              {
                  FATAL("_ATOMMASS (for species "<<species[i]<<") = "
                       <<_ATOMMASS[species[i]]<<" is not valid");
              }
              _KATOM+=tmp*_VR[i].norm2()*_ATOMMASS[species[i]]/_ATOMMASS[0];
              _PSTRESS.addnvv(2.0*tmp*_ATOMMASS[species[i]]/_ATOMMASS[0],
                              _VR[i],_VR[i]);
          }
      }
    }
    _TOTSTRESS=_PSTRESS+_VIRIAL;
    _TOTSTRESS/=_OMEGA*(1-_VACUUMRATIO); /* in eV/A^3 */
    /* end of calculating stress */
    
    dUdLAM_L=-LMAX*_TOTSTRESS[0][0]*_H[1][1]*_H[2][2];

    dUdLAM_X=-XMAX*_TOTSTRESS[0][0]*_H[1][1]*_H[2][2];
    dUdLAM_Y=-YMAX*_TOTSTRESS[1][1]*_H[2][2]*_H[0][0];
    dUdLAM_Z=-ZMAX*_TOTSTRESS[2][2]*_H[0][0]*_H[1][1];
    dUdLAM_XYZ=dUdLAM_X+dUdLAM_Y+dUdLAM_Z;
    dUdLAM=dUdLAM+dUdLAM_L+dUdLAM_XYZ;

    dEdlambda = dUdLAM;
}

void LJBONDFrame::SWITCHpotential_user(double lambda)
{
    LAMBDA = lambda;
    WALL3_X0 = LAMBDA * WALL3_X0MAX;
    LJBONDFrame::potential();
}
    
void LJBONDFrame::calcprop()
{
    MDPARALLELFrame::calcprop();

    /* extend box in x direction accorrding to lambda */ 
    _H[0][0] = _H0[0][0] + LAMBDA*2.0*WALL3_X0MAX;
}

   
#ifdef _TEST

/* Main Program Begins */
class LJBONDFrame sim;

/* The main program is defined here */
#include "main.cpp"

#endif//_TEST

