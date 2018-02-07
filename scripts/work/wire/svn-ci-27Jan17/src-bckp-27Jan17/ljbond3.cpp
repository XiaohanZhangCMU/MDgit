/*
  ljbond3.cpp
  by Wei Cai  caiwei@stanford.edu
  Last Modified : Thu Dec 16 10:05:46 PST 2010

  FUNCTION  :  MD simulation package of Lennard-Jones (LJ) potential
               LJ parameters depend on species
               Atoms can also form covalent bond with neighbors
               covalent bonds have stretching energy and bending energy

  Test cases: scripts/work/fibers/fibers.tcl
*/

#include "ljbond3.h"

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
    bindvar("BOND_K", &BOND_K, DOUBLE); /* BOND_K = E * A / L */
    bindvar("BOND_B", &BOND_B, DOUBLE); /* BOND_B = E * I / L */
    bindvar("Rcut",&LJ_RC,DOUBLE);
    bindvar("usrfile",usrfile,STRING);

    bindvar("substrate_Y",  &_SUBSTRATE_Y,DOUBLE);
    bindvar("substrate_REP",&_SUBSTRATE_REP,DOUBLE);
    bindvar("substrate_ATR",&_SUBSTRATE_ATR,DOUBLE);
    bindvar("substrate_CON",&_SUBSTRATE_CON,DOUBLE);
}

int LJBONDFrame::exec(const char *name)
{
    if(MDPARALLELFrame::exec(name)==0) return 0;
    bindcommand(name,"initLJ",initLJ());
    bindcommand(name,"makefibers",makefibers());
    bindcommand(name,"linkfibers",linkfibers());
    bindcommand(name,"writeatomeyeusr",writeatomeyeusrfile(usrfile));
    return -1;
}

void LJBONDFrame::Alloc()
{
    int size;
    MDPARALLELFrame::Alloc();

    Alloc1();
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

void LJBONDFrame::plot()
{
    int ipt,jpt,j;
    double L;
    double x1,y1,z1,x2,y2,z2,dx,dy,dz,dr;
    Vector3 sri, srj, dsrij, ri, rj;

    MDPARALLELFrame::plot();

    if(win==NULL) return;
    if(!(win->alive)) return;
    
    L=max(_H[0][0],_H[1][1]);
    L=max(L,_H[2][2])*.5;

    SHtoR();
    win->Lock();
    //win->Clear();
    
    /* Draw Bonds (from connectivity matrix) */
    for(ipt=0;ipt<_NP;ipt++)
    {
       for(j=0;j<num_bonds[ipt];j++)
       {
          jpt = bond_index[ipt*MAXNUMBOND + j];
          if(ipt>jpt) continue;

          sri=_SR[ipt];
          if(plot_map_pbc==1) sri.subint();
          ri = _H*sri;
            
          srj=_SR[jpt];
          dsrij = srj - sri;  dsrij.subint();
          srj = sri + dsrij;
          rj = _H*srj;

          if(plot_limits[0])
              if((sri.x<plot_limits[1])||(sri.x>plot_limits[2])
               ||(sri.y<plot_limits[3])||(sri.y>plot_limits[4])
               ||(sri.z<plot_limits[5])||(sri.z>plot_limits[6])
               ||(srj.x<plot_limits[1])||(srj.x>plot_limits[2])
               ||(srj.y<plot_limits[3])||(srj.y>plot_limits[4])
               ||(srj.z<plot_limits[5])||(srj.z>plot_limits[6]))
                        continue;

          x1=ri.x/L;y1=ri.y/L;z1=ri.z/L;
          x2=rj.x/L;y2=rj.y/L;z2=rj.z/L;
          dx=x2-x1;dy=y2-y1;dz=z2-z1;dr=sqrt(dx*dx+dy*dy+dz*dz);
          dx/=dr;dy/=dr;dz/=dr;
          win->DrawLine(x1+dx*atomradius[species[ipt]]/L,
                        y1+dy*atomradius[species[ipt]]/L,
                        z1+dz*atomradius[species[ipt]]/L,
                        x2-dx*atomradius[species[jpt]]/L,
                        y2-dy*atomradius[species[jpt]]/L,
                        z2-dz*atomradius[species[jpt]]/L,
                        colors[MAXCOLORS+1],bondradius/L,1);
       }
    }
    win->Unlock();
    win->Refresh();
}

/* defined in namecolor.c */
int Str2RGB(const char *name, int *r, int *g, int *b);

void LJBONDFrame::writeatomeyeusrfile(char *fname)
{
    FILE *fp;
    int i, j, ipt, jpt, r, g, b;
    
    INFO("MDFrame::writeatomeyeusrfile "<<fname);
    
    fp=fopen(fname,"w");
    if(fp==NULL)
    {
        FATAL("writeatomeyeusrfile: open file failure");
    }

    /* atom color and radius */
    Str2RGB(atomcolor[0],&r,&g,&b);
    fprintf(fp, "%g %g %g  %g\n",r/255.0, g/255.0, b/255.0, atomradius[0]);

    /* Draw Bonds (from connectivity matrix) */
    for(ipt=0;ipt<_NP;ipt++)
    {
       for(j=0;j<num_bonds[ipt];j++)
       {
          jpt = bond_index[ipt*MAXNUMBOND + j];
          if(ipt>jpt) continue;

          Str2RGB(bondcolor,&r,&g,&b);
          fprintf(fp, "%d %d  %g %g %g  %g\n", ipt, jpt, r/255.0, g/255.0, b/255.0, bondradius);
       }
    }
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
    BOND_B  = 2;
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

    _RLIST=LJ_RC+1.1;
    _SKIN=_RLIST-LJ_RC;

}

void LJBONDFrame::lennard_jones_bond()
{
    int i,j,k,ipt,jpt,kpt,bonded;
    double U,r,r2,ri6, R0,U0,K0,ri02,ri06;
    double ALJ_local,BLJ_local,Uc_local,DUDRc_local;   
    double rrij, rrik, rijdotrik, costheta;
    double invrij, invrik, tm1, tm2, tmij, tmik, dhmu;
    Vector3 sij, rij, sik, rik, fij, f3i, f3j, f3k;
    Vector3 rtij, rtik;
    
    DUMP(HIG"Lennard Jones"NOR);
        
    refreshneighborlist();
    _EPOT=0;

    for(i=0;i<_NP;i++)
    {
      /* xiaohan: this removes the interaction between substrate with fe nodes which have species = 2*/
      if (species[i] != 0) 
	continue;

      _F[i].clear(); _EPOT_IND[i]=0;}
    _VIRIAL.clear();

    for(ipt=0;ipt<_NP;ipt++)
    {
      /* xiaohan: this removes the interaction between substrate with fe nodes which have species = 2*/
      if (species[ipt] != 0) 
	continue;

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

    /* covalently bonded interaction - bond stretching */
    if(BOND_K != 0)
    {
     for(ipt=0;ipt<_NP;ipt++)
     {
      /* xiaohan: this removes the interaction between substrate with fe nodes which have species = 2*/
      if (species[ipt] != 0) 
	continue;


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
                
          U = 0.5 * BOND_K * SQR(r-BOND_R0);
          fij = rij * (-1.0*BOND_K*(1-BOND_R0/r)); /* bug fix from old ljbond.cpp */

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

    /* covalently bonded interaction - bond bending */
    if(BOND_B != 0)
    {
     for(ipt=0;ipt<_NP;ipt++)
     {
      /* xiaohan: this removes the interaction between substrate with fe nodes which have species = 2*/
      if (species[ipt] != 0) 
	continue;


       /* only applies to atoms with two bonds */
       if (num_bonds[ipt]!=2) continue;

       jpt = bond_index[ipt*MAXNUMBOND + 0];
       kpt = bond_index[ipt*MAXNUMBOND + 1];
       //INFO_Printf("ljbond: ipt=%d jpt=%d kpt=%d\n",ipt,jpt,kpt);

       sij=_SR[jpt]-_SR[ipt]; sij.subint();
       rij=_H*sij; rrij=rij.norm();
                
       sik=_SR[kpt]-_SR[ipt]; sik.subint();
       rik=_H*sik; rrik=rik.norm();

       rijdotrik = dot(rij, rik);
       costheta = rijdotrik / (rrij * rrik);

       /* adapted from sw.cpp, three-body terms */
       tm1=1.0+costheta;
       /* tm2=tm1*tm1; */
                    
       /* sw.cpp: eterm=tm2*tm3/pgam; */
       /* U = BOND_B*tm2; */ /* bug fix */
       U = BOND_B * tm1;

       /* dhmu=2.*BOND_B*tm1; */ /* bug fix */
       dhmu=BOND_B;

       rtij=rij/rrij; invrij = 1/rrij;
       rtik=rik/rrik; invrik = 1/rrik;
                   
       /* dtm1/dxij = 2*tm1*(invrij * xik/rrik - costheta*invrij * xij/rrij) */
       /* dtm1/dxik = 2*tm1*(invrik * xij/rrij - costheta*invrik * xik/rrik) */

       tmij=dhmu*( invrik-costheta*invrij);
       tmik=dhmu*( invrij-costheta*invrik);
       f3i.clear(); f3i.addnv(tmij,rtij); f3i.addnv(tmik,rtik);

       tmij=dhmu*(        costheta*invrij);
       tmik=dhmu*(-invrij                );
       f3j.clear(); f3j.addnv(tmij,rtij); f3j.addnv(tmik,rtik);

       tmij=dhmu*(-invrik                );
       tmik=dhmu*(        costheta*invrik);
       f3k.clear(); f3k.addnv(tmij,rtij); f3k.addnv(tmik,rtik);

       _F[ipt]+=f3i;
       _F[jpt]+=f3j;
       _F[kpt]+=f3k;
                    
       _VIRIAL.addnvv(1.,f3j,rij);
       _VIRIAL.addnvv(1.,f3k,rik);
                        
       _VIRIAL_IND[ipt].addnvv(0.5,f3j,rij);
       _VIRIAL_IND[jpt].addnvv(0.5,f3j,rij);
       _VIRIAL_IND[ipt].addnvv(0.5,f3k,rik);
       _VIRIAL_IND[kpt].addnvv(0.5,f3k,rik);
                        
       _EPOT_IND[ipt]+=U/3.0;
       _EPOT_IND[jpt]+=U/3.0;
       _EPOT_IND[kpt]+=U/3.0;

       _EPOT_RMV[ipt]+=U;
       _EPOT_RMV[jpt]+=U;
       _EPOT_RMV[kpt]+=U;                

       _EPOT+=U;
    }
   }
}

void LJBONDFrame::linkfibers()
{
    int fiber_type, nfiber, chainlen;
    int ipt, i, j;
    double bond_len;

    fiber_type = (int) input[0]; /* 0: loop (over PBC), 1: segment */
    nfiber     = (int) input[1];
    chainlen   = (int) input[2];
    bond_len   = input[3];

    for(ipt=0;ipt<_NP;ipt++) {
      /* xiaohan: this removes the interaction between substrate with fe nodes which have species = 2*/
      if (species[ipt] != 0) 
	continue;

      group[ipt] = 0;
    }
    if (fiber_type == 0 )
    { /* loop type */
       ERROR("fiber_type "<<fiber_type<<" not implemented!");
       return;
    }

    ipt = 0;
    for(i=0;i<nfiber;i++)
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


void LJBONDFrame::makefibers()
{
    int fiber_type, nfiber, chainlen, fix_bottom, flatten;
    int ipt, i, j;
    double bond_len,d_by_two,offset,orient_pref, bottom_coord;
    Vector3 dr, ds, ddr, dds;
    Matrix33 hinv;

    fiber_type = (int) input[0]; /* 0: loop (over PBC), 1: segment */
    nfiber     = (int) input[1];
    chainlen   = (int) input[2];
    bond_len   = input[3];

    d_by_two=0.04;
    offset=0.2;
    
    _H.clear();
    _H[0][0] = input[4];
    _H[1][1] = input[5];
    _H[2][2] = input[6];

    orient_pref = input[7];
    fix_bottom = (int)input[8];
    bottom_coord = input[9];
    flatten = (int)input[10];

    _NP = nfiber*chainlen;
    Alloc();

    hinv = _H.inv();
    ipt = 0;
    for(i=0;i<nfiber;i++)
    {
        for(j=0;j<chainlen;j++)
        {
           if(j==0)
           {
               /* randomly position the first atom */
               if (flatten == 1) _SR[ipt].set(drand48()-0.5,(drand48()-0.5)*0.05,drand48()-0.5);
               else              _SR[ipt].set(drand48()-0.5,(drand48()-0.5),     drand48()-0.5);

               if (fix_bottom) _SR[ipt].y = bottom_coord;

               /* randomly choose a growth direction */
               dr.set(drand48()-0.5,drand48()-0.5,drand48()-0.5);
               dr.y += orient_pref;
               if (flatten == 1) dr.y = 0.0;
               dr/=dr.norm();  dr*=bond_len; ds = hinv*dr;
           } 
           else 
           {
               /* randomly perturb the growth direction */
               ddr.set(drand48()-0.5,drand48()-0.5,drand48()-0.5);
               ddr/=ddr.norm();  ddr*=bond_len*0.1; dds = hinv*ddr;

               /* add a new atom in the randomly chosen direction */
               _SR[ipt] = _SR[ipt-1] + ds + dds;
           }
           species[ipt] = 0;
           ipt ++;
        }
    }
    linkfibers();
}
    
void LJBONDFrame::substrate_potential()
{
   int ipt;
   Vector3 si, ri, fi;
   double y, esub, fsub;

   /* interaction with smooth substrate */
   if ((_SUBSTRATE_REP != 0) || (_SUBSTRATE_ATR != 0))
   {
    SHtoR(); 
    for (ipt=0;ipt<_NP;ipt++)
    {
      /* xiaohan: this removes the interaction between substrate with fe nodes which have species = 2*/
      if (species[ipt] != 0) 
	continue;


        si=_SR[ipt]; 
        ri=_H*si;
        y=ri.y;
 
        if (y < _SUBSTRATE_Y) {
           esub = 0.5*_SUBSTRATE_REP*(y-_SUBSTRATE_Y)*(y-_SUBSTRATE_Y);
           fsub = _SUBSTRATE_REP*(_SUBSTRATE_Y-y); 
        } else {
           esub = _SUBSTRATE_ATR*(1.0-exp(-0.5*(y-_SUBSTRATE_Y)*(y-_SUBSTRATE_Y)));
           fsub = _SUBSTRATE_ATR*((_SUBSTRATE_Y-y)*exp(-0.5*(y-_SUBSTRATE_Y)*(y-_SUBSTRATE_Y))); 
        }
        _F[ipt].y += fsub;
        _EPOT_IND[ipt] += esub;
        _EPOT += esub;
   }
  }
  else if (_SUBSTRATE_CON != 0)
  { /* Confining (quadratic) potential */
    SHtoR(); 
    for (ipt=0;ipt<_NP;ipt++)
    {
      /* xiaohan: this removes the interaction between substrate with fe nodes which have species = 2*/
      if (species[ipt] != 0) 
	continue;

        si=_SR[ipt]; 
        ri=_H*si;
        y=ri.y;
 
        esub = 0.5*_SUBSTRATE_CON*(y-_SUBSTRATE_Y)*(y-_SUBSTRATE_Y);
        fsub = _SUBSTRATE_CON*(_SUBSTRATE_Y-y); 

        _F[ipt].y += fsub;
        _EPOT_IND[ipt] += esub;
        _EPOT += esub;
   }
  }

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
      /* xiaohan: this removes the interaction between substrate with fe nodes which have species = 2*/
      if (species[i] != 0) 
	continue;

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

    substrate_potential();

    /* end of calculating stress */
}

#ifdef _TEST

/* Main Program Begins */
class LJBONDFrame sim;

/* The main program is defined here */
#include "main.cpp"

#endif//_TEST

