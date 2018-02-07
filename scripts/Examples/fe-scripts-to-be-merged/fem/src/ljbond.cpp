/*
  ljbond.cpp
  by Wei Cai  caiwei@stanford.edu
  Last Modified : Fri Mar  2 18:51:39 2007

  FUNCTION  :  MD simulation package of Lennard-Jones (LJ) potential
               LJ parameters depend on species
               Atoms can also form covalent bond with neighbors

  Test cases: scripts/work/ar/lj2-md.script, ljbond.tcl
*/

#include "ljbond.h"

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
           bindvar(s,&(_C12[i][j]),DOUBLE);
           sprintf(s,"C6_%d%d",i,j);
           bindvar(s,&(_C6[i][j]),DOUBLE);
           sprintf(s,"C9_%d%d",i,j);
           bindvar(s,&(_C9[i][j]),DOUBLE);
       }

    bindvar("BOND_R0",&BOND_R0,DOUBLE);
    bindvar("BOND_K", &BOND_K, DOUBLE);

    bindvar("LX0", &LX0, DOUBLE);
    bindvar("LX1", &LX1, DOUBLE);
    bindvar("LY0", &LY0, DOUBLE);
    bindvar("ALKANE_SY0", &ALKANE_SY0, DOUBLE);

    bindvar("ATTRAC_EPOT",&ATTRAC_EPOT,DOUBLE);
    bindvar("ATTRAC_W",&ATTRAC_W,DOUBLE);
    bindvar("ATTRAC_K",&ATTRAC_K,DOUBLE);

    bindvar("PROBE1_EPOT",&PROBE1_EPOT,DOUBLE);
    bindvar("PROBE1_WX",&PROBE1_WX,DOUBLE);
    bindvar("PROBE1_WY",&PROBE1_WY,DOUBLE);
    bindvar("PROBE1_DY",&PROBE1_DY,DOUBLE);
    bindvar("PROBE1_SY0",&PROBE1_SY0,DOUBLE);
    bindvar("PROBE1_K",&PROBE1_K,DOUBLE);

    bindvar("PROBE2_EPOT",&PROBE2_EPOT,DOUBLE);
    bindvar("PROBE2_WX",&PROBE2_WX,DOUBLE);
    bindvar("PROBE2_WY",&PROBE2_WY,DOUBLE);
    bindvar("PROBE2_DY",&PROBE2_DY,DOUBLE);
    bindvar("PROBE2_SY0",&PROBE2_SY0,DOUBLE);
    bindvar("PROBE2_K",&PROBE2_K,DOUBLE);
    bindvar("PROBE2_KP",&PROBE2_KP,DOUBLE);

    bindvar("PROBE3_EPOT",&PROBE3_EPOT,DOUBLE);
    bindvar("PROBE3_WX0",&PROBE3_WX0,DOUBLE);
    bindvar("PROBE3_WX1",&PROBE3_WX1,DOUBLE);
    bindvar("PROBE3_WY",&PROBE3_WY,DOUBLE);
    bindvar("PROBE3_DY",&PROBE3_DY,DOUBLE);
    bindvar("PROBE3_SY0",&PROBE3_SY0,DOUBLE);
    bindvar("PROBE3_K",&PROBE3_K,DOUBLE);

    bindvar("PROBE4_EPOT",&PROBE4_EPOT,DOUBLE);
    bindvar("PROBE4_WX0",&PROBE4_WX0,DOUBLE);
    bindvar("PROBE4_WX1",&PROBE4_WX1,DOUBLE);
    bindvar("PROBE4_WY",&PROBE4_WY,DOUBLE);
    bindvar("PROBE4_DY",&PROBE4_DY,DOUBLE);
    bindvar("PROBE4_SY0",&PROBE4_SY0,DOUBLE);
    bindvar("PROBE4_K",&PROBE4_K,DOUBLE);

    bindvar("PROBEC_EPOT",&PROBEC_EPOT,DOUBLE);
    bindvar("PROBEC_WX",&PROBEC_WX,DOUBLE);
    bindvar("PROBEC_WY",&PROBEC_WY,DOUBLE);
    bindvar("PROBEC_DY",&PROBEC_DY,DOUBLE);
    bindvar("PROBEC_SY0",&PROBEC_SY0,DOUBLE);
    bindvar("PROBEC_K",&PROBEC_K,DOUBLE);

    bindvar("LJ_FACTOR", &LJ_FACTOR, DOUBLE);
    bindvar("BOND_FACTOR", &BOND_FACTOR, DOUBLE);
    bindvar("ALKANE_FACTOR", &ALKANE_FACTOR, DOUBLE);
    bindvar("LX_FACTOR", &LX_FACTOR, DOUBLE);
    bindvar("ATTRAC_K_FACTOR", &ATTRAC_K_FACTOR, DOUBLE);
    bindvar("PROBE1_WX_FACTOR",&PROBE1_WX_FACTOR, DOUBLE);
    bindvar("PROBE2_WX_FACTOR",&PROBE2_WX_FACTOR, DOUBLE);
    bindvar("PROBE3_WX_FACTOR",&PROBE3_WX_FACTOR, DOUBLE);
    bindvar("PROBE4_WX_FACTOR",&PROBE4_WX_FACTOR, DOUBLE);
    bindvar("PROBEC_WX_FACTOR",&PROBEC_WX_FACTOR, DOUBLE);

    bindvar("LJ_FACTOR_SPEC", &LJ_FACTOR_SPEC, STRING);
    bindvar("BOND_FACTOR_SPEC", &BOND_FACTOR_SPEC, STRING);
    bindvar("ALKANE_FACTOR_SPEC", &ALKANE_FACTOR_SPEC, STRING);
    bindvar("LX_FACTOR_SPEC", &LX_FACTOR_SPEC, STRING);
    bindvar("ATTRAC_K_FACTOR_SPEC", &ATTRAC_K_FACTOR_SPEC, STRING);
    bindvar("PROBE1_WX_FACTOR_SPEC",&PROBE1_WX_FACTOR_SPEC, STRING);
    bindvar("PROBE2_WX_FACTOR_SPEC",&PROBE2_WX_FACTOR_SPEC, STRING);
    bindvar("PROBE3_WX_FACTOR_SPEC",&PROBE3_WX_FACTOR_SPEC, STRING);
    bindvar("PROBE4_WX_FACTOR_SPEC",&PROBE4_WX_FACTOR_SPEC, STRING);
    bindvar("PROBEC_WX_FACTOR_SPEC",&PROBEC_WX_FACTOR_SPEC, STRING);

    bindvar("dEdLJ_FACTOR", &dEdLJ_FACTOR, DOUBLE);
    bindvar("dEdBOND_FACTOR", &dEdBOND_FACTOR, DOUBLE);
    bindvar("dEdALKANE_FACTOR", &dEdALKANE_FACTOR, DOUBLE);
    bindvar("dEdLX_FACTOR", &dEdLX_FACTOR, DOUBLE);
    bindvar("dEdATTRAC_K_FACTOR", &dEdATTRAC_K_FACTOR, DOUBLE);
    bindvar("dEdPROBE1_WX_FACTOR",&dEdPROBE1_WX_FACTOR, DOUBLE);
    bindvar("dEdPROBE2_WX_FACTOR",&dEdPROBE2_WX_FACTOR, DOUBLE);
    bindvar("dEdPROBE3_WX_FACTOR",&dEdPROBE3_WX_FACTOR, DOUBLE);
    bindvar("dEdPROBE4_WX_FACTOR",&dEdPROBE4_WX_FACTOR, DOUBLE);
    bindvar("dEdPROBEC_WX_FACTOR",&dEdPROBEC_WX_FACTOR, DOUBLE);

    bindvar("Rcut",&LJ_RC,DOUBLE);
    bindvar("usrfile",usrfile,STRING);
    bindvar("lipid_shape_spec",lipid_shape_spec,DOUBLE);
}

int LJBONDFrame::exec(const char *name)
{
    if(MDPARALLELFrame::exec(name)==0) return 0;
    bindcommand(name,"initLJ",initLJ());
    bindcommand(name,"makelipids",makelipids());
    bindcommand(name,"linklipids",linklipids());
    bindcommand(name,"makealkanes",makealkanes());
    bindcommand(name,"linkalkanes",linkalkanes());
    bindcommand(name,"makeproteins",makeproteins());
    bindcommand(name,"linkproteins",linkproteins());
    bindcommand(name,"writeatomeyeusr",writeatomeyeusrfile(usrfile));
    bindcommand(name,"update_groupID_by_position",update_groupID_by_position());
    bindcommand(name,"update_groupID_by_energy",update_groupID_by_energy());
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
    /* from ljbond2.cpp */
    int size;

    size = _NP*allocmultiple;
    Realloc(num_bonds, int,size);
    Realloc(bond_index,int,size*MAXNUMBOND);

    memset(num_bonds, 0,sizeof(int)*size);
    memset(bond_index,0,sizeof(int)*size*MAXNUMBOND);
}

void LJBONDFrame::plot()
{
    MDPARALLELFrame::plot();

    /* draw covalent bonds between atoms */
    int ipt, j, jpt;
    double L;
    double x1,y1,z1,x2,y2,z2,dx,dy,dz,dr;
    int r,g,b; double alpha; unsigned ce;
    Vector3 sri, srj, sij, ri, rj;

    if(win==NULL) return;
    if(!(win->alive)) return;
    
    L=max(_H[0][0],_H[1][1]);
    L=max(L,_H[2][2])*.5;

    SHtoR();
    win->Lock();
    //win->Clear();
    
    /* draw atom s*/
    if(BOND_K != 0)
    {
        for(ipt=0;ipt<_NP;ipt++)
        {
            sri=_SR[ipt];
            if(plot_map_pbc==1) sri.subint();
            ri = _H*sri;

            for(j=0;j<num_bonds[ipt];j++)
            {
               jpt = bond_index[ipt*MAXNUMBOND + j];
               if(ipt>jpt) continue;

               sij=_SR[jpt]-_SR[ipt];
               sij.subint();
               srj=sri+sij;
               rj = _H*srj;
                
                if(plot_limits[0])
                    if((sri.x<plot_limits[1])||(sri.x>plot_limits[2])
                       ||(sri.y<plot_limits[3])||(sri.y>plot_limits[4])
                       ||(sri.z<plot_limits[5])||(sri.z>plot_limits[6])
                       ||(srj.x<plot_limits[1])||(srj.x>plot_limits[2])
                       ||(srj.y<plot_limits[3])||(srj.y>plot_limits[4])
                       ||(srj.z<plot_limits[5])||(srj.z>plot_limits[6]))
                        continue;

                    /* only draw O-O bonds */
                    /* if((species[i]!=0)||(species[j]!=0)) continue; */
                    //INFO_Printf("atom %d %d %d %d form bond\n",i, j, i%8,j%8); 
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
    }
    
    win->Unlock();
    win->Refresh();
}

void LJBONDFrame::initvars()
{
    int i, j;
    for(i=0;i<MAXSP;i++)
       for(j=i;j<MAXSP;j++)
       {
           _C12[i][j] = 4;
           _C6 [i][j] = 4;
           _C9 [i][j] = 0;
       }

    BOND_R0 = LJ_LENGTH;
    BOND_K  = 2;
    LJ_RC = 2.37343077641 * LJ_LENGTH; /* initial value, 4th neighbor */
    _RLIST=LJ_RC+1.1;
    _SKIN=_RLIST-LJ_RC;

    strcpy(LJ_FACTOR_SPEC, "constant");
    strcpy(BOND_FACTOR_SPEC, "constant");
    strcpy(ALKANE_FACTOR_SPEC, "constant");
    strcpy(LX_FACTOR_SPEC, "constant");
    strcpy(ATTRAC_K_FACTOR_SPEC, "constant");
    strcpy(PROBE1_WX_FACTOR_SPEC, "constant");
    strcpy(PROBE2_WX_FACTOR_SPEC, "constant");
    strcpy(PROBE3_WX_FACTOR_SPEC, "constant");
    strcpy(PROBE4_WX_FACTOR_SPEC, "constant");
    strcpy(PROBEC_WX_FACTOR_SPEC, "constant");

    MDPARALLELFrame::initvars();
}

void LJBONDFrame::initLJ()
{
    int i, j;

    for(i=0;i<MAXSP;i++)
       for(j=i;j<MAXSP;j++)
       {
	     C12[i][j]=_C12[i][j]*LJ_ENERGY*POW12(LJ_LENGTH);
	     C6 [i][j]=_C6 [i][j]*LJ_ENERGY*POW6(LJ_LENGTH);

             /* 1.05 is here because sigma_SC = 1.05 * sigma
                Ref: Nanoscale, 3, 391-400 (2011)
                     J. Chem. Phys. 108, 7397-7409 (1998)
             */
	     C9 [i][j]=_C9 [i][j]*LJ_ENERGY*CUBE(CUBE((1.05*LJ_LENGTH)));

	     Uc[i][j] = (C12[i][j]*POW6(1/LJ_RC)-C6[i][j])*POW6(1/LJ_RC) + C9[i][j]*CUBE(CUBE(1/LJ_RC));
	     DUDRc[i][j] = -(12*C12[i][j]/POW6(LJ_RC)-6*C6[i][j])/POW6(LJ_RC)/LJ_RC
                           -9*C9[i][j]/CUBE(CUBE(LJ_RC))/LJ_RC;

       }
    /* symmetrize LJ parameter matrix */
    for(i=1;i<MAXSP;i++)
       for(j=0;j<i;j++)
       {
	   _C12[i][j] = _C12[j][i];
           _C6 [i][j] = _C6 [j][i];
           _C9 [i][j] = _C9 [j][i];
	    C12[i][j] =  C12[j][i];
            C6 [i][j] =  C6 [j][i];
            C9 [i][j] =  C9 [j][i];
             Uc[i][j] =   Uc[j][i];
          DUDRc[i][j] =DUDRc[j][i];
       }

    _RLIST=LJ_RC+1.1;
    _SKIN=_RLIST-LJ_RC;

}

void LJBONDFrame::SWITCHpotential_user(double lambda)
{
    _LAMBDA = lambda;

    assign_from_lambda(LJ_FACTOR_SPEC,        &LJ_FACTOR);
    assign_from_lambda(BOND_FACTOR_SPEC,      &BOND_FACTOR);
    assign_from_lambda(ALKANE_FACTOR_SPEC,    &ALKANE_FACTOR);
    assign_from_lambda(LX_FACTOR_SPEC,        &LX_FACTOR);
    assign_from_lambda(ATTRAC_K_FACTOR_SPEC,  &ATTRAC_K_FACTOR);
    assign_from_lambda(PROBE1_WX_FACTOR_SPEC, &PROBE1_WX_FACTOR);
    assign_from_lambda(PROBE2_WX_FACTOR_SPEC, &PROBE2_WX_FACTOR);
    assign_from_lambda(PROBE3_WX_FACTOR_SPEC, &PROBE3_WX_FACTOR);
    assign_from_lambda(PROBE4_WX_FACTOR_SPEC, &PROBE4_WX_FACTOR);
    assign_from_lambda(PROBEC_WX_FACTOR_SPEC, &PROBEC_WX_FACTOR);

    LJBONDFrame::lennard_jones_bond();

    dEdlambda = 0;

    add_to_dEdlambda(LJ_FACTOR_SPEC,        dEdLJ_FACTOR);
    add_to_dEdlambda(BOND_FACTOR_SPEC,      dEdBOND_FACTOR);
    add_to_dEdlambda(ALKANE_FACTOR_SPEC,    dEdALKANE_FACTOR);
    add_to_dEdlambda(LX_FACTOR_SPEC,        dEdLX_FACTOR);
    add_to_dEdlambda(ATTRAC_K_FACTOR_SPEC,  dEdATTRAC_K_FACTOR);
    add_to_dEdlambda(PROBE1_WX_FACTOR_SPEC, dEdPROBE1_WX_FACTOR);
    add_to_dEdlambda(PROBE2_WX_FACTOR_SPEC, dEdPROBE2_WX_FACTOR);
    add_to_dEdlambda(PROBE3_WX_FACTOR_SPEC, dEdPROBE3_WX_FACTOR);
    add_to_dEdlambda(PROBE4_WX_FACTOR_SPEC, dEdPROBE4_WX_FACTOR);
    add_to_dEdlambda(PROBEC_WX_FACTOR_SPEC, dEdPROBEC_WX_FACTOR);

}
    
void LJBONDFrame::calcprop()
{
    MDPARALLELFrame::calcprop();
}

void LJBONDFrame::assign_from_lambda(const char *spec, double *value)
{
    if (strcmp(spec,"constant")==0)
    {
        /* do nothing */
    }
    else if (strcmp(spec,"lambda")==0)
    {
        *value = _LAMBDA;
    }
    else if (strcmp(spec,"1-lambda")==0)
    {
        *value = 1-_LAMBDA;
    }
    else
    {
        INFO_Printf("spec = %s\n",spec);
        ERROR("unrecognized spec");
    }
}

void LJBONDFrame::add_to_dEdlambda(const char *spec, double value)
{
    if (strcmp(spec,"constant")==0)
    {
        /* do nothing */
    }
    else if (strcmp(spec,"lambda")==0)
    {
        dEdlambda += value;
    }
    else if (strcmp(spec,"1-lambda")==0)
    {
        dEdlambda -= value;
    }
    else
    {
        INFO_Printf("spec = %s\n",spec);
        ERROR("unrecognized spec");
    }
}

void LJBONDFrame::lennard_jones_bond()
{
    int i,j,ipt,jpt,k, bonded;
    double U,r,r2,ri6, R0,U0,K0,ri02,ri06;
    double C12_local,C6_local,C9_local,Uc_local,DUDRc_local;   
    Vector3 sij, rij, fij;
    Matrix33 hinv;
    
    DUMP(HIG"Lennard Jones"NOR);

    if(LX_FACTOR != 0.0)
    {
       _H[0][0] = LX0 + (LX1-LX0)*LX_FACTOR;
       _H[1][1] = LX0*LY0/_H[0][0];
    }

    refreshneighborlist();
    
    _EPOT = 0;
    LJ_EPOT = BOND_EPOT = ATTRAC_EPOT = PROBE1_EPOT = PROBE2_EPOT = 0;
    PROBE3_EPOT = PROBE4_EPOT = PROBEC_EPOT = 0;
    dEdLJ_FACTOR = dEdBOND_FACTOR = dEdALKANE_FACTOR = dEdLX_FACTOR = 0;
    dEdATTRAC_K_FACTOR = dEdPROBE1_WX_FACTOR = dEdPROBE2_WX_FACTOR = 0;
    dEdPROBE3_WX_FACTOR = dEdPROBE4_WX_FACTOR = dEdPROBEC_WX_FACTOR = 0;

    hinv = _H.inv();

    for(i=0;i<_NP;i++)
    { _F[i].clear(); _EPOT_IND[i]=0; _EPOT_RMV[i]=0; _VIRIAL_IND[i].clear(); }
    _VIRIAL.clear();

    /* LJ interaction */
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
                C12_local =   C12[species[ipt]][species[jpt]];
                C6_local  =   C6 [species[ipt]][species[jpt]];
                C9_local  =   C9 [species[ipt]][species[jpt]];
                Uc_local =    Uc[species[ipt]][species[jpt]];
		DUDRc_local = DUDRc[species[ipt]][species[jpt]];

                R0 = 0.8;  /* short range cut-off */
                if(r<R0)
                {
                  ri02  = 1./(R0*R0);
                  ri06 = ri02*ri02*ri02;
                  U0 = (7*C12_local*ri06 - 4*C6_local)*ri06 + 5.5*C9_local*ri06*ri02/R0; 
                  K0 = (6*C12_local*ri06 - 3*C6_local)*ri06*ri02 + 4.5*C9_local*ri06*ri02*ri02/R0;
                  U = U0 - K0*r2 - r*DUDRc_local + (LJ_RC*DUDRc_local-Uc_local);
                  fij = rij*(2*K0 + DUDRc_local/r);
                }
                else
		{
		    ri6 = 1./(r2*r2*r2);
		    U = (C12_local*ri6-C6_local)*ri6 + C9_local*ri6/CUBE(r)
                        -r*DUDRc_local+(LJ_RC*DUDRc_local-Uc_local);
		    fij=rij*((12*C12_local*ri6-6*C6_local)*ri6/r2 + 9*C9_local*ri6*ri6*r +DUDRc_local/r);
                }
                fij *= LJ_FACTOR;

                /* only cross LJ interaction with Alkane is scaled by ALKANE_FACTOR */
                if (((species[ipt]==3) || (species[jpt]==3)) && (species[ipt]!=species[jpt]))
                {
                    fij *= ALKANE_FACTOR;
                    _F[ipt]-=fij;
                    _F[jpt]+=fij;
                    _EPOT_IND[ipt]+=U*0.5*LJ_FACTOR*ALKANE_FACTOR;
                    _EPOT_IND[jpt]+=U*0.5*LJ_FACTOR*ALKANE_FACTOR;
                    _EPOT+=U*LJ_FACTOR*ALKANE_FACTOR;
                    LJ_EPOT+=U*LJ_FACTOR*ALKANE_FACTOR;
                    dEdLJ_FACTOR+=U*ALKANE_FACTOR;
                    dEdALKANE_FACTOR+=U*LJ_FACTOR;
                    _VIRIAL.addnvv(1.,fij,rij);

                    _VIRIAL_IND[ipt].addnvv(0.5,fij,rij);
                    _VIRIAL_IND[jpt].addnvv(0.5,fij,rij);
                } else {
                    _F[ipt]-=fij;
                    _F[jpt]+=fij;
                    _EPOT_IND[ipt]+=U*0.5*LJ_FACTOR;
                    _EPOT_IND[jpt]+=U*0.5*LJ_FACTOR;
                    _EPOT+=U*LJ_FACTOR;
                    LJ_EPOT+=U*LJ_FACTOR;
                    dEdLJ_FACTOR+=U;
                    _VIRIAL.addnvv(1.,fij,rij);

                    _VIRIAL_IND[ipt].addnvv(0.5,fij,rij);
                    _VIRIAL_IND[jpt].addnvv(0.5,fij,rij);
                }
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

          /* Note that the covalent BOND interaction is not scaled by ALKANE_FACTOR */

          fij *= BOND_FACTOR;
          _F[ipt]-=fij;
          _F[jpt]+=fij;
          _EPOT_IND[ipt]+=U*0.5*BOND_FACTOR;
          _EPOT_IND[jpt]+=U*0.5*BOND_FACTOR;
          _EPOT+=U*BOND_FACTOR;
          BOND_EPOT+=U*BOND_FACTOR;
          dEdBOND_FACTOR+=U;
          _VIRIAL.addnvv(1.,fij,rij);

          _VIRIAL_IND[ipt].addnvv(0.5,fij,rij);
          _VIRIAL_IND[jpt].addnvv(0.5,fij,rij);
       }
       //INFO_Printf("ljbond: F[%d] = (%e %e %e)\n",ipt,_F[ipt].x,_F[ipt].y,_F[ipt].z);       
    }
  }

  /* add external potential */
  Vector3 si, ri, fi;
  double x, y, e_attrac, fy_attrac, species_fact;
  double Ycut, dYcutdy, sgn_y, e_probe, fx_probe, fy_probe;

  /* Attractive potential for tail particles (repulsive for head and solvent) */ 
  if ((ATTRAC_W > 0 ) && (ATTRAC_K > 0 )) {
   for (ipt=0;ipt<_NP;ipt++)
   {
       /* Attractive (confining) potential only applied to tail particles */
       /* Repulsive  potential applied to head and solvent particles */
       if (species[ipt]==1) species_fact = 1.0;
       else if (species[ipt]==0 || species[ipt]==2) species_fact = -1.0;
       else species_fact = 0.0; /* Alkane molecules (species==3) does not feel ATTRAC potential */

       si=_SR[ipt]; si.subint();
       ri=_H*si;
       y=ri.y;

    if (fabs(y) < ATTRAC_W)
    {
       dEdATTRAC_K_FACTOR  +=        -1.0 * species_fact * ATTRAC_K * SQR( 1.0 - SQR(y/ATTRAC_W) );
       e_attrac  = -1.0 * species_fact * ATTRAC_K_FACTOR * ATTRAC_K * SQR( 1.0 - SQR(y/ATTRAC_W) );
       fy_attrac = -4.0 * species_fact * ATTRAC_K_FACTOR * ATTRAC_K *    ( 1.0 - SQR(y/ATTRAC_W) ) * y/ SQR(ATTRAC_W);

       _F[ipt].y += fy_attrac;
       _EPOT_IND[ipt] += e_attrac;
       ATTRAC_EPOT += e_attrac;
       _EPOT += e_attrac;
       _VIRIAL[1][1] += fy_attrac * y;

       _VIRIAL_IND[ipt][1][1] += fy_attrac * y;
    }
   }
  }

  /* Repulsive potential representing the probe 1 */ 
  double wx1; 
  wx1 = PROBE1_WX * PROBE1_WX_FACTOR;
  if ((wx1 > 0) && (PROBE1_K >0)) {
   for (ipt=0;ipt<_NP;ipt++)
   {
       si=_SR[ipt]; si.y-=PROBE1_SY0; si.subint();
       ri=_H*si;    x=ri.x;   y=ri.y;
  
       e_probe = fx_probe = fy_probe = 0.0;

       /* All species (head, tail, alkane) feel PROBE1 */
       if (fabs(x) < wx1) 
       {
           Ycut = 1.0 / (1.0 + exp( (fabs(y) - PROBE1_WY)/PROBE1_DY ) );
           if (y>0) sgn_y = 1.0; else if (y<0) sgn_y = -1.0; else sgn_y = 0.0;
           dYcutdy = -1.0 * SQR(Ycut) * exp ( (fabs(y) - PROBE1_WY)/PROBE1_DY ) / PROBE1_DY * sgn_y;

           e_probe  = PROBE1_K * SQR( SQR(wx1)-SQR(x) ) / SQR(wx1) * Ycut;
           fx_probe = PROBE1_K*4.0*( SQR(wx1)-SQR(x) ) * x / SQR(wx1) * Ycut;
           fy_probe = -1.0 * PROBE1_K * SQR( SQR(wx1)-SQR(x) ) / SQR(wx1) * dYcutdy;

           /* Double check this expression! */
           dEdPROBE1_WX_FACTOR += PROBE1_K * Ycut * ( 4.0*(SQR(wx1)-SQR(x)) - 2.0*SQR(SQR(wx1)-SQR(x))/SQR(wx1) ) 
                                  * PROBE1_WX / (wx1);
           _F[ipt].x += fx_probe;
           _F[ipt].y += fy_probe;
           _EPOT_IND[ipt] += e_probe;
           PROBE1_EPOT += e_probe;
           _EPOT += e_probe;

           _VIRIAL[0][0] += fx_probe * x;
           _VIRIAL[1][1] += fy_probe * y;

           _VIRIAL_IND[ipt][0][0] += fx_probe * x;
           _VIRIAL_IND[ipt][1][1] += fy_probe * y;
      }
   }
  }

  /* Repulsive potential representing the probe 2 */ 
  update_groupID_by_position();
  double wx2, dEdwf; 
  wx2 = PROBE2_WX * PROBE2_WX_FACTOR;
  if ((wx2 > 0) && (PROBE2_K >0)) {
   for (ipt=0;ipt<_NP;ipt++)
   {
       si=_SR[ipt]; si.y-=PROBE2_SY0; si.subint();
       ri=_H*si;    x=ri.x;   y=ri.y;
  
       e_probe = fx_probe = fy_probe = 0.0;

      if ( fabs(x) < _H[0][0]/4.0 ) 
      {
       if ( (((group[ipt]==2) && (x<0)) || ((group[ipt]==1) && (x>0))) 
            && ((species[ipt]==1)||(species[ipt]==2)) && (fabs(y) < PROBE2_WY+PROBE2_DY*3.0) )
       {
          Ycut = 1.0 / (1.0 + exp( (fabs(y) - PROBE2_WY)/PROBE2_DY ) );
          if (y>0) sgn_y = 1.0; else if (y<0) sgn_y = -1.0; else sgn_y = 0.0;
          dYcutdy = -1.0 * SQR(Ycut) * exp ( (fabs(y) - PROBE2_WY)/PROBE2_DY ) / PROBE2_DY * sgn_y;

          e_probe = PROBE2_K * (SQR(wx2)+0.5*SQR(wx2)*fabs(x)*PROBE2_KP) * Ycut;
          double sgn_x;
          if (x>0) sgn_x = 1.0; else if (x<0) sgn_x = -1.0; else sgn_x = 0.0;
          fx_probe = -PROBE2_K * SQR(wx2) * sgn_x *PROBE2_KP * Ycut; 
          fy_probe = -1.0 * PROBE2_K * (SQR(wx2)+0.5*SQR(wx2)*fabs(x)*PROBE2_KP) * dYcutdy;
          /* Double check this expression! */
          dEdwf = PROBE2_K * Ycut * (2.0*wx2+wx2*fabs(x)*PROBE2_KP) * PROBE2_WX;

          dEdPROBE2_WX_FACTOR += dEdwf;
          INFO_Printf("PROBE2: atom %d group = %d  x = %g  fx_probe = %g dEdwf = %g\n", ipt,group[ipt],x,fx_probe, dEdwf);
       } 
       else if (fabs(x) < wx2) 
       {
           Ycut = 1.0 / (1.0 + exp( (fabs(y) - PROBE2_WY)/PROBE2_DY ) );
           if (y>0) sgn_y = 1.0; else if (y<0) sgn_y = -1.0; else sgn_y = 0.0;
           dYcutdy = -1.0 * SQR(Ycut) * exp ( (fabs(y) - PROBE2_WY)/PROBE2_DY ) / PROBE2_DY * sgn_y;

           e_probe  = PROBE2_K * SQR( SQR(wx2)-SQR(x) ) / SQR(wx2) * Ycut;
           fx_probe = PROBE2_K*4.0*( SQR(wx2)-SQR(x) ) * x / SQR(wx2) * Ycut;
           fy_probe = -1.0 * PROBE2_K * SQR( SQR(wx2)-SQR(x) ) / SQR(wx2) * dYcutdy;

           /* Double check this expression! */
           dEdPROBE2_WX_FACTOR += PROBE2_K * Ycut * ( 4.0*(SQR(wx2)-SQR(x)) - 2.0*SQR(SQR(wx2)-SQR(x))/SQR(wx2) ) 
                                  * PROBE2_WX / (wx2);
       }

       _F[ipt].x += fx_probe;
       _F[ipt].y += fy_probe;
       _EPOT_IND[ipt] += e_probe;
       PROBE2_EPOT += e_probe;
       _EPOT += e_probe;

       _VIRIAL[0][0] += fx_probe * x;
       _VIRIAL[1][1] += fy_probe * y;

       _VIRIAL_IND[ipt][0][0] += fx_probe * x;
       _VIRIAL_IND[ipt][1][1] += fy_probe * y;
      }
   }
  }

  /* Repulsive potential representing the probe 3 */ 
  double wx3; 
  wx3 = PROBE3_WX0 + (PROBE3_WX1-PROBE3_WX0) * PROBE3_WX_FACTOR;
  if (PROBE3_K >0) {
   for (ipt=0;ipt<_NP;ipt++)
   {
       si=_SR[ipt]; si.y-=PROBE3_SY0; si.subint();
       ri=_H*si;    x=ri.x;   y=ri.y;
  
       e_probe = fx_probe = fy_probe = 0.0;

       if ( fabs(x) < _H[0][0]/4.0 ) 
       {
          Ycut = 1.0 / (1.0 + exp( (fabs(y) - PROBE3_WY)/PROBE3_DY ) );
          if (y>0) sgn_y = 1.0; else if (y<0) sgn_y = -1.0; else sgn_y = 0.0;
          dYcutdy = -1.0 * SQR(Ycut) * exp ( (fabs(y) - PROBE3_WY)/PROBE3_DY ) / PROBE3_DY * sgn_y;

          if (x<0)
          {
             if (x > -wx3)
             {
                 e_probe = 0.5* PROBE3_K * SQR(x + wx3) * Ycut;
                 fx_probe = -PROBE3_K * (x + wx3) * Ycut; 
                 fy_probe = -1.0 * PROBE3_K * SQR(x + wx3) * dYcutdy;
                 dEdPROBE3_WX_FACTOR += PROBE3_K * Ycut * (wx3 + x) * (PROBE3_WX1 - PROBE3_WX0);
             } 
          }
          else
          {
             if (x < wx3)
             {
                 e_probe = 0.5* PROBE3_K * SQR(x - wx3) * Ycut;
                 fx_probe = -PROBE3_K * (x - wx3) * Ycut; 
                 fy_probe = -1.0 * PROBE3_K * SQR(x - wx3) * dYcutdy;
                 dEdPROBE3_WX_FACTOR += PROBE3_K * Ycut * (wx3 - x) * (PROBE3_WX1 - PROBE3_WX0);
             } 
          }

          _F[ipt].x += fx_probe;
          _F[ipt].y += fy_probe;
          _EPOT_IND[ipt] += e_probe;
          PROBE3_EPOT += e_probe;
          _EPOT += e_probe;

          _VIRIAL[0][0] += fx_probe * x;
          _VIRIAL[1][1] += fy_probe * y;

          _VIRIAL_IND[ipt][0][0] += fx_probe * x;
          _VIRIAL_IND[ipt][1][1] += fy_probe * y;
       }
   }
  }

  /* Repulsive potential representing the probe 4 */ 
  update_groupID_by_energy();
  double wx4; 
  wx4 = PROBE4_WX0 + (PROBE4_WX1-PROBE4_WX0) * PROBE4_WX_FACTOR;
  if (PROBE4_K >0) {
   for (ipt=0;ipt<_NP;ipt++)
   {
       si=_SR[ipt]; si.y-=PROBE4_SY0; si.subint();
       ri=_H*si;    x=ri.x;   y=ri.y;
  
       e_probe = fx_probe = fy_probe = 0.0;

       if ( fabs(x) < _H[0][0]/4.0 ) 
       {
          Ycut = 1.0 / (1.0 + exp( (fabs(y) - PROBE4_WY)/PROBE4_DY ) );
          if (y>0) sgn_y = 1.0; else if (y<0) sgn_y = -1.0; else sgn_y = 0.0;
          dYcutdy = -1.0 * SQR(Ycut) * exp ( (fabs(y) - PROBE4_WY)/PROBE4_DY ) / PROBE4_DY * sgn_y;

          if (   ( (species[ipt]!=1) && (species[ipt]!=2) && (x<0) ) 
              || ( ((species[ipt]==1)||(species[ipt]==2)) && (group[ipt]==1) ) )
          {
             if (x > -wx4)
             {
                 e_probe = 0.5* PROBE4_K * SQR(x + wx4) * Ycut;
                 fx_probe = -PROBE4_K * (x + wx4) * Ycut; 
                 fy_probe = -1.0 * PROBE4_K * SQR(x + wx4) * dYcutdy;
                 dEdwf = PROBE4_K * Ycut * (wx4 + x) * (PROBE4_WX1 - PROBE4_WX0);
                 dEdPROBE4_WX_FACTOR += dEdwf;
                 //INFO_Printf("PROBE4: %4d %4d %8e %8e %8e %8e %8e\n", ipt, group[ipt], x, y, wx4, Ycut, dEdwf);
             } 
          }
          else if (   ( (species[ipt]!=1) && (species[ipt]!=2) && (x>0) ) 
                   || ( ((species[ipt]==1)||(species[ipt]==2)) && (group[ipt]==2) ) )
          {
             if (x < wx4)
             {
                 e_probe = 0.5* PROBE4_K * SQR(x - wx4) * Ycut;
                 fx_probe = -PROBE4_K * (x - wx4) * Ycut; 
                 fy_probe = -1.0 * PROBE4_K * SQR(x - wx4) * dYcutdy;
                 dEdwf = PROBE4_K * Ycut * (wx4 - x) * (PROBE4_WX1 - PROBE4_WX0);
                 dEdPROBE4_WX_FACTOR += dEdwf;
                 //INFO_Printf("PROBE4: %4d %4d %8e %8e %8e %8e %8e\n", ipt, group[ipt], x, y, wx4, Ycut, dEdwf);
             } 
          }

          _F[ipt].x += fx_probe;
          _F[ipt].y += fy_probe;
          _EPOT_IND[ipt] += e_probe;
          PROBE4_EPOT += e_probe;
          _EPOT += e_probe;

          _VIRIAL[0][0] += fx_probe * x;
          _VIRIAL[1][1] += fy_probe * y;

         _VIRIAL_IND[ipt][0][0] += fx_probe * x;
         _VIRIAL_IND[ipt][1][1] += fy_probe * y;
       }
   }
  }

  /* Repulsive potential representing the probe C that only acts 
     on the center of mass of lipid and solvent molecules */ 
  double wxC;  Vector3 dsj, dscom;
  wxC = PROBEC_WX * PROBEC_WX_FACTOR;
  if ((wxC > 0) && (PROBEC_K >0)) {
   for (ipt=0;ipt<_NP;ipt++)
   {
       if (species[ipt]==0) {
          /* solvent molecules */
          si=_SR[ipt]; si.y-=PROBEC_SY0; si.subint();
          ri=_H*si;    x=ri.x;   y=ri.y;
       } else if (species[ipt]==2) {
          /* head particle of lipid, find center of mass of entire lipid molecule */
          si=_SR[ipt];  dscom.clear();
          for(j=1;j<LIPID_CHAINLEN;j++)
          {
              if( species[ipt+j] != 1)
                 FATAL("atoms following head particle must be tail particles!");
              dsj=_SR[ipt+j]-_SR[ipt]; dsj.subint();
              dscom+=dsj;
          }
          dscom/=LIPID_CHAINLEN;
          si+=dscom;

          si.y-=PROBEC_SY0; si.subint();
          ri=_H*si;    x=ri.x;   y=ri.y;
       }

       e_probe = fx_probe = fy_probe = 0.0;
  
       if (fabs(x) < wxC) 
       {
           Ycut = 1.0 / (1.0 + exp( (fabs(y) - PROBEC_WY)/PROBEC_DY ) );
           if (y>0) sgn_y = 1.0; else if (y<0) sgn_y = -1.0; else sgn_y = 0.0;
           dYcutdy = -1.0 * SQR(Ycut) * exp ( (fabs(y) - PROBEC_WY)/PROBEC_DY ) / PROBEC_DY * sgn_y;

           e_probe  = PROBEC_K * SQR( SQR(wxC)-SQR(x) ) / SQR(wxC) * Ycut;
           fx_probe = PROBEC_K*4.0*( SQR(wxC)-SQR(x) ) * x / SQR(wxC) * Ycut;
           fy_probe = -1.0 * PROBEC_K * SQR( SQR(wxC)-SQR(x) ) / SQR(wxC) * dYcutdy;

           if (species[ipt]==0) {
              /* solvent molecules */
              dEdPROBEC_WX_FACTOR += PROBEC_K * Ycut * ( 4.0*(SQR(wxC)-SQR(x)) - 2.0*SQR(SQR(wxC)-SQR(x))/SQR(wxC) ) 
                                     * PROBEC_WX / (wxC);
              _F[ipt].x += fx_probe;
              _F[ipt].y += fy_probe;
              _EPOT_IND[ipt] += e_probe;
              PROBEC_EPOT += e_probe;
              _EPOT += e_probe;
              _VIRIAL[0][0] += fx_probe * x;
              _VIRIAL[1][1] += fy_probe * y;

              _VIRIAL_IND[ipt][0][0] += fx_probe * x;
              _VIRIAL_IND[ipt][1][1] += fy_probe * y;

           } else if (species[ipt]==2) {
             /* head particle of lipid */
              for(j=0;j<LIPID_CHAINLEN;j++)
              {
                  dEdPROBEC_WX_FACTOR += PROBEC_K * Ycut * ( 4.0*(SQR(wxC)-SQR(x)) - 2.0*SQR(SQR(wxC)-SQR(x))/SQR(wxC) ) 
                                         * PROBEC_WX / (wxC);
                  _F[ipt+j].x += fx_probe;
                  _F[ipt+j].y += fy_probe;
                  _EPOT_IND[ipt+j] += e_probe;
                  PROBEC_EPOT += e_probe;
                  _EPOT += e_probe;
                  _VIRIAL[0][0] += fx_probe * x;
                  _VIRIAL[1][1] += fy_probe * y;

                  _VIRIAL_IND[ipt][0][0] += fx_probe * x;
                  _VIRIAL_IND[ipt][1][1] += fy_probe * y;
              }
           }
      }
   }
  }

  /* only works for rectangular simulation box, i.e. diagonal _H matrix */
  dEdLX_FACTOR = (_VIRIAL[1][1] - _VIRIAL[0][0]) * (LX1 - LX0) / _H[0][0];

}

void LJBONDFrame::linklipids()
{
    int nsolvent, nlipid, chainlen, headid;
    int ipt, i, j;
    double bond_len;
    Vector3 dr, dr1, ds;
    Matrix33 hinv;

    nsolvent = (int) input[3];
    nlipid   = (int) input[4];
    chainlen = (int) input[5];
    headid   = (int) input[6];
    bond_len =       input[7];

    LIPID_CHAINLEN = chainlen;

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

    nsolvent = (int) input[3];
    nlipid   = (int) input[4];
    chainlen = (int) input[5];
    headid   = (int) input[6];
    bond_len =       input[7];
    nalkane  = (int) input[8];
    alkanelen= (int) input[9];
    ipt = nlipid*chainlen+nsolvent;

    LIPID_CHAINLEN = chainlen;
   
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

void LJBONDFrame::linkproteins()
{
    int nsolvent, nlipid, nalkane, chainlen, alkanelen, headid;
    int nprotein, top_len, center_len, end_len;
    int ipt, i, j;
    double bond_len;
    Vector3 dr, dr1, ds;
    Matrix33 hinv;

    nsolvent = (int) input[3];
    nlipid   = (int) input[4];
    chainlen = (int) input[5];
    headid   = (int) input[6];
    bond_len =       input[7];
    nalkane  = (int) input[8];
    alkanelen= (int) input[9];
    nprotein = (int) input[10];
    top_len  = (int) input[11];
    center_len=(int) input[12];
    end_len  = (int) input[13];

    LIPID_CHAINLEN = chainlen;

    ipt = nlipid*chainlen+nsolvent+nalkane*alkanelen;
   
    for(i=0;i<nprotein;i++)
    {
        for(j=0;j<top_len+center_len+end_len;j++)
        {
           if(j==0)
           {
               /* only one bond */
               num_bonds[ipt] = 1;
               bond_index[ipt*MAXNUMBOND + 0] = ipt+1;
               
           } else if (j==(top_len+center_len+end_len-1))
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
    double bond_len,d_by_two;
    Vector3 dr, dr1, ds;
    Matrix33 hinv;
    int make_bilayer;
    double bilayer_half_thickness;

    /* Note format change of input[] */
    _H.clear();
    _H[0][0] = input[0];
    _H[1][1] = input[1];
    _H[2][2] = input[2];

    nsolvent = (int) input[3];
    nlipid   = (int) input[4];
    chainlen = (int) input[5];
    headid   = (int) input[6];
    bond_len =       input[7];

    LIPID_CHAINLEN = chainlen;

    /* from ljbond2.cpp */
    nalkane =   (int) input[8];
    alkanelen = (int) input[9];
    d_by_two=0.025;

    make_bilayer = (int) lipid_shape_spec[0];
    bilayer_half_thickness = lipid_shape_spec[1];

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
               if (make_bilayer) {
                   _SR[ipt].set(drand48()-0.5,
                                bilayer_half_thickness/_H[1][1]*(i<nlipid/2?1:-1),
                                drand48()-0.5);
                   dr.set(0,(i<nlipid/2?-1:1),0);
               } else {
                   /* randomly position the first atom */
                   _SR[ipt].set(drand48()-0.5,drand48()-0.5,drand48()-0.5);
                   /* randomly choose a growth direction */
                   dr.set(drand48()-0.5,drand48()-0.5,drand48()-0.5);
               }
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
               dr1/=dr1.norm(); dr1*=(bond_len*0.5); dr+=dr1;
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
        if (make_bilayer) {
             _SR[ipt].set(drand48()-0.5,
                          (drand48()-0.5)*(1-2.0*bilayer_half_thickness/_H[1][1])+0.5,
                          drand48()-0.5);
             _SR[ipt].subint();
        } else {
             /* randomly position the solvent particle */
             _SR[ipt].set(drand48()-0.5,drand48()-0.5,drand48()-0.5);
        }
        species[ipt] = 0;     /* solvent particle (water) */
        ipt ++;
    }             
}

void LJBONDFrame::makealkanes()
{
    int nsolvent, nlipid, nalkane, chainlen, alkanelen, headid;
    int ipt, i, j;
    double bond_len,d_by_two,offset;
    Vector3 dr, dr1, ds;
    Matrix33 hinv;

    nsolvent = (int) input[3];
    nlipid   = (int) input[4];
    chainlen = (int) input[5];
    headid   = (int) input[6];
    bond_len =       input[7];
    nalkane  = (int) input[8];
    alkanelen= (int) input[9];

    LIPID_CHAINLEN = chainlen;

    d_by_two=0.032;
    offset=0.2;

    INFO_Printf("old NP=%d",_NP);
    _NP = nlipid*chainlen + nalkane*alkanelen + nsolvent;
    INFO_Printf("new NP=%d",_NP);
    Alloc1();


    hinv = _H.inv();
    ipt = nlipid*chainlen+nsolvent;
  
    for(i=0;i<nalkane;i++)
    //for(i=0;i<(nalkane/4);i++)
    {
        for(j=0;j<alkanelen;j++)
        {
           if(j==0)
           {
               /* randomly position the first atom */
               if (i<(nalkane/4)) 
               {
                  _SR[ipt].set(5.0/_H[0][0],d_by_two+ALKANE_SY0, 0.5-0.0973*i);
                  /* randomly choose a growth direction */
                  dr.set(1,0,0);
               }
               else if (i<(nalkane/4)*2)
               {
                  _SR[ipt].set(5.0/_H[0][0],-d_by_two+ALKANE_SY0, 0.5-0.0973*i);
                  dr.set(1,0,0);
               }
               else if (i<(nalkane/4)*3)
               {
                  _SR[ipt].set(-5.0/_H[0][0],d_by_two+ALKANE_SY0, 0.45-0.0973*i);
                  dr.set(-1,0,0);
               }
               else
               {
                  _SR[ipt].set(-5.0/_H[0][0],-d_by_two+ALKANE_SY0, 0.45-0.0973*i);
                  dr.set(-1,0,0);
               }

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
void LJBONDFrame::makeproteins()
{
    int nsolvent, nlipid, nalkane, chainlen, alkanelen, headid;
    int nprotein, top_len, center_len, end_len;
    int ipt, i, j;
    double bond_len,d_by_two,offset;
    Vector3 dr, dr1, ds;
    Matrix33 hinv;

    nsolvent = (int) input[3];
    nlipid   = (int) input[4];
    chainlen = (int) input[5];
    headid   = (int) input[6];
    bond_len =       input[7];
    nalkane  = (int) input[8];
    alkanelen= (int) input[9];
    nprotein = (int) input[10];
    top_len  = (int) input[11];
    center_len=(int) input[12];
    end_len  = (int) input[13];

    LIPID_CHAINLEN = chainlen;

    d_by_two=0.032;
    offset=0.2;

    INFO_Printf("old NP=%d",_NP);
    _NP = nlipid*chainlen + nalkane*alkanelen + nsolvent + nprotein*(top_len+center_len+end_len);
    INFO_Printf("new NP=%d",_NP);
    Alloc1();


    hinv = _H.inv();
    ipt = nlipid*chainlen+nsolvent+nalkane*alkanelen;
  
    for(i=0;i<nprotein;i++)
    {
        for(j=0;j<top_len+center_len+end_len;j++)
        {
           if(j==0)
           {
               /* randomly position the first atom */
               _SR[ipt].set(drand48()-0.5,0.3+0.3*drand48(),drand48()-0.5);

               /* randomly choose a growth direction */
               dr.set(0.1*(drand48()-0.5),-0.5*(drand48()),0.1*(drand48()-0.5));
               dr/=dr.norm();  dr*=bond_len; ds = hinv*dr;

           } else {
               /* add a new atom in the randomly chosen direction */
               _SR[ipt] = _SR[ipt-1] + ds;
           }

           if(j<top_len)
           {
               species[ipt] = 4; /* top section of protein */

           } 
           else if (j>=top_len + center_len)
           {
               species[ipt] = 6; /* end section of protein */
           } 
           else 
           {
               species[ipt] = 5; /* center section of protein */
           }

           ipt ++;
        }
    }
    linkproteins();
}


void LJBONDFrame::update_groupID_by_position()
{
   int ipt, j;
   double x, y;
   Vector3 si, ri, dsj, dscom;

   /* Only applies to lipid (head-tail) particles */
   for (ipt=0;ipt<_NP;ipt++)
   {
       if (species[ipt]==2) {
          /* head particle of lipid, find center of mass of entire lipid molecule */
          si=_SR[ipt];  dscom.clear();
          for(j=1;j<LIPID_CHAINLEN;j++)
          {
              if( species[ipt+j] != 1)
                 FATAL("atoms following head particle must be tail particles!");
              dsj=_SR[ipt+j]-_SR[ipt]; dsj.subint();
              dscom+=dsj;
          }
          dscom/=LIPID_CHAINLEN;
          si+=dscom;

          si.subint();
          ri=_H*si;    x=ri.x;   y=ri.y;
       }
       if (species[ipt]==1 || species[ipt]==2) {
          if(fabs(x*4.0) > _H[0][0])
          {
             group[ipt] == 0;
          }
          else
          {
             // allow group ID to change at every time step
             // if(group[ipt]==0)
             // to do: group id should be the one that gives the lower potential
             {
                if(x<0)
                   group[ipt]=1;
                else
                   group[ipt]=2;
             }
          }
       }
   }
}

void LJBONDFrame::update_groupID_by_energy()
{
  int ipt, j;
  Vector3 si, ri;
  double wx4, x, y, Ycut, e_probe_1, e_probe_2; 

  wx4 = PROBE4_WX0 + (PROBE4_WX1-PROBE4_WX0) * PROBE4_WX_FACTOR;

  if (PROBE4_K >0) {
   for (ipt=0;ipt<_NP;ipt++)
   {
      if( species[ipt] == 2 )
      {
          //INFO_Printf("update_groupID_by_energy: ipt = %d\n",ipt);
          si=_SR[ipt]; si.y-=PROBEC_SY0; si.subint();
          ri=_H*si;    x=ri.x;   y=ri.y;
          if (fabs(x*4.0)>_H[0][0])
          {
              for(j=0;j<LIPID_CHAINLEN;j++)
                 group[ipt+j] = 0;
          }
          else {
            e_probe_1 = e_probe_2 = 0.0;
            for(j=0;j<LIPID_CHAINLEN;j++)
            {
              //INFO_Printf("update_groupID_by_energy: ipt = %d j = %d\n",ipt,j);
              if (ipt+j>_NP)
                 FATAL("exceed maximum number of atoms");
              if( j>0 && species[ipt+j] != 1)
                 FATAL("atoms following head particle must be tail particles!");

              si=_SR[ipt+j]; si.y-=PROBEC_SY0; si.subint();
              ri=_H*si;    x=ri.x;   y=ri.y;
  
              Ycut = 1.0 / (1.0 + exp( (fabs(y) - PROBE4_WY)/PROBE4_DY ) );

              if (x > -wx4)
                 e_probe_1 += 0.5* PROBE4_K * SQR(x + wx4) * Ycut;
              if (x < wx4)
                 e_probe_2 += 0.5* PROBE4_K * SQR(x - wx4) * Ycut;
            }
            if (e_probe_1 < e_probe_2)
            {
              for(j=0;j<LIPID_CHAINLEN;j++)
                 group[ipt+j] = 1;
            }
            else {
              for(j=0;j<LIPID_CHAINLEN;j++)
                 group[ipt+j] = 2;
            }
         }
     }
  }
 }
}

void LJBONDFrame::potential()
{
    lennard_jones_bond();
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

    /* Draw Bonds (from connectivity matrix) */
    for(ipt=0;ipt<_NP;ipt++)
    {
        /* atom color and radius */
        Str2RGB(atomcolor[species[ipt]],&r,&g,&b);
        fprintf(fp, "%d  %g %g %g  %g\n",ipt, r/255.0, g/255.0, b/255.0, atomradius[species[ipt]]);
    }

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

#ifdef _TEST

/* Main Program Begins */
class LJBONDFrame sim;

/* The main program is defined here */
#include "main.cpp"

#endif//_TEST

