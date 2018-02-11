/*
  eam.cpp
  by Wei Cai  caiwei@mit.edu
  Last Modified : Wed Oct  8 22:42:59 2008

  FUNCTION  : MD++ implementation of EAM/MEAM potential
*/

#include "eam.h"

#define _CUBICSPLINE

inline double spline(double spline_coeff[NGRID][4],int ind, double qq)
{
  double a, b, c, d, qq2, qq3, f;
  a = spline_coeff[ind][0];
  b = spline_coeff[ind][1];
  c = spline_coeff[ind][2];
  d = spline_coeff[ind][3];
  qq2=qq*qq; qq3=qq2*qq;
  f=a+b*qq+c*qq2+d*qq3;
  return f;
}  

inline double spline1(double spline_coeff[NGRID][4],int ind, double qq)
{
  double a, b, c, d, qq2, qq3, fp;
  a = spline_coeff[ind][0];
  b = spline_coeff[ind][1];
  c = spline_coeff[ind][2];
  d = spline_coeff[ind][3];
  qq2=qq*qq; qq3=qq2*qq;
  fp=b+2*c*qq+3*d*qq2;
  return fp;
}


void compute_spline_coeff(int ngrid, double func[],double deriv[],double dr,double spline_coeff[NGRID][4])
{
  int ind;
  double a, b, c, d, f1, p1, A1, A2, dr2, dr3;
  
  dr2=dr*dr; dr3=dr2*dr;

  for(ind=0;ind<ngrid;ind++)
  {
    a = func[ind];
    b = deriv[ind];
    f1 = func[ind+1];
    p1 = deriv[ind+1];
    A1 = f1-a-b*dr;
    A2 = (p1-b)*dr;
    d = (A2-2*A1)/dr3;
    c = (3*A1-A2)/dr2;
    spline_coeff[ind][0] = a;
    spline_coeff[ind][1] = b;
    spline_coeff[ind][2] = c;
    spline_coeff[ind][3] = d;
  }
}

int EAMFrame::readeam()
{
    FILE *fp;
    char string[500];
    int i;
    
    LFile::SubHomeDir(potfile,potfile);
    
    fp=fopen(potfile,"r");
    if(fp==NULL)
    {
        FATAL("EAMFrame::readeam file ("<<potfile<<") not found!");
        return -1;
    }
    fgets(title3,500,fp);
    INFO_Printf("%s\n",title3);

    fgets(string,500,fp);
    sscanf(string,"%d %lf %lf",&ntyp,&rmass,&rlatt);
    fgets(string,500,fp);
    sscanf(string,"%le %le %le %le",&drar,&drhoar,&actual,&rmin);

    /*debug */
    actual*=0.995;
    
    actual2 = actual*actual;

    INFO_Printf("ntyp=%d rmass=%f rlatt=%f\n",ntyp,rmass,rlatt);
    INFO_Printf("dr=%e drho=%e rcut=%e rmin=%e\n",drar,drhoar,actual,rmin);

    _RLIST = 1.1 * actual;
    _SKIN = _RLIST - actual;

    for(i=0;i<eamgrid;i++)
    {
      fgets(string,500,fp);
      sscanf(string,"%le %le %le %le",&(rho[0][i]),&(rhop[0][i]),&(rho[1][i]),&(rhop[1][i]));
    }
    for(i=0;i<eamgrid;i++)
    {
      fgets(string,500,fp);
      sscanf(string,"%le %le %le %le",&(phi[0][i]),&(phip[0][i]),&(phi[1][i]),&(phip[1][i]));
    }
    if(eamfiletype==2 || eamfiletype==4) /* binary potential, 2:type=alloy, 4:type=fs */
    {
        for(i=0;i<eamgrid;i++) 
        {
          fgets(string,500,fp);
          sscanf(string,"%le %le",&(phix[i]),&(phipx[i]));
        }
    }
    for(i=0;i<eamgrid;i++)
    {
      fgets(string,500,fp);
      sscanf(string,"%le %le %le %le",&(frho[0][i]),&(frhop[0][i]),&(frho[1][i]),&(frhop[1][i]));
    }  
    if(eamfiletype==4) /* binary potential, 4:type=fs */
    {
        for(i=0;i<eamgrid;i++)
        {
          fgets(string,500,fp);
          sscanf(string,"%le %le %le %le",&(rho[2][i]),&(rhop[2][i]),&(rho[3][i]),&(rhop[3][i]));
        }
    }
    INFO("Finished reading EAM potential");
    
    /*   Initialize linear grid points reduntantly on each node to save
     *   transfer time.
     */
    for(i=0;i<eamgrid;i++)
    {
        rval[i] = i*drar;
        rhoval[i] = i*drhoar;
    }
    pottype = 1; /* set potential to EAM */

#ifdef _CUBICSPLINE
    /* compute spline coefficients */
    compute_spline_coeff(eamgrid,rho[0],rhop[0],drar,rho_spline[0]);
    compute_spline_coeff(eamgrid,phi[0],phip[0],drar,phi_spline[0]);
    compute_spline_coeff(eamgrid,frho[0],frhop[0],drhoar,frho_spline[0]);

    if(eamfiletype==2 || eamfiletype==4)
    {
      compute_spline_coeff(eamgrid,rho[1],rhop[1],drar,rho_spline[1]);
      compute_spline_coeff(eamgrid,phi[1],phip[1],drar,phi_spline[1]);
      compute_spline_coeff(eamgrid,phix,  phipx,  drar,phix_spline);
      compute_spline_coeff(eamgrid,frho[1],frhop[1],drhoar,frho_spline[1]);
    }

    if(eamfiletype==4)
    {
      compute_spline_coeff(eamgrid,rho[2],rhop[2],drar,rho_spline[2]);
      compute_spline_coeff(eamgrid,rho[3],rhop[3],drar,rho_spline[3]);
    }

#endif
    return 0;
}


#ifdef _TORSION_OR_BENDING
#include "../cookies/src/eam_torsion_bending.cpp"
#else

void EAMFrame::rhoeam()
{
    int i, j, jpt, idx, jdx, ind;
    Vector3 sij,rij;
    double r2ij, rmagg, qq, qr;
    double rhoi, rhoj;

    /*INFO("rhoeam");*/
    petrip = 0;
    //rhocon = 1e-10;
    rhocon = 0.0;            /* replaced by Keonwook Kang, Apr 29, 2011 */
    rhomax = eamgrid*drhoar;

    for(i=0;i<_NP;i++)
    {
        atpe3b[i]=0;
        rhotot[i]=0;
        nbst[i]=0;
    }
    for(i=0;i<_NP;i++)
    {
        _F[i].clear(); _EPOT_IND[i]=0;
        _VIRIAL_IND[i].clear();
    }
    _EPOT=0;
    _VIRIAL.clear();
    
    /*	start calculation
     *  calculate pair contribution to total energy and forces
     */

    /* do on i-particles */
    for(i=0;i<_NP;i++)
    {
        /* modify here for binary systems (0/1) */
        idx = species[i]; /* type of atom (i) */
        /* do on j-particles */
        for(j=0;j<nn[i];j++)
        {
            /* modify here for binary systems (0/1) */
            jpt=nindex[i][j];
            jdx = species[jpt]; /* type of atom (j) */
            if(i>=jpt) continue;
            sij=_SR[jpt]-_SR[i];
            sij.subint();
            rij=_H*sij;
            r2ij=rij.norm2();
            if(r2ij>actual2) continue;
            rmagg=sqrt(r2ij)-rmin;
#if 0
 if (i >=10 && i <=20) { 
	    printf("j = %d, _d_nn[%d]=%d, jpt = %d, ",j, i, nn[i], jpt);
        printf("atom[%d], j=%d, r2ij=%e, d_rmin=%e, actual2 = %e\n",i,j, r2ij, rmin, actual2);
 }
#endif
            nbst[i]++;
            ind = (int)(rmagg/drar);

            if(ind>=eamgrid)
            {
                ind=eamgrid-1;
                INFO_Printf("ind = %d r=%f in RHOEAM\n",ind, rmagg+rmin);
            }
            else if(ind<0)
            {
                ind=0;
                INFO_Printf("ind = %d r=%f in RHOEAM\n",ind, rmagg+rmin);
            }
            qq=rmagg-rval[ind];

#if 0
 if (i >=10 && i <=20) { 
      printf("atom[%d] qq=%e, rmagg=%e\n",i,qq, rmagg);
 }
#endif
            if(idx==jdx)
            {
#ifndef _CUBICSPLINE
            rhoi=rho[jdx][ind]+ qq*rhop[jdx][ind];
#else
            //rhoi = interp(rho[jdx],rhop[jdx],drar,ind,qq);
            rhoi = spline(rho_spline[jdx],ind,qq);
#endif
            rhotot[i]+=rhoi;
            rhotot[jpt]+=rhoi;
            }
            else
            {
              if (eamfiletype == 2)
              {
#ifndef _CUBICSPLINE
              rhoi=rho[jdx][ind]+ qq*rhop[jdx][ind];
              rhoj=rho[idx][ind]+ qq*rhop[idx][ind];
#else
              //rhoi = interp(rho[jdx],rhop[jdx],drar,ind,qq);
              //rhoj = interp(rho[idx],rhop[idx],drar,ind,qq);
              rhoi = spline(rho_spline[jdx],ind,qq);
              rhoj = spline(rho_spline[idx],ind,qq);
#endif
              } else if (eamfiletype == 4)
              {
#ifndef _CUBICSPLINE
              rhoi=rho[idx+2][ind]+ qq*rhop[idx+2][ind];
              rhoj=rho[jdx+2][ind]+ qq*rhop[jdx+2][ind];
#else
              rhoi = spline(rho_spline[idx+2],ind,qq);
              rhoj = spline(rho_spline[jdx+2],ind,qq);
#endif
              }
              rhotot[i]+=rhoi;
              rhotot[jpt]+=rhoj;

            }
        }
    }

    /* debug */
#if 0
    for(i=0;i<_NP;i++)
    {
       INFO_Printf("atom[%d] rhotot=%e\n",i,rhotot[i]);
    }
#endif
    
    for(i=0;i<_NP;i++)
    {
        /* modify here for binary systems (0/1) */
        idx = species[i]; /* type of atom (i) */
        if(rhotot[i]<rhocon)
        {
            rhotot[i]=rhocon;
        }
        if(rhotot[i]>rhomax)
        {
            rhotot[i]=rhomax;
        }
        ind = (int)(rhotot[i]/drhoar);
        if(ind>=eamgrid-1) ind=eamgrid-1;
        qr = rhotot[i] - rhoval[ind];

#ifndef _CUBICSPLINE
        embf[i] = frho[idx][ind] + qr*frhop[idx][ind];
        embfp[i] = frhop[idx][ind] +
          qr*(frhop[idx][ind+1]-frhop[idx][ind])/drhoar;
#else
        //embf[i] = interp(frho[idx],frhop[idx],drhoar,ind,qr);
        //embfp[i] = interp1(frho[idx],frhop[idx],drhoar,ind,qr);
        embf[i] = spline(frho_spline[idx],ind,qr);
        embfp[i] = spline1(frho_spline[idx],ind,qr);
#endif
	if (i <= 3)
	printf("i = %d,embf[i]= %e,embfp[i]=%e\n", i, embf[i], embfp[i]);
        atpe3b[i] = embf[i];
        _EPOT_IND[i]+=atpe3b[i];
        _EPOT+=atpe3b[i];
    }    

    /* debug */
#if 0
    for(i=0;i<_NP;i++)
    {
        INFO_Printf("atom[%d] embf=%e, _EPOT_IND=%e, _EPOT=%e\n",i,embf[i], _EPOT_IND[i], _EPOT);
    }
#endif
}

void EAMFrame::kraeam()
{
    int i, j, jpt, idx, jdx, ind;
    Vector3 sij,rij,fij;
    double r2ij, rmagg, pp, qq, fcp, fpp, fp, denspi, denspj;
    
    /*	start calculation
     *  calculate pair contribution to total energy and forces
     */
    
    /* do on i-particles */
    for(i=0;i<_NP;i++)
    {
        /* modify here for binary systems (0/1) */
        idx = species[i]; /* type of atom (i) */
        /* do on j-particles */
        for(j=0;j<nn[i];j++)
        {
            /* modify here for binary systems (0/1) */
            jpt=nindex[i][j];
            jdx = species[jpt]; /* type of atom (j) */
            if(i>=jpt) continue;
            sij=_SR[jpt]-_SR[i];
            sij.subint();
            rij=_H*sij;
            r2ij=rij.norm2();
            if(r2ij>actual2) continue;
            rmagg=sqrt(r2ij)-rmin;

            ind = int(rmagg/drar);
            
            if(ind>=eamgrid)
            {
                ind=eamgrid-1;
                INFO_Printf("ind = %d in RHOEAM\n",ind);
            }
            else if(ind<0)
            {
                ind=0;
                INFO_Printf("ind = %d in RHOEAM\n",ind);
            }
            qq=rmagg-rval[ind];

            if(idx==jdx)
            {
#ifndef _CUBICSPLINE
              pp = phi[idx][ind] + qq*phip[idx][ind];
              fpp = phip[idx][ind] +
                qq*(phip[idx][ind+1]-phip[idx][ind])/drar;
#else
              //pp = interp(phi[idx],phip[idx],drar,ind,qq);
              //fpp = interp1(phi[idx],phip[idx],drar,ind,qq);
              pp = spline(phi_spline[idx],ind,qq);
              fpp = spline1(phi_spline[idx],ind,qq);
#endif
            }
            else
            {
#ifndef _CUBICSPLINE
              pp = phix[ind] + qq*phipx[ind];
              fpp = phipx[ind] + qq*(phipx[ind+1]-phipx[ind])/drar;
#else
              //pp = interp(phix,phipx,drar,ind,qq);
              //fpp = interp1(phix,phipx,drar,ind,qq);
              pp = spline(phix_spline,ind,qq);
              fpp = spline1(phix_spline,ind,qq);
#endif
            }
            //INFO_Printf("phi = %20.18e\n",pp);
            if ( (idx==jdx) || (eamfiletype==2) )
            {
#ifndef _CUBICSPLINE
            denspi = rhop[idx][ind] +
                qq*(rhop[idx][ind+1]-rhop[idx][ind])/drar ;
            denspj = rhop[jdx][ind] +
                qq*(rhop[jdx][ind+1]-rhop[jdx][ind])/drar ; /* typo idx fixed to jdx ! */
#else
            //denspi = interp1(rho[idx],rhop[idx],drar,ind,qq);
            //denspj = interp1(rho[jdx],rhop[jdx],drar,ind,qq);
            denspi = spline1(rho_spline[idx],ind,qq);
            denspj = spline1(rho_spline[jdx],ind,qq);
#endif
            } else if ( (idx!=jdx) && (eamfiletype==4) ) {
#ifndef _CUBICSPLINE
            denspi = rhop[jdx+2][ind] +
                qq*(rhop[jdx+2][ind+1]-rhop[jdx+2][ind])/drar ;
            denspj = rhop[idx+2][ind] +
                qq*(rhop[idx+2][ind+1]-rhop[idx+2][ind])/drar ;
#else
            denspi = spline1(rho_spline[jdx+2],ind,qq);
            denspj = spline1(rho_spline[idx+2],ind,qq);
#endif
            }

            fcp = denspj * embfp[i] + denspi * embfp[jpt];
            fp = (fpp + fcp) / (rmagg+rmin);
            _EPOT_IND[i]+= 0.5*pp;
            _EPOT_IND[jpt]+= 0.5*pp;
            _EPOT+=pp;
//        printf("atom[%d], j =%d, pp=%e\n",i,j, pp);
            
            fij=rij*fp;
            _F[i]+=fij;
            _F[jpt]-=fij;
#if 0
            if(!(fixed[i]||fixed[jpt]))
                _VIRIAL.addnvv(-fp,rij,rij);
            else if(!(fixed[i]&&fixed[jpt]))
                _VIRIAL.addnvv(-0.5*fp,rij,rij);
#else
	        _VIRIAL_IND[i].addnvv(-.5*fp,rij,rij);
	        _VIRIAL_IND[jpt].addnvv(-.5*fp,rij,rij);
                //_VIRIAL.addnvv(-fp,rij,rij);
		
#if 0
		if (i == 0) { 
			printf("i = %d, j =%d\n", i, j);
	    for(int I = 0;I<3;I++) {
	      for(int J = 0;J<3;J++){
	        printf("_VIRIAL_i[%d][%d]=%e ",I,J, _VIRIAL_IND[i][I][J]); 
	      }
	      printf("\n");
	    }

	    for(int I = 0;I<3;I++) {
	      for(int J = 0;J<3;J++){
	        printf("_VIRIAL_jpt[%d][%d]=%e ",I,J, _VIRIAL_IND[jpt][I][J]); 
	      }
	      printf("\n");
	    }
	    }
#endif
		assert(SURFTEN==0);
                if (SURFTEN==1 && (curstep%savepropfreq)==1) AddnvvtoPtPn(_SR[jpt],rij,rij,-fp);
#endif
        }
        _VIRIAL+=_VIRIAL_IND[i];

    }
    
    /* debug */
#if 1
    for(i=0;i<_NP;i++)
    {
        INFO_Printf("atom[%d] _F=%e,%e,%e, _EPOT_IND=%e, _EPOT=%e\n",i,_F[i].x, _F[i].y, _F[i].z, _EPOT_IND[i], _EPOT);
    }
#endif

}
#endif

int EAMFrame::readMEAM()
{
    pottype = 2; /* set potential to MEAM */
    return -1;
};

void EAMFrame::rhoMEAM()
{

}

void EAMFrame::kraMEAM()
{

}

void EAMFrame::Alloc()
{
    MDPARALLELFrame::Alloc();
    int size;
    size=_NP*allocmultiple;
    
    /* Shared Memory */
    Realloc(atpe3b,double,size);
    Realloc(rhotot,double,size);
    Realloc(embf,double,size);
    Realloc(embfp,double,size);
    Realloc(nbst,int,size);
#ifdef _USECUDA
    cuda_memory_alloc();
#endif
}    

/**********************************************************/
/* Functions for MC simulations moving one atom at a time */
/**********************************************************/

void EAMFrame::MC_Alloc()
{
    int size;
    size=_NP*allocmultiple;
    
    MDPARALLELFrame::MC_Alloc();    

    /* Extra memory to allocate */
    Realloc(atpe3b1,double,size);
    Realloc(rhotot1,double,size);
    Realloc(embf1,double,size);
    Realloc(nbst1,int,size);
    Realloc(ListOfChanged,int,_NNM);

    memcpy(atpe3b1,atpe3b,sizeof(double)*size);
    memcpy(rhotot1,rhotot,sizeof(double)*size);
    memcpy(embf1,embf,sizeof(double)*size);
    memcpy(nbst1,nbst,sizeof(int)*size);

    memset(ListOfChanged,0,sizeof(int)*_NNM);
}    

void EAMFrame::MC_Update()
{
    int i,ipt;
    for (i=0;i<NumOfChanged;i++)
    {
        ipt=ListOfChanged[i];
	atpe3b1[ipt]=atpe3b[ipt];
	rhotot1[ipt]=rhotot[ipt];
	embf1[ipt]  =embf[ipt];
	nbst1[ipt]  =nbst[ipt];
        _SR_MC_SAV[ipt]   =_SR[ipt];
        _EPOT_IND1[ipt]=_EPOT_IND[ipt];
    }

    /* we can make this even faster by only
       updating the neighbor list of atoms around i */
    MC_maxd=MC_maxd_trial;
    _R[MC_atom]=_H*_SR[MC_atom];
    if (MC_maxd > _SKIN*_SKIN/4.0)
    {
      //refreshneighborlist();
      //INFO_Printf("curstep=%d, NeighborList Refreshed\n",curstep);
        NbrList_reconstruct();
        MC_maxd=0;
    }
}

void EAMFrame::MC_Recover()
{
    int i,ipt;
    for (i=0;i<NumOfChanged;i++)
    {
        ipt=ListOfChanged[i];
	atpe3b[ipt]=atpe3b1[ipt];
	rhotot[ipt]=rhotot1[ipt];
	embf[ipt]  =embf1[ipt];
	nbst[ipt]  =nbst1[ipt];
        _SR[ipt]   =_SR_MC_SAV[ipt];
        _EPOT_IND[ipt]=_EPOT_IND1[ipt];
    }
}


/********************************************************************/
/* Functions for Umbrella Sampling in MC simulations                */
/********************************************************************/

void EAMFrame::MC_AllocUMB()
{
    int size;
    size=_NP*allocmultiple;

    MDPARALLELFrame::MC_AllocUMB();    

    Realloc(atpe3bUMB,double,size);
    Realloc(rhototUMB,double,size);
    Realloc(embfUMB,double,size);
}

void EAMFrame::MC_UpdateUMB()
{
    int i, size;
    size=_NP*allocmultiple;

    //INFO_Printf("EAMFrame::MC_UpdateUMB()\n");
    MDPARALLELFrame::MC_UpdateUMB();

    /* Extra information to update */
    for (i=0;i<size;i++)
    {
	atpe3bUMB[i]=atpe3b[i];
	rhototUMB[i]=rhotot[i];
	embfUMB[i]  =embf[i];
    }
}

void EAMFrame::MC_RecoverUMB()
{
    int i;

    //INFO_Printf("EAMFrame::MC_RecoverUMB()\n");
    /* Extra information to recover */
    for (i=0;i<_NP;i++)
    {
	atpe3b[i]=atpe3bUMB[i];
        atpe3b1[i]=atpe3bUMB[i];
        rhotot[i]=rhototUMB[i];
	rhotot1[i]=rhototUMB[i];
	embf[i]=embfUMB[i];
	embf1[i]=embfUMB[i];
    }

    MDPARALLELFrame::MC_RecoverUMB();
}

/* To do: put the following function into md.cpp in the future */

/* This is not real FFS, but only to generate initial configuration for UMB */
int EAMFrame::MC_FFSprocess()
{

  if ( FFSoption > 0 )
  {
    //INFO_Printf("N_lgst_cluster = %d\n",N_lgst_cluster);
    INFO_Printf("react_coord = %d\n",react_coord);
    INFO_Printf("_Lam_array[FFSoption]=%d\n\n",_Lam_array[FFSoption]);
    if ( react_coord >= _Lam_array[FFSoption])
    {
	FFS0_check = 1;
	FFScn.write(this,zipfiles,true);
	return 1;
    }
    else
    { 
	if (FFSoption>1) 
	{
	  if ( react_coord < _Lam_array[FFSoption-2] || curstep > totalsteps+step0)
          {
   	    FFS0_check = -1;
            return -1;
          }
	}
        else 
	{
          if ( (react_coord <= lambda_A) || (curstep > totalsteps+step0) ) 
          {
	    FFS0_check = -1;
	    return -1;
	  }
	}
    }
    return 0;
  }
  else
  {
    if ( react_coord >= _Lam_array[FFSoption] && FFS0_check == 0 )
    {
	FFScn.write(this,zipfiles,true);
	FFS0_check=1;
    }

    if ( FFS0_check==1 && react_coord <= lambda_A ) FFS0_check=0;

    return 0;
  }

}




double EAMFrame::potential_energyonly(int iatom)
{
    WARNING("EAMFrame::potential_energyonly(int iatom) is not implemented!");
    return 0;
} 

double EAMFrame::potential_energyonly_before(int iatom)
{
//    INFO_Printf("EAMbeforeCall\n");
    return 0;
}

double EAMFrame::potential_energyonly_after(int iatom)
{
//    INFO_Printf("EAMafterCall\n");
    return potential_energyonly_change(iatom);
}

void EAMFrame::potential()
{
    refreshneighborlist(); // 2007.Jul.10 Chris Weinberger & Wei Cai
    if(pottype==1)
    {
#ifdef _USECUDA
        rhoeam_cuda();
        kraeam_cuda();
#else
        rhoeam();
        kraeam();
#endif

    }
    else if(pottype==2)
    {
        rhoMEAM();
        kraMEAM();
    }
    else
    {
        WARNING("Unknown pottype ("<<pottype<<")"
                " no calculation performed");
    }
}

void EAMFrame::potential_energyonly()
{
  refreshneighborlist();
  if (pottype==1)
  {
    int i, j, jpt, idx, jdx, ind;
    Vector3 sij,rij, ds, dr;
    double r2ij, rmagg, qq, qr, pp;
    double rhoi, rhoj;

    /*INFO("rhoeam");*/
    petrip = 0;
    //rhocon = 1e-10;
    rhocon = 0.0;            /* replaced by Keonwook Kang, Apr 29, 2011 */
    rhomax = eamgrid*drhoar;

    for(i=0;i<_NP;i++)
    {
        atpe3b[i]=0;
        rhotot[i]=0;
        nbst[i]=0;
    }
    for(i=0;i<_NP;i++)
    {
        _F[i].clear(); _EPOT_IND[i]=0;
    }
    _EPOT=0;
    _VIRIAL.clear();
    
    /*	start calculation
     *  calculate pair contribution to total energy and forces
     */

    /* do on i-particles */
    for(i=0;i<_NP;i++)
    {
        /* modify here for binary systems (0/1) */
        idx = species[i]; /* type of atom (i) */
        /* do on j-particles */
        for(j=0;j<nn[i];j++)
        {
            /* modify here for binary systems (0/1) */
            jpt=nindex[i][j];
            jdx = species[jpt]; /* type of atom (j) */
            if(i>=jpt) continue;
            sij=_SR[jpt]-_SR[i];
            sij.subint();
            rij=_H*sij;
            r2ij=rij.norm2();
            if(r2ij>actual2) continue;
            rmagg=sqrt(r2ij)-rmin;
            
            nbst[i]++;
            ind = (int)(rmagg/drar);

            if(ind>=eamgrid)
            {
                ind=eamgrid-1;
                INFO_Printf("ind = %d r=%f in RHOEAM\n",ind, rmagg+rmin);
            }
            else if(ind<0)
            {
                ind=0;
                INFO_Printf("ind = %d r=%f in RHOEAM\n",ind, rmagg+rmin);
            }
            qq=rmagg-rval[ind];

            if(idx==jdx)
            {
#ifndef _CUBICSPLINE
            rhoi=rho[jdx][ind]+ qq*rhop[jdx][ind];
            pp = phi[idx][ind] + qq*phip[idx][ind];
#else
            rhoi = spline(rho_spline[jdx],ind,qq);
            pp = spline(phi_spline[idx],ind,qq);
#endif
            rhotot[i]+=rhoi;
            rhotot[jpt]+=rhoi;
            }
            else
            {
#ifndef _CUBICSPLINE
              rhoi=rho[jdx][ind]+ qq*rhop[jdx][ind];
              rhoj=rho[idx][ind]+ qq*rhop[idx][ind];
              pp = phix[ind] + qq*phipx[ind];
#else
              //rhoi = interp(rho[jdx],rhop[jdx],drar,ind,qq);
              //rhoj = interp(rho[idx],rhop[idx],drar,ind,qq);
              rhoi = spline(rho_spline[jdx],ind,qq);
              rhoj = spline(rho_spline[idx],ind,qq);
              pp = spline(phix_spline,ind,qq);
#endif
              rhotot[i]+=rhoi;
              rhotot[jpt]+=rhoj;
            }
            _EPOT_IND[i]+= 0.5*pp;
            _EPOT_IND[jpt]+= 0.5*pp;
            _EPOT+=pp;
        }
    }

    /* debug */
#if 0
    for(i=0;i<_NP;i++)
    {
        INFO_Printf("atom[%d] rho=%e\n",i,rhotot[i]);
    }
#endif

    if (_ENABLE_FEXT)
    {
        SHtoR();
        for(i=0;i<_NP;i++)
            _F[i] += _Fext[i];

        if(_SR1==NULL)
        {
            for(i=0;i<_NP;i++)
                _EPOT -= dot(_Fext[i],_R[i]);
        }
        else /* if there is a reference structure */
        {
            for(i=0;i<_NP;i++)
            {
                ds = _SR[i] - _SR1[i];
                dr = _H*ds;
                _EPOT -= dot(_Fext[i],dr);
            }
        }
        _EPOT += _EPOT0_ext;
    }
    

    for(i=0;i<_NP;i++)
    {
        /* modify here for binary systems (0/1) */
        idx = species[i]; /* type of atom (i) */
        if(rhotot[i]<rhocon)
        {
            rhotot[i]=rhocon;
        }
        if(rhotot[i]>rhomax)
        {
            rhotot[i]=rhomax;
        }
        ind = (int)(rhotot[i]/drhoar);
        if(ind>=eamgrid-1) ind=eamgrid-1;
        qr = rhotot[i] - rhoval[ind];

#ifndef _CUBICSPLINE
        embf[i] = frho[idx][ind] + qr*frhop[idx][ind];
#else
        embf[i] = spline(frho_spline[idx],ind,qr);
#endif
        atpe3b[i] = embf[i];
        _EPOT_IND[i]+=atpe3b[i];
        _EPOT_RMV[i]=_EPOT_IND[i];// not accurate fix later.
        _EPOT+=atpe3b[i];
    }    
 }

}

double EAMFrame::potential_energyonly_change(int iatom)
{
/* computing the change of pair potential and embedding potential
   that involve the motion of atom i
*/
//INFO_Printf("begin potE_change\n");
 if (pottype==1) /* eam */
 {
    int i, j, jpt, idx, jdx, ind;
    Vector3 sij,rij;
    double r2ij, rmagg, qq, qr, pp;
    double rhoi, rhoj;

    int ind1;
    Vector3 sij1, rij1, dri0, dsii, drii;
    double r2ij1, rmagg1, qq1, pp1;
    double rhoi1, rhoj1;

    double Enew, Eold, deltaE;
    double rhototoldiatom, atpe3boldiatom,Eoldiatom;
    double rcut2temp;

    /*INFO("rhoeam");*/
    petrip = 0;
    //rhocon = 1e-10;
    rhocon = 0.0;            /* replaced by Keonwook Kang, Apr 29, 2011 */
    rhomax = eamgrid*drhoar;
    deltaE=0;
    rhototoldiatom=rhotot[iatom];
    atpe3boldiatom=atpe3b[iatom];
    Eoldiatom=_EPOT_IND[iatom];
    NumOfChanged=0;

//    Realloc(ListOfChanged,int,nn[iatom]); why error ?
    NumOfChanged++;
    ListOfChanged[0]=iatom;
    
    
    /*	start calculation
     *  calculate pair contribution to total energy and forces
     */

    /* do on i-particles */
        /* modify here for binary systems (0/1) */
        idx = species[iatom]; /* type of atom (i) */
        /* do on j-particles */
        i=iatom;
        dri0=(_H*_SR[i])-_R0[i];
        MC_maxd_trial=max(dri0.norm2(),MC_maxd);

        for(j=0;j<nn[iatom];j++)
        {
            /* modify here for binary systems (0/1) */
            jpt=nindex[i][j];
            jdx = species[jpt]; /* type of atom (j) */

//            if(i>=jpt) continue; // Commented out for single atom code
            if (i==jpt) continue;
            sij=_SR[jpt]-_SR[i];
            sij.subint();
            rij=_H*sij;
            r2ij=rij.norm2();

            rcut2temp=pow( (actual*1.0051+0.5*_TIMESTEP*sqrt(3)), 2);
            if( r2ij> rcut2temp ) continue;

            ListOfChanged[NumOfChanged]=jpt;
	    NumOfChanged++;

            sij1=_SR[jpt]-_SR_MC_SAV[i];
            sij1.subint();
            rij1=_H*sij1;
            r2ij1=rij1.norm2();

            rmagg=sqrt(r2ij)-rmin;
            rmagg1=sqrt(r2ij1)-rmin;

/*          Use this part if nbst is needed in other places.  
            if ( (r2ij<actual2) && (r2ij1>actual2) ) 
            {
		nbst[i]++;
		nbst[jpt]++;
	    }
	    if ( (r2ij>actual2) && (r2ij1<actual2) )
	    {
	    	nbst[i]--;
		nbst[jpt]--;
	    }
*/

            ind = (int)(rmagg/drar);
	    ind1= (int)(rmagg1/drar);

            if(ind>=eamgrid)
            {
                ind=eamgrid-1;
//                INFO_Printf("ind = %d r=%f in RHOEAM\n",ind, rmagg+rmin);
            }
            else if(ind<0)
            {
                ind=0;
//                INFO_Printf("ind = %d r=%f in RHOEAM\n",ind, rmagg+rmin);
            }

            if(ind1>=eamgrid)
            {
                ind1=eamgrid-1;
//                INFO_Printf("ind1 = %d r=%f in RHOEAM\n",ind, rmagg+rmin);
            }
            else if(ind1<0)
            {
                ind1=0;
//                INFO_Printf("ind1 = %d r=%f in RHOEAM\n",ind, rmagg+rmin);
            }

            qq=rmagg-rval[ind];
	    qq1=rmagg1-rval[ind1];

            if(idx==jdx)
            {
#ifndef _CUBICSPLINE
                rhoi=rho[jdx][ind]+ qq*rhop[jdx][ind];
	        rhoi1=rho[jdx][ind1]+qq1*rhop[jdx][ind1];
#else
                //rhoi = interp(rho[jdx],rhop[jdx],drar,ind,qq);
                rhoi = spline(rho_spline[jdx],ind,qq);
	        rhoi1= spline(rho_spline[jdx],ind1,qq1);
#endif
                rhotot[i]+=(rhoi-rhoi1);

// Compute Eold
#ifndef _CUBICSPLINE
                pp1 = phi[jdx][ind1] + qq1*phip[jdx][ind1];
                pp  = phi[jdx][ind]  + qq *phip[jdx][ind];
#else
                pp1 = spline(phi_spline[idx],ind1,qq1);
                pp  = spline(phi_spline[idx],ind,qq);
#endif
	        Eold=atpe3b[jpt];
// End of Compute Eold       
     
                rhotot[jpt]+=(rhoi-rhoi1);

// Compute Enew
                if(rhotot[jpt]<rhocon)
                {
                    rhotot[jpt]=rhocon;
                }
                if(rhotot[jpt]>rhomax)
                {
                rhotot[jpt]=rhomax;
                }
                ind = (int)(rhotot[jpt]/drhoar);
                if(ind>=eamgrid-1) ind=eamgrid-1;
                qr = rhotot[jpt] - rhoval[ind];

#ifndef _CUBICSPLINE
                embf[jpt] = frho[jdx][ind] + qr*frhop[jdx][ind];
#else
                embf[jpt] = spline(frho_spline[jdx],ind,qr);
#endif
                atpe3b[jpt]=embf[jpt];
// End of Compute Enew

                Enew=embf[jpt];
                deltaE+=((Enew-Eold)+(pp-pp1));
        	_EPOT_IND[jpt]+=(Enew-Eold+0.5*(pp-pp1));
                _EPOT_IND[i]+=0.5*(pp-pp1);
            }
            else
            {
#ifndef _CUBICSPLINE
                rhoi=rho[jdx][ind]+ qq*rhop[jdx][ind];
                rhoj=rho[idx][ind]+ qq*rhop[idx][ind];
	        rhoi1=rho[jdx][ind1]+qq1*rhop[jdx][ind1];
                rhoj1=rho[idx][ind1]+qq1*rhop[idx][ind1];
#else
                //rhoi = interp(rho[jdx],rhop[jdx],drar,ind,qq);
                //rhoj = interp(rho[idx],rhop[idx],drar,ind,qq);
                rhoi = spline(rho_spline[jdx],ind,qq);
                rhoj = spline(rho_spline[idx],ind,qq);
  	        rhoi1= spline(rho_spline[jdx],ind1,qq1);
  	        rhoj1= spline(rho_spline[idx],ind1,qq1);
#endif
                rhotot[i]+=(rhoi-rhoi1);
                rhotot[jpt]+=(rhoj-rhoj1);
            }
        }

        if(rhotot[iatom]<rhocon)
        {
            rhotot[iatom]=rhocon;
        }
        if(rhotot[iatom]>rhomax)
        {
            rhotot[iatom]=rhomax;
        }
        ind = (int)(rhotot[iatom]/drhoar);
        if(ind>=eamgrid-1) ind=eamgrid-1;
        qr = rhotot[iatom] - rhoval[ind];

#ifndef _CUBICSPLINE
        embf[iatom] = frho[idx][ind] + qr*frhop[idx][ind];
#else
        embf[iatom] = spline(frho_spline[idx],ind,qr);
#endif
        atpe3b[iatom]=embf[iatom];
	deltaE+= (embf[iatom]-atpe3boldiatom);

	if (_ENABLE_FEXT)
        {
	    dsii = _SR_MC_SAV[i] - _SR[i];
            drii = _H * dsii;
            deltaE += dot(_Fext[i],drii);
        }

        return deltaE;
 }
 else /* meam */
 {
        /* not implemented here */
        return 0;
 }
}


void EAMFrame::initvars()
{
    MDPARALLELFrame::initvars();
   
    strcpy(potfile,"eamdata");
}
void EAMFrame::initparser()
{
    MDPARALLELFrame::initparser();
    bindvar("eamgrid",&eamgrid,INT);
    bindvar("pottype",&pottype,INT);
    bindvar("eamfiletype",&eamfiletype,INT);
}

int EAMFrame::exec(const char *name)
{
    if(MDPARALLELFrame::exec(name)==0) return 0;
    bindcommand(name,"readeam",readeam());
    bindcommand(name,"readMEAM",readMEAM());

#ifdef _PARALLEL
    bindcommand(name,"Broadcast_EAM_Param",Broadcast_EAM_Param());
#endif
#ifdef _USECUDA
    bindcommand(name,"test_saxpy",test_saxpy());
    bindcommand(name,"cuda_memcpy_all",cuda_memcpy_all());
#endif
    return -1;
}
   
#ifdef _PARALLEL
void EAMFrame::Broadcast_EAM_Param()
{
    double *buf_double;
    int nparam, i, j, noffset1, noffset2;

    /* master asking slaves to call the same function */
    if(myDomain==0) Master_to_Slave("Broadcast_EAM_Param");    

    /* broadcasting */
    MPI_Bcast(&pottype,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&eamgrid,1,MPI_INT,0,MPI_COMM_WORLD);

    /* packing parameters */
    noffset1 = 12;
    noffset2 = noffset1 + 14*eamgrid;
    nparam = noffset2 + 7*4*eamgrid;

    //buf_double = (double *) malloc(nparam*sizeof(double));
    buf_double = NULL;
    Realloc(buf_double, double, nparam);
    if(myDomain==0)
    {
    buf_double[0]  = (double) ntyp;
    buf_double[1]  = rmass;
    buf_double[2]  = rlatt;
    buf_double[3]  = drar;
    buf_double[4]  = drhoar;
    buf_double[5]  = actual;
    buf_double[6]  = actual2;
    buf_double[7]  = rmin;
    buf_double[8]  = _RLIST;
    buf_double[9]  = _SKIN;
    buf_double[10] = (double) _NNM;
    buf_double[11] = (double) _NIC;

    for(i=0;i<eamgrid;i++)
    {
      buf_double[i*14 + noffset1 + 0] = rho[0][i];
      buf_double[i*14 + noffset1 + 1] = rhop[0][i];

      buf_double[i*14 + noffset1 + 2] = rho[1][i];
      buf_double[i*14 + noffset1 + 3] = rhop[1][i];

      buf_double[i*14 + noffset1 + 4] = phi[0][i];
      buf_double[i*14 + noffset1 + 5] = phip[0][i];

      buf_double[i*14 + noffset1 + 6] = phi[1][i];
      buf_double[i*14 + noffset1 + 7] = phip[1][i];

      buf_double[i*14 + noffset1 + 8] = frho[0][i];
      buf_double[i*14 + noffset1 + 9] = frhop[0][i];

      buf_double[i*14 + noffset1 +10] = frho[1][i];
      buf_double[i*14 + noffset1 +11] = frhop[1][i];

      buf_double[i*14 + noffset1 +12] = phix[i];
      buf_double[i*14 + noffset1 +13] = phipx[i];
    }
    
    for(i=0;i<eamgrid;i++)
       for(j=0;j<4;j++)
       {
          buf_double[(i*4+j)*7 + noffset2 + 0] = rho_spline[0][i][j];
          buf_double[(i*4+j)*7 + noffset2 + 1] = rho_spline[1][i][j];
          buf_double[(i*4+j)*7 + noffset2 + 2] = phi_spline[0][i][j];
          buf_double[(i*4+j)*7 + noffset2 + 3] = phi_spline[1][i][j];
          buf_double[(i*4+j)*7 + noffset2 + 4] = frho_spline[0][i][j];
          buf_double[(i*4+j)*7 + noffset2 + 5] = frho_spline[1][i][j];
          buf_double[(i*4+j)*7 + noffset2 + 6] = phix_spline[i][j];
       }
    }
    /* broadcasting */
    MPI_Bcast(buf_double,nparam,MPI_DOUBLE,0,MPI_COMM_WORLD);

    if(myDomain!=0)
    {
    /* unpacking parameters
     * the following lines can be generated from the above lines by
     * replace regexp:
     *  \(buf_double\[[0-9]+\]\)[ ]+= \([^;]+\);
     *  \2 = \1;
    */
    ntyp   = (int) buf_double[0];
    rmass  = buf_double[1];
    rlatt  = buf_double[2];
    drar   = buf_double[3];
    drhoar = buf_double[4];
    actual = buf_double[5];
    actual2= buf_double[6];
    rmin   = buf_double[7];
    _RLIST = buf_double[8];
    _SKIN  = buf_double[9]; 
    _NNM = (int) buf_double[10];
    _NIC = (int) buf_double[11];

    for(i=0;i<eamgrid;i++)
    {
      rho[0][i]  = buf_double[i*14 + noffset1 + 0];
      rhop[0][i] = buf_double[i*14 + noffset1 + 1];

      rho[1][i]  = buf_double[i*14 + noffset1 + 2];
      rhop[1][i] = buf_double[i*14 + noffset1 + 3];

      phi[0][i]  = buf_double[i*14 + noffset1 + 4];
      phip[0][i] = buf_double[i*14 + noffset1 + 5];

      phi[1][i]  = buf_double[i*14 + noffset1 + 6];
      phip[1][i] = buf_double[i*14 + noffset1 + 7];

      frho[0][i] = buf_double[i*14 + noffset1 + 8];
      frhop[0][i]= buf_double[i*14 + noffset1 + 9];

      frho[1][i] = buf_double[i*14 + noffset1 +10];
      frhop[1][i]= buf_double[i*14 + noffset1 +11];

      phix[i]    = buf_double[i*14 + noffset1 +12];
      phipx[i]   = buf_double[i*14 + noffset1 +13];
    }

    for(i=0;i<eamgrid;i++)
      for(j=0;j<4;j++)
      {
         rho_spline[0][i][j]  = buf_double[(i*4+j)*7 + noffset2 + 0];
         rho_spline[1][i][j]  = buf_double[(i*4+j)*7 + noffset2 + 1];
         phi_spline[0][i][j]  = buf_double[(i*4+j)*7 + noffset2 + 2];
         phi_spline[1][i][j]  = buf_double[(i*4+j)*7 + noffset2 + 3];
         frho_spline[0][i][j] = buf_double[(i*4+j)*7 + noffset2 + 4];
         frho_spline[1][i][j] = buf_double[(i*4+j)*7 + noffset2 + 5];
         phix_spline[i][j]    = buf_double[(i*4+j)*7 + noffset2 + 6];
      }

    }
    
    free(buf_double);

    for(i=0;i<eamgrid;i++)
    {
        rval[i] = i*drar;
        rhoval[i] = i*drhoar;
    }

}
#endif//_PARALLEL

#ifdef _TEST

/* Main Program Begins */
class EAMFrame sim;

/* The main program is defined here */
#include "main.cpp"

#endif//_TEST





























































/* old code */
#if 0

inline double interp(double func[],double deriv[],double dr,int ind,double qq)
{
    double f, a, b, c, d, f1, p1, A1, A2, dr2, dr3, qq2, qq3;
//    f = func[ind] + qq*deriv[ind];
    dr2=dr*dr; dr3=dr2*dr;
    qq2=qq*qq; qq3=qq2*qq;
    a = func[ind];
    b = deriv[ind];
    f1 = func[ind+1];
    p1 = deriv[ind+1];
    A1 = f1-a-b*dr;
    A2 = (p1-b)*dr;
    d = (A2-2*A1)/dr3;
    c = (3*A1-A2)/dr2;
    f=a+b*qq+c*qq2+d*qq3;
    return f;
}
    
    
inline double interp1(double func[],double deriv[],double dr,int ind,double qq)
{
    double fp, a, b, c, d, f1, p1, A1, A2, dr2, dr3, qq2, qq3;
//    fp = deriv[ind] + qq1/dr*(deriv[ind+1]-deriv[ind]);
    dr2=dr*dr; dr3=dr2*dr;
    qq2=qq*qq; qq3=qq2*qq;
    a = func[ind];
    b = deriv[ind];
    f1 = func[ind+1];
    p1 = deriv[ind+1];
    A1 = f1-a-b*dr;
    A2 = (p1-b)*dr;
    d = (A2-2*A1)/dr3;
    c = (3*A1-A2)/dr2;
    fp=b+2*c*qq+3*d*qq2;
    return fp;
}
#endif
