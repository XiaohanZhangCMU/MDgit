/*
  meam-lammps.cpp
  by Wei Cai  caiwei@mit.edu
  Last Modified : Wed Nov  7 15:34:37 2007

  FUNCTION  : MD++ with MEAM potential (using LAMMPS meam lib)

  LAMMPS implementation differs from Baskes's code in the two-body potential
   LAMMPS: meam_force.F 3rd order polynomial
           coefficients set up in meam_setup_done.F
           nrar drar  nr 100 -> 1000
   BASKES: phiid()

  LAMMPS does not have EPOT_IND now added

  LAMMPS change meam_force.F  strscomp = 0 -> 1

  Aug. 25 2009 Keonwook Kang
  1. Delete strscomp. Instead, use vflag_atom. 
  2. Delete eflag. Instead, use eflag_either, eflag_global, and eflag_atom
*/

#include "meam-lammps.h"

/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Greg Wagner (SNL)
------------------------------------------------------------------------- */

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define MAXLINE 1024

void MEAMFrame::Alloc()
{
    MDPARALLELFrame::Alloc();
    int size;
    size=_NP*allocmultiple;

    /* MEAM potential storage space */
    Realloc(rho,double,size);
    Realloc(rho0,double,size);
    Realloc(rho1,double,size);
    Realloc(rho2,double,size);
    Realloc(rho3,double,size);
    Realloc(frhop,double,size);
    Realloc(gamma,double,size);
    Realloc(dgamma1,double,size);
    Realloc(dgamma2,double,size);
    Realloc(dgamma3,double,size);
    Realloc(arho2b,double,size);
    Realloc(arho1,Vector3,size);
    Realloc(arho2,Vector6,size);
    Realloc(arho3,Vector10,size);
    Realloc(arho3b,Vector3,size);
    Realloc(t_ave,Vector3,size);
    Realloc(tsq_ave,Vector3,size);
    Realloc(type,int,size);
    Realloc(rtmp,Vector3,size);
    Realloc(vatom,Vector6,size);
    Realloc(num_neigh_full,int,size);
    Realloc(num_neigh_half,int,size);

    bindvar("rho",rho,DOUBLE);
    bindvar("rho0",rho0,DOUBLE);
    bindvar("rho1",rho1,DOUBLE);
    bindvar("rho2",rho2,DOUBLE);
    bindvar("rho3",rho3,DOUBLE);
    bindvar("frhop",frhop,DOUBLE);
    bindvar("gamma",gamma,DOUBLE);
    bindvar("dgamma1",dgamma1,DOUBLE);
    bindvar("dgamma2",dgamma2,DOUBLE);
    bindvar("dgamma3",dgamma3,DOUBLE);
    bindvar("arho2b",arho2b,DOUBLE);
    bindvar("arho1",arho1,DOUBLE);
    bindvar("arho2",arho2,DOUBLE);
    bindvar("arho3",arho3,DOUBLE);
    bindvar("arho3b",arho3b,DOUBLE);
    bindvar("t_ave",t_ave,DOUBLE);
    bindvar("tsq_ave",tsq_ave,DOUBLE);
    
}    

void MEAMFrame::potential()
{
    MEAM();
}

int MEAMFrame::readMEAM()
{
    int select, errorflag;
    /* single component MEAM */

#ifdef _PARALLEL    
    INFO_Printf("[%d]: readMEAM meamfile = %s element[%d]=%s element[%d]=%s\n",
                myDomain,meamfile,0,element[0],1,element[1]);
#else
    INFO_Printf("readMEAM meamfile = %s element[%d]=%s element[%d]=%s\n",
                meamfile,0,element[0],1,element[1]);
#endif
    
    LFile::SubHomeDir(meamfile,meamfile);    
    LFile::SubHomeDir(meafile, meafile);    
    read_files(meamfile,meafile);

    select=10;
    meam_setup_param_(&select,&rcut,NULL,NULL,&errorflag);
    
    meam_setup_done_(&rcut); /* Note: forced Legendre correction ! */

#ifdef _PARALLEL    
    INFO_Printf("[%d] readMEAM: rcut=%f\n",myDomain,rcut);
#else
    INFO_Printf("readMEAM: rcut=%f\n",rcut);
#endif
    
    Realloc(fmap,int,nspecies);
    for(int i=0;i<nspecies;i++) fmap[i]=i+1; /* fortran */
    
    _RLIST = 1.1*rcut;
    _SKIN = _RLIST - rcut;
    return 0;
}

void MEAMFrame::NbrList_reconstruct(int iatom)
{
    MDPARALLELFrame::NbrList_reconstruct();
}

void MEAMFrame::NbrList_translate()
{
    int i, j, jpt, n;
    /* translate nbr list to a new format */
    maxneigh=0; for(i=0;i<_NP;i++) { if (fixed[i]!=-1) maxneigh+=nn[i]; }

    //INFO_Printf("NbrList_reconstruct: maxneigh = %d\n",maxneigh);
    
    /* Allocate 2-dimensional array for neighbor lists */
    Realloc(ind_neigh_full_mem,char,_NP*_NNM*sizeof(int)+_NP*sizeof(int *));
    ind_neigh_full=(int **)(ind_neigh_full_mem+_NP*_NNM*sizeof(int));
    for(i=0;i<_NP;i++) ind_neigh_full[i]=(int *)(ind_neigh_full_mem+i*_NNM*sizeof(int));

    Realloc(ind_neigh_half_mem,char,_NP*_NNM*sizeof(int)+_NP*sizeof(int *));
    ind_neigh_half=(int **)(ind_neigh_half_mem+_NP*_NNM*sizeof(int));
    for(i=0;i<_NP;i++) ind_neigh_half[i]=(int *)(ind_neigh_half_mem+i*_NNM*sizeof(int));

    for(i=0;i<_NP;i++)
    {
        if(fixed[i]==-1) continue; /* -1 means not existing */
        
        /* prepare half neighbor list */
        n = 0;
        for(j=0;j<nn[i];j++)
        {
            jpt=nindex[i][j];
            if(fixed[jpt]==-1) continue; /* -1 means not existing */
            if(jpt>i)
            {
                ind_neigh_half[i][n]=jpt+1;  /* fortran index starts with 1 */
                ind_neigh_full[i][n]=jpt+1;  /* fortran index starts with 1 */
                n++;
            }
        }
        num_neigh_half[i] = n;        
        /* prepare full neighbor list */
        for(j=0;j<nn[i];j++)
        {
            jpt=nindex[i][j];
            if(fixed[jpt]==-1) continue;
            if(jpt<i)
            {
                ind_neigh_full[i][n]=jpt+1;  /* fortran starts with 1 */
                n++;
            }
        }
        num_neigh_full[i] = n;
    }    
}

void MEAMFrame::MEAM()
{
    int i, j, jpt;
    int eflag_either, eflag_global, eflag_atom, vflag_atom, errorflag, offset, ifort;
    Vector3 ds, ds0;
    
    refreshneighborlist();
    NbrList_translate();

    /* rearrange neighbor list into lammps format */
    Realloc(scrfcn, double,maxneigh);
    Realloc(dscrfcn,double,maxneigh);
    Realloc(fcpair, double,maxneigh);
    Realloc(rtmp, Vector3,_NP);

    SHtoR(); memcpy(rtmp,_R,sizeof(Vector3)*_NP);
    eflag_either = 1;
    eflag_global = 1;
    eflag_atom = 1;

//    INFO_Printf("MEAM: fmap[0]=%d\n",fmap[0]);
    
    for(i=0;i<_NP;i++)
    {
        /* do not skip any atom when resetting memory */
        
        _EPOT_IND[i]=0;  /* potential energy per atom */
        _F[i].clear();   /* force on every atom */
        _VIRIAL_IND[i].clear();
        _TORQUE_IND[i]=0;
        _BENDMOMENT_IND[i]=0;

        rho0[i]=0;
        arho2b[i]=0;

        arho1[i].clear();
        arho3b[i].clear();
        t_ave[i].clear();
        tsq_ave[i].clear();
        vatom[i].clear();

        for(j=0;j<6;j++) arho2[i][j]=0;
        for(j=0;j<10;j++) arho3[i][j]=0;
    }
    _EPOT=0;
    _VIRIAL.clear();
    _TORQUE = 0;
    _BENDMOMENT = 0;
    
    errorflag = 0;
    offset = 0;

    for(i=0;i<_NP;i++)
    {
        if(fixed[i]==-1) continue; /* -1 means not existing */        
        type[i] = species[i] + 1; /* fortran */        
    }
    
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]==-1) continue; /* -1 means not existing */
        /* put neighboring atoms to their nearest neighbor positions to atom i */
        ds0=_SR[i]; ds0.subint(); rtmp[i]=_H*ds0; /* needed for parallel run */
        for(j=0;j<nn[i];j++)
        {
            jpt=nindex[i][j]; /* back to C index */
            if(fixed[jpt]==-1) continue; /* -1 means not existing */            
            ds=_SR[jpt]-ds0; /* find nearest neighbor image */
            ds.subint();
            ds+=ds0;
            rtmp[jpt]=_H*ds;
        }
        ifort = i+1;
        meam_dens_init_(&ifort,&_NP,&nspecies,type,fmap,(double *)rtmp,
		    &num_neigh_half[i],ind_neigh_half[i],
                    &num_neigh_full[i],ind_neigh_full[i],
		    &scrfcn[offset],&dscrfcn[offset],&fcpair[offset],
		    rho0,&arho1[0][0],&arho2[0][0],arho2b,
		    &arho3[0][0],&arho3b[0][0],&t_ave[0][0],&tsq_ave[0][0],
                    &errorflag);
//        if(i==0)
//        {
//            INFO_Printf("[%d]: rho0[%d] = %e\n",
//                        myDomain,i,rho0[i]);
//        }

        if (errorflag) {
            FATAL("meam_dens_init library error "<<errorflag);            
        }
        offset += num_neigh_half[i];
    }

//    for(i=0;i<_NP;i++)
//        if((i>=120)&&(i<=130)) INFO_Printf("MEAM-: rho0[%d]=%e rho1[%d]=%e\n",
//                                           i,rho0[i],i,rho1[i]);

    meam_dens_final_(&_NP,&_NP,&eflag_either,&eflag_global,&eflag_atom,
                     &_EPOT,_EPOT_IND,&nspecies,type,fmap,
                     &arho1[0][0],&arho2[0][0],arho2b,&arho3[0][0],
                     &arho3b[0][0],&t_ave[0][0],&tsq_ave[0][0],
                     gamma,dgamma1,dgamma2,dgamma3,
                     rho,rho0,rho1,rho2,rho3,frhop,fixed,&errorflag);

    if (errorflag) {
        FATAL("meam_dens_final library error "<<errorflag);            
    }

    offset = 0;

    for(i=0;i<_NP;i++)
    {
        if(fixed[i]==-1) continue; /* -1 means not existing */
        /* wonder where does LAMMPS implement PBC rules? */
        /* put neighboring atoms to their nearest neighbor positions to atom i */
        ds0=_SR[i]; ds0.subint(); rtmp[i]=_H*ds0; /* needed for parallel run */
        for(j=0;j<num_neigh_full[i];j++)
        {
            jpt=ind_neigh_full[i][j]-1; /* back to C index */
            ds=_SR[jpt]-ds0; /* find nearest neighbor image */
            ds.subint();
            ds+=ds0;
            rtmp[jpt]=_H*ds;
        }
        ifort = i+1;
        vflag_atom = 1;

#ifdef _TORSION_OR_BENDING
        meam_force_(&ifort,&_NP,&_NP0,&eflag_either,$eflag_global,$eflag_atom,
                    &_EPOT,_EPOT_IND,&nspecies,type,fmap,(double *)rtmp,
                    &num_neigh_half[i],ind_neigh_half[i],
                    &num_neigh_full[i],ind_neigh_full[i],
                    &scrfcn[offset],&dscrfcn[offset],&fcpair[offset],
                    dgamma1,dgamma2,dgamma3,rho0,rho1,rho2,rho3,frhop,
                    &arho1[0][0],&arho2[0][0],arho2b,&arho3[0][0],&arho3b[0][0],
                    &t_ave[0][0],&tsq_ave[0][0],&_F[0][0],
                    &vflag_atom, &_VIRIAL_IND[0][0][0],
                    &_TORSIONSIM, _TORQUE_IND,
                    &_BENDSIM, bendspec, _BENDMOMENT_IND,
                    &errorflag);
#else
        meam_force_(&ifort,&_NP,&eflag_either,&eflag_global,&eflag_atom,
                    &vflag_atom,&_EPOT,_EPOT_IND,&nspecies,type,fmap,(double *)rtmp,
                    &num_neigh_half[i],ind_neigh_half[i],
                    &num_neigh_full[i],ind_neigh_full[i],
                    &scrfcn[offset],&dscrfcn[offset],&fcpair[offset],
                    dgamma1,dgamma2,dgamma3,rho0,rho1,rho2,rho3,frhop,
                    &arho1[0][0],&arho2[0][0],arho2b,&arho3[0][0],&arho3b[0][0],
                    &t_ave[0][0],&tsq_ave[0][0],&_F[0][0],&vatom[0][0],
                    &errorflag);
#endif

        if (errorflag) {
            FATAL("meam_force library error "<<errorflag);            
        }
        offset += num_neigh_half[i];
    }

    _TORQUE = 0;
    _BENDMOMENT = 0;
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]==-1) continue; /* -1 means not existing */        
        _VIRIAL_IND[i][0][0] = vatom[i][0];
        _VIRIAL_IND[i][1][1] = vatom[i][1];
        _VIRIAL_IND[i][2][2] = vatom[i][2];
        _VIRIAL_IND[i][0][1] = vatom[i][3]; _VIRIAL_IND[i][1][0] = vatom[i][3];
        _VIRIAL_IND[i][0][2] = vatom[i][4]; _VIRIAL_IND[i][2][0] = vatom[i][4];
        _VIRIAL_IND[i][1][2] = vatom[i][5]; _VIRIAL_IND[i][2][1] = vatom[i][5];
        //_VIRIAL_IND[i] *= -0.5; /* change of convention */
        _VIRIAL += _VIRIAL_IND[i];        
        _TORQUE += _TORQUE_IND[i];        
        _BENDMOMENT += _BENDMOMENT_IND[i];        
        //if(_TORSIONSIM&&(i<5)) INFO_Printf("TORQUE_IND[%d]=%20.12e\n",i,_TORQUE_IND[i]);
        //if(_BENDSIM&&(i<5)) INFO_Printf("BENDMOMENT_IND[%d]=%20.12e\n",i,_BENDMOMENT_IND[i]);
    }
}

void MEAMFrame::initvars()
{
    MDPARALLELFrame::initvars();

    strcpy(meamfile,"meamf");
    strcpy(meafile,"NULL");
}

void MEAMFrame::initparser()
{
    MDPARALLELFrame::initparser();
    bindvar("meamfile",meamfile,STRING);
    bindvar("meafile",meafile,STRING);
    bindvar("rcut",&rcut,DOUBLE);
}

int MEAMFrame::exec(const char *name)
{
    if(MDPARALLELFrame::exec(name)==0) return 0;
    bindcommand(name,"readMEAM",readMEAM());
    bindcommand(name,"printpairpot",printpairpot());
    
#ifdef _PARALLEL
    bindcommand(name,"Broadcast_MEAM_Param",Broadcast_MEAM_Param());
#endif
    return -1;
}


int count_words(char *line)
{
  int n = strlen(line) + 1;
  char *copy = (char *) malloc(n*sizeof(char));
  strcpy(copy,line);

  char *ptr;
  if ((ptr = strchr(copy,'#'))) *ptr = '\0';

  if (strtok(copy," \t\n\r\f") == NULL) {
    free(copy);
    return 0;
  }
  n = 1;
  while (strtok(NULL," \t\n\r\f")) n++;

  free(copy);
  return n;
}

/* ---------------------------------------------------------------- */
void MEAMFrame::read_files(char *globalfile, char *userfile)
{
  // open global meamf file on proc 0
  FILE *fp;

  fp = fopen(globalfile,"r");
  if (fp == NULL) {
      ERROR("Cannot open MEAM potential file"<<globalfile);
  }

  // allocate parameter arrays

  int params_per_line = 19;
  int nelements;

  nelements = nspecies;
  
  int *lat = new int[nelements];
  int *ielement = new int[nelements];
  int *ibar = new int[nelements];
  double *z = new double[nelements];
  double *atwt = new double[nelements];
  double *alpha = new double[nelements];
  double *b0 = new double[nelements];
  double *b1 = new double[nelements];
  double *b2 = new double[nelements];
  double *b3 = new double[nelements];
  double *alat = new double[nelements];
  double *esub = new double[nelements];
  double *asub = new double[nelements];
  double *t0 = new double[nelements];
  double *t1 = new double[nelements];
  double *t2 = new double[nelements];
  double *t3 = new double[nelements];
  double *rozero = new double[nelements];

  bool *found = new bool[nelements];
  for (int i = 0; i < nelements; i++) found[i] = false;

  // read each set of params from global MEAM file
  // one set of params can span multiple lines
  // store params if element name is in element list
  // if element name appears multiple times, only store 1st entry

  int i,n,nwords;
  char **words = new char*[params_per_line+1];
  char line[MAXLINE],*ptr;
  int eof = 0;

  int nset = 0;
  while (1) {
//    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fp);
      if (ptr == NULL) {
	eof = 1;
	fclose(fp);
      } else n = strlen(line) + 1;
//    }
//    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
//    MPI_Bcast(&n,1,MPI_INT,0,world);
//    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // strip comment, skip line if blank

    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    nwords = count_words(line);
    if (nwords == 0) continue;

    // concatenate additional lines until have params_per_line words

    while (nwords < params_per_line) {
      n = strlen(line);
//      if (comm->me == 0) {
        ptr = fgets(&line[n],MAXLINE-n,fp);
        if (ptr == NULL) {
	  eof = 1;
	  fclose(fp);
        } else n = strlen(line) + 1;
//      }
//      MPI_Bcast(&eof,1,MPI_INT,0,world);
      if (eof) break;
//      MPI_Bcast(&n,1,MPI_INT,0,world);
//      MPI_Bcast(line,n,MPI_CHAR,0,world);
      if ((ptr = strchr(line,'#'))) *ptr = '\0';
      nwords = count_words(line);
    }

    if (nwords != params_per_line)
        ERROR("Incorrect format in MEAM potential file");

    // words = ptrs to all words in line
    // strip single and double quotes from words

    nwords = 0;
    words[nwords++] = strtok(line,"' \t\n\r\f");
    while ((words[nwords++] = strtok(NULL,"' \t\n\r\f"))) continue;

    // skip if element name isn't in element list

    for (i = 0; i < nelements; i++)
      if (strcmp(words[0],element[i]) == 0) break;
    if (i == nelements) continue;
    
    // skip if element already appeared

    if (found[i] == true) continue;
    found[i] = true;
    
    // map lat string to an integer

    if (strcmp(words[1],"fcc") == 0) lat[i] = FCC;
    else if (strcmp(words[1],"bcc") == 0) lat[i] = BCC;
    else if (strcmp(words[1],"hcp") == 0) lat[i] = HCP;
    else if (strcmp(words[1],"dim") == 0) lat[i] = DIM;
    else if (strcmp(words[1],"dia") == 0) lat[i] = DIAMOND;
    else ERROR("Unrecognized lattice type in MEAM file 1");

//    INFO_Printf("i=%d  words[0]=%s element[%d]=%s found[%d]=%d lat=%d\n",
//                i,words[0],i,element[i],i,(int)found[i],lat[i]);

    // store parameters

    z[i] = atof(words[2]);
    ielement[i] = atoi(words[3]);
    atwt[i] = atof(words[4]);
    alpha[i] = atof(words[5]);
    b0[i] = atof(words[6]);
    b1[i] = atof(words[7]);
    b2[i] = atof(words[8]);
    b3[i] = atof(words[9]);
    alat[i] = atof(words[10]);
    esub[i] = atof(words[11]);
    asub[i] = atof(words[12]);
    t0[i] = atof(words[13]);
    t1[i] = atof(words[14]);
    t2[i] = atof(words[15]);
    t3[i] = atof(words[16]);
    rozero[i] = atof(words[17]);
    ibar[i] = atoi(words[18]);

    nset++;

//    INFO_Printf("[%d]: nelements=%d element[%d]=%s rozero[%d]=%e\n",
//                myDomain,nelements,i,element[i],i,rozero[i]);
  }

  // error if didn't find all elements in file

  if (nset != nelements)
    ERROR("Did not find all elements in MEAM library file");

  // pass element parameters to MEAM package

  meam_setup_global_(&nelements,lat,z,ielement,atwt,alpha,b0,b1,b2,b3,
  		     alat,esub,asub,t0,t1,t2,t3,rozero,ibar);

  // set element masses

  for (i = 0; i < nelements; i++) _ATOMMASS[i] = atwt[i];

  // clean-up memory

  delete [] words;

  delete [] lat;
  delete [] ielement;
  delete [] ibar;
  delete [] z;
  delete [] atwt;
  delete [] alpha;
  delete [] b0;
  delete [] b1;
  delete [] b2;
  delete [] b3;
  delete [] alat;
  delete [] esub;
  delete [] asub;
  delete [] t0;
  delete [] t1;
  delete [] t2;
  delete [] t3;
  delete [] rozero;
  delete [] found;

  // done if user param file is NULL

  if (strcmp(userfile,"NULL") == 0) return;

  // open user param file on proc 0

//  if (comm->me == 0) {
    fp = fopen(userfile,"r");
    if (fp == NULL) {
        ERROR("Cannot open MEAM potential file "<<userfile);
    }
//  }

  // read settings
  // pass them one at a time to MEAM package
  // match strings to list of corresponding ints

  int which;
  double value;
  int nindex,index[3];
  int maxparams = 6;
  int nparams;
  char *params[maxparams];

  eof = 0;
  while (1) {
//    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fp);
      if (ptr == NULL) {
	eof = 1;
	fclose(fp);
      } else n = strlen(line) + 1;
//    }
//    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
//    MPI_Bcast(&n,1,MPI_INT,0,world);
//    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // strip comment, skip line if blank

    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    nparams = count_words(line);
    if (nparams == 0) continue;

    // words = ptrs to all words in line

    nparams = 0;
    params[nparams++] = strtok(line,"=(), '\t\n\r\f");
    while (nparams < maxparams && 
	   (params[nparams++] = strtok(NULL,"=(), '\t\n\r\f")))
      continue;
    nparams--;

    for (which = 0; which < nkeywords; which++)
      if (strcmp(params[0],keywords[which]) == 0) break;
    if (which == nkeywords) {
        ERROR("Keyword "<<params[0]<<" in MEAM parameter file not recognized");
    }
    nindex = nparams - 2;
    for (i = 0; i < nindex; i++) index[i] = atoi(params[i+1]);

    // map lattce_meam value to an integer

    if (which == 4) {
      if (strcmp(params[nparams-1],"fcc") == 0) value = FCC;
      else if (strcmp(params[nparams-1],"bcc") == 0) value = BCC;
      else if (strcmp(params[nparams-1],"hcp") == 0) value = HCP;
      else if (strcmp(params[nparams-1],"dim") == 0) value = DIM;
      else if (strcmp(params[nparams-1],"dia") == 0) value = DIAMOND;
      else if (strcmp(params[nparams-1],"b1")  == 0) value = B1;
      else if (strcmp(params[nparams-1],"c11") == 0) value = C11;
      else if (strcmp(params[nparams-1],"l12") == 0) value = L12;
      else ERROR("Unrecognized lattice type in MEAM file 2");
    }
    else value = atof(params[nparams-1]);

    // pass single setting to MEAM package

    int errorflag = 0;
    meam_setup_param_(&which,&value,&nindex,index,&errorflag);
    if (errorflag) {
        ERROR("MEAM library error "<<errorflag);
    }
  }
}

void MEAMFrame::printpairpot()
{    
    int elti, eltj;
    double rmin, dr, rmax, r, phi, phip;
    char pairfile[200];
    FILE *fp;
    
    elti = (int) input[0] + 1;
    eltj = (int) input[1] + 1;
    rmin = input[2];
    dr   = input[3];
    rmax = input[4];

    strcpy(pairfile,"pairpot_lammps.dat");
    fp=fopen(pairfile,"w");
    if(fp==NULL)
    {
        FATAL("printpairpot: file "<<pairfile<<" open failure.");
    }
    if((rmax<rmin)||(dr<=0))
        FATAL("rmax cannot be smaller than rmin, dr must be positive");
    for(r=rmin;r<=rmax;r+=dr)
    {
        phif_(&elti,&eltj,&r,&phi,&phip);
        fprintf(fp,"%21.14e  %21.14e %21.14e\n",r,phi,phip);
    }
    fclose(fp);
    INFO("results written to "<<pairfile);
}



#ifdef _PARALLEL
void MEAMFrame::Broadcast_MEAM_Param()
{
    double *buf_double;
    int nparam;
    
    /* master asking slaves to call the same function */
    if(myDomain==0) Master_to_Slave("Broadcast_MEAM_Param");    

//    INFO_Printf("[%d]: Broadcast_MEAM_Param()\n",myDomain);

    /* send meamfile, nspecies, element0, rcut
     * call readMEAM */
    
    /* broadcasting string parameters */
//    INFO_Printf("[%d]: meamfile = %s element[%d]=%s element[%d]=%s\n",
//                myDomain,meamfile,0,element[0],1,element[1]);

    MPI_Bcast(meamfile,1000,MPI_CHAR,0,MPI_COMM_WORLD);
    MPI_Bcast(element,MAXSPECIES*10,MPI_CHAR,0,MPI_COMM_WORLD);

//    INFO_Printf("[%d]: meamfile = %s element[%d]=%s element[%d]=%s\n",
//                myDomain,meamfile,0,element[0],1,element[1]);

    /* packing parameters */    
    nparam = 4;
    buf_double = (double *) malloc(nparam*sizeof(double));

    if(myDomain==0)
    {
        buf_double[0]  = rcut;
        buf_double[1]  = (double) nspecies;
        buf_double[2]  = (double) _NNM;
        buf_double[3]  = (double) _NIC;
    }
    /* broadcasting */
    MPI_Bcast(buf_double,nparam,MPI_DOUBLE,0,MPI_COMM_WORLD);

    if(myDomain!=0)
    {
        /* unpacking parameters */
        rcut = buf_double[0];
        nspecies = (int) buf_double[1];
        _NNM = (int) buf_double[2];
        _NIC = (int) buf_double[3];
        
        readMEAM();
    }

    free(buf_double);
}
#endif//_PARALLEL




#ifdef _TEST

/* Main Program Begins */
class MEAMFrame sim;

/* The main program is defined here */
#include "main.cpp"

#endif//_
