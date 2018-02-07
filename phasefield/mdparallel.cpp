/*
  mdparallel.cpp
  by Wei Cai  caiwei@stanford.edu
  Last Modified : Thu Apr 29 09:21:25 2010

  FUNCTION  :  Easy-to-use MD simulation Framework

  Featuring :  1. Scripting input
               2. X-window display
               3. automatic simulation log
               4. convenient configuration and output file handle
               5. perfect lattice and dislocation config creator
               6. conjugate-gradient relaxation
               7. NVT MD simulation (Nose-Hoover thermostat)
               8. NPH and NPT simulation (Parrinello-Rahman + Nose-Hoover)
               9. Cell-Verlist combined neighbor list

  To Do:
      1.  debug: run_parallel does not give identical trajectory
          for runs/si-meam-lammps/relaxed2.cn
      2.  eval_parallel gives identical VIRIAL for relaxed2.cn only
          when domain skin layer is 1.4 of RLIST
      3.  Identical trajectory between serial and parallel runs on planck (mpicc and mpigpp)
          and on mc-cc (mpich, upto 16 cpus for Si NW 47580 atoms)
*/


#include "mdparallel.h"


/********************************************************************/
/* Begining of class MDFrame definition */
/********************************************************************/

#ifdef _PARALLEL
/********************************************************************/
/* Initialize script file parser */
/********************************************************************/
void MDPARALLELFrame::initvars()
{
    MDFrame::initvars();

    INFO_Printf("[%d] numDomains = %d alloc requests and status\n",myDomain,numDomains);
    inRequests  = (MPI_Request *) malloc(numDomains*sizeof(MPI_Request));
    outRequests = (MPI_Request *) malloc(numDomains*sizeof(MPI_Request));
    inStatus    = (MPI_Status  *) malloc(numDomains*sizeof(MPI_Status ));
    outStatus   = (MPI_Status  *) malloc(numDomains*sizeof(MPI_Status));
}

void MDPARALLELFrame::initparser()
{
    MDFrame::initparser();
    
    bindvar("myIX",&myIX,INT);
    bindvar("myIY",&myIY,INT);
    bindvar("myIZ",&myIZ,INT);
    bindvar("nXdoms",&nXdoms,INT);
    bindvar("nYdoms",&nYdoms,INT);
    bindvar("nZdoms",&nZdoms,INT);
    bindvar("myDomain",&myDomain,INT);
    bindvar("numDomains",&numDomains,INT);
    bindvar("myXmin",&myXmin,DOUBLE);
    bindvar("myYmin",&myYmin,DOUBLE);
    bindvar("myZmin",&myZmin,DOUBLE);    
    bindvar("myXmax",&myXmax,DOUBLE);
    bindvar("myYmax",&myYmax,DOUBLE);
    bindvar("myZmax",&myZmax,DOUBLE);
    
}    

int MDPARALLELFrame::exec(const char *name)
{
    if(MDFrame::exec(name)==0) return 0;

    bindcommand(name,"Master_to_Slave",Master_to_Slave(command));
    bindcommand(name,"Broadcast_H",{if(myDomain==0) Master_to_Slave("Broadcast_H"); Broadcast_H();});
    bindcommand(name,"Broadcast_Atoms",{if(myDomain==0) Master_to_Slave("Broadcast_Atoms"); Broadcast_Atoms();});
    bindcommand(name,"Broadcast_nebsetting",{if(myDomain==0) Master_to_Slave("Broadcast_nebsetting"); Broadcast_nebsetting();});
    bindcommand(name,"Slave_to_Master_Atoms",{if(myDomain==0) Master_to_Slave("Broadcast_Atoms");Slave_to_Master_Atoms();});
    bindcommand(name,"Slave_chdir",{if(myDomain==0) Master_to_Slave("Slave_chdir");Slave_chdir();});
    bindcommand(name,"Partition_Domains",{if(myDomain==0) Master_to_Slave("Partition_Domains");Partition_Domains();});
    bindcommand(name,"Mark_Local_Atoms",{if(myDomain==0) Master_to_Slave("Mark_Local_Atoms");Mark_Local_Atoms();});
    bindcommand(name,"eval_parallel",{if(myDomain==0) Master_to_Slave("eval_parallel");eval_parallel();});
    bindcommand(name,"run_parallel",{if(myDomain==0) Master_to_Slave("run_parallel");run_parallel();});
    bindcommand(name,"alloc_all",alloc_all());
    bindcommand(name,"quit_all",{if(myDomain==0) Master_to_Slave("quit");quit_all();});
    bindcommand(name,"nebrelax_parallel",{if(myDomain==0) Master_to_Slave("nebrelax_parallel");nebrelax_parallel();});
    bindcommand(name,"stringrelax_parallel",{if(myDomain==0) Master_to_Slave("stringrelax_parallel");stringrelax_parallel();});

    bindcommand(name,"stringrelax_parallel_2",{if(myDomain==0) Master_to_Slave("stringrelax_parallel_2");stringrelax_parallel_2();});
    bindcommand(name,"reparametrization_with_trimming",{if(myDomain==0) Master_to_Slave("reparametrization_with_trimming");reparametrization_with_trimming();});
    bindcommand(name,"allocchain_parallel",{if(myDomain==0) Master_to_Slave("allocchain_parallel");AllocChain_parallel();});
    bindcommand(name,"initRchain_parallel",{if(myDomain==0) Master_to_Slave("initRchain_parallel");initRchain_parallel();});
    bindcommand(name,"readRchain_parallel",{if(myDomain==0) Master_to_Slave("readRchain_parallel");readRchain_parallel();});
    bindcommand(name,"writeRchain_parallel",{if(myDomain==0) Master_to_Slave("writeRchain_parallel");writeRchain_parallel();});
    bindcommand(name,"copyRchaintoCN_parallel",{if(myDomain==0) Master_to_Slave("copyRchaintoCN_parallel");copyRchaintoCN_parallel();});
    bindcommand(name,"copyCNtoRchain_parallel",{if(myDomain==0) Master_to_Slave("copyCNtoRchain_parallel");copyCNtoRchain_parallel();});
    bindcommand(name,"copyCNtoRchain_masteronly",copyCNtoRchain_masteronly());
    bindcommand(name,"readcnfile_parallel",{if(myDomain==0) Master_to_Slave("readcnfile_parallel");readcnfile_parallel();});
    bindcommand(name,"writefinalcnfile_parallel",{if(myDomain==0) Master_to_Slave("writefinalcnfile_parallel");writefinalcnfile_parallel(1,false);});
    bindcommand(name,"writeatomeyecfg_parallel",{if(myDomain==0) Master_to_Slave("writeatomeyecfg_parallel");writeatomeyecfgfile_parallel();});
    return -1;
}

/* Memory Allocation */
void MDPARALLELFrame::Alloc()
{
    int size;
    size=_NP*allocmultiple;

    MDFrame::Alloc();
    
    Realloc(domainID,int,size);
    Realloc(_EPOT_IND_global,double,size);
    Realloc(_F_global,Vector3,size);
    Realloc(_VIRIAL_IND_global,Matrix33,size);
    
    bindvar("domainID",domainID,INT);
}

void MDPARALLELFrame::ParallelInit(int *pargc, char ***pargv)
{
    MPI_Init(pargc,pargv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myDomain);
    MPI_Comm_size(MPI_COMM_WORLD, &numDomains);

    if(myDomain==0)
    {
        INFO_Printf("[%d]: I am the master in charge of %d processors.\n",
                    myDomain, numDomains);
    }
    else
    {
        INFO_Printf("[%d]: Processor %d at your service\n",
                    myDomain, myDomain);
    }
}

void MDPARALLELFrame::WaitForCommand()
{
    FILE *istream, *ostream; int mypipe[2];
                                                                     
    INFO_Printf("[%d]: I am waiting for master's command.\n",myDomain);

    while(1)
    {
        if (pipe (mypipe))
        {
            fprintf (stderr, "[%d]: Pipe failed.\n",myDomain);
            return;
        }
        istream = fdopen(mypipe[0], "r");
        ostream = fdopen(mypipe[1], "w");
        MPI_Bcast(command,MAXCMDLEN,MPI_CHAR,0,MPI_COMM_WORLD);
        fprintf(ostream,"%s \n",command);
        fprintf(ostream,"\n");
        fclose(ostream);
        parse_line(istream);
        fclose(istream);
    }
}

void MDPARALLELFrame::Master_to_Slave(const char *cmd)
{
    if(command!=cmd) strcpy(command,cmd);
    MPI_Bcast(command,MAXCMDLEN,MPI_CHAR,0,MPI_COMM_WORLD);
}

void MDPARALLELFrame::alloc_all()
{
    /* issue command to all slave processors to alloc */
    sprintf(command,"NP = %d",_NP);
    Master_to_Slave(command);
}

void MDPARALLELFrame::quit_all()
{
    /* issue command to all slave processors to alloc */
    //Master_to_Slave("quit");

    if(numDomains>1)
    {
        free(domBoundX);
        free(domBoundY);
        free(domBoundZ);
    }
    
    quit();
}

void MDPARALLELFrame::Partition_Domains()
{
    Vector3 h;
    int i, j, k, i2, j2, k2, ndom, *flag;
    
    /* master asking slaves to call the same function */
    //if(myDomain==0) Master_to_Slave("Partition_Domains");

    MPI_Bcast(&nXdoms,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&nYdoms,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&nZdoms,1,MPI_INT,0,MPI_COMM_WORLD);

    if(myDomain==0)
    {
        INFO_Printf("nXdoms = %d nYdoms = %d nZdoms = %d numDomains = %d\n",
             nXdoms,nYdoms,nZdoms,numDomains);   
        if(nXdoms*nYdoms*nZdoms!=numDomains)
            FATAL("product of nXdoms nYdoms nZdoms ("<<
                  nXdoms*nYdoms*nZdoms<<" does not match numDomains");
    }

    myIZ = myDomain % nZdoms;
    myIY = ((myDomain-myIZ)/nZdoms) % nYdoms;
    myIX = ((myDomain-myIZ)/nZdoms - myIY) / nYdoms;

    if(myDomain == 0) /* verify if partition is valid */
    {
        h=_H.height();
        if ((h[0]/nXdoms) < _RLIST) FATAL("nXdoms ("<<nXdoms<<" too large, RLIST="<<_RLIST);
        if ((h[1]/nYdoms) < _RLIST) FATAL("nYdoms ("<<nYdoms<<" too large, RLIST="<<_RLIST);
        if ((h[2]/nZdoms) < _RLIST) FATAL("nZdoms ("<<nZdoms<<" too large, RLIST="<<_RLIST);
    }

    /* decide which domains are neighbors */    
    flag = (int *) malloc(numDomains*sizeof(int));
    memset(flag,0,numDomains*sizeof(int));
    for(i=-1;i<=1;i++) for(j=-1;j<=1;j++) for(k=-1;k<=1;k++)
    {
        i2=(myIX+i+nXdoms) % nXdoms;
        j2=(myIY+j+nYdoms) % nYdoms;
        k2=(myIZ+k+nZdoms) % nZdoms;        
        
        ndom = (i2*nYdoms + j2)*nZdoms + k2;
        flag[ndom] = 1;
    }    
    neighDoms = (int *) malloc(numDomains*sizeof(int));
    numNeighDoms = 0;
    for(i=0;i<numDomains;i++)
    {
        if(i==myDomain) continue;
        if(flag[i]==1)
        {
            neighDoms[numNeighDoms] = i;
            numNeighDoms++;
        }
    }
    INFO_Printf("[%d]: %d neighbor domains (",myDomain,numNeighDoms);
    for(i=0;i<numNeighDoms;i++)
        INFO_Printf("%d ",neighDoms[i]);
    INFO_Printf(")\n");
    
    free(flag);
}

void MDPARALLELFrame::Broadcast_H()
{
    /* broad cast matrix H */
    MPI_Bcast(&(_H[0][0]),9,MPI_DOUBLE,0,MPI_COMM_WORLD);
}

void MDPARALLELFrame::Broadcast_Atoms()
{
    int NP_old;
    
    /* master asking slaves to call the same function */
    //if(myDomain==0) Master_to_Slave("Broadcast_Atoms");

    /* broad cast number of atoms */
    NP_old = _NP;
    MPI_Bcast(&_NP,1,MPI_INT,0,MPI_COMM_WORLD);
    if(_NP!=NP_old) Alloc();
    
    /* broad cast _SAVEMEMORY */
    MPI_Bcast(&_SAVEMEMORY,1,MPI_INT,0,MPI_COMM_WORLD);
    
    /* broad cast matrix H */
    MPI_Bcast(&(_H[0][0]),9,MPI_DOUBLE,0,MPI_COMM_WORLD);

    /* broad cast atom positions */
    MPI_Bcast(_SR, _NP*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
    /* broad cast atom attributes */
    MPI_Bcast(fixed,  _NP,MPI_INT,0,MPI_COMM_WORLD);
    
    /* broad cast neighborlist settings */
    MPI_Bcast(&_NNM, 1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&_NIC, 1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    
    if(_SAVEMEMORY>=9) return;
    
    /* broad cast atom species */
    MPI_Bcast(species,_NP,MPI_INT,0,MPI_COMM_WORLD);
    
    if(_SAVEMEMORY>=7) return;
    
    /* broad cast atom velocities */
    MPI_Bcast(_VSR,_NP*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
    
    if(_SAVEMEMORY>=6) return;
    
    /* broad cast atom group */
    MPI_Bcast(group,  _NP,MPI_INT,0,MPI_COMM_WORLD);

}

void MDPARALLELFrame::Slave_to_Master_Atoms()
{
    int i, nsend, n, idom, *nrecv;
    double *outBuf, **inBuf;
    
    /* master asking slaves to call the same function */
    //if(myDomain==0) Master_to_Slave("Slave_to_Master_Atoms");

    MPI_Barrier(MPI_COMM_WORLD);
    
    /* determine the number of atoms to send/recv */
    nrecv = NULL;
    if(myDomain!=0) /* slave */
    {
        nsend = 0; /* find out how many atoms to send */
        for(i=0;i<_NP;i++)
        {
            if(fixed[i]==-1) continue; /* -1 means not existing */
            if(domainID[i]==myDomain)
            {
                nsend++;
            }
        }
        MPI_Isend(&nsend,1,MPI_INT,0,MSG_ATOM_LEN,MPI_COMM_WORLD,&outRequests[0]);
//        INFO_Printf("[%d]: Slave_to_Master_Atoms--(a)\n",myDomain);
    }
    else /* master */
    {
        nrecv = (int *)malloc(numDomains*sizeof(int));
        for(idom=1;idom<numDomains;idom++)
            MPI_Irecv(&nrecv[idom],1,MPI_INT,idom,MSG_ATOM_LEN,MPI_COMM_WORLD,&inRequests[idom]);
//        INFO_Printf("[%d]: Slave_to_Master_Atoms--(b)\n",myDomain);
    }
    if(myDomain!=0)
        MPI_Waitall(1, outRequests, outStatus);
    else
        MPI_Waitall(numDomains-1, inRequests+1,  inStatus+1);

//    INFO_Printf("[%d]: Slave_to_Master_Atoms--(c)\n",myDomain);
    /* send atom information */
    outBuf = NULL;  inBuf = NULL;
    if(myDomain!=0) /* slave */
    {
        outBuf = (double *) malloc(nsend*sizeof(double)*7); /* SR, VSR, ind */
        /* pack */
        n = 0;
        for(i=0;i<_NP;i++)
        {
            if(fixed[i]==-1) continue; /* -1 means not existing */
            if(domainID[i]==myDomain)
            {
                outBuf[n*7+0] = (double) i;
                outBuf[n*7+1] = _SR[i].x;
                outBuf[n*7+2] = _SR[i].y;
                outBuf[n*7+3] = _SR[i].z;
                outBuf[n*7+4] = _VSR[i].x;
                outBuf[n*7+5] = _VSR[i].y;
                outBuf[n*7+6] = _VSR[i].z;
                n++;
            }
        }
        MPI_Isend(outBuf,nsend*7,MPI_DOUBLE,0,MSG_ATOM,MPI_COMM_WORLD,&outRequests[0]);
    }
    else /* master */
    {
        inBuf = (double **) malloc(numDomains*sizeof(double *));
        for(idom=1;idom<numDomains;idom++)
        {
            inBuf[idom] = (double *) malloc(nrecv[idom]*sizeof(double)*7);
//            INFO_Printf("[%d]: Slave_to_Master_Atoms--(nrecv[%d]=%d)\n",myDomain,idom,nrecv[idom]);
        }
        for(idom=1;idom<numDomains;idom++)
            MPI_Irecv(inBuf[idom],nrecv[idom]*7,MPI_DOUBLE,idom,MSG_ATOM,MPI_COMM_WORLD,&inRequests[idom]);
    }
//    INFO_Printf("[%d]: Slave_to_Master_Atoms--(d)\n",myDomain);
    if(myDomain!=0)
    {
//        INFO_Printf("[%d]: Slave_to_Master_Atoms--(da nsend=%d)\n",myDomain,nsend);
        MPI_Waitall(1, outRequests, outStatus);
//        INFO_Printf("[%d]: Slave_to_Master_Atoms--(daa)\n",myDomain);
    }
    else
    {
//        INFO_Printf("[%d]: Slave_to_Master_Atoms--(db)\n",myDomain);
        MPI_Waitall(numDomains-1, inRequests+1,  inStatus+1);
//        INFO_Printf("[%d]: Slave_to_Master_Atoms--(dbb)\n",myDomain);
    }
//    INFO_Printf("[%d]: Slave_to_Master_Atoms--(dd)\n",myDomain);
    MPI_Barrier(MPI_COMM_WORLD);
    
    /* unpack */
    if(myDomain==0)
    {
        for(idom=1;idom<numDomains;idom++)
        {
            for(n=0;n<nrecv[idom];n++)
            {
                i   = (int) inBuf[idom][n*7+0];
                _SR[i].x  = inBuf[idom][n*7+1];
                _SR[i].y  = inBuf[idom][n*7+2];
                _SR[i].z  = inBuf[idom][n*7+3];
                _VSR[i].x = inBuf[idom][n*7+4];
                _VSR[i].y = inBuf[idom][n*7+5];
                _VSR[i].z = inBuf[idom][n*7+6];
            }
        }
    }
//    INFO_Printf("[%d]: Slave_to_Master_Atoms--(e)\n",myDomain);
    
    /* free allocated memory */
    if(myDomain!=0) /* slave */
    {
        free(outBuf);
    }
    else /* master */
    {
        for(i=1;i<numDomains;i++) free(inBuf[i]);
        free(inBuf);
        free(nrecv);
    }
//    INFO_Printf("[%d]: Slave_to_Master_Atoms--(done)\n",myDomain);   
}

void MDPARALLELFrame::Slave_chdir()
{
    char extname[200];

    MPI_Bcast(dirname,1000,MPI_CHAR,0,MPI_COMM_WORLD);
    if(myDomain!=0)
    {
        LFile::SubHomeDir(dirname, extname);
        if(chdir(extname)!=0)
        {
            ERROR("cd "<<dirname<<" failed");
        }
    }
}

void MDPARALLELFrame::potential_parallel()
{
#if 0
    if(numDomains>1)
    {
        Comm_Neighbor_Domains_Atoms();
        Mark_Local_Atoms();
    }
    potential();
    /* set forces and other properties to zero for atoms
     *  not belonging to this domain */

    if(numDomains>1)
    {
        _EPOT=0;
        for(i=0;i<_NP;i++)
        {
            if(domainID[i]!=myDomain)
            {
                _EPOT_IND[i]=0;
                _F[i].clear();
                _VIRIAL_IND[i].clear();
            }
        }
    }    
#else
   potential();
#endif
}


void MDPARALLELFrame::Master_Collect_Results()
{
    int i;
    
    MPI_Reduce(_EPOT_IND,_EPOT_IND_global,_NP,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(_F,_F_global,_NP*3,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(_VIRIAL_IND,_VIRIAL_IND_global,_NP*9,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

    memmove(_EPOT_IND,_EPOT_IND_global,_NP*sizeof(double));
    memmove(_F,_F_global,_NP*3*sizeof(double));
    memmove(_VIRIAL_IND,_VIRIAL_IND_global,_NP*9*sizeof(double));
    
    if(myDomain==0)
    {
        _EPOT=0; _VIRIAL.clear();
        for(i=0;i<_NP;i++)
        {
            _EPOT+=_EPOT_IND[i];
            _VIRIAL+=_VIRIAL_IND[i];
        }
    }
}
    
void MDPARALLELFrame::eval_parallel()
{
    /* master asking slaves to call the same function */
    //if(myDomain==0) Master_to_Slave("eval_parallel");

    Mark_Local_Atoms();

    potential_parallel();

    Master_Collect_Results();
    
    if(myDomain==0)
    {    
        calcprop();
        printResult();
    }
}

void MDPARALLELFrame::run_parallel()
{
    /* parallel Molecular Dynamics simulation
     * assuming Partition_Atoms, Broadcast_Atoms have been called
     * every CPU knows its domain boundary
     * every CPU has a copy of all atoms
     */
    int itmp;
    double pres, omega0;
    class Matrix33 prp, h0inv, h0invtran;
    int nspec, *spec;

    /* decide which algorithm to run */
    if (strcmp(ensemble_type,"NVE")==0)
        algorithm_id = 10000;
    else if (strcmp(ensemble_type,"NVT")==0)
        algorithm_id = 20000;
    else if (strcmp(ensemble_type,"NPH")==0)
        algorithm_id = 30000;
    else if (strcmp(ensemble_type,"NPT")==0)
        algorithm_id = 40000;
    else
    {
        algorithm_id = 100;
        ERROR("unknown ensemble_type("<<ensemble_type<<" use NVE\n"
            "choices: NVE, NVT, NPH, NPT");        
    }
    if (strcmp(integrator_type,"Gear6")==0)
        algorithm_id += 0;
    else if (strcmp(integrator_type,"VVerlet")==0)
        algorithm_id += 100;
    else
    {
        ERROR("unknown ensemble_type("<<ensemble_type<<" use Gear6\n"
              "choices: Gear6 VVerlet");
    }

    algorithm_id += implementation_type;
    INFO("algorithm_id = "<<algorithm_id);
    
    itmp = algorithm_id / 10000; /* 1: NVE, 2: NVT, 3: NPH, 4: NPT */

    /* master asking slaves to call the same function */
    //if(myDomain==0) Master_to_Slave("run_parallel");

    SHtoR(); /* all CPUs */
    Mark_Local_Atoms();

    /* prepare _SIGMA for Parrinello-Rahman */
    if((itmp==3)||(itmp==4)) /* NPH or NPT */
    {
        //call_potential(); /* do not need to call potential here */
        if(myDomain==0) /* master only */
        {
            pres=_EXTSTRESS.trace()/3;
            prp=_EXTSTRESS;
            prp[0][0]-=pres;
            prp[1][1]-=pres;
            prp[2][2]-=pres;
            prp*=_EXTSTRESSMUL;
            pres*=_EXTSTRESSMUL;
            pres+=_EXTPRESSADD;
            
            if(_H0.det()<1e-10)
                FATAL("you need to specify H0 by calling saveH "
                      "before using the Parrinello-Rahman method");
            
            h0inv=_H0.inv();
            h0invtran=h0inv.tran();
            omega0=_H0.det()*(1-_VACUUMRATIO);
            _SIGMA=((h0inv*prp)*h0invtran)*omega0;
            _SIGMA/=160.2e3;/* convert MPa to eV/A^3 */
        }
    }

    /* TO DO: write a function that broadcast all MD++ parameters */
    MPI_Bcast(&_SIGMA[0][0],9,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&_ATOMMASS,MAXSPECIES,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&vt2,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&_TIMESTEP,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(ensemble_type,100,MPI_CHAR,0,MPI_COMM_WORLD);
    MPI_Bcast(integrator_type,100,MPI_CHAR,0,MPI_COMM_WORLD);
    
//    INFO_Printf("[%d] run_parallel-(a)\n",myDomain);    
    nspec = 16;
    spec = (int *) malloc(nspec*sizeof(int));
    if(myDomain==0)
    {
        spec[0]  = totalsteps;
        spec[1]  = equilsteps;
        spec[2]  = 0; //usescalevelocity;
        spec[3]  = 0; //usenosehoover;
        spec[4]  = savecn;
        spec[5]  = savecnfreq;        
        spec[6]  = saveprop;
        spec[7]  = savepropfreq;
        spec[8]  = savecn;
        spec[9]  = printfreq;
        spec[10] = plotfreq;
        spec[11] = autowritegiffreq;
        spec[12] = conj_fixbox;
        spec[13] = algorithm_id;
        spec[14] = _NNM;
        spec[15] = _NIC;
    }
    MPI_Bcast(spec, nspec, MPI_INT,   0,MPI_COMM_WORLD);
    if(myDomain!=0)
    {
        totalsteps = spec[0];
        equilsteps = spec[1];
        //usescalevelocity = spec[2];
        //usenosehoover = spec[3];
        savecn = spec[4];
        savecnfreq = spec[5];
        saveprop = spec[6];
        savepropfreq = spec[7];
        savecn = spec[8];
        printfreq = spec[9];
        plotfreq = spec[10];
        autowritegiffreq = spec[11];
        conj_fixbox = spec[12];
        algorithm_id = spec[13];
        _NNM = spec[14];
        _NIC = spec[15];
    }
    free(spec);

//    INFO_Printf("[%d] run_parallel-(b) totalsteps = %d\n",myDomain,totalsteps);    
    /* main loop begins */
    for(curstep=0;curstep<totalsteps;curstep++) /* all CPUs */
    {
        while(win!=NULL) /* usually only master has a window open */
        {
            if(win->IsPaused()) sleep(1);
            else break;
        }

//        INFO_Printf("[%d] run_parallel-(c)\n",myDomain);    

        if(curstep%savepropfreq==0)
            INFO("["<<myDomain<<"]: curstep="<<curstep<<" EPOT="<<_EPOT);
        
        winplot();

//        INFO_Printf("[%d] run_parallel-(d)\n",myDomain);    

        if(curstep>=equilsteps)
        {
            if(saveprop)
            {
                if((curstep%savepropfreq)==0)
                {
                    Slave_to_Master_Atoms();/* update atoms to master */
                    //Broadcast_Atoms();
                    potential_parallel(); 
                    Master_Collect_Results();
                    if(myDomain==0) /* master only */
                    {
                        calcprop();
                        calcoutput();
                        pf.write(this);
                    }
                    potential_parallel();                    
                }
            }
            if(savecn)
                if((curstep%savecnfreq)==0&&(curstep!=0))
                {
                    Slave_to_Master_Atoms(); /* update atoms to master */
                    if(myDomain==0) /* master only */
                    {
                        intercn.write(this,zipfiles,true);
                    }
                }
        }

//        INFO_Printf("[%d] run_parallel-(e)\n",myDomain);    
        /* Integrator */
        step_parallel(); /* all CPUs */
        
//        INFO_Printf("[%d] run_parallel-(f)\n",myDomain);
        
        if(win!=NULL)
            if(win->IsAlive())
            if(autowritegiffreq>0)
                if((curstep%autowritegiffreq)==0)
                {
                    win->LockWritegif(); //lock and wait for writegif
                }
    }

    Slave_to_Master_Atoms();/* update atoms to master */
    if(saveprop)
    {
        //Broadcast_Atoms();
        potential_parallel(); /* require Broadcast */
        Master_Collect_Results();
        if(myDomain==0) /* master only */
        {
            calcprop(); 
            pf.write(this);
        }
    }
}

void MDPARALLELFrame::step_parallel()
{ /* this function is called within run_parallel()
   * so it does not require Master_to_Slave()
   */
    int i;

    /* decide which atoms are local atoms */
    //Mark_Local_Atoms(); /* other atoms fixed = -1 */

    /* do not let fixed or ghost atoms move */
    for(i=0;i<_NP;i++)
        if(domainID[i]!=myDomain)
        {
            _VSR[i].clear();
            _F[i].clear();
        }    

    /* integrator */
    switch(algorithm_id) {
    case(10000): NVE_Gear6();                     break;
    case(20000): NVT_Gear6();                     break;
    case(30000): NPH_Gear6();                     break;
    case(40000): NPT_Gear6();                     break;
    case(10100): NVE_VVerlet(curstep);            break;
    case(20100): NVT_VVerlet_Implicit(curstep);   break;
    case(20101): NVT_VVerlet_Explicit_1(curstep); break;
    case(20102): NVT_VVerlet_Explicit_2(curstep); break;
    case(30100): NPH_VVerlet(curstep);            break;
    case(40100): NPT_VVerlet_Implicit(curstep);   break;
    case(40101): NPT_VVerlet_Explicit_1(curstep); break;
    case(40102): NPT_VVerlet_Explicit_2(curstep); break;
    default: FATAL("unknown algorithm_id ("<<algorithm_id<<")"); break;
    }
  
}


void MDPARALLELFrame::Comm_Neighbor_Domains_Atoms()
{
    int i, n,  idom, yourDomain, yourIX, yourIY, yourIZ;
    int *nsend, *nrecv; double **inBuf, **outBuf;
    double dsx, dsy, dsz, sx0, sy0, sz0, sxw, syw, szw;
    Vector3 si, h;

    h=_H.height();
    dsx = fabs(1.4*_RLIST/h[0]);
    dsy = fabs(1.4*_RLIST/h[1]);
    dsz = fabs(1.4*_RLIST/h[2]);

    nsend = (int *) malloc(numNeighDoms*sizeof(int));
    nrecv = (int *) malloc(numNeighDoms*sizeof(int));
    inBuf = (double **) malloc(numNeighDoms*sizeof(double *));
    outBuf= (double **) malloc(numNeighDoms*sizeof(double *));
    
    /* decide which local atoms need to be sent to neighbors */
    for(idom=0;idom<numNeighDoms;idom++) /* all neighboring domains */
    {
        yourDomain = neighDoms[idom];
        MPI_Irecv(&nrecv[idom],1,MPI_INT,yourDomain,
                  MSG_SKIN_LEN,MPI_COMM_WORLD,&inRequests[idom]);
    }
    for(idom=0;idom<numNeighDoms;idom++) /* all neighboring domains */
    {
        yourDomain = neighDoms[idom];
        yourIZ = yourDomain % nZdoms;
        yourIY = ((yourDomain-myIZ)/nZdoms) % nYdoms;
        yourIX = ((yourDomain-myIZ)/nZdoms - myIY) / nYdoms;
        
        sx0 = (yourIX+0.5)/nXdoms-0.5; sxw = 0.5/nXdoms; /* center and width of each domain */
        sy0 = (yourIY+0.5)/nYdoms-0.5; syw = 0.5/nYdoms; 
        sz0 = (yourIZ+0.5)/nZdoms-0.5; szw = 0.5/nZdoms; 

        nsend[idom]=0;
        for(i=0;i<_NP;i++) /* first pass, count number of atoms */
        {
            if(fixed[i]==-1) continue; /* -1 means not existing */
            if(domainID[i]!=myDomain) continue; /* only examine local atoms */
            
            si=_SR[i];
            si.x-=sx0; si.y-=sy0; si.z-=sz0;
            si.subint();

            /* find atoms within the boundaries of neighboring domain */
            if((si.x >= -sxw-dsx)&&(si.x < sxw+dsx)&&
               (si.y >= -syw-dsy)&&(si.y < syw+dsy)&&
               (si.z >= -szw-dsz)&&(si.z < szw+dsz))
            { /* belong to outBuf */
                nsend[idom] ++;
            }
        }
//        INFO_Printf("[%d] nsend[%d]=%d\n",myDomain,idom,nsend[idom]);
        MPI_Isend(&nsend[idom],1,MPI_INT,yourDomain,
                  MSG_SKIN_LEN,MPI_COMM_WORLD,&outRequests[idom]);

        outBuf[idom] = (double *) malloc(nsend[idom]*7*sizeof(double));
        n = 0;
        for(i=0;i<_NP;i++) /* second pass, packing atoms */
        {
            if(fixed[i]==-1) continue; /* -1 means not existing */
            if(domainID[i]!=myDomain) continue; /* only examine local atoms */
            
            si=_SR[i]; 
            si.x-=sx0; si.y-=sy0; si.z-=sz0;
            si.subint();

            /* find atoms within the boundaries of neighboring domain */
            if((si.x >= -sxw-dsx)&&(si.x < sxw+dsx)&&
               (si.y >= -syw-dsy)&&(si.y < syw+dsy)&&
               (si.z >= -szw-dsz)&&(si.z < szw+dsz))
            { /* add to outBuf */
                outBuf[idom][n*7+0] = (double) i;
                outBuf[idom][n*7+1] = _SR[i].x;
                outBuf[idom][n*7+2] = _SR[i].y;
                outBuf[idom][n*7+3] = _SR[i].z;
                outBuf[idom][n*7+4] = _VSR[i].x;
                outBuf[idom][n*7+5] = _VSR[i].y;
                outBuf[idom][n*7+6] = _VSR[i].z;
                n++;
            }
        }
        if(n!=nsend[idom])
            FATAL("n="<<n<<" nsend["<<idom<<"]="<<nsend[idom]);
    }
    MPI_Waitall(numNeighDoms, outRequests,outStatus);
    MPI_Waitall(numNeighDoms, inRequests, inStatus );

//    for(idom=0;idom<numNeighDoms;idom++) /* all neighboring domains */
//        INFO_Printf("[%d] nrecv[%d]=%d\n",myDomain,idom,nrecv[idom]);
    
//    INFO_Printf("[%d]: Comm_Neighbor_Domains_Atoms-(a)\n",myDomain);
    /* sending the atoms in the skin layer */
    for(idom=0;idom<numNeighDoms;idom++) /* all neighboring domains */
    {
        yourDomain = neighDoms[idom];
        inBuf[idom] = (double *) malloc(nrecv[idom]*7*sizeof(double));
        MPI_Irecv(inBuf[idom],nrecv[idom]*7,MPI_DOUBLE,yourDomain,
                  MSG_SKIN,MPI_COMM_WORLD,&inRequests[idom]);
    }
//    INFO_Printf("[%d]: Comm_Neighbor_Domains_Atoms-(b)\n",myDomain);
    for(idom=0;idom<numNeighDoms;idom++) /* all neighboring domains */
    {
        yourDomain = neighDoms[idom];
//        INFO_Printf("[%d]: idom=%d yourDomain=%d\n",myDomain,idom,yourDomain);
        
        MPI_Isend(outBuf[idom],nsend[idom]*7,MPI_DOUBLE,yourDomain,
                  MSG_SKIN,MPI_COMM_WORLD,&outRequests[idom]);
    }
//    INFO_Printf("[%d]: Comm_Neighbor_Domains_Atoms-(bb)\n",myDomain);
    MPI_Waitall(numNeighDoms, outRequests,outStatus);
//    INFO_Printf("[%d]: Comm_Neighbor_Domains_Atoms-(bbb)\n",myDomain);
    MPI_Waitall(numNeighDoms, inRequests, inStatus );

//    INFO_Printf("[%d]: Comm_Neighbor_Domains_Atoms-(c)\n",myDomain);
    /* unpack */
    for(idom=0;idom<numNeighDoms;idom++) /* all neighboring domains */
    {
//        INFO_Printf("[%d]: Comm_Neighbor_Domains_Atoms-(idom=%d)\n",myDomain,idom);
        for(n=0;n<nrecv[idom];n++)
        {
            i   = (int) inBuf[idom][n*7+0];
            _SR[i].x  = inBuf[idom][n*7+1];
            _SR[i].y  = inBuf[idom][n*7+2];
            _SR[i].z  = inBuf[idom][n*7+3];
            _VSR[i].x = inBuf[idom][n*7+4];
            _VSR[i].y = inBuf[idom][n*7+5];
            _VSR[i].z = inBuf[idom][n*7+6];
        }
    }

    for(i=0;i<numNeighDoms;i++) free(inBuf [i]);
    for(i=0;i<numNeighDoms;i++) free(outBuf[i]);
    free(inBuf);
    free(outBuf);
    free(nsend);
    free(nrecv);
//    INFO_Printf("[%d]: Comm_Neighbor_Domains_Atoms-(end)\n",myDomain);
}

void MDPARALLELFrame::Mark_Local_Atoms()
{
    int i, ix, iy, iz, np0, np1;
    double dsx, dsy, dsz, sx0, sy0, sz0, sxw, syw, szw;
    Vector3 si, h;
    
    /* master asking slaves to call the same function */
    //if(myDomain==0) Master_to_Slave("Mark_Local_Atoms");

    h=_H.height();
    dsx = fabs(1.4*_RLIST/h[0]); /* factor 1.08889 from MEAM-LAMMPS */
    dsy = fabs(1.4*_RLIST/h[1]);
    dsz = fabs(1.4*_RLIST/h[2]);

    sx0 = (myIX+0.5)/nXdoms-0.5; sxw = 0.5/nXdoms; /* center and width of each domain */
    sy0 = (myIY+0.5)/nYdoms-0.5; syw = 0.5/nYdoms; 
    sz0 = (myIZ+0.5)/nZdoms-0.5; szw = 0.5/nZdoms; 

    np0=0; np1=0;
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]==-1) fixed[i]=0;
        si=_SR[i]; si.subint();
        ix = (int) floor((si.x+0.5)*nXdoms);
        iy = (int) floor((si.y+0.5)*nYdoms);
        iz = (int) floor((si.z+0.5)*nZdoms);

        /* decide which domain atom i belongs to */
        domainID[i] = (ix*nYdoms + iy)*nZdoms + iz;

        /* mark atoms beyond the boundaries of current domain
         * as "non-existing", i.e. fixed[i] = -1 */
        si.x-=sx0; si.y-=sy0; si.z-=sz0;
        si.subint();
        if((si.x < -sxw-dsx)||(si.x >= sxw+dsx)||
           (si.y < -syw-dsy)||(si.y >= syw+dsy)||
           (si.z < -szw-dsz)||(si.z >= szw+dsz))
        {
            fixed[i] = -1;            
            if(domainID[i]==myDomain)
            {
                INFO_Printf("[%d]: atom i (%e %e %e) si (%e %e %e) fixed = %d\n",
                            myDomain,_SR[i].x,_SR[i].y,_SR[i].z,si.x,si.y,si.z,fixed[i]);
                FATAL("core atoms outside skin!");
            }                
        }
        if(domainID[i]==myDomain) np0++;
        if(fixed[i]!=-1) np1++;
        
    }
    INFO_Printf("[%d]: number of core atoms = %d core + skin atoms = %d\n",
                myDomain,np0,np1);
}

void MDPARALLELFrame::Broadcast_nebsetting()
{
    /* Broadcast constrainatoms setting */
    MPI_Bcast(constrainatoms, MAXCONSTRAINATOMS, MPI_INT,   0,MPI_COMM_WORLD);
    MPI_Bcast(nebspec,        NEBSPECSIZE,       MPI_INT,   0,MPI_COMM_WORLD);
    MPI_Bcast(&_TIMESTEP,     1,                 MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&_CHAINLENGTH,  1,                 MPI_INT,   0,MPI_COMM_WORLD);
    MPI_Bcast(&totalsteps,    1,                 MPI_INT,   0,MPI_COMM_WORLD);
    MPI_Bcast(&curstep,       1,                 MPI_INT,   0,MPI_COMM_WORLD);
    MPI_Bcast(&equilsteps,    1,                 MPI_INT,   0,MPI_COMM_WORLD);
    MPI_Bcast(&printfreq,     1,                 MPI_INT,   0,MPI_COMM_WORLD);   
    MPI_Bcast(&_VACUUMRATIO,  1,                 MPI_DOUBLE,0,MPI_COMM_WORLD);

    MPI_Bcast(&_NIMAGES,      1,                 MPI_INT,   0,MPI_COMM_WORLD);   
    MPI_Bcast(&_TORSIONSIM,   1,                 MPI_INT,   0,MPI_COMM_WORLD);   
    MPI_Bcast(&_BENDSIM,      1,                 MPI_INT,   0,MPI_COMM_WORLD);   
    MPI_Bcast(&_ENABLE_SWITCH,1,                 MPI_INT,   0,MPI_COMM_WORLD);   
    MPI_Bcast(&_constrainedMD,1,                 MPI_INT,   0,MPI_COMM_WORLD);   
    MPI_Bcast(&_ENABLE_FEXT,  1,                 MPI_INT,   0,MPI_COMM_WORLD);   
    MPI_Bcast(&extforce,      EXTFORCESIZE,      MPI_DOUBLE,0,MPI_COMM_WORLD);   
    MPI_Bcast(&forcemul,      1,                 MPI_DOUBLE,0,MPI_COMM_WORLD);   

    MPI_Bcast(&conj_ftol,     1,                 MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&conj_itmax,    1,                 MPI_INT,   0,MPI_COMM_WORLD);
    MPI_Bcast(&conj_fevalmax, 1,                 MPI_INT,   0,MPI_COMM_WORLD);
    MPI_Bcast(&conj_fixbox,   1,                 MPI_INT,   0,MPI_COMM_WORLD);
    MPI_Bcast(energythreshold,ENERGYTHRESHOLDSIZE,MPI_DOUBLE,0,MPI_COMM_WORLD);
}

void MDPARALLELFrame::AllocChain_parallel()
{
    int n, size;
    
    INFO_Printf("AllocChain_parallel[%d]: \n", myDomain);

    Broadcast_nebsetting();

    /* Allocate chain */
    n = constrainatoms[0];
    size = sizeof(Vector3)*n;
    if (_Rc0==NULL) _Rc0 = (Vector3 *) malloc(size);  /* own   copy */
    if (_Rc1==NULL) _Rc1 = (Vector3 *) malloc(size);  /* left  copy */
    if (_Rc2==NULL) _Rc2 = (Vector3 *) malloc(size);  /* right copy */
    if (_Fc0==NULL) _Fc0 = (Vector3 *) malloc(size);  /* force on own copy */
    if (_Tan==NULL) _Tan = (Vector3 *) malloc(size);  /* local tangential vector */
    
    memset(_Rc0,0,size);
    memset(_Rc1,0,size);
    memset(_Rc2,0,size);
    memset(_Fc0,0,size);
    memset(_Tan,0,size);

    Realloc(_nebinterp,double,(_CHAINLENGTH+1));
    memset(_nebinterp,0,sizeof(double)*(_CHAINLENGTH+1));
}


void MDPARALLELFrame::initRchain_parallel()
{
    int i, j, n, ipt;
    double s;
    Vector3 ds;
    
    /* Broadcast atom configuration */
    Broadcast_Atoms();
    
    SHtoR(); /* all CPUs */
        
    /* Allocate _SR1 and _SR2 (start and end states) */
    if(myDomain==0)
    {
        if(_SR1 == NULL || _SR2 == NULL)
        {
            if(_SR1 == NULL)
                ERROR("config1 is not set!");
            else
                ERROR("config2 is not set!");
            return;
        }

    }
    else 
    { /* Slaves allocate _SR1, _SR2 */
        Free(_SR1); Realloc(_SR1,Vector3,_NP);
        Free(_SR2); Realloc(_SR2,Vector3,_NP);
    }        

    /* Broadcast _SR1 and _SR2 */
    MPI_Bcast(_SR1, _NP*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(_SR2, _NP*3,MPI_DOUBLE,0,MPI_COMM_WORLD);

    /* Calculate dscom12 before interpCN(s) */
    caldscom12();
    INFO_Printf("dscom12[%d] = (%f %f %f)\n",myDomain,dscom12.x,dscom12.y,dscom12.z);
    
    /* Initialize straight path */
    for(j=0;j<=_CHAINLENGTH;j++)
        _nebinterp[j] = (1.0*j)/_CHAINLENGTH;
    
    s = _nebinterp[myDomain];
    n = constrainatoms[0];
    for(i=0;i<n;i++) /* compute _Rc0 */
    {
        ipt=constrainatoms[i+1];
        ds=_SR2[ipt]-dscom12; ds-=_SR1[ipt]; ds.subint(); /* 2007/3/6 Wei Cai */
        ds*=s; ds+=_SR1[ipt];
        _Rc0[i]=_H*ds;
    }

}

void MDPARALLELFrame::allocChainVars_parallel(double **pEc, double **pEc_global, double **pFm, double **pFm_global, 
           double **pFt, double **pFt_global, double **pFs, double **pFs_global, double **pdR, double **pdR_global)
{
    /* Allocate workspace variables for nebrelax
     * Ec:       potential energy of each replica
     * Fm:       magnitude of the orthogonalized force (perpendicular to local tangent)
     * Ft:       dot product between atomic forces and the local tangent vector
     * Fs:       magnitude of spring force (along the local tangent)
     * dR:       segment length (along path) in Angstrom dR_j = |R_j - R_{j-1}|
     */
    *pEc       =(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* energy along path */
    *pEc_global=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* for communication */
    *pFm       =(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* force normal to path */
    *pFm_global=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* for communication */
    *pFt       =(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* force along path */
    *pFt_global=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* for communication */
    *pFs       =(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* spring force magnitude */
    *pFs_global=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* for communication */
    *pdR       =(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* segment length in Angstrom */
    *pdR_global=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* for communication */

    memset(*pEc,       0,sizeof(double)*(_CHAINLENGTH+1));
    memset(*pEc_global,0,sizeof(double)*(_CHAINLENGTH+1));
    memset(*pFm,       0,sizeof(double)*(_CHAINLENGTH+1));
    memset(*pFm_global,0,sizeof(double)*(_CHAINLENGTH+1));
    memset(*pFt,       0,sizeof(double)*(_CHAINLENGTH+1));
    memset(*pFt_global,0,sizeof(double)*(_CHAINLENGTH+1));
    memset(*pFs,       0,sizeof(double)*(_CHAINLENGTH+1));
    memset(*pFs_global,0,sizeof(double)*(_CHAINLENGTH+1));
    memset(*pdR,       0,sizeof(double)*(_CHAINLENGTH+1));
    memset(*pdR_global,0,sizeof(double)*(_CHAINLENGTH+1));
}

void MDPARALLELFrame::allocChainIndex_parallel(int **pleft_index,  int **pleft_index_global, int **pright_index, int **pright_index_global)
{
    /* Allocate workspace variables for stringrelax
     */
    *pleft_index        =(int *)malloc(sizeof(int)*(_CHAINLENGTH+1)); 
    *pleft_index_global =(int *)malloc(sizeof(int)*(_CHAINLENGTH+1)); 
    *pright_index       =(int *)malloc(sizeof(int)*(_CHAINLENGTH+1)); 
    *pright_index_global=(int *)malloc(sizeof(int)*(_CHAINLENGTH+1)); 

    memset(*pleft_index,        0,sizeof(int)*(_CHAINLENGTH+1));
    memset(*pleft_index_global, 0,sizeof(int)*(_CHAINLENGTH+1));
    memset(*pright_index,       0,sizeof(int)*(_CHAINLENGTH+1));
    memset(*pright_index_global,0,sizeof(int)*(_CHAINLENGTH+1));
}

void MDPARALLELFrame::calChainForce_parallel(int n, int relax_surround, double *Ec, double *Ec_global)
{
    int i, ipt;
    double s;
    Matrix33 hinv;

    /* Get atomic forces to _Fc array along the entire chain */
    s = _nebinterp[myDomain];
    hinv = _H.inv();

    /* get forces, s is reaction coordinate */
    interpCN(s); /* interpolate all atoms for a given s */

    /* copy _Rc0 to _SR */
    for(i=0;i<n;i++)
    {
        ipt=constrainatoms[i+1];
        _SR[ipt]=hinv*_Rc0[i];
    }

    /* relax surrounding atoms */
    if(relax_surround==1) relax(); // NOT WORKING !!

    /* enforces the neighborlist to be reconstructed at the very first step (important!!) */
    /* clear Verlet neighbor list => enforce neighbor list updated at every step */
    if(curstep==step0) MDFrame::clearR0();
    MDFrame::call_potential();

    /* Put potential energies into Ec array */
    memset(Ec,0,sizeof(double)*(_CHAINLENGTH+1));
    Ec[myDomain]=_EPOT;
    MPI_Allreduce(Ec,Ec_global,(_CHAINLENGTH+1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    memcpy(Ec,Ec_global,sizeof(double)*(_CHAINLENGTH+1));

    /* Put atomic forces to _Fc0 array */
    for(i=0;i<n;i++)
    {
        ipt=constrainatoms[i+1];
        _Fc0[i]=_F[ipt];
    }
}

void MDPARALLELFrame::commNeighborRc(int n)
{
    int nin, nout;
    /* receive _Rc1, _Rc2 from neighbors */
    nin = nout = 0;
    if(myDomain>0) /* receive _Rc1 from left node, send _Rc0 to left node */
    {
        MPI_Irecv(_Rc1,n*3,MPI_DOUBLE,myDomain-1,
                  MSG_NEB_ATOM,MPI_COMM_WORLD,&inRequests[nin]);
        nin ++;
        MPI_Isend(_Rc0,n*3,MPI_DOUBLE,myDomain-1,
                  MSG_NEB_ATOM,MPI_COMM_WORLD,&outRequests[nout]);
        nout ++;
    }
    if(myDomain<_CHAINLENGTH) /* receive _Rc2 from right node, send _Rc0 to right node */
    {
        MPI_Irecv(_Rc2,n*3,MPI_DOUBLE,myDomain+1,
                  MSG_NEB_ATOM,MPI_COMM_WORLD,&inRequests[nin]);
        nin ++;
        MPI_Isend(_Rc0,n*3,MPI_DOUBLE,myDomain+1,
                  MSG_NEB_ATOM,MPI_COMM_WORLD,&outRequests[nout]);
        nout ++;
    }
    MPI_Waitall(nout, outRequests,outStatus);
    MPI_Waitall(nin,  inRequests, inStatus );
}
    

void MDPARALLELFrame::calChainTang_parallel(int n, double *Ec, double *dR, double *dR_global, double *pTanMag2)
{
    int i;
    double dEmax, dEmin, r2;
    Vector3 dt;

    /* calculate tangential vector, _Tan */

    /* receive _Rc1, _Rc2 from neighbors */
    commNeighborRc(n);

    /* calculate tangential vector */
    if((myDomain>0)&&(myDomain<_CHAINLENGTH))
    {
        if (((Ec[myDomain+1]-Ec[myDomain]>0) && (Ec[myDomain-1]-Ec[myDomain]>0)) ||
            ((Ec[myDomain+1]-Ec[myDomain]<0) && (Ec[myDomain-1]-Ec[myDomain]<0))) // convex or concave profile
        {
            dEmax = max(fabs(Ec[myDomain+1]-Ec[myDomain]),fabs(Ec[myDomain-1]-Ec[myDomain]));
            dEmin = min(fabs(Ec[myDomain+1]-Ec[myDomain]),fabs(Ec[myDomain-1]-Ec[myDomain]));
            
            if (Ec[myDomain+1]>Ec[myDomain-1])
            {
                for(i=0;i<n;i++) // loop through all constrained atoms
                {
                    _Tan[i] = (_Rc2[i]-_Rc0[i])*dEmax + (_Rc0[i]-_Rc1[i])*dEmin;
                }
            }
            else
            {
                for(i=0;i<n;i++) // loop through all constrained atoms
                {
                    _Tan[i] = (_Rc2[i]-_Rc0[i])*dEmin + (_Rc0[i]-_Rc1[i])*dEmax;
                }
            }
        }
        else
        {
            if (Ec[myDomain+1]>Ec[myDomain-1]) // uphill
            {
                for(i=0;i<n;i++) // loop through all constrained atoms
                {
                    _Tan[i] = _Rc2[i]-_Rc0[i];
                }
            }
            else // downhill
            {
                for(i=0;i<n;i++) // n: number of constrained atoms
                {
                    _Tan[i] = _Rc0[i]-_Rc1[i];
                }
            }
        }
        
        *pTanMag2 = 0;        
        for(i=0;i<n;i++) // compute magnitude squared of tangential vector
            *pTanMag2 += _Tan[i].norm2();
    }
    
    /* calculate tangential vector at the left end */
    if(myDomain==0)
    {
        *pTanMag2 = 0;
        for(i=0;i<n;i++)
        {
            _Tan[i]= _Rc2[i]-_Rc0[i];
            *pTanMag2 += _Tan[i].norm2();
        }
    }   
    /* calculate tangential vector at the right end */
    if(myDomain==_CHAINLENGTH)
    {
        *pTanMag2 = 0;
        for(i=0;i<n;i++)
        {
            _Tan[i]= _Rc0[i]-_Rc1[i];
            *pTanMag2 += _Tan[i].norm2();
        }
    }

    /* calculate chain segment length dR[j] = | Rc[j] - Rc[j-1] | */
    memset(dR,0,sizeof(double)*(_CHAINLENGTH+1));
    if((myDomain>=1)&&(myDomain<=_CHAINLENGTH))
    {
        r2 = 0;
        for(i=0;i<n;i++) 
        { 
            dt=_Rc0[i]-_Rc1[i]; 
            r2+=dt.norm2(); 
        }
        dR[myDomain] = sqrt(r2);
    }
    MPI_Allreduce(dR,dR_global,(_CHAINLENGTH+1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    memcpy(dR,dR_global,sizeof(double)*(_CHAINLENGTH+1));

}

void MDPARALLELFrame::orthoChainForce_parallel(int moveleftend, int moverightend, int yesclimbimage, int EmaxDomain, double TanMag2, int n, double *Fm, double *Fm_global, double *Ft, double *Ft_global)
{
    int i;
    double fr, fm2;
    Vector3 dt;

    /* orthogonalize forces _Fc against local tangent vector _Tan */
    memset(Fm,0,sizeof(double)*(_CHAINLENGTH+1));
    memset(Ft,0,sizeof(double)*(_CHAINLENGTH+1));

    if((myDomain>=0)&&(myDomain<=_CHAINLENGTH))
    {
        fr=0;            
        for(i=0;i<n;i++)
            fr+=dot(_Fc0[i],_Tan[i]);
        Ft[myDomain] = fr/sqrt(TanMag2);
        fr /= TanMag2; // normalizing
        
        fm2=0;
        for(i=0;i<n;i++) /* orthogonalizing force against tangential vector */
        {
            dt=_Tan[i]*fr;   // tangential force vector

            if ((myDomain==0 && moveleftend==3) || (myDomain==_CHAINLENGTH && moverightend==3))
            {
               /* ends move in free fall, do not orthogonalize forces */
            } else if ((myDomain==0 && moveleftend==2) || (myDomain==_CHAINLENGTH && moverightend==2))
            {
               /* set force along local tangent direction */
               _Fc0[i] = dt;
            } else {
               /* orthogonalize forces */
               _Fc0[i] -= dt;
            }
            fm2 += _Fc0[i].norm2();

            if (yesclimbimage && (curstep>=equilsteps) && myDomain==EmaxDomain)
               _Fc0[i]-=dt;   // revert the force along tangent direction for climbing image
        }
        Fm[myDomain] = sqrt(fm2); /* store residual force magnitude */
    }
    MPI_Allreduce(Fm,Fm_global,(_CHAINLENGTH+1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    memcpy(Fm,Fm_global,sizeof(double)*(_CHAINLENGTH+1));
    MPI_Allreduce(Ft,Ft_global,(_CHAINLENGTH+1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    memcpy(Ft,Ft_global,sizeof(double)*(_CHAINLENGTH+1));

}

void MDPARALLELFrame::reparamChain_parallel(int moveleftend, int moverightend, int yesclimbimage, int EmaxDomain, int n, double *plavg, double *plavg0, double *dR, double *Ec, int *left_index, int *left_index_global, int *right_index, int *right_index_global, double energySlope, int manualcut, int manualleftend, int manualrightend)
{
    int i, jj, leftendDomain, rightendDomain, lowleft, lowright, nin, nout;
    double Ltot, new_position, alpha;
    Vector3 dt;

    /* assume that dR has already been calculated immediately before calling this function 
     *  e.g. by 
     *          calChainTang_parallel(n, Ec, dR, dR_global, &TanMag2);
     */

    leftendDomain = 0; rightendDomain = _CHAINLENGTH;

    if(EmaxDomain==0) WARNING("Max-energy domain = 0");
    if(EmaxDomain==_CHAINLENGTH) WARNING("Max-energy domain = _CHAINLENGTH");

    if (manualcut)
    {
       leftendDomain = manualleftend;
       rightendDomain = manualrightend;
    }
    else
    {
        if (moveleftend) 
        {
            /* Get the lowest energy configurations to the left of the 
             * maximum energy configuration. */
            lowleft = 0; 
            for(jj=0;jj<=EmaxDomain;jj++) {
                /* prevent very shallow region near state A */
                if(Ec[jj]<Ec[lowleft] + energySlope ) lowleft = jj; /* Wei 2/8/2014 */
            }
            /* However, we don't necessarily want to redistribute between the two lowest energy
             * configurations. We want to try to keep the left and right ends near the same energy. */
            leftendDomain = lowleft; 
            while(Ec[leftendDomain+1]-Ec[0] < 0) leftendDomain++;
        }
        if (moverightend) 
        {
            /* Get the lowest energy configurations to the left of the 
             * maximum energy configuration. */
           lowright = _CHAINLENGTH;
           for(jj=EmaxDomain;jj<=_CHAINLENGTH;jj++) {
               if(Ec[jj]<Ec[lowright]) lowright = jj;
           }
            /* However, we don't necessarily want to redistribute between the two lowest energy
             * configurations. We want to try to keep the left and right ends near the same energy. */
           rightendDomain = lowright;
           while(Ec[rightendDomain-1]-Ec[0] < 0) rightendDomain--;
        }
    }

    if (yesclimbimage)
    {
        ERROR("reparamChain_parallel: not compatible with yesclimbimage!");
    }

    /* calculate the total length from leftendDomain to rightendDomain */
    Ltot = 0.0;
    for(jj=leftendDomain+1;jj<=rightendDomain;jj++) Ltot += dR[jj];

    /* find left_index, right_index, alpha, for this replica myDomain */
    memset(left_index,        0,sizeof(int)*(_CHAINLENGTH+1));
    memset(left_index_global, 0,sizeof(int)*(_CHAINLENGTH+1));
    memset(right_index,       0,sizeof(int)*(_CHAINLENGTH+1));
    memset(right_index_global,0,sizeof(int)*(_CHAINLENGTH+1));

    new_position = (1.0*myDomain*Ltot)/_CHAINLENGTH;
    left_index[myDomain] = leftendDomain;
    while(new_position > dR[left_index[myDomain]+1] && left_index[myDomain] < _CHAINLENGTH) {
        left_index[myDomain]++;
        new_position -= dR[left_index[myDomain]];
    }
    if (left_index[myDomain] == _CHAINLENGTH) {
        left_index[myDomain] = _CHAINLENGTH-1; 
        right_index[myDomain] = _CHAINLENGTH;
        alpha = 1;
    } else { 
        right_index[myDomain] = left_index[myDomain] + 1; 
        alpha = new_position/dR[right_index[myDomain]];
    }
    MPI_Allreduce(left_index,left_index_global,(_CHAINLENGTH+1),MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    memcpy(left_index,left_index_global,sizeof(int)*(_CHAINLENGTH+1));
    MPI_Allreduce(right_index,right_index_global,(_CHAINLENGTH+1),MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    memcpy(right_index,right_index_global,sizeof(int)*(_CHAINLENGTH+1));

    /* Receive _Rc1 and _Rc2 from domains left_index and right_index */
    nin = nout = 0;
    /* receive _Rc1 from left index, send _Rc0 to those who need it as a left index */
    MPI_Irecv(_Rc1,n*3,MPI_DOUBLE,left_index[myDomain], MSG_NEB_ATOM,MPI_COMM_WORLD,&inRequests[nin]);
    nin ++;
    for(jj = 0; jj <= _CHAINLENGTH; jj++) {
	if(left_index[jj] == myDomain) { 
	    MPI_Isend(_Rc0,n*3,MPI_DOUBLE,jj, MSG_NEB_ATOM,MPI_COMM_WORLD,&outRequests[nout]);
	    nout ++;
	}
    }
    /* receive _Rc2 from right index, send _Rc0 to those who need it as a right index */
    MPI_Irecv(_Rc2,n*3,MPI_DOUBLE,right_index[myDomain], MSG_NEB_ATOM,MPI_COMM_WORLD,&inRequests[nin]);
    nin ++;
    for(jj = 0; jj <= _CHAINLENGTH; jj++) {
	if(right_index[jj] == myDomain) { 
	    MPI_Isend(_Rc0,n*3,MPI_DOUBLE,jj, MSG_NEB_ATOM,MPI_COMM_WORLD,&outRequests[nout]);
	    nout ++;
	}
    }
    MPI_Waitall(nout, outRequests,outStatus);
    MPI_Waitall(nin,  inRequests, inStatus );

    /* position _Rc0 between _Rc1 and _Rc2 */
    for(i=0;i<n;i++) { 
        _Rc0[i] = _Rc1[i]*(1.0-alpha)+_Rc2[i]*alpha;
    }

    /* update new positions of the neighbor chains */
    commNeighborRc(n);
}

void MDPARALLELFrame::nebrelax_parallel()
{/*
  * Parallel implementation of nebrelax()
  *
  * New Algorithm described in
  *
  * Improved tangent estimate in the nudged elastic band method for finding
  * minimum energy paths and saddle points, Henkelman and Johnsson,
  * Journal of Chemical Physics Vol.113 9978-9985, 2000
  *
  * relax _Rc0 by nudged elastic band method (steepest descent)
  *
  * Notes:  Assumes Broadcast_Atoms() has been called
  *         Every CPU has a copy of all atoms
  */
    int n, i, ipt, jj;
    int relax_surround, moveleftend, moverightend, yesclimbimage, EmaxDomain;
    Vector3 ds0;
    double minFm, maxFm, Fspring, spring_const, springK, TanMag2, E0, Emax;
    double *Ec, *Ec_global, *Fm, *Fm_global, *Ft, *Ft_global, *Fs, *Fs_global, *dR, *dR_global;
    FILE *fp;
    
    INFO_Printf("NEB Relax Parallel [%d]\n", myDomain);

    if( numDomains != (_CHAINLENGTH+1) )
    {
        ERROR("number of cpus("<<numDomains<<") does not match chainlength("<<_CHAINLENGTH<<")+1");
        return;
    }
    
    relax_surround= nebspec[0];
    spring_const  = (double) nebspec[1]; /* suggested spring constant, if <= 0 then automatically adjusted */
    moveleftend   = nebspec[2];
    moverightend  = nebspec[3];
    yesclimbimage = nebspec[4];

    /* if moveleftend = 0, fix the left end replica
                      = 1, move the left end replica by the spring force in the normal to the real force direction
                           (PNAS vol.104 3031-3036, 2007)
                      = 2, move the left end replica along the chain tangential direction (constrained free fall)
                      = 3, move the left end replica along the real force direction (free fall) */

    if (yesclimbimage && myDomain == 0)
        INFO_Printf("NEB: The climbing image method will be applied after %d steps.\n",equilsteps);

    /* Makes sure _SR1 and _SR2 (for configs A and B) are allocated */
    if(_SR1 == NULL || _SR2 == NULL)
    {
        if(_SR1 == NULL) ERROR("config1 is not set!");
        if(_SR2 == NULL) ERROR("config2 is not set!");
        return;
    }

    n = constrainatoms[0];
   
    fp = NULL; 
    if(myDomain==0) /* Master opens file */
    {
        fp=fopen("nebeng.out","w");
    }

    /* Allocate workspace variables for nebrelax */
    Ec = Ec_global = NULL; Fm = Fm_global = NULL; Ft = Ft_global = NULL; Fs = Fs_global = NULL; dR = dR_global = NULL; 
    allocChainVars_parallel(&Ec, &Ec_global, &Fm, &Fm_global, &Ft, &Ft_global, &Fs, &Fs_global, &dR, &dR_global);

    /* calculate center-of-mass shift from _SR1 (state A) to _SR2 (state B) */
    caldscom12();
    /* compute current constraint value (not really necessary) */
    constrain_dist=0;
    for(i=0;i<n;i++)
    {
        ipt=constrainatoms[i+1];
        ds0=(_SR2[ipt]-dscom12)-_SR1[ipt];
        constrain_dist+=ds0.norm2();
    }
    INFO_Printf("NEB[%d]: constrain_dist=%e s=%e\n", myDomain, constrain_dist, _nebinterp[myDomain]);

    if(relax_surround==1)
    {/* fix constrained atoms for relaxing surrounding atoms */
        for(i=0;i<n;i++)
        {
            ipt=constrainatoms[i+1];
            fixed[ipt]=1;
        }
    }

    /* start NEB Relax iteration */
    step0 = curstep;
    for(curstep=step0;curstep<=(step0 + totalsteps);curstep++)    
    {
        /* Get atomic forces to _Fc array along the entire chain */
        calChainForce_parallel(n, relax_surround, Ec, Ec_global);

        E0 = Ec[0];
        /* Find max-energy replica (i.e. the saddle point) */
        EmaxDomain = 0; Emax=E0;
        for(jj=1;jj<=_CHAINLENGTH;jj++)
        {
            if(Ec[jj]>Emax) { Emax=Ec[jj]; EmaxDomain=jj; }
        }

        /* calculate tangential vector, _Tan */
        calChainTang_parallel(n, Ec, dR, dR_global, &TanMag2);
         
        /* orthogonalize forces _Fc0 against local tangent vector _Tan */
        orthoChainForce_parallel(moveleftend, moverightend, yesclimbimage, EmaxDomain, 
                                 TanMag2, n, Fm, Fm_global, Ft, Ft_global);

        /* adjust spring constant if the suggested value is 0 or negative*/
        if (spring_const>0.0)
        {
            springK = spring_const;
            if (curstep==step0 && myDomain==0) INFO_Printf("use manually set spring constant k = %20.12e\n", springK);
        }
        else
        {
            minFm = maxFm = Fm[0]; // Min. and Max. magnitude of orthogonal force over jj
            for(jj=1;jj<=_CHAINLENGTH;jj++)
            {
                if (Fm[jj]<minFm) minFm = Fm[jj];
                if (Fm[jj]>maxFm) maxFm = Fm[jj];
            }
            springK = maxFm*_CHAINLENGTH;
            if (myDomain==0) INFO_Printf("automatically chosen spring constant k = %20.12e\n", springK);
        }

        /* add tangential component of spring force (not including left and right ends) */
        memset(Fs,0,sizeof(double)*(_CHAINLENGTH+1));
        if((myDomain>0)&&(myDomain<_CHAINLENGTH))
        {
            Fspring = springK*(dR[myDomain+1] - dR[myDomain]);

            Fs[myDomain] = Fspring;
            
            /* If climbing image, skip */
            if (yesclimbimage && (curstep>=equilsteps) && myDomain==EmaxDomain) continue;

            for(i=0;i<n;i++)
                _Fc0[i] += _Tan[i] * (Fspring/sqrt(TanMag2));
        }
        MPI_Allreduce(Fs,Fs_global,(_CHAINLENGTH+1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        memcpy(Fs,Fs_global,sizeof(double)*(_CHAINLENGTH+1));
        
        if(myDomain==0) /* Master print file */
        {
            if(curstep%printfreq==0)  // Mar. 08 2007 Keonwook Kang 
            {
                INFO_Printf("curstep = %d\n",curstep);
                for(jj=0;jj<=_CHAINLENGTH;jj++)
                    INFO_Printf("%20.12e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
                                _nebinterp[jj],Ec[jj]-Ec[0], Fm[jj], Ft[jj], Fs[jj], dR[jj]);
                fprintf(fp,"%d ",curstep);
                for(jj=0;jj<=_CHAINLENGTH;jj++)
                    fprintf(fp,"%20.12e %25.15e %25.15e %25.15e %25.15e %25.15e ",
                               _nebinterp[jj], Ec[jj]-Ec[0], Fm[jj], Ft[jj], Fs[jj], dR[jj]);
                fprintf(fp,"\n");
                fflush(fp); 
            }
        }

        if(curstep<(step0 + totalsteps))
        {/* move NEB chains along force */
            if ( (myDomain==0) && (moveleftend<=0) ) 
            {
                /* do nothing */
            }
            else if ( (myDomain==_CHAINLENGTH) && (moverightend<=0) ) 
            {
                /* do nothing */
            }
            else 
            {
                for(i=0;i<n;i++) _Rc0[i]+=_Fc0[i]*_TIMESTEP;
            }
        }
    }

    free(Ec); free(Ec_global); free(Fm); free(Fm_global); free(Ft); free(Ft_global); 
    free(Fs); free(Fs_global); free(dR); free(dR_global);

    if(myDomain==0) fclose(fp);

    INFO_Printf("NEB[%d]: exit\n",myDomain);
}

void MDPARALLELFrame::stringrelax_parallel()
{/*
  * Parallel implementation of stringrelax()
  * 
  * relax _Rc0 by string method 
  *
  * Notes:  stringrelax_parallel can use conjugate gradient method to relax _Rc0
  *         while stringrelax can only use steepest descent method
  *
  * Notes:  Assumes Broadcast_Atoms() has been called
  *         Every CPU has a copy of all atoms
  */
    int n, i, ipt, jj;
    int relax_surround, redistri_freq, moveleftend, moverightend, yesclimbimage, EmaxDomain;
    Vector3 ds0; Matrix33 hinv;
    double TanMag2, lavg0, lavg, E0, Emax;
    double *Ec, *Ec_global, *Fm, *Fm_global, *Ft, *Ft_global, *Fs, *Fs_global, *dR, *dR_global;
    int *left_index, *left_index_global, *right_index, *right_index_global;
    int rightendEnergyTooHigh, myEnergyTooHigh, manualcut, manualleftend, manualrightend, cgstepsmid,cgstepsrightend;
    double energySpike, energySlope;
    FILE *fp;
    
    INFO_Printf("stringrelax[%d]: String Relax Parallel\n", myDomain);

    if( numDomains != (_CHAINLENGTH+1) )
    {
        ERROR("number of cpus("<<numDomains<<") does not match chainlength("<<_CHAINLENGTH<<")+1");
        return;
    }
    
    relax_surround= nebspec[0];
    redistri_freq = nebspec[1]; 
    moveleftend   = nebspec[2];
    moverightend  = nebspec[3];
    yesclimbimage = nebspec[4];
    /* islengthconstant was nebspec[5]; deprecated */
    /* cgrelaxsteps was nebspec[6]; deprecated */
    manualcut     = nebspec[7];
    manualleftend = nebspec[8]; 
    manualrightend= nebspec[9]; 
    cgstepsmid    =nebspec[10];
    cgstepsrightend=nebspec[11];

    energySpike = energythreshold[0];
    energySlope = energythreshold[1];

    if (energySpike == 0)
       { energySpike = 0.5; WARNING("Automatically set energythreshold(0) to 0.5"); }

    if (energySlope == 0)
       { energySlope = 0.01; WARNING("Automatically set energythreshold(1) to 0.01"); }

    if (moveleftend && myDomain == 0) 
       { ERROR("string: moveleftend not yet implemented"); }

    if (yesclimbimage && myDomain == 0)
       { ERROR("string: yesclimbimage feature no longer supported. Try nebrelax."); }

    if (yesclimbimage && myDomain == 0)
       { INFO_Printf("String: The climbing image method will be applied after %d steps.\n",equilsteps); }

    if (relax_surround && myDomain == 0) 
       { ERROR("string: relax_surround feature no longer supported. Try nebrelax."); }

    /* Makes sure _SR1 and _SR2 (for configs A and B) are allocated */
    if(_SR1 == NULL || _SR2 == NULL)
    {
        if(_SR1 == NULL) ERROR("config1 is not set!");
        if(_SR2 == NULL) ERROR("config2 is not set!");
        return;
    }

    n = constrainatoms[0];
   
    fp = NULL; 
    if(myDomain==0) /* Master opens file */
    {
        fp=fopen("stringeng.out","w");
    }

    /* Allocate workspace variables for stringrelax */
    Ec = Ec_global = NULL; Fm = Fm_global = NULL; Ft = Ft_global = NULL; Fs = Fs_global = NULL; dR = dR_global = NULL; 
    allocChainVars_parallel(&Ec, &Ec_global, &Fm, &Fm_global, &Ft, &Ft_global, &Fs, &Fs_global, &dR, &dR_global);
    allocChainIndex_parallel(&left_index, &left_index_global, &right_index, &right_index_global);

    /* calculate center-of-mass shift from _SR1 (state A) to _SR2 (state B) */
    caldscom12();
    /* compute current constraint value (not really necessary) */
    constrain_dist=0;
    for(i=0;i<n;i++)
    {
        ipt=constrainatoms[i+1];
        ds0=(_SR2[ipt]-dscom12)-_SR1[ipt];
        constrain_dist+=ds0.norm2();
    }
    INFO_Printf("string[%d]: constrain_dist=%e s=%e\n", myDomain, constrain_dist, _nebinterp[myDomain]);

    /* fix constrained atoms for relaxing surrounding atoms */
    if(relax_surround==1)
    {
        for(i=0;i<n;i++)
        {
            ipt=constrainatoms[i+1];
            fixed[ipt]=1;
        }
    }

    /* store E0 which is the energy of the original leftendDomain */
    MDFrame::clearR0(); /* enforce neighbor list to be updated */
    calChainForce_parallel(n, relax_surround, Ec, Ec_global);
    E0 = Ec[0];

    if (manualcut)         
    {
       EmaxDomain = 0;
       reparamChain_parallel(moveleftend, moverightend, yesclimbimage, EmaxDomain, n, &lavg, &lavg0, dR, Ec,
                             left_index, left_index_global, right_index, right_index_global, 
                             energySlope, manualcut, manualleftend, manualrightend);
    }

    /* start stringrelax iteration */
    step0 = curstep; lavg = 0;
    for(curstep=step0;curstep<=(step0 + totalsteps);curstep++)    
    {
        /* This is in stringrelax_parallel_2, but we do not think it is necessary */
        //if(curstep==step0) MDFrame::clearR0();

        /* Get atomic forces to _Fc array along the entire chain */
        calChainForce_parallel(n, relax_surround, Ec, Ec_global);

        /* Find max-energy replica (i.e. the saddle point) */
        EmaxDomain = 0; Emax=E0;
        for(jj=1;jj<=_CHAINLENGTH;jj++)
        {
            if(Ec[jj]>Emax) { Emax=Ec[jj]; EmaxDomain=jj; }
        }

        /* calculate tangential vector, _Tan */
        calChainTang_parallel(n, Ec, dR, dR_global, &TanMag2);
         
        /* orthogonalize forces _Fc against local tangent vector _Tan */
        orthoChainForce_parallel(moveleftend, moverightend, yesclimbimage, EmaxDomain, 
                                 TanMag2, n, Fm, Fm_global, Ft, Ft_global);

        /* calculate the initial averaged chain length, lavg0 */
        if(curstep==step0)
        {
            lavg0 = 0;
            for(jj=1;jj<=_CHAINLENGTH;jj++) lavg0 += dR[jj];
            lavg0 /= _CHAINLENGTH;
        }

        if((myDomain==0) && (curstep%printfreq==0))  // Mar. 08 2007 Keonwook Kang 
        {
            INFO_Printf("curstep = %d\n",curstep);
            for(jj=0;jj<=_CHAINLENGTH;jj++)
                INFO_Printf("%8d %25.15e %25.15e %25.15e %25.15e\n", jj, Ec[jj]-Ec[0], Fm[jj], Ft[jj], dR[jj]);
            fprintf(fp,"%d ",curstep);
            for(jj=0;jj<=_CHAINLENGTH;jj++)
                fprintf(fp,"%8d %25.15e %25.15e %25.15e %25.15e ", jj, Ec[jj]-Ec[0], Fm[jj], Ft[jj], dR[jj]);
            fprintf(fp,"\n");
            fflush(fp); 
        }

        if(curstep<(step0 + totalsteps))
        {   /* move chains along force */
            if ( (myDomain==0) && (moveleftend<=0) ) 
            {
                /* do nothing */
            }
            else if ( (myDomain==_CHAINLENGTH) && (moverightend<=0) ) 
            {
                /* do nothing */
            }
            else 
            {
                rightendEnergyTooHigh = 0;
                if ( ( (Ec[_CHAINLENGTH]-Ec[0]) >= (Emax-Ec[0])*0.5 ) && (EmaxDomain >= 0.5*_CHAINLENGTH) )
                    rightendEnergyTooHigh = 1;

                myEnergyTooHigh = 0;
                if( (myDomain > 0) && (myDomain<_CHAINLENGTH) )
                    if ( (Ec[myDomain]-(Ec[myDomain-1]+Ec[myDomain+1])*0.5) > energySpike ) 
                         myEnergyTooHigh = 1;

                if ( rightendEnergyTooHigh && (myDomain==_CHAINLENGTH) )
                     myEnergyTooHigh = 1;

                if ( !rightendEnergyTooHigh ) 
                { /* update all domains by steepest descent */
                    for(i=0;i<n;i++) _Rc0[i]+=_Fc0[i]*_TIMESTEP;
                }

                if (myEnergyTooHigh)
                {
                    hinv = _H.inv();
                    for(i=0;i<n;i++)
                    {
                        ipt=constrainatoms[i+1];
                        _SR[ipt]=hinv*_Rc0[i];
                    }
                    if ( myDomain!=_CHAINLENGTH ) conj_fevalmax = cgstepsmid; 
                    else conj_fevalmax = cgstepsrightend;
                    if (curstep == step0)
                    {
                        INFO_Printf("stringrelax[%d]: conj_ftol = %g, conj_itmax = %d conj_fevalmax = %d conj_fixbox = %d\n", myDomain, conj_ftol, conj_itmax, conj_fevalmax, conj_fixbox);
                        if (conj_fevalmax == 0)
                        {
                            WARNING("cgstepmid not set! Recommended value is 1. Set it in nebspec(10).");
                        }
                        if (conj_fixbox != 1)
                        {
                            ERROR("conj_fixbox must be 1");
                        }
                    }
                    relax();
                    SHtoR();
                    /* copy _R to _Rc0 */
                    for(i=0;i<n;i++)
                    {
                        ipt=constrainatoms[i+1];
                        _Rc0[i]=_R[ipt];
                    }
                    INFO_Printf("after cgrelax: Ec[%d]-Ec[0] = %20.12e\n",myDomain,_EPOT-Ec[0]);
                }
            }

            /* reparameterize path (redistribution) */
            if(curstep%redistri_freq==0) 
            {
                 /* update new positions of the neighbor chains _Rc0 */
                 calChainTang_parallel(n, Ec, dR, dR_global, &TanMag2);
                 reparamChain_parallel(moveleftend, moverightend, yesclimbimage, EmaxDomain, n, &lavg, &lavg0, dR, Ec,
                           left_index, left_index_global, right_index, right_index_global, energySlope, 0, 0, 0);
            }

        }
    }

    free(Ec); free(Ec_global); free(Fm); free(Fm_global); free(Ft); free(Ft_global); 
    free(Fs); free(Fs_global); free(dR); free(dR_global);

    if(myDomain==0) fclose(fp);

    INFO_Printf("stringrelax[%d]: exit\n",myDomain);
}

/* a test version by William Kuykendall, temporarily kept until stringrelax_parallel is fully tested */
#include "stringrelax_parallel_2.cpp" 

int MDPARALLELFrame::readRchain_parallel()
{
    int i;
    char fname[200], fullname[200];
    char *buffer; char *pp, *q; 
    
    MPI_Bcast(incnfile, 200, MPI_CHAR,   0,MPI_COMM_WORLD);

    sprintf(fname,"%s.cpu%02d",incnfile,myDomain);
    LFile::SubHomeDir(fname, fullname);
    INFO("filename="<<fullname);
    LFile::LoadToString(fullname,buffer,0);
    
    pp=buffer;

    q=pp; pp=strchr(pp,'\n'); if(pp) *(char *)pp++=0;
    sscanf(q, "%d", &_CHAINLENGTH);

    if( numDomains > (_CHAINLENGTH+1) ) // Added by Keonwook Kang Apr/19/2010
    {
        ERROR("number of cpus("<<numDomains<<") is bigger than chainlength("<<_CHAINLENGTH<<")+1");
        return -1;
    }
    else
        _CHAINLENGTH=numDomains-1;
    INFO_Printf("[%d] : chainlength = %d\n",myDomain,_CHAINLENGTH);

    q=pp; pp=strchr(pp,'\n'); if(pp) *(char *)pp++=0;
    sscanf(q, "%d", constrainatoms);

    INFO("constrainatoms[0]="<<constrainatoms[0]);

    /* Bug fixed by Keonwook Kang Apr/20/2010
       : Read constrainatoms before initRchain_parallel() */
    for(i=0;i<constrainatoms[0];i++)
    {
            q=pp; pp=strchr(pp,'\n'); if(pp) *(char *)pp++=0;
            sscanf(q, "%d", constrainatoms+i+1);
    }

    if (_Rc0==NULL) // If _Rc0 is not allocated
    {
        INFO("Call AllocChain_parallel() from readRchain_parallel()");
        AllocChain_parallel(); initRchain_parallel();
    } else 
        WARNING("_Rc0 will be overwritten!!");

    for(i=0;i<constrainatoms[0];i++)
    {
        q=pp; pp=strchr(pp,'\n'); if(pp) *(char *)pp++=0;
        sscanf(q, "%lf %lf %lf",&(_Rc0[i].x),&(_Rc0[i].y),&(_Rc0[i].z));
    }

    Free(buffer);
    INFO_Printf("readRchain_parallel[%d] finished\n",myDomain);
    return 0;
}


int MDPARALLELFrame::writeRchain_parallel()
{
    int i;
    char fname[200], fullname[200];
    FILE *fp;
    
    MPI_Bcast(finalcnfile, 200, MPI_CHAR,   0,MPI_COMM_WORLD);

    sprintf(fname,"%s.cpu%02d",finalcnfile,myDomain);
    LFile::SubHomeDir(fname, fullname);

    fp=fopen(fullname,"w");

    fprintf(fp,"%d\n",_CHAINLENGTH);
    fprintf(fp,"%d\n",constrainatoms[0]);

    for(i=0;i<constrainatoms[0];i++)
        fprintf(fp,"%d\n",constrainatoms[i+1]);

    for(i=0;i<constrainatoms[0];i++)
        fprintf(fp,"%25.17e %25.17e %25.17e %d\n",_Rc0[i].x,_Rc0[i].y,_Rc0[i].z,group[constrainatoms[i+1]]);

    fclose(fp);
    
    INFO_Printf("writeRchain_parallel[%d] finished\n",myDomain);
    return 0;
}

void MDPARALLELFrame::copyRchaintoCN_parallel()
{
    int i, n, ipt;
    Matrix33 hinv;

    hinv = _H.inv();
    n=constrainatoms[0];
    for(i=0;i<n;i++)
    {
        ipt=constrainatoms[i+1];
        _SR[ipt]=hinv*_Rc0[i];
        //_SR[ipt].subint();                    /* 2010/04/09 Keonwook Kang */
    }
    INFO_Printf("copyRchaintoCN_parallel[%d] finished\n",myDomain);
}

void MDPARALLELFrame::copyCNtoRchain_masteronly()
{
    int i, n, ipt;

  if (myDomain == 0)
  {
    n=constrainatoms[0];
    for(i=0;i<n;i++)
    {
      ipt=constrainatoms[i+1];
      //ds=_SR[ipt]-dscom12; ds.subint(); /* 2007/3/6 Wei Cai */
      //_Rc0[i]=_H*ds;
      _Rc0[i]=_H*_SR[ipt];
    }
    INFO_Printf("copyCNtoRchain_masteronly[%d] finished\n",myDomain);
  }
}

void MDPARALLELFrame::copyCNtoRchain_parallel()
{
    int i, n, ipt;

    n=constrainatoms[0];
    for(i=0;i<n;i++)
    {
      ipt=constrainatoms[i+1];
      //ds=_SR[ipt]-dscom12; ds.subint(); /* 2007/3/6 Wei Cai */
      //_Rc0[i]=_H*ds;
      _Rc0[i]=_H*_SR[ipt];
    }
    INFO_Printf("copyCNtoRchain_parallel[%d] finished\n",myDomain);
}

int MDPARALLELFrame::readcnfile_parallel()
{
    char fname[200];

    MPI_Bcast(incnfile, 200, MPI_CHAR,   0,MPI_COMM_WORLD);
    sprintf(fname,"%s.cpu%02d",incnfile,myDomain);

    MDFrame::initcn.open(fname,LFile::O_Read);
    return MDFrame::initcn.read(this);
    INFO_Printf("readcnfile_parallel[%d] finished\n",myDomain);
}

int MDPARALLELFrame::writefinalcnfile_parallel(int zip, bool bg)
{
    char fname[200];

    MPI_Bcast(&writeall, 1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&writevelocity, 1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(finalcnfile, 200, MPI_CHAR,   0,MPI_COMM_WORLD);    
    sprintf(fname,"%s.cpu%02d",finalcnfile,myDomain);
    
    MDFrame::finalcn.open(fname);
    return MDFrame::finalcn.write(this,zip,bg);
    INFO_Printf("writefinalcnfile_parallel[%d] finished\n",myDomain);
}

void MDPARALLELFrame::writeatomeyecfgfile_parallel()
{
    char fname[200];

    MPI_Bcast(atomeyerepeat, 4,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&NCS, 1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&plot_map_pbc,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(plot_limits,10,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&plot_color_axis,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(plot_color_windows,100,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(_TOPOL, _NP,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&zipfiles, 1,MPI_INT,0,MPI_COMM_WORLD);


    if(plot_color_axis==2) calcentralsymmetry();

    MPI_Bcast(finalcnfile, 200, MPI_CHAR,   0,MPI_COMM_WORLD);    
    sprintf(fname,"%s.cpu%02d",finalcnfile,myDomain);
    
    MDFrame::writeatomeyecfgfile(fname);
    INFO_Printf("writeatomeyecfgfile_parallel[%d] finished\n",myDomain);
}
#endif


/********************************************************************/
/* End of class MDPARALLELFrame definition */
/********************************************************************/




/********************************************************************/
/* Main program starts below */
/********************************************************************/

#ifdef _TEST

#ifdef _PARALLEL
#define MDFrame MDPARALLELFrame
#endif

class TESTMD : public MDFrame
{
public:
    void potential()
    {
        /* Multi process function */
        /* a void potential */
        INFO("Empty potential");
    }
    virtual void initvars()
    {
        MDFrame::initvars();
        _RLIST=3.8;
    }
            
};

class TESTMD sim;

/* The main program is defined here */
#include "main.cpp"

#endif//_TEST
