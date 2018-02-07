/*
  md.cpp
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
    bindcommand(name,"stringrelax_parallel_1",{if(myDomain==0) Master_to_Slave("stringrelax_parallel_1");stringrelax_parallel_1();});
    bindcommand(name,"stringrelax_parallel_2",{if(myDomain==0) Master_to_Slave("stringrelax_parallel_2");stringrelax_parallel_2();});
    bindcommand(name,"stringrelax_parallel_3",{if(myDomain==0) Master_to_Slave("stringrelax_parallel_3");stringrelax_parallel_3();});
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

void MDPARALLELFrame::Master_to_Slave(char *cmd)
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
    int i;
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
    MPI_Bcast(nebspec,        10,                MPI_INT,   0,MPI_COMM_WORLD);
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
    MPI_Bcast(&extforce,      100,               MPI_DOUBLE,0,MPI_COMM_WORLD);   
    MPI_Bcast(&forcemul,      1,                 MPI_DOUBLE,0,MPI_COMM_WORLD);   

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

void MDPARALLELFrame::nebrelax_parallel()
{
    /* parallel nudged elastic band method
     * assuming Broadcast_Atoms have been called
     * every CPU has a copy of all atoms
     */
    int n, size, i, ipt, j, nin, nout;
    int moveleftend, moverightend;
    double s, dEmax, dEmin, r2, fr, fm2, minFm, maxFm, maxFs;
    double Fspring, springK, dR1norm, dR1norm2, dR2norm, dR2norm2, TanMag2;
    double *Ec, *Ec_global, *Fm, *Fm_global, *Ft, *Ft_global, *Fs, *Fs_global, *dR, *dR_global;
    Vector3 ds0, ds, dR1, dR2, tmpvec;
    Matrix33 hinv;
    FILE *fp;
    
    INFO_Printf("nebrelax[%d]: NEB Relax Parallel\n", myDomain);

    if( numDomains != (_CHAINLENGTH+1) )
    {
        ERROR("number of cpus("<<numDomains<<") does not match chainlength("<<_CHAINLENGTH<<")+1");
        return;
    }
    
    s = _nebinterp[myDomain];
    hinv = _H.inv();
    n = constrainatoms[0];
    
    /* Allocate data array along path */
    Ec=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* energy along path */
    Ec_global=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* for communication */
    Fm=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* force normal to path */
    Fm_global=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* for communication */
    Ft=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* force along path */
    Ft_global=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* for communication */
    Fs=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* spring force magnitude */
    Fs_global=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* for communication */
    dR=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* segment length in Angstrom */
    dR_global=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* for communication */

    if(myDomain==0) /* Master opens file */
    {
        fp=fopen("nebeng.out","w");
    }

    /* compute current constraint value (not really necessary) */
    constrain_dist=0;
    for(i=0;i<n;i++)
    {
        ipt=constrainatoms[i+1];
        ds0=(_SR2[ipt]-dscom12)-_SR1[ipt];
        constrain_dist+=ds0.norm2();
    }

    INFO_Printf("nebrelax[%d]: constrain_dist=%e s=%e\n", myDomain, constrain_dist, s);

    /* fix constrained atoms for relaxing surrounding atoms */
    if(nebspec[0]==1)
    {
        for(i=0;i<n;i++)
        {
            ipt=constrainatoms[i+1];
            fixed[ipt]=1;
        }
    }

    /* if moveleftend = 0, fix the left end replica
                      = 1, move the left end replica by the spring force in the normal to the real force direction
                           (PNAS vol.104 3031-3036, 2007)
                      = 2, move the left end replica along the chain tangential direction (constrained free fall)
                      = 3, move the left end replica along the real force direction (free fall) */
    moveleftend = nebspec[2]; moverightend = nebspec[3];

    /* start NEB Relax iteration */
    step0 = curstep;
    for(curstep=step0;curstep<=(step0 + totalsteps);curstep++)    
    {
        /* get forces, s is reaction coordinate */
        interpCN(s); /* interpolate all atoms for a given s */

        /* copy _Rc0 to _SR */
        for(i=0;i<n;i++)
        {
            ipt=constrainatoms[i+1];
            _SR[ipt]=hinv*_Rc0[i];
            //_SR[ipt].subint();                    /* 2010/04/09 Keonwook Kang */
        }

        /* relax surrounding atoms */
        if(nebspec[0]==1) relax(); // NOT WORKING !!

        /* enforces the neighborlist to be reconstructed at the very first step (important!!) */
        /* clear Verlet neighbor list => enforce neighbor list updated at every step */
        if(curstep==step0) MDFrame::clearR0();
        MDFrame::call_potential();

        memset(Ec,0,sizeof(double)*(_CHAINLENGTH+1));
        Ec[myDomain]=_EPOT;
        MPI_Allreduce(Ec,Ec_global,(_CHAINLENGTH+1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        memcpy(Ec,Ec_global,sizeof(double)*(_CHAINLENGTH+1));
        /* Now everybody has the entire Ec array */

        /* _Fc0 is the force on constrained atoms */
        for(i=0;i<n;i++)
        {
            ipt=constrainatoms[i+1];
            _Fc0[i]=_F[ipt];
        }
        
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
                    for(i=0;i<n;i++) // n: number of constrained atoms
                    {
                        _Tan[i] = (_Rc2[i]-_Rc0[i])*dEmax + (_Rc0[i]-_Rc1[i])*dEmin;
                    }
                }
                else
                {
                    for(i=0;i<n;i++) // n: number of constrained atoms
                    {
                        _Tan[i] = (_Rc2[i]-_Rc0[i])*dEmin + (_Rc0[i]-_Rc1[i])*dEmax;
                    }
                }
            }
            else
            {
                if (Ec[myDomain+1]>Ec[myDomain-1]) // uphill
                {
                    for(i=0;i<n;i++) // n: number of constrained atoms
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
            
            TanMag2 = 0;        
            for(i=0;i<n;i++) // compute magnitude squared of tangential vector
                TanMag2 += _Tan[i].norm2();
        }
        
        /* calculate tangential vector at the left end */
        if(myDomain==0)
        {
            TanMag2 = 0;
            for(i=0;i<n;i++)
            {
                _Tan[i]= _Rc2[i]-_Rc0[i];
                TanMag2 += _Tan[i].norm2();
            }
        }   
        /* calculate tangential vector at the right end */
        if((myDomain==_CHAINLENGTH))
        {
            TanMag2 = 0;
            for(i=0;i<n;i++)
            {
                _Tan[i]= _Rc0[i]-_Rc1[i];
                TanMag2 += _Tan[i].norm2();
            }
        }
         
        /* orthogonalize forces */
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
                tmpvec = _Fc0[i];   // store force
                _Fc0[i]-=(_Tan[i]*fr);
                fm2 += _Fc0[i].norm2();
                if((myDomain==0 && moveleftend==1) ||
                   (myDomain==_CHAINLENGTH && moverightend==1) ||
                   (myDomain==0 && moveleftend==3) ||
                   (myDomain==_CHAINLENGTH && moverightend==3))
                    _Fc0[i]=tmpvec;  // restore force
            }
            Fm[myDomain] = sqrt(fm2); /* store residual force magnitude */
        }

        //MPI_Allreduce(&Fm,&minFm,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
        MPI_Allreduce(Fm,Fm_global,(_CHAINLENGTH+1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        memcpy(Fm,Fm_global,sizeof(double)*(_CHAINLENGTH+1));
        MPI_Allreduce(Ft,Ft_global,(_CHAINLENGTH+1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        memcpy(Ft,Ft_global,sizeof(double)*(_CHAINLENGTH+1));

        /* springK: spring constant */
        if(nebspec[1]>0.0)
            springK = nebspec[1];
        else
        {
            minFm=maxFm=Fm[0];
            for(j=1;j<=_CHAINLENGTH;j++)
            {
                minFm = min(minFm,Fm[j]);
                maxFm = max(maxFm,Fm[j]);
            }
            //springK = minFm*_CHAINLENGTH;
            springK = maxFm*_CHAINLENGTH;
        }

        /* add spring force */
        memset(Fs,0,sizeof(double)*(_CHAINLENGTH+1));
        memset(dR,0,sizeof(double)*(_CHAINLENGTH+1));
        if((myDomain>0)&&(myDomain<_CHAINLENGTH))
        {
            /* abs(R_{j+1} - R_j) */
            dR2norm2=0;
            for(i=0;i<n;i++)
            {
                dR2=_Rc2[i]-_Rc0[i];
                dR2norm2+=dR2.norm2();
            }
            dR2norm=sqrt(dR2norm2);

            /* abs(R_j - R_{j-1}) */
            dR1norm2=0;
            for(i=0;i<n;i++)
            {
                dR1=_Rc0[i]-_Rc1[i];
                dR1norm2+=dR1.norm2();
            }
            dR1norm=sqrt(dR1norm2);

            dR[myDomain]=dR1norm;
            if(myDomain==_CHAINLENGTH-1)
                dR[myDomain+1]=dR2norm;
          
            Fspring = springK*(dR2norm - dR1norm);
            Fs[myDomain] = Fspring;
            
            for(i=0;i<n;i++)
                _Fc0[i] += _Tan[i] * (Fspring/sqrt(TanMag2));
        }
        MPI_Allreduce(dR,dR_global,(_CHAINLENGTH+1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        memcpy(dR,dR_global,sizeof(double)*(_CHAINLENGTH+1));

        /* calculate and orthogonalize spring force at the end */
        if((myDomain==0&&moveleftend==1) ||
           (myDomain==_CHAINLENGTH&&moverightend==1))
        {
            if(myDomain==0)
            {
                Fspring = springK; Fs[myDomain] = Fspring*sqrt(TanMag2);
            }
            else if(myDomain==_CHAINLENGTH)
            {
                Fspring = (-1.0)*springK; Fs[myDomain] = Fspring*sqrt(TanMag2);
            }

            /* store spring force in _Fc0[i] and store _Fc0[i] in _Tan[i] */
            for(i=0;i<n;i++)
            {   tmpvec  = _Fc0[i];
                _Fc0[i] = _Tan[i]*Fspring; // Now, _Fc0[i] = spring force
                _Tan[i] = tmpvec;          //      _Tan[i] = real force
            }
            /* orthogonalize spring force against real force direction */
            fr=0; TanMag2=0; 
            for(i=0;i<n;i++)
            {
                fr+=dot(_Fc0[i],_Tan[i]);
                TanMag2+=_Tan[i].norm2();
            }
            fr /= TanMag2;
            
            for(i=0;i<n;i++)
                _Fc0[i]-=(_Tan[i]*fr);
        }
        MPI_Allreduce(Fs,Fs_global,(_CHAINLENGTH+1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        memcpy(Fs,Fs_global,sizeof(double)*(_CHAINLENGTH+1));
        
        if(myDomain==0) /* Master print file */
        {
            if(curstep%printfreq==0)  // Mar. 08 2007 Keonwook Kang 
            {
                INFO_Printf("curstep = %d\n",curstep);
                for(j=0;j<=_CHAINLENGTH;j++)
                    INFO_Printf("%20.12e %25.15e %25.15e %25.15e %25.15e %25.15e\n",_nebinterp[j],Ec[j]-Ec[0], Fm[j], Ft[j], Fs[j], dR[j]);
                fprintf(fp,"%d ",curstep);
                for(j=0;j<=_CHAINLENGTH;j++)
                    fprintf(fp,"%20.12e %25.15e %25.15e %25.15e %25.15e %25.15e ",_nebinterp[j], Ec[j]-Ec[0], Fm[j], Ft[j], Fs[j], dR[j]);
                fprintf(fp,"\n");
                fflush(fp); 
            }
        }

        if(curstep<(step0 + totalsteps))
        {/* move NEB chains along force */
            if((myDomain>0)&&(myDomain<_CHAINLENGTH ||
               ((myDomain==0)&&(moveleftend>0)) ||
               ((myDomain==_CHAINLENGTH)&&(moverightend>0))))
                for(i=0;i<n;i++)
                    _Rc0[i]+=_Fc0[i]*_TIMESTEP;
        }
       
    }

    free(Ec);
    free(Ec_global);
    free(Fm);
    free(Fm_global);
    free(Ft);
    free(Ft_global);
    free(Fs);
    free(Fs_global);    
    free(dR);
    free(dR_global);

    if(myDomain==0) fclose(fp);

    INFO_Printf("nebrelax[%d]: exit\n",myDomain);

}

void MDPARALLELFrame::stringrelax_parallel()
{
    /* parallel string method
     * assuming Broadcast_Atoms have been called
     * every CPU has a copy of all atoms
     */
    int n, size, i, ipt, j, jmin, jmax, nin, nout;
    int moveleftend, moverightend, EmaxDomain;
    int yesclimbimage, islengthconstant;
    double s, dEmax, dEmin, r2, fr, fm2;
    double lavg0=0, lavg, TanMag2, E0, Emax;
    double *Ec, *Ec_global, *Fm, *Fm_global, *Ft, *Ft_global, *dR, *dR_global;
//    double *Lc, *Lc_global;
    Vector3 ds0, ds, dr, dR1, dR2, tmpvec;
    Matrix33 hinv;
    FILE *fp;
    
    INFO_Printf("stringrelax[%d]: String Relax Parallel\n", myDomain);

    if( numDomains != (_CHAINLENGTH+1) )
    {
        ERROR("number of cpus("<<numDomains<<") does not match chainlength("<<_CHAINLENGTH<<")+1");
        return;
    }
    
    s = _nebinterp[myDomain];
    hinv = _H.inv();
    n = constrainatoms[0];
    
    /* Allocate data array along path */
    Ec=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* energy along path */
    Ec_global=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* for communication */
    Fm=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* force normal to path */
    Fm_global=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* for communication */
    Ft=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* force along path */
    Ft_global=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* for communication */
    dR=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* segment length in Angstrom */
    dR_global=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* for communication */
//    Lc=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* */
//    Lc_global=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* for communication */

    if(myDomain==0) /* Master opens file */
    {
        fp=fopen("stringeng.out","w");
    }

    /* fix constrained atoms for relaxing surrounding atoms */
    if(nebspec[0]==1)
    {
        for(i=0;i<n;i++)
        {
            ipt=constrainatoms[i+1];
            fixed[ipt]=1;
        }
    }

    /* if moveleftend = 0, fix the left end replica
                      = 1, free the left end */
    moveleftend = nebspec[2]; moverightend = nebspec[3];
    yesclimbimage = nebspec[4]; islengthconstant = nebspec[5];
    if (yesclimbimage)
        INFO_Printf("String: The climbing image method will be applied after %d steps.\n",equilsteps);

    /* Initialize EmaxDomain depending on (moverightend) and (moveleftend) */
    if(moverightend && !(moveleftend))
        EmaxDomain=0;
    else if(moveleftend && !(moverightend))
        EmaxDomain=_CHAINLENGTH;
    else if(!(moveleftend) && !(moverightend))
        EmaxDomain=_CHAINLENGTH;

    /* Store the energy of the left end for future reference. 
       This was implemented for moving left end, but is also 
       applicable when the left end is fixed. Oct/04/2010 KW */
    for(i=0;i<_NP;i++) _SR[i]=_SR1[i];
    MDFrame::clearR0(); MDFrame::call_potential(); E0=_EPOT;

    /* start String Relax iteration */
    step0 = curstep;
    for(curstep=step0;curstep<=(step0 + totalsteps);curstep++)    
    {
        /* get forces, s is reaction coordinate */
        interpCN(s); /* interpolate all atoms for a given s */

        /* copy _Rc0 to _SR */
        for(i=0;i<n;i++)
        {
            ipt=constrainatoms[i+1];
            _SR[ipt]=hinv*_Rc0[i];
            //_SR[ipt].subint();                    /* 2010/04/09 Keonwook Kang */
        }

        /* relax surrounding atoms */
        if(nebspec[0]==1) relax(); // NOT WORKING !!

        /* clear neighbor list => enforce neighbor list updated at the very 1st step */
        if(curstep==step0) MDFrame::clearR0();
        MDFrame::call_potential();

        memset(Ec,0,sizeof(double)*(_CHAINLENGTH+1));
        Ec[myDomain]=_EPOT;
        MPI_Allreduce(Ec,Ec_global,(_CHAINLENGTH+1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        memcpy(Ec,Ec_global,sizeof(double)*(_CHAINLENGTH+1));
        /* Now everybody has the entire Ec array */

        /* _Fc0 is the force on constrained atoms */
        for(i=0;i<n;i++)
        {
            ipt=constrainatoms[i+1];
            _Fc0[i]=_F[ipt];
        }
        
        /* Find max-energy replica for climbing image method */
        //if ( yesclimbimage && (curstep>=equilsteps) && (curstep%nebspec[1]==0))
        //{
            Emax=Ec[0]; EmaxDomain=0;
            for(j=1;j<=_CHAINLENGTH;j++)
                if(Ec[j]>Emax)
                {
                    Emax=Ec[j]; EmaxDomain=j;
                }
            if(EmaxDomain==0 || EmaxDomain==_CHAINLENGTH)
                INFO_Printf("Warning: Max-energy domain = %d\n",EmaxDomain);
        //} 

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
    
        /* calculate tangential vector */
        if((myDomain>0)&&(myDomain<_CHAINLENGTH))
        {
            for(i=0;i<n;i++) // n: number of constrained atoms
                _Tan[i] = _Rc2[i]-_Rc1[i];
            
            TanMag2 = 0;        
            for(i=0;i<n;i++) // compute magnitude squared of tangential vector
                TanMag2 += _Tan[i].norm2();
        }
        /* calculate tangential vector at the left end */
        if(moveleftend && myDomain==0)
        {
            for(i=0;i<n;i++) // n: number of constrained atoms
                _Tan[i] = _Rc2[i]-_Rc0[i];

            TanMag2 = 0;
            for(i=0;i<n;i++) // compute magnitude squared of tangential vector
                TanMag2 += _Tan[i].norm2();
        }
        /* calculate tangential vector at the right end */
        if(moverightend && myDomain==_CHAINLENGTH)
        {
            for(i=0;i<n;i++) // n: number of constrained atoms
                _Tan[i] = _Rc0[i]-_Rc1[i];
            
            TanMag2 = 0;        
            for(i=0;i<n;i++) // compute magnitude squared of tangential vector
                TanMag2 += _Tan[i].norm2();
        }
                 
        /* orthogonalize forces */
        memset(Fm,0,sizeof(double)*(_CHAINLENGTH+1));
        memset(Ft,0,sizeof(double)*(_CHAINLENGTH+1));
        if((myDomain>0)&&(myDomain<_CHAINLENGTH))
        {
            fr=0;            
            for(i=0;i<n;i++)
                fr+=dot(_Fc0[i],_Tan[i]);
            Ft[myDomain] = fr/sqrt(TanMag2);
            fr /= TanMag2; // normalizing
            
            fm2=0;
            for(i=0;i<n;i++) /* orthogonalizing force against tangential vector */
            {
                _Fc0[i]-=(_Tan[i]*fr);
                fm2 += _Fc0[i].norm2();
                if (yesclimbimage && (curstep>=equilsteps) && myDomain==EmaxDomain)
                    _Fc0[i]-=(_Tan[i]*fr);
            }
            Fm[myDomain] = sqrt(fm2); /* store residual force magnitude */
        }
        /* orthogonalize forces at the end */
        if((moverightend && myDomain==_CHAINLENGTH)
           || (moveleftend && myDomain==0))
        {
            fr=0;            
            for(i=0;i<n;i++)
                fr+=dot(_Fc0[i],_Tan[i]);
            Ft[myDomain] = fr/sqrt(TanMag2);
            fr /= TanMag2; // normalizing
            
            fm2=0;
            for(i=0;i<n;i++) /* orthogonalizing force against tangential vector */
            {
                tmpvec = _Fc0[i];
                _Fc0[i]-=(_Tan[i]*fr);
                fm2 += _Fc0[i].norm2();
                _Fc0[i] = tmpvec;
            }
            Fm[myDomain] = sqrt(fm2); /* store residual force magnitude */
        }
        MPI_Allreduce(Fm,Fm_global,(_CHAINLENGTH+1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        memcpy(Fm,Fm_global,sizeof(double)*(_CHAINLENGTH+1));
        MPI_Allreduce(Ft,Ft_global,(_CHAINLENGTH+1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        memcpy(Ft,Ft_global,sizeof(double)*(_CHAINLENGTH+1));

        /* calculate dR */
        memset(dR,0,sizeof(double)*(_CHAINLENGTH+1));
        if((myDomain>0)&&(myDomain<=_CHAINLENGTH))
        {
            /* abs(R_j - R_{j-1}) */
            r2=0;
            for(i=0;i<n;i++)
            {
                dr=_Rc0[i]-_Rc1[i];
                r2+=dr.norm2();
            }
            dR[myDomain]=sqrt(r2);
        }
        MPI_Allreduce(dR,dR_global,(_CHAINLENGTH+1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        memcpy(dR,dR_global,sizeof(double)*(_CHAINLENGTH+1));

        /* calculate the initial averaged chain length, lavg0 */
        if(curstep==step0)
        {
            r2=0;
            for(j=1;j<=_CHAINLENGTH;j++)
                r2+=dR[j];
            
            lavg0=r2/_CHAINLENGTH;
        }
        
        if(myDomain==0) /* Master print file */
        {
            if(curstep%printfreq==0)  // Mar. 08 2007 Keonwook Kang 
            {
                INFO_Printf("curstep = %d\n",curstep);
                for(j=0;j<=_CHAINLENGTH;j++)
                    INFO_Printf("%8d %25.15e %25.15e %25.15e %25.15e\n", j, Ec[j]-Ec[0], Fm[j], Ft[j], dR[j]);
                fprintf(fp,"%d ",curstep);
                for(j=0;j<=_CHAINLENGTH;j++)
                    fprintf(fp,"%8d %25.15e %25.15e %25.15e %25.15e ", j, Ec[j]-Ec[0], Fm[j], Ft[j], dR[j]);
                fprintf(fp,"\n");
                fflush(fp); 
            }
        }

        if(curstep<(step0 + totalsteps))
        {   /* move chains along force */
            if((myDomain>0)&&(myDomain<_CHAINLENGTH))
                for(i=0;i<n;i++)
                {
		    //ipt = constrainatoms[i+1];
                    //if (fixed[ipt]==1) continue;
                    _Rc0[i]+=_Fc0[i]*_TIMESTEP;
                }

            /* move the ends along force */
            if((moverightend && myDomain==_CHAINLENGTH)
               || (moveleftend && myDomain==0))
                for(i=0;i<n;i++)
                {
		    //ipt = constrainatoms[i+1];
                    //if (fixed[ipt]==1) continue;
                    _Rc0[i]+=_Fc0[i]*_TIMESTEP;
                }

            /* update new positions of the neighbor chains */
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

            /* update dR */
            memset(dR,0,sizeof(double)*(_CHAINLENGTH+1));
            if((myDomain>0)&&(myDomain<=_CHAINLENGTH))
            {
                /* abs(R_j - R_{j-1}) */
                r2=0;
                for(i=0;i<n;i++)
                {
                    dr=_Rc0[i]-_Rc1[i];
                    r2+=dr.norm2();
                }
                dR[myDomain]=sqrt(r2);
            }
            MPI_Allreduce(dR,dR_global,(_CHAINLENGTH+1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
            memcpy(dR,dR_global,sizeof(double)*(_CHAINLENGTH+1));

            /* redistribute */
            if(curstep%nebspec[1]==0)
            {
                /* the averaged chain length, lavg */
  	        r2=0; jmin=1; jmax=_CHAINLENGTH;
                if (yesclimbimage&&(curstep>=equilsteps))
		{
                    // move only the right end
                    if (moverightend && !moveleftend)
  		        jmax=EmaxDomain;
                    // move only the left end
                    if (moveleftend && !moverightend)
 		        jmin=EmaxDomain+1;
		}
		for(j=jmin;j<=jmax;j++) r2+=dR[j];
                lavg=r2/(jmax-jmin+1);

                /* Enforce the string length varying 
                   in the following cases. */
                if ((!moverightend)&&(!moveleftend))
		{  // When the both ends are fixed
                    islengthconstant=0;
                }
                if ((moverightend&&(!moveleftend))
                     ||(moveleftend&&(!moverightend)))
		{  // When only one end is moving
                    if (yesclimbimage&&(curstep>=equilsteps))
		        islengthconstant=0;
		}
                if(islengthconstant == 1) lavg=lavg0; 

                /* redistribute */
		if(myDomain==EmaxDomain)
		{
  	            /* send _Rc0 to right node */
		    //if(myDomain<_CHAINLENGTH)
		    //    MPI_Send(_Rc0,n*3,MPI_DOUBLE,myDomain+1,
	            // 	         999,MPI_COMM_WORLD);
		    /* send _Rc0 to left node */
                    //if(myDomain>0)
		    //    MPI_Send(_Rc0,n*3,MPI_DOUBLE,myDomain-1,
		    // 	         199,MPI_COMM_WORLD);
		}
		else if(myDomain>EmaxDomain)
		{
		    /* receive _Rc1 from left node */
		    //MPI_Recv(_Rc1,n*3,MPI_DOUBLE,myDomain-1,
		    // 	     999,MPI_COMM_WORLD,&inStatus[0]);

                    if((myDomain<_CHAINLENGTH) ||
                       (myDomain==_CHAINLENGTH && moverightend))
		    {
			/* abs(R_j - R_{j-1})^2 */
			r2=0;
			for(i=0;i<n;i++)
			{
			    dr=_Rc0[i]-_Rc1[i];
			    r2+=dr.norm2();
			}

			for(i=0;i<n;i++)
			{
			    //if (fixed[ipt]==1) continue;
			    _Rc0[i]+=(_Rc0[i]-_Rc1[i])*(lavg/sqrt(r2)-1.0);
			}
                    }
		    /* send _Rc0 to right node */
		    //if(myDomain<_CHAINLENGTH)
	            // 	MPI_Send(_Rc0,n*3,MPI_DOUBLE,myDomain+1,
	            // 		 999,MPI_COMM_WORLD);
		}
		else if(myDomain<EmaxDomain)
		{
		    /* receive _Rc2 from right node */
		    //MPI_Recv(_Rc2,n*3,MPI_DOUBLE,myDomain+1,
		    // 	     199,MPI_COMM_WORLD,&inStatus[0]);

                    if((myDomain>0) ||
                       (myDomain==0 && moveleftend))
		    {
			/* abs(R_j - R_{j+1})^2 */
			r2=0;
			for(i=0;i<n;i++)
			{
			    dr=_Rc0[i]-_Rc2[i];
			    r2+=dr.norm2();
			}

			for(i=0;i<n;i++)
			{
			    //if (fixed[ipt]==1) continue;
			    _Rc0[i]+=(_Rc0[i]-_Rc2[i])*(lavg/sqrt(r2)-1.0);
			}
		    }
		    /* send _Rc0 to left node */
		    //if(myDomain>0)
		    //    MPI_Send(_Rc0,n*3,MPI_DOUBLE,myDomain-1,
	            // 	         199,MPI_COMM_WORLD);
		}
                /* update new positions of the neighbor chains */
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
            } // if(curstep%nebspec[1]==0)
        } // if(curstep<(step0 + totalsteps))
    }

    free(Ec);
    free(Ec_global);
    free(Fm);
    free(Fm_global);
    free(Ft);
    free(Ft_global);
    free(dR);
    free(dR_global);
//    free(Lc);
//    free(Lc_global);

    if(myDomain==0) fclose(fp);

    INFO_Printf("stringrelax[%d]: exit\n",myDomain);

}


/* A test version by Wei, started 4/9/2013 */
void MDPARALLELFrame::stringrelax_parallel_1()
{
    /* parallel string method
     * assuming Broadcast_Atoms have been called
     * every CPU has a copy of all atoms
     */
    int n, size, i, ipt, j, jmin, jmax, nin, nout;
    int moveleftend, moverightend, EmaxDomain;
    int yesclimbimage, islengthconstant;
    double s, dEmax, dEmin, r2, fr, fm2;
    double lavg0=0, lavg, TanMag2, E0, Emax;
    double *Ec, *Ec_global, *Fm, *Fm_global, *Ft, *Ft_global, *dR, *dR_global;
    Vector3 ds0, ds, dr, dR1, dR2, tmpvec;
    Matrix33 hinv;
    FILE *fp;
    
    INFO_Printf("stringrelax[%d]: String Relax Parallel 1\n", myDomain);

    if( numDomains != (_CHAINLENGTH+1) )
    {
        ERROR("number of cpus("<<numDomains<<") does not match chainlength("<<_CHAINLENGTH<<")+1");
        return;
    }
    
    s = _nebinterp[myDomain];
    hinv = _H.inv();
    n = constrainatoms[0];
    
    /* Allocate data array along path */
    Ec=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* energy along path */
    Ec_global=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* for communication */
    Fm=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* force normal to path */
    Fm_global=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* for communication */
    Ft=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* force along path */
    Ft_global=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* for communication */
    dR=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* segment length in Angstrom */
    dR_global=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* for communication */

    if(myDomain==0) /* Master opens file */
    {
        fp=fopen("stringeng.out","w");
    }

    /* fix constrained atoms for relaxing surrounding atoms */
    if(nebspec[0]==1)
    {
        for(i=0;i<n;i++)
        {
            ipt=constrainatoms[i+1];
            fixed[ipt]=1;
        }
    }

    /* if moveleftend = 0, fix the left end replica
                      = 1, free the left end */
    moveleftend = nebspec[2]; moverightend = nebspec[3];
    yesclimbimage = nebspec[4]; islengthconstant = nebspec[5];
    if (yesclimbimage)
        INFO_Printf("String: The climbing image method will be applied after %d steps.\n",equilsteps);

  /* Done by Keonwook?  EmaxDomain seems to be reset before it is used anyway...Commented out by WK 2013-12-06
    // Initialize EmaxDomain depending on (moverightend) and (moveleftend) 
    if(moverightend && !(moveleftend))
        EmaxDomain=0;
    else if(moveleftend && !(moverightend))
        EmaxDomain=_CHAINLENGTH;
    else if(!(moveleftend) && !(moverightend))
        EmaxDomain=_CHAINLENGTH;
  */

    /* Store the energy of the left end for future reference. 
       This was implemented for moving left end, but is also 
       applicable when the left end is fixed. Oct/04/2010 KW */
    for(i=0;i<_NP;i++) _SR[i]=_SR1[i];
    MDFrame::clearR0(); MDFrame::call_potential(); E0=_EPOT;

    /* start String Relax iteration */
    step0 = curstep;
    for(curstep=step0;curstep<=(step0 + totalsteps);curstep++)    
    {
        /* get forces, s is reaction coordinate */
        interpCN(s); /* interpolate all atoms for a given s */

        /* copy _Rc0 to _SR */
        for(i=0;i<n;i++)
        {
	    /* Original method only allowed for the relaxation of a small fraction
		of the total atoms, while the rest would maintain an interpolated
		position between A and B.  Thus, ever atom is interpolated above,
		while the atoms that are participating in the relaxation are given
		their previous positions from the real coordinates (_Rc0) below */
            ipt=constrainatoms[i+1];
            _SR[ipt]=hinv*_Rc0[i];
            //_SR[ipt].subint();                    /* 2010/04/09 Keonwook Kang */
        }

        /* relax surrounding atoms */
        if(nebspec[0]==1) relax(); // NOT WORKING !!

        /* clear neighbor list => enforce neighbor list updated at the very 1st step */
        if(curstep==step0) MDFrame::clearR0();
        MDFrame::call_potential();

        memset(Ec,0,sizeof(double)*(_CHAINLENGTH+1));
        Ec[myDomain]=_EPOT;
        MPI_Allreduce(Ec,Ec_global,(_CHAINLENGTH+1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        memcpy(Ec,Ec_global,sizeof(double)*(_CHAINLENGTH+1));
        /* Now everybody has the entire Ec array */

        /* _Fc0 is the force on constrained atoms */
        for(i=0;i<n;i++)
        {
            ipt=constrainatoms[i+1];
            _Fc0[i]=_F[ipt];
        }
        
        /* Find max-energy replica for climbing image method */
        Emax=Ec[0]; EmaxDomain=0;
        for(j=1;j<=_CHAINLENGTH;j++)
            if(Ec[j]>Emax)
            {
                Emax=Ec[j]; EmaxDomain=j;
            }
        if(EmaxDomain==0 || EmaxDomain==_CHAINLENGTH)
            INFO_Printf("Warning: Max-energy domain = %d\n",EmaxDomain);

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
    
        /* calculate tangential vector */
        if((myDomain>0)&&(myDomain<_CHAINLENGTH))
        {
            for(i=0;i<n;i++) // n: number of constrained atoms
                _Tan[i] = _Rc2[i]-_Rc1[i];
            
            TanMag2 = 0;        
            for(i=0;i<n;i++) // compute magnitude squared of tangential vector
                TanMag2 += _Tan[i].norm2();
        }
        /* calculate tangential vector at the left end */
        //if(moveleftend && myDomain==0) Better to have some version of the tangent vector than none
	if(myDomain==0)
        {
            for(i=0;i<n;i++) // n: number of constrained atoms
                _Tan[i] = _Rc2[i]-_Rc0[i];

            TanMag2 = 0;
            for(i=0;i<n;i++) // compute magnitude squared of tangential vector
                TanMag2 += _Tan[i].norm2();
        }
        /* calculate tangential vector at the right end */
        //if(moverightend && myDomain==_CHAINLENGTH) Better to have some version of the tangent vector than none
	if(myDomain==_CHAINLENGTH)
        {
            for(i=0;i<n;i++) // n: number of constrained atoms
                _Tan[i] = _Rc0[i]-_Rc1[i];
            
            TanMag2 = 0;        
            for(i=0;i<n;i++) // compute magnitude squared of tangential vector
                TanMag2 += _Tan[i].norm2();
        }
                 
        /* orthogonalize forces */
        memset(Fm,0,sizeof(double)*(_CHAINLENGTH+1));
        memset(Ft,0,sizeof(double)*(_CHAINLENGTH+1));
        if((myDomain>0)&&(myDomain<_CHAINLENGTH))
        {
            fr=0;            
            for(i=0;i<n;i++)
                fr+=dot(_Fc0[i],_Tan[i]);
            Ft[myDomain] = fr/sqrt(TanMag2);
            fr /= TanMag2; // normalizing
            
            fm2=0;
            for(i=0;i<n;i++) /* orthogonalizing force against tangential vector */
            {
                _Fc0[i]-=(_Tan[i]*fr);
                fm2 += _Fc0[i].norm2();
                if (yesclimbimage && (curstep>=equilsteps) && myDomain==EmaxDomain)
		/* We've already subtracted out the tangential component, but for the climbing
		    image we want it to climb to the highest point.  This is probably located
		    in the direction opposite where the force was pushing the configuration.
		    Thus, by subtracting the force a second time, we get a force in the opposite
		    direction as the original force, pushing the configuration up the energy hill. */
                    _Fc0[i]-=(_Tan[i]*fr);
            }
            Fm[myDomain] = sqrt(fm2); /* store residual force magnitude */
        }
        /* orthogonalize forces at the end */
        if((moverightend && myDomain==_CHAINLENGTH)
           || (moveleftend && myDomain==0))
        {
            fr=0;            
            for(i=0;i<n;i++)
                fr+=dot(_Fc0[i],_Tan[i]);
            Ft[myDomain] = fr/sqrt(TanMag2);
            fr /= TanMag2; // normalizing
            
            fm2=0;
            for(i=0;i<n;i++) /* orthogonalizing force against tangential vector */
	    /* We are not really orthogonalizing here, as the force itself is
		maintained, while the magnitude is adjusted.  This may be due to the
		fact that the tangents at the end are already not very good and
		allowing the end configuration to free fall down the energy slope
		is not necessarily a bad thing. */
            {
                tmpvec = _Fc0[i];
                _Fc0[i]-=(_Tan[i]*fr);
                fm2 += _Fc0[i].norm2();
                _Fc0[i] = tmpvec;
            }
            Fm[myDomain] = sqrt(fm2); /* store residual force magnitude */
        }
        MPI_Allreduce(Fm,Fm_global,(_CHAINLENGTH+1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        memcpy(Fm,Fm_global,sizeof(double)*(_CHAINLENGTH+1));
        MPI_Allreduce(Ft,Ft_global,(_CHAINLENGTH+1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        memcpy(Ft,Ft_global,sizeof(double)*(_CHAINLENGTH+1));

        /* calculate dR */
        memset(dR,0,sizeof(double)*(_CHAINLENGTH+1));
        if((myDomain>0)&&(myDomain<=_CHAINLENGTH))
        {
            /* abs(R_j - R_{j-1}): The distance between a configuration "D" and the one
		immediately to the left of "D"*/
            r2=0;
            for(i=0;i<n;i++)
            {
                dr=_Rc0[i]-_Rc1[i];
                r2+=dr.norm2();
            }
            dR[myDomain]=sqrt(r2);
        }
        MPI_Allreduce(dR,dR_global,(_CHAINLENGTH+1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        memcpy(dR,dR_global,sizeof(double)*(_CHAINLENGTH+1));

        /* calculate the initial averaged chain length, lavg0.  Done on first step */
        if(curstep==step0)
        {
            r2=0;
            for(j=1;j<=_CHAINLENGTH;j++)
                r2+=dR[j];
            
            lavg0=r2/_CHAINLENGTH;
        }
        
        if(myDomain==0) /* Master print file */
        {
            if(curstep%printfreq==0)  // Mar. 08 2007 Keonwook Kang 
            {
                INFO_Printf("curstep = %d\n",curstep);
                for(j=0;j<=_CHAINLENGTH;j++)
                    INFO_Printf("%8d %25.15e %25.15e %25.15e %25.15e\n", j, Ec[j]-Ec[0], Fm[j], Ft[j], dR[j]);
                fprintf(fp,"%d ",curstep);
                for(j=0;j<=_CHAINLENGTH;j++)
                    fprintf(fp,"%8d %25.15e %25.15e %25.15e %25.15e ", j, Ec[j]-Ec[0], Fm[j], Ft[j], dR[j]);
                fprintf(fp,"\n");
                fflush(fp); 
            }
        }

        if(curstep<(step0 + totalsteps))
        {   /* move chains along force */
            if((myDomain>0)&&(myDomain<_CHAINLENGTH))
                for(i=0;i<n;i++)
                {
		    //ipt = constrainatoms[i+1];
                    //if (fixed[ipt]==1) continue;
                    _Rc0[i]+=_Fc0[i]*_TIMESTEP;
                }

            /* move the ends along force */
            if((moverightend && myDomain==_CHAINLENGTH)
               || (moveleftend && myDomain==0))
                for(i=0;i<n;i++)
                {
		    //ipt = constrainatoms[i+1];
                    //if (fixed[ipt]==1) continue;
                    _Rc0[i]+=_Fc0[i]*_TIMESTEP;
                }

            /* update new positions of the neighbor chains */
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

            /* update dR */
            memset(dR,0,sizeof(double)*(_CHAINLENGTH+1));
            if((myDomain>0)&&(myDomain<=_CHAINLENGTH))
            {
                /* abs(R_j - R_{j-1}) */
                r2=0;
                for(i=0;i<n;i++)
                {
                    dr=_Rc0[i]-_Rc1[i];
                    r2+=dr.norm2();
                }
                dR[myDomain]=sqrt(r2);
            }
            MPI_Allreduce(dR,dR_global,(_CHAINLENGTH+1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
            memcpy(dR,dR_global,sizeof(double)*(_CHAINLENGTH+1));

            /* redistribute */
            if(curstep%nebspec[1]==0)
            {
                /* the average chain length, lavg, based on the section of
		    chain on the opposite side of the energy barrier from
		    the end that is free to move.  If both ends are fixed,
		    the entire energy barrier is used to get an average. */
  	        r2=0; jmin=1; jmax=_CHAINLENGTH;
                if (yesclimbimage&&(curstep>=equilsteps))
		{
                    // move only the right end
                    if (moverightend && !moveleftend)
  		        jmax=EmaxDomain;
                    // move only the left end
                    if (moveleftend && !moverightend)
 		        jmin=EmaxDomain+1;
		}
		for(j=jmin;j<=jmax;j++) r2+=dR[j];
                lavg=r2/(jmax-jmin+1);

	//	// Cai Note: These two sections changing islengthconstant "Need to change" //	//
                /* Enforce the string length varying 
                   in the following cases. */
                if ((!moverightend)&&(!moveleftend))
		{  // When the both ends are fixed
                    islengthconstant=0;
                }
                if ((moverightend&&(!moveleftend))
                     ||(moveleftend&&(!moverightend)))
		{  // When only one end is moving
                    if (yesclimbimage&&(curstep>=equilsteps))
		        islengthconstant=0;
		}
                //if(islengthconstant == 1) lavg=lavg0; 
                if(islengthconstant == 1)
                {
                    if (Ec[_CHAINLENGTH]<Ec[0]) 
                    {
			/* Make lavg shorter to pull the right end up
			    Also, lavg0 must also be shrunk to avoid oscillations
			*/
                        lavg = lavg * 0.9 + lavg0 * 0.1;
                        // Reduce the length of lavg0. 08/28/2013 Keonwook Kang
                        // The lower bound of lavg0 = 1 (Angstrom). 
                        lavg0 = lavg0 * 0.9 + 0.1; 
                    }
		         /* I believe this was also done to avoid oscillations */
                    else lavg0 = lavg0 * 0.9 + lavg * 0.1;
                }

                /* update new positions of the neighbor chains */
		// I don't think any changes have taken place in _Rc yet, so this doesn't do anything
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

              //if ( Ec[5] < Ec[0] + 0.01 ) /* the condition that A need to shift to the right */
              // If the chain is folded many times,
              if ( fabs(Ec[2]-Ec[1])<0.1 && 
                   fabs(Ec[3]-Ec[2])<0.1 && 
                   fabs(Ec[4]-Ec[2])<0.1 && 
                   fabs(Ec[5]-Ec[2])<0.1 ) 
              {
                    /* Move config along chain according to length */
                    double alpha;
                    alpha = 1.0 - (myDomain*1.0)/(_CHAINLENGTH*1.0);
                    if (alpha>1) alpha = 1; // How would we even get to this case?
                    for(i=0;i<n;i++)
                    {
			/* Moves Domain 0 to the position of Domain 1 and shifts the
			    rest of the domains to the right by proportionally less
			    until reaching the final domain, which stays in the same
			    spot
			*/
    	                _Rc0[i] = _Rc0[i]*(1.0-alpha)+_Rc2[i]*alpha;
    	            } 
              }
              else /* length reparametrization */
              {
                if(myDomain > 0)
                {
		/* If the current domain is not in the position it would occupy if
		    the entire chain was evenly spaced, it is shifted to the left
		    (top) or right (bottom) so it goes to the correct position.  
		    This is divided into two sections because if the configuration
		    is moving left, it is moved along the path between the current
		    position and the left position (similarly with the right
		    position for moving right).  This avoids problems where the
		    tangent or other path may move the configuration into some
		    high energy position.  Theoretically we are safer moving in
		    a direct line between two current configurations than trying
		    a move in some other direction.
		*/
                    /* Move config along chain according to length */
                    double length_from_A, alpha;
                    length_from_A = 0;
                    for(i=1;i<=myDomain;i++) length_from_A += dR[i];
                    if (length_from_A > myDomain*lavg)
                    { /* need to move left to shorten path */
                      alpha = 1.0 - (length_from_A - myDomain*lavg)/dR[myDomain];
                      if (alpha<0) alpha = 0;
                      for(i=0;i<n;i++)
		      {
                         _Rc0[i] = _Rc1[i]*(1.0-alpha)+_Rc0[i]*alpha;
	              } 
                    } 
                    else 
                    { /* need to move right to lengthen path */
                       if (myDomain<_CHAINLENGTH)
                       {
                          alpha = (myDomain*lavg - length_from_A)/dR[myDomain+1];
                          if (alpha>1) alpha = 1;
                          for(i=0;i<n;i++)
		          {
    	                     _Rc0[i] = _Rc0[i]*(1.0-alpha)+_Rc2[i]*alpha;
    	                  } 
	              } 
                    }
                 }
              }
            } // if(curstep%nebspec[1]==0)
        } // if(curstep<(step0 + totalsteps))
    }

    free(Ec);
    free(Ec_global);
    free(Fm);
    free(Fm_global);
    free(Ft);
    free(Ft_global);
    free(dR);
    free(dR_global);
//    free(Lc);
//    free(Lc_global);

    if(myDomain==0) fclose(fp);

    INFO_Printf("stringrelax[%d]: exit\n",myDomain);

}

/* A test version by William, started 12/06/2013 */
void MDPARALLELFrame::stringrelax_parallel_2()
{
    /* parallel string method
     * assuming Broadcast_Atoms have been called
     * every CPU has a copy of all atoms
     */
    int n, size, i, ipt, j, jmin, jmax, nin, nout, pullup = 50;
    int moveleftend, moverightend, EmaxDomain, Bnegative = 0;
    int yesclimbimage, islengthconstant, cleanupflag, count;
    double s, dEmax, dEmin, r2, fr, fm2, alpha;
    double lavg0=0, lavg, TanMag2, E0, Emax;
    double *Ec, *Ec_global, *Fm, *Fm_global, *Ft, *Ft_global, *dR, *dR_global;
    Vector3 ds0, ds, dr, dR1, dR2, tmpvec;
    Matrix33 hinv;
    FILE *fp;
    int *ind_left, *ind_left_global, *ind_right, *ind_right_global;
    int lowleft, lowright, leftendDomain, rightendDomain;
    double new_position, quarter_chain, new_totalR;
    int left_index, right_index;
    int cgrelaxsteps;
 
    INFO_Printf("stringrelax[%d]: String Relax Parallel 2\n", myDomain);

    if( numDomains != (_CHAINLENGTH+1) )
    {
        ERROR("number of cpus("<<numDomains<<") does not match chainlength("<<_CHAINLENGTH<<")+1");
        return;
    }
    
    s = _nebinterp[myDomain];
    hinv = _H.inv();
    n = constrainatoms[0];
    
    /* Allocate data array along path */
    Ec=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* energy along path */
    _Ec=Ec;
    Ec_global=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* for communication */
    Fm=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* force normal to path */
    Fm_global=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* for communication */
    Ft=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* force along path */
    Ft_global=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* for communication */
    dR=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* segment length in Angstrom */
    dR_global=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* for communication */

    ind_left=(int *)malloc(sizeof(int)*(_CHAINLENGTH+1)); /* configuration to the left */
    ind_left_global=(int *)malloc(sizeof(int)*(_CHAINLENGTH+1)); /* for communication */
    ind_right=(int *)malloc(sizeof(int)*(_CHAINLENGTH+1)); /* configuration to the right */
    ind_right_global=(int *)malloc(sizeof(int)*(_CHAINLENGTH+1)); /* for communication */
    memset(ind_left,0,sizeof(int)*(_CHAINLENGTH+1));
    memset(ind_right,0,sizeof(int)*(_CHAINLENGTH+1));
    memset(ind_left_global,0,sizeof(int)*(_CHAINLENGTH+1));
    memset(ind_right_global,0,sizeof(int)*(_CHAINLENGTH+1));

    if(myDomain==0) /* Master opens file */
    {
        fp=fopen("stringeng.out","w");
    }

    /* fix constrained atoms for relaxing surrounding atoms */
    if(nebspec[0]==1)
    {
        for(i=0;i<n;i++)
        {
            ipt=constrainatoms[i+1];
            fixed[ipt]=1;
        }
    }

    /* if moveleftend = 0, fix the left end replica
                      = 1, free the left end */
    moveleftend = nebspec[2]; moverightend = nebspec[3];
    yesclimbimage = nebspec[4]; islengthconstant = nebspec[5];
    cgrelaxsteps = nebspec[6];  if (cgrelaxsteps<=0) cgrelaxsteps = 10;
    if (yesclimbimage)
        INFO_Printf("String: The climbing image method will be applied after %d steps.\n",equilsteps);

  /* Done by Keonwook?  EmaxDomain seems to be reset before it is used anyway...Commented out by WK 2013-12-06
    // Initialize EmaxDomain depending on (moverightend) and (moveleftend) 
    if(moverightend && !(moveleftend))
        EmaxDomain=0;
    else if(moveleftend && !(moverightend))
        EmaxDomain=_CHAINLENGTH;
    else if(!(moveleftend) && !(moverightend))
        EmaxDomain=_CHAINLENGTH;
  */

    /* Store the energy of the left end for future reference. 
       This was implemented for moving left end, but is also 
       applicable when the left end is fixed. Oct/04/2010 KW */
    for(i=0;i<_NP;i++) _SR[i]=_SR1[i];
    MDFrame::clearR0(); MDFrame::call_potential(); E0=_EPOT;

    /* start String Relax iteration */
    step0 = curstep;
    for(curstep=step0;curstep<=(step0 + totalsteps);curstep++)    
    {
        /* get forces, s is reaction coordinate */
        interpCN(s); /* interpolate all atoms for a given s */

        /* copy _Rc0 to _SR */
        for(i=0;i<n;i++)
        {
	    /* Original method only allowed for the relaxation of a small fraction
		of the total atoms, while the rest would maintain an interpolated
		position between A and B.  Thus, ever atom is interpolated above,
		while the atoms that are participating in the relaxation are given
		their previous positions from the real coordinates (_Rc0) below */
            ipt=constrainatoms[i+1];
            _SR[ipt]=hinv*_Rc0[i];
            //_SR[ipt].subint();                    /* 2010/04/09 Keonwook Kang */
        }

        /* relax surrounding atoms */
        if(nebspec[0]==1) relax(); // NOT WORKING !!

        /* clear neighbor list => enforce neighbor list updated at the very 1st step */
        if(curstep==step0) MDFrame::clearR0();
	/* Get energy and forces */
        MDFrame::call_potential();
        memset(Ec,0,sizeof(double)*(_CHAINLENGTH+1));
        Ec[myDomain]=_EPOT;
        MPI_Allreduce(Ec,Ec_global,(_CHAINLENGTH+1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        memcpy(Ec,Ec_global,sizeof(double)*(_CHAINLENGTH+1));
        /* Now everybody has the entire Ec array */

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

        /* _Fc0 is the force on constrained atoms */
        for(i=0;i<n;i++)
        {
            ipt=constrainatoms[i+1];
            _Fc0[i]=_F[ipt];
        }
        
        /* Find max-energy replica for climbing image method */
        Emax=Ec[0]; EmaxDomain=0;
        for(j=1;j<=_CHAINLENGTH;j++)
            if(Ec[j]>Emax)
            {
                Emax=Ec[j]; EmaxDomain=j;
            }
        if(EmaxDomain==0 || EmaxDomain==_CHAINLENGTH)
            INFO_Printf("Warning: Max-energy domain = %d\n",EmaxDomain);
    
        /* calculate tangential vector */
        if((myDomain>0)&&(myDomain<_CHAINLENGTH))
        {
            for(i=0;i<n;i++) // n: number of constrained atoms
                _Tan[i] = _Rc2[i]-_Rc1[i];
            
            TanMag2 = 0;        
            for(i=0;i<n;i++) // compute magnitude squared of tangential vector
                TanMag2 += _Tan[i].norm2();
        }
        /* calculate tangential vector at the left end */
        //if(moveleftend && myDomain==0) Better to have some version of the tangent vector than none
	if(myDomain==0)
        {
            for(i=0;i<n;i++) // n: number of constrained atoms
                _Tan[i] = _Rc2[i]-_Rc0[i];

            TanMag2 = 0;
            for(i=0;i<n;i++) // compute magnitude squared of tangential vector
                TanMag2 += _Tan[i].norm2();
        }
        /* calculate tangential vector at the right end */
        //if(moverightend && myDomain==_CHAINLENGTH) Better to have some version of the tangent vector than none
	if(myDomain==_CHAINLENGTH)
        {
            for(i=0;i<n;i++) // n: number of constrained atoms
                _Tan[i] = _Rc0[i]-_Rc1[i];
            
            TanMag2 = 0;        
            for(i=0;i<n;i++) // compute magnitude squared of tangential vector
                TanMag2 += _Tan[i].norm2();
        }
                 
        /* orthogonalize forces */
        memset(Fm,0,sizeof(double)*(_CHAINLENGTH+1));
        memset(Ft,0,sizeof(double)*(_CHAINLENGTH+1));
        if((myDomain>0)&&(myDomain<_CHAINLENGTH))
        {
            fr=0;            
            for(i=0;i<n;i++)
                fr+=dot(_Fc0[i],_Tan[i]);
            Ft[myDomain] = fr/sqrt(TanMag2);
            fr /= TanMag2; // normalizing
            
            fm2=0;
            for(i=0;i<n;i++) /* orthogonalizing force against tangential vector */
            {
                _Fc0[i]-=(_Tan[i]*fr);
                fm2 += _Fc0[i].norm2();
                if (yesclimbimage && (curstep>=equilsteps) && myDomain==EmaxDomain)
		/* We've already subtracted out the tangential component, but for the climbing
		    image we want it to climb to the highest point.  This is probably located
		    in the direction opposite where the force was pushing the configuration.
		    Thus, by subtracting the force a second time, we get a force in the opposite
		    direction as the original force, pushing the configuration up the energy hill. */
                    _Fc0[i]-=(_Tan[i]*fr);
            }
            Fm[myDomain] = sqrt(fm2); /* store residual force magnitude */
        }
        /* orthogonalize forces at the end */
        if((moverightend && myDomain==_CHAINLENGTH)
           || (moveleftend && myDomain==0))
        {
            fr=0;            
            for(i=0;i<n;i++)
                fr+=dot(_Fc0[i],_Tan[i]);
            Ft[myDomain] = fr/sqrt(TanMag2);
            fr /= TanMag2; // normalizing
            
            fm2=0;
            for(i=0;i<n;i++) /* orthogonalizing force against tangential vector */
	    /* We are not really orthogonalizing here, as the force itself is
		maintained, while the magnitude is adjusted.  This may be due to the
		fact that the tangents at the end are already not very good and
		allowing the end configuration to free fall down the energy slope
		is not necessarily a bad thing. */
            {
                tmpvec = _Fc0[i];
                _Fc0[i]-=(_Tan[i]*fr);
                fm2 += _Fc0[i].norm2();
                _Fc0[i] = tmpvec;
            }
            Fm[myDomain] = sqrt(fm2); /* store residual force magnitude */
        }
        MPI_Allreduce(Fm,Fm_global,(_CHAINLENGTH+1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        memcpy(Fm,Fm_global,sizeof(double)*(_CHAINLENGTH+1));
        MPI_Allreduce(Ft,Ft_global,(_CHAINLENGTH+1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        memcpy(Ft,Ft_global,sizeof(double)*(_CHAINLENGTH+1));

        /* calculate dR */
        memset(dR,0,sizeof(double)*(_CHAINLENGTH+1));
        if((myDomain>0)&&(myDomain<=_CHAINLENGTH))
        {
            /* abs(R_j - R_{j-1}): The distance between a configuration "D" and the one
		immediately to the left of "D"*/
            r2=0;
            for(i=0;i<n;i++)
            {
                dr=_Rc0[i]-_Rc1[i];
                r2+=dr.norm2();
            }
            dR[myDomain]=sqrt(r2);
        }
        MPI_Allreduce(dR,dR_global,(_CHAINLENGTH+1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        memcpy(dR,dR_global,sizeof(double)*(_CHAINLENGTH+1));

        /* calculate the initial averaged chain length, lavg0.  Done on first step */
        if(curstep==step0)
        {
            r2=0;
            for(j=1;j<=_CHAINLENGTH;j++)
                r2+=dR[j];
            
            lavg0=r2/_CHAINLENGTH;
        }
        
        if(myDomain==0) /* Master print file */
        {
            if(curstep%printfreq==0)  // Mar. 08 2007 Keonwook Kang 
            {
                INFO_Printf("curstep = %d\n",curstep);
                for(j=0;j<=_CHAINLENGTH;j++)
                    INFO_Printf("%8d %25.15e %25.15e %25.15e %25.15e\n", j, Ec[j]-Ec[0], Fm[j], Ft[j], dR[j]);
                fprintf(fp,"%d ",curstep);
                for(j=0;j<=_CHAINLENGTH;j++)
                    fprintf(fp,"%8d %25.15e %25.15e %25.15e %25.15e ", j, Ec[j]-Ec[0], Fm[j], Ft[j], dR[j]);
                fprintf(fp,"\n");
                fflush(fp); 
            }
        }

        if(curstep<(step0 + totalsteps))
        {   /* move chains along force */
          /* do not allow other replicas to move if energy of right end is too high */
          if( ((Ec[_CHAINLENGTH]-Ec[0]) < (Emax-Ec[0])*0.5) || (EmaxDomain < 0.5*_CHAINLENGTH) ) {

            if((myDomain>0)&&(myDomain<_CHAINLENGTH))
                for(i=0;i<n;i++)
                {
		    //ipt = constrainatoms[i+1];
                    //if (fixed[ipt]==1) continue;
                    _Rc0[i]+=_Fc0[i]*_TIMESTEP;
                }
            /* move the ends along force */
            if((moverightend && myDomain==_CHAINLENGTH)
               || (moveleftend && myDomain==0))
                for(i=0;i<n;i++)
                {
		    //ipt = constrainatoms[i+1];
                    //if (fixed[ipt]==1) continue;
                    _Rc0[i]+=_Fc0[i]*_TIMESTEP;
                }

         	//if(curstep%nebspec[1]==0) reparametrization_with_trimming();
           }
           else {

            /* only move right end, no trimming */
            if (moverightend && myDomain==_CHAINLENGTH)
            { /* do CGRelax */
               /* copy _Rc0 to _SR */
               for(i=0;i<n;i++)
               {
                  ipt=constrainatoms[i+1];
                  _SR[ipt]=hinv*_Rc0[i];
               }
               conj_ftol = 1e-4; conj_itmax = 10; conj_fevalmax = 10;
               conj_fixbox = 1;
               relax();
               SHtoR();
               /* copy _R to _Rc0 */
               for(i=0;i<n;i++)
               {
                   ipt=constrainatoms[i+1];
                   _Rc0[i]=_R[ipt];
               }
               //Ec[myDomain]=_EPOT;
               INFO_Printf("Ec[%d]-Ec[0] = %20.12e\n",myDomain,_EPOT-Ec[0]);
            }

	    //if(curstep%nebspec[1]==0) reparametrization_with_trimming();

           }

#if 0
            for(i=0;i<n;i++)
            {
               ipt=constrainatoms[i+1];
               _SR[ipt]=hinv*_Rc0[i];
            }
            MDFrame::call_potential();
            INFO_Printf("[%d] _EPOT = %20.12e  Ec[%d] = %20.12e",myDomain,_EPOT,myDomain,Ec[myDomain]);

            /* need to modifiy this condition to avoid one more call_potential */
            if (_EPOT > Ec[myDomain]) /* if energy increases after steepest descent step, do CGRelax */
#endif
          if( (myDomain > 0) && (myDomain<_CHAINLENGTH) )
          {
            if ( (Ec[myDomain]-(Ec[myDomain-1]+Ec[myDomain+1])*0.5) > 0.5 ) 
            { /* do CGRelax */
               INFO_Printf("[%d] _EPOT = %20.12e  Ec[%d] = %20.12e  Need relax",myDomain,_EPOT,myDomain,Ec[myDomain]);
               /* copy _Rc0 to _SR */
               for(i=0;i<n;i++)
               {
                  ipt=constrainatoms[i+1];
                  _SR[ipt]=hinv*_Rc0[i];
               }
               conj_ftol = 1e-4; conj_itmax = cgrelaxsteps; conj_fevalmax = cgrelaxsteps;
               conj_fixbox = 1;
               relax();
               SHtoR();
               /* copy _R to _Rc0 */
               for(i=0;i<n;i++)
               {
                   ipt=constrainatoms[i+1];
                   _Rc0[i]=_R[ipt];
               }
               //Ec[myDomain]=_EPOT;
               INFO_Printf("Ec[%d]-Ec[0] = %20.12e\n",myDomain,_EPOT-Ec[0]);
            }
          }
         
	   if(curstep%nebspec[1]==0) reparametrization_with_trimming();

	   //if( curstep%nebspec[1]==0)  { /* pasting the content of reparameterization with trimming */
           //} /* end of reparameterization_with_trimming block */
	    
        } // if(curstep<(step0 + totalsteps))
    }

    free(Ec);
    free(Ec_global);
    free(Fm);
    free(Fm_global);
    free(Ft);
    free(Ft_global);
    free(dR);
    free(dR_global);
    free(ind_left);
    free(ind_left_global);
    free(ind_right);
    free(ind_right_global);

    if(myDomain==0) fclose(fp);

    INFO_Printf("stringrelax[%d]: exit\n",myDomain);

}

/* A test version by Wei, using CGRelax, started 2/8/2014 */
void MDPARALLELFrame::stringrelax_parallel_3()
{
    /* parallel string method
     * assuming Broadcast_Atoms have been called
     * every CPU has a copy of all atoms
     */
    int n, size, i, ipt, j, jmin, jmax, nin, nout, pullup = 50;
    int moveleftend, moverightend, EmaxDomain, Bnegative = 0;
    int yesclimbimage, islengthconstant, cleanupflag, count;
    double s, dEmax, dEmin, r2, fr, fm2, alpha;
    double lavg0=0, lavg, TanMag2, E0, Emax;
    double *Ec, *Ec_global, *Fm, *Fm_global, *Ft, *Ft_global, *dR, *dR_global;
    Vector3 ds0, ds, dr, dR1, dR2, tmpvec;
    Matrix33 hinv;
    FILE *fp;
    int *ind_left, *ind_left_global, *ind_right, *ind_right_global;
    int lowleft, lowright, leftendDomain, rightendDomain;
    double new_position, quarter_chain, new_totalR;
    int left_index, right_index;
    int cgrelaxsteps;
    
    INFO_Printf("stringrelax[%d]: String Relax Parallel 3\n", myDomain);

    if( numDomains != (_CHAINLENGTH+1) )
    {
        ERROR("number of cpus("<<numDomains<<") does not match chainlength("<<_CHAINLENGTH<<")+1");
        return;
    }
    
    s = _nebinterp[myDomain];
    hinv = _H.inv();
    n = constrainatoms[0];
    
    /* Allocate data array along path */
    Ec=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* energy along path */
    _Ec=Ec;
    Ec_global=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* for communication */
    Fm=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* force normal to path */
    Fm_global=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* for communication */
    Ft=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* force along path */
    Ft_global=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* for communication */
    dR=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* segment length in Angstrom */
    dR_global=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* for communication */

    ind_left=(int *)malloc(sizeof(int)*(_CHAINLENGTH+1)); /* configuration to the left */
    ind_left_global=(int *)malloc(sizeof(int)*(_CHAINLENGTH+1)); /* for communication */
    ind_right=(int *)malloc(sizeof(int)*(_CHAINLENGTH+1)); /* configuration to the right */
    ind_right_global=(int *)malloc(sizeof(int)*(_CHAINLENGTH+1)); /* for communication */
    memset(ind_left,0,sizeof(int)*(_CHAINLENGTH+1));
    memset(ind_right,0,sizeof(int)*(_CHAINLENGTH+1));
    memset(ind_left_global,0,sizeof(int)*(_CHAINLENGTH+1));
    memset(ind_right_global,0,sizeof(int)*(_CHAINLENGTH+1));

    if(myDomain==0) /* Master opens file */
    {
        fp=fopen("stringeng.out","w");
    }

    /* fix constrained atoms for relaxing surrounding atoms */
    if(nebspec[0]==1)
    {
        for(i=0;i<n;i++)
        {
            ipt=constrainatoms[i+1];
            fixed[ipt]=1;
        }
    }

    /* if moveleftend = 0, fix the left end replica
                      = 1, free the left end */
    moveleftend = nebspec[2]; moverightend = nebspec[3];
    yesclimbimage = nebspec[4]; islengthconstant = nebspec[5]; 
    cgrelaxsteps = nebspec[6];  if (cgrelaxsteps<=0) cgrelaxsteps = 10;
    if (yesclimbimage)
        //INFO_Printf("String: The climbing image method will be applied after %d steps.\n",equilsteps);
        INFO_Printf("String: yesclimbimage is specified, but is not implemented in this version.\n");

    /* Store the energy of the left end for future reference. 
       This was implemented for moving left end, but is also 
       applicable when the left end is fixed. Oct/04/2010 KW */
    for(i=0;i<_NP;i++) _SR[i]=_SR1[i];
    MDFrame::clearR0(); MDFrame::call_potential(); E0=_EPOT;

    /* start String Relax iteration */
    step0 = curstep;
    for(curstep=step0;curstep<=(step0 + totalsteps);curstep++)    
    {
        /* get forces, s is reaction coordinate */
        interpCN(s); /* interpolate all atoms for a given s */

        /* copy _Rc0 to _SR */
        for(i=0;i<n;i++)
        {
	    /* Original method only allowed for the relaxation of a small fraction
		of the total atoms, while the rest would maintain an interpolated
		position between A and B.  Thus, ever atom is interpolated above,
		while the atoms that are participating in the relaxation are given
		their previous positions from the real coordinates (_Rc0) below */
            ipt=constrainatoms[i+1];
            _SR[ipt]=hinv*_Rc0[i];
            //_SR[ipt].subint();                    /* 2010/04/09 Keonwook Kang */
        }

        /* relax surrounding atoms */
        if(nebspec[0]==1) relax(); // NOT WORKING !!

        /* clear neighbor list => enforce neighbor list updated at the very 1st step */
        if(curstep==step0) MDFrame::clearR0();
	/* Get energy and forces */
        MDFrame::call_potential();
        memset(Ec,0,sizeof(double)*(_CHAINLENGTH+1));
        Ec[myDomain]=_EPOT;
        MPI_Allreduce(Ec,Ec_global,(_CHAINLENGTH+1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        memcpy(Ec,Ec_global,sizeof(double)*(_CHAINLENGTH+1));
        /* Now everybody has the entire Ec array */

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

        /* _Fc0 is the force on constrained atoms */
        for(i=0;i<n;i++)
        {
            ipt=constrainatoms[i+1];
            _Fc0[i]=_F[ipt];
        }
        
        /* Find max-energy replica for climbing image method */
        Emax=Ec[0]; EmaxDomain=0;
        for(j=1;j<=_CHAINLENGTH;j++)
            if(Ec[j]>Emax)
            {
                Emax=Ec[j]; EmaxDomain=j;
            }
        if(EmaxDomain==0 || EmaxDomain==_CHAINLENGTH)
            INFO_Printf("Warning: Max-energy domain = %d\n",EmaxDomain);
    
        /* calculate tangential vector */
        if((myDomain>0)&&(myDomain<_CHAINLENGTH))
        {
            for(i=0;i<n;i++) // n: number of constrained atoms
                _Tan[i] = _Rc2[i]-_Rc1[i];
            
            TanMag2 = 0;        
            for(i=0;i<n;i++) // compute magnitude squared of tangential vector
                TanMag2 += _Tan[i].norm2();
        }
        /* calculate tangential vector at the left end */
        //if(moveleftend && myDomain==0) Better to have some version of the tangent vector than none
	if(myDomain==0)
        {
            for(i=0;i<n;i++) // n: number of constrained atoms
                _Tan[i] = _Rc2[i]-_Rc0[i];

            TanMag2 = 0;
            for(i=0;i<n;i++) // compute magnitude squared of tangential vector
                TanMag2 += _Tan[i].norm2();
        }
        /* calculate tangential vector at the right end */
        //if(moverightend && myDomain==_CHAINLENGTH) Better to have some version of the tangent vector than none
	if(myDomain==_CHAINLENGTH)
        {
            for(i=0;i<n;i++) // n: number of constrained atoms
                _Tan[i] = _Rc0[i]-_Rc1[i];
            
            TanMag2 = 0;        
            for(i=0;i<n;i++) // compute magnitude squared of tangential vector
                TanMag2 += _Tan[i].norm2();
        }
                 
        /* orthogonalize forces */
        memset(Fm,0,sizeof(double)*(_CHAINLENGTH+1));
        memset(Ft,0,sizeof(double)*(_CHAINLENGTH+1));
        if((myDomain>0)&&(myDomain<_CHAINLENGTH))
        {
            fr=0;            
            for(i=0;i<n;i++)
                fr+=dot(_Fc0[i],_Tan[i]);
            Ft[myDomain] = fr/sqrt(TanMag2);
            fr /= TanMag2; // normalizing
            
            fm2=0;
            for(i=0;i<n;i++) /* orthogonalizing force against tangential vector */
            {
                _Fc0[i]-=(_Tan[i]*fr);
                fm2 += _Fc0[i].norm2();
                if (yesclimbimage && (curstep>=equilsteps) && myDomain==EmaxDomain)
		/* We've already subtracted out the tangential component, but for the climbing
		    image we want it to climb to the highest point.  This is probably located
		    in the direction opposite where the force was pushing the configuration.
		    Thus, by subtracting the force a second time, we get a force in the opposite
		    direction as the original force, pushing the configuration up the energy hill. */
                    _Fc0[i]-=(_Tan[i]*fr);
            }
            Fm[myDomain] = sqrt(fm2); /* store residual force magnitude */
        }
        /* orthogonalize forces at the end */
        if((moverightend && myDomain==_CHAINLENGTH)
           || (moveleftend && myDomain==0))
        {
            fr=0;            
            for(i=0;i<n;i++)
                fr+=dot(_Fc0[i],_Tan[i]);
            Ft[myDomain] = fr/sqrt(TanMag2);
            fr /= TanMag2; // normalizing
            
            fm2=0;
            for(i=0;i<n;i++) /* orthogonalizing force against tangential vector */
	    /* We are not really orthogonalizing here, as the force itself is
		maintained, while the magnitude is adjusted.  This may be due to the
		fact that the tangents at the end are already not very good and
		allowing the end configuration to free fall down the energy slope
		is not necessarily a bad thing. */
            {
                tmpvec = _Fc0[i];
                _Fc0[i]-=(_Tan[i]*fr);
                fm2 += _Fc0[i].norm2();
                _Fc0[i] = tmpvec;
            }
            Fm[myDomain] = sqrt(fm2); /* store residual force magnitude */
        }
        MPI_Allreduce(Fm,Fm_global,(_CHAINLENGTH+1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        memcpy(Fm,Fm_global,sizeof(double)*(_CHAINLENGTH+1));
        MPI_Allreduce(Ft,Ft_global,(_CHAINLENGTH+1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        memcpy(Ft,Ft_global,sizeof(double)*(_CHAINLENGTH+1));

        /* calculate dR */
        memset(dR,0,sizeof(double)*(_CHAINLENGTH+1));
        if((myDomain>0)&&(myDomain<=_CHAINLENGTH))
        {
            /* abs(R_j - R_{j-1}): The distance between a configuration "D" and the one
		immediately to the left of "D"*/
            r2=0;
            for(i=0;i<n;i++)
            {
                dr=_Rc0[i]-_Rc1[i];
                r2+=dr.norm2();
            }
            dR[myDomain]=sqrt(r2);
        }
        MPI_Allreduce(dR,dR_global,(_CHAINLENGTH+1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        memcpy(dR,dR_global,sizeof(double)*(_CHAINLENGTH+1));

        /* calculate the initial averaged chain length, lavg0.  Done on first step */
        if(curstep==step0)
        {
            r2=0;
            for(j=1;j<=_CHAINLENGTH;j++)
                r2+=dR[j];
            
            lavg0=r2/_CHAINLENGTH;
        }
        
        if(myDomain==0) /* Master print file */
        {
            if(curstep%printfreq==0)  // Mar. 08 2007 Keonwook Kang 
            {
                INFO_Printf("curstep = %d\n",curstep);
                for(j=0;j<=_CHAINLENGTH;j++)
                    INFO_Printf("%8d %25.15e %25.15e %25.15e %25.15e\n", j, Ec[j]-Ec[0], Fm[j], Ft[j], dR[j]);
                fprintf(fp,"%d ",curstep);
                for(j=0;j<=_CHAINLENGTH;j++)
                    fprintf(fp,"%8d %25.15e %25.15e %25.15e %25.15e ", j, Ec[j]-Ec[0], Fm[j], Ft[j], dR[j]);
                fprintf(fp,"\n");
                fflush(fp); 
            }
        }

        if(curstep<(step0 + totalsteps))
        {   /* no longer move chains along force, use CGRelax at every step */

            if( (myDomain>0) ) /* need to refine this condition */
            { /* do CGRelax */
               INFO_Printf("[%d] _EPOT = %20.12e  Ec[%d] = %20.12e  Need relax",myDomain,_EPOT,myDomain,Ec[myDomain]);
               /* copy _Rc0 to _SR */
               for(i=0;i<n;i++)
               {
                  ipt=constrainatoms[i+1];
                  _SR[ipt]=hinv*_Rc0[i];
               }
               conj_ftol = 1e-4; conj_itmax = cgrelaxsteps; conj_fevalmax = cgrelaxsteps;
               conj_fixbox = 1;
               relax();
               SHtoR();
               /* copy _R to _Rc0 */
               for(i=0;i<n;i++)
               {
                   ipt=constrainatoms[i+1];
                   _Rc0[i]=_R[ipt];
               }
               //Ec[myDomain]=_EPOT;
               INFO_Printf("Ec[%d]-Ec[0] = %20.12e\n",myDomain,_EPOT-Ec[0]);
            }

            Ec[myDomain]=_EPOT;
            MPI_Allreduce(Ec,Ec_global,(_CHAINLENGTH+1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
            memcpy(Ec,Ec_global,sizeof(double)*(_CHAINLENGTH+1));
         
	   if(curstep%nebspec[1]==0) reparametrization_with_trimming();

	   //if( curstep%nebspec[1]==0)  { /* pasting the content of reparameterization with trimming */
           //} /* end of reparameterization_with_trimming block */
	    
        } // if(curstep<(step0 + totalsteps))
    }

    free(Ec);
    free(Ec_global);
    free(Fm);
    free(Fm_global);
    free(Ft);
    free(Ft_global);
    free(dR);
    free(dR_global);
    free(ind_left);
    free(ind_left_global);
    free(ind_right);
    free(ind_right_global);

    if(myDomain==0) fclose(fp);

    INFO_Printf("stringrelax[%d]: exit\n",myDomain);

}

void MDPARALLELFrame::reparametrization_with_trimming() 
{
    int n, size, i, ipt, j, jmin, jmax, nin, nout, right_index;
    int moveleftend, moverightend, EmaxDomain, left_index;
    int yesclimbimage, islengthconstant, new_endDomain;
    int lowleft, lowright, leftendDomain, rightendDomain;
    double s, r2, alpha, new_position, quarter_chain, new_totalR;
    double dEmax, dEmin, E0, Emax;
    double *Ec, *Ec_global, *dR, *dR_global, *Ap, *Ap_global;
    int *ind_left, *ind_left_global, *ind_right, *ind_right_global;
    Vector3 ds0, ds, dr, dR1, dR2, tmpvec;
    Matrix33 hinv;
    
    INFO_Printf("RP[%d]: Reparametrization With Trimming\n", myDomain);

    if( numDomains != (_CHAINLENGTH+1) )
    {
        ERROR("number of cpus("<<numDomains<<") does not match chainlength("<<_CHAINLENGTH<<")+1");
        return;
    }

//INFO_Printf("RP0\n");
    s = _nebinterp[myDomain];
//INFO_Printf("RP0.1\n");
    hinv = _H.inv();
//INFO_Printf("RP0.2\n");
    n = constrainatoms[0];
//INFO_Printf("RP0.3\n");
    
    /* Allocate data array along path */
    //Ec=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* energy along path */
    Ec=_Ec;
    Ec_global=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* for communication */
    dR=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* segment length in Angstrom */
    dR_global=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* for communication */
    ind_left=(int *)malloc(sizeof(int)*(_CHAINLENGTH+1)); /* configuration to the left */
    ind_left_global=(int *)malloc(sizeof(int)*(_CHAINLENGTH+1)); /* for communication */
    ind_right=(int *)malloc(sizeof(int)*(_CHAINLENGTH+1)); /* configuration to the right */
    ind_right_global=(int *)malloc(sizeof(int)*(_CHAINLENGTH+1)); /* for communication */
    memset(ind_left,0,sizeof(int)*(_CHAINLENGTH+1));
    memset(ind_right,0,sizeof(int)*(_CHAINLENGTH+1));
    memset(ind_left_global,0,sizeof(int)*(_CHAINLENGTH+1));
    memset(ind_right_global,0,sizeof(int)*(_CHAINLENGTH+1));

#ifdef _DEBUG
    Ap=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* alpha along path */
    Ap_global=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* for communication */
//INFO_Printf("RP0.4\n");
#endif


#ifdef _DEBUG   //////////NEW
    /* fix constrained atoms for relaxing surrounding atoms */
    if(nebspec[0]==1)
    {
        for(i=0;i<n;i++)
        {
            ipt=constrainatoms[i+1];
            fixed[ipt]=1;
        }
    }
#endif

//INFO_Printf("RP0.5\n");
    /* if moveleftend = 0, fix the left end replica
                      = 1, free the left end */
    moveleftend = nebspec[2]; moverightend = nebspec[3];
    yesclimbimage = nebspec[4]; islengthconstant = nebspec[5];

#if 0  /* no longer needed since Ec = _Ec,  Wei 2/8/2014 */
#ifdef _DEBUG  //////////////NEW
//INFO_Printf("RP0.6\n");
    for(i=0;i<_NP;i++) _SR[i]=_SR1[i];
//    MDFrame::clearR0(); MDFrame::call_potential(); E0=_EPOT;
//INFO_Printf("RP0.7: Energy: %g, Domain: %d\n", E0, myDomain);
        /* get forces, s is reaction coordinate */
        interpCN(s); /* interpolate all atoms for a given s */
//INFO_Printf("RP0.8\n");
        /* copy _Rc0 to _SR */
        for(i=0;i<n;i++)
        {
            ipt=constrainatoms[i+1];
            _SR[ipt]=hinv*_Rc0[i];
            //_SR[ipt].subint();                    /* 2010/04/09 Keonwook Kang */
	}
//INFO_Printf("RP0.9\n");
        /* clear neighbor list => enforce neighbor list updated at the very 1st step */
//	MDFrame::clearR0(); 
#endif

	MDFrame::call_potential();
//INFO_Printf("RP0.10\n");

    memset(Ec,0,sizeof(double)*(_CHAINLENGTH+1));

#ifdef _DEBUG
    memset(Ap,0,sizeof(double)*(_CHAINLENGTH+1));
#endif

    Ec[myDomain]=_EPOT;

#ifdef _DEBUG
INFO_Printf("Energy = %g\n", _EPOT);
#endif

    MPI_Allreduce(Ec,Ec_global,(_CHAINLENGTH+1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    memcpy(Ec,Ec_global,sizeof(double)*(_CHAINLENGTH+1));

#endif /* no longer needed since Ec = _Ec,  Wei 2/8/2014 */

    /* Now everybody has the entire Ec array */

    /* Find max-energy replica for climbing image method */
    Emax=Ec[0]; EmaxDomain=0;
    for(j=1;j<=_CHAINLENGTH;j++)
        if(Ec[j]>Emax)
        {
            Emax=Ec[j]; EmaxDomain=j;
        }
    if(EmaxDomain==0 || EmaxDomain==_CHAINLENGTH)
        INFO_Printf("Warning: Max-energy domain = %d\n",EmaxDomain);
    // Get the lowest energy configurations to the left and right of the
    // maximum energy configuration.
    lowleft = 0; lowright = _CHAINLENGTH;
    for(j=0;j<=EmaxDomain;j++) {
	//if(Ec[j]<Ec[lowleft]) lowleft = j;
        //prevent very shallow region near state A
	if(Ec[j]<Ec[lowleft] + 0.01 ) lowleft = j; /* Wei 2/8/2014 */
    }
    for(j=EmaxDomain;j<=_CHAINLENGTH;j++) {
	if(Ec[j]<Ec[lowright]) lowright = j;
    }
    // We don't necessarily want to redistribute between the two lowest configurations,
    // however.  We want to try to keep the left and right ends near the same energy.
    leftendDomain = lowleft; rightendDomain = lowright;
    while(Ec[leftendDomain+1]-Ec[0] < 0) leftendDomain++;
    while(Ec[rightendDomain-1]-Ec[0] < 0) rightendDomain--;

/*	for(i = 0; i <= _CHAINLENGTH; i++) {
//		ind_left[i] = 0;
		ind_left_global[i] = 0;
//		ind_right[i] = 0;
		ind_right_global[i] = 0;
	}
*/

//INFO_Printf("RP1\n");


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


    /* calculate dR */
    memset(dR,0,sizeof(double)*(_CHAINLENGTH+1));
    if((myDomain>0)&&(myDomain<=_CHAINLENGTH))
    {
        /* abs(R_j - R_{j-1}): The distance between a configuration "D" and the one
	    immediately to the left of "D" */
        r2=0;
        for(i=0;i<n;i++)
        {
            dr=_Rc0[i]-_Rc1[i];
            r2+=dr.norm2();
        }
        dR[myDomain]=sqrt(r2);
    }
    MPI_Allreduce(dR,dR_global,(_CHAINLENGTH+1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    memcpy(dR,dR_global,sizeof(double)*(_CHAINLENGTH+1));

#ifdef _DEBUG
INFO_Printf("RP2\n");
    for(i = 0; i <= _CHAINLENGTH; i++) {
	INFO_Printf("Configuration %d Energy: %g dR: %g \n",i,Ec[i]-Ec[0],dR[i]);
    }
INFO_Printf("moverightend = %d, moveleftend = %d, Energy Difference = %g\n", moverightend,moveleftend,Ec[0]-Ec[_CHAINLENGTH]);
#endif

    /* Redistribution! */
//INFO_Printf("RP2.3\n");
    new_totalR = 0.0;
    for(j=leftendDomain+1;j<=rightendDomain;j++) {
	new_totalR+=dR[j];

#ifdef _DEBUG
	INFO_Printf("Config %d dR: %g, TotalR: %g\n",j,dR[j],new_totalR);
#endif

    }
    new_position = myDomain*new_totalR/_CHAINLENGTH;
    left_index = leftendDomain;
    // The second condition on this while loop should be unnecessary...
    while(new_position > dR[left_index+1] && left_index < _CHAINLENGTH) {
	left_index++;
	new_position -= dR[left_index];
    }


 
//INFO_Printf("RP3\n");
    if (left_index == _CHAINLENGTH) {
	left_index = _CHAINLENGTH-1; 
	right_index = _CHAINLENGTH;
	new_position += dR[right_index];
    } else { right_index = left_index+1; }
    ind_left[myDomain] = left_index; ind_right[myDomain] = right_index;
    MPI_Allreduce(ind_left,ind_left_global,(_CHAINLENGTH+1),MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    memcpy(ind_left,ind_left_global,sizeof(int)*(_CHAINLENGTH+1));
    MPI_Allreduce(ind_right,ind_right_global,(_CHAINLENGTH+1),MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    memcpy(ind_right,ind_right_global,sizeof(int)*(_CHAINLENGTH+1));
    /* Now everyone knows what everyone needs */

#ifdef _DEBUG
INFO_Printf("RP3.1: Config 23: Left: %d; Right: %d\n",ind_left[23],ind_right[23]);
#endif

//	for(i = 0; i <= _CHAINLENGTH; i++) {
//		INFO_Printf("Config %d Left: %d; Right: %d\n",i,ind_left[i],ind_right[i]);
//	}
    /* receive _Rc1, _Rc2 from neighbors */
    nin = nout = 0;
//INFO_Printf("RP3.11\n");
    /* receive _Rc1 from left index, send _Rc0 to those who need it as a left index */
    MPI_Irecv(_Rc1,n*3,MPI_DOUBLE,ind_left[myDomain],
		MSG_NEB_ATOM,MPI_COMM_WORLD,&inRequests[nin]);
//INFO_Printf("RP3.12\n");
    nin ++;
//INFO_Printf("RP3.13\n");
    for(i = 0; i <= _CHAINLENGTH; i++) {
	if(ind_left[i] == myDomain) { 
//INFO_Printf("RP3.14\n");
	    MPI_Isend(_Rc0,n*3,MPI_DOUBLE,i,
			MSG_NEB_ATOM,MPI_COMM_WORLD,&outRequests[nout]);
	    nout ++;
	}
    }
//INFO_Printf("RP3.2\n", myDomain);
    /* receive _Rc2 from right index, send _Rc0 to those who need it as a right index */
    MPI_Irecv(_Rc2,n*3,MPI_DOUBLE,ind_right[myDomain],
		MSG_NEB_ATOM,MPI_COMM_WORLD,&inRequests[nin]);
    nin ++;
    for(i = 0; i <= _CHAINLENGTH; i++) {
	if(ind_right[i] == myDomain) { 
	    MPI_Isend(_Rc0,n*3,MPI_DOUBLE,i,
			MSG_NEB_ATOM,MPI_COMM_WORLD,&outRequests[nout]);
	    nout ++;
	}
    }
//INFO_Printf("RP3.3\n", myDomain);
    MPI_Waitall(nout, outRequests,outStatus);
    MPI_Waitall(nin,  inRequests, inStatus );
    /* Now we need to move everything */
//    new_position += dR[left_index];
    alpha = new_position/dR[right_index];

#ifdef _DEBUG
INFO_Printf("RP4: alpha = %g\n", alpha);
    Ap[myDomain] = alpha;
    MPI_Allreduce(Ap,Ap_global,(_CHAINLENGTH+1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    memcpy(Ap,Ap_global,sizeof(double)*(_CHAINLENGTH+1));
	for(i = 0; i <= _CHAINLENGTH; i++) {
		INFO_Printf("Config %d Left: %d; Right: %d; Alpha: %g \n",i,ind_left[i],ind_right[i],Ap[i]);
	}
#endif



    for(i=0;i<n;i++) { 
	_Rc0[i] = _Rc1[i]*(1.0-alpha)+_Rc2[i]*alpha;
    }

	/* This only works when we are applying this to every atom */
        /* copy _Rc0 to _SR */
        for(i=0;i<n;i++)
        {
            ipt=constrainatoms[i+1];
            _SR[ipt]=hinv*_Rc0[i];
            //_SR[ipt].subint();                    /* 2010/04/09 Keonwook Kang */
	}

#ifdef _DEBUG
        /* clear neighbor list => enforce neighbor list updated at the very 1st step */
	MDFrame::clearR0(); MDFrame::call_potential();
    memset(Ec,0,sizeof(double)*(_CHAINLENGTH+1));
    memset(Ec_global,0,sizeof(double)*(_CHAINLENGTH+1));
    Ec[myDomain]=_EPOT;
INFO_Printf("Energy = %g\n", _EPOT);
    MPI_Allreduce(Ec,Ec_global,(_CHAINLENGTH+1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    memcpy(Ec,Ec_global,sizeof(double)*(_CHAINLENGTH+1));
    /* Now everybody has the entire Ec array */
    for(i = 0; i <= _CHAINLENGTH; i++) {
	INFO_Printf("Configuration %d Energy: %g\n",i,Ec[i]-Ec[0]);
    }
INFO_Printf("RP5\n");
#endif



#ifdef _DEBUG
/* call_potential first */
/* put this at the end of the function to print out eng to screen */
       
    if(myDomain==0) /* Master print to screen */
    {
	for(j=0;j<=_CHAINLENGTH;j++) INFO_Printf("%8d %25.15e\n", j, Ec[j]-Ec[0]);
    }

    /* Everyone print to screen */
    INFO_Printf("CPU[%d] needs config %d (left) and %d (right)\n", myDomain, ind_left[myDomain], ind_right[myDomain]);
#endif



    //free(Ec);  /* Wei 2/8/2014 */
    free(Ec_global);
    free(dR);
    free(dR_global);
    free(ind_left);
    free(ind_left_global);
    free(ind_right);
    free(ind_right_global);
}


int MDPARALLELFrame::readRchain_parallel()
{
    int i;
    char fname[200], fullname[200];
    char *buffer; char *pp, *q; 
    FILE *fp;
    
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
        INFO("Warning: _Rc0 will be overwritten!!");

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
