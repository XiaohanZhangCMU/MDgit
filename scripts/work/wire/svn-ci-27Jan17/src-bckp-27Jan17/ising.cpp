// Last Modified : Fri Sep  5 18:24:11 2008

#include "ising.h"

void IsingFrame::initparser()
{
    char s[100];

    /* Parser */
    bindvar("myname",myname,STRING);
    bindvar("command",command,STRING);
    bindvar("input",input,DOUBLE);
    bindvar("NX",&_NX,INT);
    bindvar("NY",&_NY,INT);
    bindvar("NZ",&_NZ,INT);
    bindvar("nsample",&_NSAMPLE,INT);
    bindvar("nsuccess",&_NSUCCESS,INT);

    /* Simulation setting */
    bindvar("kBT",&kBT,DOUBLE);
    bindvar("J",&J,DOUBLE);
    bindvar("H",&H,DOUBLE);
    
    bindvar("totalsteps",&totalsteps,INT);
    bindvar("savepropfreq",&savepropfreq,INT);
    bindvar("randseed",&_RANDSEED,INT);

    /* Forward Flux Sampling (FFS) parameters */
    bindvar("curstep",&curstep,INT);
    bindvar("continue_curstep",&continue_curstep,INT);
    bindvar("N_lgst_cluster",&N_lgst_cluster,INT);
    bindvar("N_lgst_avg",&N_lgst_avg,INT);
    bindvar("saveFFScn",&saveFFScn,INT);
    bindvar("FFSoption",&FFSoption,INT);
    bindvar("FFSpruning",&FFSpruning,INT);
    bindvar("FFScn_weight",&FFScn_weight,DOUBLE);
    bindvar("FFS_Pp",&FFS_Pp,DOUBLE);
    bindvar("FFS0_check",&FFS0_check,INT);
    bindvar("lambda_A",&lambda_A,INT);
    bindvar("lambda_B",&lambda_B,INT);
    bindvar("FFScurstep",&FFScurstep,INT);
    bindvar("FFSautoend",&FFSautoend,INT);

    /* Umbrella Sampling parameters */
    bindvar("YES_UMB",&YES_UMB,INT);
    bindvar("UMB_K",&UMB_K,DOUBLE);
    bindvar("UMB_equilstep",&UMB_equilstep,INT);
    bindvar("UMB_curstep",&UMB_curstep,INT);
    bindvar("n_center",&n_center,INT);
    bindvar("delta_n",&delta_n,INT);
    bindvar("Kinetic",&Kinetic,INT);
    bindvar("Kinetic_Time",&Kinetic_Time,DOUBLE);
    bindvar("Kinetic_Swip",&Kinetic_Swip,INT);

    /* File input and output */
    bindvar("savecn",&savecn,INT);
    bindvar("saveprop",&saveprop,INT);
    bindvar("savecnfreq",&savecnfreq,INT);
    bindvar("savepropfreq",&savepropfreq,INT);
    bindvar("printfreq",&printfreq,INT);
    bindvar("calcryfreq",&calcryfreq,INT); 
    bindvar("calcryfreq1",&calcryfreq1,INT);
    bindvar("filecounter",&filecounter,INT);
    bindvar("FFSfilecounter",&FFSfilecounter,INT);
    bindvar("incnfile",incnfile,STRING);
    bindvar("finalcnfile",finalcnfile,STRING);
    bindvar("outpropfile",outpropfile,STRING);
    bindvar("intercnfile",intercnfile,STRING);
    bindvar("FFScnfile",FFScnfile,STRING);
    bindvar("zipfiles",&zipfiles,INT);

    /* Visualization */
    bindvar("win_width",&win_width,INT);
    bindvar("win_height",&win_height,INT);
    bindvar("plotfreq",&plotfreq,INT);
    bindvar("atomradius",atomradius,DOUBLE);

    bindvar("backgroundcolor",backgroundcolor,STRING);
    bindvar("rotateangles",rotateangles,DOUBLE);
    for(int i=0;i<MAXSPECIES;i++)
    {
        sprintf(s,"atomcolor%d",i);
        bindvar(s,atomcolor[i],STRING);
    }
    bindvar("atomcolor",atomcolor[0],STRING);
    for(int i=0;i<MAXCOLORS;i++)
    {
        sprintf(s,"color%02d",i);
        bindvar(s,colornames[i],STRING);
    }
}    

int IsingFrame::exec(const char *name)
{
    if(Organizer::exec(name)==0) return 0;

    /* Parser */
    bindcommand(name,"runcommand",runcommand());

    bindcommand(name,"initspin",initspin());    
    bindcommand(name,"initpath",initpath((int)input[0]));

    /* Monte Carlo simulation */
    bindcommand(name,"MCrun",MCrun());    
    bindcommand(name,"srand",srand((unsigned int)_RANDSEED));
    bindcommand(name,"srand48",srand48((unsigned int)_RANDSEED));
    bindcommand(name,"srandbytime",srand((unsigned int)time(NULL)));
    bindcommand(name,"srand48bytime",srand48((unsigned int)time(NULL)));

    bindcommand(name,"copyPathtoS",copyPathtoS((int)input[0],(int)input[1]));
    bindcommand(name,"SampleMCpath",SampleMCpath(input[0],input[1],(int)input[2],(int)input[3]));
    bindcommand(name,"WalkonChain",WalkonChain(input[0],input[1],(int)input[2]));
    bindcommand(name,"ComputeSuccessRate",_NSUCCESS=ComputeSuccessRate(_NSAMPLE,input[0],input[1],
                                                                       (int)input[2]));
    bindcommand(name,"AnalyzeConfigs",AnalyzeConfigs((int)input[0],(int)input[1],_NSAMPLE,
                                                     input[2],input[3],(int)input[4]));
    /* FFS implementation */
    bindcommand(name,"FFSprocess",FFSprocess());
    bindcommand(name,"calcrystalorder",calcrystalorder());
    bindcommand(name,"allocDFS",allocDFS());    
    bindcommand(name,"assign_Lam",assign_Lam());
            
    /* File input and output */
    bindcommand(name,"writecn",writefinalcnfile());
    bindcommand(name,"readcn",readcn());
    bindcommand(name,"openintercnfile",openintercnfile());
    bindcommand(name,"openFFScnfile",openFFScnfile());
    bindcommand(name,"writeintercn",writeintercnfile());
    bindcommand(name,"writeFFScn",writeFFScnfile());
    bindcommand(name,"setfilecounter",setfilecounter());
    bindcommand(name,"setFFSfilecounter",setFFSfilecounter());
    bindcommand(name,"openpropfile",openpropfile());
    bindcommand(name,"writepath",writepath((int)input[0]));
    bindcommand(name,"readpath",readpath((int)input[0]));

    /* Visualization */
    bindcommand(name,"openwin",openwin());
    bindcommand(name,"plot",plot());
    bindcommand(name,"alloccolors",alloccolors());
    bindcommand(name,"alloccolorsX",alloccolorsX());
    bindcommand(name,"rotate",rotate());
    bindcommand(name,"saverot",saverot());
    bindcommand(name,"reversergb",win->reversergb());
    bindcommand(name,"wintogglepause",wintogglepause());
    
    return -1;
}

void IsingFrame::runcommand()
{
    char extcomm[400];
    
    LFile::SubHomeDir(command,command);
    INFO("run:  "<<command);
    if(nolog)
        system(command);
    else
    {
        if (ncom==0) system("rm -f B*.log");
        strcpy(extcomm,command);
        strcpy(extcomm+strlen(extcomm)," > B");
        sprintf(extcomm+strlen(extcomm),"%d.log",ncom);
        system(extcomm); 
        sprintf(extcomm,"echo -- '\n'B%d.log: \"%s\"'\n' >> B.log",ncom,command);
        system(extcomm);
        sprintf(extcomm,"cat B%d.log >> B.log",ncom);
        system(extcomm);
        ncom++;
    }
}

void IsingFrame::initvars()
{
    initparser();
    strcpy(command,"echo Hello World");
    
    strcpy(incnfile,"spin.cn");
    strcpy(intercnfile,"inter.cn");
    strcpy(FFScnfile,"FFS.cn");
    strcpy(finalcnfile,"final.cn");
    strcpy(outpropfile,"prop.out");

    strcpy(bondcolor,"red");
    strcpy(highlightcolor,"purple");

    for(int i=0;i<MAXSPECIES;i++)
    {
        atomradius[i]=0.1;
        strcpy(atomcolor[i],"orange");
    }
    for(int i=0;i<MAXCOLORS;i++)
    {
        strcpy(colornames[i],"gray50");
    }

}

#ifdef _USETCL
/********************************************************************/
/* Initialize Tcl parser */
/********************************************************************/
int IsingFrame::Tcl_AppInit(Tcl_Interp *interp)
{
    /*
     * Tcl_Init reads init.tcl from the Tcl script library.
     */
    if (Tcl_Init(interp) == TCL_ERROR) {
        return TCL_ERROR;
    }
    /*
     * Register application-specific commands.
     */
    Tcl_CreateCommand(interp, "MD++", MD_Cmd,
                      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);

    Tcl_CreateCommand(interp, "MD++_Set", MD_SetVar,
                      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);

    Tcl_CreateCommand(interp, "MD++_Get", MD_GetVar,
                      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);

    /*
      Random_Init(interp);
      Blob_Init(interp);
    */
    return TCL_OK;
}

int Tcl_Main_Wrapper(int argc, char *argv[])
{
    int i, pid;
    char oldargv[1000], newargv[1000], cmd[1000], c, *p, *filename;
    FILE *oldfile, *newfile;

    pid=getpid();
    for(i=1;i<argc;i++)
    {
        strcpy(oldargv,argv[i]);
        /* identify suffix .script */
        p = strrchr(oldargv,'.');
        if(p!=NULL) p = strstr(p,".script");
        if(p!=NULL)
        {
            /* identify file name */
            filename = strrchr(oldargv,'/');
            if(filename==NULL)
                filename = oldargv;
            else
                filename ++;

            p[0] = 0; sprintf(newargv,"/tmp/%s.tmp%d.%s",filename,pid,"tcl");
            p[0] = '.';

            INFO_Printf("replace: %s -> %s\n",oldargv,newargv);

            strcpy(argv[i],newargv);

            /* create new file */
            oldfile = fopen(oldargv,"r");
            newfile = fopen(newargv,"w");

            if(oldfile==NULL) ERROR("old input file open failure");
            if(newfile==NULL) ERROR("new input file open failure");
            fprintf(newfile,"# -*-shell-script-*-\n");
            fprintf(newfile,"#TCL script file created from %s\n\n",oldargv);
            fprintf(newfile,"MD++ {\n");
            for(;;)
            {
                c=getc(oldfile);
                if(c==EOF) break;
                putc(c,newfile);
            }
            fprintf(newfile,"\n}\n");
            fclose(oldfile);
            fclose(newfile);
        }
    }
    INFO_Printf("\n");

    setenv("TCL_LIBRARY","/usr/share/tcl8.4",0);
    Tcl_Main(argc, argv, Tcl_AppInit);

    /* remove tmp files */
    sprintf(cmd,"rm -f /tmp/*.tmp%d.tcl\n",pid);
    //INFO_Printf(cmd);
    system(cmd);

    return TCL_OK;
}
#endif



/* Memory Allocation */
void IsingFrame::Alloc()
{
    Ntot = _NX*_NY*_NZ;

    Realloc(_S,int,Ntot);
    memset(_S,-1,sizeof(int)*Ntot);

    Realloc(_S0,int,Ntot);
    memset(_S0,-1,sizeof(int)*Ntot);    
}

void IsingFrame::copyStoS0()
{
    memcpy(_S0,_S,sizeof(int)*Ntot);
}

void IsingFrame::copyS0toS()
{
    memcpy(_S,_S0,sizeof(int)*Ntot);
}


/* initialization */
void IsingFrame::initspin()
{
    int i, j, k, size;

    size = _NX*_NY*_NZ;
    
    Alloc();

    if (input[0] == -1)
    {
        memset(_S,-1,sizeof(int)*size);
    }
    else if (input[0] == 1)
    {
        memset(_S,1,sizeof(int)*size);
    }
    else if (input[0] == 0)
    {
        /* randomize */
        for(i=0;i<_NX;i++)
            for(j=0;j<_NY;j++)
                for(k=0;k<_NZ;k++)
                {
                    _S[i*_NY*_NZ+j*_NZ+k]=2*(drand48()>0.5)-1;
                }
    }
}

/********************************************************************/
/* Monte Carlo Simulation */
/********************************************************************/
void IsingFrame::MCrun()
{
#define MAXNN 6
    int i, j, k, ip, im, jp, jm, kp, km;
    int snn, acc, iter, key, Ntemp, Ntemp2;
    double Rp[2*MAXNN+1], Rm[2*MAXNN+1], R;
    int Sum; double avg_temp;
    int N_lgst_temp, umb_deltan;
    double umb_temp, kinetic_flag, kinetic_step0;
    int kinetic0;
    int Stot0;

    N_lgst_temp = 0;
    kinetic_flag = 0;
    kinetic0 = 0;
    Stot0 = 0;
    
    Sum=0;
    /* precompute flip rate */
    for(i=-MAXNN;i<=MAXNN;i++)
    {
        Rp[i+MAXNN] = exp(-2*(J*i+H)/kBT);
        Rm[i+MAXNN] = exp( 2*(J*i+H)/kBT);
    }

    INFO_Printf("Ntot = %d\n",Ntot);
    /* compute total spin */
    Stot = 0;
    for(i=0;i<Ntot;i++) Stot+=_S[i];

    /* continue step implmented Sep/01/2008, Seunghwa */
    if(continue_curstep) step0 = curstep;
    else step0 = 0;

//   printf("babo1\n");
    /* setup umbrella sampling arrays */
    if ( YES_UMB==1) {
	allocUMB(); 
	copyStoS0(); Stot0=Stot;
        calcrystalorder();
        N_lgst_temp = N_lgst_cluster;
        kinetic_flag=0;
    }
//    printf("babo2\n");
   
    /* FFS algorithm, works only if automatic end is initiated */
    if (FFSautoend==1) {
        Ntemp = (int) floor(0.5*Ntot);
        Ntemp2= (int) _Lam_array[FFSoption*2];
        if ( Ntemp < Ntemp2 ) {
            _Lam_array[FFSoption]=Ntemp;
        } else {
            _Lam_array[FFSoption]=Ntemp2;
        }
    } 

    for(curstep=step0;curstep<(step0+totalsteps);curstep++)
    { 
        key = 0;
        for(iter=0;iter<Ntot;iter++)
        {
            while(win!=NULL)
            {
                if(win->IsPaused()) sleep(1);
                else break;
            }

        /* simulation step */
        i = (int) floor(drand48()*_NX);
        j = (int) floor(drand48()*_NY);

        ip = (i+1)%_NX;
        im = (i-1+_NX)%_NX;
        jp = (j+1)%_NY;
        jm = (j-1+_NY)%_NY;

        if (_NZ==1)
        {
            snn = _S[ip*_NY+j]+_S[im*_NY+j]+_S[i*_NY+jp]+_S[i*_NY+jm];
            if (_S[i*_NY+j]>0)
                R = Rp[snn+MAXNN];
            else
                R = Rm[snn+MAXNN];
            acc = (drand48()<=R);
            if (acc)
            {
                _S[i*_NY+j] = -_S[i*_NY+j];
                Stot += 2*_S[i*_NY+j];
            }
        }
        else if (_NZ>1)
        {
            k = (int) floor(drand48()*_NZ);
            kp = (k+1)%_NZ;
            km = (k-1+_NZ)%_NZ;
            snn = _S[ip*_NY*_NZ+j*_NZ+k] + _S[im*_NY*_NZ+j*_NZ+k]
                + _S[i*_NY*_NZ+jp*_NZ+k] + _S[i*_NY*_NZ+jm*_NZ+k]
                + _S[i*_NY*_NZ+j*_NZ+kp] + _S[i*_NY*_NZ+j*_NZ+km];

            if (_S[i*_NY*_NZ+j*_NZ+k]>0)
                R = Rp[snn+MAXNN];
            else
                R = Rm[snn+MAXNN];
            acc = (drand48()<=R);
            if (acc)
            {
                _S[i*_NY*_NZ+j*_NZ+k] = -_S[i*_NY*_NZ+j*_NZ+k];
                Stot += 2*_S[i*_NY*_NZ+j*_NZ+k];
            }
        }
        UMB_curstep=curstep;
	if((curstep%calcryfreq==0)&&((iter+1)%calcryfreq1==0))
        {
	    calcrystalorder();

        if (YES_UMB==1 && UMB_K==0 ) {
//            printf("N_lgst_cluster=%d\n",N_lgst_cluster);
            if(N_lgst_cluster<Narray[0])
	    {
//            printf("op1 curstep=%d iter=%d\n",curstep,iter);
		copyS0toS(); Stot=Stot0;
		Occurence[0]++;
                Probability[0]++;
	    } 
            else if (N_lgst_cluster>Narray[2*delta_n])
            {
//            printf("op2\n");
                copyS0toS(); Stot=Stot0;
		Occurence[2*delta_n]++;
                Probability[2*delta_n]++;
            }
	    else
            {
//            printf("op3\n");
		umb_deltan=N_lgst_cluster-Narray[0];
		Occurence[umb_deltan]++;
                Probability[umb_deltan]++;
                copyStoS0(); Stot0=Stot;
                N_lgst_temp=N_lgst_cluster;
            }
	}

        if ( YES_UMB ==1 && UMB_K > 0 && kinetic_flag == 0 ) 
        {
	    umb_temp=(N_lgst_cluster-N_lgst_temp)*(N_lgst_cluster+N_lgst_temp-2*n_center);
	    R=exp(-0.5*UMB_K*(umb_temp)/kBT);
            acc = (drand48()<=R);
//            printf("N_cur = %d, N_pre = %d\n",N_lgst_cluster,N_lgst_temp);
//            printf("acc = %d, R = %f\n",acc,R);
            if (acc) {
                 umb_deltan=(N_lgst_cluster-Narray[0]);
                 umb_temp=exp(0.5*UMB_K*pow((N_lgst_cluster-n_center),2)/kBT);
                 copyStoS0(); Stot0=Stot;
                 if (umb_deltan>=0 && umb_deltan <= 2*delta_n && Kinetic < 1 && curstep > UMB_equilstep) {
		    Occurence[umb_deltan]++;
                    Probability[umb_deltan]+=umb_temp;
		}
                N_lgst_temp=N_lgst_cluster;
            }
            else 
            {
		umb_deltan=(N_lgst_temp-Narray[0]);
 	 	umb_temp=exp(0.5*UMB_K*pow((N_lgst_temp-n_center),2)/kBT); 
                copyS0toS(); Stot=Stot0;
                if (umb_deltan>=0 && umb_deltan <= 2*delta_n && Kinetic < 1 && curstep > UMB_equilstep ) {
		   Occurence[umb_deltan]++;
                   Probability[umb_deltan]+=umb_temp;
		}
            }
        }


        if ( Kinetic ==1 && kinetic_flag == 1 ) {
	    kinetic0++;
//	    printf("kinetic0=%d, N_clus=%d\n",kinetic0, N_lgst_cluster);
	    n_array[kinetic0]=n_array[kinetic0]+(N_lgst_cluster-n_center)*(N_lgst_cluster-n_center);
	    if ( kinetic0 > (int) 1.0*Ntot/calcryfreq1*Kinetic_Time ) kinetic_flag=0;
        }

        if ( Kinetic == 1 && N_lgst_cluster == n_center && kinetic_flag == 0 && curstep < totalsteps-(int)Kinetic_Time-1 )
  	{
//	printf("kinetic_initiated, curstep=%d\n N_clus=%d \n",curstep,N_lgst_cluster);
	    kinetic_step0=curstep;
	    kinetic_flag=1;
            Kinetic_Swip++;
	    kinetic0=0;
	}
        
        }

        if((curstep%printfreq==0)&&((iter+1)%Ntot==0))
        {
            INFO_Printf("curstep = %6d Stot = %6d (%6d/%6d up) Nc = %6d Check =%d\n",curstep,Stot,(Stot+Ntot)/2,Ntot,N_lgst_cluster,FFS0_check);
     	}

        if( saveFFScn == 1 && ((iter+1)%calcryfreq1==0) )
        {
	    FFScurstep=curstep;
            key=FFSprocess(); 
            if (FFSoption==0 && N_lgst_cluster > lambda_B)
		initspin();
        }
        
	if ((saveprop)&&(iter==Ntot-1))
	{
	    if((curstep%savepropfreq)==0) pf.write(this);
	}
            if (key==1) break;
        if((savecn)&&(iter==Ntot-1))
            if((curstep%savecnfreq)==0&&(curstep!=0))
                intercn.write(this);
        
        if(iter==Ntot-1) {
            Sum=Sum+N_lgst_cluster; 
            winplot();
        }
        }
        if (key==1) break;
    }
    avg_temp=Sum/totalsteps;
    N_lgst_avg=(int) avg_temp;
    if ( YES_UMB == 1 && Kinetic < 1) printhist();
    if ( Kinetic ==1 ) kineticprint();
}

/* compute energy of current config */
double IsingFrame::Energy()
{
    int i, j, k, ip, im, jp, jm, kp, km;
    int snn, snntot;
    
    Stot = 0;
    for(i=0;i<Ntot;i++)
        Stot+=_S[i];
    
    if(_NZ==1)
    {
        snntot = 0;
        for(i=0;i<_NX;i++)
          for(j=0;j<_NY;j++)
          {
              ip = (i+1)%_NX;
              im = (i-1+_NX)%_NX;
              jp = (j+1)%_NY;
              jm = (j-1+_NY)%_NY;
              
              snn = _S[ip*_NY+j]+_S[im*_NY+j]+_S[i*_NY+jp]+_S[i*_NY+jm];
              snntot += snn*_S[i*_NY+j];
          }
    }
    else
    {
        snntot = 0;
        for(i=0;i<_NX;i++)
          for(j=0;j<_NY;j++)
              for(k=0;k<_NZ;k++)
              {
                  ip = (i+1)%_NX;
                  im = (i-1+_NX)%_NX;
                  jp = (j+1)%_NY;
                  jm = (j-1+_NY)%_NY;
                  kp = (k+1)%_NZ;
                  km = (k-1+_NZ)%_NZ;
              
                  snn = _S[ip*_NY*_NZ+j*_NZ+k] + _S[im*_NY*_NZ+j*_NZ+k]
                      + _S[i*_NY*_NZ+jp*_NZ+k] + _S[i*_NY*_NZ+jm*_NZ+k]
                      + _S[i*_NY*_NZ+j*_NZ+kp] + _S[i*_NY*_NZ+j*_NZ+km];
                  
                  snntot += snn*_S[i*_NY*_NZ+j*_NZ+k];
          }
    }
    E = -snntot*J*0.5-Stot*H;
    return E;
}


/********************************************************************/
/* Path sampling */
/********************************************************************/
void Path::Alloc(int nc,int nx,int ny,int nz)
{
    maxlen = nc;
    _NX = nx; _NY = ny; _NZ = nz; 
    Ntot = _NX*_NY*_NZ;
    
    Realloc(_S0,int,Ntot);
    memset(_S0,-1,sizeof(int)*Ntot);

    Realloc(ix,int,maxlen);
    Realloc(iy,int,maxlen);
    Realloc(iz,int,maxlen);
    Realloc(newS,int,maxlen);
    Realloc(newStot,int,maxlen);    
}

void Path::SetS0(int *s)
{
    memcpy(_S0,s,sizeof(int)*Ntot);
}

void Path::Free()
{
    free(_S0);
    free(ix);
    free(iy);
    free(iz);
    free(newS);
    free(newStot);
}

void Path::Append(int i,int j,int k,int s)
{
    if(len==maxlen)
    {
        /* expand memory allocation */
        maxlen = maxlen + blocksize;
        Realloc(ix,int,maxlen);
        Realloc(iy,int,maxlen);
        Realloc(iz,int,maxlen);
        Realloc(newS,int,maxlen);
        Realloc(newStot,int,maxlen);    
    }
    else if(len>maxlen)
    {
        ERROR("Path: memory out of bound");
    }

    /* Append entry at the end of list */
    ix[len]=i;
    iy[len]=j;
    iz[len]=k;
    newS[len]=s;
    /* assuming that the new spin always flips sign */
    if(len>0) newStot[len]=newStot[len-1]+s*2;
    else newStot[len]=S0tot+s*2;

    len ++;
}

/* reading and writing paths */
int IsingFrame::readpath(int pn)
{
    class Path *p;
    int i, j, nc, nx, ny, nz, ci, cj, ck, cs, newStot;
    char *buffer; char *pp, *q; 
    char fname[300];

    p = NULL;
    if(pn==1)      p=&pathA;
    else if(pn==2) p=&pathB;
    else if(pn==3) p=&pathC;
    else ERROR("initpath: pn="<<pn<<" not allowed");

    /* read file: incnfile */
    LFile::SubHomeDir(incnfile,fname);
    LFile::LoadToString(fname,buffer,0);

    pp=buffer;
    
    q=pp; pp=strchr(pp,'\n'); if(pp) *(char *)pp++=0;
    sscanf(q, "%d", &nc);

    q=pp; pp=strchr(pp,'\n'); if(pp) *(char *)pp++=0;
    sscanf(q, "%d", &nx);

    q=pp; pp=strchr(pp,'\n'); if(pp) *(char *)pp++=0;
    sscanf(q, "%d", &ny);

    q=pp; pp=strchr(pp,'\n'); if(pp) *(char *)pp++=0;
    sscanf(q, "%d", &nz);

    if((nc>p->maxlen)|(nx!=p->_NX)|(ny!=p->_NY)|(nz!=p->_NZ))
        p->Alloc(nc,nx,ny,nz);

    /* read matrix _S0 */
    newStot = 0;
    for(i=0;i<p->Ntot;i++)
    {
        q=pp; pp=strchr(pp,'\n'); if(pp) *(char *)pp++=0;
        sscanf(q, "%d", p->_S0+i);
        newStot += p->_S0[i];
    }

    /* read incremental changes */
    p->len = 0;
    for(j=0;j<p->len;j++)
    {
        q=pp; pp=strchr(pp,'\n'); if(pp) *(char *)pp++=0;
        sscanf(q, "%d",&ci);
        q=pp; pp=strchr(pp,'\n'); if(pp) *(char *)pp++=0;
        sscanf(q, "%d",&cj);
        q=pp; pp=strchr(pp,'\n'); if(pp) *(char *)pp++=0;
        sscanf(q, "%d",&ck);
        q=pp; pp=strchr(pp,'\n'); if(pp) *(char *)pp++=0;
        sscanf(q, "%d",&cs);
        
        p->Append(ci,cj,ck,cs);

        newStot += p->newS[i]*2; /* assuming spin always flips in each step of path*/
        p->newStot[i] = newStot;        
    }
    Free(buffer);
    DUMP("readpath finished");
    return 0;
}


int IsingFrame::writepath(int pn)
{
    class Path *p;
    int i;
    FILE *fp;

    p = NULL;
    if(pn==1)      p=&pathA;
    else if(pn==2) p=&pathB;
    else if(pn==3) p=&pathC;
    else ERROR("writepath: pn="<<pn<<" not allowed");

    /* read file: finalcnfile */
    fp=fopen(finalcnfile,"w");

    fprintf(fp,"%d\n",p->len);
    fprintf(fp,"%d\n",p->_NX);
    fprintf(fp,"%d\n",p->_NY);
    fprintf(fp,"%d\n",p->_NZ);

    /* write spin matrix */
    for(i=0;i<p->Ntot;i++)
    {
        fprintf(fp,"%d\n", p->_S0[i]);
    }

    for(i=0;i<p->len;i++)
    {
        fprintf(fp,"%d\n%d\n%d\n%d\n",p->ix[i],p->iy[i],p->iz[i],p->newS[i]);
    }
    DUMP("writepath finished");
    return 0;
}

void IsingFrame::initpath(int pn)
{
    class Path *p;
    int i, j, k, ix, iy, iz, *a;
    
    p = NULL;
    if(pn==1)      p=&pathA;
    else if(pn==2) p=&pathB;
    else if(pn==3) p=&pathC;
    else ERROR("initpath: pn="<<pn<<" not allowed");

    /* create a rand path from all spin down (-1) to all spin up (+1) */
    p->Alloc(Ntot,_NX,_NY,_NZ);

    /* set initial spin to point down (-1) */
    memset(p->_S0,-1,sizeof(int)*Ntot);
    p->S0tot = -Ntot;
    
    /* create random permutation array */
    a = (int*) malloc(sizeof(int)*Ntot);
    for(i=0;i<Ntot;i++) a[i]=i;
    for(i=0;i<Ntot;i++)
    {
        j=(int)floor((Ntot-i)*drand48());
        k = a[i];
        a[i] = a[i+j];
        a[i+j] = k;
    }

    /* start builing the path from the beginning */
    p->len = 0;
    for(i=0;i<Ntot;i++)
    {
        iz = a[i]%_NZ;
        iy = ((a[i]-iz)/_NZ)%_NY;
        ix = (a[i]-iz-iy*_NZ)/_NX;
        p->Append(ix,iy,iz,1);
    }

    free(a);
}

void IsingFrame::copyPathtoS(int pn,int nstep)
{
    class Path *p;
    int i, id;

    p = NULL;
    if(pn==1)      p=&pathA;
    else if(pn==2) p=&pathB;
    else if(pn==3) p=&pathC;
    else ERROR("copyPathtoS: pn="<<pn<<" not allowed");

    if((_NX!=p->_NX)|(_NY!=p->_NY)|(_NZ!=p->_NZ))
        ERROR("copyPathtoS: dimension mismatch!");
    
    memcpy(_S,p->_S0,sizeof(int)*Ntot);

    for(i=0;i<nstep;i++)
    {
        id = (p->ix[i])*_NY*_NZ+(p->iy[i])*_NZ+(p->iz[i]);
        _S[id] = p->newS[i];
    }
    //winplot();
}


/* sampling of Monte Carlo path */
void IsingFrame::SampleMCpath(double smin, double smax,int nmax,int pn)
{  /* starting from _S (current spin) configuartion
    *  do MC simulation, record spin flips into path
    *  terminate if fraction of spin up < smin, or > smax
    *  or if maximum iteration exceeds nmax
    */
    class Path *p;
    int i, j, ip, im, jp, jm, snn, acc, iter;
    double Rp[9], Rm[9], R, frac_spin_up;

    p = NULL;
    if(pn==1)      p=&pathA;
    else if(pn==2) p=&pathB;
    else if(pn==3) p=&pathC;
    else ERROR("SampleMCpath: pn="<<pn<<" not allowed");

    //INFO("SampleMCpath");
    
    p->SetS0(_S);
    p->len = 0;

    /* precompute flip rate */
    for(i=-4;i<=4;i++)
    {
        Rp[i+4] = exp(-2*J*(i+H)/kBT);
        Rm[i+4] = exp( 2*J*(i+H)/kBT);
    }

    /* compute total spin */
    Stot = 0;
    for(i=0;i<Ntot;i++)
        Stot+=_S[i];
    p->S0tot = Stot;
    
     for(iter=0;iter<nmax;iter++)
     {
         while(win!=NULL)
         {
             if(win->IsPaused()) sleep(1);
             else break;
         }

         /* simulation step */
         i = (int) floor(drand48()*_NX);
         j = (int) floor(drand48()*_NY);

         ip = (i+1)%_NX;
         im = (i-1+_NX)%_NX;
         jp = (j+1)%_NY;
         jm = (j-1+_NY)%_NY;
         
         snn = _S[ip*_NY+j]+_S[im*_NY+j]+_S[i*_NY+jp]+_S[i*_NY+jm];
         if (_S[i*_NY+j]>0)
             R = Rp[snn+4];
         else
             R = Rm[snn+4];
         
         acc = (drand48()<=R);
         if (acc)
         {
             _S[i*_NY+j] = -_S[i*_NY+j];
             Stot += 2*_S[i*_NY+j];
             p->Append(i,j,0,_S[i*_NY+j]);
         }

         frac_spin_up = (Stot+Ntot)/2.0/Ntot;

         if((frac_spin_up<smin)|(frac_spin_up>smax))
             break;
         
//         if(iter%printfreq==0)
//            INFO_Printf("iter = %6d Stot = %6d (%d/%d up)\n",
//                        iter,Stot,(Stot+Ntot)/2,Ntot);

         //if(iter%plotfreq==0) winplot();
     }
//     INFO_Printf("SampleMCpath: Stot %d - %d\n",Stot,p->newStot[p->len-1]);     
//     INFO_Printf("SampleMCpath: num spin up %d\n",(Stot+Ntot)/2);     
}

void IsingFrame::copypath(int pdes, int psrc)
{
    /* copy path psrc to pdes */
    class Path *pd, *ps;

    pd = NULL;
    if(pdes==1)      pd=&pathA;
    else if(pdes==2) pd=&pathB;
    else if(pdes==3) pd=&pathC;
    else ERROR("copypath: pdes="<<pdes<<" not allowed");
    
    ps = NULL;
    if(psrc==1)      ps=&pathA;
    else if(psrc==2) ps=&pathB;
    else if(psrc==3) ps=&pathC;
    else ERROR("copypath: psrc="<<psrc<<" not allowed");

    pd->Alloc(ps->maxlen,ps->_NX,ps->_NY,ps->_NZ);

    pd->len = ps->len;
    pd->Ntot = ps->Ntot;
    pd->S0tot = ps->S0tot;
    
    memcpy(pd->_S0,ps->_S0,sizeof(int)*(ps->Ntot));
    memcpy(pd->ix,ps->ix,sizeof(int)*(ps->maxlen));
    memcpy(pd->iy,ps->iy,sizeof(int)*(ps->maxlen));
    memcpy(pd->iz,ps->iz,sizeof(int)*(ps->maxlen));
    memcpy(pd->newS,ps->newS,sizeof(int)*(ps->maxlen));
    memcpy(pd->newStot,ps->newStot,sizeof(int)*(ps->maxlen));

}

    
void IsingFrame::cutpath(int pn, int istart, int newlen)
{
    /* only retain segments from istart to iend-1 */
    class Path *p;
    int i, id;
    
    p = NULL;
    if(pn==1)      p=&pathA;
    else if(pn==2) p=&pathB;
    else if(pn==3) p=&pathC;
    else ERROR("cutpath: pn="<<pn<<" not allowed");

    if((istart<0)||(newlen<0)||(newlen>p->len))
        ERROR("cutpath: istart="<<istart<<"  newlen="<<newlen);

    /* walk along the path to istart-1 */
    for(i=0;i<istart;i++)
    {
        id = (p->ix[i])*_NY*_NZ+(p->iy[i])*_NZ+(p->iz[i]);
        p->_S0[id] = p->newS[i];
        p->S0tot = p->newStot[i];
    }

    /* move path segment up front */
    memmove(p->ix,p->ix+istart,sizeof(int)*(newlen));
    memmove(p->iy,p->iy+istart,sizeof(int)*(newlen));
    memmove(p->iz,p->iz+istart,sizeof(int)*(newlen));
    memmove(p->newS,p->newS+istart,sizeof(int)*(newlen));
    memmove(p->newStot,p->newStot+istart,sizeof(int)*(newlen));
            

    memset(p->ix+newlen,0,sizeof(int)*(p->maxlen-newlen));
    memset(p->iy+newlen,0,sizeof(int)*(p->maxlen-newlen));
    memset(p->iz+newlen,0,sizeof(int)*(p->maxlen-newlen));
    memset(p->newS+newlen,0,sizeof(int)*(p->maxlen-newlen));
    memset(p->newStot+newlen,0,sizeof(int)*(p->maxlen-newlen));
    
    p->len = newlen;
}

void IsingFrame::mergepath(int phead, int ptail)
{
    /* merge ptail to the end of phead */
    class Path *ph, *pt;
    int i;
    
    ph = NULL;
    if(phead==1)      ph=&pathA;
    else if(phead==2) ph=&pathB;
    else if(phead==3) ph=&pathC;
    else ERROR("mergepath: phead="<<phead<<" not allowed");
    
    pt = NULL;
    if(ptail==1)      pt=&pathA;
    else if(ptail==2) pt=&pathB;
    else if(ptail==3) pt=&pathC;
    else ERROR("mergepath: ptail="<<ptail<<" not allowed");

    if((ph->_NX!=pt->_NX)||(ph->_NY!=pt->_NY)||(ph->_NZ!=pt->_NZ))
        ERROR("mergepath: dimension mismatch");

    if(ph->newStot[ph->len-1]!=pt->S0tot)
    {
        INFO_Printf("ph Stot_end = %d pt S0tot = %d",ph->newStot[ph->len-1],pt->S0tot);
        ERROR("mergepath: total spin mismatch");
    }
    for(i=0;i<pt->len;i++)
        ph->Append(pt->ix[i],pt->iy[i],pt->iz[i],pt->newS[i]);

}

void IsingFrame::reversepath(int pn)
{
    /* reverse direction of path */
    class Path *p;
    int i, k, id, len;

    p = NULL;
    if(pn==1)      p=&pathA;
    else if(pn==2) p=&pathB;
    else if(pn==3) p=&pathC;
    else ERROR("cutpath: pn="<<pn<<" not allowed");
    
    /* walk along the path to the end */
    len = p->len;
    for(i=0;i<len;i++)
    {
        id = (p->ix[i])*_NY*_NZ+(p->iy[i])*_NZ+(p->iz[i]);
        p->_S0[id] = p->newS[i];
        p->S0tot = p->newStot[i];        
    }

    //k=0;
    //for(i=0;i<p->Ntot;i++) k+=p->_S0[i];
    //INFO_Printf("reversepath: S0 = %d %d\n",k,p->S0tot);
    
    /* swap the arrays */
    for(i=0;i<len/2;i++)
    {
        k=p->ix[i]; p->ix[i]=p->ix[len-i-1]; p->ix[len-i-1]=k;
        k=p->iy[i]; p->iy[i]=p->iy[len-i-1]; p->iy[len-i-1]=k;
        k=p->iz[i]; p->iz[i]=p->iz[len-i-1]; p->iz[len-i-1]=k;
        k=p->newS[i]; p->newS[i]=p->newS[len-i-1]; p->newS[len-i-1]=k;
    }

    /* change the sign of newS */
    for(i=0;i<len;i++)
        p->newS[i] = -p->newS[i];

    /* update newStot */
    if(len>0) p->newStot[0]=p->S0tot+p->newS[0]*2;
    for(i=1;i<len;i++)
    {
        p->newStot[i]=p->newStot[i-1]+p->newS[i]*2;
    }
    //INFO_Printf("reversepath: Stot_end = %d\n",p->newStot[len-1]);
}

int IsingFrame::ComputeSuccessRate(int nsample,double smin,double smax,int nmax)
{
    int nsucc, isample;
    double frac_spin_up;

    copyStoS0();
    nsucc = 0;

    for(isample=0;isample<nsample;isample++)
    {
            
        SampleMCpath(smin, smax, nmax, 3); /* resample path C starting from S */
        frac_spin_up = (pathC.newStot[pathC.len-1]+Ntot)/2.0/Ntot;
            
        if(frac_spin_up<=smin)
        {
            /* failure */
        }
        else if(frac_spin_up>=smax)
        {
            /* success */
            nsucc ++;
        }
        else
        {
            /* cannot tell */
            WARNING("ComputeSuccessRate: nmax too short to tell");
        }
        copyS0toS();
     }
     INFO_Printf("ComputeSuccessRate: nsample = %d nsuccess = %d sratio = %.5f\n",
                 nsample,nsucc,(1.0*nsucc)/nsample);
     return nsucc;
}
    
void IsingFrame::AnalyzeConfigs(int nstart,int nend,int nsample,double smin, double smax,int nmax)
{
    int iconfig, isample, nsucc;
    double frac_spin_up;
    
    initpath(3); /* initialize path C */
    for(iconfig=nstart;iconfig<=nend;iconfig++)
    {
        /* read in cn (spin) file */
        sprintf(incnfile,"inter%04d.cn",iconfig);
        readcn();  copyStoS0();
        Energy();
        INFO_Printf("Analyzing %s  E = %e  Stot = %d Nup = %d\n",incnfile,E,Stot,(Stot+Ntot)/2);

        curstep = 0; winplot();
        
        nsucc = 0;
        for(isample=0;isample<nsample;isample++)
        {
            /* read in cn (spin) file again */
            if(isample>0) copyS0toS();
            
            SampleMCpath(smin, smax, nmax, 3); /* resample path C starting from S */
            frac_spin_up = (pathC.newStot[pathC.len-1]+Ntot)/2.0/Ntot;
            
            if(frac_spin_up<=smin)
            {
                /* failure */
            }
            else if(frac_spin_up>=smax)
            {
                /* success */
                nsucc ++;
            }
            else
            {
                /* cannot tell */
                WARNING("AnalyzeConfig: nmax too short to tell");
            }
            if(isample%printfreq==0) INFO_Printf("nsucc = %d / %d\n",nsucc,isample+1);
        }
        INFO_Printf("Analyzing %s  E = %e  Stot = %d nsucc = %d nsample = %d sratio = %.5f\n",
                    incnfile,E,Stot,nsucc,nsample,(1.0*nsucc)/nsample);
    }
}


void IsingFrame::WalkonChain(double smin,double smax,int nmax)
{
    int len, slide, dslide;
    double frac_spin_up;
    
    /* randomize path A */
    initpath(1);
    initpath(2);
    initpath(3);

    len = pathA.len;
    slide = (int) floor(drand48()*len);
    for(curstep=0;curstep<totalsteps;curstep++)
    {
        //slide = 3500;
        //slide = 2000;
    

        copyPathtoS(1, slide);             /* move S along path A */
        winplot();                 /* plot S before path sampling */
        if(savecn && (curstep%savecnfreq==0) && (curstep!=0) )
        {
            _NSUCCESS = ComputeSuccessRate(_NSAMPLE,smin,smax,nmax);
            intercn.write(this);   /* save saddle conf as intercn */
        }

        SampleMCpath(smin, smax, nmax, 3); /* resample path C */
        //winplot();                 /* plot S after path sampling */

        copypath(2, 1);               /* save path A to path B */
    
        frac_spin_up = (pathC.newStot[pathC.len-1]+Ntot)/2.0/Ntot;

        //INFO_Printf("Num spin up = %d\n",(pathC.newStot[pathC.len-1]+Ntot)/2);
    
        INFO_Printf("WalkonChain: curstep = %6d slide = %6d  len = %7d",
                    curstep,slide,len);
        if(frac_spin_up<=smin)
        { /* reach state A (failure) */
            //INFO("spin down");
            INFO_Printf(" (-) down %.3f\n",frac_spin_up);
            reversepath(1);           /* reverse path A */
            cutpath(1, 0, len-slide); /* truncate path A */
            mergepath(1, 3);          /* merge path A and C */
            reversepath(1);           /* reverse path A */
            slide = pathC.len;
            //dslide = (int) floor(drand48()*(len-slide)*0.3);
            dslide = (int) floor(drand48()*Ntot);
            slide += dslide;
            len = pathA.len;
            if(slide>len) slide-=dslide;
        }
        else if(frac_spin_up>=smax)
        { /* reach state B (success) */
            //INFO("spin up");
            INFO_Printf(" (+) up   %.3f\n",frac_spin_up);            
            cutpath(1, 0, slide);  /* truncate path A from 0 to slide */
            mergepath(1, 3);       /* merge path A and C */
            //slide -= (int) floor(drand48()*slide*0.3);
            dslide = (int) floor(drand48()*Ntot);
            slide -= dslide;
            len = pathA.len;
            if(slide<0) slide+=dslide;
        }
        else
        {
            //ERROR("WalkonChain: frac_spin_up="<<frac_spin_up);
            /* maximum iteration exceeded, try again */
            INFO_Printf(" try again %.3f\n",frac_spin_up);            
            slide = (int) floor(drand48()*len);
        }

        if(savecn && (curstep%savecnfreq==0) && (curstep!=0) )
        {
            INFO("save inter-path-A.pcn");
            strcpy(finalcnfile,"inter-path-A.pcn");
            writepath(1);
        }
    }       
}


/********************************************************************/
/* Forward Flux Sampling (FFS) algorithm                            */
/********************************************************************/
int IsingFrame::FFSprocess()
{
    int i;
    if (FFSoption==0) 
    {
	if (N_lgst_cluster >= _Lam_array[0] && FFS0_check == 0 )
	{
	    FFScn.write(this,zipfiles,true);
	    FFS0_check=1;
	}
	if (FFS0_check==1 && N_lgst_cluster <= lambda_A )
	    FFS0_check=0;
        return 0;
    }
    else
    {
	if (N_lgst_cluster >= _Lam_array[FFSoption])
	{
	    FFS0_check = 1;
	    FFScn.write(this,zipfiles,true);
	    return 1;
	}
	
	if ( FFSpruning > 0 )
	{
	    for (i=FFSoption;i>1;i--)
	    {
		if ( N_lgst_cluster<=_Lam_array[i-2] && _Lam_check[i-2]==0 )
		{
		    if (drand48()>FFS_Pp)
		    {
			FFScn_weight=FFScn_weight/(1-FFS_Pp);
			_Lam_check[i-2]=1;
		    }
		    else
		    {
		 	FFS0_check = -1;
		    }
		}
		if ( N_lgst_cluster <= lambda_A ) FFS0_check=-1;
	    }
	    if ( N_lgst_cluster <= lambda_A ) FFS0_check=-1;
	}
	else
	{
	    if ( FFScurstep > step0 + totalsteps || N_lgst_cluster <= lambda_A )
		FFS0_check = -1;
	}
   	
	if (FFS0_check<0) return 1;
	else return 0;
    }
}

void IsingFrame::calcrystalorder()
{
    int i, Nc_i,j,k0;
    int N_cluster, Size_cluster[Ntot];
    N_cluster = 0;
    N_lgst_cluster = 0;

    for (i=0;i<Ntot;i++)
    {
	DFS_mark[i]=0;
        Size_cluster[i]=0;
    }
    k0=0;
    for (i=0;i<Ntot;i++)
    {
	if (_S[i]>0)
	{
	    N_cluster++;
	    N_cluster_temp=0;
	    if (N_cluster == 1)
	        Nc_i=i;
	    DFS(i);
	    if (N_cluster_temp > N_lgst_cluster)
	    {
	   	N_lgst_cluster = N_cluster_temp;
	   	Nc_i=i;
	    }
/************************UMBs**/
	    if (YES_UMB==1 && UMB_K < 0 && N_cluster_temp > 0 )
            {
	        for (j=0;j<=2*delta_n;j++) if (Narray[j]==N_cluster_temp) {
                if ( UMB_curstep > UMB_equilstep ) { 
		Occurence[j]++;
		Probability[j]++;
		}
		}
//                if (N_cluster_temp>Narray[2*delta_n]) printf("n_clu=%d\n",N_cluster_temp);
	    }
/************************UMBe**/
	}
/*************************UMBs**/
	else 
 	{ 
            k0++;
        }
/*************************UMBe**/
    }
//    printf("Ntot=%d, S_=%d",Ntot,k0);
/*************************UMBs**/
    if (YES_UMB == 1 && UMB_K < 0 ) 
    { 
	for (j=0;j<=2*delta_n;j++) if (Narray[j]==0) {
	if (UMB_curstep > UMB_equilstep ) {
		Occurence[j]=Occurence[j]+k0; 
		Probability[j]=Probability[j]+k0;
	}
	}
    } 
/*************************UMBe**/
}

void IsingFrame::allocDFS()
{
    Realloc(DFS_mark,int,Ntot);
}

void IsingFrame::allocUMB()
{
    int i, no;
    no= (int) round(1.0*Ntot/calcryfreq1*Kinetic_Time);
    Realloc(Narray,int,2*delta_n+1);
    Realloc(Occurence,int,2*delta_n+1);
    Realloc(Probability,double,2*delta_n+1);
    Realloc(time_array,double,no+1);
    Realloc(n_array,int,no+1);
    for (i=0;i<2*delta_n+1;i++)
    {
	Narray[i]=n_center-delta_n+i;
        Occurence[i]=0;
        Probability[i]=0;
    }
    for (i=0;i<no+1;i++)
    {
	time_array[i]=1.0*i*calcryfreq1/Ntot;
	n_array[i]=0;
    }
} 

void IsingFrame::printhist()
{
    FILE *fp, *fp1, *fp2;
    int i;
    fp  = fopen("Narray.txt","w");
    fp1 = fopen("Freq.txt","w");
    fp2 = fopen("Prob.txt","w");
//    printf("n_center=%d, delta_n=%d\n",n_center,delta_n);
    for (i=0;i<2*delta_n+1;i++)
    {
       fprintf(fp,"%d\n",Narray[i]);
       if (UMB_K >= 0 ) {
       fprintf(fp1,"%d\n",Occurence[i]);
       fprintf(fp2,"%f\n",Probability[i]);
       } else {
       fprintf(fp1,"%f\n",Probability[i]);
       fprintf(fp2,"%f\n",Probability[i]);
       }
    }
    fclose(fp);
    fclose(fp1);
    fclose(fp2);
}
void IsingFrame::kineticprint()
{
    FILE *fp, *fp1, *fp2;
    int i;
    fp = fopen("time_array.txt","w");
    fp1= fopen("n_array.txt","w");
    fp2= fopen("tot_swip.txt","w");
    for (i=0;i<(int)round(1.0*Ntot/calcryfreq1*Kinetic_Time+1);i++) {
	fprintf(fp,"%f\n",time_array[i]);
	fprintf(fp1,"%d\n",n_array[i]);
    }
	fprintf(fp2,"%d %d\n",Ntot,Kinetic_Swip);
    fclose(fp);fclose(fp1);fclose(fp2);
}

void IsingFrame::assign_Lam()
{
    int i, n;
    n = (int) round(input[0]);
    Realloc(_Lam_array,int,n);
    Realloc(_Lam_check,int,n);
    for (i=1;i<=n;i++)
    {
	_Lam_array[i-1]=(int)round(input[i]);
	_Lam_check[i-1]=0;
    }
    FFS0_check=1;
}

void IsingFrame::DFS(int i)
{
    int j, ip, im, jp, jm, kp, km, ind[6];
    int icur, jcur, kcur;
    div_t divresult;

    if (DFS_mark[i]==1 || _S[i]<0)
	return;
    else
    {
	if (_NZ<2) 
	{
		DFS_mark[i]=1;
		N_cluster_temp++;
	        divresult = div(i,_NY);
	        icur = divresult.quot;
	        jcur = divresult.rem;
	        ip = (icur+1)%_NX;
		im = (icur-1+_NX)%_NX;
		jp = (jcur+1)%_NY;
		jm = (jcur-1+_NY)%_NY;
	        ind[0]=ip*_NY+jcur;
	        ind[1]=im*_NY+jcur;
	        ind[2]=icur*_NY+jp;
	        ind[3]=icur*_NY+jm;
	        for (j=0;j<4;j++)
		{
		    if (_S[ind[j]] < 0 || DFS_mark[ind[j]] ==1 )
			continue;
		    DFS(ind[j]);
	        }
	}
	else
	{
		DFS_mark[i]=1;
		N_cluster_temp++;
		divresult = div(i,_NZ);
		icur = divresult.quot;
		kcur = divresult.rem;
		divresult = div(icur,_NY);
		icur = divresult.quot;
		jcur = divresult.rem;
		ip = (icur+1)%_NX;
		im = (icur-1+_NZ)%_NX;
		jp = (jcur+1)%_NY;
		jm = (jcur-1+_NY)%_NY;
		kp = (kcur+1)%_NZ;
		km = (kcur-1+_NZ)%_NZ; 
		ind[0]=ip*_NY*_NZ+jcur*_NY+kcur;
		ind[1]=im*_NY*_NZ+jcur*_NY+kcur;
		ind[2]=icur*_NY*_NZ+jp*_NY+kcur;
		ind[3]=icur*_NY*_NZ+jm*_NY+kcur;
		ind[4]=icur*_NY*_NZ+jcur*_NY+kp;
		ind[5]=icur*_NY*_NZ+jcur*_NY+km;
		for (j=0;j<6;j++)
		{
			if (_S[ind[j]] < 0 || DFS_mark[ind[j]]==1)
				continue;
			DFS(ind[j]);
		}
	}
 
    }

}

/********************************************************************/
/* File input and output */
/********************************************************************/
int IsingFrame::readcn()
{
    initcn.open(incnfile,LFile::O_Read);
    return initcn.read(this);
}

int IsingFrame::setfilecounter()
{
    intercn.setcount(filecounter);
    return 0;
}

int IsingFrame::setFFSfilecounter()
{
    FFScn.setcount(FFSfilecounter);
    return 0;
}

int IsingFrame::writefinalcnfile(int zip, bool bg)
{
    finalcn.open(finalcnfile);
    return finalcn.write(this,zip,bg);
}

void IsingFrame::saveintercn(int step)
{
    if(savecn)
        if((step%savecnfreq)==0&&(step!=0))
            writeintercnfile();
}

int IsingFrame::writeintercnfile(int zip,bool bg)
{
    return intercn.write(this,zip,bg);
}

int IsingFrame::writeFFScnfile(int zip, bool bg)
{
    return FFScn.write(this,zip,bg);
}

int IsingFrame::openintercnfile()
{
    return intercn.open(intercnfile);
}

int IsingFrame::openFFScnfile()
{
    return FFScn.open(FFScnfile);
}

int IsingFrame::openpropfile()
{
    return pf.open(outpropfile);
}

/* Configuration Files */
char * CNFile::describe()
{
    static char tmp[500];
    sprintf(tmp,"%s","spin configurations, Format:\n"
            "_NX _NY _NZ; s(1,1,1); s(1,1,2); ...; s(_NX,_NY,_NZ)");
    return tmp;
}

int CNFile::writeblock(void *p)
{
    int i;
    IsingFrame &d=*((IsingFrame *)p);
    f->Printf("%10d\n%10d\n%10d\n",d._NX,d._NY,d._NZ);
    for(i=0;i<d.Ntot;i++)
    {
        f->Printf("%3d\n",d._S[i]);
    }
    f->Printf("%10d\n",d._NSAMPLE);
    f->Printf("%10d\n",d._NSUCCESS);
    return 0;
}

int CNFile::readblock(void *p)
{
    int i;
    char *buffer, *pp, *q;
    IsingFrame &d=*((IsingFrame *)p);
    
    LFile::LoadToString(fname,buffer,0);

    pp=buffer;
    sscanf(pp, "%d %d %d", &d._NX, &d._NY, &d._NZ);
    INFO("readblock: NX="<<d._NX<<" NY="<<d._NY<<" NZ="<<d._NZ);    
    
    d.Alloc();
    
    pp=strchr(pp, '\n');
    pp++;
    pp=strchr(pp, '\n');
    pp++;
    pp=strchr(pp, '\n');
    pp++;

    d.Stot = 0;
    for(i=0;i<d.Ntot;i++)
    {
        q=pp;
        pp=strchr(pp,'\n');
        if(pp) *(char *)pp++=0;
        sscanf(q, "%d",&(d._S[i]));
        d.Stot += d._S[i];
    }
    q=pp; pp=strchr(pp,'\n'); 
    if(pp) 
    {
       *(char *)pp++=0;
       sscanf(q, "%d",&(d._NSAMPLE));
       q=pp; pp=strchr(pp,'\n'); 
       if(pp)
       {
          *(char *)pp++=0;
          sscanf(q, "%d",&(d._NSUCCESS));
       }
    }

    Free(buffer);
    DUMP("readblock finished");
    return 0;
}

/* Output Property Files */
char * PropFile::describe()
{
    static char tmp[500];
    sprintf(tmp,"%s","Format:\n"
            "curstep total_spin num_spin_up ");
    return tmp;
}

int PropFile::writeentry(void *p)
{
    char tmp[500];
    IsingFrame &d=*((IsingFrame *)p);
    sprintf(tmp,"%10d %10d %10.0f %10d %10d %2d\n",
            d.curstep,d.Stot,(d.Stot+d.Ntot)/2.0,d.Ntot,d.N_lgst_cluster,d.FFS0_check);
    *f<<tmp;
    return 0;
}

/********************************************************************/
/* Visualization */
/********************************************************************/

void IsingFrame::openwin()
{
    openwindow(win_width,win_height,dirname);
}

int IsingFrame::openwindow(int w,int h,const char *n)
{
    if(win!=NULL)
        if(win->alive) return -1;
    win=new YWindow(w,h,n,true,true,false);
    if(win==NULL) return -1;
    else
    {
        win->setinterval(100);
        win->Run();
        return 0;
    }
}

void IsingFrame::closewindow() {delete(win);}

void IsingFrame::winplot()
{
    if(win!=NULL)
        if(win->alive)
            if((curstep%plotfreq)==0)
            {
                plot();
            }
}

void IsingFrame::winplot(int step)
{
    if(win!=NULL)
        if(win->alive)
            if((step%plotfreq)==0)
            {
                plot();
            }
}        

void IsingFrame::wintogglepause()
{
    if(win!=NULL)
        if(win->IsAlive())
            win->TogglePause();
}

#include "namecolor.c"
void IsingFrame::alloccolors()
{
    int r,g,b;
    if(win!=NULL)
    {
        if(win->alive)
        {
            for(int i=0;i<MAXCOLORS;i++)
            {
                Str2RGB(colornames[i],&r,&g,&b);
                colors[i]=win->AllocShortRGBColor(r,g,b);
            }
            
            Str2RGB(atomcolor[0],&r,&g,&b);
            colors[MAXCOLORS+0]=win->AllocShortRGBColor(r,g,b);
            Str2RGB(bondcolor,&r,&g,&b);
            colors[MAXCOLORS+1]=win->AllocShortRGBColor(r,g,b);
            Str2RGB(highlightcolor,&r,&g,&b);
            colors[MAXCOLORS+2]=win->AllocShortRGBColor(r,g,b);
            Str2RGB(fixatomcolor,&r,&g,&b);
            colors[MAXCOLORS+3]=win->AllocShortRGBColor(r,g,b);
            Str2RGB(backgroundcolor,&r,&g,&b);
            win->bgcolor=win->AllocShortRGBColor(r,g,b);

            for(int i=0;i<MAXSPECIES;i++)
            {
                Str2RGB(atomcolor[i],&r,&g,&b);
                colors[MAXCOLORS+5+i]=win->AllocShortRGBColor(r,g,b);
            }
        }
        else
        {
            WARNING("No window to allocate color for!");
        }
    }
    else
    {
        WARNING("No window to allocate color for!");
    }
}
void IsingFrame::alloccolorsX()
{
    if(win!=NULL)
    {
        if(win->alive)
        {
            colors[0]=win->AllocNamedColor(atomcolor[0]);
            colors[1]=win->AllocNamedColor(bondcolor);
            colors[2]=win->AllocNamedColor(highlightcolor);
            colors[3]=win->AllocNamedColor(fixatomcolor);
            win->bgcolor=win->AllocNamedColor(backgroundcolor);
        }
        else
        {
            WARNING("No window to allocate color for!");
        }
    }
    else
    {
        WARNING("No window to allocate color for!");
    }
}

void IsingFrame::rotate()
{
    if(win!=NULL)
    {
        if(win->alive)
        {
            win->horizontalRot(rotateangles[0]*M_PI/180);
            win->verticalRot(rotateangles[1]*M_PI/180);
            win->spinRot(rotateangles[2]*M_PI/180);
            if(rotateangles[3]!=0) win->zoom(rotateangles[3]);
            if(rotateangles[4]!=0) win->project(rotateangles[4]);
        }
        else
        {
            WARNING("No window to rotate for!");
        }
    }
    else
    {
        WARNING("No window to rotate for!");
    }
}

void IsingFrame::saverot()
{
    if(win!=NULL)
    {
        if(win->alive)
        {
            //win->saveRot(); win->saveScale();
            win->saveView();
        }
        else
        {
            WARNING("No window to rotate for!");
        }
    }
    else
    {
        WARNING("No window to rotate for!");
    }
}


//#include "colormap.h"
  
void IsingFrame::plot()
{
    int i, j, k, spin;
    unsigned long c;
    
    if(win==NULL) return;
    if(!(win->alive)) return;
    
    win->Lock();
    win->Clear();
    
    /* draw atom s*/
    for(i=0;i<_NX;i++)
        for(j=0;j<_NY;j++)
            for(k=0;k<_NZ;k++)
            {
                spin = _S[i*_NY*_NZ+j*_NZ+k];

                if((_NZ>1)&&(spin<0)) continue;
                
                if(spin<0) c=colors[0]; else c=colors[1];

                win->DrawPoint((2.0*i)/_NX-1,1-(2.0*j)/_NY,(2.0*k)/_NZ-1,
                               atomradius[0]/_NX,c);
            }
    
    win->Unlock();
    win->Refresh();

}




/***********************************************************/
/*  Main Program Begins */

#ifdef _TEST

class IsingFrame sim;

/* The main program is defined here */
#include "main.cpp"

#endif//_TEST
