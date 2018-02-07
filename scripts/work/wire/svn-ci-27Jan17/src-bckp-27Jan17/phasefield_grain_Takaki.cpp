// Following Fan 1997 formulation 
// Last Modified : Mon March 31, 2014
// For grain coarsening phase field simulation

#include "phasefield_grain_Takaki.h"

void PhaseFieldFrame::initparser()
{
    char s[100];
    int i, j;

    /* Parser */
    bindvar("myname",myname,STRING);
    bindvar("command",command,STRING);
    bindvar("input",input,DOUBLE);
    bindvar("NX",&_NX,INT);
    bindvar("NY",&_NY,INT);
    bindvar("NZ",&_NZ,INT);
    bindvar("num_fields",&num_fields,INT);
    bindvar("ROT_MATRIX",&(_ROT_MATRIX[0][0]),DOUBLE);

    /* Simulation setting */
    bindvar("totalsteps",&totalsteps,INT);
    bindvar("savepropfreq",&savepropfreq,INT);
    bindvar("randseed",&_RANDSEED,INT);
    bindvar("gridsize",&gridsize,DOUBLE);
    bindvar("timestep",&timestep,DOUBLE);
    bindvar("dtmax",&dtmax,DOUBLE);
    bindvar("dphimax",&dphimax,DOUBLE);
    bindvar("alpha",&alpha,DOUBLE);
    bindvar("beta",&beta,DOUBLE);
    bindvar("gamma",&gamma,DOUBLE);
    bindvar("kappa",&kappa,DOUBLE);
    bindvar("L",&L,DOUBLE);
    /* Simulation variables */
    bindvar("F",&_F,DOUBLE);
    bindvar("G",&_G,DOUBLE);

    /* File input and output */
    bindvar("savecn",&savecn,INT);
    bindvar("saveprop",&saveprop,INT);
    bindvar("savecnfreq",&savecnfreq,INT);
    bindvar("savepropfreq",&savepropfreq,INT);
    bindvar("printfreq",&printfreq,INT);
    bindvar("filecounter",&filecounter,INT);
    bindvar("incnfile",incnfile,STRING);
    bindvar("finalcnfile",finalcnfile,STRING);
    bindvar("outpropfile",outpropfile,STRING);
    bindvar("intercnfile",intercnfile,STRING);
    bindvar("zipfiles",&zipfiles,INT);

    /* Visualization */
    bindvar("win_width",&win_width,INT);
    bindvar("win_height",&win_height,INT);
    bindvar("plotfreq",&plotfreq,INT);
    bindvar("atomradius",atomradius,DOUBLE);
    bindvar("plot_threshold",&plot_threshold,DOUBLE);

    bindvar("backgroundcolor",backgroundcolor,STRING);
    bindvar("rotateangles",rotateangles,DOUBLE);

    bindvar("ROT_MATRIX_11",&(_ROT_MATRIX[0][0]),DOUBLE);
    bindvar("ROT_MATRIX_12",&(_ROT_MATRIX[0][1]),DOUBLE);
    bindvar("ROT_MATRIX_13",&(_ROT_MATRIX[0][2]),DOUBLE);
    bindvar("ROT_MATRIX_21",&(_ROT_MATRIX[1][0]),DOUBLE);
    bindvar("ROT_MATRIX_22",&(_ROT_MATRIX[1][1]),DOUBLE);
    bindvar("ROT_MATRIX_23",&(_ROT_MATRIX[1][2]),DOUBLE);
    bindvar("ROT_MATRIX_31",&(_ROT_MATRIX[2][0]),DOUBLE);
    bindvar("ROT_MATRIX_32",&(_ROT_MATRIX[2][1]),DOUBLE);
    bindvar("ROT_MATRIX_33",&(_ROT_MATRIX[2][2]),DOUBLE);

    for(i=0;i<MAXNUMFIELDS;i++)
    {
        sprintf(s,"atomcolor%d",i);
        bindvar(s,atomcolor[i],STRING);
    }
    bindvar("atomcolor",atomcolor[0],STRING);
    for(i=0;i<MAXCOLORS;i++)
    {
        sprintf(s,"color%02d",i);
        bindvar(s,colornames[i],STRING);
    }
}    

int PhaseFieldFrame::exec(const char *name)
{
    if(Organizer::exec(name)==0) return 0;

    /* Parser */
    bindcommand(name,"runcommand",runcommand());

    bindcommand(name,"initphasefield",initphasefield());    

    /* Monte Carlo simulation */
    bindcommand(name,"eval",FreeEnergy());    
#ifdef _USECUDA
    bindcommand(name,"cuda_eval",cuda_FreeEnergy_single());    
#endif
    bindcommand(name,"run",run());    
    bindcommand(name,"srand",srand((unsigned int)_RANDSEED));
    bindcommand(name,"srand48",srand48((unsigned int)_RANDSEED));
    bindcommand(name,"srandbytime",srand((unsigned int)time(NULL)));
    bindcommand(name,"srand48bytime",srand48((unsigned int)time(NULL)));

    /* File input and output */
    bindcommand(name,"writecn",writefinalcnfile());
    bindcommand(name,"readcn",readcn());
    bindcommand(name,"openintercnfile",openintercnfile());
    bindcommand(name,"writeintercn",writeintercnfile());
    bindcommand(name,"setfilecounter",setfilecounter());
    bindcommand(name,"openpropfile",openpropfile());

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

void PhaseFieldFrame::runcommand()
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

void PhaseFieldFrame::initvars()
{
    int i, j;
    initparser();
    strcpy(command,"echo Hello World");
    
    strcpy(incnfile,"in.cn");
    strcpy(intercnfile,"inter.cn");
    strcpy(finalcnfile,"final.cn");
    strcpy(outpropfile,"prop.out");

    strcpy(bondcolor,"red");
    strcpy(highlightcolor,"purple");

    for(i=0;i<MAXNUMFIELDS;i++)
    {
        atomradius[i]=0.1;
        strcpy(atomcolor[i],"orange");
    }
    for(i=0;i<MAXCOLORS;i++)
    {
        strcpy(colornames[i],"gray50");
    }
    for(i=0;i<MAXNUMFIELDS;i++)
    {
        plot_threshold[2*i+0] = -0.05;
        plot_threshold[2*i+1] = 0.05;
    }

    _ROT_MATRIX[0][0] = 1.0; _ROT_MATRIX[0][1] = 0.0; _ROT_MATRIX[0][2] = 0.0; 
    _ROT_MATRIX[1][0] = 0.0; _ROT_MATRIX[1][1] = 1.0; _ROT_MATRIX[1][2] = 0.0; 
    _ROT_MATRIX[2][0] = 0.0; _ROT_MATRIX[2][1] = 0.0; _ROT_MATRIX[2][2] = 1.0; 
}

#ifdef _USETCL
/********************************************************************/
/* Initialize Tcl parser */
/********************************************************************/
int PhaseFieldFrame::Tcl_AppInit(Tcl_Interp *interp)
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
void PhaseFieldFrame::Alloc()
{
    char s[100];
    int n;

    Ntot = _NX*_NY*_NZ;

    if (num_fields <= 0)
      ERROR("Alloc: num_fields="<<num_fields);

    INFO("Alloc: num_fields="<<num_fields);

#ifndef _USEFLOAT
    for(n=0;n<num_fields;n++)
    {
      Realloc(_PHI[n],double,Ntot);
      memset(_PHI[n],0,sizeof(double)*Ntot);

      Realloc(_dPHIdt[n],double,Ntot);
      memset(_dPHIdt[n],0,sizeof(double)*Ntot);
      sprintf(s,"dPHIdt_%d",n); bindvar(s,_dPHIdt[n],DOUBLE);

      Realloc(_dPHIdx[n],double,Ntot);
      memset(_dPHIdx[n],0,sizeof(double)*Ntot);

      Realloc(_dPHIdy[n],double,Ntot);
      memset(_dPHIdy[n],0,sizeof(double)*Ntot);

      Realloc(_dPHIdz[n],double,Ntot);
      memset(_dPHIdz[n],0,sizeof(double)*Ntot);

        
      Realloc(_dFdPHIi[n],double,Ntot);
      memset(_dFdPHIi[n],0,sizeof(double)*Ntot);

      Realloc(_d2PHI[n],double,Ntot);
      memset(_d2PHI[n],0,sizeof(double)*Ntot);
       
    }


    
#else
    for(n=0;n<num_fields;n++)
    {
      Realloc(_PHI[n],float,Ntot);
      memset(_PHI[n],0,sizeof(float)*Ntot);

      Realloc(_dPHIdt[n],float,Ntot);
      memset(_dPHIdt[n],0,sizeof(float)*Ntot);

      Realloc(_dPHIdx[n],float,Ntot);
      memset(_dPHIdx[n],0,sizeof(float)*Ntot);

      Realloc(_dPHIdy[n],float,Ntot);
      memset(_dPHIdy[n],0,sizeof(float)*Ntot);

      Realloc(_dPHIdz[n],float,Ntot);
      memset(_dPHIdz[n],0,sizeof(float)*Ntot);

      Realloc(_d2PHI[n],float,Ntot);
      memset(_d2PHI[n],0,sizeof(float)*Ntot);

    }
#endif
}

/* initialization */
void PhaseFieldFrame::initphasefield()
{
    int i, j, k, n, size, ind;

    size = _NX*_NY*_NZ;
    
    Alloc();
    /* Initialize m matrix */
      
   
    if (input[0] == 1)
    {
        /* randomize */
        for(n=0;n<num_fields;n++)
          for(i=0;i<_NX;i++)
            for(j=0;j<_NY;j++)
                for(k=0;k<_NZ;k++)
                {
                    _PHI[n][i*_NY*_NZ+j*_NZ+k]=((drand48()-0.5)/500);
                }
    }
    else
    {
       ERROR("initphasefield: unknown type "<<input[0]);
    }
}


/********************************************************************/
/* Phase Field Simulation */
/********************************************************************/
void PhaseFieldFrame::run()
{
    int n, i, ind;
    double phase_volume[MAXNUMFIELDS];

    for(curstep=step0;curstep<(step0+totalsteps);curstep++)
    { 
        /* compute derivative of free energy functional */
        FreeEnergy();

        /* compute relative volume of the phases */
        for(n=0;n<num_fields;n++)
        {
             phase_volume[n] = 0;
             for (ind=0; ind<_NX*_NY*_NZ; ind++)
                    phase_volume[n] += _PHI[n][ind];
                    phase_volume[n] /= (_NX*_NY*_NZ);
        }

        
        /* adjust time step at every step */
        timestep = dtmax/fmax(_G*dtmax/dphimax,1);
       
        
        /* time integrator */
        for(n=0;n<num_fields;n++)
        {
             for (ind=0; ind<_NX*_NY*_NZ; ind++)
                 _PHI[n][ind] += _dPHIdt[n][ind] * timestep;
        }

        
        /* print information */
        if(curstep%printfreq==0)
        {
            if (num_fields == 1)
              INFO_Printf("curstep = %6d F = %18.12e G = %12.6e timestep = %8.3e rel_vol = %5.3g%%\n",
                          curstep,_F, _G, timestep, phase_volume[0]*100);
            else if (num_fields == 2)
              INFO_Printf("curstep = %6d F = %20.12e rel_vol = (%5.3g,%5.3g)%% \n",curstep,_F,
                           phase_volume[0]*100, phase_volume[1]*100);
            else if (num_fields == 3)
              INFO_Printf("curstep = %6d F = %20.12e G = %20.12e timestep = %5.3g \n"
                          " rel_vol = (%5.3g,%5.3g,%5.3g)%% M01=%5.3g M02=%5.3g \n",
                           curstep,_F, _G, timestep,
                           phase_volume[0]*100, phase_volume[1]*100, phase_volume[2]*100);
            else
              INFO_Printf("curstep = %6d F = %20.12e\n",curstep,_F);
     	}

	if (saveprop)
	{
	    if((curstep%savepropfreq)==0) pf.write(this);
	}

        if(savecn)
            if((curstep%savecnfreq)==0&&(curstep!=0))
                intercn.write(this);
        
        winplot();
    }
}

/* compute free energy of current phase field config */
double PhaseFieldFrame::FreeEnergy()
{
      return FreeEnergy_multi();
}

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////


double PhaseFieldFrame::FreeEnergy_multi() /* multi-phase field model (MPF) */
{
    int n, p, i, j, k, ip, im, jp, jm, kp, km, ind;
    double h2, f, g, dphidx, dphidy, dphidz, d2phi;
    double grad_term;

    _F = 0;

    if (gridsize<=0)
      FATAL("FreeEnergy: gridesize = "<<gridsize);

    
    /* Zero out time derivative of phase fields */
    for(n=0;n<num_fields;n++)
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {
            _dPHIdt[n][ind] = 0;
            _dFdPHIi[n][ind] = 0;
        }

    /* compute spatial gradients */
    h2 = gridsize*gridsize;
    for(n=0;n<num_fields;n++)
    {
        for(i=0;i<_NX;i++)
        {
            ip = (i+1)%_NX;
            im = (i-1+_NX)%_NX;
            for (j=0;j<_NY;j++)
            {
                jp = (j+1)%_NY;
                jm = (j-1+_NY)%_NY;
                for (k=0;k<_NZ;k++)
                {
                    kp = (k+1)%_NZ;
                    km = (k-1+_NZ)%_NZ;
                    
                    ind = i*_NY*_NZ+j*_NZ+k;
                    d2phi = ( _PHI[n][ip*_NY*_NZ+j*_NZ+k] + _PHI[n][im*_NY*_NZ+j*_NZ+k]
                             +_PHI[n][i*_NY*_NZ+jp*_NZ+k] + _PHI[n][i*_NY*_NZ+jm*_NZ+k]
                             +_PHI[n][i*_NY*_NZ+j*_NZ+kp] + _PHI[n][i*_NY*_NZ+j*_NZ+km]
                             -_PHI[n][ind] * 6.0 ) / h2;
                    
                    _d2PHI[n][ind] = d2phi;
                    
                    dphidx = (_PHI[n][ip*_NY*_NZ+j*_NZ+k] - _PHI[n][im*_NY*_NZ+j*_NZ+k]) / (2.0*gridsize);
                    dphidy = (_PHI[n][i*_NY*_NZ+jp*_NZ+k] - _PHI[n][i*_NY*_NZ+jm*_NZ+k]) / (2.0*gridsize);
                    dphidz = (_PHI[n][i*_NY*_NZ+j*_NZ+kp] - _PHI[n][i*_NY*_NZ+j*_NZ+km]) / (2.0*gridsize);
                    
                    _dPHIdx[n][ind] = dphidx;
                    _dPHIdy[n][ind] = dphidy;
                    _dPHIdz[n][ind] = dphidz;
                    
                }
            }
        }
    }

   
          
   for (ind=0; ind<_NX*_NY*_NZ; ind++)
   {
            
       for(n=0;n<num_fields;n++)
       {
           grad_term = SQR(_dPHIdx[n][ind]) + SQR(_dPHIdy[n][ind]) + SQR(_dPHIdz[n][ind]) ;
           _F += SQR(_PHI[n][ind])*(-alpha/2) + SQR(SQR(_PHI[n][ind]))*(beta/4) + grad_term*(kappa/2);                       
           for(p=0;p<num_fields;p++)
           {
               if ( p!=n )
               {
                   _F += gamma*SQR(_PHI[n][ind])*SQR(_PHI[p][ind]);
               }
           }
           
           _dFdPHIi[n][ind] = -alpha*_PHI[n][ind] + beta*CUBE(_PHI[n][ind]) - kappa*_d2PHI[n][ind]; 
           for(p=0;p<num_fields;p++)
           {
               if ( p!=n )
               { 
                   _dFdPHIi[n][ind] = _dFdPHIi[n][ind] + 2*gamma*_PHI[n][ind]*SQR(_PHI[p][ind]);
               }
           }
        }
	}
            
    for(n=0;n<num_fields;n++)
    {          
	  for (ind=0; ind<_NX*_NY*_NZ; ind++)
          {
           	_dPHIdt[n][ind] = -L*_dFdPHIi[n][ind];
           }
    }    
       
    /* calculate the max value of dphidt */
    double Gmiddle = 0;
    for (ind=0; ind<_NX*_NY*_NZ; ind++)   
    {                
        if ( Gmiddle <= fmax(fmax(fabs(_dPHIdt[0][ind]), fabs(_dPHIdt[1][ind])),fabs(_dPHIdt[2][ind])) )
            Gmiddle = fmax(fmax(fabs(_dPHIdt[0][ind]), fabs(_dPHIdt[1][ind])),fabs(_dPHIdt[2][ind]));
    }
    _G = Gmiddle;
        
       
# if 0
        INFO_Printf("_F = %20.12e\n", _F);
        INFO_Printf("_G = %20.12e\n", _G);
        write_array(_PHI[0],"PHI_0.dat");
        write_array(_PHI[1],"PHI_1.dat");
        write_array(_PHI[2],"PHI_2.dat");  
        sleep(sleepseconds);
# endif
    

    return _F*CUBE(gridsize);        

}



/********************************************************************/
/* File input and output */
/********************************************************************/
int PhaseFieldFrame::readcn()
{
    initcn.open(incnfile,LFile::O_Read);
    return initcn.read(this);
}

int PhaseFieldFrame::setfilecounter()
{
    intercn.setcount(filecounter);
    return 0;
}

int PhaseFieldFrame::writefinalcnfile(int zip, bool bg)
{
    finalcn.open(finalcnfile);
    return finalcn.write(this,zip,bg);
}

void PhaseFieldFrame::saveintercn(int step)
{
    if(savecn)
        if((step%savecnfreq)==0&&(step!=0))
            writeintercnfile();
}

int PhaseFieldFrame::writeintercnfile(int zip,bool bg)
{
    return intercn.write(this,zip,bg);
}

int PhaseFieldFrame::openintercnfile()
{
    return intercn.open(intercnfile);
}

int PhaseFieldFrame::openpropfile()
{
    return pf.open(outpropfile);
}

/* Configuration Files */
char * CNFile::describe()
{
    static char tmp[500];
    sprintf(tmp,"%s","phase field configurations, Format:\n"
            "num_fields _NX _NY _NZ; Phi[0](0,0,0); s(0,0,1); ...; s(_NX-1,_NY-1,_NZ-1)");
    return tmp;
}

int CNFile::writeblock(void *p)
{
    int n,i;
    PhaseFieldFrame &d=*((PhaseFieldFrame *)p);
    f->Printf("%10d\n%10d\n%10d\n%10d\n",d.num_fields,d._NX,d._NY,d._NZ);
    for(n=0;n<d.num_fields;n++)
      for(i=0;i<d.Ntot;i++)
      {
        f->Printf("%20.12e\n",d._PHI[n][i]);
      }
    return 0;
}

int CNFile::readblock(void *p)
{
    int n,i;
    char *buffer, *pp, *q;
    PhaseFieldFrame &d=*((PhaseFieldFrame *)p);
    
    LFile::LoadToString(fname,buffer,0);

    pp=buffer;
    sscanf(pp, "%d %d %d %d", &d.num_fields, &d._NX, &d._NY, &d._NZ);
    INFO("readblock: num_fields="<<d.num_fields<<"NX="<<d._NX<<" NY="<<d._NY<<" NZ="<<d._NZ);    
    
    d.Alloc();
    
    pp=strchr(pp, '\n');
    pp++;
    pp=strchr(pp, '\n');
    pp++;
    pp=strchr(pp, '\n');
    pp++;
    pp=strchr(pp, '\n');
    pp++;

    for(n=0;n<d.num_fields;n++)
      for(i=0;i<d.Ntot;i++)
      {
        q=pp;
        pp=strchr(pp,'\n');
        if(pp) *(char *)pp++=0;
#ifndef _USEFLOAT
        sscanf(q, "%le",&(d._PHI[n][i]));
#else
        sscanf(q, "%e",&(d._PHI[n][i]));
#endif
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
            "curstep free_energy ");
    return tmp;
}

int PropFile::writeentry(void *p)
{
    char tmp[500];
    PhaseFieldFrame &d=*((PhaseFieldFrame *)p);
    sprintf(tmp,"%10d %20.12e %20.12e\n",
            d.curstep,d._F,d._G);
    *f<<tmp;
    return 0;
}

void PhaseFieldFrame::write_array(double * data, const char *fname)
{
    FILE *fp;
    int i;
    fp = fopen(fname,"w");
    for(i=0;i<_NX*_NY*_NZ;i++)
       fprintf(fp,"%20.12e\n",data[i]);
    fclose(fp);
}

/********************************************************************/
/* Visualization */
/********************************************************************/

void PhaseFieldFrame::openwin()
{
    openwindow(win_width,win_height,dirname);
}

int PhaseFieldFrame::openwindow(int w,int h,const char *n)
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

void PhaseFieldFrame::closewindow() {delete(win);}

void PhaseFieldFrame::winplot()
{
    if(win!=NULL)
        if(win->alive)
            if((curstep%plotfreq)==0)
            {
                plot();
            }
}

void PhaseFieldFrame::winplot(int step)
{
    if(win!=NULL)
        if(win->alive)
            if((step%plotfreq)==0)
            {
                plot();
            }
}        

void PhaseFieldFrame::wintogglepause()
{
    if(win!=NULL)
        if(win->IsAlive())
            win->TogglePause();
}

#include "namecolor.c"
void PhaseFieldFrame::alloccolors()
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

            for(int i=0;i<MAXNUMFIELDS;i++)
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
void PhaseFieldFrame::alloccolorsX()
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

void PhaseFieldFrame::rotate()
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

void PhaseFieldFrame::saverot()
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
  
void PhaseFieldFrame::plot()
{
    int n, i, j, k, NMAX;
    double val;
    unsigned long c;
    Vector3 s1, s2, vr1, vr2, sri, srj, ri, rj;
    
    if(win==NULL) return;
    if(!(win->alive)) return;
    
    win->Lock();
    win->Clear();

    NMAX = max(max(_NX, _NY), _NZ);
    
    /* draw atom s*/
    for(n=0;n<num_fields;n++)
      for(i=0;i<_NX;i++)
        for(j=0;j<_NY;j++)
            for(k=0;k<_NZ;k++)
            {
                val = _PHI[n][i*_NY*_NZ+j*_NZ+k];

                //INFO_Printf("plot: phi[%d] %d,%d,%d = %20.12e \n", n, i, j, k, val);

                /* make these threshold value changable variables later */
                if (val<plot_threshold[2*n+0] || val>plot_threshold[2*n+1]) continue;

                //INFO_Printf("plot: draw point %d,%d,%d with color %d\n", i, j, k, n);
               
                c = colors[n]; 

                win->DrawPoint((2.0*i-_NX)/NMAX,-(2.0*j-_NY)/NMAX,(2.0*k-_NZ)/NMAX,
                               atomradius[0]/NMAX,c);
            }
    
    /* draw frame */
#define drawsline(a,b,c,d,e,f,g,h,i) \
        win->DrawLine(a,b,c,d,e,f,g,h,i);
    
    drawsline(-1,-1,-1,-1,-1, 1,colors[MAXCOLORS+2],0,0);
    drawsline(-1,-1, 1,-1, 1, 1,colors[MAXCOLORS+2],0,0);
    drawsline(-1, 1, 1,-1, 1,-1,colors[MAXCOLORS+2],0,0);
    drawsline(-1, 1,-1,-1,-1,-1,colors[MAXCOLORS+2],0,0);
    drawsline( 1,-1,-1, 1,-1, 1,colors[MAXCOLORS+2],0,0);
    drawsline( 1,-1, 1, 1, 1, 1,colors[MAXCOLORS+2],0,0);
    drawsline( 1, 1, 1, 1, 1,-1,colors[MAXCOLORS+2],0,0);
    drawsline( 1, 1,-1, 1,-1,-1,colors[MAXCOLORS+2],0,0);
    drawsline(-1,-1,-1, 1,-1,-1,colors[MAXCOLORS+2],0,0);
    drawsline(-1,-1, 1, 1,-1, 1,colors[MAXCOLORS+2],0,0);
    drawsline(-1, 1, 1, 1, 1, 1,colors[MAXCOLORS+2],0,0);
    drawsline(-1, 1,-1, 1, 1,-1,colors[MAXCOLORS+2],0,0);
    
    win->Unlock();
    win->Refresh();

}




/***********************************************************/
/*  Main Program Begins */

#ifdef _TEST

class PhaseFieldFrame sim;

/* The main program is defined here */
#include "main.cpp"

#endif//_TEST
