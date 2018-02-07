// Following Folch 2003 formulation 
// Last Modified : Mon July 8, 2014

#include <iostream>
#include "phasefield.h"

#ifdef _STK_MPI
#include "phasefield_mpi.cpp"
#include <mpi.h>
#else
#ifdef _OPENMP
#include "phasefield_omp.cpp"
#endif
#endif


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
    bindvar("model_type",&model_type,INT);
    bindvar("ROT_MATRIX",&(_ROT_MATRIX[0][0]),DOUBLE);

    /* Simulation setting */
    bindvar("totalsteps",&totalsteps,INT);
    bindvar("savepropfreq",&savepropfreq,INT);
    bindvar("randseed",&_RANDSEED,INT);
    bindvar("gridsize",&gridsize,DOUBLE);
    bindvar("timestep",&timestep,DOUBLE);
    bindvar("dtmax",&dtmax,DOUBLE);
    bindvar("dphimax",&dphimax,DOUBLE);
    bindvar("dynamics_type",&dynamics_type,INT);
    bindvar("Mob_GL",&Mob_GL,DOUBLE);
    bindvar("Mob_D",&Mob_D,DOUBLE);
    bindvar("Mob_SV",&Mob_SV,DOUBLE);
    bindvar("Mob_LS",&Mob_LS,DOUBLE);
    bindvar("Mob_LV",&Mob_LV,DOUBLE);
    bindvar("Fext_x",&Fext_x,DOUBLE);
    bindvar("Fext_y",&Fext_y,DOUBLE);
    bindvar("Penalty_3phase",&Penalty_3phase,DOUBLE);
    bindvar("stop_flag_x",&stop_flag_x,INT);
    bindvar("stop_flag_y",&stop_flag_y,INT);
    bindvar("stop_flag_solid",&stop_flag_solid,INT);
    bindvar("solid_obj",&solid_obj,DOUBLE);
    bindvar("com_x_obj",&com_x_obj,DOUBLE);
    bindvar("com_y_obj",&com_y_obj,DOUBLE);
    bindvar("save_NW",&save_NW,INT);
    bindvar("Penalty_NW",&Penalty_NW,DOUBLE);
    bindvar("vol_incr0",&vol_incr0,DOUBLE);
    bindvar("x_incr0",&x_incr0,DOUBLE);
    bindvar("y_incr0",&y_incr0,DOUBLE);
    bindvar("create_vapor",&create_vapor,INT);
    bindvar("createfreq",&createfreq,INT);
    bindvar("shiftheight",&shiftheight,INT);

    /* Simulation variables */
    bindvar("F",&_F,DOUBLE);
    bindvar("Fraw",&_Fraw,DOUBLE);
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
#ifdef _STK_MPI
    bindvar("myDomain",&myDomain,INT);
    bindvar("mySize",&mySize,INT);
#endif
    
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
       for(j=0;j<MAXNUMFIELDS;j++)
       {
           sprintf(s,"EPS_%d%d",i,j);
           bindvar(s,&(_EPS[i][j]),DOUBLE);

           sprintf(s,"EPS1_%d%d",i,j);
           bindvar(s,&(_EPS1[i][j]),DOUBLE);
           sprintf(s,"EPS2_%d%d",i,j);
           bindvar(s,&(_EPS2[i][j]),DOUBLE);
           sprintf(s,"EPS3_%d%d",i,j);
           bindvar(s,&(_EPS3[i][j]),DOUBLE);

           sprintf(s,"U_%d%d",i,j);
           bindvar(s,&(_U[i][j]),DOUBLE);
           sprintf(s,"M_%d%d",i,j);
           bindvar(s,&(_M[i][j]),DOUBLE);
       }
    for(i=0;i<MAXNUMFIELDS;i++)
    {
        sprintf(s,"atomcolor%d",i);
        bindvar(s,atomcolor[i],STRING);
        sprintf(s,"MU_%d",i);
        bindvar(s,&(_MU[i]),DOUBLE);
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
        phase_volume[i] = 0;
    }
    for(i=0;i<MAXCOLORS;i++)
    {
        strcpy(colornames[i],"gray50");
    }
    for(i=0;i<MAXNUMFIELDS;i++)
    {
        for(j=0;j<MAXNUMFIELDS;j++)
        {
           if (i!=j) 
           _U[i][j] = _EPS[i][j] = 0.0;
           else
           _U[i][j] = _EPS[i][j] = 1.0;

           _EPS1[i][j] = _EPS2[i][j] = _EPS3[i][j] = 0.0;
        }
        plot_threshold[2*i+0] = 0.45;
        plot_threshold[2*i+1] = 0.55;
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

/* For MPI implementation based on the StencilToolkit */
#ifdef _STK_MPI
    /* should specify the number of cores here */
    _node = new Node3D(2, 2, 2);
    _idx = new Index3D(*_node, _NZ, _NY, _NX);    
    MPI_Comm_size(MPI_COMM_WORLD, &mySize);
    MPI_Comm_rank(MPI_COMM_WORLD, &myDomain);
    
    data_size = Ntot*num_fields/mySize;
    phase_size = Ntot/mySize;
    total_size = Ntot*num_fields;
    
    Realloc(sendbuf, REAL_T, data_size);
    memset(sendbuf, 0, sizeof(REAL_T) * data_size);
    Realloc(sendbuf_ind, int, data_size*3);
    memset(sendbuf_ind, 0, sizeof(int) * data_size*3);
    
    if (myDomain == 0)
    {
        Realloc(rbuf, REAL_T, total_size);
        memset(rbuf, 0, sizeof(REAL_T) * total_size);
        Realloc(rbuf_ind, int, total_size*3);
        memset(rbuf_ind, 0, sizeof(int) * total_size*3);
    }
#endif

    for(int n = 0; n < num_fields; n++)
    {
        Realloc(_PHI[n], REAL_T, Ntot);
        memset(_PHI[n], 0, sizeof(REAL_T) * Ntot);
        
#ifdef _STK_MPI                
        s_PHI[n] = new Grid3D<REAL_T>(*_idx, 3, 3, 3);        
        s_dPHIdt[n] = new Grid3D<REAL_T>(*_idx);
        s_dPHIdt0[n] = new Grid3D<REAL_T>(*_idx);
        s_dPHIdtCH[n] = new Grid3D<REAL_T>(*_idx, 1, 1, 1);
        s_d2dPHIdt0[n] = new Grid3D<REAL_T>(*_idx);
        s_dPHIdx[n] = new Grid3D<REAL_T>(*_idx);
        s_dPHIdy[n] = new Grid3D<REAL_T>(*_idx);
        s_dPHIdz[n] = new Grid3D<REAL_T>(*_idx);
        s_dFdPHI[n] = new Grid3D<REAL_T>(*_idx, 1, 1, 1);
        s_d2PHI[n] = new Grid3D<REAL_T>(*_idx);
        s_EPS_loc[n] = new Grid3D<REAL_T>(*_idx, 1, 1, 1);
        s_tmp_x[n] = new Grid3D<REAL_T>(*_idx, 1, 1, 1);
        s_tmp_y[n] = new Grid3D<REAL_T>(*_idx, 1, 1, 1);
        s_tmp_z[n] = new Grid3D<REAL_T>(*_idx, 1, 1, 1);
        s_d2dFdPHI[n] = new Grid3D<REAL_T>(*_idx);
#else
        Realloc(_dPHIdt[n], REAL_T, Ntot);
        memset(_dPHIdt[n], 0, sizeof(REAL_T) * Ntot);
        sprintf(s, "dPHIdt_%d", n); bindvar(s, _dPHIdt[n], REAL_S);
        
        Realloc(_dPHIdt0[n], REAL_T, Ntot);
        memset(_dPHIdt0[n], 0, sizeof(REAL_T) * Ntot);
        sprintf(s, "dPHIdt0_%d", n); bindvar(s, _dPHIdt0[n], REAL_S);
        
        Realloc(_dPHIdtCH[n], REAL_T, Ntot);
        memset(_dPHIdtCH[n], 0, sizeof(REAL_T) * Ntot);
        
        Realloc(_d2dPHIdt0[n], REAL_T, Ntot);
        memset(_d2dPHIdt0[n], 0, sizeof(REAL_T) * Ntot);
        
        Realloc(_dPHIdx[n], REAL_T, Ntot);
        memset(_dPHIdx[n], 0, sizeof(REAL_T) * Ntot);
        
        Realloc(_dPHIdy[n], REAL_T, Ntot);
        memset(_dPHIdy[n], 0, sizeof(REAL_T) * Ntot);
        
        Realloc(_dPHIdz[n], REAL_T, Ntot);
        memset(_dPHIdz[n], 0, sizeof(REAL_T) * Ntot);
        
        Realloc(_dFdPHI[n], REAL_T, Ntot);
        memset(_dFdPHI[n], 0, sizeof(REAL_T) * Ntot);
        sprintf(s, "dFdPHI_%d", n); bindvar(s, _dFdPHI[n], REAL_S);
        
        Realloc(_dFdPHIi[n], REAL_T, Ntot);
        memset(_dFdPHIi[n], 0, sizeof(REAL_T) * Ntot);
        
        Realloc(_d2PHI[n], REAL_T, Ntot);
        memset(_d2PHI[n], 0, sizeof(REAL_T) * Ntot);
        
        Realloc(_EPS_loc[n], REAL_T, Ntot);
        memset(_EPS_loc[n], 0, sizeof(REAL_T) * Ntot);
        
        Realloc(tmp_x[n], REAL_T, Ntot);
        memset(tmp_x[n], 0, sizeof(REAL_T) * Ntot);
        
        Realloc(tmp_y[n], REAL_T, Ntot);
        memset(tmp_y[n], 0, sizeof(REAL_T) * Ntot);
        
        Realloc(tmp_z[n], REAL_T, Ntot);
        memset(tmp_z[n], 0, sizeof(REAL_T) * Ntot);
        
        Realloc(_d2dFdPHI[n], REAL_T, Ntot);
        memset(_d2dFdPHI[n], 0, sizeof(REAL_T) * Ntot);
#endif
    }
    
#ifdef _STK_MPI
    s_mprime = new Grid3D<REAL_T>(*_idx);
    s_K_S = new Grid3D<REAL_T>(*_idx);
    s_K_L = new Grid3D<REAL_T>(*_idx);
    s_K_SV = new Grid3D<REAL_T>(*_idx);
    s_K_LS = new Grid3D<REAL_T>(*_idx);
    s_K_LV = new Grid3D<REAL_T>(*_idx);
    s_K_other = new Grid3D<REAL_T>(*_idx);
    s_K_D = new Grid3D<REAL_T>(*_idx);
    s_matrix_x = new Grid3D<REAL_T>(*_idx);
    s_matrix_y = new Grid3D<REAL_T>(*_idx);
    s_NW_orig = new Grid3D<REAL_T>(*_idx);
    s_dGLdmuL = new Grid3D<REAL_T>(*_idx);
    s_dGLdmuV = new Grid3D<REAL_T>(*_idx);
    s_dGVdmuL = new Grid3D<REAL_T>(*_idx);
    s_dGVdmuV = new Grid3D<REAL_T>(*_idx);
    s_dGSdmuL = new Grid3D<REAL_T>(*_idx);
    s_dGSdmuV = new Grid3D<REAL_T>(*_idx);
        
#else
    Realloc(_mprime, REAL_T, Ntot);
    memset(_mprime, 0, sizeof(REAL_T) * Ntot);
    
    Realloc(_K_S, REAL_T, Ntot);
    memset(_K_S, 0, sizeof(REAL_T) * Ntot);
    
    Realloc(_K_L, REAL_T, Ntot);
    memset(_K_L, 0, sizeof(REAL_T) * Ntot);
    
    Realloc(_K_SV, REAL_T, Ntot);
    memset(_K_SV, 0, sizeof(REAL_T) * Ntot);
    
    Realloc(_K_LS, REAL_T, Ntot);
    memset(_K_LS, 0, sizeof(REAL_T) * Ntot);
    
    Realloc(_K_LV, REAL_T, Ntot);
    memset(_K_LV, 0, sizeof(REAL_T) * Ntot);
    
    Realloc(_K_other, REAL_T, Ntot);
    memset(_K_other, 0, sizeof(REAL_T) * Ntot);
    
    Realloc(_K_D, REAL_T, Ntot);
    memset(_K_D, 0, sizeof(REAL_T) * Ntot);
    
    Realloc(matrix_x, REAL_T, Ntot);
    memset(matrix_x, 0, sizeof(REAL_T) * Ntot);
    
    Realloc(matrix_y, REAL_T, Ntot);
    memset(matrix_y, 0, sizeof(REAL_T) * Ntot);
    
    Realloc(_dGLdmuL,REAL_T,Ntot);
    memset(_dGLdmuL,0,sizeof(REAL_T)*Ntot);
    
    Realloc(_dGLdmuV,REAL_T,Ntot);
    memset(_dGLdmuV,0,sizeof(REAL_T)*Ntot);
    
    Realloc(_dGVdmuL,REAL_T,Ntot);
    memset(_dGVdmuL,0,sizeof(REAL_T)*Ntot);
    
    Realloc(_dGVdmuV,REAL_T,Ntot);
    memset(_dGVdmuV,0,sizeof(REAL_T)*Ntot);
    
    Realloc(_dGSdmuL,REAL_T,Ntot);
    memset(_dGSdmuL,0,sizeof(REAL_T)*Ntot);
    
    Realloc(_dGSdmuV,REAL_T,Ntot);
    memset(_dGSdmuV,0,sizeof(REAL_T)*Ntot);
    
    Realloc(_NW_orig,REAL_T,Ntot);
    memset(_NW_orig,0,sizeof(REAL_T)*Ntot);

#endif
    
#ifdef _USECUDA
    cuda_memory_alloc();
#endif
    
}    
    
/* initialization */
void PhaseFieldFrame::initphasefield()
{
    int i, j, k, n, size, ind;

    size = _NX*_NY*_NZ;
    
    Alloc();
    /* Initialize m matrix */
   
    for (i=0;i<3;i++)
    {
    	for (j=0;j<3;j++)
        {
            _M[i][j] = _MU[i]-_MU[j];
        }
    }

    if (input[0] == 0)
    {
        for(n=0;n<num_fields;n++)
        {
           memset(_PHI[n],0,sizeof(int)*size);
        }
    }
    else if (input[0] == 1)
    {
        for(n=0;n<num_fields;n++)
        {
           memset(_PHI[n],1,sizeof(int)*size);
        }
    }
    else if (input[0] == -1)
    {
        /* randomize */
        for(n=0;n<num_fields;n++)
          for(i=0;i<_NX;i++)
            for(j=0;j<_NY;j++)
                for(k=0;k<_NZ;k++)
                {
                    _PHI[n][i*_NY*_NZ+j*_NZ+k]=(drand48());
                }
    }
    else if (input[0] == 2)
    {
        /* cube */
        for(n=0;n<num_fields;n++)
          for(i=0;i<_NX;i++)
            for(j=0;j<_NY;j++)
                for(k=0;k<_NZ;k++)
                {
                    if ((i+1>_NX/4)&&(i+1<_NX*3/4)&&(j+1>_NY/4)&&(j+1<_NY*3/4)&&(k+1>_NZ/4)&&(k+1<_NZ*3/4))
                       _PHI[n][i*_NY*_NZ+j*_NZ+k]=1.0-1e-6;
                    else
                       _PHI[n][i*_NY*_NZ+j*_NZ+k]=1e-6;
                }
    }
    else if (input[0] == 3)
    {
        /* three beams */
        for(n=0;n<num_fields;n++)
          for(i=0;i<_NX;i++)
            for(j=0;j<_NY;j++)
                for(k=0;k<_NZ;k++)
                {
                    if ( ((i+1>_NX/4)&&(i+1<_NX*3/4)&&(j+1>_NY/4)&&(j+1<_NY*3/4)) ||
                         ((j+1>_NY/4)&&(j+1<_NY*3/4)&&(k+1>_NZ/4)&&(k+1<_NZ*3/4)) ||
                         ((i+1>_NX/4)&&(i+1<_NX*3/4)&&(k+1>_NZ/4)&&(k+1<_NZ*3/4)) )
                       _PHI[n][i*_NY*_NZ+j*_NZ+k]=1.0-1e-6;
                    else
                       _PHI[n][i*_NY*_NZ+j*_NZ+k]=1e-6;
                }
    }
    else if (input[0] == 4)
    {
        /* sphere */
        for(n=0;n<num_fields;n++)
          for(i=0;i<_NX;i++)
            for(j=0;j<_NY;j++)
                for(k=0;k<_NZ;k++)
                {
                    if (SQR(i+1-_NX/2)+SQR(j+1-_NY/2)+SQR(k+1-_NZ/2)<SQR(max(max(_NX,_NY),_NZ)*input[1]))
                       _PHI[n][i*_NY*_NZ+j*_NZ+k]=1.0-1e-6;
                    else
                       _PHI[n][i*_NY*_NZ+j*_NZ+k]=1e-6;
                }
    }
    else if (input[0] == 5)
    {
        /* spherical liquid at planar solid-vapor interface */
        if (num_fields != 3)
        {
           ERROR("initphasefield: input = 5 only works for n = 3!");
           return;
        }
        for(i=0;i<_NX;i++)
          for(j=0;j<_NY;j++)
              for(k=0;k<_NZ;k++)
              {
                 if (SQR(i+1-_NX/2)+SQR(j+1-_NY/2)+SQR(k+1-_NZ*input[2])<SQR(max(_NX,_NY)*input[1]))
                    _PHI[0][i*_NY*_NZ+j*_NZ+k] = 1.0-1e-6;
                 else
                    _PHI[0][i*_NY*_NZ+j*_NZ+k] = 0.0+1e-6;

                  if (k+1 < _NZ*input[2])
                    _PHI[1][i*_NY*_NZ+j*_NZ+k] = 1.0-1e-6 - _PHI[0][i*_NY*_NZ+j*_NZ+k];
                  else
                    _PHI[1][i*_NY*_NZ+j*_NZ+k] = 0.0+1e-6;

                 _PHI[2][i*_NY*_NZ+j*_NZ+k] = 1.0 - _PHI[0][i*_NY*_NZ+j*_NZ+k] 
                                                  - _PHI[1][i*_NY*_NZ+j*_NZ+k];
              }
    }
    else if (input[0] == 6)
    {
        /* cylindrical liquid at planar solid-vapor interface */
        if (num_fields != 3)
        {
           ERROR("initphasefield: input = 6 only works for n = 3!");
           return;
        }
        for(i=0;i<_NX;i++)
          for(j=0;j<_NY;j++)
              for(k=0;k<_NZ;k++)
              {
                 if (SQR(j+1-_NY/2)+SQR(k+1-_NZ/2)<SQR(max(max(_NX,_NY),_NZ)*input[1]))
                    _PHI[0][i*_NY*_NZ+j*_NZ+k] = 1.0-1e-6;
                 else
                    _PHI[0][i*_NY*_NZ+j*_NZ+k] = 0.0+1e-6;

                 if (k+1 < _NZ/2)
                    _PHI[1][i*_NY*_NZ+j*_NZ+k] = 1.0-1e-6 - _PHI[0][i*_NY*_NZ+j*_NZ+k];
                 else
                    _PHI[1][i*_NY*_NZ+j*_NZ+k] = 0.0+1e-6;


                 _PHI[2][i*_NY*_NZ+j*_NZ+k] = 1.0 - _PHI[0][i*_NY*_NZ+j*_NZ+k] 
                                                  - _PHI[1][i*_NY*_NZ+j*_NZ+k];
              }
    }
    else if (input[0] == 7)
    {
        /* spherical liquid cap on top of planar solid-vapor interface 
         *  input = [ relative_radius relative_sphere_center cos_contact_angle ]
         *  suggested value: input = [ 7 0.35  0.125  0.74 ]
         */
        if (num_fields != 3)
        {
           ERROR("initphasefield: input = 7 only works for n = 3!");
           return;
        }
        double radius, Zc, cos_theta;

        radius    = max(_NX,_NY)*input[1];
        Zc        = _NZ*input[2];
        cos_theta = input[3];

        for(i=0;i<_NX;i++)
          for(j=0;j<_NY;j++)
              for(k=0;k<_NZ;k++)
              {
                 if (SQR(i+1-_NX/2)+SQR(j+1-_NY/2)+SQR(k+1-Zc)<SQR(max(_NX,_NY)*input[1]) && 
                     k+1-Zc-radius*cos_theta>=0)
                    _PHI[0][i*_NY*_NZ+j*_NZ+k] = 1.0-1.1e-6;
                 else
                    _PHI[0][i*_NY*_NZ+j*_NZ+k] = 0.0+1.0e-6;

                 if (k+1 < Zc+radius*cos_theta)
                    _PHI[1][i*_NY*_NZ+j*_NZ+k] = 1.0-1.15e-6 - _PHI[0][i*_NY*_NZ+j*_NZ+k];
                 else
                    _PHI[1][i*_NY*_NZ+j*_NZ+k] = 0.0+2e-6;

                 _PHI[2][i*_NY*_NZ+j*_NZ+k] = 1.0 - _PHI[0][i*_NY*_NZ+j*_NZ+k]
                                                  - _PHI[1][i*_NY*_NZ+j*_NZ+k];
              }
    }
    else if (input[0] == 8)
    {
        /* cylinder along z direction */
          for(i=0;i<_NX;i++)
            for(j=0;j<_NY;j++)
                for(k=0;k<_NZ;k++)
                {
                    if (SQR(i+1-_NX/2)+SQR(j+1-_NY/2)<SQR(max(_NX,_NY)*input[1]))
                       _PHI[1][i*_NY*_NZ+j*_NZ+k]=1.0-1e-6+0.01*(drand48()-0.5);
                    else
                       _PHI[2][i*_NY*_NZ+j*_NZ+k]=1.0-1e-6;
                   _PHI[0][i*_NY*_NZ+j*_NZ+k] = 1.0 - _PHI[1][i*_NY*_NZ+j*_NZ+k]
                                                  - _PHI[2][i*_NY*_NZ+j*_NZ+k];
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
    int n, i, j, k, ind;
    
#ifdef _STK_MPI
        
    const int z_begin = _idx->get_range_x().begin();
    const int z_end   = _idx->get_range_x().end();

    const int y_begin = _idx->get_range_y().begin();
    const int y_end   = _idx->get_range_y().end();

    const int x_begin = _idx->get_range_z().begin();
    const int x_end   = _idx->get_range_z().end();
#else
    const int z_begin = 0;
    const int z_end   = _NZ;

    const int y_begin = 0;
    const int y_end   = _NY;

    const int x_begin = 0;
    const int x_end   = _NX;
#endif


#ifdef _STK_MPI
    /* assign phase field values to each node */
    for (n = 0; n < num_fields; n++)
        for (i = x_begin; i < x_end; i++)
            for (j = y_begin; j < y_end; j++)
                for (k = z_begin; k < z_end; k++)
                    (*s_PHI[n])(k, j, i) = _PHI[n][i*_NY*_NZ+j*_NZ+k];
#endif


    if (save_NW == 1) 
    {
        for (i = x_begin; i < x_end; i++)
            for (j = y_begin; j < y_end; j++)
                for (k = z_begin; k < z_end; k++)
                {
#ifdef _STK_MPI
                    (*s_NW_orig)(k, j, i) = _PHI[1][i*_NY*_NZ+j*_NZ+k];
#else
                    ind = i*_NY*_NZ+j*_NZ+k;
                    _NW_orig[ind] = _PHI[1][i*_NY*_NZ+j*_NZ+k];
#endif    
                }
    }

    for (int i = x_begin; i < x_end; i++)
	{
		for (int j = y_begin; j < y_end; j++)
		{
			for (int k = z_begin; k < z_end; k++)
			{
#ifdef _STK_MPI
				/* calculate the position matrix for external force term */
				(*s_matrix_x)(k, j, i) = (i - _NX / 2)*gridsize;
				(*s_matrix_y)(k, j, i) = (j - _NY / 2)*gridsize;
#else                                
                        ind = i*_NY*_NZ+j*_NZ+k;
                		matrix_x[ind] = (i - _NX/2)*gridsize;
                		matrix_y[ind] = (j - _NY/2)*gridsize;
#endif                
			}
		}
	}

#ifdef _USECUDA
    cuda_init();
#endif
    
    for(curstep=step0;curstep<(step0+totalsteps);curstep++)
    {   
        
        liquid_volume = 0; com_x = 0; com_y = 0;
        
if ( model_type != 20 )
{
        for (i = x_begin; i < x_end; i++)
            for (j = y_begin; j < y_end; j++)
                for (k = z_begin; k < z_end; k++)
                {
#ifdef _STK_MPI
                    double phi_0 = (*s_PHI[0])(k, j, i);
                    liquid_volume += tanh((phi_0-0.5)/0.10)/2 + 0.5; 
#else
                    int idx = i*_NY*_NZ + j*_NZ + k;
                    liquid_volume += tanh((_PHI[0][idx]-0.5)/0.10)/2 + 0.5;      
#endif
                }  
        
#ifdef _STK_MPI
        _node->reduce(STK_SUM, STK_DOUBLE, &liquid_volume);
#endif         
        
}
        /* compute derivative of free energy functional */
        /* calculate dphidt and time step */
        FreeEnergy();
        
        
if ( model_type != 20 )
{
         /* compute relative volume of the phases */
        for(n = 0; n < num_fields; n++)
        {
            phase_volume[n] = 0;
            for (i = x_begin; i < x_end; i++)
                for (j = y_begin; j < y_end; j++)
                    for (k = z_begin; k < z_end; k++)
                    {
#ifdef _STK_MPI
                        phase_volume[n] += (*s_PHI[n])(k, j, i);
#else
                        phase_volume[n] += _PHI[n][i*_NY*_NZ+j*_NZ+k];
#endif
                    }

#ifdef _STK_MPI
            _node->reduce(STK_SUM, STK_DOUBLE, &(phase_volume[n]));
#endif
            phase_volume[n] /= (_NX * _NY * _NZ);
        }

    for (i = x_begin; i < x_end; i++)
        for (j = y_begin; j < y_end; j++)
            for (k = z_begin; k < z_end; k++)
            {
#ifdef _STK_MPI
                double phi_0 = (*s_PHI[0])(k, j, i);
                com_x += (tanh((phi_0-0.5)/0.10)/2 + 0.5)*(*s_matrix_x)(k,j,i)/liquid_volume;
                com_y += (tanh((phi_0-0.5)/0.10)/2 + 0.5)*(*s_matrix_y)(k,j,i)/liquid_volume;
#else
                int idx = i*_NY*_NZ + j*_NZ + k;
                com_x += (tanh((_PHI[0][idx]-0.5)/0.10)/2 + 0.5)*matrix_x[idx]/liquid_volume;
                com_y += (tanh((_PHI[0][idx]-0.5)/0.10)/2 + 0.5)*matrix_y[idx]/liquid_volume;
#endif
            }
    
#ifdef _STK_MPI
    _node->reduce(STK_SUM, STK_DOUBLE, &com_x);
    _node->reduce(STK_SUM, STK_DOUBLE, &com_y);
#endif
    
     /* time integrator */
        for(n = 0; n < num_fields; n++)
        {
            for (i = x_begin; i < x_end; i++)
                for (j = y_begin; j < y_end; j++)
                    for (k = z_begin; k < z_end; k++)
                    {
#ifdef _STK_MPI
                        (*s_PHI[n])(k, j, i) += (*s_dPHIdt[n])(k, j, i) * timestep;
                        _PHI[n][(i * _NY * _NZ) + (j * _NZ) + k] = (*s_PHI[n])(k, j, i);
#else
                        int idx = i*_NY*_NZ + j*_NZ + k;
                        _PHI[n][idx] += _dPHIdt[n][idx]*timestep;
#endif
                    }
        }

}
else if (model_type == 20)
    {
    #ifdef _USECUDA
        cuda_integrator();
    #endif
    }

        
        /* if the x center of mass position reaches the designated value, break */
     if (stop_flag_x == 1) 
        {   
            if (fabs(com_x - com_x_obj)<= 0.25) break;
        }

        /* if the y center of mass position reaches the designated value, break */
      if (stop_flag_y == 1) 
        {   
            if (fabs(com_y - com_y_obj)<= 0.25) break;
        }
         
      if (stop_flag_solid == 1)       
        {
            if (fabs(phase_volume[1] - solid_obj) <= 0.0005) break;
        } 

        /* shift the simulation box and create vapor phase */
        if (create_vapor == 1)
        {
           if (curstep%createfreq==0)
           {
               
#ifdef _STK_MPI
               /* pack all data into an ordered 1D array */
               MPI_Barrier(MPI_COMM_WORLD);
               
               int i_pack = 0;
               
               for (n = 0; n < num_fields; n++)
                   for (i = x_begin; i < x_end; i++)
                       for (j = y_begin; j < y_end; j++)
                           for (k = z_begin; k < z_end; k++)
                           {
                               sendbuf[i_pack] = (*s_PHI[n])(k, j, i);
                               
                               sendbuf_ind[i_pack*3] = i;
                               sendbuf_ind[i_pack*3+1] = j;
                               sendbuf_ind[i_pack*3+2] = k;
                               
                               i_pack++;
                           }
               
               MPI_Barrier(MPI_COMM_WORLD);
               MPI_Gather(sendbuf, data_size, MPI_DOUBLE, rbuf, data_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
               MPI_Gather(sendbuf_ind, data_size*3, MPI_INT, rbuf_ind, data_size*3, MPI_INT, 0, MPI_COMM_WORLD);
               
               if (myDomain == 0)
               {
                   /* rearrange order of the phase field value */
                   for (int m = 0; m < mySize; m++)
                       for (int n = 0; n < num_fields; n++)
                           for (int i= 0; i < phase_size; i++)
                           {
                               int ind = m*data_size + n*phase_size + i;
                               _PHI[n][rbuf_ind[ind*3]*_NY*_NZ + rbuf_ind[ind*3+1]*_NZ + rbuf_ind[ind*3+2]] = rbuf[ind];
                           }
               }
               
               MPI_Bcast(&(_PHI[0][0]),total_size,MPI_DOUBLE,0,MPI_COMM_WORLD);
               
               MPI_Barrier(MPI_COMM_WORLD);
#endif 

               int droplet_height = 0; 
               /* shift the phase field downwards */
               for(i=0;i<_NX;i++)
                   for(j=0;j<_NY;j++)
                       for(k=0;k<_NZ;k++)
                       {   
                           if((_PHI[0][i*_NY*_NZ+j*_NZ+k]>0.5) && (k > droplet_height))
                           {
                                droplet_height = k;
                            }
                       }
#ifndef _STK_MPI
               INFO_Printf("curstep = %6d droplet_height = %d \n",curstep,droplet_height);
#else
               
               if (myDomain == 0)
               {
                   INFO_Printf("curstep = %6d droplet_height = %d \n",curstep,droplet_height);
               }
#endif

	       if (_NZ - droplet_height <= 10)
	       {               
               /* shift the phase field downwards */
               for(i=0;i<_NX;i++)
                   for(j=0;j<_NY;j++)
                       for(k=0;k<_NZ-shiftheight;k++)
                       {   
                           _PHI[0][i*_NY*_NZ+j*_NZ+k] = _PHI[0][i*_NY*_NZ+j*_NZ+k+shiftheight];
                           _PHI[1][i*_NY*_NZ+j*_NZ+k] = _PHI[1][i*_NY*_NZ+j*_NZ+k+shiftheight];
                           _PHI[2][i*_NY*_NZ+j*_NZ+k] = _PHI[2][i*_NY*_NZ+j*_NZ+k+shiftheight];
                       }
               
               /* create vapor phase */
                for(i=0;i<_NX;i++)
                   for(j=0;j<_NY;j++)
                       for(k=_NZ-shiftheight;k<_NZ;k++)
                       {
                           _PHI[2][i*_NY*_NZ+j*_NZ+k] = 1 - 1e-6;
                           _PHI[0][i*_NY*_NZ+j*_NZ+k] = 2e-6;
                           _PHI[1][i*_NY*_NZ+j*_NZ+k] = 1.0 - _PHI[0][i*_NY*_NZ+j*_NZ+k] - _PHI[2][i*_NY*_NZ+j*_NZ+k];
                       }
               }
#ifdef _STK_MPI               
               /* reassign phase field values back to each node */
                for (n = 0; n < num_fields; n++)
                   for (i = x_begin; i < x_end; i++)
                       for (j = y_begin; j < y_end; j++)
                           for (k = z_begin; k < z_end; k++)
                               (*s_PHI[n])(k, j, i) = _PHI[n][i*_NY*_NZ+j*_NZ+k];
#endif
               
           }
        }
        
    
        
        /* print information */
        if(curstep%printfreq==0)
        {
#ifdef _STK_MPI
	if (myDomain == 0)
	{
            if (num_fields == 1)
              INFO_Printf("curstep = %6d F = %18.12e  Fraw = %18.12e G = %12.6e timestep = %8.3e rel_vol = %5.3g%%\n",
                          curstep,_F,_Fraw, _G, timestep, phase_volume[0]*100);
            else if (num_fields == 2)
              INFO_Printf("curstep = %6d F = %20.12e Fraw = %18.12e rel_vol = (%5.3g,%5.3g)%% \n",curstep,_F,_Fraw,
                           phase_volume[0]*100, phase_volume[1]*100);
            else if (num_fields == 3) {
              INFO_Printf("curstep = %6d F = %20.12e Fraw = %20.12e G = %20.12e timestep = %5.3g \n"
                          " rel_vol = (%5.3g,%5.3g,%5.3g)%% M01=%5.3g M02=%5.3g COM_x=%5.3g COM_y=%5.3g\n",
                           curstep,_F, _Fraw, _G, timestep,
                           phase_volume[0]*100, phase_volume[1]*100, phase_volume[2]*100,_M[0][1],_M[0][2],com_x, com_y);        
	    }
	    else
              INFO_Printf("curstep = %6d F = %20.12e\n",curstep,_F);
	}
#else
            if (num_fields == 1)
              INFO_Printf("curstep = %6d F = %18.12e  Fraw = %18.12e G = %12.6e timestep = %8.3e rel_vol = %5.3g%%\n",
                          curstep,_F,_Fraw, _G, timestep, phase_volume[0]*100);
            else if (num_fields == 2)
              INFO_Printf("curstep = %6d F = %20.12e Fraw = %18.12e rel_vol = (%5.3g,%5.3g)%% \n",curstep,_F,_Fraw,
                           phase_volume[0]*100, phase_volume[1]*100);
            else if (num_fields == 3) {
              INFO_Printf("curstep = %6d F = %20.12e Fraw = %20.12e G = %20.12e timestep = %5.3g \n"
                          " rel_vol = (%5.3g,%5.3g,%5.3g)%% M01=%5.3g M02=%5.3g COM_x=%5.3g COM_y=%5.3g\n",
                           curstep,_F, _Fraw, _G, timestep,
                           phase_volume[0]*100, phase_volume[1]*100, phase_volume[2]*100,_M[0][1],_M[0][2],com_x, com_y);        
	    }
	    else
              INFO_Printf("curstep = %6d F = %20.12e\n",curstep,_F);
#endif
     	}

	if (saveprop)
	{
#ifdef _STK_MPI
	if (myDomain == 0)
	{
	    if((curstep%savepropfreq)==0) pf.write(this);
	}
#else
	    if((curstep%savepropfreq)==0) pf.write(this);
#endif
	}

    if(savecn)
	{
#ifdef _STK_MPI
            if((curstep%savecnfreq)==0&&(curstep!=0))
            {   
    
                int i_pack = 0;

                for (n = 0; n < num_fields; n++)
                for (i = x_begin; i < x_end; i++)
                    for (j = y_begin; j < y_end; j++)
                    for (k = z_begin; k < z_end; k++)
                    {        
                        sendbuf[i_pack] = (*s_PHI[n])(k, j, i);                   
                            sendbuf_ind[i_pack*3] = i;
                            sendbuf_ind[i_pack*3+1] = j;
                            sendbuf_ind[i_pack*3+2] = k;                   
                            i_pack++;
                    }

                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Gather(sendbuf, data_size, MPI_DOUBLE, rbuf, data_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                MPI_Gather(sendbuf_ind, data_size*3, MPI_INT, rbuf_ind, data_size*3, MPI_INT, 0, MPI_COMM_WORLD);
                     

                if (myDomain == 0)
                {
                    for (int m = 0; m < mySize; m++)
                        for (int n = 0; n < num_fields; n++)
                            for (int i= 0; i < phase_size; i++)
                        {
                            int ind = m*data_size + n*phase_size + i;
                            _PHI[n][rbuf_ind[ind*3]*_NY*_NZ + rbuf_ind[ind*3+1]*_NZ + rbuf_ind[ind*3+2]] = rbuf[ind];
                        }
    		     }
                MPI_Barrier(MPI_COMM_WORLD);
              
                if (myDomain == 0)
                {
                  intercn.open(intercnfile);
               	  intercn.write(this);
                }
                MPI_Barrier(MPI_COMM_WORLD);
        }
#else            
	    if((curstep%savecnfreq)==0&&(curstep!=0))
            {
                intercn.open(intercnfile);
                intercn.write(this);
            }
#endif
	}
#ifndef _STK_MPI
        winplot();
#endif
    }


#ifdef _STK_MPI      
    /* pack all data into an ordered 1D array */
    MPI_Barrier(MPI_COMM_WORLD);
    
    int i_pack = 0;
    
    for (n = 0; n < num_fields; n++)
        for (i = x_begin; i < x_end; i++)
            for (j = y_begin; j < y_end; j++)
                for (k = z_begin; k < z_end; k++)
                {        
                    sendbuf[i_pack] = (*s_PHI[n])(k, j, i);
                    
                    sendbuf_ind[i_pack*3] = i;
                    sendbuf_ind[i_pack*3+1] = j;
                    sendbuf_ind[i_pack*3+2] = k;
                    
                    i_pack++;
                }
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(sendbuf, data_size, MPI_DOUBLE, rbuf, data_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(sendbuf_ind, data_size*3, MPI_INT, rbuf_ind, data_size*3, MPI_INT, 0, MPI_COMM_WORLD);
    
    if (myDomain == 0)
    {
        /* rearrange order of the phase field value */
        for (int m = 0; m < mySize; m++)
            for (int n = 0; n < num_fields; n++)
                for (int i= 0; i < phase_size; i++)
            {
                int ind = m*data_size + n*phase_size + i;
                _PHI[n][rbuf_ind[ind*3]*_NY*_NZ + rbuf_ind[ind*3+1]*_NZ + rbuf_ind[ind*3+2]] = rbuf[ind];
            }
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif 
}



/* compute free energy of current phase field config */
double PhaseFieldFrame::FreeEnergy()
{
    if (model_type == 0) 
    {
        /* single phase field model */
#ifdef _STK_MPI
        return FreeEnergy_single_mpi();
#elif _OPENMP
        return FreeEnergy_single_omp();
#else
        return FreeEnergy_single();
#endif
    }
    else if (model_type == 1)
    {
        /* multi phase field model */
#ifdef _STK_MPI
        return FreeEnergy_multi_mpi();
#elif _OPENMP
        return FreeEnergy_multi_omp();
#else
        return FreeEnergy_multi();
#endif
    }
#ifdef _USECUDA
    else if (model_type == 10) 
    {
        /* single phase field model */
        return cuda_FreeEnergy_single();
    }
    else if (model_type == 20)
    {
        /* multi-phase field model */
        return cuda_FreeEnergy_multi();
    }
#endif
    else
    {
        FATAL("FreeEnergy: unknown model_type " << model_type);
        return 0;
    }
}

double PhaseFieldFrame::FreeEnergy_single() /* traditional (single) phase field model */
{
    int n, i, j, k, ip, im, jp, jm, kp, km, ind;
    double h2, f, g, phi, dphidx, dphidy, dphidz, d2phi, dphi_SQR;
    double absdphi, EPS_loc, EPS_prefactor, triplesum;
    double nx_lab, ny_lab, nz_lab, nx_cryst, ny_cryst, nz_cryst;
    double dnx_lab_dphidx, dny_lab_dphidx, dnz_lab_dphidx; 
    double dnx_lab_dphidy, dny_lab_dphidy, dnz_lab_dphidy; 
    double dnx_lab_dphidz, dny_lab_dphidz, dnz_lab_dphidz; 
    double dEPS_dnx_lab, dEPS_dny_lab, dEPS_dnz_lab;
    double dEPS_dnx_cryst, dEPS_dny_cryst, dEPS_dnz_cryst;
    double dEPS_dphidx, dEPS_dphidy, dEPS_dphidz;
    double avg_dPHIdt0;

    _F = 0;
    _Fraw = 0;

    if (gridsize<=0)
      FATAL("FreeEnergy: gridesize = "<<gridsize);

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
            phi = _PHI[n][ind];
            d2phi = (_PHI[n][ip*_NY*_NZ+j*_NZ+k] + _PHI[n][im*_NY*_NZ+j*_NZ+k]
                    +_PHI[n][i*_NY*_NZ+jp*_NZ+k] + _PHI[n][i*_NY*_NZ+jm*_NZ+k]
                    +_PHI[n][i*_NY*_NZ+j*_NZ+kp] + _PHI[n][i*_NY*_NZ+j*_NZ+km]
                    -_PHI[n][i*_NY*_NZ+j*_NZ+k] * 6.0 ) / h2;

            _d2PHI[n][ind] = d2phi;

            dphidx = (_PHI[n][ip*_NY*_NZ+j*_NZ+k] - _PHI[n][im*_NY*_NZ+j*_NZ+k]) / (2.0*gridsize);
            dphidy = (_PHI[n][i*_NY*_NZ+jp*_NZ+k] - _PHI[n][i*_NY*_NZ+jm*_NZ+k]) / (2.0*gridsize);
            dphidz = (_PHI[n][i*_NY*_NZ+j*_NZ+kp] - _PHI[n][i*_NY*_NZ+j*_NZ+km]) / (2.0*gridsize);

            _dPHIdx[n][ind] = dphidx;
            _dPHIdy[n][ind] = dphidy;
            _dPHIdz[n][ind] = dphidz;

            f = phi*phi*(1-phi)*(1-phi);
            g = 2.0*phi*(1-phi)*(1-2.0*phi);

            _F += f * _U[n][n]; 
            _dFdPHI[n][ind] = g * _U[n][n];

            if ((_EPS1[n][n] == 0)&&(_EPS2[n][n] == 0)&&(_EPS3[n][n] == 0))
            { /* isotropic interface energy */
              EPS_loc = _EPS[n][n];
              EPS_prefactor = 0.5*SQR(EPS_loc);

              _F += (dphidx*dphidx + dphidy*dphidy + dphidz*dphidz) * EPS_prefactor;
              _dFdPHI[n][ind] += - 2.0 * d2phi * EPS_prefactor;

            }
            else
            { /* cubic anisotropic interface energy */
              dphi_SQR   = SQR(dphidx) + SQR(dphidy) + SQR(dphidz);
              absdphi = sqrt(dphi_SQR + 1e-24);
              nx_lab = dphidx / absdphi;
              ny_lab = dphidy / absdphi;
              nz_lab = dphidz / absdphi;

#define _R _ROT_MATRIX 
              nx_cryst = _R[0][0]*nx_lab + _R[0][1]*ny_lab + _R[0][2]*nz_lab;
              ny_cryst = _R[1][0]*nx_lab + _R[1][1]*ny_lab + _R[1][2]*nz_lab;
              nz_cryst = _R[2][0]*nx_lab + _R[2][1]*ny_lab + _R[2][2]*nz_lab;

              dnx_lab_dphidx = 1./absdphi - dphidx * dphidx / CUBE(absdphi);
              dny_lab_dphidx =            - dphidy * dphidx / CUBE(absdphi);
              dnz_lab_dphidx =            - dphidz * dphidx / CUBE(absdphi);
   
              dnx_lab_dphidy =            - dphidx * dphidy / CUBE(absdphi);
              dny_lab_dphidy = 1./absdphi - dphidy * dphidy / CUBE(absdphi);
              dnz_lab_dphidy =            - dphidz * dphidy / CUBE(absdphi);
   
              dnx_lab_dphidz =            - dphidx * dphidz / CUBE(absdphi);
              dny_lab_dphidz =            - dphidy * dphidz / CUBE(absdphi);
              dnz_lab_dphidz = 1./absdphi - dphidz * dphidz / CUBE(absdphi);
   
              triplesum = SQR(nx_cryst)*SQR(ny_cryst) + SQR(ny_cryst)*SQR(nz_cryst) + SQR(nz_cryst)*SQR(nx_cryst);
              EPS_loc   = _EPS[n][n] + _EPS1[n][n]*triplesum + _EPS2[n][n]*SQR(nx_cryst*ny_cryst*nz_cryst)
                                     + _EPS3[n][n]*SQR(triplesum);
              _EPS_loc[n][ind]   = EPS_loc;

              _F += 0.5 * SQR(EPS_loc) * dphi_SQR;
              /* dFdPHI += -d/dx ( SQR(EPS_loc) * dphidx + dphi_SQR * EPS_loc * d(EPS_loc)/d(dphidx) )
                           -d/dy ( SQR(EPS_loc) * dphidy + dphi_SQR * EPS_loc * d(EPS_loc)/d(dphidy) )
                           -d/dz ( SQR(EPS_loc) * dphidz + dphi_SQR * EPS_loc * d(EPS_loc)/d(dphidz) ) ;
              */

              dEPS_dnx_cryst = 2.*_EPS1[n][n]*nx_cryst*(SQR(ny_cryst)+SQR(nz_cryst))
                             + 2.*_EPS2[n][n]*nx_cryst*SQR(ny_cryst)*SQR(nz_cryst)
                             + 4.*_EPS3[n][n]*nx_cryst*(SQR(ny_cryst)+SQR(nz_cryst))*triplesum;
              dEPS_dny_cryst = 2.*_EPS1[n][n]*ny_cryst*(SQR(nz_cryst)+SQR(nx_cryst))
                             + 2.*_EPS2[n][n]*ny_cryst*SQR(nz_cryst)*SQR(nx_cryst)
                             + 4.*_EPS3[n][n]*ny_cryst*(SQR(nz_cryst)+SQR(nx_cryst))*triplesum;
              dEPS_dnz_cryst = 2.*_EPS1[n][n]*nz_cryst*(SQR(nx_cryst)+SQR(ny_cryst))
                             + 2.*_EPS2[n][n]*nz_cryst*SQR(nx_cryst)*SQR(ny_cryst)
                             + 4.*_EPS3[n][n]*nz_cryst*(SQR(nx_cryst)+SQR(ny_cryst))*triplesum;
  
              dEPS_dnx_lab = _R[0][0]*dEPS_dnx_cryst + _R[1][0]*dEPS_dny_cryst + _R[2][0]*dEPS_dnz_cryst;
              dEPS_dny_lab = _R[0][1]*dEPS_dnx_cryst + _R[1][1]*dEPS_dny_cryst + _R[2][1]*dEPS_dnz_cryst;
              dEPS_dnz_lab = _R[0][2]*dEPS_dnx_cryst + _R[1][2]*dEPS_dny_cryst + _R[2][2]*dEPS_dnz_cryst;
 
              dEPS_dphidx = dEPS_dnx_lab*dnx_lab_dphidx + dEPS_dny_lab*dny_lab_dphidx + dEPS_dnz_lab*dnz_lab_dphidx;
              dEPS_dphidy = dEPS_dnx_lab*dnx_lab_dphidy + dEPS_dny_lab*dny_lab_dphidy + dEPS_dnz_lab*dnz_lab_dphidy;
              dEPS_dphidz = dEPS_dnx_lab*dnx_lab_dphidz + dEPS_dny_lab*dny_lab_dphidz + dEPS_dnz_lab*dnz_lab_dphidz;

              //tmp_x[n][ind] = SQR(EPS_loc)*dphidx + dphi_SQR*EPS_loc*dEPS_dphidx;
              //tmp_y[n][ind] = SQR(EPS_loc)*dphidy + dphi_SQR*EPS_loc*dEPS_dphidy;
              //tmp_z[n][ind] = SQR(EPS_loc)*dphidz + dphi_SQR*EPS_loc*dEPS_dphidz;

              /* new scheme for computing variational derivative, 2012/04/27  */
              tmp_x[n][ind] = dphi_SQR*EPS_loc*dEPS_dphidx;
              tmp_y[n][ind] = dphi_SQR*EPS_loc*dEPS_dphidy;
              tmp_z[n][ind] = dphi_SQR*EPS_loc*dEPS_dphidz;

              _dFdPHI[n][ind] += - 1.0 * SQR(EPS_loc) * d2phi;
            }
          }
        }
     }

     if ((_EPS1[n][n] != 0)||(_EPS2[n][n] != 0)||(_EPS3[n][n] != 0))
     { /* cubic anisotropic interface energy */
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
            _dFdPHI[n][ind] +=
                  - (tmp_x[n][ip*_NY*_NZ+j*_NZ+k]-tmp_x[n][im*_NY*_NZ+j*_NZ+k]) / (2.0*gridsize)
                  - (tmp_y[n][i*_NY*_NZ+jp*_NZ+k]-tmp_y[n][i*_NY*_NZ+jm*_NZ+k]) / (2.0*gridsize)
                  - (tmp_z[n][i*_NY*_NZ+j*_NZ+kp]-tmp_z[n][i*_NY*_NZ+j*_NZ+km]) / (2.0*gridsize) ;

            /* new scheme for computing variational derivative, 2012/04/27  */
            _dFdPHI[n][ind] +=
                - (SQR(_EPS_loc[n][ip*_NY*_NZ+j*_NZ+k])-SQR(_EPS_loc[n][im*_NY*_NZ+j*_NZ+k])) / (2.0*gridsize) * _dPHIdx[n][ind]
                - (SQR(_EPS_loc[n][i*_NY*_NZ+jp*_NZ+k])-SQR(_EPS_loc[n][i*_NY*_NZ+jm*_NZ+k])) / (2.0*gridsize) * _dPHIdy[n][ind]
                - (SQR(_EPS_loc[n][i*_NY*_NZ+j*_NZ+kp])-SQR(_EPS_loc[n][i*_NY*_NZ+j*_NZ+km])) / (2.0*gridsize) * _dPHIdz[n][ind];
          }
        }
      }
     }

     _F *= CUBE(gridsize);
     _Fraw = _F;

     if (dynamics_type == 1) /* Cahn-Hillard equation */
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
            _d2dFdPHI[n][ind] = 
                    (_dFdPHI[n][ip*_NY*_NZ+j*_NZ+k] + _dFdPHI[n][im*_NY*_NZ+j*_NZ+k]
                    +_dFdPHI[n][i*_NY*_NZ+jp*_NZ+k] + _dFdPHI[n][i*_NY*_NZ+jm*_NZ+k]
                    +_dFdPHI[n][i*_NY*_NZ+j*_NZ+kp] + _dFdPHI[n][i*_NY*_NZ+j*_NZ+km]
                    -_dFdPHI[n][ind] * 6.0 ) / h2;
          }
        }
     }

     if ((dynamics_type == 0) || (dynamics_type == 2)) /* Cahn-Hillard equation */
     {
              for(ind=0;ind<_NX*_NY*_NZ;ind++)
                 _dPHIdt0[n][ind] = -1.0 * Mob_GL * _dFdPHI[n][ind];

     } 
     else if (dynamics_type == 1) /* Ginzburg-Landau equation */
     {
              for(ind=0;ind<_NX*_NY*_NZ;ind++)
                 _dPHIdt0[n][ind] = -1.0 * Mob_D * _d2dFdPHI[n][ind];
     }

     if ((dynamics_type == 0) || (dynamics_type == 1)) /* no additional constraints */
     {
       /* 0: Ginzburg-Landau equation: not conserved */
       /* 1: Cahn-Hilliard equation: conserved */
       for(ind=0;ind<_NX*_NY*_NZ;ind++)
              _dPHIdt[n][ind] = _dPHIdt0[n][ind];
     }
     else if (dynamics_type == 2) 
     {
        /* Ginzburg-Landau equation but with constant volume constraint */ 
        avg_dPHIdt0 = 0;

        for(ind=0;ind<_NX*_NY*_NZ;ind++)
              avg_dPHIdt0 += _dPHIdt0[n][ind];

        avg_dPHIdt0 /= (_NX*_NY*_NZ);

        for(ind=0;ind<_NX*_NY*_NZ;ind++)
              _dPHIdt[n][ind] = _dPHIdt0[n][ind] - avg_dPHIdt0;
     }
     else
     {
        ERROR("unknown dynamics_type = "<<dynamics_type);
     }

    }/* end of for(n=0;n<num_fields;n++) */

        /* calculate the max value of dphidt */
        double Gmiddle = 0;
        for(n=0;n<num_fields;n++)
        {
          for(ind=0;ind<_NX*_NY*_NZ;ind++)
          {
               if ( Gmiddle <= fabs(_dPHIdt[n][ind]) )  Gmiddle = fabs(_dPHIdt[n][ind]);
          }
        }
        _G = Gmiddle;
        
        /* adjust time step at every step */
        timestep = dtmax/fmax(_G*dtmax/dphimax,1);
    
    return _F;
}

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////


double PhaseFieldFrame::FreeEnergy_multi() /* multi-phase field model (MPF) */
{
    int n, p, i, j, k, ip, ip2, ip3, im, im2, im3, jp, jp2, jp3, jm, jm2, jm3, kp, kp2, kp3, km, km2, km3, ind;
    double h2, f, g, phi, dphidx, dphidy, dphidz, d2phi, dphi_SQR, dmm, dmu1, fp;
    double absdphi, EPS_loc, EPS_prefactor, triplesum;
    double nx_lab, ny_lab, nz_lab, nx_cryst, ny_cryst, nz_cryst;
    double dnx_lab_dphidx, dny_lab_dphidx, dnz_lab_dphidx; 
    double dnx_lab_dphidy, dny_lab_dphidy, dnz_lab_dphidy; 
    double dnx_lab_dphidz, dny_lab_dphidz, dnz_lab_dphidz; 
    double dEPS_dnx_lab, dEPS_dny_lab, dEPS_dnz_lab;
    double dEPS_dnx_cryst, dEPS_dny_cryst, dEPS_dnz_cryst;
    double dEPS_dphidx, dEPS_dphidy, dEPS_dphidz;
    double eps1, eps2, eps3, dEPS_dphidx_loc, dEPS_dphidy_loc, dEPS_dphidz_loc;
    double grad_loc_x, grad_loc_y, grad_loc_z, grad_loc_SQR;
    double grad_term, avg_dPHIdt0;
    double A_11, A_12, A_13, A_14, A_21, A_22, A_23, A_24, A_31, A_32, A_33, A_34, A_41, A_42, A_43, A_44;
    double B_1, B_2, B_3, B_4;
    double dmuL, dmuV, dfX, dfY, divider;
    double vol_incr, x_incr, y_incr;
    
    _F = 0;
    _Fraw = 0;

    if (gridsize<=0)
      FATAL("FreeEnergy: gridesize = "<<gridsize);
   
    if (num_fields!=3)
      FATAL("FreeEnergy: num_field = "<<num_fields);
    
    /* Zero out time derivative of phase fields */
    for(n=0;n<num_fields;n++)
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {
            _dPHIdt0[n][ind] = 0;
            _dFdPHIi[n][ind] = 0;
            _dPHIdtCH[n][ind] = 0;
        }
    
    /* Reset liquid chemical potential, will be updated by constraint */
    if (dynamics_type == 2 || dynamics_type == 3 || dynamics_type == 8)
    {
        _M[0][1] = 0; _M[0][2] = _M[1][2];
        _M[1][0] = 0; _M[2][0] = _M[2][1];
        _MU[0] = _MU[1];
    }

    if (dynamics_type == 4 || dynamics_type == 5 || dynamics_type == 6 || dynamics_type == 7)
    {
        // _MU[1] = 0;
        _MU[0] = _MU[1];
        _MU[2] = _MU[1];
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
                    phi = _PHI[n][ind];
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
     
 
    /* Compute Free energy and time derivative of phase fields */
    
    /* assign the interface mobility and assigan the gradient matrix for external force term*/
    
    for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {
            _K_D[ind] = Mob_D;
        }          
    
    for (i=0;i<_NX;i++)
    {
        ip  = (i+1)%_NX;
        ip2 = (i+2)%_NX;
        ip3 = (i+3)%_NX;
        im  = (i-1+_NX)%_NX;
        im2 = (i-2+_NX)%_NX;
        im3 = (i-3+_NX)%_NX;
        for (j=0;j<_NY;j++)
        {
            jp  = (j+1)%_NY;
            jp2 = (j+2)%_NY;
            jp3 = (j+3)%_NY;
            jm  = (j-1+_NY)%_NY;
            jm2 = (j-2+_NY)%_NY;
            jm3 = (j-3+_NY)%_NY;
            for (k=0;k<_NZ;k++)
            {
                kp  = (k+1)%_NZ;
                kp2 = (k+2)%_NZ;
                kp3 = (k+3)%_NZ;
                km  = (k-1+_NZ)%_NZ;
                km2 = (k-2+_NZ)%_NZ; 
                km3 = (k-3+_NZ)%_NZ;
                ind = i*_NY*_NZ+j*_NZ+k;
                
                /* define the solid phase boundary */
                if (     (((_PHI[1][ind]-0.5)*(_PHI[1][ip*_NY*_NZ+j*_NZ+k]-0.5)<=0)
                       || ((_PHI[1][ind]-0.5)*(_PHI[1][im*_NY*_NZ+j*_NZ+k]-0.5)<=0)
                       || ((_PHI[1][ind]-0.5)*(_PHI[1][ip2*_NY*_NZ+j*_NZ+k]-0.5)<=0)
                       || ((_PHI[1][ind]-0.5)*(_PHI[1][im2*_NY*_NZ+j*_NZ+k]-0.5)<=0)
                       || ((_PHI[1][ind]-0.5)*(_PHI[1][ip3*_NY*_NZ+j*_NZ+k]-0.5)<=0)
                       || ((_PHI[1][ind]-0.5)*(_PHI[1][im3*_NY*_NZ+j*_NZ+k]-0.5)<=0)
                       || ((_PHI[1][ind]-0.5)*(_PHI[1][i*_NY*_NZ+jp*_NZ+k]-0.5)<=0)
                       || ((_PHI[1][ind]-0.5)*(_PHI[1][i*_NY*_NZ+jm*_NZ+k]-0.5)<=0)
                       || ((_PHI[1][ind]-0.5)*(_PHI[1][i*_NY*_NZ+jp2*_NZ+k]-0.5)<=0)
                       || ((_PHI[1][ind]-0.5)*(_PHI[1][i*_NY*_NZ+jm2*_NZ+k]-0.5)<=0)
     	               || ((_PHI[1][ind]-0.5)*(_PHI[1][i*_NY*_NZ+jp3*_NZ+k]-0.5)<=0)
                       || ((_PHI[1][ind]-0.5)*(_PHI[1][i*_NY*_NZ+jm3*_NZ+k]-0.5)<=0)
                       || ((_PHI[1][ind]-0.5)*(_PHI[1][i*_NY*_NZ+j*_NZ+kp]-0.5)<=0)
                       || ((_PHI[1][ind]-0.5)*(_PHI[1][i*_NY*_NZ+j*_NZ+km]-0.5)<=0)
                       || ((_PHI[1][ind]-0.5)*(_PHI[1][i*_NY*_NZ+j*_NZ+kp2]-0.5)<=0)
                       || ((_PHI[1][ind]-0.5)*(_PHI[1][i*_NY*_NZ+j*_NZ+km2]-0.5)<=0)
                       || ((_PHI[1][ind]-0.5)*(_PHI[1][i*_NY*_NZ+j*_NZ+kp3]-0.5)<=0)
                       || ((_PHI[1][ind]-0.5)*(_PHI[1][i*_NY*_NZ+j*_NZ+km3]-0.5)<=0)) )
                { 
                    _K_S[ind] = 1;
                }
                else
                {
                    _K_S[ind] = 0;
                }
                
                /* define the liquid phase boundary */
                if (        (((_PHI[0][ind]-0.5)*(_PHI[0][ip*_NY*_NZ+j*_NZ+k]-0.5)<=0)
                          || ((_PHI[0][ind]-0.5)*(_PHI[0][im*_NY*_NZ+j*_NZ+k]-0.5)<=0)
                          || ((_PHI[0][ind]-0.5)*(_PHI[0][ip2*_NY*_NZ+j*_NZ+k]-0.5)<=0)
                          || ((_PHI[0][ind]-0.5)*(_PHI[0][im2*_NY*_NZ+j*_NZ+k]-0.5)<=0)
      	                  || ((_PHI[0][ind]-0.5)*(_PHI[0][ip3*_NY*_NZ+j*_NZ+k]-0.5)<=0)
                          || ((_PHI[0][ind]-0.5)*(_PHI[0][im3*_NY*_NZ+j*_NZ+k]-0.5)<=0)
                          || ((_PHI[0][ind]-0.5)*(_PHI[0][i*_NY*_NZ+jp*_NZ+k]-0.5)<=0)
                          || ((_PHI[0][ind]-0.5)*(_PHI[0][i*_NY*_NZ+jm*_NZ+k]-0.5)<=0)
                          || ((_PHI[0][ind]-0.5)*(_PHI[0][i*_NY*_NZ+jp2*_NZ+k]-0.5)<=0)
                          || ((_PHI[0][ind]-0.5)*(_PHI[0][i*_NY*_NZ+jm2*_NZ+k]-0.5)<=0)
                          || ((_PHI[0][ind]-0.5)*(_PHI[0][i*_NY*_NZ+jp3*_NZ+k]-0.5)<=0)
                          || ((_PHI[0][ind]-0.5)*(_PHI[0][i*_NY*_NZ+jm3*_NZ+k]-0.5)<=0)
                          || ((_PHI[0][ind]-0.5)*(_PHI[0][i*_NY*_NZ+j*_NZ+kp]-0.5)<=0)
                          || ((_PHI[0][ind]-0.5)*(_PHI[0][i*_NY*_NZ+j*_NZ+km]-0.5)<=0)
                          || ((_PHI[0][ind]-0.5)*(_PHI[0][i*_NY*_NZ+j*_NZ+kp2]-0.5)<=0)
                          || ((_PHI[0][ind]-0.5)*(_PHI[0][i*_NY*_NZ+j*_NZ+km2]-0.5)<=0) 
                          || ((_PHI[0][ind]-0.5)*(_PHI[0][i*_NY*_NZ+j*_NZ+kp3]-0.5)<=0)
                          || ((_PHI[0][ind]-0.5)*(_PHI[0][i*_NY*_NZ+j*_NZ+km3]-0.5)<=0)))
                { 
                    _K_L[ind] = 1;
                }
                else
                {
                    _K_L[ind] = 0;
                }
                /* define the l-s interface */
                if (( _K_L[ind]==1 ) && ( _K_S[ind]==1 )) 
                {
                    _K_LS[ind] = 1;
                }
                else
                {
                    _K_LS[ind] = 0;
                }
                /* define the s-v interface */
                if (( _K_L[ind]==0 ) && ( _K_S[ind]==1 ))
                {
                    _K_SV[ind] = 1;
                }
                else
                {
                    _K_SV[ind] = 0;
                }
                /* define the l-v interface */
                if (( _K_L[ind]==1 ) && ( _K_S[ind]==0 ))
                {
                    _K_LV[ind] = 1;
                }
                else
                {
                    _K_LV[ind] = 0;
                }
                
                _K_other[ind] = 1 - _K_LV[ind] - _K_SV[ind] - _K_LS[ind];
                
            }
        }
    }
      
    
    
    for(n=0;n<num_fields;n++)
    {
          for (ind=0; ind<_NX*_NY*_NZ; ind++)
          {
                grad_term = SQR(_dPHIdx[n][ind]) + SQR(_dPHIdy[n][ind]) + SQR(_dPHIdz[n][ind]) ;
                
                f = SQR(_PHI[n][ind]) * SQR(1-_PHI[n][ind]);
                //fp = _PHI[n][ind];
                fp = tanh((_PHI[n][ind]-0.5)/0.10)/2 + 0.5;
                _F += f * _U[n][n] + fp * _MU[n];
                _Fraw += f * _U[n][n];
                /* the influence of the external force (only for liquid) */ 
                if ( n == 0 )
                { 
                   _F += Fext_x * fp * matrix_x[ind] + Fext_y * fp * matrix_y[ind]; 
                }
                 
                // Penalty term for three phase co-existing    
                _F += Penalty_3phase/num_fields*SQR(_PHI[0][ind])*SQR(_PHI[1][ind])*SQR(_PHI[2][ind]);
                 
                // Penalty term for deviation from original NW shape
                _F += Penalty_NW/num_fields*SQR(_PHI[1][ind]-_NW_orig[ind]);
               
            if ( ((_EPS1[n][n] == 0)&&(_EPS2[n][n] == 0)&&(_EPS3[n][n] == 0)) )
            { /* isotropic interface energy */
                EPS_loc = 1.0;
                EPS_prefactor = SQR(EPS_loc);
                
                _F += _EPS[n][n] * EPS_prefactor * grad_term;
                _Fraw += _EPS[n][n] * EPS_prefactor * grad_term;

                /* the equation of motion should follow the Steinbach 1999 formulation */
                _dFdPHIi[n][ind] = _U[n][n]*(4.0*_PHI[n][ind]*_PHI[n][ind]*_PHI[n][ind]
                               + 2.0*_PHI[n][ind] - 6.0*_PHI[n][ind]*_PHI[n][ind]) 
                               + _MU[n]*(5.0-5.0*SQR(tanh(10*_PHI[n][ind]-5.0)))
                               - 2.0*_EPS[n][n]*EPS_prefactor*_d2PHI[n][ind]; 
              
                /* with the choice of C(phi) = phi */  
                /* _dFdPHIi[n][ind] = _U[n][n]*(4.0*_PHI[n][ind]*_PHI[n][ind]*_PHI[n][ind]
                               + 2.0*_PHI[n][ind] - 6.0*_PHI[n][ind]*_PHI[n][ind]) 
                               + _MU[n] - 2.0*_EPS[n][n]*EPS_prefactor*_d2PHI[n][ind]; */
                
                /* with the external force */
                if ( n==0 )
                { _dFdPHIi[0][ind] += Fext_x * (5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))*matrix_x[ind]
                                    + Fext_y * (5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))*matrix_y[ind];
                }
                
            }
            else
            {   /* anisotropic interface energy */                          
                /* the EPS term will be computed and added below */
                /* surface normal defined by solid phase gradient */
                dphidx = _dPHIdx[n][ind];
                dphidy = _dPHIdy[n][ind];
                dphidz = _dPHIdz[n][ind];

                dphi_SQR   = SQR(dphidx) + SQR(dphidy) + SQR(dphidz);
                absdphi = sqrt(dphi_SQR + 1e-24);

                nx_lab = dphidx / absdphi;
                ny_lab = dphidy / absdphi;
                nz_lab = dphidz / absdphi;

                nx_cryst = _R[0][0]*nx_lab + _R[0][1]*ny_lab + _R[0][2]*nz_lab;
                ny_cryst = _R[1][0]*nx_lab + _R[1][1]*ny_lab + _R[1][2]*nz_lab;
                nz_cryst = _R[2][0]*nx_lab + _R[2][1]*ny_lab + _R[2][2]*nz_lab;

                eps1 = _EPS1[n][n]; /* for solid(1) phase only */
                eps2 = _EPS2[n][n]; /* for solid(1) phase only */
                eps3 = _EPS3[n][n]; /* for solid(1) phase only */

                triplesum = SQR(nx_cryst)*SQR(ny_cryst) + SQR(ny_cryst)*SQR(nz_cryst) + SQR(nz_cryst)*SQR(nx_cryst);

                /* Notice the difference from the Matlab implementation */
                EPS_loc   = 1.0 + eps1*triplesum + eps2*SQR(nx_cryst*ny_cryst*nz_cryst) + eps3*SQR(triplesum);
                _EPS_loc[n][ind] = EPS_loc;

                /* Free energy for anisotropic surface energy */                
                _F += _EPS[n][n] * SQR(_EPS_loc[n][ind]) * grad_term;
                _Fraw += _EPS[n][n] * SQR(_EPS_loc[n][ind]) * grad_term;
                
                dnx_lab_dphidx = 1./absdphi - dphidx * dphidx / CUBE(absdphi);
                dny_lab_dphidx =            - dphidy * dphidx / CUBE(absdphi);
                dnz_lab_dphidx =            - dphidz * dphidx / CUBE(absdphi);
   
                dnx_lab_dphidy =            - dphidx * dphidy / CUBE(absdphi);
                dny_lab_dphidy = 1./absdphi - dphidy * dphidy / CUBE(absdphi);
                dnz_lab_dphidy =            - dphidz * dphidy / CUBE(absdphi);
   
                dnx_lab_dphidz =            - dphidx * dphidz / CUBE(absdphi);
                dny_lab_dphidz =            - dphidy * dphidz / CUBE(absdphi);
                dnz_lab_dphidz = 1./absdphi - dphidz * dphidz / CUBE(absdphi);
   
                dEPS_dnx_cryst = 2.*eps1*nx_cryst*(SQR(ny_cryst)+SQR(nz_cryst))
                               + 2.*eps2*nx_cryst*SQR(ny_cryst)*SQR(nz_cryst)
                               + 4.*eps3*nx_cryst*(SQR(ny_cryst)+SQR(nz_cryst))*triplesum;
                dEPS_dny_cryst = 2.*eps1*ny_cryst*(SQR(nz_cryst)+SQR(nx_cryst))
                               + 2.*eps2*ny_cryst*SQR(nz_cryst)*SQR(nx_cryst)
                               + 4.*eps3*ny_cryst*(SQR(nz_cryst)+SQR(nx_cryst))*triplesum;
                dEPS_dnz_cryst = 2.*eps1*nz_cryst*(SQR(nx_cryst)+SQR(ny_cryst))
                               + 2.*eps2*nz_cryst*SQR(nx_cryst)*SQR(ny_cryst)
                               + 4.*eps3*nz_cryst*(SQR(nx_cryst)+SQR(ny_cryst))*triplesum;
    
                dEPS_dnx_lab = _R[0][0]*dEPS_dnx_cryst + _R[1][0]*dEPS_dny_cryst + _R[2][0]*dEPS_dnz_cryst;
                dEPS_dny_lab = _R[0][1]*dEPS_dnx_cryst + _R[1][1]*dEPS_dny_cryst + _R[2][1]*dEPS_dnz_cryst;
                dEPS_dnz_lab = _R[0][2]*dEPS_dnx_cryst + _R[1][2]*dEPS_dny_cryst + _R[2][2]*dEPS_dnz_cryst;
 
                dEPS_dphidx = dEPS_dnx_lab*dnx_lab_dphidx + dEPS_dny_lab*dny_lab_dphidx + dEPS_dnz_lab*dnz_lab_dphidx;
                dEPS_dphidy = dEPS_dnx_lab*dnx_lab_dphidy + dEPS_dny_lab*dny_lab_dphidy + dEPS_dnz_lab*dnz_lab_dphidy;
                dEPS_dphidz = dEPS_dnx_lab*dnx_lab_dphidz + dEPS_dny_lab*dny_lab_dphidz + dEPS_dnz_lab*dnz_lab_dphidz;

                                                      
              
                dEPS_dphidx_loc = dEPS_dphidx;
                dEPS_dphidy_loc = dEPS_dphidy;
                dEPS_dphidz_loc = dEPS_dphidz;
                

  
                grad_loc_x = _dPHIdx[n][ind];
                grad_loc_y = _dPHIdy[n][ind];
                grad_loc_z = _dPHIdz[n][ind];
                grad_loc_SQR = SQR(grad_loc_x) + SQR(grad_loc_y) + SQR(grad_loc_z);
  
                /* new scheme for computing variational derivative, 2012/04/27  */
                tmp_x[n][ind] = grad_loc_SQR*EPS_loc*dEPS_dphidx_loc;
                tmp_y[n][ind] = grad_loc_SQR*EPS_loc*dEPS_dphidy_loc;
                tmp_z[n][ind] = grad_loc_SQR*EPS_loc*dEPS_dphidz_loc;
                
                _dFdPHIi[n][ind] += _U[n][n]*(4*_PHI[n][ind]*_PHI[n][ind]*_PHI[n][ind]
                                  + 2*_PHI[n][ind] - 6*_PHI[n][ind]*_PHI[n][ind]) 
                                  + _MU[n]*(5.0-5.0*SQR(tanh(10*_PHI[n][ind]-5.0))); 

                /* for the choice of C(phi) = phi */
                /* _dFdPHIi[n][ind] += _U[n][n]*(4*_PHI[n][ind]*_PHI[n][ind]*_PHI[n][ind]
                                  + 2*_PHI[n][ind] - 6*_PHI[n][ind]*_PHI[n][ind]) 
                                  + _MU[n]; */

                _dFdPHIi[n][ind] += (-2.0) * _EPS[n][n]*SQR(EPS_loc) *(_d2PHI[n][ind]);
                
                // with the external force term
                if ( n==0 )
                { _dFdPHIi[0][ind] += Fext_x * (5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))*matrix_x[ind]
                    + Fext_y * (5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))*matrix_y[ind];
                }
            }

         }

          
       if ((_EPS1[n][n] != 0)||(_EPS2[n][n] != 0)||(_EPS3[n][n] != 0))
       { /* cubic anisotropic interface energy */
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
              _dFdPHIi[n][ind] += _EPS[n][n]*(-2.0)*(tmp_x[n][ip*_NY*_NZ+j*_NZ+k]-tmp_x[n][im*_NY*_NZ+j*_NZ+k]) / (2.0*gridsize)
                                + _EPS[n][n]*(-2.0)*(tmp_y[n][i*_NY*_NZ+jp*_NZ+k]-tmp_y[n][i*_NY*_NZ+jm*_NZ+k]) / (2.0*gridsize)
                                + _EPS[n][n]*(-2.0)*(tmp_z[n][i*_NY*_NZ+j*_NZ+kp]-tmp_z[n][i*_NY*_NZ+j*_NZ+km]) / (2.0*gridsize);

              /* new scheme for computing variational derivative, 2012/04/27  */
              grad_loc_x = _dPHIdx[n][ind];
              grad_loc_y = _dPHIdy[n][ind];
              grad_loc_z = _dPHIdz[n][ind];


              _dFdPHIi[n][ind] += _EPS[n][n]*(-4.0)*_EPS_loc[n][ind]*(_EPS_loc[n][ip*_NY*_NZ+j*_NZ+k]
                                - _EPS_loc[n][im*_NY*_NZ+j*_NZ+k])/ (2.0*gridsize) * grad_loc_x
                                + _EPS[n][n]*(-4.0)*_EPS_loc[n][ind]*(_EPS_loc[n][i*_NY*_NZ+jp*_NZ+k]
                                - _EPS_loc[n][i*_NY*_NZ+jm*_NZ+k])/ (2.0*gridsize) * grad_loc_y                                 
                                + _EPS[n][n]*(-4.0)*_EPS_loc[n][ind]*(_EPS_loc[n][i*_NY*_NZ+j*_NZ+kp]
                                - _EPS_loc[n][i*_NY*_NZ+j*_NZ+km])/ (2.0*gridsize) * grad_loc_z;                
            
            }
          }
        }
      }
    }

    /* variational derivative due to the triple junction penalty term */
    for (ind=0; ind<_NX*_NY*_NZ; ind++)
    {  
        _dFdPHIi[0][ind] += 2*Penalty_3phase*_PHI[0][ind]*SQR(_PHI[1][ind])*SQR(_PHI[2][ind]);
        _dFdPHIi[1][ind] += 2*Penalty_3phase*_PHI[1][ind]*SQR(_PHI[0][ind])*SQR(_PHI[2][ind]);
        _dFdPHIi[2][ind] += 2*Penalty_3phase*_PHI[2][ind]*SQR(_PHI[0][ind])*SQR(_PHI[1][ind]);
        
        /* variational derivative due to the NW shape penalty term */
        _dFdPHIi[1][ind] += 2*Penalty_NW*(_PHI[1][ind]-_NW_orig[ind]);
    }
    
    /* equation of motion based on Steinback 1999 formulation */
    for (ind=0; ind<_NX*_NY*_NZ; ind++)
    {
        /* G-L equation of motion */
        _dPHIdt0[0][ind] = - _K_SV[ind]*(Mob_GL*(_dFdPHIi[0][ind] - _dFdPHIi[1][ind]) 
                           + Mob_GL*(_dFdPHIi[0][ind] - _dFdPHIi[2][ind])) 
                           - _K_LV[ind]*(Mob_GL*((_dFdPHIi[0][ind]+_dFdPHIi[2][ind])/2 - _dFdPHIi[1][ind]) 
                           + Mob_LV*(_dFdPHIi[0][ind] - _dFdPHIi[2][ind])) 
                           - _K_LS[ind]*(Mob_GL*((_dFdPHIi[0][ind]+_dFdPHIi[1][ind])/2 - _dFdPHIi[2][ind])
                           + Mob_LS*(_dFdPHIi[0][ind] - _dFdPHIi[1][ind])) 
                           - _K_other[ind]*(Mob_GL*(_dFdPHIi[0][ind] - _dFdPHIi[1][ind])
                           + Mob_GL*(_dFdPHIi[0][ind] - _dFdPHIi[2][ind]));
                        
        _dPHIdt0[1][ind] = - _K_LV[ind]*(Mob_GL*(_dFdPHIi[1][ind] - _dFdPHIi[2][ind])
                           + Mob_GL*(_dFdPHIi[1][ind] - _dFdPHIi[0][ind])) 
                           - _K_SV[ind]*(Mob_GL*((_dFdPHIi[1][ind]+_dFdPHIi[2][ind])/2 - _dFdPHIi[0][ind])
                           + Mob_SV*(_dFdPHIi[1][ind] - _dFdPHIi[2][ind]))
                           - _K_LS[ind]*(Mob_GL*((_dFdPHIi[0][ind]+_dFdPHIi[1][ind])/2 - _dFdPHIi[2][ind])
                           + Mob_LS*(_dFdPHIi[1][ind] - _dFdPHIi[0][ind])) 
                           - _K_other[ind]*(Mob_GL*(_dFdPHIi[1][ind] - _dFdPHIi[0][ind])
                           + Mob_GL*(_dFdPHIi[1][ind] - _dFdPHIi[2][ind]));
                    
        _dPHIdt0[2][ind] = - _K_LS[ind]*(Mob_GL*(_dFdPHIi[2][ind] - _dFdPHIi[1][ind])
                           + Mob_GL*(_dFdPHIi[2][ind] - _dFdPHIi[0][ind])) 
                           - _K_LV[ind]*(Mob_GL*((_dFdPHIi[0][ind]+_dFdPHIi[2][ind])/2 - _dFdPHIi[1][ind])
                           + Mob_LV*(_dFdPHIi[2][ind] - _dFdPHIi[0][ind]))
                           - _K_SV[ind]*(Mob_GL*((_dFdPHIi[1][ind]+_dFdPHIi[2][ind])/2 - _dFdPHIi[0][ind])
                           + Mob_SV*(_dFdPHIi[2][ind] - _dFdPHIi[1][ind])) 
                           - _K_other[ind]*(Mob_GL*(_dFdPHIi[2][ind] - _dFdPHIi[0][ind])
                           + Mob_GL*(_dFdPHIi[2][ind] - _dFdPHIi[1][ind]));
        
     	/* C-H equation of motion */	    
        _dPHIdtCH[0][ind] = - _K_D[ind]*(_dFdPHIi[0][ind] - _dFdPHIi[1][ind])
                            - _K_D[ind]*(_dFdPHIi[0][ind] - _dFdPHIi[2][ind]);
                           
        _dPHIdtCH[1][ind] = - _K_D[ind]*(_dFdPHIi[1][ind] - _dFdPHIi[0][ind])
                            - _K_D[ind]*(_dFdPHIi[1][ind] - _dFdPHIi[2][ind]);
                           
        _dPHIdtCH[2][ind] = - _K_D[ind]*(_dFdPHIi[2][ind] - _dFdPHIi[0][ind])
                            - _K_D[ind]*(_dFdPHIi[2][ind] - _dFdPHIi[1][ind]);    
		
    }                        

    _F *= CUBE(gridsize) ;
    _Fraw *= CUBE(gridsize) ;
    
    /* calculate the maximum change in the phase feild per step */
    double Gmiddle = 0;
    if (dynamics_type != 1)
    {
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {
            if ( Gmiddle <= fmax(fmax(fabs(_dPHIdt0[0][ind]), fabs(_dPHIdt0[1][ind])),fabs(_dPHIdt0[2][ind])) )
                 Gmiddle = fmax(fmax(fabs(_dPHIdt0[0][ind]), fabs(_dPHIdt0[1][ind])),fabs(_dPHIdt0[2][ind]));
        }
    }
    else if (dynamics_type == 1)
    {
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {
            if ( Gmiddle <= fmax(fmax(fabs(_dPHIdtCH[0][ind]), fabs(_dPHIdtCH[1][ind])),fabs(_dPHIdtCH[2][ind])) )
                 Gmiddle = fmax(fmax(fabs(_dPHIdtCH[0][ind]), fabs(_dPHIdtCH[1][ind])),fabs(_dPHIdtCH[2][ind]));
        }
    }
    _G = Gmiddle;

    /* adjust time step at every step */
    timestep = dtmax/fmax(_G*dtmax/dphimax,1);
    
    vol_incr = vol_incr0/timestep;
    x_incr = x_incr0*liquid_volume/timestep;
    y_incr = y_incr0*liquid_volume/timestep;
    
    
    if (dynamics_type != 2 && dynamics_type != 3 && dynamics_type != 4 && dynamics_type != 5 && dynamics_type != 6 && dynamics_type != 7 && dynamics_type != 8) 
    {
        for(n=0;n<num_fields;n++)
        {
         if (dynamics_type == 1) /* Cahn-Hillard equation */
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

                _d2dPHIdt0[n][i*_NY*_NZ+j*_NZ+k] = 
                        (_dPHIdtCH[n][ip*_NY*_NZ+j*_NZ+k] + _dPHIdtCH[n][im*_NY*_NZ+j*_NZ+k]
                        +_dPHIdtCH[n][i*_NY*_NZ+jp*_NZ+k] + _dPHIdtCH[n][i*_NY*_NZ+jm*_NZ+k]
                        +_dPHIdtCH[n][i*_NY*_NZ+j*_NZ+kp] + _dPHIdtCH[n][i*_NY*_NZ+j*_NZ+km]
                        -_dPHIdtCH[n][i*_NY*_NZ+j*_NZ+k] * 6.0 ) / h2;
              }
            }
          }
        }
        
        if (dynamics_type == 0) /* no additional constraints */
        {
        /* 0: Ginzburg-Landau equation: not conserved */
            for (ind=0; ind<_NX*_NY*_NZ; ind++)
                _dPHIdt[n][ind] = _dPHIdt0[n][ind];
        }
        else if (dynamics_type == 1) /* no additional constraints */
        {
        /* 1: Cahn-Hilliard equation: conserved */
            for (ind=0; ind<_NX*_NY*_NZ; ind++)
                _dPHIdt[n][ind] = -_d2dPHIdt0[n][ind];
        }
      }
    }
    
    else if (dynamics_type == 2) 
    {
        /* Update _mprime */
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {
            _mprime[ind] = - _K_LS[ind]*(Mob_GL/2+Mob_LS)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
                           - _K_SV[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
                           - _K_LV[ind]*(Mob_GL/2+Mob_LV)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
                           - _K_other[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
        }       
        /* Update DM */	
        dmm = 0; dmu1 = 0;
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
                    dmm += _mprime[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
        for (ind=0; ind<_NX*_NY*_NZ; ind++)    
                    dmu1 -= ((5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0))) * _dPHIdt0[0][ind]) / dmm;
        
        /* Update gradient due to the change of m */
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {           
            _dPHIdt[0][ind] = _dPHIdt0[0][ind] + dmu1* _mprime[ind];
                    
            _dPHIdt[1][ind] = _dPHIdt0[1][ind] + dmu1
                            * (_K_LS[ind]*(-Mob_GL/2+Mob_LS)+_K_LV[ind]*Mob_GL+_K_SV[ind]*Mob_GL+_K_other[ind]*Mob_GL)
                            *(5.0-5.0*SQR(tanh(10.0*_PHI[0][ind]-5.0)));
                    
            _dPHIdt[2][ind] = _dPHIdt0[2][ind] + dmu1
                            * (_K_LS[ind]*Mob_GL+_K_LV[ind]*(-Mob_GL/2+Mob_LV)+_K_SV[ind]*Mob_GL+_K_other[ind]*Mob_GL)
                            *(5.0-5.0*SQR(tanh(10.0*_PHI[0][ind]-5.0)));                     
        }
        /* update chemical potential */
        _MU[0] = _MU[1] + dmu1;
       
        for (i=0;i<3;i++) {
            for (j=0;j<3;j++) {
                 _M[i][j] = _MU[i]-_MU[j];
            }
        }
    }
    
    else if (dynamics_type == 3)
    {       
        /* Update _mprime */
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {
            _mprime[ind] = - _K_LS[ind]*(Mob_GL/2+Mob_LS)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_SV[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_LV[ind]*(Mob_GL/2+Mob_LV)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_other[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
        }       
        /* Update DM */	
        dmm = 0; dmu1 = 0;
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
            dmm += _mprime[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
        for (ind=0; ind<_NX*_NY*_NZ; ind++)    
            dmu1 -= ((5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0))) * _dPHIdt0[0][ind]) / dmm;
        
        /* Update gradient due to the change of m */
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {           
            _dPHIdt[0][ind] = _dPHIdt0[0][ind] + dmu1* _mprime[ind];
            
            _dPHIdt[1][ind] = _dPHIdt0[1][ind] + dmu1
            * (_K_LS[ind]*(-Mob_GL/2+Mob_LS)+_K_LV[ind]*Mob_GL+_K_SV[ind]*Mob_GL+_K_other[ind]*Mob_GL)
            *(5.0-5.0*SQR(tanh(10.0*_PHI[0][ind]-5.0)));
            
            _dPHIdt[2][ind] = _dPHIdt0[2][ind] + dmu1
            * (_K_LS[ind]*Mob_GL+_K_LV[ind]*(-Mob_GL/2+Mob_LV)+_K_SV[ind]*Mob_GL+_K_other[ind]*Mob_GL)
            *(5.0-5.0*SQR(tanh(10.0*_PHI[0][ind]-5.0))); 
        }
        /* update chemical potential */
        _MU[0] = _MU[1] + dmu1;
        
        for (i=0;i<3;i++) {
            for (j=0;j<3;j++) {
                _M[i][j] = _MU[i]-_MU[j];
            }
        }
        
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
                            
			_d2dPHIdt0[n][i*_NY*_NZ+j*_NZ+k] = 
                 	   (_dPHIdtCH[n][ip*_NY*_NZ+j*_NZ+k] + _dPHIdtCH[n][im*_NY*_NZ+j*_NZ+k]
                   	 + _dPHIdtCH[n][i*_NY*_NZ+jp*_NZ+k] + _dPHIdtCH[n][i*_NY*_NZ+jm*_NZ+k]
                	    +_dPHIdtCH[n][i*_NY*_NZ+j*_NZ+kp] + _dPHIdtCH[n][i*_NY*_NZ+j*_NZ+km]
                   	 - _dPHIdtCH[n][i*_NY*_NZ+j*_NZ+k] * 6.0 ) / h2;
                        }
                    }
                }
            }
        
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {
            _dPHIdt[0][ind] = _dPHIdt[0][ind] - _d2dPHIdt0[0][ind];                      
            _dPHIdt[1][ind] = _dPHIdt[1][ind] - _d2dPHIdt0[1][ind];                      
            _dPHIdt[2][ind] = _dPHIdt[2][ind] - _d2dPHIdt0[2][ind];  
        }
    }
    
    else if (dynamics_type == 4) /* constrain all phase volume */
    {    
        _MU[1] = 0;
        /* Update _mprime */
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {
            _dGLdmuL[ind] = - _K_LS[ind]*(Mob_GL/2+Mob_LS)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_SV[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_LV[ind]*(Mob_GL/2+Mob_LV)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_other[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
           
            _dGLdmuV[ind] = - _K_LS[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_SV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_LV[ind]*(Mob_GL/2-Mob_LV)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_other[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)));
        
            _dGVdmuL[ind] = - _K_LS[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_SV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_LV[ind]*(Mob_GL/2-Mob_LV)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_other[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));       
             
            _dGVdmuV[ind] = - _K_LS[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_SV[ind]*(Mob_GL/2+Mob_SV)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_LV[ind]*(Mob_GL/2+Mob_LV)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_other[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)));

            _dGSdmuL[ind] = - _K_LS[ind]*(Mob_GL/2-Mob_LS)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_SV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_LV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_other[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));

            _dGSdmuV[ind] = - _K_LS[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_SV[ind]*(Mob_GL/2-Mob_SV)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_LV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_other[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)));
        }  
            /* initialize the parameters */
            A_11=0; A_12=0; A_21=0; A_22=0; B_1=0; B_2=0; dmuL=0; dmuV=0;
            
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {   
            A_11 += _dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
            A_12 += _dGLdmuV[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
            A_21 += _dGVdmuL[ind]*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)));
            A_22 += _dGVdmuV[ind]*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0))); 
            B_1  += -(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))*_dPHIdt0[0][ind];
            B_2  += -(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))*_dPHIdt0[2][ind];
        }
        
        B_2 += vol_incr;
        
        dmuL = (A_12*B_2 - A_22*B_1)/(A_12*A_21-A_11*A_22);
        dmuV = -(A_11*B_2 - A_21*B_1)/(A_12*A_21-A_11*A_22);
        
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {           
            _dPHIdt[0][ind] = _dPHIdt0[0][ind] + dmuL*_dGLdmuL[ind] + dmuV*_dGLdmuV[ind];
            
            _dPHIdt[1][ind] = _dPHIdt0[1][ind] + dmuL*_dGSdmuL[ind] + dmuV*_dGSdmuV[ind];
            
            _dPHIdt[2][ind] = _dPHIdt0[2][ind] + dmuL*_dGVdmuL[ind] + dmuV*_dGVdmuV[ind];

        }
    
        /* update chemical potential */
        _MU[0] = _MU[1] + dmuL;
        _MU[2] = _MU[1] + dmuV;
        
        for (i=0;i<3;i++) {
            for (j=0;j<3;j++) {
                _M[i][j] = _MU[i]-_MU[j];
            }
        }
    }
    
    else if (dynamics_type == 5) // constrain y-center of mass and all phase volumes
    {   
        /* Update _mprime */
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {
            _dGLdmuL[ind] = - _K_LS[ind]*(Mob_GL/2+Mob_LS)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_SV[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_LV[ind]*(Mob_GL/2+Mob_LV)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_other[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
            
            _dGLdmuV[ind] = - _K_LS[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_SV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_LV[ind]*(Mob_GL/2-Mob_LV)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_other[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)));
            
            _dGVdmuL[ind] = - _K_LS[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_SV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_LV[ind]*(Mob_GL/2-Mob_LV)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_other[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));       
            
            _dGVdmuV[ind] = - _K_LS[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_SV[ind]*(Mob_GL/2+Mob_SV)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_LV[ind]*(Mob_GL/2+Mob_LV)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_other[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)));
            
            _dGSdmuL[ind] = - _K_LS[ind]*(Mob_GL/2-Mob_LS)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_SV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_LV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_other[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
            
            _dGSdmuV[ind] = - _K_LS[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_SV[ind]*(Mob_GL/2-Mob_SV)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_LV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_other[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)));
        }  
        /* initialize the parameters */
        A_11=0; A_12=0; A_13=0; A_21=0; A_22=0; A_23=0; A_31=0; A_32=0; A_33=0;
        B_1=0; B_2=0; B_3=0; dmuL=0; dmuV=0; dfY=0; divider=0; _MU[1]=0;
        
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {   
            A_11 += _dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
            A_12 += _dGLdmuV[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
            A_13 += matrix_y[ind]*_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
            A_21 += _dGVdmuL[ind]*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)));
            A_22 += _dGVdmuV[ind]*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0))); 
            A_23 += matrix_y[ind]*_dGVdmuL[ind]*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)));
            A_31 += matrix_y[ind]*_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
            A_32 += matrix_y[ind]*_dGLdmuV[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
            A_33 += SQR(matrix_y[ind])*_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
            B_1 += -(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))*_dPHIdt0[0][ind];
            B_2 += -(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))*_dPHIdt0[2][ind];
            B_3 += -matrix_y[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))*_dPHIdt0[0][ind];
        }

        B_2 += vol_incr;
        B_3 += y_incr;
        
        divider = A_11*A_22*A_33-A_11*A_23*A_32-A_12*A_21*A_33+A_12*A_23*A_31+A_13*A_21*A_32-A_13*A_22*A_31;
        dmuL =  (A_12*A_23*B_3-A_12*A_33*B_2-A_13*A_22*B_3+A_13*A_32*B_2+A_22*A_33*B_1-A_23*A_32*B_1)/divider;
        dmuV = -(A_11*A_23*B_3-A_11*A_33*B_2-A_13*A_21*B_3+A_13*A_31*B_2+A_21*A_33*B_1-A_23*A_31*B_1)/divider;
        dfY  =  (A_11*A_22*B_3-A_11*A_32*B_2-A_12*A_21*B_3+A_12*A_31*B_2+A_21*A_32*B_1-A_22*A_31*B_1)/divider;


        
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {           
            _dPHIdt[0][ind] = _dPHIdt0[0][ind] + dmuL*_dGLdmuL[ind] + dmuV*_dGLdmuV[ind] + dfY*_dGLdmuL[ind]*matrix_y[ind];
            
            _dPHIdt[1][ind] = _dPHIdt0[1][ind] + dmuL*_dGSdmuL[ind] + dmuV*_dGSdmuV[ind] + dfY*_dGSdmuL[ind]*matrix_y[ind];
            
            _dPHIdt[2][ind] = _dPHIdt0[2][ind] + dmuL*_dGVdmuL[ind] + dmuV*_dGVdmuV[ind] + dfY*_dGVdmuL[ind]*matrix_y[ind];
            
        }
        
        /* update chemical potential */
        _MU[0] = _MU[1] + dmuL;
        _MU[2] = _MU[1] + dmuV;
        Fext_y = 0.0 + dfY;
        
        for (i=0;i<3;i++) {
            for (j=0;j<3;j++) {
                _M[i][j] = _MU[i]-_MU[j];
            }
        }
    }

    else if (dynamics_type == 6) // constrain x-center of mass and all phase volumes
    {   
        /* Update _mprime */
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {
            _dGLdmuL[ind] = - _K_LS[ind]*(Mob_GL/2+Mob_LS)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_SV[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_LV[ind]*(Mob_GL/2+Mob_LV)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_other[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
            
            _dGLdmuV[ind] = - _K_LS[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_SV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_LV[ind]*(Mob_GL/2-Mob_LV)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_other[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)));
            
            _dGVdmuL[ind] = - _K_LS[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_SV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_LV[ind]*(Mob_GL/2-Mob_LV)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_other[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));       
            
            _dGVdmuV[ind] = - _K_LS[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_SV[ind]*(Mob_GL/2+Mob_SV)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_LV[ind]*(Mob_GL/2+Mob_LV)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_other[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)));
            
            _dGSdmuL[ind] = - _K_LS[ind]*(Mob_GL/2-Mob_LS)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_SV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_LV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_other[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
            
            _dGSdmuV[ind] = - _K_LS[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_SV[ind]*(Mob_GL/2-Mob_SV)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_LV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_other[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)));
        }  
        /* initialize the parameters */
        A_11=0; A_12=0; A_13=0; A_21=0; A_22=0; A_23=0; A_31=0; A_32=0; A_33=0;
        B_1=0; B_2=0; B_3=0; dmuL=0; dmuV=0; dfX=0; divider=0; _MU[1]=0;
        
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {   
            A_11 += _dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
            A_12 += _dGLdmuV[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
            A_13 += matrix_x[ind]*_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
            A_21 += _dGVdmuL[ind]*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)));
            A_22 += _dGVdmuV[ind]*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0))); 
            A_23 += matrix_x[ind]*_dGVdmuL[ind]*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)));
            A_31 += matrix_x[ind]*_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
            A_32 += matrix_x[ind]*_dGLdmuV[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
            A_33 += SQR(matrix_x[ind])*_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
            B_1 += -(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))*_dPHIdt0[0][ind];
            B_2 += -(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))*_dPHIdt0[2][ind];
            B_3 += -matrix_x[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))*_dPHIdt0[0][ind];
        }

        B_2 += vol_incr;
        B_3 += x_incr;
        
        divider = A_11*A_22*A_33-A_11*A_23*A_32-A_12*A_21*A_33+A_12*A_23*A_31+A_13*A_21*A_32-A_13*A_22*A_31;
        dmuL =  (A_12*A_23*B_3-A_12*A_33*B_2-A_13*A_22*B_3+A_13*A_32*B_2+A_22*A_33*B_1-A_23*A_32*B_1)/divider;
        dmuV = -(A_11*A_23*B_3-A_11*A_33*B_2-A_13*A_21*B_3+A_13*A_31*B_2+A_21*A_33*B_1-A_23*A_31*B_1)/divider;
        dfX  =  (A_11*A_22*B_3-A_11*A_32*B_2-A_12*A_21*B_3+A_12*A_31*B_2+A_21*A_32*B_1-A_22*A_31*B_1)/divider;


        
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {           
            _dPHIdt[0][ind] = _dPHIdt0[0][ind] + dmuL*_dGLdmuL[ind] + dmuV*_dGLdmuV[ind] + dfX*_dGLdmuL[ind]*matrix_x[ind];
            
            _dPHIdt[1][ind] = _dPHIdt0[1][ind] + dmuL*_dGSdmuL[ind] + dmuV*_dGSdmuV[ind]+ dfX*_dGSdmuL[ind]*matrix_x[ind];
            
            _dPHIdt[2][ind] = _dPHIdt0[2][ind] + dmuL*_dGVdmuL[ind] + dmuV*_dGVdmuV[ind]+ dfX*_dGVdmuL[ind]*matrix_x[ind];
            
        }
        
        /* update chemical potential */
        _MU[0] = _MU[1] + dmuL;
        _MU[2] = _MU[1] + dmuV;
        Fext_x = 0.0 + dfX;
        
        for (i=0;i<3;i++) {
            for (j=0;j<3;j++) {
                _M[i][j] = _MU[i]-_MU[j];
            }
        }
    }
    
    else if (dynamics_type == 7) // constrain center of mass and all phase volumes 
    {   
        /* Update _mprime */
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {
            _dGLdmuL[ind] = - _K_LS[ind]*(Mob_GL/2+Mob_LS)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_SV[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_LV[ind]*(Mob_GL/2+Mob_LV)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_other[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
            
            _dGLdmuV[ind] = - _K_LS[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_SV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_LV[ind]*(Mob_GL/2-Mob_LV)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_other[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)));
            
            _dGVdmuL[ind] = - _K_LS[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_SV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_LV[ind]*(Mob_GL/2-Mob_LV)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_other[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));       
            
            _dGVdmuV[ind] = - _K_LS[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_SV[ind]*(Mob_GL/2+Mob_SV)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_LV[ind]*(Mob_GL/2+Mob_LV)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_other[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)));
            
            _dGSdmuL[ind] = - _K_LS[ind]*(Mob_GL/2-Mob_LS)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_SV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_LV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_other[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
            
            _dGSdmuV[ind] = - _K_LS[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_SV[ind]*(Mob_GL/2-Mob_SV)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_LV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_other[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)));
        }  
        /* initialize the parameters */
        A_11=0; A_12=0; A_13=0; A_14=0; A_21=0; A_22=0; A_23=0; A_24=0; A_31=0; A_32=0; A_33=0; A_34=0; A_41=0; A_42=0; A_43=0; A_44=0;
        B_1=0; B_2=0; B_3=0; B_4=0; dmuL=0; dmuV=0; dfX=0; dfY=0; divider=0; _MU[1]=0;
        
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {   
            A_11 += _dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
            A_12 += _dGLdmuV[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
            A_13 += matrix_x[ind]*_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
            A_14 += matrix_y[ind]*_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
            A_21 += _dGVdmuL[ind]*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)));
            A_22 += _dGVdmuV[ind]*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0))); 
            A_23 += matrix_x[ind]*_dGVdmuL[ind]*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)));
            A_24 += matrix_y[ind]*_dGVdmuL[ind]*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)));
            A_31 += matrix_x[ind]*_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
            A_32 += matrix_x[ind]*_dGLdmuV[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
            A_33 += SQR(matrix_x[ind])*_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
            A_34 += matrix_x[ind]*matrix_y[ind]*_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
            A_41 += matrix_y[ind]*_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
            A_42 += matrix_y[ind]*_dGLdmuV[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
            A_43 += matrix_x[ind]*matrix_y[ind]*_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
            A_44 += SQR(matrix_y[ind])*_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
            B_1 += -(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))*_dPHIdt0[0][ind];
            B_2 += -(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))*_dPHIdt0[2][ind];
            B_3 += -matrix_x[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))*_dPHIdt0[0][ind];
            B_4 += -matrix_y[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))*_dPHIdt0[0][ind];
        }
        
        B_2 += vol_incr;
        B_3 += x_incr;
        B_4 += y_incr;
        
        divider = A_11*A_22*A_33*A_44-A_11*A_22*A_34*A_43-A_11*A_23*A_32*A_44+A_11*A_23*A_34*A_42+A_11*A_24*A_32*A_43
                - A_11*A_24*A_33*A_42-A_12*A_21*A_33*A_44+A_12*A_21*A_34*A_43+A_12*A_23*A_31*A_44-A_12*A_23*A_34*A_41
                - A_12*A_24*A_31*A_43+A_12*A_24*A_33*A_41+A_13*A_21*A_32*A_44-A_13*A_21*A_34*A_42-A_13*A_22*A_31*A_44
                + A_13*A_22*A_34*A_41+A_13*A_24*A_31*A_42-A_13*A_24*A_32*A_41-A_14*A_21*A_32*A_43+A_14*A_21*A_33*A_42
                + A_14*A_22*A_31*A_43-A_14*A_22*A_33*A_41-A_14*A_23*A_31*A_42+A_14*A_23*A_32*A_41;
        dmuL = -(A_12*A_23*A_34*B_4-A_12*A_23*A_44*B_3-A_12*A_24*A_33*B_4+A_12*A_24*A_43*B_3
                +A_12*A_33*A_44*B_2-A_12*A_34*A_43*B_2-A_13*A_22*A_34*B_4+A_13*A_22*A_44*B_3
                +A_13*A_24*A_32*B_4-A_13*A_24*A_42*B_3-A_13*A_32*A_44*B_2+A_13*A_34*A_42*B_2
                +A_14*A_22*A_33*B_4-A_14*A_22*A_43*B_3-A_14*A_23*A_32*B_4+A_14*A_23*A_42*B_3
                +A_14*A_32*A_43*B_2-A_14*A_33*A_42*B_2-A_22*A_33*A_44*B_1+A_22*A_34*A_43*B_1
                +A_23*A_32*A_44*B_1-A_23*A_34*A_42*B_1-A_24*A_32*A_43*B_1+A_24*A_33*A_42*B_1)/divider;
        dmuV =  (A_11*A_23*A_34*B_4-A_11*A_23*A_44*B_3-A_11*A_24*A_33*B_4+A_11*A_24*A_43*B_3
                +A_11*A_33*A_44*B_2-A_11*A_34*A_43*B_2-A_13*A_21*A_34*B_4+A_13*A_21*A_44*B_3
                +A_13*A_24*A_31*B_4-A_13*A_24*A_41*B_3-A_13*A_31*A_44*B_2+A_13*A_34*A_41*B_2
                +A_14*A_21*A_33*B_4-A_14*A_21*A_43*B_3-A_14*A_23*A_31*B_4+A_14*A_23*A_41*B_3
                +A_14*A_31*A_43*B_2-A_14*A_33*A_41*B_2-A_21*A_33*A_44*B_1+A_21*A_34*A_43*B_1
                +A_23*A_31*A_44*B_1-A_23*A_34*A_41*B_1-A_24*A_31*A_43*B_1+A_24*A_33*A_41*B_1)/divider;
        dfX  = -(A_11*A_22*A_34*B_4-A_11*A_22*A_44*B_3-A_11*A_24*A_32*B_4+A_11*A_24*A_42*B_3
                +A_11*A_32*A_44*B_2-A_11*A_34*A_42*B_2-A_12*A_21*A_34*B_4+A_12*A_21*A_44*B_3
                +A_12*A_24*A_31*B_4-A_12*A_24*A_41*B_3-A_12*A_31*A_44*B_2+A_12*A_34*A_41*B_2
                +A_14*A_21*A_32*B_4-A_14*A_21*A_42*B_3-A_14*A_22*A_31*B_4+A_14*A_22*A_41*B_3
                +A_14*A_31*A_42*B_2-A_14*A_32*A_41*B_2-A_21*A_32*A_44*B_1+A_21*A_34*A_42*B_1
                +A_22*A_31*A_44*B_1-A_22*A_34*A_41*B_1-A_24*A_31*A_42*B_1+A_24*A_32*A_41*B_1)/divider;
        dfY  =  (A_11*A_22*A_33*B_4-A_11*A_22*A_43*B_3-A_11*A_23*A_32*B_4+A_11*A_23*A_42*B_3
                +A_11*A_32*A_43*B_2-A_11*A_33*A_42*B_2-A_12*A_21*A_33*B_4+A_12*A_21*A_43*B_3
                +A_12*A_23*A_31*B_4-A_12*A_23*A_41*B_3-A_12*A_31*A_43*B_2+A_12*A_33*A_41*B_2
                +A_13*A_21*A_32*B_4-A_13*A_21*A_42*B_3-A_13*A_22*A_31*B_4+A_13*A_22*A_41*B_3
                +A_13*A_31*A_42*B_2-A_13*A_32*A_41*B_2-A_21*A_32*A_43*B_1+A_21*A_33*A_42*B_1
                +A_22*A_31*A_43*B_1-A_22*A_33*A_41*B_1-A_23*A_31*A_42*B_1+A_23*A_32*A_41*B_1)/divider;
        
        
        
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {           
            _dPHIdt[0][ind] = _dPHIdt0[0][ind] + dmuL*_dGLdmuL[ind] + dmuV*_dGLdmuV[ind] 
                            + dfX*_dGLdmuL[ind]*matrix_x[ind] + dfY*_dGLdmuL[ind]*matrix_y[ind];
            
            _dPHIdt[1][ind] = _dPHIdt0[1][ind] + dmuL*_dGSdmuL[ind] + dmuV*_dGSdmuV[ind]
                            + dfX*_dGSdmuL[ind]*matrix_x[ind] + dfY*_dGSdmuL[ind]*matrix_y[ind];
            
            _dPHIdt[2][ind] = _dPHIdt0[2][ind] + dmuL*_dGVdmuL[ind] + dmuV*_dGVdmuV[ind]
                            + dfX*_dGVdmuL[ind]*matrix_x[ind] + dfY*_dGVdmuL[ind]*matrix_y[ind];
            
        }
        
        /* update chemical potential */
        _MU[0] = _MU[1] + dmuL;
        _MU[2] = _MU[1] + dmuV;
        Fext_x = 0.0 + dfX;
        Fext_y = 0.0 + dfY;
        
        for (i=0;i<3;i++) {
            for (j=0;j<3;j++) {
                _M[i][j] = _MU[i]-_MU[j];
            }
        }
    }
    
    else if (dynamics_type == 8)  // constrain center of mass and droplet volume
    {   
        /* Update _mprime */
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {
            _dGLdmuL[ind] = - _K_LS[ind]*(Mob_GL/2+Mob_LS)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_SV[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_LV[ind]*(Mob_GL/2+Mob_LV)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_other[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
            
            _dGLdmuV[ind] = - _K_LS[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_SV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_LV[ind]*(Mob_GL/2-Mob_LV)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_other[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)));
            
            _dGVdmuL[ind] = - _K_LS[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_SV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_LV[ind]*(Mob_GL/2-Mob_LV)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_other[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));       
            
            _dGVdmuV[ind] = - _K_LS[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_SV[ind]*(Mob_GL/2+Mob_SV)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_LV[ind]*(Mob_GL/2+Mob_LV)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_other[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)));
            
            _dGSdmuL[ind] = - _K_LS[ind]*(Mob_GL/2-Mob_LS)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_SV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_LV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))
            - _K_other[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
            
            _dGSdmuV[ind] = - _K_LS[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_SV[ind]*(Mob_GL/2-Mob_SV)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_LV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)))
            - _K_other[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*_PHI[2][ind]-5.0)));
        }  
        /* initialize the parameters */
        A_11=0; A_12=0; A_13=0; A_21=0; A_22=0; A_23=0; A_31=0; A_32=0; A_33=0;
        B_1=0; B_2=0; B_3=0; dmuL=0; dfX=0; dfY=0; divider=0; _MU[1]=0;
        
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {   
            A_11 += _dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
            A_12 += matrix_x[ind]*_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
            A_13 += matrix_y[ind]*_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));            
            A_21 += matrix_x[ind]*_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
            A_22 += SQR(matrix_x[ind])*_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
            A_23 += matrix_x[ind]*matrix_y[ind]*_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
            A_31 += matrix_y[ind]*_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
            A_32 += matrix_x[ind]*matrix_y[ind]*_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));
            A_33 += SQR(matrix_y[ind])*_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)));            
            B_1 += -(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))*_dPHIdt0[0][ind];
            B_2 += -matrix_x[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))*_dPHIdt0[0][ind];
            B_3 += -matrix_y[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))*_dPHIdt0[0][ind];
        }
        
        divider = A_11*A_22*A_33-A_11*A_23*A_32-A_12*A_21*A_33+A_12*A_23*A_31+A_13*A_21*A_32-A_13*A_22*A_31;
        dmuL =  (A_12*A_23*B_3-A_12*A_33*B_2-A_13*A_22*B_3+A_13*A_32*B_2+A_22*A_33*B_1-A_23*A_32*B_1)/divider;
        dfX  = -(A_11*A_23*B_3-A_11*A_33*B_2-A_13*A_21*B_3+A_13*A_31*B_2+A_21*A_33*B_1-A_23*A_31*B_1)/divider;
        dfY  =  (A_11*A_22*B_3-A_11*A_32*B_2-A_12*A_21*B_3+A_12*A_31*B_2+A_21*A_32*B_1-A_22*A_31*B_1)/divider;
        
        
        
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {           
            _dPHIdt[0][ind] = _dPHIdt0[0][ind] + dmuL*_dGLdmuL[ind] 
                            + dfX*_dGLdmuL[ind]*matrix_x[ind] + dfY*_dGLdmuL[ind]*matrix_y[ind];
            
            _dPHIdt[1][ind] = _dPHIdt0[1][ind] + dmuL*_dGSdmuL[ind] 
                            + dfX*_dGSdmuL[ind]*matrix_x[ind]+ dfY*_dGSdmuL[ind]*matrix_y[ind];
            
            _dPHIdt[2][ind] = _dPHIdt0[2][ind] + dmuL*_dGVdmuL[ind] 
                            + dfX*_dGVdmuL[ind]*matrix_x[ind]+ dfY*_dGVdmuL[ind]*matrix_y[ind];
            
        }
        
        /* update chemical potential */
        _MU[0] = _MU[1] + dmuL;
        Fext_x = 0.0 + dfX;
        Fext_y = 0.0 + dfY;
        
        for (i=0;i<3;i++) {
            for (j=0;j<3;j++) {
                _M[i][j] = _MU[i]-_MU[j];
            }
        }
    }    
    
    else
    {
        ERROR("unknown dynamics_type = "<<dynamics_type);
    }
       

      
#if 0
        INFO_Printf("_F = %20.12e\n", _F);
        INFO_Printf("_G = %20.12e\n", _G);
	INFO_Printf("A = %20.12e\n", A);
	INFO_Printf("B = %20.12e\n", B);
	INFO_Printf("C = %20.12e\n", C);
	INFO_Printf("D = %20.12e\n", D);
	INFO_Printf("E = %20.12e\n", E);
	INFO_Printf("FF = %20.12e\n", FF);
	INFO_Printf("A_h = %20.12e\n", A_h);
	INFO_Printf("B_h = %20.12e\n", B_h);
	INFO_Printf("C_h = %20.12e\n", C_h);
	INFO_Printf("E_h = %20.12e\n", E_h);
	INFO_Printf("A_hh = %20.12e\n", A_hh);
	INFO_Printf("dmuL = %20.12e\n", dmuL);
	INFO_Printf("dmuV = %20.12e\n", dmuV);
	INFO_Printf("dfy = %20.12e\n", dfY);
        INFO_Printf("divider = %20.12e\n", divider);
        write_array(matrix_y,"matrix_y.dat");		
//#else
        INFO_Printf("dmu1 = %20.12e _mprime[0] = %20.12e dmm = %20.12e\n", dmu1,_mprime[0],dmm);
        INFO_Printf("_dPHIdt0[0] = %20.12e _dPHIdt[0] = %20.12e \n", _dPHIdt0[0][0],_dPHIdt[0][0]);
        //INFO_Printf("_dPHIdx[0][0] = %20.12e \n", _dPHIdx[0][0]);
        write_array(_mprime,"_mprime.dat");
        write_array(_PHI[0],"PHI_0.dat");
        write_array(_PHI[1],"PHI_1.dat");
        write_array(_PHI[2],"PHI_2.dat");
        write_array(_d2PHI[0],"d2PHI_0.dat");
        write_array(_d2PHI[1],"d2PHI_1.dat");
        write_array(_d2PHI[2],"d2PHI_2.dat");
        write_array(_dPHIdt0[0],"dPHIdt0_0.dat");
        write_array(_dPHIdt0[1],"dPHIdt0_1.dat");
        write_array(_dPHIdt0[2],"dPHIdt0_2.dat");
        write_array(_dPHIdt[0],"dPHIdt_0.dat");
        write_array(_dPHIdx[0],"dPHIdx_0.dat"); 
        write_array(_dPHIdy[0],"dPHIdy_0.dat");
        write_array(_dPHIdz[0],"dPHIdz_0.dat");
        write_array(_dPHIdx[1],"dPHIdx_1.dat"); 
        write_array(_dPHIdy[1],"dPHIdy_1.dat");
        write_array(_dPHIdz[1],"dPHIdz_1.dat");
        write_array(_dPHIdx[2],"dPHIdx_2.dat"); 
        write_array(_dPHIdy[2],"dPHIdy_2.dat");
        write_array(_dPHIdz[2],"dPHIdz_2.dat");
        write_array(_K_S,"K_S.dat");
        write_array(_K_L,"K_L.dat");
        write_array(_K_LS,"K_LS.dat");
        write_array(_K_LV,"K_LV.dat");
        write_array(_K_SV,"K_SV.dat");
        write_array(_K_other,"K_other.dat");
        write_array(_dFdPHIi[0],"dFdPHIi_0.dat");  
        write_array(_dFdPHIi[1],"dFdPHIi_1.dat"); 
        write_array(_dFdPHIi[2],"dFdPHIi_2.dat"); 
        sleep(sleepseconds);
#endif
    

    return _F;        

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
#ifdef _STK_MPI
	if (myDomain == 0)
	{
        	finalcn.open(finalcnfile);
        	finalcn.write(this,zip,bg);
	}
	MPI_Barrier(MPI_COMM_WORLD);
#else
	finalcn.open(finalcnfile);
    finalcn.write(this,zip,bg);
#endif
	return 0;
}

void PhaseFieldFrame::saveintercn(int step)
{
	 if(savecn)
       	 if((step%savecnfreq)==0&&(step!=0))
	 	writeintercnfile();
}

int PhaseFieldFrame::writeintercnfile(int zip,bool bg)
{
#ifdef _STK_MPI
	if (myDomain == 0)
	{    
        printf("myDomain == %d/n",myDomain);
	    strcpy(intercnfile,"inter.cn");
        intercn.open(intercnfile);
	    intercn.write(this,zip,bg);
	}
	MPI_Barrier(MPI_COMM_WORLD);		
#else   
    intercn.open(intercnfile);
	intercn.write(this,zip,bg);

#endif
	return 0;
}

/*setup for writing prop.out file */
int PhaseFieldFrame::openpropfile()
{
#ifdef _STK_MPI
    if (myDomain == 0)
    {
	    return pf.open(outpropfile);
    }
    else
    {
	return 0;
    }
#else
    return pf.open(outpropfile);
#endif
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
    int n, i;   
    PhaseFieldFrame &d=*((PhaseFieldFrame *)p);
    f->Printf("%10d\n%10d\n%10d\n%10d\n",d.num_fields,d._NX,d._NY,d._NZ);   
        for (int n = 0; n< d.num_fields; n++)
        for (int i = 0; i < d.Ntot; i++)
            f->Printf("%20.12e\n", d._PHI[n][i]);
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
    
    //d.Alloc();
    
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
    sprintf(tmp,"%10d %20.12e %20.12e %20.12e\n",
            d.curstep,d._F,d._Fraw,d._G);
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
