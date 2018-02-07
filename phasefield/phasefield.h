// Last Modified : Tue July 8 18:03:42 2014

#ifndef _PHASEFIELD_H
#define _PHASEFIELD_H

#include "mdparallel.h"

/* Physical constants */
#define EV   1.6021772E-19     /* (J)    */
#define BOLZ 1.380658E-23      /* (J/K)  */
#define KB   0.8617336E-4      /* (eV/K) */
#ifndef max
#define max(a,b) (((a)>(b))?(a):(b))
#endif
#ifndef min
#define min(a,b) (((a)<(b))?(a):(b))
#endif

#define MAXNUMFIELDS 3         /* maximum number of phase fields */
#define MAXCOLORS  10          /* maximum number of colors in display */

/* Random number generator for platforms not supporting drand48 (e.g. cygwin) */
#ifdef _NODRAND48
#define drand48 (1.0/RAND_MAX)*rand
#endif

#define Realloc(p,t,s) {p=(t *)realloc(p,sizeof(t)*s);}

#ifdef _STK_MPI
#include <StencilToolkit/Grid3D.hpp>
using namespace StencilToolkit;
#endif

#ifdef __USEFLOAT
#define REAL_T float
#define REAL_S FLOAT
#else
#define REAL_T double
#define REAL_S DOUBLE
#endif


/* Phase Field simulation framework */
class PhaseFieldFrame : public Organizer
{
public:
    int _NX, _NY, _NZ; /* number of spins in X, Y, Z direction */
                       /* _NZ = 1 for 2D problem */

    int Ntot;          /* Ntot = _NX*_NY*_NZ */
    int num_fields; 

    int model_type;    /* 0: traditional (single) phase field model, 1: multi phase field model (MPF) */
    
#ifdef _STK_MPI
    Index3D *_idx;
    Node3D *_node;
    int myDomain;
    int mySize;
    int data_size;
    int phase_size;
    int total_size;
#endif

#ifdef _USECUDA
    double *d_PHI, *d_dFdPHIi, *d_d2PHI, *d_dPHIdt0, *d_dPHIdt, *d_abs_dPHIdt, *d_dPHIdtCH;
    double *d_tmp_x, *d_tmp_y, *d_tmp_z, *d_EPS_loc, *d_dPHIdx, *d_dPHIdy, *d_dPHIdz;
    double *d_K_S, *d_K_L, *d_K_LS, *d_K_SV, *d_K_LV, *d_K_other, *d_K_D, *d_NW_orig;
    double *d_fden, *d_fden_raw, *d_matrix_x, *d_matrix_y, *d_mprime;
    double *d_U, *d_EPS, *d_EPS1, *d_EPS2, *d_EPS3, *d_R, *d_MU;
    double *d_UP, *d_EPSP, *d_EPS1P, *d_EPS2P, *d_EPS3P;
    double *d_dGLdmuL, *d_dGLdmuV, *d_dGVdmuL, *d_dGVdmuV, *d_dGSdmuL, *d_dGSdmuV, *d_Cprime0, *d_Cprime2, *d_CPHI;
#endif
    
#ifdef _STK_MPI
    
    Grid3D<REAL_T> *s_PHI[MAXNUMFIELDS]; /* phase fields */
    Grid3D<REAL_T> *s_dPHIdt[MAXNUMFIELDS], *s_dPHIdt0[MAXNUMFIELDS], *s_dPHIdtCH[MAXNUMFIELDS]; /* time derivatives */
    Grid3D<REAL_T> *s_d2dPHIdt0[MAXNUMFIELDS];
    Grid3D<REAL_T> *s_dPHIdx[MAXNUMFIELDS], *s_dPHIdy[MAXNUMFIELDS], *s_dPHIdz[MAXNUMFIELDS]; /* spatial derivatives */
    Grid3D<REAL_T> *s_d2PHI[MAXNUMFIELDS], *s_dFdPHI[MAXNUMFIELDS]; /* intermediate fields */
    Grid3D<REAL_T> *s_EPS_loc[MAXNUMFIELDS];
    Grid3D<REAL_T> *s_tmp_x[MAXNUMFIELDS], *s_tmp_y[MAXNUMFIELDS], *s_tmp_z[MAXNUMFIELDS]; 
    Grid3D<REAL_T> *s_d2dFdPHI[MAXNUMFIELDS]; /* intermediate fields */
    Grid3D<REAL_T> *s_mprime;
    Grid3D<REAL_T> *s_K_S, *s_K_L, *s_K_LS, *s_K_LV, *s_K_SV, *s_K_other; /* Phase region */
    Grid3D<REAL_T> *s_K_D;
    Grid3D<REAL_T> *s_matrix_x, *s_matrix_y;
    Grid3D<REAL_T> *s_NW_orig, *s_dGLdmuL, *s_dGLdmuV, *s_dGVdmuL, *s_dGVdmuV, *s_dGSdmuL, *s_dGSdmuV;

    REAL_T *_PHI[MAXNUMFIELDS]; /* phase fields */
    REAL_T *_dPHIdt[MAXNUMFIELDS], *_dPHIdt0[MAXNUMFIELDS], *_dPHIdtCH[MAXNUMFIELDS]; /* time derivatives */
    REAL_T *_d2dPHIdt0[MAXNUMFIELDS];
    REAL_T *_dPHIdx[MAXNUMFIELDS], *_dPHIdy[MAXNUMFIELDS], *_dPHIdz[MAXNUMFIELDS]; /* spatial derivatives */
    REAL_T *_d2PHI[MAXNUMFIELDS], *_dFdPHI[MAXNUMFIELDS], *_dFdPHIi[MAXNUMFIELDS]; /* intermediate fields */
    REAL_T *_EPS_loc[MAXNUMFIELDS];
    REAL_T *tmp_x[MAXNUMFIELDS], *tmp_y[MAXNUMFIELDS], *tmp_z[MAXNUMFIELDS]; 
    REAL_T *_d2dFdPHI[MAXNUMFIELDS]; /* intermediate fields */
    REAL_T *_mprime;
    REAL_T *_K_S, *_K_L, *_K_LS, *_K_LV, *_K_SV, *_K_other; /* Phase region */
    REAL_T *_K_D;
    REAL_T *matrix_x, *matrix_y;
    REAL_T *_NW_orig, *_dGLdmuL, *_dGLdmuV, *_dGVdmuL, *_dGVdmuV, *_dGSdmuL, *_dGSdmuV;
    REAL_T *sendbuf, *rbuf;
    int *sendbuf_ind, *rbuf_ind;
#else
#ifdef _USEFLOAT
    float *_PHI[MAXNUMFIELDS];     /* phase fields */
    float *_dPHIdt[MAXNUMFIELDS], *_dPHIdt0[MAXNUMFIELDS], *_dPHIdtCH[MAXNUMFIELDS]; /* time derivatives */
    float *_d2dPHIdt0[MAXNUMFIELDS]; 
    float *_dPHIdx[MAXNUMFIELDS], *_dPHIdy[MAXNUMFIELDS], *_dPHIdz[MAXNUMFIELDS];  /* spatial derivatives */
    float *_d2PHI[MAXNUMFIELDS], *_dFdPHI[MAXNUMFIELDS], *_dFdPHIi[MAXNUMFIELDS]; /* intermediate fields */
    float *_EPS_loc[MAXNUMFIELDS];
    float *tmp_x[MAXNUMFIELDS], *tmp_y[MAXNUMFIELDS], *tmp_z[MAXNUMFIELDS]; 
    float *_d2dFdPHI[MAXNUMFIELDS];                      /* intermediate fields */
    float *_mprime;
    float *_dGLdmuL, *_dGLdmuV, *_dGVdmuL, *_dGVdmuV, *_dGSdmuL, *_dGSdmuV;
    float *_K_S, *_K_L, *_K_LS, *_K_LV, *_K_SV, *_K_other;  /* Phase region */
    float *_K_D;
    float *matrix_x, *matrix_y;
    float *_NW_orig;
#else
    double *_PHI[MAXNUMFIELDS];     /* phase fields */
    double *_dPHIdt[MAXNUMFIELDS], *_dPHIdt0[MAXNUMFIELDS], *_dPHIdtCH[MAXNUMFIELDS]; /* time derivatives */
    double *_d2dPHIdt0[MAXNUMFIELDS]; 
    double *_dPHIdx[MAXNUMFIELDS], *_dPHIdy[MAXNUMFIELDS], *_dPHIdz[MAXNUMFIELDS];  /* spatial derivatives */
    double *_d2PHI[MAXNUMFIELDS], *_dFdPHI[MAXNUMFIELDS], *_dFdPHIi[MAXNUMFIELDS]; /* intermediate fields */
    double *_EPS_loc[MAXNUMFIELDS];
    double *tmp_x[MAXNUMFIELDS], *tmp_y[MAXNUMFIELDS], *tmp_z[MAXNUMFIELDS]; 
    double *_d2dFdPHI[MAXNUMFIELDS];                      /* intermediate fields */
    double *_mprime;
    double *_dGLdmuL, *_dGLdmuV, *_dGVdmuL, *_dGVdmuV, *_dGSdmuL, *_dGSdmuV;
    double *_K_S, *_K_L, *_K_LS, *_K_LV, *_K_SV, *_K_other;  /* Phase region */
    double *_K_D;
    double *matrix_x, *matrix_y;
    double *_NW_orig;
#endif
#endif

    double gridsize;
    double _EPS [MAXNUMFIELDS][MAXNUMFIELDS];  /* parameters for the Free Energy Functional */
    double _EPS1[MAXNUMFIELDS][MAXNUMFIELDS];  /* parameters for cubic anisotropic energy   */
    double _EPS2[MAXNUMFIELDS][MAXNUMFIELDS];  /* parameters for cubic anisotropic energy   */
    double _EPS3[MAXNUMFIELDS][MAXNUMFIELDS];  /* parameters for cubic anisotropic energy   */
    double _U   [MAXNUMFIELDS][MAXNUMFIELDS];  /* parameters for the Free Energy Functional */
    double _EPSP [MAXNUMFIELDS][MAXNUMFIELDS];  /* parameters for the Free Energy Functional */
    double _EPS1P[MAXNUMFIELDS][MAXNUMFIELDS];  /* parameters for cubic anisotropic energy   */
    double _EPS2P[MAXNUMFIELDS][MAXNUMFIELDS];  /* parameters for cubic anisotropic energy   */
    double _EPS3P[MAXNUMFIELDS][MAXNUMFIELDS];  /* parameters for cubic anisotropic energy   */
    double _UP   [MAXNUMFIELDS][MAXNUMFIELDS];  /* parameters for the Free Energy Functional */ 

   
    double _M   [MAXNUMFIELDS][MAXNUMFIELDS];
    //double _DM  [MAXNUMFIELDS][MAXNUMFIELDS]; /* disabled due to change of equation of motion */
    double _MU  [MAXNUMFIELDS]; 
    double _VOL [MAXNUMFIELDS]; 

    double _F;          /* current energy (call Energy() to update */
    double _E;
    double _Fraw;
    double _G;          /* maximum time derivative */
            
    int save_NW;
    double Mob_GL;
    double Mob_SV;
    double Mob_LS;
    double Mob_LV;
    double Mob_D;
    double Fext_x;
    double Fext_y;
    double Penalty_3phase;      /* factor for the three phase co-existing pentaly function */ 
    double Penalty_NW;  /* factor for the pentaly function accounts for the solid deviating from the original shape */
    double phase_volume[MAXNUMFIELDS];
    double liquid_volume, com_x, com_y, com_x0, com_y0;
    double com_y_incr, com_y_old;
    
    class Matrix33 _ROT_MATRIX;
   

    /* simulation setting */
    int totalsteps, elasteps, curstep, continue_curstep, step0;
    int dynamics_type;
    int mobility_type;
    int stop_flag_x;
    int stop_flag_y;
    double com_x_obj;
    double com_y_obj;
    int stop_flag_solid;
    double solid_obj;
    double L_vol_incr0;
    double vol_incr0;
    double x_incr0, y_incr0;
    int create_vapor, createfreq, shiftheight;
    int interfaceheight; 
   
    double timestep;
    double timestep_ela;
    double dtmax;
    double dphimax;
    int _RANDSEED;

    /* Visualization */
    YWindow *win;                   /* Display window object */

    int win_width, win_height;      /* size of the window */
    double rotateangles[5];         /* rotation angles and scaling parameter */
    int plotfreq;                   /* frequency of re-plotting */
    double atomradius[MAXNUMFIELDS];  /* size of phase field sample points */
    double bondradius, bondlength;  /* thickness and length of bonds */

    char atomcolor[MAXNUMFIELDS][30]; /* color of phase field sample points */
    char bondcolor[30];             /* color of bonds */    
    char fixatomcolor[30];          /* color of fixed atoms, fixed[i]=1 */
    char highlightcolor[30];        /* color of highlighted atoms */
    char backgroundcolor[30];       /* background color of window */
    char colornames[MAXCOLORS][30]; /* names of allocated colors */
    unsigned colors[MAXCOLORS+15];  /* value of allocated colors */
    double plot_threshold[MAXNUMFIELDS*2];

    /* File input and output */
    double input[20000];                    /* general input variable in script file */
    class CNFile initcn,intercn,finalcn;    /* configuration file objects */
    class PropFile pf;                      /* property file objects */
    char incnfile[200], finalcnfile[200];   /* configuration file name */
    char intercnfile[200];                  /* configuration file name */
    char outpropfile[200];                  /* property file name */
    char epseigfile[200];                   /* eigen strain file name */
    char epsfile[200];                      /* strain file name */
    char sigfile[200];                      /* stress file name */
    char dispfile[200];                     /* displacement file name */
    char Kfile[200];                        /* kinetic file name */
    char myname[200];
    char command[1000]; int ncom;           /* shell command */
    int  savecn, savecnfreq;                /* frequency of saving cn files */
    int  saveprop, savepropfreq;            /* frequency of saving prop files */
    int  printfreq;                         /* frequency of printing on scren */
    int  filecounter;                       /* number of current file */
    int  zipfiles;
   
#ifdef _ELASTIC 
    /* Elasticity */
    double _Cstiff[3][3][3][3];
    double *_udisp[3];
    double *_eps_eigen[3][3];
    double *_eps[3][3];
    double *_sig[3][3];
    double *_sigt[3][3];
    double *_dEdu[3];
    double *_dEdphi;
#endif

    double *_K;
    
    /* Constructor (set initial values of member variables) */
    PhaseFieldFrame(): _NX(100), _NY(100), _NZ(1), num_fields(1), model_type(0), gridsize(0), _F(0), _G(0), _Fraw(0),
                       save_NW(0), Mob_GL(1.0), Mob_SV(1.0), Mob_LS(1.0), Mob_LV(1.0), Mob_D(1.0), 
                       Fext_x(0.0), Fext_y(0.0), Penalty_3phase(0.5), Penalty_NW(0), stop_flag_x(0), stop_flag_y(0), com_x_obj(0), com_y_obj(0), stop_flag_solid(0), solid_obj(0), L_vol_incr0(0), vol_incr0(0), x_incr0(0), y_incr0(0), liquid_volume(0), com_x(0.0), com_y(0.0), com_x0(0.0), com_y0(0.0), com_y_incr(0.0), com_y_old(0.0), create_vapor(0), createfreq(100), shiftheight(0), interfaceheight(0),
               
                  /* simulation steps */
                  totalsteps(100), curstep(0), continue_curstep(0), step0(0), dynamics_type(0),
                  timestep(1e-6), dtmax(2.0), dphimax(0.001), _RANDSEED(12345),
                  
                  /* Visualization */
                  win(0),win_width(350),win_height(350),
                  plotfreq(1),bondradius(0.1),bondlength(0),
                  
                  /* File input and output */
                  initcn(AUXFile::BLOCK),intercn(AUXFile::SERIES),
                  finalcn(AUXFile::BLOCK),pf(AUXFile::ENTRIES),
                  ncom(0),savecn(0),savecnfreq(100),
                  saveprop(0),savepropfreq(100),
                  printfreq(100), filecounter(1),zipfiles(0)
    { 
        /* Input/output control */
        input[0]=0;
        
    };

    virtual ~PhaseFieldFrame() {delete(win);}

    /* Parser */    
    void initparser();
    virtual int exec(const char *name);            
    void initvars();
    void runcommand();        /* run shell command */    

#ifdef _USETCL    
    virtual int Tcl_AppInit(Tcl_Interp *interp);
#endif

    void Alloc();
    void initphasefield();        /* initialize _PHI */
    
    /* Phase Field simulation */
    void run();                    /* run the simulation */
    double FreeEnergy();           /* compute free energy of current phase field config */
    double FreeEnergy_single();    /* traditional (single) phase field model */
    double FreeEnergy_multi();     /* multi-phase field model (MPF) */
    double FreeEnergy_elastic();   /* phase field elasticity model */
  
#ifdef _USECUDA 
    void cuda_init();
    void cuda_memory_alloc();
    void cuda_integrator();
#ifdef _USEFLOAT
    float  cuda_FreeEnergy_single();    /* traditional (single) phase field model */
    float  cuda_FreeEnergy_multi();
#else
    double cuda_FreeEnergy_single();    /* traditional (single) phase field model */
    double cuda_FreeEnergy_multi();
#endif
#endif

#ifdef _OPENMP
    double FreeEnergy_single_omp();     /* traditional (single) phase field model */
    double FreeEnergy_multi_omp();     /* traditional (single) phase field model */
    void Phase_evolution();
#endif

#ifdef _STK_MPI
    double FreeEnergy_single_mpi();	/* traditional (single) phase field model */
    double FreeEnergy_multi_mpi();	/* multi-phase field model (MPF) */
#endif

    /* File input and output */
    virtual void saveintercn(int);   /* save intermediate cn files */
    int readcn();
    int setfilecounter();
    int openpropfile();
    int writefinalcnfile(int zip=1,bool bg=true);
    int writeintercnfile(int zip=1,bool bg=false);
    void write_array(double *,const char *);

    void setK();
#ifdef _ELASTIC
    void setdisp();
    void setinitialphi();
    void setepseig();

    int readepseig(const char *fname);
    int readdisp(const char *fname);
    int readK(const char *fname);
    void writeeps(char *fname);
    void writesig(char *fname);
    void writedisp(char *fname);
#endif
    
    /* Visualization */
    void winplot();        
    void winplot(int);
    void openwin();
    void plot();
    int  openwindow(int w,int h,const char *n);
    void closewindow();
    void wintogglepause();
    void alloccolors();
    void alloccolorsX();
    void rotate();
    void saverot();
    
};    

/* TCL driver program */
#ifdef _USETCL

int MD_Cmd(ClientData clientData,Tcl_Interp *interp,int argc, const char *argv[]);
int MD_SetVar(ClientData clientData,Tcl_Interp *interp,int argc, const char *argv[]);
int MD_GetVar(ClientData clientData,Tcl_Interp *interp,int argc, const char *argv[]);
int Tcl_Main_Wrapper(int argc, char *argv[]);

#define Tcl_Parser_Init(sim)                                         \
int Tcl_AppInit(Tcl_Interp *interp)                                  \
{                                                                    \
    return sim.Tcl_AppInit(interp);                                  \
}                                                                    \
                                                                     \
int MD_Cmd(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[])\
{                                                                    \
    FILE *istream, *ostream; int mypipe[2];                          \
                                                                     \
    if (pipe (mypipe))                                               \
    {                                                                \
        fprintf (stderr, "Pipe failed.\n");                          \
        return EXIT_FAILURE;                                         \
    }                                                                \
    istream = fdopen(mypipe[0], "r");                                \
    ostream = fdopen(mypipe[1], "w");                                \
    for(int i=1;i<argc;i++)                                          \
        fprintf(ostream,"%s ",(char *)argv[i]);                      \
    fprintf(ostream,"\n");                                           \
    fclose(ostream);                                                 \
    sim.parse_line(istream);                                         \
    fclose(istream);                                                 \
    return TCL_OK;                                                   \
}                                                                    \
                                                                     \
int MD_SetVar(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[])\
{                                                                    \
    int i, curn;                                                     \
    curn = sim.identify(argv[1]);                                    \
    if(curn<0)                                                       \
    {                                                                \
        WARNING("Unrecognized variable "<<argv[1]);                  \
        return TCL_OK;                                               \
    }                                                                \
    for(i=2;i<argc;i++)                                              \
    {                                                                \
        switch(sim.vartype[curn])                                    \
        {                                                            \
        case(INT): sscanf(argv[i],"%d",(int *)sim.varptr[curn]+i-2+sim.shift); break; \
        case(LONG): sscanf(argv[i],"%ld",(long *)sim.varptr[curn]+i-2+sim.shift); break;\
        case(DOUBLE): sscanf(argv[i],"%lf",(double *)sim.varptr[curn]+i-2+sim.shift); break;\
        case(STRING): strcpy((char *)sim.varptr[curn]+i-2+sim.shift,argv[i]); break;   \
        default: FATAL("unknown vartype ("<<sim.vartype[curn]<<")"); \
        }                                                            \
    }                                                                \
    return TCL_OK;                                                   \
}                                                                    \
                                                                     \
int MD_GetVar(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[])\
{                                                                    \
    int i, curn, offset;  char buf[1000];                            \
    curn = sim.identify(argv[1]);                                    \
    if(curn<0)                                                       \
    {                                                                \
        WARNING("Unrecognized variable "<<argv[1]);                  \
        return TCL_OK;                                               \
    }                                                                \
    for(i=1;i<argc;i++)                                              \
    {                                                                \
        if(i==1) offset=0;                                           \
        else sscanf(argv[i],"%d",&offset);                           \
        if(i>2) sprintf(buf," ");                                    \
        switch(sim.vartype[curn])                                    \
        {                                                            \
        case(INT): sprintf(buf,"%d",*((int *)sim.varptr[curn]+offset+sim.shift)); break; \
        case(LONG): sprintf(buf,"%ld",*((long *)sim.varptr[curn]+offset+sim.shift)); break;\
        case(DOUBLE): sprintf(buf,"%.16g",*((double *)sim.varptr[curn]+offset+sim.shift)); break;\
        case(STRING): sprintf(buf,"%s",(char *)sim.varptr[curn]+offset+sim.shift); break;\
        default: FATAL("unknown vartype ("<<sim.vartype[curn]<<")"); \
        }                                                            \
    }                                                                \
    Tcl_SetResult(interp, buf, TCL_VOLATILE);                        \
    return TCL_OK;                                                   \
}                                                                    

#endif

#endif // _PHASEFIELD_H

