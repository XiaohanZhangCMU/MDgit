// Last Modified : Fri Sep  5 18:03:42 2008

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

#define MAXNUMFIELDS 60         /* maximum number of phase fields */
#define MAXCOLORS  10          /* maximum number of colors in display */

/* Random number generator for platforms not supporting drand48 (e.g. cygwin) */
#ifdef _NODRAND48
#define drand48 (1.0/RAND_MAX)*rand
#endif

#define Realloc(p,t,s) {p=(t *)realloc(p,sizeof(t)*s);}


/* Phase Field simulation framework */
class PhaseFieldFrame : public Organizer
{
public:
    int _NX, _NY, _NZ; /* number of spins in X, Y, Z direction */
                       /* _NZ = 1 for 2D problem */

    int Ntot;          /* Ntot = _NX*_NY*_NZ */
    int num_fields; 

#ifdef _USEFLOAT
    float *_PHI[MAXNUMFIELDS];     /* phase fields */
    float *_dPHIdt[MAXNUMFIELDS]; /* time derivatives */
    float *_dPHIdx[MAXNUMFIELDS], *_dPHIdy[MAXNUMFIELDS], *_dPHIdz[MAXNUMFIELDS];  /* spatial derivatives */
    float *_d2PHI[MAXNUMFIELDS], *_dFdPHIi[MAXNUMFIELDS]; /* intermediate fields */
#else
    double *_PHI[MAXNUMFIELDS];     /* phase fields */
    double *_dPHIdt[MAXNUMFIELDS], *_dPHIdtCH[MAXNUMFIELDS]; /* time derivatives */
    double *_dPHIdx[MAXNUMFIELDS], *_dPHIdy[MAXNUMFIELDS], *_dPHIdz[MAXNUMFIELDS];  /* spatial derivatives */
    double *_d2PHI[MAXNUMFIELDS], *_dFdPHIi[MAXNUMFIELDS]; /* intermediate fields */
#endif

    double gridsize;

    
    double _F;          /* current energy (call Energy() to update */
    double _G;          /* maximum time derivative */
    
    double alpha;
    double beta;
    double gamma;
    double kappa;
    double L;    

    class Matrix33 _ROT_MATRIX;
   

    /* simulation setting */
    int totalsteps, curstep, continue_curstep, step0;

    double timestep;
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
    char myname[200];
    char command[1000]; int ncom;           /* shell command */
    int  savecn, savecnfreq;                /* frequency of saving cn files */
    int  saveprop, savepropfreq;            /* frequency of saving prop files */
    int  printfreq;                         /* frequency of printing on scren */
    int  filecounter;                       /* number of current file */
    int  zipfiles;
    

    /* Constructor (set initial values of member variables) */
    PhaseFieldFrame(): _NX(100), _NY(100), _NZ(1), num_fields(1), gridsize(0), _F(0), _G(0),
                       alpha(1.0), beta(1.0), gamma(1.0), kappa(2.0), L(1.0),
               
                  /* simulation steps */
                  totalsteps(100), curstep(0), continue_curstep(0), step0(0),
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
  
#ifdef _USECUDA 
#ifdef _USEFLOAT
    float  cuda_FreeEnergy_single();    /* traditional (single) phase field model */
#else
    double cuda_FreeEnergy_single();    /* traditional (single) phase field model */
#endif
#endif

#ifdef _OPENMP
    double FreeEnergy_single_omp();     /* traditional (single) phase field model */
    double FreeEnergy_multi_omp();     /* traditional (single) phase field model */
#endif

    /* File input and output */
    virtual void saveintercn(int);   /* save intermediate cn files */
    int readcn();
    int setfilecounter();
    int openintercnfile();
    int openpropfile();
    int writefinalcnfile(int zip=1,bool bg=true);
    int writeintercnfile(int zip=1,bool bg=false);
    void write_array(double *,const char *);
    
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

