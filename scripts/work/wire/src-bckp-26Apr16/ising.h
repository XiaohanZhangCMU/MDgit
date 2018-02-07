// Last Modified : Fri Sep  5 18:03:42 2008

#ifndef _ISING_H
#define _ISING_H

#ifdef _USETCL
#include "tcl.h"
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "general.h"
#include "organizer.h"
#include "display.h"

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

#define MAXSPECIES 10          /* maximum number of species */
#define MAXCOLORS  10          /* maximum number of colors in display */

/* Random number generator for platforms not supporting drand48 (e.g. cygwin) */
#ifdef _NODRAND48
#define drand48 (1.0/RAND_MAX)*rand
#endif

#define Realloc(p,t,s) {p=(t *)realloc(p,sizeof(t)*s);}


/* Input/output files */
class CNFile : public AUXFile /* configuration(cn) file */
{
public:
    CNFile(int t):AUXFile(t){};
    virtual char * describe();
private:
    virtual int writeblock(void *p);
    virtual int readblock(void *p);
};

class PropFile : public AUXFile/* Output property file */
{
public:
    PropFile(int t):AUXFile(t){};
    virtual char * describe();
    virtual int writeentry(void *p);
};

class Path
{
    enum { blocksize = 10000 };
public:
    int maxlen;
    int len;
    int _NX, _NY, _NZ, Ntot;
    int *_S0, S0tot;
    int *ix, *iy, *iz, *newS, *newStot;

    Path(int nc=0,int nx=100,int ny=100,int nz=1):
        maxlen(0),len(nc),_NX(nx),_NY(ny),_NZ(nz),
        _S0(0),ix(0),iy(0),iz(0),newS(0),newStot(0) {};
    
    ~Path(){ Free(); };

    void Alloc(int nc,int nx,int ny,int nz);
    void Free();

    void Append(int i,int j,int k,int s);
    void SetS0(int *s);

};

/* MD simulation framework */
class IsingFrame : public Organizer
{
public:
    int _NX, _NY, _NZ; /* number of spins in X, Y, Z direction */
                       /* _NZ = 1 for 2D problem */

    int *_S, *_S0;     /* spin array */
    int Stot;          /* sum of all spin */
    int Ntot;          /* total number of spins */

    double kBT, J, H;  /* parameters for the Hamiltonian */
    double E;          /* current energy (call Energy() to update */
    
    /* simulation setting */
    int totalsteps, curstep, continue_curstep, step0;

    double MC_spin, MC_accept, MC_accept_tot, MC_accept_ratio;
    int _RANDSEED;
    int _NSAMPLE, _NSUCCESS;  /* committor information for current spin configuration */

    /* paths */
    class Path pathA, pathB, pathC;
    
    /* Forward Flux Sampling (FFS) parameters */
    int N_lgst_cluster;/* total number of spins in the largest cluster */ 
    int N_lgst_avg;
    int N_cluster_temp;
    int *DFS_mark;     /* DFS algorithm marker */
    int *_Lam_array, *_Lam_check; /* lambda arrays for FFS */
    
    int     saveFFScn;  /* configuration save if saveFFScn=1 */
    int     FFSoption;	/* 0 for q/T, >1 for P(i+1/i) */
    int	    FFSpruning;	/* 0 for non, 1 for pruning */
    double  FFScn_weight; /* default = 1, pruning factor */
    double  FFS_Pp;	/* Pruning probability */
    int	    FFS0_check;	/* check for each FFS run */
    int	    lambda_A;	/* region A if lambda < lambda_A */
    int     lambda_B;   /* region B if lambda > lambda_B */ 
    int     FFScurstep; /* step marker for time-limited FFS */ 
    int     FFSautoend; /* autoend mark */ 

    /* Umbrella Sampling parameters */
    int YES_UMB; /* set to be 0 as a default, 1 if UMB is running */
    double UMB_K; /* k of the bias function 0.5k(n-n_c)^2 */
    int UMB_equilstep; /* equilibration step before UMB */
    int UMB_curstep; /* curstep check for UMB */
    int n_center; /* n_c of bias function */
    int delta_n; /* the width of the sampling */
    int *Narray; /* the array of number including all n in the window */
    int *Occurence; /* row occurence data */
    double *Probability; /* occurence times the weighting factor */

    int Kinetic; /* set to be 0 as a default, 1 if kinetic running */
    double Kinetic_Time;
    int Kinetic_Swip;
    double *time_array;
    int *n_array;
       
    /* Visualization */
    YWindow *win;                   /* Display window object */

    int win_width, win_height;      /* size of the window */
    double rotateangles[5];         /* rotation angles and scaling parameter */
    int plotfreq;                   /* frequency of re-plotting */
    double atomradius[MAXSPECIES];  /* size of atoms */
    double bondradius, bondlength;  /* thickness and length of bonds */

    char atomcolor[MAXSPECIES][30]; /* color of atom species */
    char bondcolor[30];             /* color of bonds */    
    char fixatomcolor[30];          /* color of fixed atoms, fixed[i]=1 */
    char highlightcolor[30];        /* color of highlighted atoms */
    char backgroundcolor[30];       /* background color of window */
    char colornames[MAXCOLORS][30]; /* names of allocated colors */
    unsigned colors[MAXCOLORS+15];  /* value of allocated colors */

    /* File input and output */
    double input[20000];                    /* general input variable in script file */
    class CNFile initcn,intercn,FFScn, finalcn;    /* configuration file objects */
    class PropFile pf;                      /* property file objects */
    char incnfile[200], finalcnfile[200];   /* configuration file name */
    char intercnfile[200];                  /* configuration file name */
    char FFScnfile[200];                    /* cn file for FFS */
    char outpropfile[200];                  /* property file name */
    char myname[200];
    char pathfile[200];                     /* path file */
    char command[1000]; int ncom;           /* shell command */
    int  savecn, savecnfreq;                /* frequency of saving cn files */
    int  saveprop, savepropfreq;            /* frequency of saving prop files */
    int  printfreq, calcryfreq, calcryfreq1;/* frequency of printing on scren */
    int  filecounter;                       /* number of current file */
    int  FFSfilecounter;
    int  zipfiles;
    

    /* Constructor (set initial values of member variables) */
    IsingFrame(): _NX(100), _NY(100), _NZ(1), _S(0),
               
                  /* simulation steps */
                  totalsteps(100), curstep(0), continue_curstep(0), step0(0),
                  
                  /* Monte Carlo simulation */
                  MC_spin(0),MC_accept(0),MC_accept_tot(0), MC_accept_ratio(0),
                  _RANDSEED(12345),
                  _NSAMPLE(100),_NSUCCESS(0),

                  /* Path sampling */
                  pathA(0,_NX,_NY,_NZ),
                  pathB(0,_NX,_NY,_NZ),
                  pathC(0,_NX,_NY,_NZ),
                  
		  /* FFS parameters */ 	
                   N_lgst_cluster(0), N_lgst_avg(0), N_cluster_temp(0),
                  DFS_mark(0), _Lam_array(0), _Lam_check(0), saveFFScn(0),
                  FFSoption(0), FFSpruning(0), FFScn_weight(1), FFS_Pp(0),
                  FFS0_check(1), lambda_A(0),lambda_B(0),FFScurstep(0),
                  FFSautoend(0),
 
                  /* Umbrella Sampling parameters */
                  YES_UMB(0), UMB_K(0), UMB_equilstep(0), UMB_curstep(0), n_center(0), delta_n(0), Narray(0), Occurence(0), 
        	  Probability(0), Kinetic(0), Kinetic_Time(10), Kinetic_Swip(0), time_array(0), n_array(0),
                  
                  /* Visualization */
                  win(0),win_width(350),win_height(350),
                  plotfreq(1),bondradius(0.1),bondlength(0),
                  
                  /* File input and output */
                  initcn(AUXFile::BLOCK),intercn(AUXFile::SERIES),
                  FFScn(AUXFile::SERIES), finalcn(AUXFile::BLOCK),pf(AUXFile::ENTRIES),
                  ncom(0),savecn(0),savecnfreq(100),
                  saveprop(0),savepropfreq(100),
                  printfreq(100), calcryfreq(100), calcryfreq1(1), filecounter(1),FFSfilecounter(1),zipfiles(0)
    { 
        /* Input/output control */
        input[0]=0;
        
    };

    virtual ~IsingFrame() {delete(win);}

    /* Parser */    
    void initparser();
    virtual int exec(const char *name);            
    void initvars();
    void runcommand();        /* run shell command */    

#ifdef _USETCL    
    virtual int Tcl_AppInit(Tcl_Interp *interp);
#endif

    void Alloc();
    void initspin();               /* initialize _S */
    void initpath(int pn);         /* initialize path C */
    void copypath(int pdes, int psrc); /* copy path psrc to pdes */
    void cutpath(int pn, int istart, int newlen);
    void reversepath(int pn);
    void mergepath(int phead,int ptail);
    
    /* Monte Carlo simulation */
    void MCrun();              /* run the simulation */
    double Energy();           /* compute energy of current spin config */
   
    /* Path sampling */
    void SampleMCpath(double smin, double smax, int nmax, int pn);
    void WalkonChain(double smin,double smax,int nmax);
    int ComputeSuccessRate(int nsample,double smin,double smax,int nmax);
    void AnalyzeConfigs(int nstart,int nend,int nsample,double smin, double smax,int nmax);
    void copyPathtoS(int pn, int nstep); /* copy a state along path to S */
    void copyStoS0();
    void copyS0toS();
    
    /* Forward Flux Sampling (FFS) implementation */
    int FFSprocess();
    void calcrystalorder();
    void allocDFS();
    void assign_Lam(); 
    void DFS(int);
  
    /* Umbrella Sampling implementation */
    void allocUMB();
    void printhist();
    void kineticprint();
    
    /* File input and output */
    virtual void saveintercn(int);   /* save intermediate cn files */
    int readcn();
    int setfilecounter();
    int setFFSfilecounter();
    int openintercnfile();
    int openFFScnfile();
    int openpropfile();
    int writefinalcnfile(int zip=1,bool bg=true);
    int writeintercnfile(int zip=1,bool bg=false);
    int writeFFScnfile(int zip=1,bool bg=false); 
    int writepath(int pn);                 /* write pathC into a file */
    int readpath(int pn);                  /* read  pathC from a file */
    
    /* Interface with Fortran code (Relax) */
    void writefortraninifile(char *fname);
    void fortranrelax();

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

#endif // _ISING_H

