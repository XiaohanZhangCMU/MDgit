/*
  kmc.h
  by Wei Cai  caiwei@stanford.edu, Eunseok Lee  euniv@stanford.edu
  Last Modified : Wed Jan 20 16:00:06 2010

  FUNCTION  :  Easy-to-use KMC simulation Framefork

  Featuring :  1. Scripting input
               2. X-window display
               3. automatic simulation log
               4. convenient configuration and output file handle
               
  This is a single CPU code.
*/


#ifndef _KMC_H
#define _KMC_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef _USETCL
#include "tcl.h"
#endif

#include "general.h"
#include "filecls.h"
#include "organizer.h"
#include "display.h"
#include "linalg3.h"

/* Internal units
   Energy:  eV
   Length:  m

*/

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
class PropFile : public AUXFile /* Output property file */
{
public:
    PropFile(int t):AUXFile(t){};
    virtual char * describe();
    virtual int writeentry(void *p);
};

class CNFile : public AUXFile /* configuration(cn) file */
{
public:
    CNFile(int t):AUXFile(t){};
    virtual char * describe();
private:
    virtual int writeblock(void *p);
    virtual int readblock(void *p);
};


/* KMC simulation framework */
class KMCFrame : public Organizer
{
public:
    int _NP;                   /* total number of particles */

    class Vector3 *_R, *_SR;   /* real and scaled coordinates of atoms */
    class Vector3 *_R0;        /* place to save old coordinates of atoms */

    double Diffusivity;        /* diffusion coefficient to be computed */
    double treal;              /* real time at current step */

    int *fixed;                /* 1 if atom is fixed */
    int *species;              /* species type of atom */
    int *group;                /* group ID */
    int nspecies;              /* number of atom species */
    char element[MAXSPECIES][10]; /* species names */

    class Matrix33 _H, _H0;     /* simulation box vectors and copy */
    double _EPOT, *_EPOT_IND;   /* total potential energy and contribution per atom */
    double *_TOPOL;
    
    /* Kinetic Monte Carlo parameters */
    double _ATOMMASS[MAXSPECIES];/* in (g/mol) */
    double _OMEGA;               /* box volume (in A^3) */
    double _TIMESTEP;            /* in picosecond (ps) */
    double _T;                   /* desired and current temperature */

    /* Neighbor list */
    double _RLIST,_SKIN;        /* potential cut off radius (in A) */
    int _NIC;                   /* maximum number of atoms per cell */
    int _NNM;                   /* maximum number of neighbors per atom */
    int *nn;                    /* nn[i]: number of neighbors of atom i */
    int **nindex;               /* nindex[i][j]: the jth neighbor of atom i */
    int ****celllist;
    char *cell_mem, *nindex_mem;
    bool firsttime;
       
    /* Kinetic Monte Carlo simulation */
    int totalsteps, equilsteps, curstep;  /* simulation steps */
    int _RANDSEED;
    int select_atom, select_dir;
    
    /* Configuration manipulation */
    class Vector3 *storedr;         /* atomic displacements */
    class Matrix33 dH;              /* change of matrix H (cell size) */
    char  crystalstructure[30];     /* crystal structure for makecrystal() */
    double latticeconst[3];         /* lattice constant  for makecrystal() */
    double latticesize[3][4];       /* lattice size for makecrystal() */
    
    /* File input and output */
    double input[20000];                    /* general input variable in script file */
    char  output[10000];                    /* general output */
    char  output_fmt[10000];                /* output format */
    class CNFile initcn,intercn,finalcn;    /* configuration file objects */
    class PropFile pf;                      /* property file objects */
    char incnfile[200], finalcnfile[200];   /* configuration file name */
    char intercnfile[200];                  /* configuration file name */
    char outpropfile[200];                  /* property file name */
    char myname[200];                       /* name of the simulation */
    char potfile[100];                      /* name of potential file */    
    char command[1000]; int ncom;           /* shell command */
    int  savecn, savecnfreq;                /* frequency of saving cn files */
    int  saveprop, savepropfreq;            /* frequency of saving prop files */
    int  printfreq;                         /* frequency of printing on scren */
    int  filecounter;                       /* number of current file */
    int allocmultiple;                      /* allocate more memory than current number of atoms */

    /* Interface to AtomEye and other viewers */
    char atomeyepath[100], atomeyeexe[100];
    int atomeyerepeat[4];

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

    int plot_highlight_atoms[10000];/* indices of atoms to be highlighed */    
    double plot_limits[10];         /* set x,y,z, limits of plotting regime */
    int    plot_atom_info;          /* print information when atom is clicked */
    int    plot_map_pbc;            /* map all atoms in primary cell when plotting */
    double plot_color_windows[100]; /* only plot atoms whose properties fall into specified windows */
    double plot_color_bar[10];      /* how to color atoms according to their properties */
    int    plot_color_axis;          /* which atomic property specifies color */

    int autowritegiffreq;           /* frequency of outputing .GIF graphics file */
    double *color_ind;              /* property linked to atom color */

    /* Constructor (set initial values of member variables) */
    KMCFrame(): _NP(0),

                _R(0),_SR(0),_R0(0),treal(0.0),
                
               fixed(0),species(0),group(0),nspecies(1),
                _EPOT(0),_EPOT_IND(0), _TOPOL(0),

               /* Molecular Dynamics parameters */
               _OMEGA(0), _TIMESTEP(0), _T(0),
               
               /* Neighbor list */
               _RLIST(0),_SKIN(0),_NIC(100),_NNM(60),nn(0),nindex(0),
               celllist(0), cell_mem(0), nindex_mem(0), firsttime(true),

               /* KMC simulation steps */
               totalsteps(100), equilsteps(0), curstep(0), _RANDSEED(12345),
               select_atom(0), select_dir(0),
                
               /* File input and output */
               initcn(AUXFile::BLOCK),intercn(AUXFile::SERIES),
               finalcn(AUXFile::BLOCK),pf(AUXFile::ENTRIES),
               ncom(0),savecn(0),savecnfreq(100),saveprop(0),savepropfreq(100),
               printfreq(100),filecounter(1),allocmultiple(1),

               /* Visualization */
               win(0),win_width(350),win_height(350),
               plotfreq(1),bondradius(0.1),bondlength(0),plot_atom_info(1),plot_map_pbc(0),
               plot_color_axis(0),autowritegiffreq(0),color_ind(0)
        
    {
        /* Input/output control */
        input[0]=0;
        output[0]=0;
        output_fmt[0]=0;
        
        /* Configuration manipulation */
        latticeconst[0]=latticeconst[1]=latticeconst[2]=1.0;
        sprintf(crystalstructure,"simple-cubic");

        /* Plot settings */
        plot_highlight_atoms[0]=0;
        plot_limits[0]=0;
        plot_color_windows[0]=0;
        plot_color_bar[0]=0;

        
    };

    virtual ~KMCFrame() { delete(win);}
            
    /* Parser */    
    virtual void initparser();
    virtual int exec(const char *name);            
    virtual void initvars();
    void runcommand();        /* run shell command */
#ifdef _USETCL    
    virtual int Tcl_AppInit(Tcl_Interp *interp);
#endif
    
    /* Coordinate transformation */
    void SHtoR();
    void RHtoS();
    void RtoR0();
    void R0toR();
    bool Bond(int I, int J) const;
    
    virtual void Alloc();

    /* Potentials */
    virtual void potential();
    virtual void ext_potential();
    virtual double ext_potential(int iatom);
    virtual void int_potential();
    virtual double int_potential(int iatom);
    virtual double int_potential(int iatom, int jatom);

    /* Neighbor list */
    void NbrList_reconstruct();  /* use cell list combined with Verlist */
    void NbrList_print();        /* for debug use */
    bool NbrList_needrefresh();
    void NbrList_init(int, int);
    void NbrList_initcell(int, int, int, int);
    void NbrList_free();
    void NbrList_freecell();
    void NbrList_refresh();
#define refreshneighborlist NbrList_refresh
            
    /* KMC simulation */
    virtual int kmc_step();

    virtual void runkmc();           /* run the simulation */
    virtual void calcprop();         /* calculate physical properties */
    virtual void calcoutput();       /* calculate output string */

    void randomposition();

    void eval();
    void multieval();                /* for testing/timing purposes */

    /* Configuration manipulation */
    void makecrystal();         /* create perfect crystal structure */
    
    void scaleH();              /* hydrostatic scale of simulation box */
    void setH();
    void saveH();
    void restoreH();            /* restore H from H0 */
    
    void shiftbox();            /* shift PBC box */
    void redefinepbc();
    void extendbox();           /* make multiple copies of the box*/

    void moveatom();
    void movegroup();
    void pbcshiftatom();        /* shift atoms in PBC box */
    
    void fixatoms_by_ID();      /* fix a set of atoms whose ID are specified */
    void fixatoms_by_position();/* fix all atoms whose position falls within a specified regime */
    void fixallatoms();
    void freeallatoms();
    void reversefixedatoms();

    void setfixedatomsspecies();
    void setfixedatomsgroup();
    void reversespecies();      /* species: 0->1 , 1->0 */
    void movefixedatoms();      /* move atoms that are marked as fixed */
    void removefixedatoms();    /* remove atoms that are marked as fixed */
    void markremovefixedatoms();/* mark fixed atoms for removal (fix=-1) */
    void removeellipsoid();     /* remove all atoms within an ellipsoid */    
    void removerectbox();       /* remove all atoms within a parallelpepid */
    
    /* File input and output */
    virtual void saveintercn(int);   /* save intermediate cn files */
    int readcn();
    int setfilecounter();
    int openintercnfile();
    int openpropfile();
        
    int writefinalcnfile(int zip=1,bool bg=true);
    int writeintercnfile(int zip=1,bool bg=false);

    /* Interface with AtomEye and other viewers */
    void writeatomeyecfgfile(char *fname);
    void convertCNtoCFG();
    void atomeye();
    
    /* Print out energies of current configuration */
    void writeENERGY(char *fname);
    void writePOSITION(char *fname);
    void GnuPlotHistogram();

    /* Visualization */
    virtual void winplot();        
    virtual void winplot(int);
    virtual void openwin();
    virtual void plot();
    int     openwindow(int w,int h,const char *n);
    void    closewindow();
    void    wintogglepause();
    void    alloccolors();
    void    alloccolorsX();
    void    rotate();
    void    saverot();

    /* set variables from arg list */
    void set_dirname(char *);
    void set_randseed(char *);
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

#endif //_KMC_H
