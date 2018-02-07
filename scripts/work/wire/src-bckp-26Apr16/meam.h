/*
  meam.h
  by Wei Cai  caiwei@mit.edu
  Last Modified : Fri Feb 16 09:58:48 2007

  FUNCTION  : MD++ implementation of MEAM potential
*/

#ifndef _MEAM_H
#define _MEAM_H

#include "mdparallel.h"

class Matrix333
{
public:
    double element[3][3][3];
    void clear() { memset((double *)element,0,sizeof(double)*27);}
    Matrix333(){clear();}
};

class MEAMFrame : public MDPARALLELFrame
{
public:
    /* MEAM potential parameters */
    double zsmeam, alphas, betas[4], esubs, asubs, ts[4], rozros;
    double rcutmeam, cmin, cmax, repuls, attrac, legend;
    int ibarr, noscr;
    double res, alat, ielement, rcutmeam2;
    double sconst, scnres, xzbl, xzbl0, hmeam;
    char elt[10], lat[10];
    int enable_zbl_fdimer, enable_square_rscrn;
    
    /* MEAM potential storage space */
    double *atpe2b, *atpe3b, *rhotot, *embf;
    double *c8a, *dang1, *dang2, *cg8c ;

    Vector3   *tav, *rhsq, *a8b, *ag;
    Matrix33  *b8c;
    Matrix333 *d8d;

    double **scrnab;  char *scrnab_mem;
    double **dscrnab; char *dscrnab_mem;
    
    MEAMFrame() :zsmeam(0),alphas(0),esubs(0),asubs(0),
                 rozros(0),rcutmeam(0),cmin(2.0),cmax(2.8),
                 repuls(0),attrac(0),legend(0),ibarr(0),noscr(0),
                 res(0),alat(0),ielement(1),rcutmeam2(0),
                 sconst(1.3),scnres(1.0e-6),xzbl(-3),xzbl0(-1),hmeam(1.0e-5),
                 enable_zbl_fdimer(1),enable_square_rscrn(1),                 
                 atpe2b(0),atpe3b(0),rhotot(0),embf(0),
                 c8a(0),dang1(0),dang2(0),cg8c(0),
                 tav(0),rhsq(0),a8b(0),ag(0),b8c(0),d8d(0),
                 scrnab(0),scrnab_mem(0),dscrnab(0),dscrnab_mem(0)
                 
    {};

    virtual void potential ();
    void screen();
    void dscreen();
    void dscrfor();
    void rhoMEAM();
    void kraMEAM();
    int readMEAM();

    /* MEAM utility functions */
    inline double rhof(double r,double abc,double re,double rozero);
    inline double frhoi(double rhotp,double asub,double esub);
    inline double dfrhoi(double rho,double asub,double esub);
    inline double rscrn(double);
    inline double erose(double,double,double,double,double,double);
    inline double zbar(int,double,char *,double,double,double);

    double phiid(double rmagg, int i);
    double bar(double rho0,double A,int ibar,double z,double *dang1,double *dang2);
    double zbl(double,double charge1,double charge2);
    Matrix33 dcmij(Vector3 rr,double rs);
    double dscrn(double,double,double,double,double);

    virtual void Alloc();
    virtual void initvars();
    
    virtual int exec(const char *nam);
    virtual void initparser();
    
    void printpairpot();

#ifdef _PARALLEL
    void Broadcast_MEAM_Param();
#endif
};

#define drhof(abc,re,rho) (-(abc)*(rho)/(re))
#define xcut(x) SQR(1.0-(x)*(x)*(x)*(x))

#endif // _MEAM_H

