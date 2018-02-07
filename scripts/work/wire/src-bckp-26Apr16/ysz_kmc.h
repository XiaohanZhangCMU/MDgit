/*
  ysz.h
  by Wei Cai  caiwei@stanford.edu, Eunseok Lee  euniv@stanford.edu
  Last Modified : Mon Dec 28 15:43:31 2009

  FUNCTION  :  Easy-to-use KMC simulation Framework for YSZ

  This is a single CPU code.
*/


#ifndef _YSZ_H
#define _YSZ_H

#include "kmc.h"

#define LINKLEN    12
#define NUMJUMPDIR 6
#define JUMPUNIT   2

/* Input/output files */
class CorrFileShort : public AUXFile /* correlation fnction file */
{
public:
    CorrFileShort(int t):AUXFile(t){};
    virtual char * describe();
private:
    virtual int writeblock(void *p);
};

class CorrFileLong : public AUXFile /* correlation fnction file */
{
public:
    CorrFileLong(int t):AUXFile(t){};
    virtual char * describe();
private:
    virtual int writeblock(void *p);
};

class ACFile : public AUXFile /* response to AC external potential */
{
public:
    ACFile(int t):AUXFile(t){};
    virtual char * describe();
private:
    virtual int writeblock(void *p);
};

class DiffFile : public AUXFile /* response to AC external potential */
{
public:
    DiffFile(int t):AUXFile(t){};
    virtual char * describe();
private:
    virtual int writeblock(void *p);
};

class HistFile : public AUXFile /* response to AC external potential */
{
public:
    HistFile(int t):AUXFile(t){};
    virtual char * describe();
private:
    virtual int writeblock(void *p);
};

class AtomEyeFile : public AUXFile /* store atomeyefile */
{
public:
    AtomEyeFile(int t):AUXFile(t){};
    virtual char * describe();
private:
    virtual int writeblock(void *p);
};

class YSZFrame : public KMCFrame
{
public:
    /* number of smallest length units in X,Y,Z direction
     * in terms of this unit, Cations occupy even sites, Anions occupy odd sites
     */    
    int NCellX, NCellY, NCellZ; 
    int NCation, NAnion;    
    
    double epsilon_0, epsilon_r; /* permitivity */
    double trialfreq;            /* frequency prefactor for kmc */
    double VCoulomb;             /* prefractor to compute Coulomb potentials */

    /* public arrays and pointers */
    int *_RI;
    int *linklist, *anionlist, *cationlist, *vaclist, *vaclist0, *Ylist, *nnYlist, *snYlist;
    double *Ebarrier,*EbarrierN; /* look up table for every possible jump */
    double *CoulombE, *CatPot;   /* look up table for the entire volume */
    int *rx, *ry, *rz;           /* store the vacanies' real position */
    int *netcharge;              /* net charge on each plane */
    double *V_sc;                /* voltage on a plane */
    double Ebarrier_tab[64];  /* computed from ab initio */
    
    int total_cation, total_anion;
    int num_vac, num_jump, num_loop, num_Y;
    int count, count1, count2, count3;    

    char linklistfile[1000], ebarrierfile[1000], potfile[1000], Ydist[1000], catpotfile[1000];
    char catspeciesfile[1000], nnYlistfile[1000], snYlistfile[1000];
    
    /* external potential profile*/
    /* U = extU0 * sin(extw0*t)*/
    double extU0;             /* amplitude */
    double extw0;             /* frequency */
    double timelimit;      /* timelimit allowed for small external potential change */

    /* data */
    double *enposwext;
    double extdt;
    int ncycle, nextdt;

    /*AC response file */
    class ACFile af;
    char acfile[200];
    int saveac, saveacfreq;
    
    /* correlation file */
    class CorrFileShort cfs;       /* correlation file object */
    char corrfiles[200];
    class CorrFileLong cfl;       /* correlation file object */
    char corrfilel[200];
    int savecorr, savecorrfreq;

    /* diffusion coefficient file */
    class DiffFile df;
    char difffile[200];
    int savediff, savedifffreq;

    /* energy barrier histogram file & data */
    class HistFile hf;
    char histfile[200];
    int savehist, savehistfreq;
    double *histEb0, *histEb;
    
    /* atomeye file */
    class AtomEyeFile mf;
    char atomeyefile[200];
    int saveatomeye, saveatomeyefreq;
    
    /* correlation data */
    double enpos[3];
    double *corrs, *corrdxs, *corrdVs, *corrtimes;
    int iids,eids,ncorrsteps,ncorrbins,corrbinsizes;
    double maxcorrtimes, tdurs, corrdts;
    double *corrl, *corrdxl, *corrdVl, *corrtimel;
    int iidl,eidl,ncorrstepl,ncorrbinl,corrbinsizel;
    double maxcorrtimel, tdurl, corrdtl;
    
    /* interaction turn on/off*/
    int CTL_Ebarrier, CTL_CoulombE, CTL_CatPot;

    /* short range interaction constant*/
    double EnnY, EsnY;
    double Ebconst;

    /* save data for repeating */
    int savecontdata, savecontfreq;
    
    /* cluster exapnsion method */
    char cemfilesdir[1000];
    char clusterCC_file[1000];
    char clusterAA_file[1000];
    char clusterCA_file[1000];
    char clusterCCC_file[1000];
    char clusterCCA_file[1000];
    char clusterCAA_file[1000];
    char clusterAAA_file[1000];
    char eci_opt_file[1000];
    int *cCC, *cAA, *cCA, *cCCC, *cCCA, *cCAA, *cAAA, *cluster_spin, *cluster_phi, *n_cluster_phi;
    int *map_AA, *map_CC, *map_CA, *map_CCC, *map_CCA, *map_CAA, *map_AAA;
    int *map_An_clusters, *map_An_AA, *map_An_CA, *map_An_AAA, *map_An_CAA, *map_An_CCA;
    int nAA, nCC, nCA, nCCC, nCCA, nCAA, nAAA;
    double cluster_rcut, cluster_dcut;
    double *eci_opt;
    int neci;
    int *clusters_to_use;
    int nonzeroeci;
    int maxsize_cluster_buffer;

    /* avoid transient period */
    int ntransient;

    /* avoid too much energy change in Eb correction */
    double ebmaxratio;
    double ebaddratio;
    char ebcorrtable_file[1000];
    char ebxcorrtable_file[1000];
    double *ebcorrtable, *ebxcorrtable;
    int nebcorrtable;

    double EbZZ_Car, EbZY_Car, EbYY_Car;
    
    YSZFrame():NCellX(0),NCellY(0),NCellZ(0),
               epsilon_0(8.854e-12),epsilon_r(29), trialfreq(1e13), VCoulomb(1),
               Ebarrier(0),EbarrierN(0),CoulombE(0),CatPot(0),rx(0),ry(0),rz(0),
               netcharge(0),V_sc(0),
               total_cation(0),total_anion(0),
               num_vac(0),num_jump(0),num_loop(0),num_Y(0),
               count(0),count1(0),count2(0),count3(0),
               extU0(0), extw0(0), timelimit(0),
               enposwext(0), extdt(0), ncycle(0), nextdt(1),
               af(AUXFile::SERIES),saveac(0),saveacfreq(100),
               cfs(AUXFile::SERIES),cfl(AUXFile::SERIES), savecorr(0),savecorrfreq(100),
               df(AUXFile::SERIES),savediff(0),savedifffreq(100),
               hf(AUXFile::SERIES),savehist(0),savehistfreq(100),
               mf(AUXFile::SERIES),saveatomeye(0),saveatomeyefreq(100),
               corrs(0),corrdxs(0),corrdVs(0),corrtimes(0),
               iids(0),eids(0),ncorrsteps(1),ncorrbins(0),corrbinsizes(1),
               maxcorrtimes(0),tdurs(0), corrdts(0),
               corrl(0),corrdxl(0),corrdVl(0),corrtimel(0),
               iidl(0),eidl(0),ncorrstepl(1),ncorrbinl(0),corrbinsizel(1),
               maxcorrtimel(0),tdurl(0), corrdtl(0),
               CTL_Ebarrier(1), CTL_CoulombE(1), CTL_CatPot(1),
               EnnY(0), EsnY(0), Ebconst(0.0),
               savecontdata(1), savecontfreq(100),
               cCC(0), cAA(0), cCA(0), cCCC(0), cCCA(0), cCAA(0), cAAA(0),
               cluster_spin(0), cluster_phi(0), n_cluster_phi(0),
               map_AA(0), map_CC(0), map_CA(0), map_CCC(0), map_CCA(0), map_CAA(0), map_AAA(0),
               map_An_clusters(0), map_An_AA(0), map_An_CA(0), map_An_AAA(0), map_An_CAA(0), map_An_CCA(0),
               nAA(0), nCC(0), nCA(0), nCCC(0), nCCA(0), nCAA(0), nAAA(0),
               cluster_rcut(1.0), cluster_dcut(1.0),
               eci_opt(0), neci(0), clusters_to_use(0), nonzeroeci(0), maxsize_cluster_buffer(500000),
               ntransient(1),
               ebmaxratio(2.0),ebaddratio(0.085),
               ebcorrtable(0), ebxcorrtable(0), nebcorrtable(0),
               EbZZ_Car(0.58), EbZY_Car(1.29), EbYY_Car(1.86)
    {};

    virtual void initparser();
    virtual int  exec(const char *);
    virtual void initvars();
    virtual void Alloc();
    
    /* YSZ initialization functions */
    void AllocLinkList();
    void BuildLinkList();
    int ReadLinkList();
    void BuildEnergyBarrierList();
    int  ReadEnergyBarrierList();
    void BuildCoulombEList();
    void BuildEwald2DList();
    void BuildCatPotList();
    void BuildCatSpeciesList();
    void BuildNNYList();
    void BuildSNYList();
    void InitV_sc();
    int  ReadCoulombEList();
    int  ReadCatPotList();
    int  ReadCatSpeciesList();
    int  ReadNNYList();
    int  ReadSNYList();
    void CreateConfig_Vacancy_Only();
    
    void PrintLinkList();
    void PrintCationList();
    void PrintAnionList();

    void MigCellCationSites_unsorted(int,int,int [6]);
    void MigCellCationSites_sorted  (int,int,int [6]);
    void Add_Yttrium();
    void Add_Yttrium_Manually();
    void Add_Yttrium_Random();
    void Add_Yttrium_On_Surface();
    void Add_Yttrium_On_Surface_X();
    int ChecknY(int);
    void Add_Vacancy();
    void Add_Vacancy_On_Surface();
    
    void CreateDatabase();
    void FindMenergy();
    void kmcrun();
    void kmcruneq();
    void kmcrunAC();
    void kmcrun_continued();
    double CalInteract(int);
    double PDiffInteract(int,int,int,int);
    int Barrier(int,int,int);
    double mymod(int,int);
    double mymaxabs(int,int,int);
    void mysort(int,int *);
    
    int opencorrfiles();
    int opencorrfilel();
    int openacfile();
    int opendifffile();
    int openhistfile();
    int openatomeyefile();
    void CorrFcnX(double);
    void CorrFcnXn(int,int);
    void CorrFcnXV(double,double);

    void run_diffusion_1d();

    int set_corr_filecounter(int);
    int set_diff_filecounter(int);
    int set_hist_filecounter(int);
    int set_atomeye_filecounter(int);

    int ReadClusters();
    int obtain_AA(int,int,int,int,int,int);
    int obtain_CC(int,int,int,int,int,int);
    int obtain_CA(int,int,int,int,int,int);
    double cal_energy_by_CEM_initial();
    double cal_energy_by_CEM();
    double cal_energy_by_CEM_update(int,int);
    void mysort2(int *);
};

extern "C" void CoulombPotential(int lx, int ly, int lz, double *V);
extern "C" void Ewald2DPotential(int lx, int ly, int lz, double *V);

#endif // _YSZ_H

