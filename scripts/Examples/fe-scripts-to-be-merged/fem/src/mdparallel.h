/*
  mdparallel.h
  by Wei Cai  caiwei@stanford.edu
  Last Modified : Tue Aug  7 15:01:00 2007

  FUNCTION  :  Easy-to-use MD simulation Framefork

*/


#ifndef _MDPARALLEL_H
#define _MDPARALLEL_H

#include "md.h"

#ifndef _PARALLEL
class MDPARALLELFrame : public MDFrame { };
#else
#include "mpi.h"

enum {MSG_ATOM_LEN, MSG_ATOM, MSG_SKIN_LEN, MSG_SKIN, MSG_NEB_ATOM};

/* MD simulation framework */
class MDPARALLELFrame : public MDFrame
{
public:

    /* Domain decomposition */
    int myIX, myIY, myIZ;           /* IX, IY, IZ of current domain */
    int nXdoms, nYdoms, nZdoms;     /* nXdoms * nYdoms * nZdoms = total number of processes */
    int myDomain, numDomains;       /* domain ID of current process */
    double *domBoundX, *domBoundY, *domBoundZ; /* domain boundaries */
    double myXmin, myYmin, myZmin;  /* boundary of my domain */
    double myXmax, myYmax, myZmax;
    int *domainID, *globalID;       /* domain and global ID for each atom */
    double *_EPOT_IND_global;
    Vector3 *_F_global;
    Matrix33 *_VIRIAL_IND_global;

    int numNeighDoms, *neighDoms;
    MPI_Request *inRequests, *outRequests;
    MPI_Status  *inStatus,   *outStatus;

    /* Parallel Nudged Elastic Band Method */
    Vector3 *_Rc0, *_Rc1, *_Rc2, *_Fc0, *_Tan;
    double *_Ec;
    
    /* Constructor (set initial values of member variables) */
    MDPARALLELFrame(): myIX(0),myIY(0),myIZ(0),
                       nXdoms(1),nYdoms(1),nZdoms(1),
                       myDomain(0),numDomains(1),
                       domBoundX(0),domBoundY(0),domBoundZ(0),
                       myXmin(-10),myYmin(-10),myZmin(-10),
                       myXmax(+10),myYmax(+10),myZmax(+10),
                       domainID(0),globalID(0),
                       _EPOT_IND_global(0),_F_global(0),_VIRIAL_IND_global(0),
                       numNeighDoms(0),neighDoms(0),
                       inRequests(0),outRequests(0),inStatus(0),outStatus(0),
                       _Rc0(0),_Rc1(0),_Rc2(0),_Fc0(0),_Tan(0),_Ec(0)
    {
        MDFrame::MDFrame();
    }

    //virtual void call_potential() { potential_parallel(); }
    virtual void call_potential() { MDFrame::call_potential(); }
    
    virtual void initvars();
    virtual void initparser();
    virtual int exec(const char *name);
    virtual void Alloc();
    
    void ParallelInit(int *pargc, char ***pargv);
    void WaitForCommand();
    void Master_to_Slave(char *cmd);
    void Slave_to_Master_Atoms();
    void Slave_chdir();
    void Broadcast_H();
    void Broadcast_Atoms();
    void Master_Collect_Results();
    
    void Partition_Domains();
    void Comm_Neighbor_Domains_Atoms();
    void Mark_Local_Atoms();
    void potential_parallel();
    void eval_parallel();
    void run_parallel();
    void step_parallel();
    
    void alloc_all();
    void quit_all();

    /* Parallel Nudged Elastic Band Method */
    void nebrelax_parallel();
    void stringrelax_parallel();
    void stringrelax_parallel_1();
    void stringrelax_parallel_2();
    void stringrelax_parallel_3();
    void reparametrization_with_trimming();
    void AllocChain_parallel();
    void Broadcast_nebsetting();
    void initRchain_parallel();
    int  readRchain_parallel();
    int  writeRchain_parallel();
    void copyRchaintoCN_parallel();
    void copyCNtoRchain_parallel();
    void copyCNtoRchain_masteronly();
    int  readcnfile_parallel();
    int  writefinalcnfile_parallel(int, bool);
    void writeatomeyecfgfile_parallel();

};    
#endif //_PARALLEL

#endif //_MD_H
