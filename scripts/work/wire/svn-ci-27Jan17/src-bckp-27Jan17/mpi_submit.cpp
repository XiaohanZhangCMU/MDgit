/* 
   Utility to launch a series of serial jobs to the queue of a parallel
   cluster as one package.

   As a test case, compile your codes by

     make sw build=R SYS=wcr
     make mpi_submit build=R SYS=wcr_mpich

   You can then submit your PBS file by

     qsub scripts/work/mpi_submit/testsubmit.pbs

   Your PBS file may look like the following

--------------------------------------------------------
#!/bin/bash
#PBS -N mpi_submit
#PBS -j oe
#PBS -l nodes=2:ppn=8,walltime=1:00:00
#PBS -V

ncpu=`cat $PBS_NODEFILE | wc -w`
mpiexec -np $ncpu bin/mpi_submit_wcr_mpich bin/sw_wcr scripts/work/mpi_submit/si-test.tcl 0
--------------------------------------------------------

This is equivalent to submit 16 serial MD++ jobs:


bin/sw_wcr scripts/work/mpi_submit/si-test.tcl 0 0
bin/sw_wcr scripts/work/mpi_submit/si-test.tcl 0 1
bin/sw_wcr scripts/work/mpi_submit/si-test.tcl 0 2
...

bin/sw_wcr scripts/work/mpi_submit/si-test.tcl 0 15


The two arguments follwing the tcl file name can be used to
inform tcl to do different things in each job.

*/

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

void submit_jobs(int *pargc, char ***pargv)
{
    int myCPU, nCPUs, myArg, i;
    char command[100];

    MPI_Init(pargc,pargv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myCPU);
    MPI_Comm_size(MPI_COMM_WORLD, &nCPUs);

    printf("[%d] argc = %d\n", myCPU, *pargc);

    if (*pargc == 3)
    {
       for(i=1;i<=2;i++)
          printf("[%d] myArg[%d] = %s\n",myCPU, i, (*pargv)[i]);
       sprintf(command,"%s %s %d", (*pargv)[1], (*pargv)[2], myCPU);
       printf("[%d] command = %s\n", myCPU, command);
    }
    else if (*pargc == 4)
    {
       for(i=1;i<=3;i++) 
          printf("[%d] myArg[%d] = %s\n",myCPU, i, (*pargv)[i]);
       sprintf(command,"%s %s %s %d", (*pargv)[1], (*pargv)[2], (*pargv)[3], myCPU);
       printf("[%d] command = %s\n", myCPU, command);
    }
    else
    {
       printf("number of arguments need to be either 2 or 3!");
       sprintf(command," ");
    }

    system(command);    
}

int main(int argc, char *argv[])
{
    submit_jobs(&argc, &argv);

    MPI_Finalize();
    return 0;
}

