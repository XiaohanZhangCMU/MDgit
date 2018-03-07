#!/bin/bash



sw_mc2 ../scripts/disl_nuc_hetero.tcl 1 0.120 001 1


sw_mc2 ../scripts/disl_nuc_hetero.tcl 1 0.110 001 1


sw_mc2 ../scripts/disl_nuc_hetero.tcl 1 0.100 001 1


sw_mc2 ../scripts/disl_nuc_hetero.tcl 1 0.090 001 1


sw_mc2 ../scripts/disl_nuc_hetero.tcl 1 0.080 001 1


sw_mc2 ../scripts/disl_nuc_hetero.tcl 1 0.070 001 1


sw_mc2 ../scripts/disl_nuc_hetero.tcl 1 0.060 001 1


#sleep 1
#mpirun -np $ncpu sw_mc2_mpich ../scripts/x_Si_disl/w_disl_nuc_hetero.tcl 4 0.012  001 1 

wait
#wipe $PBS_NODEFILE

