#!/bin/bash

nwin=24

cat > scripts/LLNL2010/PBS/ising.pbs << FIN
#!/bin/bash
#PBS -N PNAS_US_${strain_fmt}_${T}
#PBS -j oe
#PBS -l nodes=1:ppn=24,walltime=72:00:00
#PBS -V
cd /home/xzhang11/Planet/Libs/MD++.svn

declare -i f 
f=0
while (( \$f < $nwin ))
do
 bin/ising_mc2 scripts/LLNL2010/Day1_ising_US.tcl \$f 6 8 0.06 1.7 &
 sleep 17
 let f+=1
done

wait
FIN

qsub scripts/LLNL2010/PBS/ising.pbs

