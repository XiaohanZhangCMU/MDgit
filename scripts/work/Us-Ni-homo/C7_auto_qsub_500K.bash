#!/bin/bash

T=500
declare -i strain
for strain in 90 92 95 100 105 110 115 
do 
 strain_fmt=`printf "%03d" $strain`
 echo "submit US for strain = ${strain_fmt}/1000  Temperature = ${T}"

 if (( $strain <= 110 ))
 then
  nwin=24
 else 
  if (( $strain <= 120 ))
  then
   nwin=18
  else 
   if (( $strain <= 130 ))
   then
    nwin=16
   else
    nwin=14
   fi
  fi
 fi

 cat > scripts/work/PNAShomo/PBS/C7_${strain_fmt}_${T}.pbs << FIN
#!/bin/bash
#PBS -N PNAS_US_${strain_fmt}_${T}
#PBS -j oe
#PBS -l nodes=1:ppn=24,walltime=72:00:00
#PBS -V
cd $PBS_O_WORKDIR
echo $PWD

declare -i f 
f=0
while (( \$f < $nwin ))
do
 bin1/eam_mc2 scripts/work/PNAShomo/C7_run_US.tcl \$f $nwin ${strain_fmt} ${T} &
 sleep 17
 let f+=1
done

wait
FIN

# qsub scripts/work/PNAShomo/PBS/C7_${strain_fmt}_$T.pbs

done

