#!/bin/bash

T=400
declare -i strain
for strain in 90 94 100 108 112 118
do 
 strain_fmt=`printf "%03d" $strain`
 echo "submit US for strain = ${strain_fmt}/1000  Temperature = ${T}"

 if (( $strain <= 100 ))
 then
  nwin=24
 else 
  if (( $strain <= 120 ))
  then
   nwin=20
  else 
   if (( $strain <= 180 ))
   then
    nwin=10
   else
    nwin=11
   fi
  fi
 fi

 cat > scripts/work/SiHomo-skin/PBS/status_7_umb/C7_auto_${strain_fmt}_${T}.pbs << FIN
#!/bin/bash
#PBS -N SiHomo_US_${strain_fmt}_${T}
#PBS -j oe
#PBS -l nodes=1:ppn=24,walltime=999:00:00
#PBS -V

cd /home/xzhang11/Planet/Libs/MD++UMB.svn3

declare -i f 
f=0
while (( \$f < $nwin ))
do
 bin/sw_mc2 scripts/work/SiHomo-skin/C7_run_US.tcl \$f $nwin ${strain_fmt} ${T} 50  &
 sleep 17
 let f+=1
done

wait
FIN

qsub scripts/work/SiHomo-skin/PBS/C7_auto_${strain_fmt}_${T}.pbs

done

