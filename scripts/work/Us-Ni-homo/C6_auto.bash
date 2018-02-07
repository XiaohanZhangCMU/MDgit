#!/bin/bash

for strain in 090 092 095 100 105 110 115 120 125 130 135 140 145
 do echo "bin/eam_mc2 scripts/work/PNAShomo/C6_prepare_US.tcl  $strain  0"
 #do bin/eam_mc2 scripts/work/PNAShomo/C6_prepare_US.tcl  $strain  0
done

for strain in 090 092 095 100 105 110 115 120 125 130 
 do echo "bin/eam_mc2 scripts/work/PNAShomo/C6_prepare_US.tcl  $strain  100"
 #do bin/eam_mc2 scripts/work/PNAShomo/C6_prepare_US.tcl  $strain  100
done

for strain in 090 092 095 100 105 110 115 120 125 
 do echo "bin/eam_mc2 scripts/work/PNAShomo/C6_prepare_US.tcl  $strain  200"
 #do bin/eam_mc2 scripts/work/PNAShomo/C6_prepare_US.tcl  $strain  200
done

for strain in 090 092 095 100 105 110 115 120 125
 do echo "bin/eam_mc2 scripts/work/PNAShomo/C6_prepare_US.tcl  $strain  300"
 #do bin/eam_mc2 scripts/work/PNAShomo/C6_prepare_US.tcl  $strain  300
done

for strain in 090 092 095 100 105 110 115 120 
 do echo "bin/eam_mc2 scripts/work/PNAShomo/C6_prepare_US.tcl  $strain  400"
 #do bin/eam_mc2 scripts/work/PNAShomo/C6_prepare_US.tcl  $strain  400
done

for strain in 090 092 095 100 105 110 115 
 do echo "bin/eam_mc2 scripts/work/PNAShomo/C6_prepare_US.tcl  $strain  500"
 #do bin/eam_mc2 scripts/work/PNAShomo/C6_prepare_US.tcl  $strain  500
done


