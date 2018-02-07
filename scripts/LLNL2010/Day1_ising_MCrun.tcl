# Ising 2D Monte Carlo simulation
#-------------------------------------------

set graphic "On" 
source "scripts/Examples/Tcl/startup.tcl"

# default value is h=0.06 and kBT=1.7
if { $argc==0} {
        set h 0.06
        set t 1.7 
        set totsteps 10000
} elseif { $argc > 0 } {
        set h [ lindex $argv 0 ]
	set t [ lindex $argv 1 ]
        set totsteps [ lindex $argv 2 ]
}

# make a directory "runs/2D_ISING" if it does not exist
set check [ glob -nocomplain runs/2D_ISING ]
set key [ llength $check ]
if { $key < 1 } {
exec mkdir runs/2D_ISING
}

#*******************************************
# Definition of procedures
#*******************************************
proc init_ising { h t } {
#MD++ setnolog
MD++ setoverwrite
MD++ dirname = "runs/2D_ISING/MCrun_${h}_${t}"
MD++ srand48bytime 
MD++ srandbytime
#--------------------------------------------
# Initializtion
#
MD++ NX = 100  NY = 100 input = -1  initspin
#-----------------------------------------
# Physical parameters
MD++ J = 1     # interaction strength
MD++ H = ${h}  # applied field
MD++ kBT = ${t} # temperature
MD++ zipfiles = 1
} 

#--------------------------------------------
proc openwindow { } { MD++ {
# Plot Configuration
atomradius = 1.0 win_width = 480 win_height = 480 rotateangles = [ 0 0 0 1.7]
color00 = white color01 = blue backgroundcolor = gray80
openwin alloccolors rotate saverot plot
} }
#end of proc openwindow

proc run_setup { } {
MD++ {
	saveprop = 1 savepropfreq=1 plotfreq = 100
        savecn   = 0 calcryfreq=1 calcryfreq1=250
	openpropfile
      }
}	

proc compute_M { kBT } {
set M [ expr pow(1.0-1.0/pow(sinh(2.0/$kBT),4.0),1.0/8.0) ];
return $M
}

proc equil_run { totsteps } {
   MD++ {
        allocDFS
        calcryfreq1 = 250
	calcryfreq = 1 
   }
   MD++  totalsteps = $totsteps MCrun
}

proc datafile_process { filename index frac fracend operation } {
   set fp [ open $filename r ]
   set data [ split [read $fp] \n]
   close $fp
   set NUM [ llength $data ]
   set NUM [ expr $NUM-1 ]
   set Nini [ expr round($NUM*$frac) ]
   set Nend [ expr round($NUM*$fracend) ]
#   puts "total line $NUM \n"
   set Sum 0
   set k 0
   set Var 0
   set Std 0
   for { set i $Nini } { $i < $Nend } { incr i 1 } {
       set k [expr $k+1]
       set data1 [ lindex $data $i ]
       split $data1
       set datum [ lindex $data1 $index ]
       if { $i == $Nini } {
        set MAX [lindex $data1 $index]
        set MIN [lindex $data1 $index]
       } elseif { $i > $Nini } {
        set MAX [ expr ($MAX>$datum)?$MAX:$datum ]
        set MIN [ expr ($MIN<$datum)?$MIN:$datum ]
       }
#       puts "$datum"
       set Sum [expr $Sum+$datum]
   }
   set Ave [ expr $Sum/$k]
#   puts "$Sum $Ave"
#   puts "$MAX $MIN"

    if { [ string match "*STD*" $operation ] || [ string match "*VAR*" $operation ] } {
        for { set i $Nini } { $i < $Nend } { incr i 1 } {
          set data1 [ lindex $data $i ]
          split $data1
          set datum [ lindex $data1 $index ]
          set Var [expr $Var+($datum-$Ave)*($datum-$Ave) ]
        }
        set Var [ expr $Var/$k ]
        set Std [ expr sqrt($Var) ]
   }
   split $operation
   set Nvar [ llength $operation ]
   for { set i 0 } { $i < $Nvar } { incr i 1 } {
     set var [ lindex $operation $i ]
     if { $var=="SUM" } {
        lappend LIST $Sum
     } elseif { $var=="AVE" } {
        lappend LIST $Ave
     } elseif { $var=="MAX" } {
        lappend LIST $MAX
     } elseif { $var=="MIN" } {
        lappend LIST $MIN
     } elseif { $var=="VAR" } {
        lappend LIST $Var
     } elseif { $var=="STD" } {
        lappend LIST $Std
     }
   }
#   puts "$Std $Var" 
#   puts "$LIST"
   return $LIST
}



#*******************************************
# Main program starts here
#*******************************************
# status 0:prepare liquid by MD (NPT) simulation
#        1:MC equilibration
#        2:compute free energy by MC switching simulation
#
# read in status from command line argument

  # Monte Carlo simulation (brute force)
   
  if { $graphic=="On" } {
  MD++ setnolog 
  MD++ printfreq = 100
  } 

  init_ising $h $t

  if { $graphic=="On" } {
  openwindow
  }   
  run_setup
  equil_run $totsteps

  set Stotave [ datafile_process "prop.out" 1 0.2 1 "AVE" ]
  set Stotstd [ datafile_process "prop.out" 1 0.2 1 "STD" ]
  set Mag_fromMC [ expr 1.0*$Stotave/10000 ]
  set deltaM_MC [ expr 1.0*$Stotstd/10000 ]
  set Nc_ave [ datafile_process "prop.out" 4 0.2 1 "AVE" ]
  set Mag_analytic [ compute_M $t ]
  puts "Magnetization from Monte Carlo Simulation   =$Mag_fromMC +- $deltaM_MC"
  puts "Magnetization from Analytic Onsager Solution=$Mag_analytic"
  puts "Average Lagest Cluster Size = $Nc_ave"

MD++ quit

