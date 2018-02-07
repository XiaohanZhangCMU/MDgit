# -*-shell-script-*-
# Ising 2D Monte Carlo simulation
#-------------------------------------------

source "scripts/Examples/Tcl/startup.tcl"

#*******************************************
# Definition of procedures
#*******************************************
proc init_ising { n size } {
MD++ setnolog
MD++ setoverwrite
MD++ dirname = "runs/ising3d-$n"
MD++ srand48bytime 

#--------------------------------------------
# Initializtion
#
MD++ NX = $size  NY = $size   NZ = $size  input = -1  initspin
#incnfile = ../ising2d/s20.cn readcn
#--------------------------------------------
# Physical parameters
MD++ J = 1     # interaction strength
#MD++ H = 0.065  kBT = 1.8  
MD++ H = 1.20 kBT = 1.8
#MD++ H = 0 kBT = 2.25  # roughening temperature
#MD++ H = 0 kBT = 3.0
#MD++ H = 0 kBT = 4.51149 # critical temperature
} 

#--------------------------------------------
proc openwindow { } { MD++ {
# Plot Configuration
#
atomradius = 1.0 win_width = 600 win_height = 600 rotateangles = [ 0 0 0 1.2]
color00 = white color01 = blue backgroundcolor = gray80
openwin alloccolors rotate saverot plot
} }
#end of proc openwindow





#*******************************************
# Main program starts here
#*******************************************
# status 0:prepare liquid by MD (NPT) simulation
#        1:MC equilibration
#        2:compute free energy by MC switching simulation
#
# read in status from command line argument
if { $argc == 0 } {
 set status 0
} elseif { $argc > 0 } {
 set status [lindex $argv 0]
}
puts "status = $status"

if { $argc <= 1 } {
 set n 0
} elseif { $argc > 1 } {
 set n [lindex $argv 1]
}
puts "n = $n"


if { $status == 0 } {
  # Monte Carlo simulation (brute force)
  set size 50
  init_ising $n $size
  openwindow
  #MD++ sleep quit
  MD++ totalsteps = 1000001 finalcnfile = s$size.cn writecn 
  MD++ plotfreq = 20 #1000 
  MD++ saveprop = 1 savepropfreq = 100 printfreq = 1000
  MD++ savecn = 1   savecnfreq   = 10000 openintercnfile
  MD++ srand48bytime openpropfile MCrun
  MD++ finalcnfile = s$size-a.cn writecn
  
  exitmd

} else {

  puts "unknown status = $status"
  exitmd

}
