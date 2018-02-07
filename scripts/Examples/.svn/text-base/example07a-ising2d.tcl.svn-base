# -*-shell-script-*-
# Ising 2D Monte Carlo simulation
#-------------------------------------------

source "scripts/Examples/Tcl/startup.tcl"

#*******************************************
# Definition of procedures
#*******************************************
proc init_ising { n } {
MD++ setnolog
MD++ setoverwrite
MD++ dirname = "runs/ising2d-$n"
MD++ srand48bytime 

#--------------------------------------------
# Initializtion
#
MD++ NX = 100  NY = 100 input = -1  initspin
#incnfile = ../ising2d/s100.cn readcn
#--------------------------------------------
# Physical parameters
MD++ J = 1     # interaction strength
MD++ H = 0.065  # applied field
#MD++ H = 0.06  # applied field
#MD++ kBT = 1.5 # temperature
MD++ kBT = 1.8 # temperature
} 

#--------------------------------------------
proc openwindow { } { MD++ {
# Plot Configuration
#
atomradius = 1.0 win_width = 480 win_height = 480 rotateangles = [ 0 0 0 1.7]
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
  init_ising $n
  openwindow
  MD++ totalsteps = 1000001 finalcnfile = s100.cn writecn 
  MD++ plotfreq = 1000 
  MD++ saveprop = 1 savepropfreq = 100 printfreq = 1000
  MD++ savecn = 1   savecnfreq   = 10000 openintercnfile
  MD++ srand48bytime openpropfile MCrun
  MD++ finalcnfile = s100-a.cn writecn
  
  exitmd

} elseif { $status == 1 } {
  # Path Sampling (manual test)
  init_ising $n
  openwindow
  MD++ input =  1  initpath finalcnfile = path0.pcn writepath
  MD++ input =  3  initpath
  MD++ input = \[ 1 3500 \] copyPathtoS 
  MD++ input = \[ 0.03 0.5 10000000 1 \] 
  MD++ plotfreq = 1000 printfreq = 1000 SampleMCpath 
  MD++ finalcnfile = path1.pcn input = 1  writepath 
  MD++ nsample = 100 ComputeSuccessRate 
  exitmd

} elseif { $status == 2 } {
  # Path Sampling (automatic)
  init_ising $n
  openwindow
  MD++ input = \[ 0.03 0.5 100000000 \] plotfreq = 1 printfreq = 10000
  MD++ totalsteps = 1000   savecn = 1  savecnfreq = 10  openintercnfile
  MD++ saveprop = 1 savepropfreq = 1 printfreq = 1000 openpropfile
#  MD++ nsample = 100 
  MD++ WalkonChain

  MD++ finalcnfile = pathA.pcn input =  1  writepath
  MD++ finalcnfile = pathB.pcn input =  2  writepath
  MD++ finalcnfile = pathC.pcn input =  3  writepath

  exitmd

} elseif { $status == 3 } {
  # Saddle points analysis
  init_ising $n
  openwindow
  MD++ input = \[ 1 5 0.03 0.5 100000000 \] nsample = 100 plotfreq = 10000 printfreq = 10000
  MD++ AnalyzeConfigs

  exitmd

} else {

  puts "unknown status = $status"
  exitmd

}
