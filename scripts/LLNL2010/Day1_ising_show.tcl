# -*-shell-script-*-
# Ising 2D Monte Carlo simulation
#-------------------------------------------

source "scripts/Examples/Tcl/startup.tcl"

#*******************************************
# Definition of procedures
#*******************************************
proc init_ising { } {
MD++ setnolog
MD++ setoverwrite
MD++ srand48bytime 

#--------------------------------------------
# Initializtion
#
MD++ NX = 32  NY = 32 NZ =32 input = -1  initspin
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
# read in status from command line argument
if { $argc == 0 } {
 set status 0
} elseif { $argc > 0 } {
 set status [lindex $argv 0]
}
  init_ising 
  MD++ N_lgst_index = 1
  MD++ incnfile=$status readcn
  MD++ allocDFS
  MD++ calcrystalorder zipfiles = 0 
  openwindow
  set Nlgst [MD++_Get N_lgst_cluster ]
  puts "Nlgst_cluster=$Nlgst"
  MD++ finalcnfile=temp.cn writecn
  MD++ sleep 
