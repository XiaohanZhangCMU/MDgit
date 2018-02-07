# -*-shell-script-*-
# relax buckball via conjugate gradient method and run MD
# to compile:
#  make rebo build=R SYS=intel
# to run:
#  bin/rebo_gcc scripts/Examples/example08-buckyball.tcl 0
#  bin/rebo_gcc scripts/Examples/example08-buckyball.tcl 1

source "scripts/Examples/Tcl/startup.tcl"

#*******************************************
# Definition of procedures
#*******************************************
proc initmd { } {
MD++ setnolog
MD++ setoverwrite
MD++ dirname = runs/buckyball
}

#------------------------------------------------------------
proc read_rebo { } { MD++ {
rebofile = "~/Codes/MD++/potentials/AIREBO/CH.airebo" 
rcut = 3.0  # should be the biggest cut-off radius of all interactions
readREBO
} }

#--------------------------------------------
proc relax_fixbox { } { MD++ {
# Conjugate-Gradient relaxation
conj_ftol = 1e-7 conj_itmax = 1000 conj_fevalmax = 3000
conj_fixbox = 1
relax
} }
#end of proc relax_fixbox

#--------------------------------------------
proc relax_freebox { } { MD++ {
# Conjugate-Gradient relaxation
conj_ftol = 1e-7 conj_itmax = 1000 conj_fevalmax = 3000
conj_fixbox = 0  conj_fixboxvec = [ 1 1 1   1 1 1   1 1 0 ]
relax
} }
#end of proc relax_fixbox

proc setup_window { } { MD++ {
#------------------------------------------------------------
#colors for Central symmetry view
color00 = "red" color01 = "blue" color02 = "green"
color03 = "magenta" color04 = "cyan" color05 = "purple"
color06 = "gray80" color07 = "white" color08 = "orange"
#--------------------------------------------
# Plot Configuration
#
atomradius = [0.4 0.28] bondradius = 0.15 bondlength = 1.9
win_width=600 win_height=600
atomcolor = cyan highlightcolor = purple  bondcolor = red
fixatomcolor = yellow backgroundcolor = gray70
plot_color_axis = 0
plot_color_windows = [ 3
                       -10.0 -0.90      1  
                       -0.90 -0.50      5
                       -0.50  10.0      4 
                       1
                     ]
#plot_atom_info = 1 # reduced coordinates of atoms
#plot_atom_info = 2 # real coordinates of atoms
plot_atom_info = 3 # energy of atoms
plotfreq = 1
rotateangles = [ 0 0 0 1.2 ]
} }


proc openwindow { } {
setup_window
MD++ openwin alloccolors rotate saverot eval plot
}

#--------------------------------------------
proc exitmd { } { MD++ quit }
#end of proc exitmd

#--------------------------------------------
proc setup_md { } { MD++ {     
equilsteps = 0  timestep = 0.001 # (ps)
atommass = [12.0107 1.00794] # (g/mol) for Carbon and Hydrogen
DOUBLE_T = 1
saveprop = 1 savepropfreq = 100 openpropfile #run
savecn = 1 savecnfreq = 10000 openintercnfile
savecfg = 1 savecfgfreq = 2500
plotfreq = 10 printfreq = 10
vt2 = 3e26  #1e28 2e28 5e28
wallmass = 2e3     # atommass * NP = 14380
boxdamp = 1e-3     # optimal damping for 216 atoms and wallmass 1e-3
saveH # Use current H as reference (H0), needed for specifying stress
fixboxvec  = [ 1 1 1 
               1 1 1
               1 1 0 ]
stress = [ 0 0 0
           0 0 0
           0 0 0 ]
output_fmt = "curstep EPOT KATOM Tinst TSTRESSinMPa_zz" 
ensemble_type = "NVT" integrator_type = "VVerlet" implementation_type = 1
T_OBJ = 300 #(in K) # Desired Temperature
zerorot = "all" # zero out angular momentum when initvelocity
writeall = 1 # write velocity and species information into file
} }





#*******************************************
# Main program starts here
#*******************************************
# status 0:
#        1:
#        2:
#
# read in status from command line argument
if { $argc == 0 } {
 set status 0
} elseif { $argc > 0 } {
 set status [lindex $argv 0]
}
puts "status = $status"

set myname  [ MD++_Get "myname"]



if { $status == 0 } {
  # read buckyball structure and relax

  initmd
  read_rebo
 
  MD++ incnfile = "~/Codes/MD++/structures/Examples/buckyball.cn" readcn
  MD++ nspecies = 1  element0 = "C" freeallatoms

  openwindow

  # debug
  MD++_PrintArray EPOT_IND "eV" 0 60 "%12.6e"
  MD++ eval

  relax_fixbox
  MD++ finalcnfile = "buckyball-relaxed.cn" writecn

  MD++ sleep

  exitmd

} elseif { $status == 1 } {
  # read buckyball structure and run MD

  initmd
  read_rebo
 
  MD++ incnfile = "buckyball-relaxed.cn" readcn
  MD++ nspecies = 1  element0 = "C" freeallatoms

  openwindow

  setup_md
  MD++ T_OBJ = 600  initvelocity 
  MD++ totalsteps = 100000   run 

  MD++ finalcnfile = equil.cn writeall = 1 writecn

  MD++ sleep

  exitmd

} else {
        
 puts "unknown status = $status"
 exitmd 

} 

