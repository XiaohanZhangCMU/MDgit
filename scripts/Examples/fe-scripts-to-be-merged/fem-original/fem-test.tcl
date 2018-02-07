# -*-shell-script-*-
# FEM mode of elastomer
#
# compile as: 
#  make fem build=R SYS=intel
#
# run as (in md++.pbs):
#  bin/fem_intel scripts/work/fem/fem-test.tcl 0

source "scripts/Examples/Tcl/startup.tcl"

#*******************************************
# Definition of procedures
#*******************************************
proc initmd { status } {
MD++ setnolog
MD++ setoverwrite
MD++ dirname = runs/fem-test-$status
MD++ NNM = 200
}

#--------------------------------------------
proc makeperfcrystal { NX NY NZ a } { 
# Create Perfect Lattice Configuration
#
MD++ crystalstructure = diamond-cubic
MD++ latticeconst = $a
MD++ element0 = Silicon
MD++ latticesize = \[ 1 1 1  $NX   -1  -1  2 $NY 0.5 -0.5 0 1 \]
MD++ makecrystal #finalcnfile = perf.cn writecn #eval
}

#--------------------------------------------
proc make_membrane_2d { NX NY } { 
# make a 2d membrane
#
#MD++ input = \[3 2 0 0 $radius 0 \]
#MD++ makecylinder #finalcnfile = cylinder.cn writecn
} 

#-------------------------------------------------------------
proc setup_window { } { MD++ {
#Plot Settings
#
atomradius = 0.1 bondradius = 0.01 bondlength = 0 #for Si
atomcolor = orange highlightcolor = purple  bondcolor = red backgroundcolor = gray70
#plot_color_bar = [ 1 -4.85 -4.50 ]  highlightcolor = red
plotfreq = 10  rotateangles = [ 0 0 0 1.25 ] #[ 0 -90 -90 1.5 ]
} }

#-------------------------------------------------------------
proc openwindow { } {
#Configure and open a plot
#
setup_window
MD++ openwin  alloccolors rotate saverot eval plot
}

#-------------------------------------------------------------
proc relax_fixbox { } { MD++ {
#Conjugate-Gradient relaxation
#
conj_ftol = 1e-7 conj_itmax = 1000 conj_fevalmax = 10000
conj_fixbox = 1 #conj_monitor = 1 conj_summary = 1
relax 
} }

#-------------------------------------------------------------
proc relax_freebox { } { MD++ {
#Conjugate-Gradient relaxation
#
conj_ftol = 1e-7 conj_itmax = 1000 conj_fevalmax = 10000
conj_fixbox = 0 #conj_monitor = 1 conj_summary = 1
relax 
} }

#-------------------------------------------------------------
proc setup_md { } { MD++ {
#MD settings (without running MD)
#           
equilsteps = 0  totalsteps = 100 timestep = 0.0001 # (ps)
atommass = 28.0855 # (g/mol)
atomTcpl = 200.0 boxTcpl = 20.0
DOUBLE_T = 0
srand48bytime
initvelocity totalsteps = 10000 saveprop = 0
saveprop = 1 savepropfreq = 10 openpropfile
savecn = 1  savecnfreq = 100
writeall = 1
savecfg = 1 savecfgfreq = 100
ensemble_type = "NVT" integrator_type = "VVerlet"
implementation_type = 0 vt2=1e28
totalsteps = 20000
output_fmt = "curstep EPOT KATOM Tinst HELM HELMP TSTRESS_xx TSTRESS_yy TSTRESS_zz"
#output_fmt = "curstep EPOT KATOM Tinst HELM HELMP TSTRESS_xx TSTRESS_yy TSTRESS_zz"
plotfreq = 50 #autowritegiffreq = 10
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

if { $argc <= 1 } {
 set n 0
} elseif { $argc > 1 } {
 set n [lindex $argv 1]
}
puts "n = $n"

if { $argc <= 2 } {
 set m 1
} elseif { $argc > 2 } {
 set m [lindex $argv 2]
}
puts "m = $m"


set myname  [ MD++_Get "myname"]
puts "myname = $myname"

if { $status == 0 } {
  # create 2d membrane, relax

  initmd $status 
 
  # size of supercell
  #set NX 10
  #set NY 8
  #make_membrane_2d $NX $NY 

  MD++ atommass = 1.0
  MD++ fem_coeff_file = "../../scripts/work/fem/CPE4-Q4.dat"  read_fem_coeff

  #MD++ elements_file = "../../scripts/work/fem/membrane2d-1element.ele"  read_elements
  #MD++ incnfile = "../../scripts/work/fem/membrane2d-1element.cn" readcn

  MD++ elements_file = "../../scripts/work/fem/membrane2d-5element.ele"  read_elements
  MD++ incnfile = "../../scripts/work/fem/membrane2d-5element.cn" readcn

  MD++ RtoRref

  openwindow
  #setup_window

  MD++ eval

  if { 0 } { # 1element
    MD++ R(3) = 1.2   #x(1)
    #MD++ R(6) = 1.2   #x(2)
    MD++ R(7)  = 1.2   #y(2)
    MD++ R(9) = 0.1   #x(3)
    #MD++ R(10) = 1.2   #y(3)
    MD++ RHtoS
  }

  if { 1 } { # 5element
    MD++ R(15) = 3.0   #x(5)
    MD++ R(18) = 3.0   #x(6)
    MD++ R(16) = 0.5   #y(5)
    MD++ R(19) = 1.5   #y(6)
    MD++ input = \[ 4 0 11 5 6 \] fixatoms_by_ID
    MD++ RHtoS
  }

  MD++ eval  plot

  #MD++ finalcnfile = "membrane2d-init.cn" writecn
  #MD++ finalcnfile = "membrane2d-init.cn" writeatomeyecfg

  relax_fixbox
  MD++ plot
  MD++ writeall = 1
  MD++ finalcnfile = "membrane2d-relaxed.cn" writecn
  MD++ finalcnfile = "membrane2d-relaxed.cfg" writeatomeyecfg

  MD++ sleep
  exitmd 

} elseif { $status == 1 } {

  exitmd 

} elseif { $status == 8 } {
  # Molecular Dynamics simulation
  initmd $status 

#  MD++ incnfile = "../al-pillar-disl-0/pillar-relaxed.cn" readcn
#
#  #openwindow
#  setup_window
#
#  setup_md
#  MD++ T_OBJ = 300 totalsteps = 500000 initvelocity run
#  MD++ finalcnfile = "pillar-deformed.cn" writecn
#  MD++ eval

  exitmd

} elseif { $status == 10 } {
  # view
  initmd view $n $m

  MD++ incnfile = "../al-pillar-disl-0/pillar-relaxed.cn" readcn

  openwindow
  #setup_window

  MD++ sleep

  exitmd

} else {
        
 puts "unknown status = $status"
 exitmd 

} 



