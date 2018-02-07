# -*-shell-script-*-
# MD code of Stinger-Weber Silicon

#*******************************************
# Definition of procedures
#*******************************************
proc initmd { T } { MD++ {
setnolog
setoverwrite }
MD++ dirname = runs/si-example-$T
}

proc makecrystal { nx ny nz } { 
#Create Perfect Lattice Configuration
#
MD++ crystalstructure = diamond-cubic
MD++ latticeconst = 5.4309529817532409 #(A) for Si
MD++ element0 = Silicon
MD++ latticesize  = \[ 1 0 0 $nx 0 1 0 $ny 0 0 1 $nz \]
MD++ makecrystal  finalcnfile = perf.cn writecn
MD++ finalcnfile = perf.cfg
}

#-------------------------------------------------------------
proc setup_window { } { MD++ {
#Plot Settings
#
atomradius = 0.67 bondradius = 0.3 bondlength = 2.8285 #for Si
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

# read in status from command line argument
if { $argc == 0 } {
 set status 0
} elseif { $argc > 0 } {
 set status [lindex $argv 0]
} 
puts "status = $status"

if { $argc <=1 } {
  set T_OBJ 300
} elseif { $argc > 1 } {
 set T_OBJ [lindex $argv 1]
}
puts "T_OBJ = $T_OBJ"

initmd $T_OBJ

if { $status == 0 } {
  #Create crystal, relax, run MD/NVT simulation

  makecrystal  4 4 4

  openwindow

  relax_fixbox
  MD++ finalcnfile = relaxed.cn writecn
  MD++ finalcnfile = relaxed.lammps writeLAMMPS
  MD++ finalcnfile = relaxed.cfg writeatomeyecfg
  setup_md
  MD++ T_OBJ = $T_OBJ
  MD++ run 

  MD++ finalcnfile = si100.cn writecn
  MD++ finalcnfile = si100.cfg writeatomeyecfg  

  MD++ sleep quit

} elseif { $status == 1} {
  #Read in relaxed, NVT

  MD++ incnfile = relaxed.cn readcn
  openwindow
  setup_md
  MD++ T_OBJ = $T_OBJ
  MD++ run 
  MD++ finalcnfile = si100.cn writecn

  MD++ sleep quit

} elseif { $status == 2} {
  #Visualization of final structure

  MD++ incnfile = si100.cn readcn
  openwindow

  MD++ sleep quit

} else { 
  puts "unknown argument. status = $status" 
  MD++ quit 
}

MD++ quit

