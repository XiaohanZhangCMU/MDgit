# -*-shell-script-*-
# FEM mode of elastomer
#
# compile as: 
#  make fem build=R SYS=intel
#
# run as (in md++.pbs):
#  bin/fem_intel scripts/work/fem/fem-test.tcl 0


source "$::env(MDPLUS_DIR)/scripts/Examples/Tcl/startup.tcl"
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
#plotfreq = 10  rotateangles = [ 0 0 0 1.25 ] #[ 0 -90 -90 1.5 ]
plotfreq = 10  rotateangles = [ 0 0 0 1 ]
#win_width = 800 win_height = 800
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
 set pbid 1
} elseif { $argc > 1 } {
 set pbid [lindex $argv 1]
}
puts "pbid = $pbid"

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

  if { $pbid == 1 } { # 5element
    MD++ fem_coeff_file = "../../scripts/work/fem/membrane2d-5element-CPE4-Q4.dat"  read_fem_coeff
    MD++ elements_file = "../../scripts/work/fem/membrane2d-5element.ele"  read_elements
    MD++ incnfile = "../../scripts/work/fem/membrane2d-5element.cn" readcn
    MD++ map23 = \[ 0 1 \] 
    MD++ RtoRref
    #setup_window
    #openwindow

    MD++ test_saxpy
    MD++ sleep 
    MD++ eval

     MD++ R(15) = 3.0   #x(5)
     MD++ R(18) = 3.0   #x(6)
     MD++ R(16) = 0   #y(5)
#     MD++ R(17) = 0   #z(5)
     MD++ R(19) = 1   #y(6)
#     MD++ R(20) = 1   #z(6)
     MD++ input = \[ 4 0 11 5 6 \] fixatoms_by_ID
     MD++ RHtoS
  }

  if { $pbid == 2 } { # 8x5 elements

    MD++ fem_coeff_file = "../../scripts/work/fem/membrane2d-8x5element-CPE4-Q4.dat"  read_fem_coeff
    MD++ elements_file = "../../scripts/work/fem/membrane2d-8x5element.ele"  read_elements
    MD++ incnfile = "../../scripts/work/fem/membrane2d-8x5element.cn" readcn
    MD++ fem_bdy_nodes_file = "../../scripts/work/fem/membrane2d-8x5element.bdy" read_bdy_nodes
    MD++ map23 = \[ 0 1 \] 

    set h1 [ MD++_Get H_11 ]
    set h2 [ MD++_Get H_22 ]
    set h3 [ MD++_Get H_33 ]
    set i1 [ MD++_Get map23(0) ]
    set i2 [ MD++_Get map23(1) ]

    puts "$h1 $h2 $h3 $i1 $i2"
    MD++ RtoRref    
    #setup_window
    openwindow

#    puts bdy_group
    set n_bdy_nodes [MD++_Get n_bdy_nodes]
    set n_bdy_tags [MD++_Get n_bdy_tags]
    puts "n_bdy_nodes = $n_bdy_nodes"
    puts "n_bdy_tags = $n_bdy_tags"

    if { 1 } {
    set curr_pos 10.0    
    set input_temp {}
    for { set niter 0 } { $niter <= $n_bdy_nodes } { incr niter 1 } {
	set bnd [ MD++_Get bNds($niter) ]
	set bid [ MD++_Get bTags($niter) ]
	set xb  [ MD++_Get bXs($niter) ]
	set yb  [ MD++_Get bYs($niter) ]
	set zb  [ MD++_Get bZs($niter) ]
	
	if { $bid == 2 } {	    
	    MD++ R([expr $bnd * 3 ]) = $curr_pos

	    puts "yb = $yb, h2 = $h2"
#            MD++ R([expr $bnd * 3 + 1 ]) = [expr $yb*$h2]
#            MD++ R([expr $bnd * 3 + 2 ]) = $zb
	    MD++ fixed($bnd) = 1
	}
	if { $bid == 0 } {	
	    MD++ fixed($bnd) = 1    
	}
    }
    }
     MD++ RHtoS
  }

  if { $pbid == 3 } {


    MD++ fem_coeff_file = "../../scripts/work/fem/membrane2d-8x5element-CPE4-Q4.dat"  read_fem_coeff
    MD++ elements_file = "../../scripts/work/fem/membrane2d-8x5element.ele"  read_elements
    MD++ incnfile = "../../scripts/work/fem/membrane2d-8x5element.cn" readcn
    MD++ fem_bdy_nodes_file = "../../scripts/work/fem/membrane2d-8x5element.bdy" read_bdy_nodes
    MD++ map23 = \[ 0 1 \] 

    MD++ RtoRref    
    setup_window
    openwindow

#    MD++ sleep

#    puts bdy_group
    set n_bdy_nodes [MD++_Get n_bdy_nodes]
    set n_bdy_tags [MD++_Get n_bdy_tags]
    puts "n_bdy_nodes = $n_bdy_nodes"
    puts "n_bdy_tags = $n_bdy_tags"

    if { 1 } {
    for { set niter 0 } { $niter < $n_bdy_nodes } { incr niter 1 } {
	puts "niter = $niter"

	set bnd [ MD++_Get bNds($niter) ]
	set bid [ MD++_Get bTags($niter) ]
	set xb  [ MD++_Get bXs($niter) ]
	set yb  [ MD++_Get bYs($niter) ]
	set zb  [ MD++_Get bZs($niter) ]
	
	puts "bnd = $bnd"
	if { $bid == 0 } {
	    puts "$bnd is fixed"
	    MD++ fixed($bnd) = 1
	    
	}
	if { $bid == 2 } {
	    puts "$bnd is fixed"
	    MD++ fixed($bnd) = 1
	}
    }

   set Lx0 [ MD++_Get H_11 ]
   set Ly0 [ MD++_Get H_22 ]
   set Lz0 [ MD++_Get H_33 ]

   set epsilon 0.5
   set Lx_fix  [expr $Lx0*(1+$epsilon) ] 
   MD++ H_11 = $Lx_fix
#   MD++ fixallatoms
   }
#     MD++ RHtoS
  }

  MD++ eval  plot

  MD++ { conj_fixboxvec = [ 1 0 1
                            0 0 1
			    1 1 1 ] }

#  MD++ input = \[ 1 1 0.5 \] changeH_keepR
  MD++ relax 
 # relax_fixbox
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



