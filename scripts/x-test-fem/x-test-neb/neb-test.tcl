# -*-shell-script-*-
# compile as: 
# make xfem build=R SYS=mc2
# bin/xfem_mpich scripts/ox-test-fem/x-test-neb/neb-test.tcl 0 [1-4] [1-5]


source "$::env(MDPLUS_DIR)/scripts/Examples/Tcl/startup.tcl"
#*******************************************
# Definition of procedures
#*******************************************
proc initmd { status pbid meshid} {
MD++ setnolog
MD++ setoverwrite
if {$status == 20} {
   MD++ dirname = runs/x-test-fem-view-$pbid-$meshid
} else {
   MD++ dirname = runs/x-test-fem-$status-$pbid-$meshid
}
MD++ NNM = 400
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

proc setup_window { } { MD++ {
#Plot Settings
#
atomradius = 0.01 bondradius = 0.01 bondlength = 0 #for Si
atomcolor = orange highlightcolor = purple  bondcolor = red backgroundcolor = gray70
#plot_color_bar = [ 1 -4.85 -4.50 ]  highlightcolor = red
#plotfreq = 10  rotateangles = [ 0 0 0 1.25 ] #[ 0 -90 -90 1.5 ]
plotfreq = 10  rotateangles = [ 0 0 0 1 ]
#win_width = 800 win_height = 800
} }

#-------------------------------------------------------------
proc openwindow { } {
#Configure and open a plot
setup_window
MD++ openwin  alloccolors rotate saverot eval plot
}

#-------------------------------------------------------------
proc relax_fixbox { } { MD++ {
#Conjugate-Gradient relaxation
#
conj_ftol = 1e-7 conj_itmax = 1000 conj_fevalmax = 10000000
conj_fixbox = 1 #conj_monitor = 1 conj_summary = 1
relax 
} }

#-------------------------------------------------------------
proc relax_freebox { } { MD++ {
#Conjugate-Gradient relaxation
#
conj_ftol = 1e-7 conj_itmax = 1000 conj_fevalmax = 10000
conj_fixbox = 0 
#conj_monitor = 1 conj_summary = 1
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
 set meshid 1
} elseif { $argc > 2 } {
 set meshid [lindex $argv 2]
}
puts "meshid = $meshid"

set myname  [ MD++_Get "myname"]
puts "myname = $myname"


  initmd $status $pbid $meshid
  MD++ atommass = 1.0

  set meshfolder "../../scripts/x-test-fem/fe-meshes/M"
  set meshfolder $meshfolder$meshid
  set datafile [ glob -directory $meshfolder *.dat ]
  set fecnfile [ glob -directory $meshfolder *.cn  ]
  set elefile  [ glob -directory $meshfolder *.ele ]
  set bdyfile  [ glob -directory $meshfolder *.bdy ]

  MD++ fem_coeff_file = $datafile read_fem_coeff
  MD++ elements_file = $elefile  read_elements
  MD++ incnfile = $fecnfile  readcn
  MD++ fem_bdy_nodes_file = $bdyfile  read_bdy_nodes

  if { $meshid == 1 || $meshid == 2 || $meshid == 5  || $meshid == 7 || $meshid == 8 } {
    MD++ map23 = \[ 0 1 \] 
  } elseif { $meshid == 3 || $meshid == 4 || $meshid == 6 } {
    MD++ map23 = \[ 0 2 \] 
  } else {
    puts "mesh id is not valid. must be 1 -- 7"
    MD++ quit
  }

  # n_bdy_nodes have duplicates in it. all corner nodes are double counted
  set n_bdy_nodes [MD++_Get n_bdy_nodes]
  set n_bdy_tags [MD++_Get n_bdy_tags]
  set Lx0 [ MD++_Get H_11 ]
  set Ly0 [ MD++_Get H_22 ]
  set Lz0 [ MD++_Get H_33 ]


if { $status == 0 } {
  puts "n_bdy_nodes = $n_bdy_nodes -----------------------"
  MD++ RtoRref
  #setup_window
  openwindow

  if { $pbid == 1 } { # 5element    	  
    for { set niter 0 } { $niter <= $n_bdy_nodes } { incr niter 1 } {
	set bnd [ MD++_Get bNds($niter) ]
	set bid [ MD++_Get bTags($niter) ]
	set xb  [ MD++_Get bXs($niter) ]
	set yb  [ MD++_Get bYs($niter) ]
	set zb  [ MD++_Get bZs($niter) ]
	
	if { $bid == 2 } {	    
	    MD++ R([expr $bnd * 3 ]) = [ expr 2 * $xb ]
	    MD++ fixed($bnd) = 1
	    set res [ MD++_Get R([expr $bnd*3]) ]
	    puts "R = $res"
	}
	if { $bid == 0 } {	
	    MD++ fixed($bnd) = 1    
	} 
    }      
    MD++ RHtoS

  } elseif { $pbid == 2 || $pbid == 3 } { 	  
    for { set niter 0 } { $niter < $n_bdy_nodes } { incr niter 1 } {
	set bnd [ MD++_Get bNds($niter) ]
	set bid [ MD++_Get bTags($niter) ]
	puts "bid = $bid"
	if { $bid == 0 || $bid == 2 } {
	    puts "bnd = $bnd is fixed"
	    MD++ fixed($bnd) = 1	    
	}
    }
    if { $pbid == 2 } {
       set epsilon 0.5
    } elseif { $pbid == 3 } {
       set epsilon -0.5
    }
    set Lx_fix  [expr $Lx0*(1+$epsilon) ] 
    puts "Lx = $Lx_fix"
    MD++ H_11 = $Lx_fix
#   You cannot put RHtoS here, since it will scale the membrane w.r.t the box.
#   and then no elastic energy can be stored 
#   MD++ RHtoS 

  } elseif { $pbid == 4 } {
    for { set niter 0 } { $niter < $n_bdy_nodes } { incr niter 1 } {
	set bnd [ MD++_Get bNds($niter) ]
	set bid [ MD++_Get bTags($niter) ]
	set xb  [ MD++_Get bXs($niter) ]
	set yb  [ MD++_Get bYs($niter) ]
	set zb  [ MD++_Get bZs($niter) ]

	# puts "bid = $bid"
	# perturb to get two modes
	set modeII 1
	set perturb 0.0001
	#set perturb 0.0

	if {$modeII ==0} {
        if { $bid == 1 || $bid == 3 } {		
	   set PI 3.1415926
	   set newR [expr $yb + $perturb*cos($PI/$Lx0 * $xb)]
	   set oldR [MD++_Get R([expr $bnd*3+1])]
	   set barlength 5
	   set pert [expr $perturb*cos($PI/$barlength * $xb)]
	   MD++ R([expr $bnd*3 + 1]) = [expr $yb + $perturb*cos($PI/$barlength * $xb)]
#           puts "yb = $yb, xb = $xb, oldR = $oldR, newR = $newR, pert = $pert"
        } 
        }
	if { $bid == 0 || $bid == 2 } {
	    puts "bnd = $bnd is fixed"
	    MD++ fixed($bnd) = 1	    
	}
    }
   if {$modeII == 0} {
       MD++ RHtoS
   }
#    MD++ sleep

#perturb = 0.0001
#-0.0125: no curvature
#-0.015: mode I curved

#perturb = 0
#-0.3 -- 0.0345: rubbed no curvature
#-0.5: mode I curved

    set epsilon -0.015
    set Lx_fix  [expr $Lx0*(1+$epsilon) ] 
    MD++ H_11 = $Lx_fix


  } else {
	 puts "pbid value not valid"
  }

  MD++ { conj_fixboxvec = [ 1 0 1
                            0 0 1
			    1 1 1 ] }
  #MD++ relax
  relax_fixbox
  MD++ plot
  MD++ writeall = 1
  MD++ finalcnfile = "NEBfinal.cn" writecn
  MD++ sleep
  exitmd 

} elseif { $status == 4 } {
  MD++ RtoRref
  
  set ncpu [ MD++_Get numDomains ]
  if { $ncpu < 1 } { set ncpu 1 }
  puts "ncup = $ncpu"
  
  MD++ Slave_chdir
  set initNEBfile ../x-test-fem-0-$pbid-$meshid/NEBinit.cn
  set finalNEBfile ../x-test-fem-0-$pbid-$meshid/NEBfinal.cn 
  MD++ incnfile = $initNEBfile readcn saveH
  MD++ incnfile = $finalNEBfile readcn restoreH SHtoR setconfig2
  MD++ eval
  MD++ incnfile = $initNEBfile readcn setconfig1
  MD++ eval

  # rename previous result files in order to not overwrite them
  set check [ glob -nocomplain "stringrelax.chain.cn.cpu00" ]
  set exist [ llength $check ]
  if { $exist > 0 } {
     for { set i 0 } { $i < $ncpu } { incr i 1 } {
         set ival [ format %02d $i ]
         exec cp stringrelax.chain.cn.cpu$ival Prev_stringrelax.chain.cn.cpu$ival
     }
     exec cp stringeng.out Prev_stringeng.out
  }

  set chainlength [expr $ncpu - 1]
  MD++ chainlength = $chainlength  
  MD++ timestep = 0.01 printfreq = 1

  MD++ nebspec = \[ 0 1 0 1 0 0 1 \] totalsteps = 1000 equilsteps = 10000

  MD++ eval 

  MD++ Broadcast_Atoms
  MD++ Broadcast_FEM_Param

  MD++ fixallatoms  constrain_fixedatoms freeallatoms
  MD++ allocchain_parallel initRchain_parallel

  if { $exist > 0 } {
    MD++ incnfile = stringrelax.chain.cn readRchain_parallel
  }

  MD++ eval
  MD++ stringrelax_parallel_2
  MD++ finalcnfile = stringrelax.chain.cn writeRchain_parallel

  if { $ncpu == 1 }   {
        MD++ quit
  } else {
        MD++ quit_all
  }
 

} elseif { $status == 20 } {
  MD++ setnolog
  initmd $status $pbid $meshid 

  set chain_no 0
  set total_no 7
  MD++ chainlength = $total_no

  MD++ incnfile = "../x-test-fem-0-4-7/NEBinit.cn" readcn  saveH
  MD++ incnfile = "../x-test-fem-4-4-7/stringrelax.chain.cn" 
  MD++ input = $chain_no  readRchain_parallel_toCN  RHtoS

  #MD++ clearR0 refreshnnlist eval quit

#  setup_window

  MD++ zipfiles = 0
  for { set iter 0 } { $iter <= $total_no } { incr iter 1 } {
     # set chain number
     set chain_no $iter
     MD++ input = $chain_no  readRchain_parallel_toCN  RHtoS
     MD++ NCS = 8 eval calcentralsymmetry 
     
     MD++ finalcnfile = "chain_no_${chain_no}.cfg" writecn
     MD++ finalcnfile = "chain_no_${chain_no}.cfg" writeatomeyecfg
     MD++ sleepseconds = 1  sleep
  }
  MD++ sleepseconds = 100  sleep
  exitmd
 
} else {        
  puts "unknown status = $status"
  exitmd 
} 



