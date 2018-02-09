# -*-shell-script-*-
# compile:
# make fem build=R SYS=mc2
# bin/fem_mpich scripts/work/fem/neb-test.tcl $status $pbid $meshid $equationType

source "$::env(MDPLUS_DIR)/scripts/Examples/Tcl/startup.tcl"
#*******************************************
# Definition of procedures
#*******************************************
proc initmd { status pbid meshid eqnType } {
  MD++ setnolog
  MD++ setoverwrite
  if { $status==19 || $status==20 || $status==21 || $status==22 || $status==23 || $status==24 } {
    MD++ dirname = runs/fem/fem-view-$status-$pbid-$meshid-$eqnType
  } else {
    MD++ dirname = runs/fem/fem-$status-$pbid-$meshid-$eqnType
  }
  exec cp -r ../../../scripts/work/fem/neb-test.tcl .
  MD++ NNM = 400
}

#--------------------------------------------

proc setup_window { } { MD++ {
#Plot Settings
  atomradius = 0.001 bondradius = 0.01 bondlength = 0 #for Si
  atomcolor = orange highlightcolor = purple  bondcolor = red backgroundcolor = gray70
  color00 = "orange"  color01 = "green" color02 = "purple"
  color03 = "magenta" color04 = "cyan"   color05 = "purple"
  color06 = "gray80"  color07 = "white"
  #plot_color_bar = [ 1 -4.85 -4.50 ]  highlightcolor = red
  plotfreq = 10  rotateangles = [ 0 0 0 6 ]
  #plot_atom_info = 1 # reduced coordinates of atoms
  plot_atom_info = 2 # real coordinates of atoms       
  win_width = 1000 win_height = 800
  plot_map_pbc = 1
} }

#-------------------------------------------------------------
proc openwindow { } {
#Configure and open a plot
  setup_window
  MD++ openwin  alloccolors rotate saverot plot
}

#-------------------------------------------------------------
proc relax_fixbox { } { MD++ {
#Conjugate-Gradient relaxation
  conj_ftol = 1e-2 conj_itmax = 2600 conj_fevalmax =30 # 3300000
  conj_fixbox = 1 #conj_monitor = 1 conj_summary = 1
  relax 
} }

#-------------------------------------------------------------
proc relax_freebox { } { MD++ {
#Conjugate-Gradient relaxation
  conj_ftol = 1e-7 conj_itmax = 1000 conj_fevalmax = 10000
  conj_fixbox = 0 
  #conj_monitor = 1 conj_summary = 1
  relax 
} }

#*******************************************
# Main program starts here
#*******************************************
# status 0:
#        1:
#        2:
# eqnType 0: beam. 1: snapband. 2: islam
if { $argc <= 3 } {
  puts "Need to specify status, pbid, meshid and eqnType"
  MD++ quit
} elseif { $argc > 3 } {
  set status [lindex $argv 0]
  set pbid [lindex $argv 1]
  set meshid [lindex $argv 2]
  set eqnType [lindex $argv 3]
}
puts "status = $status"
puts "pbid = $pbid"
puts "meshid = $meshid"
puts "eqnType = $eqnType"

set USEGPU 0
set myname [ MD++_Get "myname" ]
if { [ string match "*femcpu*" $myname ] } {
  set USEGPU 0
} elseif { [ string match "*femgpu*" $myname ] } {
  puts "femgpu, set USEGPU =1"
  set USEGPU 1
} else { 
  puts "USEGPU set default zero"
}
if { $USEGPU == 1 } {
  MD++ test_saxpy
}

initmd $status $pbid $meshid $eqnType
MD++ EquationType = $eqnType
set meshfolder "../../../scripts/work/fem/fe-meshes/M${meshid}"

set fecnfile [ glob -directory $meshfolder *.cn  ]
set elefile  [ glob -directory $meshfolder *.ele ]
set bdyfile  [ glob -directory $meshfolder *.bdy ]
set length [string length $elefile]
set filename [string range $elefile 0 [expr $length-5] ]

MD++ elements_file = $elefile  read_elements
puts "I am here 1"
MD++ fem_coeff_file = ${filename}-C3D8-Q8.dat read_fem_coeff
#MD++ fem_coeff_file = ${filename}-CPE4-Q4.dat read_fem_coeff
puts "I am here 2"
MD++ incnfile = $fecnfile  readcn
puts "I am here 3"
MD++ fem_bdy_nodes_file = $bdyfile  read_bdy_nodes

puts "I am here 4"

# 2-d mesh xy
#16. islam, fancy, expanded
#12. snap band 3d.
if { $meshid==1 || $meshid==2 || $meshid==5  || $meshid==7 || $meshid==8 || $meshid==11  || $meshid==14 || $meshid==15 || $meshid==16 || $meshid==19 || $meshid==20 || $meshid==30 } {
  MD++ map23 = \[ 0 1 \] 
#xz
} elseif { $meshid==3 || $meshid==4 || $meshid==6 } {
  MD++ map23 = \[ 0 2 \] 
# 3-d mesh
} elseif { $meshid==9 || $meshid==10 || $meshid==12 || $meshid==13 || $meshid==17 || $meshid==18 || $meshid==21 || $meshid==22 || $meshid==23  } {
  MD++ map23 = \[-1 -1 \]
} else {
  puts "mesh id is not valid."
  MD++ quit
}

# n_bdy_nodes have duplicates in it. all corner nodes are double counted
set n_bdy_nodes [MD++_Get n_bdy_nodes]
set n_bdy_tags [MD++_Get n_bdy_tags]
set Lx0 [ MD++_Get H_11 ]
set Ly0 [ MD++_Get H_22 ]
set Lz0 [ MD++_Get H_33 ]
puts "n_bdy_nodes = $n_bdy_nodes"


if { $status == 0 } {
  MD++ RtoRref
  if { $USEGPU == 1 } {
    MD++ cuda_memcpy_all
  }
#  setup_window
#  openwindow

  if { $pbid == 1 } {     	  
    for { set niter 0 } { $niter <= $n_bdy_nodes } { incr niter 1 } {
      set bnd [ MD++_Get bNds($niter) ]
      set bid [ MD++_Get bTags($niter) ]
      set xb  [ MD++_Get bXs($niter) ]
      set yb  [ MD++_Get bYs($niter) ]
      set zb  [ MD++_Get bZs($niter) ]
      
      if { $bid == 1 } {	    
        MD++ R([expr $bnd * 3 ]) = [ expr 3.6 * $xb ]
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

  } elseif { $pbid == 4 } {
    for { set niter 0 } { $niter < $n_bdy_nodes } { incr niter 1 } {
      set bnd [ MD++_Get bNds($niter) ]
      set bid [ MD++_Get bTags($niter) ]
      set xb  [ MD++_Get bXs($niter) ]
      set yb  [ MD++_Get bYs($niter) ]
      set zb  [ MD++_Get bZs($niter) ]

      # perturb to get two modes
      set modeII 0
      set perturb 0.0001
      if {$modeII ==0} {
        if { $bid == 1 || $bid == 3 } {		
	  set PI 3.1415926
	  set oldR [MD++_Get R([expr $bnd*3+1])]
          set barlength 5
	  set pert [expr $perturb*cos($PI/$barlength * $xb)]
	  MD++ R([expr $bnd*3+1])=[expr $yb+$perturb*cos($PI/$barlength*$xb)]
#	  MD++ R([expr $bnd*3+1])=[expr $yb+$perturb*(floor(rand()*2)-1)]
        } 
      }
      if { $bid == 0 || $bid == 2 } {
	MD++ fixed($bnd) = 1	    
      }
    }
    if { $modeII == 0 || $modeII == 1 } {
      MD++ RHtoS
    }
#    MD++ sleep

#perturb = 0.0001, -0.0125: no curvature, -0.015: mode I curved

#perturb = 0, -0.3 -- 0.0345: rubbed no curvature, -0.5: mode I curved
    if { $modeII == 0 } { 
      set epsilon -0.035
    }
    if { $modeII == 1 } {
      set epsilon -0.035
    }
    set Lx_fix  [expr $Lx0*(1+$epsilon) ] 
    MD++ H_11 = $Lx_fix

# An equivalent pid as pid=4. Create buckled 3d beam
  } elseif { $pbid == 5 } {
    for { set niter 0 } { $niter < $n_bdy_nodes } { incr niter 1 } {
      set bnd [ MD++_Get bNds($niter) ]
      set bid [ MD++_Get bTags($niter) ]
      set xb  [ MD++_Get bXs($niter) ]
      set yb  [ MD++_Get bYs($niter) ]
      set zb  [ MD++_Get bZs($niter) ]

      set modeII 0
      set perturb 0.1
      if {$modeII ==0} {
        if { $bid == 4 || $bid == 5 } {		
          set PI 3.1415926
          set newR [expr $yb + $perturb*cos($PI/$Lx0 * $xb)]
          set oldR [MD++_Get R([expr $bnd*3+1])]
          # cos perturb to get determined mode 1 upward/downward
	  set barlength 5
	  set pert [expr $perturb*cos($PI/$barlength * $xb)]
	  MD++ R([expr $bnd*3+1])=[expr $yb+$perturb*cos($PI/$barlength*$xb)]
        } 
      }
      if { $bid == 0 || $bid == 1 } {
         MD++ fixed($bnd) = 1	    
      }
    }
    if { $modeII == 0 || $modeII == 1 } {
      MD++ RHtoS
    }

    if { $modeII == 0 } { 
        set epsilon -0.05
    }
    if { $modeII == 1 } {
        set epsilon -0.05
    }
    set Lx_fix  [expr $Lx0*(1+$epsilon) ] 
    MD++ H_11 = $Lx_fix

# Create 3-d snapband in one direction using eigenstrain formulation
  } elseif { $pbid == 6 } {
    for { set niter 0 } { $niter < $n_bdy_nodes } { incr niter 1 } {
      set bnd [ MD++_Get bNds($niter) ]
      set bid [ MD++_Get bTags($niter) ]
      set xb  [ MD++_Get bXs($niter) ]
      set yb  [ MD++_Get bYs($niter) ]
      set zb  [ MD++_Get bZs($niter) ]
      
# perturb to get two modes
      set modeII 0
      set perturb 0.1
      if { $bid == 0 || $bid == 1 } {
         MD++ fixed($bnd) = 1	    
      }
    }
    if { $modeII == 0 || $modeII == 1 } {
      MD++ RHtoS
    }
# Create 3-d snapband in either x or y direction using eigenstrain formulation
# Do not fix the whole x=xmin and x=max domain, but only fix the corners of that surfaces.
# use snapband executable instead of fem
  } elseif { $pbid == 7 } {
    set modeII 0
    set barlength 8
    set PI 3.1415926
    set perturb 2
    set NP [ MD++_Get "NP" ] 
    # MD++ incnfile = "x-curve.cn" readcn
    # MD++ incnfile = "y-curve-large-152.cn" readcn
    # MD++ incnfile = "xonlycurve.cn" readcn

    set x_eigen_strain 1
    set y_eigen_strain 1.38
puts "I am here 5"
    
    MD++ x_eigen_strain = $x_eigen_strain
    MD++ y_eigen_strain = $y_eigen_strain
    # they are scaled coordinates
    MD++ x_eigen_zbound_min = -10
    MD++ x_eigen_zbound_max =  0.0025
    MD++ y_eigen_zbound_min = -10
    MD++ y_eigen_zbound_max = -0.0025

puts "I am here 6"

    MD++  conj_ftol = 1e-2 conj_itmax = 260 conj_fevalmax = 330000
    MD++ conj_fixbox = 1 
     
    if { 1 } {
puts "I am here 7"
#      relax_fixbox
puts "I am here 7.1"

      if { $x_eigen_strain == 1 } { 
         MD++ finalcnfile = "y-curve-1.cn" writecn
      } elseif { $y_eigen_strain == 1 } { 
         MD++ finalcnfile = "x-curve-1.cn" writecn
      } 
    } else { 
      if { $x_eigen_strain == 1 } { 
         MD++ incnfile = "y-curve-1.cn" readcn
      } elseif { $y_eigen_strain == 1 } { 
         MD++ incnfile = "x-curve-1.cn" readcn
      } 
    } 
puts "I am here 8"
     
    MD++ x_eigen_strain = 0.72
    MD++ y_eigen_strain = 1.38

    MD++  conj_ftol = 1e-2 conj_itmax = 2600 conj_fevalmax = 33000

  MD++ conj_fevalmax = 1000 relax quit
  MD++ eval 
    MD++ finalcnfile = "curve-2.cn" writecn

    set perturb_2_curved_band 1

    if { $perturb_2_curved_band == 1 } {

    if { 0 } {
      for { set iter 0 } { $iter < $NP } { incr iter 1 } {
        set x [ MD++_GetVector R $iter x ] 
        set y [ MD++_GetVector R $iter y ]
        set z [ MD++_GetVector R $iter z ]
        MD++ R([expr $iter*3+2])=[expr $z+$perturb*cos($PI/$barlength*($x)) ]
      }
    } 
    if { 0 } {
      for { set niter 0 } { $niter < $NP } { incr niter 1 } {
        set x [ MD++_GetVector R $niter x ] 
        set y [ MD++_GetVector R $niter y ]
        set z [ MD++_GetVector R $niter z ]
        MD++ R([expr $niter*3+2]) = [expr $z - 2*cos($PI/3.0 * $y)]
      }
    }
    }
    if { $modeII == 0 || $modeII == 1 } {
      MD++ RHtoS
    }
  } elseif { $pbid == 8 } {
# 2d islamic cubes. 	  
    if { 0 } {
      for { set niter 0 } { $niter < $n_bdy_nodes } { incr niter 1 } {
        set bnd [ MD++_Get bNds($niter) ]
        set bid [ MD++_Get bTags($niter) ]
        if { $meshid == 14 || $meshid == 16 } {
          if { $bnd == 12 || $bnd == 11   } {     
            # MD++ fixed($bnd) = 1	    
          }
        }
        if { $meshid == 15 } {
          if { $bid == 0 || $bid == 2 } {
            puts "bnd = $bnd is fixed"
            MD++ fixed($bnd) = 1
          } 
        }              
      }
      
      set epsilon  0
      set Lx_fix [expr $Lx0 * (1+$epsilon)]
      set Ly_fix [expr $Ly0 * (1+$epsilon)]
      MD++ H_11 = $Lx_fix
      MD++ H_22 = $Ly_fix

    } else {
      MD++ saveH	       
      #set initfile ../fem-0-$pbid-$meshid/NEBinit-1.cn
      set initfile ../../../scripts/work/fem/fe-meshes/M15/membrane2d-islam-element.cn	       
      MD++ incnfile = $initfile readcn		     
      MD++ SHtoR
      MD++ restoreH
      MD++ plot
    }
#2d islamic fancy pattern
  } elseif { $pbid == 9 } {

    MD++ saveH	       
    set initfile ../fem-0-$pbid-$meshid/NEBinit-3.cn
    MD++ incnfile = $initfile readcn		     
    MD++ restoreH
    MD++ plot          
    MD++ sleep

  } elseif { $pbid == 10 } {
# 2d islamic fancy pattern, apply force 
    for { set niter 0 } { $niter < $n_bdy_nodes } { incr niter 1 } {
      set bnd [ MD++_Get bNds($niter) ]
      set bid [ MD++_Get bTags($niter) ]	
      if { $bid == 0 } {     
          MD++ fixed($bnd) = 1
          MD++ input = 1 setfixedatomsgroup		    
          MD++ freeallatoms
      }
      if { $bid == 2 } {
          MD++ fixed($bnd) = 1
          MD++ input = 3 setfixedatomsgroup
          MD++ freeallatoms
      }             
      if { $bid == 1 } {
          MD++ fixed($bnd) = 1
          MD++ input = 2 setfixedatomsgroup
          MD++ freeallatoms
      }             
      if { $bid == 3 } {
          MD++ fixed($bnd) = 1
          MD++ input = 4 setfixedatomsgroup
          MD++ freeallatoms
      }             
    }
    set magnitude 0.0
    MD++ extforce = \[ 4 -$magnitude 0 0  0 0 0  $magnitude 0 0  0 0.0 0 \]
  } elseif { $pbid == 11 } {
#1: apply force
    MD++ fixed(0)  = 1
    MD++ input = 1 setfixedatomsgroup 
    MD++ freeallatoms
    
    MD++ fixed(14) = 1
    MD++ input = 2 setfixedatomsgroup
    MD++ freeallatoms
    
    MD++ fixed(51) = 1
    MD++ input = 3 setfixedatomsgroup
    MD++ freeallatoms
    
    MD++ fixed(59) = 1
    MD++ input = 4 setfixedatomsgroup
    MD++ freeallatoms
    
    set magnitude 0.1
    MD++ extforce = \[ 4 -$magnitude -$magnitude 0   $magnitude -$magnitude 0     -$magnitude $magnitude 0     $magnitude  $magnitude 0 \]

    #for 2d truss elements. 
  } elseif { $pbid == 12 } {     	  

    MD++ fixed(11) = 1
    MD++ input = 1 setfixedatomsgroup
    MD++ freeallatoms
    
    MD++ fixed(0) = 1
#   MD++ fixed(3) = 1
    for { set niter 0 } { $niter <= $n_bdy_nodes } { incr niter 1 } {
      set bnd [ MD++_Get bNds($niter) ]
      set bid [ MD++_Get bTags($niter) ]
      set xb  [ MD++_Get bXs($niter) ]
      set yb  [ MD++_Get bYs($niter) ]
      set zb  [ MD++_Get bZs($niter) ]
	    
      if { $bid == 0 } {	    
# MD++ fixed($bnd) = 1
        puts "$bnd is fixed"
      }
    } 
    MD++ RHtoS
  } else {
    puts "pbid value not valid"
  }

#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
  if { 0 } {
    set epsilon 0
    for { set niter 0 } { $niter < 10000000000 } { incr niter 1 } {
      #set epsilon [expr $epsilon + 0.00002]
      #set Lx_fix  [expr $Lx0*(1+$epsilon) ]
      set Ly_fix  [expr $Ly0*(1+$epsilon*0.8) ]
      #MD++ H_11 = $Lx_fix     
       MD++ H_22 = $Ly_fix
      set magnitude [expr 0.5]
      #MD++ extforce = \[ 4 -$magnitude 0 0 0 -$magnitude 0 $magnitude 0 0  0 $magnitude 0 \]
      MD++ extforce = \[ 4 -$magnitude 0 0 0 0 0 $magnitude 0 0 0 0 0 \]
      set temp [expr $niter % 50 ]
      if { $temp == 0 } {
         MD++ finalcnfile = NEBinit${niter}.cn writecn
         MD++ finalcnfile = NEBinit${niter}.cfg  writeatomeyecfg
      }
      MD++ relax
      MD++ plot
    }
  }
  if { 0 } {
    set magnitude 0
    for { set niter 0 } { $niter < 100000 } { incr niter 1 } {
      set magnitude [expr $magnitude + 0.002 ]
      MD++ extforce = \[ 4 -$magnitude 0 0 0 0 0 $magnitude 0 0 0 0 0 \]
      relax_fixbox
      MD++ plot
     }
  }
  relax_fixbox
  
  MD++ plot
  MD++ writeall = 1
  MD++ finalcnfile = "NEBinit-2.cn" writecn
  MD++ finalcnfile = "NEBinit-2.cfg"  writeatomeyecfg
  MD++ sleepseconds = 600
  MD++ sleep
  exitmd 

# same as status 4. neb2d. used for islamic pattern. 
} elseif { $status == 1 } {
  MD++ RtoRref
  
  set initNEBfile ../../scripts/fem/fe-meshes/M15/membrane2d-islam-element.cn		
  #set initNEBfile ../fem-0-$pbid-$meshid/NEBinit-1.cn
  #set finalNEBfile ../fem-0-$pbid-$meshid/NEBinit-2.cn 
  set finalNEBfile ../../scripts/fem/fe-meshes/M16/membrane2d-islam-element.cn
  
  MD++ incnfile = $initNEBfile readcn saveH
  MD++ incnfile = $finalNEBfile readcn restoreH SHtoR setconfig2
  MD++ incnfile = $initNEBfile readcn setconfig1
  
  MD++ fixallatoms  constrain_fixedatoms freeallatoms
  
  MD++ {  chainlength = 24 allocchain   totalsteps = 500000
#500
#chainlength = 20 allocchain   totalsteps = 500
    timestep = 0.00001 printfreq = 500
    initRchain
#incnfile = neb.chain.500 readRchain

    nebspec = [ 0  #0:interpolate surrounding atoms, 1:relax surrounding atoms
		1  #redistribution frequency
	    ]
    nebrelax
    finalcnfile = neb.chain.500 writeRchain
  }
} elseif { $status == 4 } {
  MD++ RtoRref
  	
  set initNEBfile ../fem-0-$pbid-$meshid/NEBinit-1.cn
  set finalNEBfile ../fem-0-$pbid-$meshid/NEBinit-2.cn 
  
  MD++ incnfile = $initNEBfile readcn saveH
  MD++ incnfile = $finalNEBfile readcn restoreH SHtoR setconfig2
  MD++ incnfile = $initNEBfile readcn setconfig1
  
  MD++ fixallatoms  constrain_fixedatoms freeallatoms
  
  if {1} {
    for { set niter 0 } { $niter < $n_bdy_nodes } { incr niter 1 } {
      set bnd [ MD++_Get bNds($niter) ]
      set bid [ MD++_Get bTags($niter) ]
      if { $bid == 0 || $bid == 1 } {
        puts "bnd = $bnd is fixed"
        MD++ fixed($bnd) = 1	    
      }
    }
  }
  
  MD++ {  chainlength = 24 allocchain   totalsteps = 5
#500
#chainlength = 20 allocchain   totalsteps = 500
    timestep = 0.00001 printfreq = 2
    initRchain
#incnfile = neb.chain.500 readRchain
    nebspec = [ 0  #0:interpolate surrounding atoms, 1:relax surrounding atoms
		1  #redistribution frequency
	    ]
    nebrelax
    finalcnfile = neb.chain.500 writeRchain
  }
} elseif { $status == 5 } {
  MD++ RtoRref
  
  set ncpu [ MD++_Get numDomains ]
  if { $ncpu < 1 } { set ncpu 1 }
    puts "ncup = $ncpu"
    
    MD++ Slave_chdir
    #set initNEBfile ../fem-0-$pbid-$meshid/NEBinit.cn
    #set finalNEBfile ../fem-0-$pbid-$meshid/NEBfinal.cn 
    set initNEBfile ../fem-0-$pbid-$meshid/stateA.cn
    set finalNEBfile ../fem-0-$pbid-$meshid/stateB.cn 
    MD++ incnfile = $initNEBfile readcn saveH
    MD++ incnfile = $finalNEBfile readcn restoreH SHtoR setconfig2
    MD++ incnfile = $initNEBfile readcn setconfig1

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
    MD++ timestep = 0.00001 printfreq = 1
	    
    MD++ nebspec = \[ 0 1 0 1 0 0 1 \] totalsteps = 10000 equilsteps = 10000
    MD++ Broadcast_Atoms
    MD++ Broadcast_FEM_Param
    MD++ eval_parallel
    MD++ fixallatoms  constrain_fixedatoms freeallatoms
    if {0} {
      for { set niter 0 } { $niter < $n_bdy_nodes } { incr niter 1 } {
        set bnd [ MD++_Get bNds($niter) ]
        set bid [ MD++_Get bTags($niter) ]
        if { $bid == 0 || $bid == 1 } {
          puts "bnd = $bnd is fixed"
          MD++ fixed($bnd) = 1	    
        }
      }
      MD++ Broadcast_Atoms
    }

    MD++ allocchain_parallel initRchain_parallel

#   if { $exist > 0 } {
#     MD++ incnfile = stringrelax.chain.cn readRchain_parallel
#   }
    MD++ stringrelax_parallel_2
#   MD++ stringrelax_parallel
    MD++ finalcnfile = stringrelax.chain.cn writeRchain_parallel
    if { $ncpu == 1 }   {
       MD++ quit
    } else {
       MD++ quit_all
    }
			    
# Start from an earlier chain of states. 
} elseif { $status == 6 } {
  MD++ RtoRref
  
  puts "meshfolder = $meshfolder"
  set initNEBfile ../fem-0-$pbid-$meshid/NEBinit-1.cn
  set finalNEBfile ../fem-0-$pbid-$meshid/NEBinit-3.cn 
  
  MD++ incnfile = $initNEBfile readcn saveH
  MD++ incnfile = $finalNEBfile readcn restoreH SHtoR setconfig2
  MD++ incnfile = $initNEBfile readcn setconfig1
  
  MD++ fixallatoms  constrain_fixedatoms freeallatoms

  for { set niter 0 } { $niter < $n_bdy_nodes } { incr niter 1 } {
    set bnd [ MD++_Get bNds($niter) ]
    set bid [ MD++_Get bTags($niter) ]
    if { $bid == 0 || $bid == 1 } {
      puts "bnd = $bnd is fixed"
      MD++ fixed($bnd) = 1	    
    }
  }

  MD++ readrchainspec = 48
  MD++ {  chainlength = 48 allocchain   totalsteps = 1500
    timestep = 0.00001 printfreq = 2
    incnfile = ../../scripts/fem/fe-meshes/M11/membrane2d_neb.chain.500 readRchain
  }
  MD++ { 
    nebspec = [ 0 1] 
    nebrelax
    finalcnfile = neb.chain.1000 writeRchain
  }
#same as status = 4, neb for 3d.
} elseif { $status == 7 } {
  MD++ RtoRref
  #set initNEBfile ../fem-0-$pbid-$meshid/NEBinit-1.cn
  #set finalNEBfile ../fem-0-$pbid-$meshid/NEBinit-2.cn 
  
  set initNEBfile ../fem-0-$pbid-$meshid/stateA.cn
  set finalNEBfile ../fem-0-$pbid-$meshid/stateB.cn 
  
  MD++ incnfile = $initNEBfile readcn saveH
  MD++ incnfile = $finalNEBfile readcn restoreH SHtoR setconfig2
  MD++ incnfile = $initNEBfile readcn setconfig1
  
  MD++ fixallatoms  constrain_fixedatoms freeallatoms
  
  if {0} {
    for { set niter 0 } { $niter < $n_bdy_nodes } { incr niter 1 } {
      set bnd [ MD++_Get bNds($niter) ]
      set bid [ MD++_Get bTags($niter) ]
      if { $bid == 0 || $bid == 1 } {
        puts "bnd = $bnd is fixed"
        MD++ fixed($bnd) = 1	    
      }
    }
  }

  MD++ {  chainlength = 48 allocchain   totalsteps = 500000
#500
#chainlength = 20 allocchain   totalsteps = 500
    timestep = 0.00001 printfreq = 2
    initRchain
#incnfile = neb.chain.500 readRchain
    nebspec = [ 0  #0:interpolate surrounding atoms, 1:relax surrounding atoms
		    1  #redistribution frequency
	      ]
    nebrelax
    finalcnfile = neb.chain.500000 writeRchain
  }

#same as status = 7, neb for 3d, pbid = 7, mesh13. .
} elseif { $status == 8 } {
  MD++ RtoRref
  #set initNEBfile ../fem-0-$pbid-$meshid/NEBinit-1.cn
  #set finalNEBfile ../fem-0-$pbid-$meshid/NEBinit-2.cn 
  
  set initNEBfile ../fem-0-$pbid-$meshid/stateA.cn
  set finalNEBfile ../fem-0-$pbid-$meshid/stateB.cn 
  
  MD++ incnfile = $initNEBfile readcn saveH
  MD++ incnfile = $finalNEBfile readcn restoreH SHtoR setconfig2
  MD++ incnfile = $initNEBfile readcn setconfig1
  
  MD++ fixallatoms  constrain_fixedatoms freeallatoms
  
  for { set niter 0 } { $niter < $n_bdy_nodes } { incr niter 1 } {
      set bnd [ MD++_Get bNds($niter) ]
      set bid [ MD++_Get bTags($niter) ]
      set xb  [ MD++_Get bXs($niter) ]
      set yb  [ MD++_Get bYs($niter) ]
      set zb  [ MD++_Get bZs($niter) ]
  
      set xmax 8
      set ymax 8
      set zmax 0.125
      set zloc [expr abs(abs($zb) - $zmax)]
      set yloc [expr abs(abs($yb) - $ymax)]
      set xloc [expr abs(abs($xb) - $xmax)]
  
      if { $xloc < 0.001 && $yloc < 0.001 & $zloc < 0.001 } {
  	    puts "bnd = $bnd is fixed"
  	    #MD++ fixed($bnd) = 1	    
  	}
  }
	
  MD++ {  chainlength = 24 allocchain   totalsteps = 50000
#500
#chainlength = 20 allocchain   totalsteps = 500
    timestep = 0.00001 printfreq = 2

    initRchain
#incnfile = neb.chain.500 readRchain
    nebspec = [ 0  #0:interpolate surrounding atoms, 1:relax surrounding atoms
		1  #redistribution frequency
	    ]
    nebrelax
    finalcnfile = neb.chain.500 writeRchain
  }

# Start from an earlier chain of states. same as status = 6, specifically for islamic 2d.
} elseif { $status == 9 } {
  MD++ RtoRref
  
  puts "meshfolder = $meshfolder"
  set initNEBfile ../fem-0-$pbid-$meshid/NEBinit-1.cn
  set finalNEBfile ../fem-0-$pbid-$meshid/NEBinit-2.cn 
  
  MD++ incnfile = $initNEBfile readcn saveH
  MD++ incnfile = $finalNEBfile readcn restoreH SHtoR setconfig2
  MD++ incnfile = $initNEBfile readcn setconfig1
  
  MD++ fixallatoms  constrain_fixedatoms freeallatoms
  MD++ readrchainspec = 20
  MD++ {  chainlength = 20 allocchain   totalsteps = 1500
       timestep = 0.00001 printfreq = 2
       incnfile = ../fem-0-10-15/NEBinit-neb.chain.500 readRchain
  }
      
  MD++ { 
       nebspec = [ 0 1] 
       nebrelax
       finalcnfile = neb.chain.1500 writeRchain
  }
} elseif { $status == 19 } {

  MD++ setnolog
  MD++ RtoRref
  set initNEBfile ../fem-0-$pbid-$meshid/NEBinit.cn
  set finalNEBfile ../fem-0-$pbid-$meshid/NEBfinal.cn 
  MD++ incnfile = $initNEBfile readcn saveH
  MD++ incnfile = $finalNEBfile readcn restoreH SHtoR setconfig2
  MD++ incnfile = $initNEBfile readcn setconfig1
  MD++ fixallatoms  constrain_fixedatoms freeallatoms
  
  for { set niter 0 } { $niter < $n_bdy_nodes } { incr niter 1 } {
    set bnd [ MD++_Get bNds($niter) ]
    set bid [ MD++_Get bTags($niter) ]
    if { $bid == 0 || $bid == 2 } {
      puts "bnd = $bnd is fixed"
      MD++ fixed($bnd) = 1	    
    }
  }

  MD++ {  chainlength = 24 allocchain   totalsteps = 200
#chainlength = 20 allocchain   totalsteps = 500
    timestep = 0.00001 printfreq = 2
    incnfile = "../fem-4-4-7/neb.chain.500" readRchain 
  }
  set chain_no 0
  set total_no 23
  MD++ chainlength = [expr $total_no+1]
  MD++ input = $chain_no  copyRchaintoCN # RHtoS
  #MD++ clearR0 refreshnnlist eval quit

  setup_window
  openwindow
  MD++ zipfiles = 0
  for { set iter 0 } { $iter <= $total_no } { incr iter 1 } {
    # set chain number
    set chain_no $iter
    MD++ input = $chain_no  copyRchaintoCN #RHtoS
    # MD++ NCS = 8 eval calcentralsymmetry 
    MD++ plot
    MD++ finalcnfile = "chain_no_${chain_no}.cn" writecn
    MD++ finalcnfile = "chain_no_${chain_no}.cfg" writeatomeyecfg
    MD++ sleepseconds = 1  sleep
  }
  MD++ sleepseconds = 100  sleep
  exitmd
} elseif { $status == 20 } {
  MD++ setnolog
  set chain_no 0
  set total_no 23
  MD++ chainlength = $total_no
#  MD++ incnfile = "../fem-0-4-7/NEBinit.cn" readcn  saveH
#  MD++ incnfile = "../fem-5-4-7/stringrelax.chain.cn" 
  MD++ incnfile = "../fem-0-$pbid-$meshid/NEBinit.cn" readcn  saveH
  MD++ incnfile = "../fem-5-$pbid-$meshid/stringrelax.chain.cn" 
  MD++ input = $chain_no  readRchain_parallel_toCN  RHtoS
  #MD++ clearR0 refreshnnlist eval quit
  setup_window
  openwindow
  MD++ zipfiles = 0
  for { set iter 0 } { $iter <= $total_no } { incr iter 1 } {
  # set chain number
      set chain_no $iter
      MD++ input = $chain_no  readRchain_parallel_toCN  RHtoS
      # MD++ NCS = 8 eval calcentralsymmetry 
      MD++ plot
      MD++ finalcnfile = "chain_no_${chain_no}.cn" writecn
      MD++ finalcnfile = "chain_no_${chain_no}.cfg" writeatomeyecfg
      MD++ sleepseconds = 1  sleep
  }
  exitmd
} elseif { $status == 21 } {
  MD++ setnolog
  set chain_no 9
  MD++ RtoRref
  MD++ incnfile = "chain_no_${chain_no}.cn" readcn
  for { set niter 0 } { $niter < $n_bdy_nodes } { incr niter 1 } {
    set bnd [ MD++_Get bNds($niter) ]
    set bid [ MD++_Get bTags($niter) ]
    if { $bid == 0 || $bid == 2 } {
      puts "bnd = $bnd is fixed"
      MD++ fixed($bnd) = 1	    
    }
  }
  
  setup_window
  openwindow
  MD++ eval
  relax_fixbox
  MD++ plot
  relax_fixbox
  MD++ sleepseconds = 100  sleep
  exitmd

} elseif { $status == 22 } {
  MD++ setnolog
  set chain_no 1
  MD++ RtoRref
  MD++ incnfile = "../fem-view-5-8/chain_no_${chain_no}.cn" readcn
  setup_window
  openwindow
  MD++ eval
  MD++ sleep
  exitmd

# visualization of single cpu nebrelax.
} elseif { $status == 23 } {
  MD++ setnolog
  MD++ incnfile = "../29July16-curve-paths-5/fem-0-$pbid-$meshid/NEBinit-1.cn" readcn saveH setconfig1
  MD++ incnfile = "../29July16-curve-paths-5/fem-0-$pbid-$meshid/NEBinit-2.cn" readcn restoreH SHtoR setconfig2
  if {0} {
    setup_window
    MD++ NCS = 8 eval calcentralsymmetry
    MD++ incnfile = "" readcn restoreH SHtoR
    MD++ finalcnfile = "" writeatomeyecfg
    openwindow
    MD++ sleep
    exitmd
  }
  set total_no 48
  MD++ zipfiles = 0
  for { set iter 0 } { $iter <= $total_no } { incr iter 1 } {
#   MD++ incnfile ="../fem-7-$pbid-$meshid/neb.chain.500" readRchain
    MD++ incnfile ="../29July16-curve-paths-5/fem-6-4-8/neb.chain.1000" readRchain
     # set chain number
    set chain_no $iter
    MD++ input = $iter  copyRchaintoCN eval
    if { $iter == 0 } {
      openwindow
    }
    MD++ eval plot 
    MD++ finalcnfile = "chain_no_${chain_no}.cn"  writecn
    MD++ finalcnfile = "chain_no_${chain_no}.cfg" writeatomeyecfg
    MD++ sleepseconds = 1  sleep
  }
  MD++ sleepseconds = 100  sleep
  exitmd

# visualization of single cpu nebrelax. 3d snapband. islamic meshid=15,pbid = 10; status =8. pbid = 7. meshid = 13. 
} elseif { $status == 24 } {
  MD++ setnolog
  MD++ incnfile = "../fem-0-$pbid-$meshid/stateA.cn" readcn saveH setconfig1
  openwindow
  MD++ plot 
  MD++ sleepseconds = 1 sleep
  MD++ incnfile = "../fem-0-$pbid-$meshid/stateB.cn" readcn restoreH SHtoR setconfig2
  openwindow 
  set total_no 24
  MD++ zipfiles = 0
  for { set iter 0 } { $iter <= $total_no } { incr iter 1 } {
    MD++ incnfile ="../fem-8-$pbid-$meshid/neb.chain.500" readRchain
     # set chain number
    set chain_no $iter

    MD++ input = $iter  copyRchaintoCN eval
    #MD++ input = \[ 1 1 2 \] changeH_keepR
    #MD++ input = \[ 2 2 2 \] changeH_keepR
    #MD++ input = \[ 3 3 2 \] changeH_keepR
    
    if { $iter == 0 } {
    	openwindow
    }
    MD++ eval plot 
    MD++ finalcnfile = "chain_no_${chain_no}.cn"  writecn
    MD++ finalcnfile = "chain_no_${chain_no}.cfg" writeatomeyecfg
    MD++ sleepseconds = 2  sleep
  }
  MD++  sleepseconds = 100 sleep
  exitmd
} else {        
  puts "unknown status = $status"
  exitmd 
} 
