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
  conj_ftol = 5e-4 conj_itmax = 1000 conj_fevalmax = 130000
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

proc Fatal { mesg } {
  puts $mesg
  MD++ quit
}

#*******************************************
# Main program starts here
#*******************************************
# status 0:
#        1:
#        2:
# eqnType 0: beam. 1: snapband. 2: islam
if { $argc <= 3 } {
  Fatal "Need to specify status, pbid, meshid and eqnType"
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

puts "elefile = $elefile"
puts "filename = $filename"

MD++ elements_file = $elefile  read_elements
MD++ fem_coeff_file = ${filename}-C3D8-Q8.dat read_fem_coeff
#MD++ fem_coeff_file = ${filename}-CPE4-Q4.dat read_fem_coeff
MD++ incnfile = $fecnfile  readcn
MD++ fem_bdy_nodes_file = $bdyfile  read_bdy_nodes

# meshid
# 1: 2d beam. 
# 2: islam, fancy, expanded
# 3: snap band 3d.
if { $meshid==1 || $meshid == 2 } { 
  MD++ map23 = \[ 0 1 \] 
} elseif { $meshid==3 } { 
  MD++ map23 = \[-1 -1 \]
} else {
  Fatal "mesh id is not valid."
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
  } else { 
    setup_window
    openwindow
  }

  if { $pbid == 1 } {     	  
    MD++ writeall = 1
    set PI 3.1415926
    set barlength 2.5
    set epsilon -0.035
    for { set sAB 0 } { $sAB < 3 } { incr sAB 1 } {  
      set sgn [ expr $sAB * 2 - 1 ] 
      for { set niter 0 } { $niter < $n_bdy_nodes } { incr niter 1 } {
        set bnd [ MD++_Get bNds($niter) ]
        set bid [ MD++_Get bTags($niter) ]
        set xb  [ MD++_Get bXs($niter) ]
        set yb  [ MD++_Get bYs($niter) ]
        # perturb to get two modes
        set perturb [expr 0.01 * $sgn] 
        if { $bid == 0 || $bid == 2 } {
          MD++ fixed($bnd) = 1	    
        }
        if { $bid == 1 || $bid == 3 } {		
          MD++ R([expr $bnd*3+1])=[expr $yb+$perturb*cos($PI/$barlength*$xb)]
	  if { $sAB == 2 } { 
            MD++ R([expr $bnd*3+1])=[expr $yb+$perturb*0.1*sin($PI/$barlength*$xb)]
          } 
        } 
      }
      MD++ RHtoS
      set Lx_fix  [expr $Lx0*(1+$epsilon) ] 
      MD++ H_11 = $Lx_fix
      MD++ plot
      relax_fixbox
      MD++ finalcnfile = "NEBinit-$sAB.cn" writecn
      MD++ finalcnfile = "NEBinit-$sAB.cfg" writeatomeyecfg
      MD++ incnfile = $fecnfile  readcn 
      MD++ SHtoR
    }
    MD++ sleep
    exitmd

  } elseif { $pbid == 2 } { 

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

    set magnitude 0
    for { set niter 0 } { $niter < 10000 } { incr niter 1 } {
      set magnitude [expr $magnitude + 0.002 ]
      MD++ extforce = \[ 4 -$magnitude 0 0 0 0 0 $magnitude 0 0 0 0 0 \]
      relax_fixbox
      MD++ plot
    }
    relax_fixbox
    
    MD++ plot
    MD++ writeall = 1
    MD++ finalcnfile = "NEBinit-2.cn" writecn
    MD++ finalcnfile = "NEBinit-2.cfg"  writeatomeyecfg

# Create 3-d snapband in either x or y direction using eigenstrain formulation
# Do not fix the whole x=xmin and x=max domain, but only fix the corners of that surfaces.
  } elseif { $pbid == 3 } {
    set modeII 0
    set barlength 8
    set PI 3.1415926
    set perturb 2
    set NP [ MD++_Get "NP" ] 
    set x_curve  0
    set y_curve  1
    set ic_preparing 0

    if { $x_curve == 1 } { 
      MD++ x_eigen_strain = 0.73
      MD++ y_eigen_strain = 1.29
    } else {
      MD++ x_eigen_strain = 0.83
      MD++ y_eigen_strain = 1.29
    }
    # they are scaled coordinates
    MD++ x_eigen_zbound_min = -10
    MD++ x_eigen_zbound_max =  0
    MD++ y_eigen_zbound_min = -10
    MD++ y_eigen_zbound_max =  0

     
    if { ${ic_preparing} == 1 } {
      MD++  conj_ftol = 1e-2 conj_itmax = 260 conj_fevalmax = 330000
      MD++ conj_fixbox = 1 relax
      if { $x_curve == 1 } { 
         MD++ finalcnfile = "x-curve-1.cn" writecn
      } elseif { $y_curve == 1 } { 
         MD++ finalcnfile = "y-curve-1.cn" writecn
      } 
      exitmd
    } else { 
      if { $x_curve == 1 } { 
         MD++ incnfile = "x-curve-1.cn" readcn
#         MD++ plot sleep 
      } elseif { $y_curve == 1 } { 
         MD++ incnfile = "y-curve-1.cn" readcn
      } 
    } 
    if { $USEGPU == 1 } {
      MD++ cuda_memcpy_all
    }
    MD++ plot
     
    MD++ x_eigen_strain = 0.7323
    MD++ y_eigen_strain = 1.42

#    MD++ relaxation_algorithm =  PRPLUS
    MD++ conj_fixbox = 1 
    if { $x_curve == 1} { 
      MD++  conj_ftol = 1e-5 conj_itmax = 2600 conj_fevalmax = 130000 relax 
    } else {
      MD++  conj_ftol = 1e-5 conj_itmax = 2600 conj_fevalmax = 130000 relax 
    }
    MD++ eval 
    if { $y_curve == 1 } { 
      MD++ finalcnfile = "NEBinit-0.cn" writecn
    } else { 
      MD++ finalcnfile = "NEBinit-1.cn" writecn
    }

    MD++ plot
    MD++ writeall = 1
    MD++ sleepseconds = 600
    MD++ sleep
    exitmd 

  } else {
    puts "pbid value not valid"
  }

} elseif { $status == 1 } {
  	
  if { $pbid == 1 } { 

    for { set sAB 0 } { $sAB < 3 } { incr sAB 1 } {  
      MD++ RtoRref 
      if { $sAB == 0 } {
         set A 0
	 set B 1
      } elseif { $sAB == 1 } { 
         set A 0
	 set B 2
      } elseif { $sAB == 2 } {
         set A 2 
	 set B 1
      } 
      set initNEBfile ../fem-0-$pbid-$meshid-$eqnType/NEBinit-$A.cn
      set finalNEBfile ../fem-0-$pbid-$meshid-$eqnType/NEBinit-$B.cn 

      MD++ incnfile = $initNEBfile readcn  saveH
      MD++ incnfile = $finalNEBfile readcn restoreH SHtoR setconfig2
      MD++ incnfile = $initNEBfile readcn setconfig1
      MD++ fixallatoms  constrain_fixedatoms freeallatoms
      
      for { set niter 0 } { $niter < $n_bdy_nodes } { incr niter 1 } {
        set bnd [ MD++_Get bNds($niter) ]
        set bid [ MD++_Get bTags($niter) ]
        if { $bid == 0 || $bid == 2 } {
          MD++ fixed($bnd) = 1	    
        }
      }

      MD++ {  chainlength = 24 allocchain   totalsteps = 10
        timestep = 0.00001 printfreq = 2
        initRchain
        nebspec = [ 0  #0:interpolate surrounding atoms, 1:relax surrounding atoms
            	1  #redistribution frequency
                ]
        nebrelax
      }
      set ttlstps [ MD++_Get totalsteps ]
      exec mv nebeng.out nebeng_${A}_${B}_${ttlstps}.out
      MD++ finalcnfile = neb_chain_${A}_${B}.${ttlstps} writeRchain
    } 

# same as status 4. neb2d. used for islamic pattern. 
  } elseif { $pbid == 2 } {
    MD++ RtoRref 
    set initNEBfile ../../../scripts/work/fem/fe-meshes/M2/membrane2d-islam-element.cn		
    set finalNEBfile ../../../scripts/work/fem/fe-meshes/M22/membrane2d-islam-element.cn
    MD++ incnfile = $initNEBfile readcn saveH
    MD++ incnfile = $finalNEBfile readcn restoreH SHtoR setconfig2
    MD++ incnfile = $initNEBfile readcn setconfig1
    MD++ fixallatoms  constrain_fixedatoms freeallatoms
    MD++ {  chainlength = 24 allocchain   totalsteps = 500000
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
  } elseif { $pbid == 3 } {
    MD++ RtoRref 
    MD++ x_eigen_strain = 0.7323
    MD++ y_eigen_strain = 1.42
    MD++ x_eigen_zbound_min = -10
    MD++ x_eigen_zbound_max =  0
    MD++ y_eigen_zbound_min = -10
    MD++ y_eigen_zbound_max =  0

    set initNEBfile ../fem-0-$pbid-$meshid-$eqnType/NEBinit-0.cn
#    set finalNEBfile ../fem-0-$pbid-$meshid-$eqnType/NEBinit-1.cn 
    set finalNEBfile ../fem-view-23-$pbid-$meshid-$eqnType/chain_no_23.cn 

    MD++ incnfile = $initNEBfile readcn saveH
    #setup_window
    #openwindow
    #MD++ eval plot sleep 

    MD++ incnfile = $finalNEBfile readcn restoreH SHtoR setconfig2


    MD++ incnfile = $initNEBfile readcn setconfig1
    MD++ eval sleep
    MD++ fixallatoms  constrain_fixedatoms freeallatoms
    if { $USEGPU == 1 } {
      MD++ cuda_memcpy_all
    }
    MD++ {  chainlength = 48 allocchain   totalsteps = 5000
      timestep = 0.0001 printfreq = 100
      initRchain
      #incnfile = neb.chain.500 readRchain
      nebspec = [ 0  #0:interpolate surrounding atoms, 1:relax surrounding atoms
          	1  #redistribution frequency
              ]
      nebrelax
      finalcnfile = neb.chain.5000 writeRchain
    }
    MD++ {  chainlength = 48 allocchain   totalsteps = 50000
      timestep = 0.00001 printfreq = 100
      initRchain
      incnfile = neb.chain.5000 readRchain
      nebspec = [ 0  #0:interpolate surrounding atoms, 1:relax surrounding atoms
          	1  #redistribution frequency
              ]
      nebrelax
      finalcnfile = neb.chain.50000 writeRchain
    }
  }
} elseif { $status == 2 } {

  # First mannually paste two neb_chain_${A}_${B}.totalsteps paths together
  # which are created from status = 1, then perturb the path and nebrelax() 
  set ttlstps 10
  set bfld "../fem-1-$pbid-$meshid-$eqnType/"
  set filein_1 [open "${bfld}neb_chain_0_2.${ttlstps}" r]
  set filein_2 [open "${bfld}neb_chain_2_1.${ttlstps}" r]
  set data_1 [ split [read $filein_1 ] "\n" ]
  set data_2 [ split [read $filein_2 ] "\n" ]
  close $filein_1
  close $filein_2

  set fileout [ open init_path.${ttlstps} w ]
  # read in the straight initial path. perturb each node
  set nchain_1 [ lindex $data_1 0 ] 
  set nchain_2 [ lindex $data_2 0 ] 
  # remove the middle connecting state
  set nchain   [ expr $nchain_1+$nchain_2 ]
  set nnodes   [ lindex $data_1 1 ]   
  puts "nchain =  $nchain"
  puts "nnodes =  $nnodes"
  puts $fileout $nchain
  puts $fileout $nnodes


  for { set i 0 } { $i < [expr $nnodes] } { incr i 1 } {
     set nd [ lindex $data_1 [expr $i+2] ]
     puts $fileout $nd
  }

  for { set chain 0 } { $chain <= $nchain_1-1 } { incr chain 1 } {  
    for { set node 0 } { $node < $nnodes } { incr node 1 } {
      # read in rchain 
      set data1 [ lindex $data_1 [expr 2+$nnodes+$chain*$nnodes+$node] ]
      if { $data1 == "" } { 
        break
      }  
      set xold [ lindex $data1 0 ]
      set yold [ lindex $data1 1 ]
      set zold [ lindex $data1 2 ] 
      puts $fileout "$xold $yold $zold"
    }
  }

  for { set chain 0 } { $chain <= $nchain_2 } { incr chain 1 } {  
    for { set node 0 } { $node < $nnodes } { incr node 1 } {
      # read in rchain 
      set data2 [ lindex $data_2 [expr 2+$nnodes+$chain*$nnodes+$node] ]
      if { $data2 == "" } { 
        break
      }  
      set xold [ lindex $data2 0 ]
      set yold [ lindex $data2 1 ]
      set zold [ lindex $data2 2 ] 
      puts $fileout "$xold $yold $zold"
    }
  }

  close $fileout
#  MD++ sleep 

  # Mannually perturb initial path 0 -- 2 -- 1 by adding random modes to it
  MD++ RtoRref
  MD++ saveH

  set PI 3.1415926
  set perturb 0.0005
#  set perturb 0.0
  set fileout [open init_path.${ttlstps} r] 
  set data [ split [read $fileout ] \n]
  close $fileout

  # ntrials > 1 does not work. Need to figure this out 
  set ntrials 1 
  set Nend [ llength $data ]
  for { set itrial 0 } { $itrial < $ntrials } { incr itrial 1 } { 
    set fileout trial_path-${itrial}.chain.${ttlstps}
    set fp [ open $fileout w ]
    # read in the straight initial path. perturb each node
    set nchain [lindex $data 0 ] 
    set nnodes [lindex $data 1 ]   
    puts "nchain =  $nchain"
    puts "nnodes =  $nnodes"
    puts $fp $nchain
    puts $fp $nnodes

    for { set i 0 } { $i < [expr $nnodes] } { incr i 1 } {
       set data1 [ lindex $data [expr $i+2] ]
       puts $fp $data1
    }

    for { set chain 0 } { $chain <= $nchain } { incr chain 1 } {  
      for { set node 0 } { $node < $nnodes } { incr node 1 } {
        # read in rchain 
        set data1 [ lindex $data [expr 2+$nnodes+$chain*$nnodes+$node] ]
        if { $data1 == "" } { 
           break
        }  
        set xref [ MD++_Get Rref([expr $node*3]) ]
        set xold [ lindex $data1 0 ]
        set yold [ lindex $data1 1 ]
        set zold [ lindex $data1 2 ] 
               
        set y $yold
        set barlength 2.5
        if { $chain != 0 && $chain != $nchain } {
           set nmodes 10
           for {set j 0 } { $j < $nmodes } { incr j 1 } {
               set freq [expr int(rand()*3)+1 ]
               set mag [expr rand()*$perturb]
               set y [expr $y+$mag*sin(2*$PI*$freq/$barlength*$xref) ]
           }
        }
	       
        puts $fp "$xold $y $zold"
      }
    }
    close $fp

    # do neb relaxation with initial chain
    MD++ RtoRref
    MD++ restoreH
    MD++ RHtoS

    set initNEBfile ../fem-0-$pbid-$meshid-$eqnType/NEBinit-0.cn
    set finalNEBfile ../fem-0-$pbid-$meshid-$eqnType/NEBinit-1.cn
    MD++ incnfile = $initNEBfile readcn saveH
    MD++ incnfile = $finalNEBfile readcn restoreH SHtoR setconfig2
    MD++ incnfile = $initNEBfile readcn setconfig1

    MD++ fixallatoms constrain_fixedatoms freeallatoms
    for { set niter 0 } { $niter < $n_bdy_nodes } { incr niter 1 } {
        set bnd [ MD++_Get bNds($niter) ]
        set bid [ MD++_Get bTags($niter) ]
        if { $bid == 0 || $bid == 1 } {
             puts "bnd = $bnd is fixed"
             MD++ fixed($bnd) = 1	    
        }
    }
  
    MD++ { chainlength = 48 allocchain timestep = 0.00001 printfreq = 2 }
    MD++ readrchainspec = $nchain
    MD++ incnfile = $fileout  readRchain
    set totalsteps 5000
    set chainout "path-${itrial}-neb.chain.$totalsteps"
    MD++ finalcnfile = $chainout
    MD++ totalsteps = $totalsteps
    MD++ { nebspec = [ 0 1 ]  nebrelax writeRchain }
    file mkdir "TRIAL-PATH-${itrial}"
#   exec mv $fileout nebeng.out  "path-${itrial}-neb.chain.${totalsteps}" gp_rawdata.out strain12.out strain11.out strain22.out energy.out "TRIAL-PATH-${itrial}/"
    exec mv $fileout nebeng.out "path-${itrial}-neb.chain.${totalsteps}" "TRIAL-PATH-${itrial}/"
  }
  MD++ quit

#MPI stringrelax_parallel
} elseif { $status == 3 } {
  MD++ RtoRref
  set ncpu [ MD++_Get numDomains ]
  if { $ncpu < 1 } { set ncpu 1 }
  puts "ncup = $ncpu"
    
  MD++ Slave_chdir
  set initNEBfile ../fem-0-$pbid-$meshid-$eqnType/NEBinit-0.cn
  set finalNEBfile ../fem-0-$pbid-$meshid-$eqnType/NEBinit-1.cn 
  MD++ incnfile = $initNEBfile readcn saveH
  MD++ incnfile = $finalNEBfile readcn restoreH SHtoR setconfig2
  MD++ incnfile = $initNEBfile readcn setconfig1

  if { $pbid == 3 } {
      MD++ x_eigen_strain = 0.7323
      MD++ y_eigen_strain = 1.42
      MD++ x_eigen_zbound_min = -10
      MD++ x_eigen_zbound_max =  0
      MD++ y_eigen_zbound_min = -10
      MD++ y_eigen_zbound_max =  0
  }

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
  MD++ timestep = 0.0001 printfreq = 10

  MD++ nebspec = \[ 0 1 0 1 0 0 1 \] totalsteps = 10000 equilsteps = 10000
  MD++ fixallatoms  constrain_fixedatoms freeallatoms

  MD++ allocchain_parallel initRchain_parallel
  if { $exist > 0 } {
    MD++ incnfile = stringrelax.chain.cn readRchain_parallel
  }

  MD++ Broadcast_Atoms
  MD++ Broadcast_FEM_Param
  MD++ Broadcast_nebsetting
  MD++ eval_parallel
  MD++ stringrelax_parallel_2

  MD++ finalcnfile = stringrelax.chain.cn writeRchain_parallel
  if { $ncpu == 1 }   {
    MD++ quit
  } else {
    MD++ quit_all
  }
# Read in chain_no_${iter} from status 24, then re-do nebrelax
} elseif { $status == 8 } {
  MD++ RtoRref
  set ncpu [ MD++_Get numDomains ]
  if { $ncpu < 1 } { set ncpu 1 }
  puts "ncup = $ncpu"
    
  MD++ Slave_chdir
  set initNEBfile ../fem-view-23-$pbid-$meshid-$eqnType/chain_no_6.cn
  set finalNEBfile ../fem-view-23-$pbid-$meshid-$eqnType/chain_no_10.cn 
  MD++ incnfile = $initNEBfile readcn saveH
  MD++ incnfile = $finalNEBfile readcn restoreH SHtoR setconfig2
  MD++ incnfile = $initNEBfile readcn setconfig1

  if { $pbid == 3 } {
      MD++ x_eigen_strain = 0.7323
      MD++ y_eigen_strain = 1.42
      MD++ x_eigen_zbound_min = -10
      MD++ x_eigen_zbound_max =  0
      MD++ y_eigen_zbound_min = -10
      MD++ y_eigen_zbound_max =  0
  }

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
  MD++ timestep = 0.0001 printfreq = 10

  MD++ nebspec = \[ 0 1 0 1 0 0 1 \] totalsteps = 10000 equilsteps = 10000
  MD++ fixallatoms  constrain_fixedatoms freeallatoms

  MD++ allocchain_parallel initRchain_parallel
  if { $exist > 0 } {
    MD++ incnfile = stringrelax.chain.cn readRchain_parallel
  }

  MD++ Broadcast_Atoms
  MD++ Broadcast_FEM_Param
  MD++ Broadcast_nebsetting
  MD++ eval_parallel
  MD++ stringrelax_parallel_2

  MD++ finalcnfile = stringrelax.chain.cn writeRchain_parallel
  if { $ncpu == 1 }   {
    MD++ quit
  } else {
    MD++ quit_all
  }

} elseif { $status == 9 } {
  MD++ RtoRref 
  MD++ x_eigen_strain = 0.7323
  MD++ y_eigen_strain = 1.42
  MD++ x_eigen_zbound_min = -10
  MD++ x_eigen_zbound_max =  0
  MD++ y_eigen_zbound_min = -10
  MD++ y_eigen_zbound_max =  0

  set initNEBfile ../fem-view-23-$pbid-$meshid-$eqnType/chain_no_29.cn
  set finalNEBfile ../fem-view-23-$pbid-$meshid-$eqnType/chain_no_35.cn

  MD++ incnfile = $initNEBfile readcn saveH
  #setup_window
  #openwindow
  #MD++ eval plot sleep 

  MD++ incnfile = $finalNEBfile readcn restoreH SHtoR setconfig2
  MD++ incnfile = $initNEBfile readcn setconfig1
  MD++ fixallatoms  constrain_fixedatoms freeallatoms
  if { $USEGPU == 1 } {
    MD++ cuda_memcpy_all
  }
  MD++ {  chainlength = 8 allocchain   totalsteps = 50000
    timestep = 0.0001 printfreq = 100
    initRchain
    #incnfile = neb.chain.500 readRchain
    nebspec = [ 0  #0:interpolate surrounding atoms, 1:relax surrounding atoms
        	1  #redistribution frequency
            ]
    nebrelax
    finalcnfile = neb_29_35.chain.500 writeRchain
  }

# Manipulate several paths and put them together as an new initial path
} elseif { $status == 10 } {

  set bfld9 "../fem-9-$pbid-$meshid-$eqnType"
  set bfld1 "../fem-1-$pbid-$meshid-$eqnType"
  set chain_file_1 [open "${bfld9}-s13-e20/neb_13_20.chain.500" r]
  set chain_file_2 [open "${bfld9}-s29-e35/neb_29_35.chain.500" r]
  set chain_file_3 [open "${bfld1}/neb.chain.50000" r]
  set path_1 [ split [read ${chain_file_1} ] "\n" ]
  set path_2 [ split [read ${chain_file_2} ] "\n" ]
  set path_3 [ split [read ${chain_file_3} ] "\n" ]
  close ${chain_file_1}
  close ${chain_file_2}
  close ${chain_file_3}

  #replace subchains of 13-20 and 29-35 of path_3 with path_1 and path_2
  set fileout [ open neb.chain.adjusted w ]
  set nchain [ lindex ${path_3} 0 ]
  set nnodes [ lindex ${path_3} 1 ]
  puts $fileout $nchain
  puts $fileout $nnodes
  
  for { set i 0 } { $i < [expr $nnodes] } { incr i 1 } {
     set nd [ lindex ${path_3} [expr $i+2] ]
     puts $fileout $nd
  }

  set perturb 0.001
  for { set chain 0 } { $chain <= $nchain } { incr chain 1 } {  
    for { set node 0 } { $node < $nnodes } { incr node 1 } {
      if { $chain >= 16 && $chain <= 18 } { 
        set d1 [ lindex ${path_3} [expr 2+$nnodes+(15-0)*$nnodes+$node] ]
        set x [ lindex $d1 0 ]; set y [ lindex $d1 1 ]; set z [ lindex $d1 2 ];
        set x [ expr $x + rand()*$perturb ]
        set y [ expr $y + rand()*$perturb ]
        set z [ expr $z + rand()*$perturb ]
        puts $fileout "$x $y $z"
      } elseif { $chain >= 31 && $chain <= 32 } { 
        set d1 [ lindex ${path_3} [expr 2+$nnodes+(30-0)*$nnodes+$node] ]
        set x [ lindex $d1 0 ]; set y [ lindex $d1 1 ]; set z [ lindex $d1 2 ];
        set x [ expr $x + rand()*$perturb ]
        set y [ expr $y + rand()*$perturb ]
        set z [ expr $z + rand()*$perturb ]
        puts $fileout "$x $y $z"
      } else { 
        set d1 [ lindex ${path_3} [expr 2+$nnodes+($chain-0)*$nnodes+$node] ]
        set x [ lindex $d1 0 ]; set y [ lindex $d1 1 ]; set z [ lindex $d1 2 ];
        puts $fileout "$x $y $z"
      }
      if { $d1 == "" } { 
        break
      }  
    }
  }
  close $fileout

  MD++ RtoRref 
  MD++ x_eigen_strain = 0.7323
  MD++ y_eigen_strain = 1.42
  MD++ x_eigen_zbound_min = -10
  MD++ x_eigen_zbound_max =  0
  MD++ y_eigen_zbound_min = -10
  MD++ y_eigen_zbound_max =  0

  set initNEBfile ../fem-0-$pbid-$meshid-$eqnType/NEBinit-0.cn
  set finalNEBfile ../fem-0-$pbid-$meshid-$eqnType/NEBinit-1.cn 
  MD++ incnfile = $initNEBfile readcn saveH
  MD++ incnfile = $finalNEBfile readcn restoreH SHtoR setconfig2
  MD++ incnfile = $initNEBfile readcn setconfig1
  MD++ fixallatoms  constrain_fixedatoms freeallatoms
  if { $USEGPU == 1 } {
    MD++ cuda_memcpy_all
  }
  MD++ {  chainlength = 48 allocchain   totalsteps = 5000
    timestep = 0.0001 printfreq = 100
    initRchain
    incnfile = neb.chain.adjusted readRchain
    nebspec = [ 0  #0:interpolate surrounding atoms, 1:relax surrounding atoms
        	1  #redistribution frequency
            ]
    nebrelax
    finalcnfile = neb.chain.5000 writeRchain
  }

} elseif { $status == 11} {
  MD++ setnolog
  MD++ RtoRref 
  if { $pbid == 3 } {
    MD++ x_eigen_strain = 0.7323
    MD++ y_eigen_strain = 1.42
    MD++ x_eigen_zbound_min = -10
    MD++ x_eigen_zbound_max =  0
    MD++ y_eigen_zbound_min = -10
    MD++ y_eigen_zbound_max =  0
  }
  set initNEBfile ../fem-0-$pbid-$meshid-$eqnType/NEBinit-0.cn
  set finalNEBfile ../fem-0-$pbid-$meshid-$eqnType/NEBinit-1.cn 
  MD++ incnfile = $initNEBfile readcn saveH
  MD++ incnfile = $finalNEBfile readcn restoreH SHtoR setconfig2
  MD++ incnfile = $initNEBfile readcn setconfig1
  MD++ fixallatoms  constrain_fixedatoms freeallatoms
  setup_window
#  openwindow
  if { $USEGPU == 1 } {
    MD++ cuda_memcpy_all
  }
  MD++ {  chainlength = 48 allocchain   totalsteps = 5000
    timestep = 0.0001 printfreq = 100
    initRchain
    nebspec = [ 0  #0:interpolate surrounding atoms, 1:relax surrounding atoms
          	1  #redistribution frequency
              ]
  }

  set nchain [MD++_Get chainlength]
  set perturb 0.0
  set nnodes [ MD++_Get NP ]
  MD++ incnfile ="../fem-1-$pbid-$meshid-$eqnType/neb.chain.50000" readRchain
  for { set iter 0 } { $iter <= $nchain } { incr iter 1 } {
    if {$iter >=16 && $iter <=18  } { 
      MD++ input = $iter  copyRchaintoCN
      for { set node 0 } { $node < $nnodes } { incr node 1 } {
        set xold [ MD++_Get SR(3*$node)  ]
        set yold [ MD++_Get SR(3*$node+1)  ]
        set zold [ MD++_Get SR(3*$node+2)  ]
        set x [ expr $xold +$perturb*rand()]
        set y [ expr $yold +$perturb*rand()]
        set z [ expr $zold +$perturb*rand()]
        MD++ SR([expr $node*3]) = $x
        MD++ SR([expr $node*3+1]) = $y
        MD++ SR([expr $node*3+2]) = $z
      }
      MD++ SHtoR
      MD++ input = $iter copyCNtoRchain 
      #MD++ conj_fevalmax = 3000 conj_fixbox = 1 relax
      #MD++ sleepseconds = 2 sleep
    } elseif { $iter >=31 && $iter <=32 } {
      MD++ input = $iter  copyRchaintoCN
      for { set node 0 } { $node < $nnodes } { incr node 1 } {
        set xold [ MD++_Get SR(3*$node)  ]
        set yold [ MD++_Get SR(3*$node+1)  ]
        set zold [ MD++_Get SR(3*$node+2)  ]
        set x [ expr $xold +$perturb*rand()]
        set y [ expr $yold +$perturb*rand()]
        set z [ expr $zold +$perturb*rand()]
        MD++ SR([expr $node*3]) = $x
        MD++ SR([expr $node*3+1]) = $y
        MD++ SR([expr $node*3+2]) = $z
      }
      MD++ SHtoR
      MD++ input = $iter  copyCNtoRchain 
      #MD++ conj_fevalmax = 3000 conj_fixbox = 1 relax
      #MD++ sleepseconds = 2 sleep
    } else { 
      MD++ input = $iter  copyRchaintoCN
      MD++ SHtoR
      MD++ eval plot 
      MD++ input = $iter  copyCNtoRchain 
    }
  }
#  MD++ finalcnfile = neb.chain.adjusted writeRchain

  MD++ {
#    incnfile = neb.chain.adjusted readRchain
    nebrelax
    finalcnfile = neb.chain.5000 writeRchain
  }

  exitmd


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
  MD++ incnfile = "../fem-0-$pbid-$meshid-$eqnType/NEBinit.cn" readcn  saveH
  MD++ incnfile = "../fem-5-$pbid-$meshid-$eqnType/stringrelax.chain.cn" 
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
  MD++ RtoRref 
  if { $pbid == 3 } {
    MD++ x_eigen_strain = 0.7323
    MD++ y_eigen_strain = 1.42
    MD++ x_eigen_zbound_min = -10
    MD++ x_eigen_zbound_max =  0
    MD++ y_eigen_zbound_min = -10
    MD++ y_eigen_zbound_max =  0
  }

  MD++ incnfile = "../fem-0-$pbid-$meshid-$eqnType/NEBinit-0.cn" readcn saveH setconfig1
  MD++ incnfile = "../fem-0-$pbid-$meshid-$eqnType/NEBinit-1.cn" readcn restoreH SHtoR setconfig2
  if {0} {
    setup_window
    MD++ NCS = 8 eval calcentralsymmetry
    MD++ incnfile = "" readcn restoreH SHtoR
    MD++ finalcnfile = "" writeatomeyecfg
    MD++ sleep
    exitmd
  }
#  setup_window
#  openwindow
#  MD++ incnfile = "chain_23.cn" readcn
#  MD++ plot sleep 

  set total_no 48
  MD++ zipfiles = 0
  for { set iter 0 } { $iter <= $total_no } { incr iter 1 } {
    MD++ incnfile ="../snap_saved/fem-10-$pbid-$meshid-$eqnType/neb.chain.5000" readRchain
    set chain_no $iter
    MD++ input = $iter  copyRchaintoCN SHtoR eval
    if { $iter == 230 } { 
      MD++ SHtoR

      if { $USEGPU == 1 } {
        MD++ cuda_memcpy_all
      } 

      MD++ conj_ftol = 5e-1 conj_itmax = 1000 conj_fevalmax = 130000
      MD++ relax
      MD++ finalcnfile = "chain_23.cn" writecn
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
  set chain_no 0
  set total_no 24
  MD++ chainlength = $total_no
  MD++ incnfile = "../fem-0-$pbid-$meshid-$eqnType/NEBinit-0.cn" readcn saveH setconfig1
  MD++ incnfile = "../fem-0-$pbid-$meshid-$eqnType/NEBinit-1.cn" readcn restoreH SHtoR setconfig2

  MD++ incnfile = "../fem-3-$pbid-$meshid-$eqnType/stringrelax.chain.cn" 
  MD++ input = $chain_no  readRchain_parallel_toCN  RHtoS
  setup_window
  openwindow
  MD++ zipfiles = 0
  for { set iter 0 } { $iter <= $total_no } { incr iter 1 } {
      set chain_no $iter
      MD++ input = $chain_no  readRchain_parallel_toCN  RHtoS
      MD++ plot
      MD++ finalcnfile = "chain_no_${chain_no}.cn" writecn
      MD++ finalcnfile = "chain_no_${chain_no}.cfg" writeatomeyecfg
      MD++ sleepseconds = 1  sleep
  }
  exitmd

} else {        
  puts "unknown status = $status"
  exitmd 
} 
