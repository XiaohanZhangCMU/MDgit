# -*-shell-script-*-
# compute energy barrier of homogeneous dislocation nucleation in BCC W
# relax all other stress components (except yz) to zero
#
# make fs build=R SYS=mc2
# sw_mc2_mpich scripts/disl_nuc_hetero.tcl 0 0 001 1     -generates 0.01 0.02 0.03 strain perfect states A
# sw_mc2_mpich scripts/disl_nuc_hetero.tcl 1 0.030 001 1 -generates 0.03 state B.
# mpirun -np $ncpu sw_mc2_mpich scripts/disl_nuc_hetero.tcl 4 0.030 001 1 ----string_relax_parallel
source "$::env(MDPLUS_DIR)/scripts/Examples/Tcl/startup.tcl"

#*******************************************
# Definition of procedures
#*******************************************
proc initmd { n } {
#MD++ setnolog
MD++ setoverwrite
MD++ dirname = runs/frankMEP/
#max number neigbors
set myname [ MD++_Get "myname" ]
readmeam-according-to-myname $myname
#max number neigbors
exec cp ../../scripts/work/frankMEP/disl_nuc_hetero.tcl .
exec cp ../../scripts/work/frankMEP/global_search.py .
MD++ NNM = 200
}

proc readpot { } { MD++ {
#--------------------------------------------
#Read in potential file
#
potfile = $::env(MDPLUS_DIR)/potentials/w_pot readpot
} }

#------------------------------------------------------------
proc readmeam-lammps { } { 
#Read in MEAM potential (Baskes format)
MD++ meamfile = "~/Planet/Libs/MD++UMB.svn3/potentials/MEAMDATA/meamf"
#MD++ meafile = "~/Planet/LIbs/MD++UMB.svn3/potentials/MEAMDATA/AuSi2nn.meam" 
MD++ nspecies = 2  element0 = "Siz" element1 = "Ge" 
MD++ {
rcut = 4.5  readMEAM
NNM = 300
} }

# make sure the coordinate is right hand sided.

proc make_perfect_crystal { nx ny nz } {
    MD++ crystalstructure = diamond-cubic latticeconst =  5.4309529817532409 #(A) for Si
    MD++ latticesize = \[  1 1 0  $nx  -1 1 0  $ny  0 0 1  $nz \]
    MD++ makecrystal #finalcnfile = perf.cn writecn #eval
}

#0: Si. 1 : Ge.
proc set_all_atoms_species { id } {
    MD++ fixedallatoms
    MD++ input = $id setfixedatomsspecies
    MD++ freeallatoms
}

proc make_dislocation_loop { flag epsilon } { 
    set store 1

    # xiaohan
    # create structure B'
    # 1) perfect dislocation on shuffle plane {001}
    #    b = a/2[110]
    # 2) partial dislocation on glide plane, {001}
    #    b = a/2[011], b1=a/6[121], b2= a/6[-1,11]
    # H11 gives dimension in xx
    set nu 0.28
    set a 5.4309529817532409

    set bx  [expr 0.5/sqrt(2)]
    set by  [expr 0.5/sqrt(2)]
    set bz  [expr -0.5 ] 

#     set bx  0.5
#     set by  0.5
#     set bz  -0.7072

    set Lx [MD++_Get H_11]
    set Ly [MD++_Get H_22]
    set Lz [MD++_Get H_33]

    set espratio [expr 1 ]
    if { $flag == 001 } {
       #first generate the perfect lattice and visualize with ovito
       #pick one atoms above slip plane, one below. both on top face
       #put their real cooardinate here, 


       if { $epsilon <= 0.030 } {
	       set w [expr $Ly * 0.495]
       } elseif { $epsilon == 0.0500 } {
	       set w [expr $Ly * 0.4875]
       } elseif { $epsilon == 0.0520 } {
	       set w [expr $Ly * 0.486]
       } elseif { $epsilon == 0.0540 } {
	       set w [expr $Ly * 0.486]
       } elseif { $epsilon == 0.0560 } {
               set w [expr $Ly * 0.485 ] 
       } elseif { $epsilon == 0.0580 } {
               set w [expr $Ly * 0.483 ] 
       } elseif { $epsilon == 0.0600 } {
               set w [expr $Ly * 0.482 ]
       } elseif { $epsilon == 0.0620 } {
              set w [expr $Ly * 0.481 ] 
       } elseif { $epsilon == 0.0630 } {
              set w [expr $Ly * 0.480 ] 
       } elseif { $epsilon == 0.0640 } {
              set w [expr $Ly * 0.479 ] 
       } elseif { $epsilon == 0.0700 } {
               set w [expr $Ly * 0.481 ]
       } elseif { $epsilon == 0.0800 } {
               set w [expr $Ly * 0.479 ]
       } else {
              set w [expr $Ly * 0.419]
       }

       set w [expr $w * 0.65]
       set h [expr $w * $espratio]
       set angle 0.9553

#       set RDy  0
#       set RDx  0.7072
       set RDx 0.70749
       set RDz -1

       #calculate the middle point of the loop, which should sit on the shuffle/glide plane
       set P0y [expr  0.0  ]
       set P0x [expr  0.5*(-0.2898-0.2702)  * $Lx]
       set P0z [expr  0.5*(0.2525 +0.2611)  * $Lz]
       
       set P1y $P0y
       set P1x [expr $P0x + $RDx * $h]
       set P1z [expr $P0z + $RDz * $h]
       
       set y1 [expr $P0y - $w]
       set x1 [expr $P0x ]
       set z1 [expr $P0z ]

       set y2 [expr $P0y + $w]
       set x2 [expr $P0x ]
       set z2 [expr $P0z ]

       set y3 [expr $P1y + 0.55*$w]
       set x3 [expr $P1x ]
       set z3 [expr $P1z ]

       set y4 [expr $P1y - 0.55*$w]
       set x4 [expr $P1x ]
       set z4 [expr $P1z ]
    } else {
       puts "unknown flag = $flag, must be 001 or 110"
       return
    } 

    puts "x1 = $x1"
    puts "y1 = $y1"
    puts "z1 = $z1"

    puts "x2 = $x2"
    puts "y2 = $y2"
    puts "z2 = $z2"

    puts "x3 = $x3"
    puts "y3 = $y3"
    puts "z3 = $z3"

    puts "x4 = $x4"
    puts "y4 = $y4"
    puts "z4 = $z4"

    puts "w = $w"
    puts "h = $h"

    puts "Lx = $Lx"
    puts "Ly = $Ly"
    puts "Lz = $Lz"
    
    MD++ input = \[ 1 $store $nu $a $bx $by $bz 4 $x1 $y1 $z1 $x2 $y2 $z2 $x3 $y3 $z3 $x4 $y4 $z4 \] 
    MD++ makedislpolygon
}

proc readmeam-according-to-myname { myname  } {
 if { [ string match "*baskes*" $myname ] } {
  puts "readmeam-baskes"
  readmeam-baskes
 } elseif { [ string match "*lammps*" $myname ] } {
  puts "readmeam-lammps"
  readmeam-lammps
 } elseif { [ string match "*meam*" $myname ] } {
  puts "readmeam"
  readmeam
 } elseif { [ string match "*eam*" $myname ] } {
  puts "readeam"
  readeam
 } else {
  puts "not an eam potential, not reading any files"
 }
}


#--------------------------------------------
proc relax_fixbox { } { MD++ {
# Conjugate-Gradient relaxation
conj_ftol = 2e-6 conj_itmax = 1000 conj_fevalmax = 1000
conj_fixbox = 1
relax
} }
#end of proc relax_fixbox

#--------------------------------------------
proc relax_freebox { } { MD++ {
# Conjugate-Gradient relaxation
conj_ftol = 1e-4 conj_itmax = 1000 conj_fevalmax = 1000
conj_fixbox = 0
conj_fixboxvec = [ 0 1 1
                   1 1 1
                   1 1 0 ]
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
atomradius = [1.0 0.78] bondradius = 0.3 bondlength = 2.8285 #for Si
win_width=600 win_height=600
#atomradius = 0.9 bondradius = 0.3 bondlength = 0 #2.8285 #for Si
atomcolor = orange highlightcolor = purple  bondcolor = red
fixatomcolor = red backgroundcolor = gray70
#atomcolor = lightgrey highlightcolor = purple  bondcolor = darkgrey
plot_color_windows = [ 0
                 -5   -4   0
                 -5.5 -5   8   #color00 = red
                 -6  -5.5  1
               -6.55 -6    2
               -6.72 -6.55 3
                 -6.5 -6   6   #color06 = gray80
                 -7.5 -6.5 4
                ]

#xiaohan, if you want to plot the dislocations, uncomment the following

plot_color_axis = 2  NCS = 4
plot_color_windows = [ 0
                       0.6 9.4   1  
                       9.4  10   5
                       10 20   6
                       20 50   8
                       0  0.6  4
                     ]

#plot_limits = [ 1 -10 10 -10 10.0 0.05 10 ]


#
#xiaohan

#plot_color_windows = 5
#plot_color_windows = 0 
plot_atom_info = 1 # reduced coordinates of atoms
#plot_limits = [ 1 -10 10  -10 10 0 10 ]
#plot_atom_info = 2 # real coordinates of atoms
#plot_atom_info = 3 # energy of atoms
#plot_highlight = [ 0 0 1 2 3 4 5 6 7 8 9 ]
plotfreq = 10
#

rotateangles = [ -0 90 0 1.2 ]

#rotateangles = [ 0 0 0 1.7 ]
#rotateangles = [ 0 -90 0 1.7 ]
#openwin alloccolors rotate saverot plot
#plot_color_axis = 0 input = [ -8 -3 10] GnuPlotHistogram
#plot_color_axis = 2 input = [ 0.6 50 50 ] GnuPlotHistogram
} }

proc openwindow { } { 
setup_window
MD++ openwin alloccolors rotate saverot eval plot
}

#--------------------------------------------
proc exitmd { } { MD++ quit }
#end of proc exitmd
#--------------------------------------------

#--------------------------------------------
proc setup_md { } { MD++ {     
T_OBJ = 300 #Kelvin #add by xiaohan

equilsteps = 0  totalsteps = 5000 timestep = 0.0001 # (ps)
atommass = 28.0855 # (g/mol)
DOUBLE_T = 1
saveprop = 1 savepropfreq = 100 openpropfile #run
savecn = 1 savecnfreq = 10000 openintercnfile
plotfreq = 100 printfreq = 100
#ensemble_type = "NPH" integrator_type = "Gear6" implementation_type = 0
ensemble_type = "NVE" integrator_type = "VVerlet" implementation_type = 0
vt2 = 1e28  #1e28 2e28 5e28
wallmass = 2e3     # atommass * NP = 14380
boxdamp = 1e-3     # optimal damping for 216 atoms and wallmass 1e-3
saveH # Use current H as reference (H0), needed for specifying stress
fixboxvec = [ 0 0 1
              1 0 1
              0 0 0 ]
output_fmt = "curstep EPOT KATOM Tinst HELM HELMP TSTRESS_xx TSTRESS_yy TSTRESS_zz H_11 H_22 H_33" 
} }
#end of proc setup_md

proc myRand {min max} { return [expr int(rand()*($max-$min+1)) + $min] }

#*******************************************
# Main program starts here
#*******************************************
# status 0:
#        1:
#        2:
#
# read in status from command line argument
if { $argc <= 0 } {
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
 set flag 0
} elseif { $argc > 2 } {
 set flag [lindex $argv 2]
}
puts "flag = $flag"

if { $argc <= 3 } {
 set opt1 0
} elseif { $argc > 3 } {
 set opt1 [lindex $argv 3]
}
puts "opt1 = $opt1"

if { $argc <= 4 } {
 set opt2 0
} elseif { $argc > 4 } {
 set opt2 [lindex $argv 4]
}
puts "opt2 = $opt2"

if { $argc <= 5 } {
 set opt3 0
} elseif { $argc > 5 } {
 set opt3 [lindex $argv 5]
}
puts "opt3 = $opt3"

if { $argc <= 6 } {
 set opt4 0
} elseif { $argc > 6 } {
 set opt4 [lindex $argv 6]
}
puts "opt4 = $opt4"

if { $argc <= 7 } {
 set opt5 0
} elseif { $argc > 7 } {
 set opt5 [lindex $argv 7]
}
puts "opt5 = $opt5"

if { $argc <= 8 } {
 set opt6 0
} elseif { $argc > 8 } {
 set opt6 [lindex $argv 8]
}
puts "opt6 = $opt6"

if { $argc <= 9 } {
 set opt7 0
} elseif { $argc > 9 } {
 set opt7 [lindex $argv 9]
}
puts "opt7 = $opt7"

if { $argc <= 10 } {
 set opt8 0
} elseif { $argc > 10 } {
 set opt8 [lindex $argv 10]
}
puts "opt8 = $opt8"

if { $argc <= 11 } {
 set opt9 0
} elseif { $argc > 11 } {
 set opt9 [lindex $argv 11]
}
puts "opt9 = $opt9"


if { $status == 0 } {
  MD++ setnolog
  initmd $n
  readpot

  make_perfect_crystal 12 12 16
#  make_perfect_crystal 42 21 42

  MD++ finalcnfile = "w-perf.cn" writecn
  MD++ finalcnfile = "w-perf.cfg" writeatomeyecfg

  MD++ eval saveH conj_fixboxvec =  \[ 0 1 0  0 0 0   1 1 0 \]
  MD++ conj_fixbox = 1 conj_ftol = 2e-2 conj_itmax = 1000 conj_fevalmax = 2000
  MD++ relax finalcnfile = "0K_perfect.cn" writecn
  MD++ eval

  
  set dighole 0
  set make_free_surface 1  
  set surface_reconstruct 1
  
  if { $dighole } {
     MD++ input = \[ 1 -10 -0.25  -10 10   0.25 10 \] fixatoms_by_position
     MD++ removefixedatoms
     MD++ freeallatoms      
  }
   #end of dig hole

  if {$make_free_surface} {
    MD++ vacuumratio = [expr 1.0-1.0/(1.0 + 0.2)]
    MD++ conj_fixbox = 0 
    if { $flag == 001 } {
       # create 001 surface

#MD++ conj_ftol = 2e-6 conj_itmax = 1000 conj_fevalmax = 2

       MD++ input = \[ 3 3 0.2 \] changeH_keepR relax finalcnfile = "0K_surf${flag}.cn" writecn
       set strain_data_file "strain_surf_001.dat"
       set stress_data_file "stress_surf_001.dat"
       MD++ { conj_fixboxvec = [ 0 1 1
                                 1 1 1
                                 1 1 1 ] }
    }
    MD++ conj_fixbox = 1
  }
   # end of make free surface


#  move top and bot layer atoms to form dimers
#  Top layer
#  0.2440 -- 0.2381
if { $surface_reconstruct } {
 set ny  [ MD++_Get latticesize(7) ]
 for { set i 0 } { $i < $ny } { incr i 1 } {
    set ymin [ expr -0.5006+1.0/$ny*$i ]
    set ymax [ expr -0.5006+1.0/$ny*($i+0.5) ]
    MD++ input = \[ 1 -10 10 $ymin $ymax 0.403   10 \]  fixatoms_by_position
    MD++ input = \[ 1 -10 10 $ymin $ymax  -10  -0.416 \]  fixatoms_by_position
 }
 MD++ input = 1  setfixedatomsgroup  freeallatoms

 for { set i 0 } { $i < $ny } { incr i 1 } {
    set ymin [ expr -0.5006+1.0/$ny*($i+0.5) ]
    set ymax [ expr -0.5006+1.0/$ny*($i+1) ]
    MD++ input = \[ 1 -10 10 $ymin $ymax 0.403   10 \]  fixatoms_by_position
    MD++ input = \[ 1 -10 10 $ymin $ymax -10 -0.416 \]  fixatoms_by_position
 }
 MD++ input = 2  setfixedatomsgroup  freeallatoms

 MD++ input = \[ 1  0  0.8 0  1 \] movegroup  
 MD++ input = \[ 1  0 -0.8 0  2 \] movegroup  
}
 #end of surf reconstruct

  relax_fixbox
  MD++ finalcnfile = "0K_0.0_relaxed_surf${flag}.cn" writecn
  MD++ finalcnfile = "0K_0.0_relaxed_surf${flag}.cfg" writeatomeyecfg

  set epsilon 0.06
  set H11_0 [ MD++_Get H_11 ]; set H22_0 [ MD++_Get H_22 ]; set H33_0 [ MD++_Get H_33 ]
  MD++ saveH
  set H11_fix [ expr $H11_0*(1.0-$epsilon) ]
  MD++ H_11 = $H11_fix
  MD++ conj_ftol = 1e-4 conj_itmax = 3800 conj_fevalmax = 6000
  MD++ conj_fixbox = 1  relax
  MD++ SHtoR
  MD++ eval
  set fp [ open "EPOT_1.dat" a+ ]; puts $fp [ MD++_Get EPOT ]; close $fp
  MD++ finalcnfile = "0K_${epsilon}_strained_${n}.cn" writecn
  MD++ finalcnfile = "0K_${epsilon}_strained_${n}.cfg" writeatomeyecfg

  exitmd

# Read cn file from global_search.py, a treched configuration
# Relax with a little help of applied strain to close up the trench
# Evaluate energy and write (append) to trench_energy.out at line ($n)
# 

} elseif { $status == 1 } {
  MD++ setnolog
  initmd $status
  set sliceid 1

  set A $flag
  set B $opt1
  set C $opt2
  set D $opt3

  set x0 $opt4
  set y0 $opt5
  set z0 $opt6

  set icdx $opt7
  set icdy $opt8
  set icdz $opt9
    
  MD++ incnfile = "trench_${n}.cn" readcn
  MD++ relax
  set strain_data_file "strain_surf_001.dat"
  set stress_data_file "stress_surf_001.dat"

  set fp [ open "nucleus-$n.dat" r ]
  set tol [expr 0.05]

  set file_data [ read $fp ]
  close $fp
  set data [split $file_data "\n"];

  set heterogeneous 1
  if { $heterogeneous == 1 } { 
    set NP [ MD++_Get NP ] 
    set Nnuc [ llength $data ]
    set Nnuc [ expr $Nnuc -1 ]
    puts "Nnuc = $Nnuc"
    #puts "data = $data"

    for { set i 0 } { $i < $NP } { incr i 1 } {    
      set xi [ MD++_Get SR([expr 3*$i])   ]
      set yi [ MD++_Get SR([expr 3*$i+1]) ]
      set zi [ MD++_Get SR([expr 3*$i+2]) ]
      if { [expr abs($xi-$x0)] < $icdx && [expr abs($yi-$y0)] < $icdy && [expr abs($zi-$z0)] < $icdz } { 
        set count [expr 0 ]
        foreach line $data {
            if { $count < $Nnuc } { 
              set linedata [ split $line " " ]
              lassign $linedata x y z
              set dx [ expr $xi - $x ]
              set dy [ expr $yi - $y ]
              set dz [ expr $zi - $z ]
              if { [expr abs($dx) ] < $tol && [expr abs($dz) ] < $tol && [expr abs($dz)]<$tol } {
                 set tmp [expr $A * $xi + $B * $yi + $C * $zi + $D ]
                 if { $tmp >0 } {
                   MD++ input = \[ 1 $i \] fixatoms_by_ID
                 }
              }
            }
            set count [expr $count +1]
        }  
      }
    }
    MD++ input = 1 setfixedatomsgroup freeallatoms

    for { set i 0 } { $i < $NP } { incr i 1 } {    
      set xi [ MD++_Get SR([expr 3*$i])   ]
      set yi [ MD++_Get SR([expr 3*$i+1]) ]
      set zi [ MD++_Get SR([expr 3*$i+2]) ]
      if { [expr abs($xi-$x0)] < $icdx && [expr abs($yi-$y0)] < $icdy && [expr abs($zi-$z0)] < $icdz } { 
        set count [expr 0 ]
        foreach line $data {
            if { $count < $Nnuc } { 
              set linedata [ split $line " " ]
              lassign $linedata x y z
              set dx [ expr $xi - $x ]
              set dy [ expr $yi - $y ]
              set dz [ expr $zi - $z ]
              if { [expr abs($dx) ] < $tol && [expr abs($dz) ] < $tol && [expr abs($dz)]<$tol } {
                 set tmp [expr $A * $xi + $B * $yi + $C * $zi + $D ]
                 if { $tmp < 0 } {
                   MD++ input = \[ 1 $i \] fixatoms_by_ID
                 }
              }
            }
            set count [expr $count +1]
        }  
      }
    }

    MD++ input = 2 setfixedatomsgroup freeallatoms
    
    set mag  1
    set magx [ expr $A*$mag ]
    set magy [ expr $B*$mag ]
    set magz [ expr $C*$mag ]

    MD++ input = \[ 1 -$magx -$magy -$magz 1 \] movegroup
    MD++ input = \[ 1  $magx  $magy  $magz 2 \] movegroup
  }

  set epsilon 0.06
  set H11_0 [ MD++_Get H_11 ]; set H22_0 [ MD++_Get H_22 ]; set H33_0 [ MD++_Get H_33 ]
  MD++ saveH
  set H11_fix [ expr $H11_0*(1.0-$epsilon) ]
  MD++ H_11 = $H11_fix
  MD++ conj_ftol = 1e-4 conj_itmax = 3800 conj_fevalmax = 6000
  MD++ conj_fixbox = 1  relax
  MD++ SHtoR

  MD++ finalcnfile = "0K_${epsilon}_relaxed_${n}.cn" writecn
  MD++ finalcnfile = "0K_${epsilon}_relaxed_${n}.cfg" writeatomeyecfg
  MD++ eval 
  set fp [ open "EPOT_2.dat" a+ ]; puts $fp [ MD++_Get EPOT ]; close $fp

  exitmd

} elseif { $status == 2 } { 

  initmd $status
  MD++ setnolog
  #MD++ incnfile = "0K_0.0_relaxed_surf${flag}.cn" readcn
  make_perfect_crystal 12 12 16
  MD++ relax eval

  if { 1 } { 
    set EPOT_0 [MD++_Get EPOT]
    puts "EPOT_0 = $EPOT_0"
    set NP [MD++_Get NP]
    set n [myRand 0 $NP]
    MD++ input=\[ 1 $n \] fixatoms_by_ID
    MD++ removefixedatoms
    MD++ freeallatoms
    MD++ relax
    MD++ eval 
    set EPOT_1 [MD++_Get EPOT]
    puts "EPOT_1 = $EPOT_1"
  }


  MD++ sleep 
  MD++ freeallatoms

  set fp [ open "numvoids.txt" ]
  set file_data [ read $fp ]
  close $fp
  set data [ split $file_data "\n" ]
  set Maxiters [ llength $data ] 
  set voidFrac 0.01
  set sum 0
  for { set iter 0 } { $iter < $Maxiters } { incr iter 1 } {
    set voidFrac [expr $voidFrac * 2 ]
    set NP [MD++_Get NP]
    #set N_voids [expr int($NP*$voidFrac/100.0)]
    set N_voids [ lindex $data $iter ]
    puts "N_voids = $N_voids"
#    setup_window
#    openwindow

    MD++ eval
    set epot0 [ MD++_Get EPOT ]
    set fp [ open "VOIDS_EPOT_$iter.dat" a+ ]; puts $fp [ MD++_Get EPOT ]; 
    close $fp

    for { set x 0 } { $x < $N_voids } { incr x 1 } {
      set NP [MD++_Get NP]
      set n [myRand 0 $NP]
      MD++ input=\[ 1 $n \] fixatoms_by_ID
      MD++ removefixedatoms
    }

    MD++ freeallatoms
    MD++ eval
    set epot1 [ MD++_Get EPOT ]
    relax_fixbox
    MD++ eval
    set epot2 [ MD++_Get EPOT ]

    set tmp [ expr ($epot1-$epot0)/$N_voids ]
    set sum [ expr $sum + ($epot1-$epot0)/$N_voids ]
    puts "energy difference / # of voids = $tmp"

    set fp [ open "VOIDS_EPOT.dat" a+ ]; 
    puts $fp "${N_voids} $epot0 $epot1 $epot2 $sum"
    close $fp
  }
  set sum [ expr $sum/$Maxiters ]
  puts "energy of a void is $sum"
 
} elseif { $status == 3 } { 


  MD++ setnolog
  initmd $status
  set sliceid 1

  set A $flag
  set B $opt1
  set C $opt2
  set D $opt3

  set x0 $opt4
  set y0 $opt5
  set z0 $opt6

  set icdx $opt7
  set icdy $opt8
  set icdz $opt9
    
  MD++ incnfile = "trench_${n}.cn" readcn
  MD++ relax
  MD++ eval
  set fp [ open "EPOT_1.dat" a+ ]; puts $fp [ MD++_Get EPOT ]; close $fp

  #openwindow
  
  set strain_data_file "strain_surf_001.dat"
  set stress_data_file "stress_surf_001.dat"

  set fp [ open "nucleus-$n.dat" r ]
  set tol [expr 0.05]

  set file_data [ read $fp ]
  close $fp
  set data [split $file_data "\n"];

  set heterogeneous 1
  if { $heterogeneous == 1 } { 
  set NP [ MD++_Get NP ] 
  set Nnuc [ llength $data ]
  set Nnuc [ expr $Nnuc -1 ]
  puts "Nnuc = $Nnuc"
  #puts "data = $data"

  for { set i 0 } { $i < $NP } { incr i 1 } {    
    set xi [ MD++_Get SR([expr 3*$i])   ]
    set yi [ MD++_Get SR([expr 3*$i+1]) ]
    set zi [ MD++_Get SR([expr 3*$i+2]) ]
    if { [expr abs($xi-$x0)] < $icdx && [expr abs($yi-$y0)] < $icdy && [expr abs($zi-$z0)] < $icdz } { 
      set count [expr 0 ]
      foreach line $data {
	  if { $count < $Nnuc } { 
            set linedata [ split $line " " ]
	    lassign $linedata x y z
#	    puts "n = $n, x = $x, y = $y, z = $z, count = $count"
            set dx [ expr $xi - $x ]
            set dy [ expr $yi - $y ]
            set dz [ expr $zi - $z ]
            if { [expr abs($dx) ] < $tol && [expr abs($dz) ] < $tol && [expr abs($dz)]<$tol } {
               set tmp [expr $A * $xi + $B * $yi + $C * $zi + $D ]
               if { $tmp >0 } {
                 MD++ input = \[ 1 $i \] fixatoms_by_ID
               }
            }
	  }
	  set count [expr $count +1]
      }  
    }
  }
  MD++ input = 1 setfixedatomsgroup freeallatoms

  for { set i 0 } { $i < $NP } { incr i 1 } {    
    set xi [ MD++_Get SR([expr 3*$i])   ]
    set yi [ MD++_Get SR([expr 3*$i+1]) ]
    set zi [ MD++_Get SR([expr 3*$i+2]) ]
    if { [expr abs($xi-$x0)] < $icdx && [expr abs($yi-$y0)] < $icdy && [expr abs($zi-$z0)] < $icdz } { 
      set count [expr 0 ]
      foreach line $data {
	  if { $count < $Nnuc } { 
            set linedata [ split $line " " ]
	    lassign $linedata x y z
#	    puts "n = $n, x = $x, y = $y, z = $z, count = $count"
            set dx [ expr $xi - $x ]
            set dy [ expr $yi - $y ]
            set dz [ expr $zi - $z ]
            if { [expr abs($dx) ] < $tol && [expr abs($dz) ] < $tol && [expr abs($dz)]<$tol } {
               set tmp [expr $A * $xi + $B * $yi + $C * $zi + $D ]
               if { $tmp < 0 } {
                 MD++ input = \[ 1 $i \] fixatoms_by_ID
               }
            }
	  }
	  set count [expr $count +1]
      }  
    }
  }

  MD++ input = 2 setfixedatomsgroup freeallatoms
  
  set mag  1
  set magx [ expr $A*$mag ]
  set magy [ expr $B*$mag ]
  set magz [ expr $C*$mag ]
#  set magx [ expr $mag ]
#  set magy [ expr 0 ]
#  set magz [ expr 0 ]

  MD++ input = \[ 1 -$magx -$magy -$magz 1 \] movegroup
  MD++ input = \[ 1  $magx  $magy  $magz 2 \] movegroup
  }

  #MD++ plot 
  MD++ finalcnfile = "0K_unstrained_${n}.cn" writecn
  MD++ finalcnfile = "0K_unstrained_${n}.cfg" writeatomeyecfg

  set H11_0 [ MD++_Get H_11 ] ; set H22_0 [ MD++_Get H_22 ] ; set H33_0 [ MD++_Get H_33 ]

  MD++ saveH

  set maxitereps 0
  set epsilon 0.05
  for { set itereps 0 } { $itereps <= $maxitereps } { incr itereps 1 } {

      set H11_fix [ expr $H11_0*(1.0-$epsilon) ]
      MD++ H_11 = $H11_fix
      MD++ conj_ftol = 1e-4 conj_itmax = 3800 conj_fevalmax = 6000
      MD++ conj_fixbox = 1  relax
      MD++ SHtoR

      MD++ finalcnfile = "0K_${epsilon}_strained_${n}.cn" writecn
      MD++ finalcnfile = "0K_${epsilon}_strained_${n}.cfg" writeatomeyecfg

      MD++ eval 
      set fp [ open "EPOT_2.dat" a+ ]; puts $fp [ MD++_Get EPOT ]; close $fp

      set H11_fix [ MD++_Get H_11 ] ; set H22_new [ MD++_Get H_22 ] ; set H33_new [ MD++_Get H_33 ] 
      set H31_cur [ MD++_Get H_31 ] ; set H12_cur [ MD++_Get H_12 ] ; set H23_new [ MD++_Get H_23 ]

      set sig_xx [ MD++_Get TSTRESSinMPa_xx ] ; set sig_yy [ MD++_Get TSTRESSinMPa_yy ]
      set sig_zz [ MD++_Get TSTRESSinMPa_zz ] ; set sig_xy [ MD++_Get TSTRESSinMPa_xy ]
      set sig_xz [ MD++_Get TSTRESSinMPa_xz ] ; set sig_yz [ MD++_Get TSTRESSinMPa_yz ]
      set e_xx [ expr ($H11_fix-$H11_0)/$H11_0 ] ; set e_yy [ expr ($H22_new-$H22_0)/$H22_0 ]
      set e_zz [ expr ($H33_new-$H33_0)/$H33_0 ] ; set e_yz [ expr ($H23_new)/$H33_0 ]
      set e_xy [ expr ($H12_cur)/$H22_0 ] ;        set e_xz [ expr ($H31_cur)/$H33_0 ]

      set fp [ open $strain_data_file a+ ] ; puts $fp "$e_xx $e_yy $e_zz $e_xy $e_xz $e_yz" ; close $fp
      set fp [ open $stress_data_file a+ ] ; puts $fp "$sig_xx $sig_yy $sig_zz $sig_xy $sig_xz $sig_yz" ; close $fp

      set epsilon [ expr $epsilon+0.002 ]
      set epsilon [format "%.4f" $epsilon]

  }
  exitmd

} elseif { $status == 4 } { 

  MD++ setnolog
  initmd $status
    
  MD++ incnfile = "newcfg.cn" readcn
  MD++ relax eval
  set H11_0 [ MD++_Get H_11 ]; set H22_0 [ MD++_Get H_22 ]; set H33_0 [ MD++_Get H_33 ]
  MD++ saveH
  set H11_fix [ expr $H11_0*(1.0-$epsilon) ]
  MD++ H_11 = $H11_fix
  MD++ conj_ftol = 1e-4 conj_itmax = 3800 conj_fevalmax = 6000
  MD++ conj_fixbox = 1  relax
  MD++ SHtoR

  MD++ finalcnfile = "0K_${epsilon}_strained_${n}.cn" writecn
  MD++ finalcnfile = "0K_${epsilon}_strained_${n}.cfg" writeatomeyecfg
  MD++ eval 
  set fp [ open "EPOT_2.dat" a+ ]; puts $fp [ MD++_Get EPOT ]; close $fp
  exitmd

} elseif { $status == 100 } {
  initmd $status
  MD++ incnfile = "${n}.cn" readcn
  MD++ finalcnfile = "${n}.cfg" writeatomeyecfg
} elseif { $status == 1000 } {
  MD++ incnfile = "runs/frankMEP/0K_0.05_relaxed_18.cn" readcn
  MD++ relax

#  for { set itereps 0 } { $itereps <= 99 } { incr itereps 1 } {
#     MD++ incnfile = "runs/frankMEP/data_saved_2/0K_-0.040_strained_${itereps}.cn" readcn
#     MD++ eval
#     set fp [ open "EPOT.dat" a+ ]; puts $fp [ MD++_Get EPOT ]; close $fp
#  }
 MD++ eval
  MD++ incnfile = "runs/frankMEP/data_saved_3/0K_0.05_relaxed_0.cn" readcn
MD++ eval
} else {
 puts "unknown status = $status"
 exitmd 

} 

