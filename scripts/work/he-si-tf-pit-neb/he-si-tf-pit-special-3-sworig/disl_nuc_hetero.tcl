# -*-shell-script-*-
# compute energy barrier of shuffle-glide dislocation nucleation from silicon (001) surface with a pit
#
# make sworig build=R SYS=mc2
# sworig_mc2_mpich scripts/disl_nuc_hetero.tcl 0 0 001 1     -generates strained perfect thin film --> state A.
# sworig_mc2_mpich scripts/disl_nuc_hetero.tcl 1 0.030 001 1 -generates shuffle glide dislcoaiton -->  state B.
# mpirun -np $ncpu sworig_mc2_mpich scripts/disl_nuc_hetero.tcl 4 0.030 001 1 ----string_relax_parallel.
# sworig_mc2_mpich scripts/disl_nuc_hetero.tcl 20 0.030 001 1 ---- visualize relaxed chain.
source "$::env(MDPLUS_DIR)/scripts/Examples/Tcl/startup.tcl"

#*******************************************
# Definition of procedures
#*******************************************
proc initmd { n } {
#MD++ setnolog
MD++ setoverwrite
#relative to the  directory of lauching 
MD++ dirname = runs/he-si-tf-pit-neb/he-si-tf-pit-special-3-sworig
#max number neigbors
MD++ NNM = 200
MD++ nspecies = 1  element0 = "Siz"
}

proc readpot { } { MD++ {
#--------------------------------------------
#Read in potential file
#
potfile = $::env(MDPLUS_DIR)/potentials/w_pot readpot
} }

# make sure the coordinate is right hand sided.

proc make_perfect_crystal { nx ny nz } {
    MD++ crystalstructure = diamond-cubic latticeconst =  5.4309529817532409 #(A) for Si
    MD++ latticesize = \[  1 1 0  $nx  -1 1 0  $ny  0 0 1  $nz \]
    MD++ makecrystal #finalcnfile = perf.cn writecn #eval
}


#make a partial dislocation on one shuffle plane
proc make_glide_dislocation_loop_1 { flag epsilon } { 
    set store 1
    set nu 0.28
    set a 5.4309529817532409


    #set bx  [expr 1.0* sqrt(2.0) /12.0 ]
    #set by  [expr 1.0* sqrt(2.0) / 4.0 ]
    #set bz  [expr 1.0/6.0 ]

     set bx [ expr 1.0*sqrt(2)/6.0 ] 
     set by [ expr 0 ] 
     set bz [ expr 1.0/3.0]

    set Lx [MD++_Get H_11]
    set Ly [MD++_Get H_22]
    set Lz [MD++_Get H_33]

    set espratio [expr 0.7 ]
    if { $flag == 001 } {
       #first generate the perfect lattice and visualize with ovito
       #pick one atoms above slip plane, one below. both on top face
       #put their real cooardinate here, 


       if { $epsilon <= 0.030 } {
	       set w [expr $Ly * 0.495]
       } elseif { $epsilon == 0.0500 } {
	       set w [expr $Ly * 0.486]
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
              set w [expr $Ly * 0.486 ]
       }

       set w [expr $w * 0.65 ]
       set h [expr $w * $espratio]
#       set w [expr $w * 40]
       set angle 0.9553

#       set RDy  0
#       set RDx  0.7072
       set RDx -0.70749
       set RDz -1

       #calculate the middle point of the loop, which should sit on the shuffle/glide plane


       set P0y [expr  0.0 * $Ly ]
       set P0x [expr  0.5*(0.2711+0.2921) * $Lx]
       set P0z [expr  0.5*(0.2615+0.2529) * $Lz ]
       
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






#make another partial dislocation on one shuffle plane
proc make_glide_dislocation_loop_2 { flag epsilon } { 
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


    #set bx  [expr 1.0* sqrt(2.0) /12.0 ]
    #set by  [expr 1.0* sqrt(2.0) / 4.0 ]
    #set bz  [expr 1.0/6.0 ]

     set bx [ expr 1.0*sqrt(2)/6.0 ] 
     set by [ expr 0 ] 
     set bz [ expr 1.0/3.0]

    set Lx [MD++_Get H_11]
    set Ly [MD++_Get H_22]
    set Lz [MD++_Get H_33]

    set espratio [expr 0.7 ]
    if { $flag == 001 } {
       #first generate the perfect lattice and visualize with ovito
       #pick one atoms above slip plane, one below. both on top face
       #put their real cooardinate here, 


       if { $epsilon <= 0.030 } {
	       set w [expr $Ly * 0.495]
       } elseif { $epsilon == 0.0500 } {
	       set w [expr $Ly * 0.486]
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
              set w [expr $Ly * 0.486 ]
       }

       set w [expr $w * 0.65 ]
       set h [expr $w * $espratio]
 #      set w [expr $w * 40]
       set angle 0.9553

#       set RDy  0
#       set RDx  0.7072
       set RDx -0.70749
       set RDz -1

       #calculate the middle point of the loop, which should sit on the shuffle/glide plane


       set P0y [expr  0.0* $Ly  ]
       set P0x [expr  0.5*(0.2507+0.2717) * $Lx]
       set P0z [expr  0.5*(0.2752+0.2686) * $Lz ]
       
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



#make the middle row of atoms between the two shuffle plane to rotate.
proc make_glide_dislocation_loop_3 { flag epsilon } { 
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


    #set bx  [expr 1.0* sqrt(2.0) /12.0 ]
    #set by  [expr 1.0* sqrt(2.0) / 4.0 ]
    #set bz  [expr 1.0/6.0 ]

     set fac 0.1
     set bx [ expr sqrt(2) * $fac ] 
     set by [ expr 0 ] 
     set bz [ expr -1.0 * $fac]

    set Lx [MD++_Get H_11]
    set Ly [MD++_Get H_22]
    set Lz [MD++_Get H_33]

    set espratio [expr 0.7 ]
    if { $flag == 001 } {
       #first generate the perfect lattice and visualize with ovito
       #pick one atoms above slip plane, one below. both on top face
       #put their real cooardinate here, 


       if { $epsilon <= 0.030 } {
	       set w [expr $Ly * 0.495]
       } elseif { $epsilon == 0.0500 } {
	       set w [expr $Ly * 0.486]
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
              set w [expr $Ly * 0.486 ]
       }

       set w [expr $w * 0.65 ]
       set h [expr $w * $espratio]
   #    set w [expr $w * 40]
       set angle 0.9553

#       set RDy  0
#       set RDx  0.7072
       set RDx -0.70749
       set RDz -1

       #calculate the middle point of the loop, which should sit on the shuffle/glide plane


       set P0y [expr  0.0*$Ly  ]
       set P0x [expr  0.25*( 0.2507 + 0.2717 + 0.2711 + 0.2921 ) * $Lx]
       set P0z [expr  0.25*( 0.2752 + 0.2686 + 0.2615 + 0.2529 ) * $Lz ]
       
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

#plot_limits = [ 1 -10 10 -0.25 0.25 0.05 10 ]


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
 set flag 0
} elseif { $argc > 2 } {
 set flag [lindex $argv 2]
}
puts "flag = $flag"

if { $argc <= 3 } {
 set opt 0
} elseif { $argc > 3 } {
 set opt [lindex $argv 3]
}
puts "opt = $opt"


if { $status == 0 } {
  MD++ setnolog
  initmd $n
  readpot

#  make_perfect_crystal 30 15 42
  make_perfect_crystal 12 12  21

  MD++ finalcnfile = "w-perf.cn" writecn
  MD++ finalcnfile = "w-perf.cfg" writeatomeyecfg

  MD++ eval saveH conj_fixboxvec =  \[ 0 1 0  0 0 0   1 1 0 \]
  MD++ conj_fixbox = 1 conj_ftol = 2e-2 conj_itmax = 1000 conj_fevalmax = 2000
  MD++ relax finalcnfile = "0K_perfect.cn" writecn

  MD++ eval
# setup_window
# openwindow

  
  set dighole 0
  set make_free_surface 1  
  set surface_reconstruct 1
  
  if { $dighole } {
     MD++ input = \[ 1 -0.25 0.25  -0.25 0.25   0.4 10 \] fixatoms_by_position
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

       MD++ input = \[ 3 3 0.6 \] changeH_keepR relax  relax finalcnfile = "0K_surf${flag}.cn" writecn
       set strain_data_file "strain_surf_001.dat"
       set stress_data_file "stress_surf_001.dat"
       MD++ { conj_fixboxvec = [ 0 1 1
                                 1 1 1
                                 1 1 1 ] }
    }
    MD++ conj_fixbox = 1
  }
   # end of make free surface


set NP [ MD++_Get "NP" ]
set xholemin [ expr -0.5 + 0.2083 ]
set yholemin [ expr -0.5 + 0.2083 ]
set xholemax [ expr 0.5 - 0.2083 ]
set yholemax [ expr 0.5 - 0.2083 ]
set zholemax [ expr 0.5*5.0/6.0 - 0.14 ]

for { set ip 0 } { $ip < $NP } { incr ip 1 } {
     set sxi [ MD++_GetVector SR $ip x ] 
     set syi [ MD++_GetVector SR $ip y ]
     set szi [ MD++_GetVector SR $ip z ]
     if { $sxi < $xholemax && $sxi > $xholemin && $syi < $yholemax && $syi > $yholemin && $szi > $zholemax } {
        MD++ fixed($ip) = 1
     }
}

MD++ removefixedatoms
MD++ freeallatoms

MD++ plot
#MD++ sleep 

#  move top and bot layer atoms to form dimers
#  Top layer

#0.403 -0.416


set zholemin_tol [expr $zholemax-0.02]
set zholemax_tol [expr $zholemax+0.02]
if { $surface_reconstruct } {
 set ny  [ MD++_Get latticesize(7) ]
 for { set i 0 } { $i < $ny } { incr i 1 } {
    set ymin [ expr -0.5006+1.0/$ny*$i ]
    set ymax [ expr -0.5006+1.0/$ny*($i+0.5) ]
    MD++ input = \[ 1 -10 10 $ymin $ymax 0.3025   10 \]  fixatoms_by_position
    MD++ input = \[ 1 -10 10 $ymin $ymax  -10  -0.310 \]  fixatoms_by_position
    MD++ input = \[ 1 -10 10 $ymin $ymax  $zholemin_tol  $zholemax_tol \]  fixatoms_by_position
 }
 MD++ input = 1  setfixedatomsgroup  freeallatoms

 for { set i 0 } { $i < $ny } { incr i 1 } {
    set ymin [ expr -0.5006+1.0/$ny*($i+0.5) ]
    set ymax [ expr -0.5006+1.0/$ny*($i+1) ]
    MD++ input = \[ 1 -10 10 $ymin $ymax 0.3025   10 \]  fixatoms_by_position
    MD++ input = \[ 1 -10 10 $ymin $ymax -10 -0.310 \]  fixatoms_by_position
    MD++ input = \[ 1 -10 10 $ymin $ymax $zholemin_tol $zholemax_tol \]  fixatoms_by_position
 }
 MD++ input = 2  setfixedatomsgroup  freeallatoms

 MD++ input = \[ 1  0  0.8 0  1 \] movegroup  
 MD++ input = \[ 1  0 -0.8 0  2 \] movegroup  
}
 #end of surf reconstruct


  relax_fixbox

  set H11_0 [ MD++_Get H_11 ] ; set H22_0 [ MD++_Get H_22 ] ; set H33_0 [ MD++_Get H_33 ]

  set epsilon 0.0
  set C11 520000  
  set C44 160000
  set factor 0.7

  set maxitereps 60
  set maxiter    100
  for { set itereps 0 } { $itereps <= $maxitereps } { incr itereps 1 } {

      puts "\n\n epsilon=$epsilon"

      MD++ eval plot

      if { 1 } {

	  set H11_fix [ expr $H11_0*(1.0-$epsilon) ]
          MD++ H_11 = $H11_fix
	  
	  if { $flag == 001 } {
             MD++ { conj_fixboxvec = [ 0 1 1
                                       0 1 1
                                       1 1 1 ] }
          } else {
             puts "unknown flag = $flag, must be 001"
             exitmd
          }

#MD++ conj_ftol = 2e-6 conj_itmax = 1000 conj_fevalmax = 2
          MD++ relax
          MD++ eval

          set H11_fix [ MD++_Get H_11 ] ; set H22_new [ MD++_Get H_22 ] ; set H33_new [ MD++_Get H_33 ] 
          set H31_cur [ MD++_Get H_31 ] ; set H12_cur [ MD++_Get H_12 ] ; set H23_new [ MD++_Get H_23 ]

          set sig_xx [ MD++_Get TSTRESSinMPa_xx ] ; set sig_yy [ MD++_Get TSTRESSinMPa_yy ]
          set sig_zz [ MD++_Get TSTRESSinMPa_zz ] ; set sig_xy [ MD++_Get TSTRESSinMPa_xy ]
          set sig_xz [ MD++_Get TSTRESSinMPa_xz ] ; set sig_yz [ MD++_Get TSTRESSinMPa_yz ]

          puts "iter = $itereps"
	  puts "epsilon = $epsilon"
          puts "sig_xx = $sig_xx sig_yy = $sig_yy sig_zz = $sig_zz"

      } else {

          set bx  [expr 0.5]
	  set by  [expr 0.5]
	  set bz  [expr -0.5*sqrt(2) ] 

          set nx  [expr 1.0/sqrt(2)]
	  set ny  [expr 0]
	  set nz  [expr 1.0/sqrt(2)] 
	 
	  set scaleStress [expr 3500.0*$epsilon/0.051]
	  # calculate stress from burgers and n. 
	  set targetStressMpa_xx [expr ($bx*$nx+$bx*$nx) * $scaleStress]; 
	  set targetStressMpa_yy [expr ($by*$ny+$by*$ny) * $scaleStress];
	  set targetStressMpa_zz [expr ($bz*$nz+$bz*$nz) * $scaleStress];
	  set targetStressMpa_xy [expr ($bx*$ny+$by*$nx) * $scaleStress];
	  set targetStressMpa_xz [expr ($bx*$nz+$bz*$nx) * $scaleStress];
	  set targetStressMpa_yz [expr ($by*$nz+$bz*$ny) * $scaleStress];
	  
	  puts "targetStress = $targetStressMpa_xx,  $targetStressMpa_yy,  $targetStressMpa_zz,  $targetStressMpa_xy,  $targetStressMpa_xz,  $targetStressMpa_yz"

          # adjust other strains in small steps
          for { set iter 0 } { $iter <= $maxiter } { incr iter 1 } {

          set sig_xx [ MD++_Get TSTRESSinMPa_xx ] ; set sig_yy [ MD++_Get TSTRESSinMPa_yy ]
          set sig_zz [ MD++_Get TSTRESSinMPa_zz ] ; set sig_xy [ MD++_Get TSTRESSinMPa_xy ]
          set sig_xz [ MD++_Get TSTRESSinMPa_xz ] ; set sig_yz [ MD++_Get TSTRESSinMPa_yz ]

          set dsig_xx [expr -1.0*( $targetStressMpa_xx -  $sig_xx )] ; 
	  set dsig_yy [expr -1.0*( $targetStressMpa_yy -  $sig_yy )]
          set dsig_zz [expr -1.0*( $targetStressMpa_zz -  $sig_zz )] ;
	  set dsig_xy [expr -1.0*( $targetStressMpa_xy -  $sig_xy )]
          set dsig_xz [expr -1.0*( $targetStressMpa_xz -  $sig_xz )] ; 
	  set dsig_yz [expr -1.0*( $targetStressMpa_yz -  $sig_yz )]

	  #Here e is actually gamma. 
          set e_xx [ expr $dsig_xx / $C11 ] ; set e_yy [ expr $dsig_yy / $C11 ] ; set e_zz [ expr $dsig_zz / $C11 ]
          set e_xy [ expr $dsig_xy / $C44 ] ; set e_xz [ expr $dsig_xz / $C44 ] ; set e_yz [ expr $dsig_yz / $C44 ]

          #puts "iter = $iter"
          #puts "sig_xx = $sig_xx sig_yy = $sig_yy sig_zz = $sig_zz"
          #puts "sig_xy = $sig_xy sig_xz = $sig_xz sig_yz = $sig_yz"
          #puts "e_xx = $e_xx e_yy = $e_yy e_zz = $e_zz"
          #puts "e_xy = $e_xy e_xz = $e_xz e_yz = $epsilon"

          set H11_cur [ MD++_Get H_11 ] ; set H22_cur [ MD++_Get H_22 ] ; set H33_cur [ MD++_Get H_33 ] 
          set H13_cur [ MD++_Get H_13 ] ; set H12_cur [ MD++_Get H_12 ] ; set H23_cur [ MD++_Get H_23 ]

          set H11_new [ expr ${H11_cur}*(1.0+$e_xx*$factor) ] ;        MD++ H_11 = ${H11_new}
          set H22_new [ expr ${H22_cur}*(1.0+$e_yy*$factor) ] ;        MD++ H_22 = ${H22_new}
          set H33_new [ expr ${H33_cur}*(1.0+$e_zz*$factor) ] ;        MD++ H_33 = ${H33_new}
          set H13_new [ expr ${H13_cur} + ${H33_cur}*$e_xz*$factor ] ; MD++ H_13 = ${H13_new}
          set H12_new [ expr ${H12_cur} + ${H22_cur}*$e_xy*$factor ] ; MD++ H_12 = ${H12_new}
          set H23_new [ expr ${H23_cur} + ${H33_cur}*$e_yz*$factor ] ; MD++ H_23 = ${H23_new}

#          set H23_new ${H23_cur}

          if { $iter == [expr $maxiter ] } {
              relax_fixbox		  
#              MD++ conj_ftol = 2e-6 conj_fixbox = 1 relax
          }
          MD++ eval
     }
     set SIGXY [ expr int($sig_xy) ]
      relax_fixbox
#     MD++ conj_ftol = 2e-6 conj_fixbox = 1 relax
     }

     relax_fixbox

     MD++ finalcnfile = "0K_${epsilon}_relaxed_surf${flag}.cn" writecn
     MD++ finalcnfile = "0K_${epsilon}_relaxed_surf${flag}.cfg" writeatomeyecfg

     set e_xx [ expr ($H11_fix-$H11_0)/$H11_0 ] ; set e_yy [ expr ($H22_new-$H22_0)/$H22_0 ]
     set e_zz [ expr ($H33_new-$H33_0)/$H33_0 ] ; set e_yz [ expr ($H23_new)/$H33_0 ]
     set e_xy [ expr ($H12_cur)/$H22_0 ] ;        set e_xz [ expr ($H31_cur)/$H33_0 ]

     set fp [ open $strain_data_file a+ ] ; puts $fp "$e_xx $e_yy $e_zz $e_xy $e_xz $e_yz" ; close $fp
     set fp [ open $stress_data_file a+ ] ; puts $fp "$sig_xx $sig_yy $sig_zz $sig_xy $sig_xz $sig_yz" ; close $fp

     set epsilon [ expr $epsilon+0.001 ]
     set epsilon [format "%.4f" $epsilon]
  }

#  MD++ sleep
  exitmd

# To plot shear stress strain curve
#    cd runs/w-disl-nuc-hetero-0
#    octave
#         load stress.dat ; load strain.dat
#         plot(-strain(:,2),stress(:,2)/1000,'o-');
#         xlabel('\epsilon_{yy}');  ylabel('\sigma_{yy}  (GPa)');
  
} elseif { $status == 1 } {
  # prepare and relax dislocation structure
  MD++ setnolog
  initmd 0
#  readpot

  set epsilon $n

  if { $flag == 001 } {
    #set epsilon 0.04 
    #set epsilon 0.074
    set FEVALMAX 20
  } elseif { $flag == 110 } {
    #set epsilon 0.04
    #set epsilon 0.074
    set FEVALMAX 20
  } else {
     puts "unknown flag = $flag, must be 001 or 110"
     exitmd
  }

# shouldn`t here be reading from state A?

#  make_perfect_crystal 12 12 16 
#  make_dislocation_loop $flag $epsilon
 

#  MD++ vacuumratio = [expr 1.0-1.0/(1.0 + 0.2)]
#  if { $flag == 001 } {
#     # create 001 surface
#     MD++ input = \[ 3 3 0.2 \] changeH_keepR 
#  } elseif { $flag == 110 } {
#     # create 110 surface
#     MD++ input = \[ 1 1 0.2 \] changeH_keepR 
#  } else {
#     puts "unknown flag = $flag, must be 001 or 110"
#     exitmd
#  }



#make dislocation after reconstru but before the straining
  MD++ incnfile = "../w-disl-nuc-hetero-0/0K_0.0_relaxed_surf${flag}.cn" readcn

  setup_window
  openwindow
# MD++ sleep 
#  make_shuffle_dislocation_loop $flag $epsilon
  make_glide_dislocation_loop_1 $flag $epsilon
  make_glide_dislocation_loop_2 $flag $epsilon
  make_glide_dislocation_loop_3 $flag $epsilon


  MD++ incnfile = "../w-disl-nuc-hetero-0/0K_${epsilon}_relaxed_surf${flag}.cn" readcn
  MD++ commit_storedr
   
  MD++ finalcnfile = "w-loop-init.cn" writecn
  MD++ finalcnfile = "w-loop-init.cfg" writeatomeyecfg

  MD++ eval plot

#  MD++ sleep 


#  MD++ eval plot 

#  MD++ conj_ftol = 1e-4 conj_itmax = 30 conj_fevalmax = $FEVALMAX #30 #90 #60 #30
#  MD++ conj_fixbox = 1  relax
#  MD++ finalcnfile = "w-loop-relaxed-surf${flag}.cn" writecn
#  MD++ finalcnfile = "w-loop-relaxed-surf${flag}.cfg" writeatomeyecfg

  MD++ eval plot 

  # this is to prepare the structure for NEB calculation
  #MD++ conj_ftol = 1e-4 conj_itmax = 10 conj_fevalmax = 30
 # MD++ finalcnfile = "w-loop-compression-surf${flag}-eps${epsilon}.cn"  writecn
 # MD++ finalcnfile = "w-loop-compression-surf${flag}-eps${epsilon}.cfg" writeatomeyecfg
#  MD++ sleep



  # this is to confirm the structure is beyond critical
  MD++ conj_ftol = 1e-4 conj_itmax = 3800 conj_fevalmax = 50

  MD++ conj_fixbox = 1  relax

#  MD++ conj_ftol = 1e-5 conj_itmax = 3800 conj_fevalmax = 30
#  MD++ conj_fixbox = 1 relax

  MD++ finalcnfile = "w-loop-compression-surf${flag}-eps${epsilon}.cn"  writecn
  MD++ finalcnfile = "w-loop-compression-surf${flag}-eps${epsilon}.cfg" writeatomeyecfg

  MD++ eval plot

  MD++ sleep
  exitmd
  
} elseif { $status == 2 } {
  # find minimum energy path for dislocation loop nucleation
  #MD++ setnolog

 puts "I am here 1"
  initmd $flag-$n
#  readpot
puts "I am here 2"
  set epsilon $n

  MD++ incnfile = "../w-disl-nuc-hetero-0/0K_${epsilon}_relaxed_surf${flag}.cn" readcn  saveH
  MD++ incnfile = "../w-disl-nuc-hetero-0/w-loop-compression-surf${flag}-eps${epsilon}.cn" readcn   restoreH SHtoR  setconfig2
  MD++ eval
  MD++ incnfile = "../w-disl-nuc-hetero-0/0K_${epsilon}_relaxed_surf${flag}.cn" readcn  setconfig1
  MD++ eval 

  #setup_window
  #openwindow

  MD++ fixallatoms constrain_fixedatoms freeallatoms
  MD++ chainlength = 25 totalsteps = 1000

  MD++ timestep = 0.01 printfreq = 2
  MD++ allocchain initRchain
puts "I am here 3"
  set maxiter 20
  for { set iter 0 } { $iter <= $maxiter } { incr iter 1 } {
     #MD++ { nebspec = [ 0 1 0 1 0  0 ] nebrelax }
     ## [ relax_surround k moveleftend moverightend yesclimbimage ]
     #MD++ finalcnfile = nebrelax.chain.cn.iter${iter} writeRchain
     #file rename nebeng.out nebeng.out.iter${iter}
puts " I am here 4"
     MD++ { nebspec = [ 0 1 0 1 1 0 ] stringrelax }
     ## [ relax_surround redistr_freq moveleftend moverightend yesclimbimage islengthconstant ]
     MD++ finalcnfile = stringrelax.chain.cn.iter${iter} writeRchain
     file rename stringeng.out stringeng.out.iter${iter}
  }

  #MD++ sleep
  exitmd
  
} elseif { $status == 3 } {

  # using parallel run to find minimum energy path for dislocation loop nucleation
  #MD++ setnolog
  set ncpu [ MD++_Get numDomains ]
  if { $ncpu < 1 } { set ncpu 1 }
  puts "ncpu = $ncpu"

  initmd $flag-$n.$ncpu

  set epsilon $n

  MD++ Slave_chdir #Slave change to same directory as Master (dirname)
#  readpot

  MD++ incnfile = "../w-disl-nuc-hetero-0/0K_${epsilon}_relaxed_surf${flag}.cn" readcn  saveH
  MD++ incnfile = "../w-disl-nuc-hetero-0/w-loop-compression-surf${flag}.cn" readcn   restoreH SHtoR  setconfig2
  MD++ eval
  MD++ incnfile = "../w-disl-nuc-hetero-0/0K_${epsilon}_relaxed_surf${flag}.cn" readcn  setconfig1
  MD++ eval 

  #setup_window
  #openwindow

  set chainlength [expr $ncpu - 1]
  MD++ chainlength = $chainlength  
  MD++ timestep = 0.01 printfreq = 10

  if { $opt == 0 } {
   if { $n >= 0.09 } {
    MD++ { nebspec = [ 0 1 0 1 0 0 ] totalsteps = 200 equilsteps = 2000 }
   } else {
    MD++ { nebspec = [ 0 1 0 1 0 0 ] totalsteps = 200 equilsteps = 2000 }
   }
  } else {
    MD++ { nebspec = [ 0 1 0 1 0 1 ] totalsteps = 5000 equilsteps = 2000 }
    ## [ relax_surround redistr_freq moveleftend moverightend yesclimbimage islengthconstant ]
  }
  MD++ Broadcast_FS_Param

  MD++ fixallatoms constrain_fixedatoms freeallatoms
  MD++ allocchain_parallel initRchain_parallel

  if { $opt == 1 } {
    MD++ incnfile = "neb.chain.1000" readRchain_parallel
  }

#  MD++ incnfile = "../w-disl-nuc-hetero-$flag-$n.$ncpu/neb.chain.500" readRchain_parallel

#  MD++ { nebspec = [ 0 1 0 1 0 0 ] totalsteps = 20000 }
#  MD++ stringrelax_parallel
#  MD++ finalcnfile = neb.chain.500 writeRchain_parallel

#  MD++ incnfile = "neb.chain.500" readRchain_parallel
#  MD++ { nebspec = [ 0 1 0 1 0 1 ] totalsteps = 20000 }

  MD++ stringrelax_parallel

  if { $opt == 0 } {
    MD++ finalcnfile = neb.chain.1000 writeRchain_parallel
  } else {
    MD++ finalcnfile = neb.chain.2000 writeRchain_parallel
  }

  if { $ncpu == 1 }   {
        MD++ quit
  } else {
        MD++ quit_all
  }
  
} elseif { $status == 4 } {
  # using parallel run to find minimum energy path for dislocation loop nucleation
  # two step relaxation (first with length constant, then with length free)
  #MD++ setnolog
  set ncpu [ MD++_Get numDomains ]
  if { $ncpu < 1 } { set ncpu 1 }
  puts "ncpu = $ncpu"
  #set ncpu 24

  initmd $flag-$n.$ncpu

  set epsilon $n

  MD++ Slave_chdir #Slave change to same directory as Master (dirname)
#  readpot

  MD++ incnfile = "../w-disl-nuc-hetero-0/0K_${epsilon}_relaxed_surf${flag}.cn" readcn  saveH
  MD++ incnfile = "../w-disl-nuc-hetero-0/w-loop-compression-surf${flag}-eps${epsilon}.cn" readcn   restoreH SHtoR  setconfig2
  MD++ eval
  MD++ incnfile = "../w-disl-nuc-hetero-0/0K_${epsilon}_relaxed_surf${flag}.cn" readcn  setconfig1
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

#  setup_window
  #openwindow

  set chainlength [expr $ncpu - 1]
  MD++ chainlength = $chainlength  
  #MD++ timestep = 0.01 printfreq = 10
  MD++ timestep = 0.01 printfreq = 10

  # First relaxation
  if { 0 } {
   if { $n >= 0.2 } {
    MD++ { nebspec = [ 0 1 0 1 0 1 ] totalsteps = 200 equilsteps = 2000 }
  ## [ relax_surround redistr_freq moveleftend moverightend yesclimbimage islengthconstant ]
   } else {
    MD++ { nebspec = [ 0 1 0 1 0 1 ] totalsteps = 1000 equilsteps = 2000 }
   }
  } else {
    MD++ nebspec = \[ 0 1 0 1 0 0 1 \] totalsteps = 100000 equilsteps = 10000
  }   

  MD++ Broadcast_FS_Param

  MD++ fixallatoms constrain_fixedatoms freeallatoms
  MD++ allocchain_parallel initRchain_parallel

  if { $exist > 0 } {
    MD++ incnfile = stringrelax.chain.cn readRchain_parallel
  }

  MD++ eval
  #MD++ stringrelax_parallel
  MD++ stringrelax_parallel_2
  MD++ finalcnfile = stringrelax.chain.cn writeRchain_parallel

  if { 0 } {
  exec cp stringeng.out stringeng_step1.out

  MD++ nebspec = \[ 0 1 0 1 0 1 \] totalsteps = 20000 equilsteps = 100000
  MD++ Broadcast_nebsetting
  MD++ stringrelax_parallel_1
  #we adopt William's revised formulation
  MD++ finalcnfile = stringrelax.chain.cn writeRchain_parallel
  exec cp stringeng.out stringeng_step2.out
  }

  if { $ncpu == 1 }   {
        MD++ quit
  } else {
        MD++ quit_all
  }
  
} elseif { $status == 5 } {
  # anneal surface with MD
  MD++ setnolog
  set ncpu [ MD++_Get numDomains ]
  if { $ncpu < 1 } { set ncpu 1 }
  puts "ncpu = $ncpu"

  initmd $flag-$n.$ncpu

  set epsilon $n

  readpot

  MD++ incnfile = "../w-disl-nuc-hetero-0/0K_${epsilon}_relaxed_surf${flag}.cn" readcn  saveH
  MD++ eval

  setup_window
  MD++ plot_color_windows = 4
#  openwindow

  setup_md
  MD++  T_OBJ = 100   initvelocity  totalsteps = 1000 fixbox = 1  run

  MD++ eval
  relax_fixbox 

  MD++ sleep

  exitmd 
} elseif { $status == 10 } {
  # visualization (from stringrelax serial runs)
  MD++ setnolog
  initmd "view" 
  readpot

#  set epsilon 0.09
  set epsilon $n

  MD++ incnfile = "../w-disl-nuc-hetero-0/0K_${epsilon}_relaxed_surf${flag}.cn" readcn  saveH
  MD++ incnfile = "../w-disl-nuc-hetero-0/w-loop-compression-surf${flag}.cn" readcn restoreH SHtoR setconfig2
  MD++ eval
  MD++ incnfile = "../w-disl-nuc-hetero-0/0K_${epsilon}_relaxed_surf${flag}.cn" readcn  setconfig1
  MD++ eval 

  if { 0 } {
    setup_window
    MD++ NCS = 8 eval calcentralsymmetry 
    MD++ incnfile = "../w-disl-nuc-hetero0/w-loop-compression-surf${flag}.cn" readcn   restoreH SHtoR 
    MD++ finalcnfile = "w-loop-compression-surf${flag}.cfg" writeatomeyecfg
    openwindow
    MD++ sleep
    exitmd
  }

  MD++ zipfiles = 0
  for { set iter 0 } { $iter <= 25 } { incr iter 1 } {
     MD++ incnfile ="../w-disl-nuc-hetero-${flag}-${epsilon}/stringrelax.chain.cn.iter0" readRchain
     MD++ finalcnfile = stringrelax.chain.cn.iter0 writeRchain

     # set chain number
     set chain_no $iter
     set total_no 25

     MD++ input = $iter  copyRchaintoCN eval

     if { $iter == 0 } {
         openwindow
     }
     MD++ eval plot 
     # write cn or cfg file
     MD++ finalcnfile = "chain_no_${chain_no}.cn"  writecn
     MD++ finalcnfile = "chain_no_${chain_no}.cfg" writeatomeyecfg
     MD++ sleepseconds = 1  sleep

  }
  MD++ sleepseconds = 100  sleep
  exitmd

} elseif { $status == 20 } {
  # visualization (from stringrelax_parallel runs)
  MD++ setnolog
  initmd "view-$n" 
#  readpot

  #set epsilon 0.084
  set epsilon $n
  set chain_no 0
  #set chain_no 8
  set total_no 23
  MD++ chainlength = $total_no

  MD++ incnfile = "../w-disl-nuc-hetero-0/0K_${epsilon}_relaxed_surf${flag}.cn" readcn  saveH
  MD++ incnfile = "../w-disl-nuc-hetero-${flag}-${epsilon}.24/stringrelax.chain.cn" 
  MD++ input = $chain_no  readRchain_parallel_toCN  RHtoS

  #MD++ clearR0 refreshnnlist eval quit

#  setup_window

  MD++ zipfiles = 0
  for { set iter 0 } { $iter <= $total_no } { incr iter 1 } {
     # set chain number
     set chain_no $iter
     MD++ input = $chain_no  readRchain_parallel_toCN  RHtoS
     MD++ NCS = 8 eval calcentralsymmetry 

     if { $iter == 0 } {
#         openwindow
     }
     MD++ eval plot 

     if { $iter == 1} {
     #    relax_fixbox
     #    MD++ eval plot sleepseconds = 100 sleep
     }
     # write cn or cfg file
     #MD++ finalcnfile = "chain_no_${chain_no}.cn"  writecn
     MD++ finalcnfile = "chain_no_${chain_no}.cfg" writeatomeyecfg
     MD++ sleepseconds = 1  sleep

  }
  MD++ sleepseconds = 100  sleep
  exitmd

  # This is to adjust Ec[0] for a relaxed path
} elseif { $status == 23 } { 

  set ncpu 24
  #MD++ setnolog
  initmd $flag-$n.$ncpu
#  readpot

  set epsilon $n
  set chain_no 0
  set total_no 23
  MD++ chainlength = $total_no

  set fp [ open "EPOT.dat" a+ ]

  MD++ incnfile = "../w-disl-nuc-hetero-0/0K_${epsilon}_relaxed_surf${flag}.cn" readcn 
  MD++ eval
  set EPOT_StateA [ MD++_Get EPOT ]
  puts $fp [ MD++_Get EPOT ]
  puts "EPOT_StateA = $EPOT_StateA"
  #set chainfile "../w-disl-nuc-hetero-${flag}-${epsilon}.$ncpu.saved/stringrelax.chain.cn" 
  set chainfile "../w-disl-nuc-hetero-${flag}-${epsilon}.$ncpu/stringrelax.chain.cn" 
  MD++ incnfile = $chainfile
  for { set iter 0 } { $iter <= $total_no } { incr iter 1 } {
     set chain_no $iter
     MD++ input = $chain_no  readRchain_parallel_toCN  RHtoS
     MD++ eval
     set EPOT_Chain0 [ MD++_Get EPOT ]
     puts "EPOT_Chain0 = $EPOT_Chain0"
     puts $fp [ MD++_Get EPOT ]
     break
  }
  puts $fp [expr $EPOT_StateA - $EPOT_Chain0 ]
  close $fp
} elseif { $status == 21 } {
  MD++ incnfile = "runs/w-disl-nuc-hetero-0/w-loop-compression-surf001-eps0.0600.cn" readcn 
  MD++ finalcnfile = "final.cfg" writeatomeyecfg
} else {
        
 puts "unknown status = $status"
 exitmd 

} 

