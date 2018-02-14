# -*-shell-script-*-
# dislocation structure in Nickel
# make eam build=R SYS=mc2 SPEC+=-DMAXCONSTRAINATOMS=800001

source "$::env(MDPLUS_DIR)/scripts/Examples/Tcl/startup.tcl"

#*******************************************
# Definition of procedures
#*******************************************
proc initmd { status { n 0 } } {
MD++ setnolog
MD++ setoverwrite
MD++ dirname = runs/ni-cs/ni-Rao-screw-30x30x20-$status
MD++ NNM = 200
exec cp -r $::env(MDPLUS_DIR)/scripts/work/NiCrossSlip/ni_Rao_screw_cross_slip_Yanming_v1.tcl .  
                                                                            
}

proc readpot { } { MD++ {
#--------------------------------------------
#Read in potential file
#potfile = "~/Codes/MD++.svn/potentials/EAMDATA/eamdata.Ni.LAMMPS.jea.eam.fs" 
#eamgrid = 1000 readeam
potfile = /home/xzhang11/MD++/potentials/EAMDATA/eamdata.Ni.Rao99 
eamgrid = 5000 readeam
} }

proc make_perfect_crystal {  } { 
#------------------------------------------------------------
#Create Perfect Lattice Configuration
#
MD++ {
    latticestructure = face-centered-cubic
    latticeconst = 3.52 #(A) for Ni
    element0 = "Ni"
    latticesize = [  1 -1  0  30    #(x)
                     1  1  1  30    #(y)
                    -1 -1  2  20  ] #(z)            
    #latticesize = [  1 -1  0  40    #(x)
    #                 1  1  1  40    #(y)
    #                -1 -1  2  20  ] #(z)            
    #latticesize = [  1 -1  0  40    #(x)
    #                 1  1  1  25    #(y)
    #                -1 -1  2  20  ] #(z)            
    makecn finalcnfile = perf.cn writecn #refreshnnlist eval quit
}
    MD++ makecrystal #finalcnfile = perf.cn writecn #eval

    # remove top and bottom layers of atoms
    MD++ input = \[ 1 -1 1 -1 -0.40 -1 1 \] fixatoms_by_position
    MD++ input = \[ 1 -1 1  0.40  1 -1 1 \] fixatoms_by_position
    MD++ removefixedatoms    
        
    # remove front and back layers of atoms
    #MD++ input = \[ 1 -1 1 -1 1 -1 -0.40 \] fixatoms_by_position
    #MD++ input = \[ 1 -1 1 -1 1 0.40  1  \] fixatoms_by_position
    #MD++ removefixedatoms    

    # fix top and bottom layers of atoms
    MD++ input = \[ 1  0.00 0.38 -1 -0.37 -1 1 \] fixatoms_by_position
    MD++ input = \[ 1 -0.38 0.00  0.37  1 -1 1 \] fixatoms_by_position
    MD++ input = \[ 10 \] setfixedatomsgroup    
    MD++ freeallatoms
}
        
proc make_screw_dipole { store } { 
#------------------------------------------------------------
#Introduce screw dislocation dipole
#
    set bx [expr 0.5 / [MD++_Get latticesize(3) ] ]
    set by 0.0
    set bz 0.0
    
    # dislocation position
    set z0 [expr 0.111 /   [MD++_Get latticesize(3)] ]
    set y0 [expr -0.5+0.125/[MD++_Get latticesize(7)] ]
    set y1 [expr $y0 + 0.5]
    MD++ input = \[ 1 2 $bx $by $bz  $z0 $y0 $y1 0.305 -10 10 -10 10 1 0 0 0 $store \]
    MD++ makedipole    
}

proc make_glide_partials { store } {
#
    set a 3.52
    set bx [expr 0.5*0.5*sqrt(2)]
    set by 0.0
    set bz 0.0
    set Lx [MD++_Get H_11]
    set Ly [MD++_Get H_22]
    set d  30.0
    set nu 0.3
    set x1 [expr $Lx]
    set y1 [expr $Ly*0.125/[MD++_Get latticesize(7)] ]
    set z1 [expr $d/2]
    set x2 $x1
    set y2 $y1
    set z2 0
    set x3 [expr -$x1]
    set y3 $y1
    set z3 0
    set x4 [expr -$x1]
    set y4 $y1
    set z4 $z1
    MD++ input = \[ 1 1 $nu $a  $bx $by $bz  4 $x1 $y1 $z1 $x2 $y2 $z2 $x3 $y3 $z3 $x4 $y4 $z4 \]
    MD++ makedislpolygon

    set x1 [expr $Lx]
    set y1 [expr $Ly*0.125/[MD++_Get latticesize(7)] ]
    set z1 [expr -$d/2]
    set x2 $x1
    set y2 $y1
    set z2 0
    set x3 [expr -$x1]
    set y3 $y1
    set z3 0
    set x4 [expr -$x1]
    set y4 $y1
    set z4 $z1
    MD++ input = \[ 1 $store $nu $a  $bx $by $bz  4 $x1 $y1 $z1 $x2 $y2 $z2 $x3 $y3 $z3 $x4 $y4 $z4 \]
    MD++ makedislpolygon

}

proc make_cross_partials { store } {
#
    set a 3.52
    set bx [expr 0.5*0.5*sqrt(2)]
    set by 0.0
    set bz 0.0
    set Lx [MD++_Get H_11]
    set Ly [MD++_Get H_22]
    set d  30.0
    set nu 0.3
    set y0 [expr $Ly*0.125/[MD++_Get latticesize(7)] ]
    set theta [expr acos(1.0/3.0)]
    set x1 [expr $Lx]
    set y1 [expr $y0+$d*sin($theta)/2] 
    set z1 [expr $d*cos($theta)/2]
    set x2 $x1
    set y2 $y0
    set z2 0
    set x3 [expr -$x1]
    set y3 $y0
    set z3 0
    set x4 [expr -$x1]
    set y4 $y1
    set z4 $z1
    MD++ input = \[ 1 1 $nu $a  $bx $by $bz  4 $x1 $y1 $z1 $x2 $y2 $z2 $x3 $y3 $z3 $x4 $y4 $z4 \]
    MD++ makedislpolygon

    set x1 [expr $Lx]
    set y1 [expr $y0-$d*sin($theta)/2] 
    set z1 [expr -$d*cos($theta)/2]
    set x2 $x1
    set y2 $y0
    set z2 0
    set x3 [expr -$x1]
    set y3 $y0
    set z3 0
    set x4 [expr -$x1]
    set y4 $y1
    set z4 $z1
    MD++ input = \[ 1 $store $nu $a  $bx $by $bz  4 $x1 $y1 $z1 $x2 $y2 $z2 $x3 $y3 $z3 $x4 $y4 $z4 \]
    MD++ makedislpolygon
}

proc make_cross_slip { store } {
#
    set a 3.52
    set bx [expr 0.5*0.5*sqrt(2)]
    set by 0.0
    set bz 0.0
    set Lx [MD++_Get H_11]
    set Ly [MD++_Get H_22]
    set d  30.0
    set nu 0.3
# make glide partials
    set x1 [expr $Lx]
    set y1 [expr $Ly*0.125/[MD++_Get latticesize(7)] ]
    set z1 [expr $d/2]
    set x2 $x1
    set y2 $y1
    set z2 0
    set x3 [expr $Lx/4]
    set y3 $y1
    set z3 0
    set x4 $x3
    set y4 $y1
    set z4 $z1
    MD++ input = \[ 1 1 $nu $a  $bx $by $bz  4 $x1 $y1 $z1 $x2 $y2 $z2 $x3 $y3 $z3 $x4 $y4 $z4 \]
    MD++ makedislpolygon
    set x1 [expr -$Lx/4]
    set y1 [expr $Ly*0.125/[MD++_Get latticesize(7)] ]
    set z1 [expr $d/2]
    set x2 $x1
    set y2 $y1
    set z2 0
    set x3 [expr -$Lx]
    set y3 $y1
    set z3 0
    set x4 $x3
    set y4 $y1
    set z4 $z1
    MD++ input = \[ 1 1 $nu $a  $bx $by $bz  4 $x1 $y1 $z1 $x2 $y2 $z2 $x3 $y3 $z3 $x4 $y4 $z4 \]
    MD++ makedislpolygon
    set x1 [expr $Lx]
    set y1 [expr $Ly*0.125/[MD++_Get latticesize(7)] ]
    set z1 [expr -$d/2]
    set x2 $x1
    set y2 $y1
    set z2 0
    set x3 [expr $Lx/4]
    set y3 $y1
    set z3 0
    set x4 $x3
    set y4 $y1
    set z4 $z1
    MD++ input = \[ 1 1 $nu $a  $bx $by $bz  4 $x1 $y1 $z1 $x2 $y2 $z2 $x3 $y3 $z3 $x4 $y4 $z4 \]
    MD++ makedislpolygon
    set x1 [expr -$Lx/4]
    set y1 [expr $Ly*0.125/[MD++_Get latticesize(7)] ]
    set z1 [expr -$d/2]
    set x2 $x1
    set y2 $y1
    set z2 0
    set x3 [expr -$Lx]
    set y3 $y1
    set z3 0
    set x4 $x3
    set y4 $y1
    set z4 $z1
    MD++ input = \[ 1 1 $nu $a  $bx $by $bz  4 $x1 $y1 $z1 $x2 $y2 $z2 $x3 $y3 $z3 $x4 $y4 $z4 \]
    MD++ makedislpolygon

# make cross slip partials
    set y0 [expr $Ly*0.125/[MD++_Get latticesize(7)] ]
    set theta [expr acos(1.0/3.0)]
    set x1 [expr $Lx/4]
    set y1 [expr $y0+$d*sin($theta)/2] 
    set z1 [expr $d*cos($theta)/2]
    set x2 $x1
    set y2 $y0
    set z2 0
    set x3 [expr -$x1]
    set y3 $y0
    set z3 0
    set x4 [expr -$x1]
    set y4 $y1
    set z4 $z1
    MD++ input = \[ 1 1 $nu $a  $bx $by $bz  4 $x1 $y1 $z1 $x2 $y2 $z2 $x3 $y3 $z3 $x4 $y4 $z4 \]
    MD++ makedislpolygon

    set y0 [expr $Ly*0.125/[MD++_Get latticesize(7)] ]
    set theta [expr acos(1.0/3.0)]
    set x1 [expr $Lx/4]
    set y1 [expr $y0-$d*sin($theta)/2] 
    set z1 [expr -$d*cos($theta)/2]
    set x2 $x1
    set y2 $y0
    set z2 0
    set x3 [expr -$x1]
    set y3 $y0
    set z3 0
    set x4 [expr -$x1]
    set y4 $y1
    set z4 $z1
    MD++ input = \[ 1 $store $nu $a  $bx $by $bz  4 $x1 $y1 $z1 $x2 $y2 $z2 $x3 $y3 $z3 $x4 $y4 $z4 \]
    MD++ makedislpolygon
}

proc make_glide_partials_Fleischer { store } {
#
    set a 3.52
    set bx [expr 0.5*0.5*sqrt(2)]
    set by 0.0
    set bz 0.0
    set Lx [MD++_Get H_11]
    set Ly [MD++_Get H_22]
    set d  30.0
    set nu 0.3
    set x1 [expr $Lx]
    set y1 [expr $Ly*0.125/[MD++_Get latticesize(7)] ]
#    set z1 [expr -$d]
    set z1 [expr $d]
    set x2 $x1
    set y2 $y1
    set z2 0
    set x3 [expr -$x1]
    set y3 $y1
    set z3 0
    set x4 [expr -$x1]
    set y4 $y1
    set z4 $z1
    MD++ input = \[ 1 $store $nu $a  $bx $by $bz  4 $x1 $y1 $z1 $x2 $y2 $z2 $x3 $y3 $z3 $x4 $y4 $z4 \]
    MD++ makedislpolygon
}

proc make_cross_partials_Fleischer { store } {
#
    set a 3.52
    set bx [expr 0.5*0.5*sqrt(2)]
    set by 0.0
    set bz 0.0
    set Lx [MD++_Get H_11]
    set Ly [MD++_Get H_22]
    set d  30.0
    set nu 0.3
    set y0 [expr $Ly*0.125/[MD++_Get latticesize(7)] ]
    set theta [expr acos(1.0/3.0)]
    set x1 [expr $Lx]
    set y1 [expr $y0+$d*sin($theta)] 
    set z1 [expr $d*cos($theta)]
    set x2 $x1
    set y2 $y0
    set z2 0
    set x3 [expr -$x1]
    set y3 $y0
    set z3 0
    set x4 [expr -$x1]
    set y4 $y1
    set z4 $z1
    MD++ input = \[ 1 $store $nu $a  $bx $by $bz  4 $x1 $y1 $z1 $x2 $y2 $z2 $x3 $y3 $z3 $x4 $y4 $z4 \]
    MD++ makedislpolygon
}

proc make_cross_slip_Fleischer { store } {
#
    set a 3.52
    set bx [expr 0.5*0.5*sqrt(2)]
    set by 0.0
    set bz 0.0
    set Lx [MD++_Get H_11]
    set Ly [MD++_Get H_22]
    set d  30.0
    set nu 0.3

# Made differently than standard cross slip
# One glide partial with one side (A) coinciding with dislocation
# One cross slip partial with one side (A) coinciding with dislocation

# make glide partial Fleischer
    set x1 [expr $Lx]
    set y1 [expr $Ly*0.125/[MD++_Get latticesize(7)] ]
#    set z1 [expr -$d]
    set z1 [expr $d]
    set x2 $x1
    set y2 $y1
    set z2 0
    set x3 [expr -$Lx]
    set y3 $y1
    set z3 0
    set x4 $x3
    set y4 $y1
    set z4 $z1
    MD++ input = \[ 1 1 $nu $a  $bx $by $bz  4 $x1 $y1 $z1 $x2 $y2 $z2 $x3 $y3 $z3 $x4 $y4 $z4 \]
    MD++ makedislpolygon
# make cross slip partial Fleischer
    set y0 [expr $Ly*0.125/[MD++_Get latticesize(7)] ]
    set theta [expr acos(1.0/3.0)]
    set x1 [expr $Lx/4]
    set y1 [expr $y0+$d*sin($theta)] 
    set z1 [expr $d*cos($theta)]
    set x2 $x1
    set y2 $y0
    set z2 0
    set x3 [expr -$x1]
    set y3 $y0
    set z3 0
    set x4 [expr -$x1]
    set y4 $y1
    set z4 $z1
    MD++ input = \[ 1 $store $nu $a  $bx $by $bz  4 $x1 $y1 $z1 $x2 $y2 $z2 $x3 $y3 $z3 $x4 $y4 $z4 \]
    MD++ makedislpolygon
}


proc mark_groups { } {
#------------------------------------------------------------
#Mark surface atoms for applying stress later
#
  MD++ NCS = 12 eval calcentralsymmetry freeallatoms
  
  #identify top surface atom
  MD++ input = \[ 1 -10 10 0 10 -10 10 10 20 \] fixatoms_by_pos_topol
  set  ntop     [MD++_Get NPfixed]
  MD++ input = 1 setfixedatomsgroup freeallatoms
  
  #identify bottom surface atom
  MD++ input = \[ 1 -10 10 -10 0 -10 10 10 20 \] fixatoms_by_pos_topol
  set  nbot     [MD++_Get NPfixed]
  MD++ input = 2 setfixedatomsgroup freeallatoms

  set natoms [list $ntop $nbot]
  return $natoms
}

#--------------------------------------------
proc relax_fixbox { } { MD++ {
# Conjugate-Gradient relaxation
conj_ftol = 1e-4 conj_itmax = 1000 conj_fevalmax = 1000
conj_fixbox = 1
relax
} }
#end of proc relax_fixbox


#--------------------------------------------
proc relax_freebox { } { MD++ {
# Conjugate-Gradient relaxation
conj_ftol = 1e-4 conj_itmax = 1000 conj_fevalmax = 1000 #10
conj_fixbox = 0
conj_fixboxvec = [ 0 0 0
                   1 0 1
                   1 0 0 ]
relax
} }
#end of proc relax_freebox

proc setup_window { } { MD++ {
#------------------------------------------------------------
#colors for Central symmetry view
color00 = "red" color01 = "blue" color02 = "green"
color03 = "magenta" color04 = "cyan" color05 = "purple"
color06 = "gray80" color07 = "white" color08 = "orange"
#-------------------------------------------------------------
#Plot Configuration
#Plot Configuration
atomradius = 1.0 bondradius = 0.3 bondlength = 0
atomcolor = cyan highlightcolor = purple backgroundcolor = gray
bondcolor = red fixatomcolor = yellow
color00 = "orange"  color01 = "red"    color02 = "green"
color03 = "magenta" color04 = "cyan"   color05 = "purple"
color06 = "gray80"  color07 = "white"  color08 = "blue"
#
plot_color_axis = 2  #2: use CSD (default 0: use local energy)
#
plot_color_windows = [ 3
                       1  3  8 #color08 = blue
                       3 20  6 #color06 = gray80
                      20 100 0 #color00 = orange
                       1
                      ]
#plot_limits = [ 1 -10 10 -10 10 -0.37 0.37 ]
plot_atom_info = 2  plotfreq = 10 
rotateangles = [ 0 0 0 1.3 ]
#
win_width = 600 win_height = 600
#openwin alloccolors rotate saverot refreshnnlist eval plot
} }

proc openwindow { } { 
setup_window
MD++ openwin alloccolors rotate saverot refreshnnlist eval plot
}

#--------------------------------------------
proc exitmd { } { MD++ quit }
#end of proc exitmd
#--------------------------------------------

#--------------------------------------------
proc setup_md { } { MD++ {     
equilsteps = 0  timestep = 0.001 # (ps)
atommass = 58.71 # (g/mol)
DOUBLE_T = 1
saveprop = 1 savepropfreq = 1000 openpropfile
savecn = 1 savecnfreq = 1000 openintercnfile
savecfg = 1 savecfgfreq = 10000
plotfreq = 10 printfreq = 10
randseed = 12345 srand48  #randomize random number generator
#vt2 = 1e28  #1e28 2e28 5e28
NHMass = 2e-3
wallmass = 2e3     # atommass * NP = 14380
boxdamp = 1e-3     # optimal damping for 216 atoms and wallmass 1e-3
saveH # Use current H as reference (H0), needed for specifying stress
fixbox  =  1
output_fmt = "curstep EPOT KATOM Tinst HELM HELMP TSTRESSinMPa_xx TSTRESSinMPa_yy TSTRESSinMPa_zz TSTRESSinMPa_xy TSTRESSinMPa_xz TSTRESSinMPa_yz H_11 H_22 H_33" 
ensemble_type = "NVT" integrator_type = "VVerlet" implementation_type = 1
writeall = 1
} }
#end of proc setup_md
#-------------------------------------------
#-------------------------------------------
proc swap_velocities { NP v_ratio } { 
puts "Perturb velocities ..."
  set maxAtomsPerturbed  [expr ($NP*1.0)/$v_ratio]
  set vx [MD++_Get VSR(0)]
  set vy [MD++_Get VSR(1)]
  set vz [MD++_Get VSR(2)]
  for { set i 0 } { $i < $maxAtomsPerturbed } { incr i } {
     set atom [expr { int(rand()*$NP) } ]
     set ind [expr $atom*3]
     set vx0 [MD++_Get VSR($ind)]
     MD++_Set VSR($ind) $vx
     set ind [expr $atom*3+1]
     set vy0 [MD++_Get VSR($ind)]
     MD++_Set VSR($ind) $vy
     set ind [expr $atom*3+2]
     set vz0 [MD++_Get VSR($ind)]
     MD++_Set VSR($ind) $vz
     if { $i == 0 } {
       puts "vx = $vx, vx0 = $vx0"
       puts "vy = $vy, vy0 = $vy0"
       puts "vz = $vz, vz0 = $vz0"
       puts "..."
     }
     set vx $vx0
     set vy $vy0
     set vz $vz0
  }
}
#-------------------------------------------
#--------------------------------------------
#proc datafile_process { filename index frac fracend operation } {
#   puts "total line $NUM \n"
#   set Nini [ expr round($NUM*$frac) ] 
#   set Nend [ expr round($NUM*$fracend) ]

proc datafile_process { filename index Nlast operation } {
   set fp [ open $filename r ]
   set data [ split [read $fp] \n]
   close $fp
   set NUM [ expr [ llength $data ] - 1 ]

   set Nend $NUM
   set Nini [ expr $Nend - $Nlast - 1 ]

   set Sum 0
   set k   0
   set Var 0
   set Std 0
   for { set i $Nini } { $i < $Nend } { incr i 1 } {
       set k [expr $k+1]
       set data1 [ lindex $data $i ]
       split $data1 
       set datum [ lindex $data1 $index ]
       if { $i == $Nini } { 
	set MAX [lindex $data1 $index]
 	set MIN [lindex $data1 $index]
       } elseif { $i > $Nini } {
        set MAX [ expr ($MAX>$datum)?$MAX:$datum ]
        set MIN [ expr ($MIN<$datum)?$MIN:$datum ]
       }
       set Sum [expr $Sum+$datum]
   }
   set Ave [ expr $Sum/$k]

   if { [ string match "*STD*" $operation ] || [ string match "*VAR*" $operation ] } {
      for { set i $Nini } { $i < $Nend } { incr i 1 } {
          set data1 [ lindex $data $i ]
          split $data1
          set datum [ lindex $data1 $index ]  
	  set Var [expr $Var+($datum-$Ave)*($datum-$Ave) ]
      }
      set Var [ expr $Var/$k ]
      set Std [ expr sqrt($Var) ]
   }
   split $operation
   set Nvar [ llength $operation ]
   for { set i 0 } { $i < $Nvar } { incr i 1 } {
      set var [ lindex $operation $i ]
      if { $var=="SUM" } { 
          lappend LIST $Sum
      } elseif { $var=="AVE" } {
          lappend LIST $Ave
      } elseif { $var=="MAX" } {
          lappend LIST $MAX
      } elseif { $var=="MIN" } {
          lappend LIST $MIN
      } elseif { $var=="VAR" } {
          lappend LIST $Var
      } elseif { $var=="STD" } {
          lappend LIST $Std
      }
   }
   return $LIST
}

#-------------------------------------------
proc fix_atom_strips { } {
# fix atoms to prevent annihilation of double kink
  MD++ input = \[ 0 \] 
 MD++ fixatoms_by_ID

 MD++ input = \[ 0 \]
 MD++ fixatoms_by_ID
}

#-------------------------------------------
proc set_constrainatoms { nconstr } {
# set constrain atoms for constr-relaxation or NEB relax

  if { $nconstr == 183 } {

  MD++ constrainatoms = \[ 0 \]


  } elseif { $nconstr == 327 } {

  MD++ constrainatoms = \[ 0 \]

  } elseif { $nconstr == 500 } {

  MD++ constrainatoms = \[ 0 \]

  } else {
    puts "unrecongnized nconstr ($nconstr)!"
    exit
  }
}

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
 set n2 0
} elseif { $argc > 2 } {
 set n2 [lindex $argv 2]
}
puts "n2 = $n2"

if { $argc <= 3 } {
 set n3 0
} elseif { $argc > 3 } {
 set n3 [lindex $argv 3]
}
puts "n3 = $n3"

if { $argc <= 4 } {
 set n4 0
} elseif { $argc > 4 } {
 set n4 [lindex $argv 4]
}
puts "n4 = $n4"

if { $argc <= 5 } {
 set n5 0
} elseif { $argc > 5 } {
 set n5 [lindex $argv 5]
}
puts "n5 = $n5"

if { $argc <= 6 } {
 set n6 0
} elseif { $argc > 6 } {
 set n6 [lindex $argv 6]
}
puts "n6 = $n6"

if { $argc <= 7 } {
 set n7 0
} elseif { $argc > 7 } {
 set n7 [lindex $argv 7]
}
puts "n7 = $n7"

set USEGPU 0
set myname [ MD++_Get "myname" ]
if { [ string match "*eamcpu*" $myname ] } {
  set USEGPU 0
} elseif { [ string match "*eamgpu*" $myname ] } {
  puts "eamgpu, set USEGPU =1"
  set USEGPU 1
} else { 
  puts "USEGPU set default zero"
}
if { $USEGPU == 1 } {
  MD++ test_saxpy
}
if { $status == -1 } {
  MD++ setnolog
  initmd $status
  readpot

  make_perfect_crystal 

  MD++ input = \[ 1 1 2 \] changeH_keepR
  if { $USEGPU == 1 } {
    MD++ cuda_memcpy_all
  }
  MD++ eval quit
#  MD++ conj_fevalmax = 100 relax eval quit



} elseif { $status == 0 } {
  # prepare and relax dislocation structure on glide plane
  MD++ setnolog
  initmd $status
  readpot

  make_perfect_crystal 


  make_screw_dipole 1
  make_glide_partials 0
  #make_cross_partials 0
  #make_cross_slip 0
  
  #setup_window
  #openwindow
  MD++ finalcnfile = "ni-screw-gp-init.cfg" writeatomeyecfg
  
  MD++ conj_ftol = 1e-4 conj_fixbox = 1 conj_fevalmax = 100
  for { set i 1 } { $i <= 100 } { incr i } { 
      MD++ relax
      MD++ eval
      MD++ finalcnfile = "ni-screw-gp-relax-$i.cn" writecn
      MD++ finalcnfile = "ni-screw-gp-relax-$i.cfg" writeatomeyecfg
  }

  MD++ finalcnfile = "ni-screw-gp-relaxed.cn" writecn
  MD++ eval

  #MD++ sleep
  exitmd
  
} elseif { $status == 0.5 } {
  # prepare and relax dislocation structure on glide plane
  MD++ setnolog
  initmd $status
  readpot

  make_perfect_crystal 

  make_screw_dipole 1

  make_glide_partials_Fleischer 0
  #make_cross_partials 0
  #make_cross_slip 0
  
  setup_window
  #openwindow
  MD++ eval
  MD++ finalcnfile = "ni-screw-gp-init.cfg" writeatomeyecfg
  
  MD++ conj_ftol = 1e-4 conj_fixbox = 1 conj_fevalmax = 100
  for { set i 1 } { $i <= 100 } { incr i } { 
      MD++ relax
      MD++ eval
      MD++ finalcnfile = "ni-screw-gp-relax-$i.cn" writecn
      MD++ finalcnfile = "ni-screw-gp-relax-$i.cfg" writeatomeyecfg
  }

  MD++ finalcnfile = "ni-screw-gp-relaxed.cn" writecn
  MD++ eval

  #MD++ sleep
  exitmd

} elseif { $status == 1 } {
  # prepare and relax dislocation structure on cross-slip plane
  MD++ setnolog
  initmd $status
  readpot

  make_perfect_crystal 

  make_screw_dipole 1

  make_cross_partials 0
 
  setup_window 
  #openwindow
  MD++ eval
  MD++ finalcnfile = "ni-screw-cp-init.cfg" writeatomeyecfg
  
  MD++ conj_ftol = 1e-4 conj_fixbox = 1 conj_fevalmax = 100
  for { set i 1 } { $i <= 100 } { incr i } { 
      MD++ relax
      MD++ eval
      MD++ finalcnfile = "ni-screw-cp-relax-$i.cn" writecn
      MD++ finalcnfile = "ni-screw-cp-relax-$i.cfg" writeatomeyecfg
  }

  MD++ finalcnfile = "ni-screw-cp-relaxed.cn" writecn
  MD++ eval

  #MD++ sleep
  exitmd
  
} elseif { $status == 1.5 } {
  # prepare and relax dislocation structure
  MD++ setnolog
  initmd $status
  readpot

  make_perfect_crystal 

  make_screw_dipole 1

  make_cross_partials_Fleischer 0
 
  setup_window 
  #openwindow
  MD++ eval
  MD++ finalcnfile = "ni-screw-cp-init.cfg" writeatomeyecfg
  
  MD++ conj_ftol = 1e-4 conj_fixbox = 1 conj_fevalmax = 100
  for { set i 1 } { $i <= 100 } { incr i } { 
      MD++ relax
      MD++ eval
      MD++ finalcnfile = "ni-screw-cp-relax-$i.cn" writecn
      MD++ finalcnfile = "ni-screw-cp-relax-$i.cfg" writeatomeyecfg
  }

  MD++ finalcnfile = "ni-screw-cp-relaxed.cn" writecn
  MD++ eval

  #MD++ sleep
  exitmd
  
} elseif { $status == 2 } {
  # prepare and relax dislocation structure
  MD++ setnolog
  initmd $status
  readpot

  make_perfect_crystal 

  make_screw_dipole 1

  make_cross_slip 0
 
  setup_window 
  #openwindow
  MD++ eval
  MD++ finalcnfile = "ni-cross-slip-init.cfg" writeatomeyecfg
  
  MD++ conj_ftol = 1e-4 conj_fixbox = 1 conj_fevalmax = 20
  for { set i 1 } { $i <= 100 } { incr i } { 
      MD++ relax
      MD++ eval
      MD++ finalcnfile = "ni-cross-slip-relax-$i.cn" writecn
      MD++ finalcnfile = "ni-cross-slip-relax-$i.cfg" writeatomeyecfg
  }

  MD++ finalcnfile = "ni-cross-slip-relaxed.cn" writecn
  MD++ eval

  #MD++ sleep
  exitmd
  
} elseif { $status == 2.5 } {
  # prepare and relax dislocation structure
  # Fleischer Mechanism
  MD++ setnolog
  initmd $status
  readpot

  make_perfect_crystal 

  make_screw_dipole 1

#  make_cross_slip 0

  make_cross_slip_Fleischer 0
 
  setup_window 
  #openwindow
  MD++ eval
  MD++ finalcnfile = "ni-cross-slip-init.cfg" writeatomeyecfg
  
  MD++ conj_ftol = 1e-4 conj_fixbox = 1 conj_fevalmax = 20
  for { set i 1 } { $i <= 100 } { incr i } { 
      MD++ relax
      MD++ eval
      MD++ finalcnfile = "ni-cross-slip-relax-$i.cn" writecn
      MD++ finalcnfile = "ni-cross-slip-relax-$i.cfg" writeatomeyecfg
  }

  MD++ finalcnfile = "ni-cross-slip-relaxed.cn" writecn
  MD++ eval

  #MD++ sleep
  exitmd
  
} elseif { $status == 3 } {
  # visualize dislocation structure
  #  compute order parameter
  MD++ setnolog
  initmd $status
  readpot

  MD++ incnfile = "../ni-Rao-screw-30x30x20-0/ni-screw-gp-relax-50.cn" readcn setconfig2

  MD++ incnfile = "../ni-Rao-screw-30x30x20-2/ni-cross-slip-relax-12.cn" readcn
  
  MD++ latticeconst = 3.52 #(A) for Ni
  set a0 3.52
  set rc [expr $a0*sqrt(2.0)/2.0*1.2] 
  set cutoff 0.2
  MD++ UMB_order_param = 30  Rc_for_QLM = $rc
  MD++ L_for_QLM = 1 allocQLM

  set nx 0
  set ny [expr 1.0/3.0]
  set nz [expr -sqrt(8.0)/3.0]
  MD++ SHtoR
  set NP [MD++_Get "NP"]

  puts "set atoms on cross slip plane to group 950129 ..."
  for { set i 0 } { $i < $NP } { incr i } { 
     set ind [expr $i*3]
     set rx [MD++_Get R($ind)]
     set ind [expr $i*3+1]
     set ry [MD++_Get R($ind)]
     set ind [expr $i*3+2]
     set rz [MD++_Get R($ind)]
     set rdotn [expr $rx*$nx + $ry*$ny + $rz*$nz]
     set temp [expr $rdotn/($a0*sqrt(3.0)/4.0)]
     if { $temp > -0.5 && $temp < 0.5 } {
         MD++_Set group($i) 950129
     }
  }

  openwindow

  MD++ QLM_cutoff = $cutoff
  MD++ refreshnnlist caldislocationorder

  MD++ plot_color_axis = 9  #9: use dislocation_order
  MD++ plot_color_windows = \[ 2 $cutoff 10 8  -10 -$cutoff  0  1 \]
  MD++ plot

  MD++ finalcnfile = cross-slip-order-param.cfg writeatomeyecfg

  MD++ sleep
  exitmd
  
} elseif { $status == 3.5 } {
  # visualize Fleischer dislocation structure
  #  compute order parameter
  MD++ setnolog
  initmd $status
  readpot

  MD++ incnfile = "../ni-Rao-screw-30x30x20-0.5/ni-screw-gp-relax-50.cn" readcn setconfig2

  MD++ incnfile = "../ni-Rao-screw-30x30x20-2.5/ni-cross-slip-relax-12.cn" readcn
  
  MD++ latticeconst = 3.52 #(A) for Ni
  set a0 3.52
  set rc [expr $a0*sqrt(2.0)/2.0*1.2] 
  set cutoff 0.2
  MD++ UMB_order_param = 30  Rc_for_QLM = $rc
  MD++ L_for_QLM = 1 allocQLM

  set nx 0
  set ny [expr 1.0/3.0]
  set nz [expr -sqrt(8.0)/3.0]
  MD++ SHtoR
  set NP [MD++_Get "NP"]

  puts "set atoms on cross slip plane to group 950129 ..."
  for { set i 0 } { $i < $NP } { incr i } { 
     set ind [expr $i*3]
     set rx [MD++_Get R($ind)]
     set ind [expr $i*3+1]
     set ry [MD++_Get R($ind)]
     set ind [expr $i*3+2]
     set rz [MD++_Get R($ind)]
     set rdotn [expr $rx*$nx + $ry*$ny + $rz*$nz]
     set temp [expr $rdotn/($a0*sqrt(3.0)/4.0)]
     if { $temp > -0.5 && $temp < 0.5 } {
         MD++_Set group($i) 950129
     }
  }

  openwindow

  MD++ QLM_cutoff = $cutoff
  MD++ refreshnnlist caldislocationorder

  MD++ plot_color_axis = 9  #9: use dislocation_order
  MD++ plot_color_windows = \[ 2 $cutoff 10 8  -10 -$cutoff  0  1 \]
  MD++ plot

  MD++ finalcnfile = cross-slip-order-param.cfg writeatomeyecfg

  MD++ sleep
  exitmd
  
} elseif { $status == 4 } {
  # visualize dislocation structure
  #  select atoms on one plane to be analysed for cross slip (status == 4)
  MD++ setnolog
  initmd $status
  readpot

  # set atoms on cross slip plane to group 950129

  #MD++ incnfile = "../ni-Rao-screw-30x30x20-0/ni-screw-gp-relax-50.cn" readcn
  MD++ incnfile = "../ni-Rao-screw-30x30x20-2/ni-cross-slip-relax-50.cn" readcn
  
  MD++ latticeconst = 3.52 #(A) for Ni
  set a0 3.52
  set nx 0
  set ny [expr 1.0/3.0]
  set nz [expr -sqrt(8.0)/3.0]
  MD++ SHtoR
  set NP [MD++_Get "NP"]

  puts "find atoms on cross slip plane..."
  for { set i 0 } { $i < $NP } { incr i } { 
     set ind [expr $i*3]
     set rx [MD++_Get R($ind)]
     set ind [expr $i*3+1]
     set ry [MD++_Get R($ind)]
     set ind [expr $i*3+2]
     set rz [MD++_Get R($ind)]
     set rdotn [expr $rx*$nx + $ry*$ny + $rz*$nz]
     set temp [expr $rdotn/($a0*sqrt(3.0)/4.0)]
     if { $temp > -0.5 && $temp < 0.5 } {
         MD++_Set fixed($i)  1
     }
  }
  
  openwindow
  MD++ eval

  MD++ sleep
  exitmd
  
} elseif { $status == 4.5 } {
  # visualize dislocation structure
  #  select atoms on one plane to be analysed for cross slip (status == 3)
  MD++ setnolog
  initmd $status
  readpot

  # set atoms on cross slip plane to group 950129

  #MD++ incnfile = "../ni-Rao-screw-30x30x20-0.5/ni-screw-gp-relax-50.cn" readcn
  MD++ incnfile = "../ni-Rao-screw-30x30x20-2.5/ni-cross-slip-relax-50.cn" readcn
  #MD++ incnfile = "../ni-Rao-screw-30x30x20-2.5/ni-cross-slip-relax-5.cn" readcn
  
  MD++ latticeconst = 3.52 #(A) for Ni
  set a0 3.52
  set nx 0
  set ny [expr 1.0/3.0]
  set nz [expr -sqrt(8.0)/3.0]
  MD++ SHtoR
  set NP [MD++_Get "NP"]

  puts "find atoms on cross slip plane..."
  for { set i 0 } { $i < $NP } { incr i } { 
     set ind [expr $i*3]
     set rx [MD++_Get R($ind)]
     set ind [expr $i*3+1]
     set ry [MD++_Get R($ind)]
     set ind [expr $i*3+2]
     set rz [MD++_Get R($ind)]
     set rdotn [expr $rx*$nx + $ry*$ny + $rz*$nz]
     set temp [expr $rdotn/($a0*sqrt(3.0)/4.0)]
     if { $temp > -0.5 && $temp < 0.5 } {
         MD++_Set fixed($i)  1
     }
  }
  
  openwindow
  MD++ eval

  MD++ sleep
  exitmd
  
} elseif { $status == 5 } {
 # relax screw dislocation structure A with stresses

 # glide plane normal: y-axis
 # dislocation line/Burger vector: x-axis
 
 # Schmid stress on glide plane: sigma_xy (must be 0)

 # Escaig stress on glide plane:      sigma_yz (applied by Fext)
 # Escaig stress on cross slip plane: sigma_zz (applied by H)
 # Escaig stress on cross slip plane: sigma_yy (applied by Fext)
 # Schmid stress on cross slip plane: sigma_xz (applied by H)

  set sigma_xx   0.0
  set sigma_xy   0.0
  set sigma_yz   $n2
  set sigma_zz   $n3
  set sigma_yy   $n4
  set sigma_xz   $n5

  set modulus_xx 450e3 ; #MPa
  set modulus_zz 450e3 ; #MPa
  set modulus_xz 100e3 ; #MPa

  #MD++ setnolog
  MD++ setoverwrite
  file mkdir "runs/ni-Rao-screw-30x30x20-$status"
  MD++ dirname = runs/ni-Rao-screw-30x30x20-$status/stress_${n2}_${n3}_${n4}_${n5}
  MD++ NNM = 200

  readpot
  MD++ incnfile = "../../ni-Rao-screw-30x30x20-0/perf.cn" readcn  writeall = 1
  set NP0 [MD++_Get "NP"]

  MD++ incnfile = "../../ni-Rao-screw-30x30x20-0/ni-screw-gp-relaxed.cn" readcn  writeall = 1

  # test new state A at Schmid stress 400 (does not relax well)
  #MD++ incnfile = "../../ni-Rao-screw-30x30x20-0/chain_no_7.cn" readcn freeallatoms writeall = 1

  setup_window
  #openwindow

  set NP [MD++_Get "NP"]
  set vacuumratio  [expr 1.0 - ($NP*1.0)/$NP0]
  puts "vacuumratio = $vacuumratio"
  MD++ vacuumratio = $vacuumratio

  #Top surface atoms selection
  MD++ plot_color_axis = 2 eval # 2 means CSD
  MD++ input = \[ 1 -10 10 0.2 10 -10 10 10 100 \] fixatoms_by_pos_topol
  MD++ input = 1 setfixedatomsgroup 
  set N1 [MD++_Get "NPfixed"]
  MD++ freeallatoms
  #Bottom surface atoms selection
  MD++ input = \[ 1 -10 10 -10 -0.2 -10 10 10 100 \] fixatoms_by_pos_topol
  MD++ input = 2 setfixedatomsgroup 
  set N2 [MD++_Get "NPfixed"]
  MD++ freeallatoms

  # the number of atoms in each group should equal to the following
  #set N1 4800
  #set N2 4800

  set lx [MD++_Get "H_11"]
  set lz [MD++_Get "H_33"]

  #Compute Fext in terms of eV/A
  set f1x [expr -${sigma_xy}*$lx*$lz/$N1*0.6242E-05]
  set f1y [expr -${sigma_yy}*$lx*$lz/$N1*0.6242E-05]
  set f1z [expr -${sigma_yz}*$lx*$lz/$N1*0.6242E-05]

  set f2x [expr  ${sigma_xy}*$lx*$lz/$N1*0.6242E-05]
  set f2y [expr  ${sigma_yy}*$lx*$lz/$N1*0.6242E-05]
  set f2z [expr  ${sigma_yz}*$lx*$lz/$N1*0.6242E-05]

  MD++ extforce = \[ 2   $f1x $f1y $f1z  $f2x $f2y $f2z \]

  #MD++ stress = \[ 0 0 ${sigma_xz} 0 0 0 ${sigma_xz} 0 ${sigma_zz} \]  stressmul = 1.0
  MD++ eval

  # adjust box shape for sigma_xz and sigma_zz
  set maxiter 20
  set factor 1.0
  for { set iter 1 } { $iter <= $maxiter } { incr iter 1 } {
      if { $iter < [expr $maxiter/2] } {
              MD++ conj_ftol = 1.0
      } else {
              MD++ conj_ftol = 1e-4 
      }
      MD++ conj_itmax = 1000 conj_fevalmax = 1000
      # fix box
      MD++ conj_fixbox = 1  conj_fixboxvec = \[ 0 0 0   1 0 1   1 0 0 \]
      MD++ relax  eval
      set sig_xx [ MD++_Get TSTRESSinMPa_xx ] ; set sig_zz [ MD++_Get TSTRESSinMPa_zz ] ; 
      set sig_xz [ MD++_Get TSTRESSinMPa_xz ] ; 
      set e_xx [ expr ($sig_xx-$sigma_xx) / $modulus_xx ] ; 
      set e_zz [ expr ($sig_zz-$sigma_zz) / $modulus_zz ]
      set e_xz [ expr ($sig_xz-$sigma_xz) / $modulus_xz ] ; 

      puts "iter = $iter"
      puts "sig_xx = $sig_xx sig_zz = $sig_zz sig_xz = $sig_xz" 

      set H11_cur [ MD++_Get H_11 ] ; set H13_cur [ MD++_Get H_13 ] ; set H33_cur [ MD++_Get H_33 ] 

      set H11_new [ expr ${H11_cur}*(1.0+$e_xx*$factor) ] ;   MD++ H_11 = ${H11_new}
      set H33_new [ expr ${H33_cur}*(1.0+$e_zz*$factor) ] ;   MD++ H_33 = ${H33_new}
      set H13_new [ expr ${H13_cur} + ${H33_cur}*$e_xz*$factor ] ; MD++ H_13 = ${H13_new}

      MD++ finalcnfile = "ni-screw-gp-relax-$iter.cn" writecn
  }

  MD++ eval

  MD++ finalcnfile = "ni-screw-gp-relax-${n2}_${n3}_${n4}_${n5}.cn"  writecn
  MD++ finalcnfile = "ni-screw-gp-relax-${n2}_${n3}_${n4}_${n5}.cfg" writeatomeyecfg

  exitmd

} elseif { $status == 5.5 } {
 # relax screw dislocation structure A with stresses for Fleischer

 # glide plane normal: y-axis
 # dislocation line/Burger vector: x-axis
 
 # Schmid stress on glide plane: sigma_xy (must be 0)

 # Escaig stress on glide plane:      sigma_yz (applied by Fext)
 # Escaig stress on cross slip plane: sigma_zz (applied by H)
 # Escaig stress on cross slip plane: sigma_yy (applied by Fext)
 # Schmid stress on cross slip plane: sigma_xz (applied by H)

  set sigma_xx   0.0
  set sigma_xy   0.0
  set sigma_yz   $n2
  set sigma_zz   $n3
  set sigma_yy   $n4
  set sigma_xz   $n5

  set modulus_xx 450e3 ; #MPa
  set modulus_zz 450e3 ; #MPa
  set modulus_xz 100e3 ; #MPa

  #MD++ setnolog
  MD++ setoverwrite
  file mkdir "runs/ni-Rao-screw-30x30x20-$status"
  MD++ dirname = runs/ni-Rao-screw-30x30x20-$status/stress_${n2}_${n3}_${n4}_${n5}
  MD++ NNM = 200

  readpot
  MD++ incnfile = "../../ni-Rao-screw-30x30x20-0.5/perf.cn" readcn  writeall = 1
  set NP0 [MD++_Get "NP"]

  MD++ incnfile = "../../ni-Rao-screw-30x30x20-0.5/ni-screw-gp-relaxed.cn" readcn  writeall = 1

  # test new state A at Schmid stress 400 (does not relax well)
  #MD++ incnfile = "../../ni-Rao-screw-30x30x20-0.5/chain_no_7.cn" readcn freeallatoms writeall = 1

  setup_window
  #openwindow

  set NP [MD++_Get "NP"]
  set vacuumratio  [expr 1.0 - ($NP*1.0)/$NP0]
  puts "vacuumratio = $vacuumratio"
  MD++ vacuumratio = $vacuumratio

  #Top surface atoms selection
  MD++ plot_color_axis = 2 eval # 2 means CSD
  MD++ input = \[ 1 -10 10 0.2 10 -10 10 10 100 \] fixatoms_by_pos_topol
  MD++ input = 1 setfixedatomsgroup 
  set N1 [MD++_Get "NPfixed"]
  MD++ freeallatoms
  #Bottom surface atoms selection
  MD++ input = \[ 1 -10 10 -10 -0.2 -10 10 10 100 \] fixatoms_by_pos_topol
  MD++ input = 2 setfixedatomsgroup 
  set N2 [MD++_Get "NPfixed"]
  MD++ freeallatoms

  # the number of atoms in each group should equal to the following
  #set N1 4800
  #set N2 4800

  set lx [MD++_Get "H_11"]
  set lz [MD++_Get "H_33"]

  #Compute Fext in terms of eV/A
  set f1x [expr -${sigma_xy}*$lx*$lz/$N1*0.6242E-05]
  set f1y [expr -${sigma_yy}*$lx*$lz/$N1*0.6242E-05]
  set f1z [expr -${sigma_yz}*$lx*$lz/$N1*0.6242E-05]

  set f2x [expr  ${sigma_xy}*$lx*$lz/$N1*0.6242E-05]
  set f2y [expr  ${sigma_yy}*$lx*$lz/$N1*0.6242E-05]
  set f2z [expr  ${sigma_yz}*$lx*$lz/$N1*0.6242E-05]

  MD++ extforce = \[ 2   $f1x $f1y $f1z  $f2x $f2y $f2z \]

  #MD++ stress = \[ 0 0 ${sigma_xz} 0 0 0 ${sigma_xz} 0 ${sigma_zz} \]  stressmul = 1.0
  MD++ eval

  # adjust box shape for sigma_xz and sigma_zz
  set maxiter 20
  set factor 1.0
  for { set iter 1 } { $iter <= $maxiter } { incr iter 1 } {
      if { $iter < [expr $maxiter/2] } {
              MD++ conj_ftol = 1.0
      } else {
              MD++ conj_ftol = 1e-4 
      }
      MD++ conj_itmax = 1000 conj_fevalmax = 1000
      # fix box
      MD++ conj_fixbox = 1  conj_fixboxvec = \[ 0 0 0   1 0 1   1 0 0 \]
      MD++ relax  eval
      set sig_xx [ MD++_Get TSTRESSinMPa_xx ] ; set sig_zz [ MD++_Get TSTRESSinMPa_zz ] ; 
      set sig_xz [ MD++_Get TSTRESSinMPa_xz ] ; 
      set e_xx [ expr ($sig_xx-$sigma_xx) / $modulus_xx ] ; 
      set e_zz [ expr ($sig_zz-$sigma_zz) / $modulus_zz ]
      set e_xz [ expr ($sig_xz-$sigma_xz) / $modulus_xz ] ; 

      puts "iter = $iter"
      puts "sig_xx = $sig_xx sig_zz = $sig_zz sig_xz = $sig_xz" 

      set H11_cur [ MD++_Get H_11 ] ; set H13_cur [ MD++_Get H_13 ] ; set H33_cur [ MD++_Get H_33 ] 

      set H11_new [ expr ${H11_cur}*(1.0+$e_xx*$factor) ] ;   MD++ H_11 = ${H11_new}
      set H33_new [ expr ${H33_cur}*(1.0+$e_zz*$factor) ] ;   MD++ H_33 = ${H33_new}
      set H13_new [ expr ${H13_cur} + ${H33_cur}*$e_xz*$factor ] ; MD++ H_13 = ${H13_new}

      MD++ finalcnfile = "ni-screw-gp-relax-$iter.cn" writecn
  }

  MD++ eval

  MD++ finalcnfile = "ni-screw-gp-relax-${n2}_${n3}_${n4}_${n5}.cn"  writecn
  MD++ finalcnfile = "ni-screw-gp-relax-${n2}_${n3}_${n4}_${n5}.cfg" writeatomeyecfg

  exitmd

} elseif { $status == 6 } {
 # relax screw dislocation structure B_FE' with stresses

 # glide plane normal: y-axis
 # dislocation line/Burger vector: x-axis
 
 # Schmid stress on glide plane: sigma_xy (must be 0)

 # Escaig stress on glide plane:      sigma_yz (applied by Fext)
 # Escaig stress on cross slip plane: sigma_zz (applied by H)
 # Escaig stress on cross slip plane: sigma_yy (applied by Fext)
 # Schmid stress on cross slip plane: sigma_xz (applied by H)

  set sigma_xy   0.0
  set sigma_yz   $n2
  set sigma_zz   $n3
  set sigma_yy   $n4
  set sigma_xz   $n5

  set modulus_xz 100e3 ; #MPa

  # Try to choose a good starting point for the dislocation configuration
  set iselect 40

 if { $sigma_yy != 0 } {
  if { $sigma_yy <= 300 } {
     set iselect 40
     #set iselect 20
  } elseif { $sigma_yy <= 500 } {
     set iselect 12
  } else {
     #set iselect 9
     set iselect 5
  }
 }
  
 if { $sigma_yz != 0 } {
  if { $sigma_yz <= 300 } {
      if { $iselect > 40 } {
        set iselect 40
        #set iselect 20
      }
  } elseif { $sigma_yz <= 500 } {
      if { $iselect > 12 } {
        set iselect 12
      }
  } else {
      if { $iselect > 5 } {
        #set iselect 9
	set iselect 5
      }
  }
 }

 if { $sigma_xz != 0 } {
  if { $sigma_xz <= 300 } {
      if { $iselect > 40 } {
        set iselect 40
        #set iselect 20
      }
  } elseif { $sigma_xz <= 500 } {
      if { $iselect > 12 } {
        set iselect 12
      }
  } else {
      if { $iselect > 5 } {
        #set iselect 9
	set iselect 5
      }
  }
 }
 
   puts "iselect = $iselect"

  #MD++ setnolog
  MD++ setoverwrite
  # cannot create folder, create manually
  file mkdir "runs/ni-Rao-screw-30x30x20-$status"
  MD++ dirname = runs/ni-Rao-screw-30x30x20-$status/stress_${n2}_${n3}_${n4}_${n5}
  MD++ NNM = 200

  readpot

  # Read perfect crystal (for computing vacuumratio)
  MD++ incnfile = "../../ni-Rao-screw-30x30x20-0/perf.cn" readcn  writeall = 1
  set NP0 [MD++_Get "NP"]

  # Read state A (for its box and surrounding atoms, save coordinates to R0)
  #MD++ incnfile = "../../ni-Rao-screw-30x30x20-5/stress_${n2}_${n3}_${n4}_${n5}/ni-screw-gp-relax-${n2}_${n3}_${n4}_${n5}.cn" readcn saveH
  # use new A
  MD++ incnfile = "../../ni-Rao-screw-30x30x20-5/stress_${n2}_${n3}_${n4}_${n5}/ni-screw-gp-relax-1.cn" readcn saveH
  MD++ SHtoR RtoR0 writeall = 1
  
  # Read state B' (for its coordinates around disl)
  MD++ incnfile = "../../ni-Rao-screw-30x30x20-2/ni-cross-slip-relax-${iselect}.cn" readcn
  #MD++ incnfile = "../../ni-Rao-screw-30x30x20-0/chain_no_23.cn" readcn 
  MD++ restoreH  
  MD++ freeallatoms input = \[ 1 2 0 0 20 1 -10 10 \] makecylinder input = 3 setfixedatomsgroup 
  MD++ freeallatoms
  MD++ SHtoR input = \[ 1  3 \] R0toR_by_group RHtoS
  #MD++ R0toR RHtoS

  setup_window
  #openwindow

  set NP [MD++_Get "NP"]
  set vacuumratio  [expr 1.0 - ($NP*1.0)/$NP0]
  puts "vacuumratio = $vacuumratio"
  MD++ vacuumratio = $vacuumratio

  #Top surface atoms selection
  MD++ plot_color_axis = 2 eval # 2 means CSD
  MD++ input = \[ 1 -10 10 0.2 10 -10 10 10 100 \] fixatoms_by_pos_topol
  MD++ input = 1 setfixedatomsgroup 
  set N1 [MD++_Get "NPfixed"]
  MD++ freeallatoms
  #Bottom surface atoms selection
  MD++ input = \[ 1 -10 10 -10 -0.2 -10 10 10 100 \] fixatoms_by_pos_topol
  MD++ input = 2 setfixedatomsgroup 
  set N2 [MD++_Get "NPfixed"]
  MD++ freeallatoms

  # the number of atoms in each group should equal to the following
  #set N1 4800
  #set N2 4800

  set lx [MD++_Get "H_11"]
  set lz [MD++_Get "H_33"]

  #Compute Fext in terms of eV/A
  set f1x [expr -${sigma_xy}*$lx*$lz/$N1*0.6242E-05]
  set f1y [expr -${sigma_yy}*$lx*$lz/$N1*0.6242E-05]
  set f1z [expr -${sigma_yz}*$lx*$lz/$N1*0.6242E-05]

  set f2x [expr  ${sigma_xy}*$lx*$lz/$N1*0.6242E-05]
  set f2y [expr  ${sigma_yy}*$lx*$lz/$N1*0.6242E-05]
  set f2z [expr  ${sigma_yz}*$lx*$lz/$N1*0.6242E-05]

  MD++ extforce = \[ 2   $f1x $f1y $f1z  $f2x $f2y $f2z \]

  #MD++ stress = \[ 0 0 ${sigma_xz} 0 0 0 ${sigma_xz} 0 ${sigma_zz} \]  stressmul = 1.0
  MD++ eval
  for { set i 1 } { $i <= 20 } { incr i } {
      MD++ conj_ftol = 1e-4 conj_itmax = 1000 conj_fevalmax = 20
      # fix box
      MD++ conj_fixbox = 1  
      MD++ conj_fixboxvec = \[ 0 0 0   1 0 1   1 0 0 \]
      MD++ relax
      MD++ eval
      MD++ finalcnfile = "ni-cross-slip-relax-$i.cn" writecn
      MD++ finalcnfile = "ni-cross-slip-relax-$i.cfg" writeatomeyecfg
  }

  MD++ finalcnfile = "ni-cross-slip-relax-${n2}_${n3}_${n4}_${n5}.cn"  writecn
  MD++ finalcnfile = "ni-cross-slip-relax-${n2}_${n3}_${n4}_${n5}.cfg" writeatomeyecfg

  exitmd

} elseif { $status == 6.5 } {
 # relax screw dislocation structure B' with stresses
 # Fleischer Mechanism

 # glide plane normal: y-axis
 # dislocation line/Burger vector: x-axis
 
 # Schmid stress on glide plane: sigma_xy (must be 0)

 # Escaig stress on glide plane:      sigma_yz (applied by Fext)
 # Escaig stress on cross slip plane: sigma_zz (applied by H)
 # Escaig stress on cross slip plane: sigma_yy (applied by Fext)
 # Schmid stress on cross slip plane: sigma_xz (applied by H)

  set sigma_xy   0.0
  set sigma_yz   $n2
  set sigma_zz   $n3
  set sigma_yy   $n4
  set sigma_xz   $n5

  set modulus_xz 100e3 ; #MPa

 # Try to choose a good starting point for the dislocation configuration
 if { $sigma_yy != 0 } {
  if { $sigma_yy <= 300 } {
     set iselect 40
     #set iselect 20
  } elseif { $sigma_yy <= 500 } {
     set iselect 12
  } else {
     #set iselect 9
     set iselect 5
  }
 }
  
 if { $sigma_yz != 0 } {
  if { $sigma_yz <= 300 } {
     set iselect 40
     #set iselect 20
  } elseif { $sigma_yz <= 500 } {
     set iselect 12
  } else {
     #set iselect 9
	 set iselect 5
  }
 }

 if { $sigma_xz != 0 } {
  if { $sigma_xz <= 300 } {
     set iselect 40
     #set iselect 20
  } elseif { $sigma_xz <= 500 } {
     set iselect 12
  } else {
     #set iselect 9
     set iselect 5
  }
 }

# Dislocation on cross-slip plane is not visible in window by 5th CG step, so choosing a smaller iselect, at least initially
set iselect 1

  #MD++ setnolog
  MD++ setoverwrite
  # cannot create folder, create manually
  file mkdir "runs/ni-Rao-screw-30x30x20-$status"
  MD++ dirname = runs/ni-Rao-screw-30x30x20-$status/stress_${n2}_${n3}_${n4}_${n5}
  MD++ NNM = 200

  readpot

  # Read perfect crystal (for computing vacuumratio)
  MD++ incnfile = "../../ni-Rao-screw-30x30x20-0.5/perf.cn" readcn  writeall = 1
  set NP0 [MD++_Get "NP"]

  # Read state A (for its box and surrounding atoms, save coordinates to R0)
  #MD++ incnfile = "../../ni-Rao-screw-30x30x20-5.5/stress_${n2}_${n3}_${n4}_${n5}/ni-screw-gp-relax-${n2}_${n3}_${n4}_${n5}.cn" readcn saveH
  # use new A
  MD++ incnfile = "../../ni-Rao-screw-30x30x20-5.5/stress_${n2}_${n3}_${n4}_${n5}/ni-screw-gp-relax-1.cn" readcn saveH
  MD++ SHtoR RtoR0 writeall = 1
  
  # Read state B' (for its coordinates around disl)
  MD++ incnfile = "../../ni-Rao-screw-30x30x20-2.5/ni-cross-slip-relax-${iselect}.cn" readcn
  #MD++ incnfile = "../../ni-Rao-screw-30x30x20-0.5/chain_no_23.cn" readcn 
  MD++ restoreH  
  MD++ freeallatoms input = \[ 1 2 0 0 20 1 -10 10 \] makecylinder input = 3 setfixedatomsgroup 
  MD++ freeallatoms
  MD++ SHtoR input = \[ 1  3 \] R0toR_by_group RHtoS
  #MD++ R0toR RHtoS

  setup_window
  #openwindow

  set NP [MD++_Get "NP"]
  set vacuumratio  [expr 1.0 - ($NP*1.0)/$NP0]
  puts "vacuumratio = $vacuumratio"
  MD++ vacuumratio = $vacuumratio

  #Top surface atoms selection
  MD++ plot_color_axis = 2 eval # 2 means CSD
  MD++ input = \[ 1 -10 10 0.2 10 -10 10 10 100 \] fixatoms_by_pos_topol
  MD++ input = 1 setfixedatomsgroup 
  set N1 [MD++_Get "NPfixed"]
  MD++ freeallatoms
  #Bottom surface atoms selection
  MD++ input = \[ 1 -10 10 -10 -0.2 -10 10 10 100 \] fixatoms_by_pos_topol
  MD++ input = 2 setfixedatomsgroup 
  set N2 [MD++_Get "NPfixed"]
  MD++ freeallatoms

  # the number of atoms in each group should equal to the following
  #set N1 4800
  #set N2 4800

  set lx [MD++_Get "H_11"]
  set lz [MD++_Get "H_33"]

  #Compute Fext in terms of eV/A
  set f1x [expr -${sigma_xy}*$lx*$lz/$N1*0.6242E-05]
  set f1y [expr -${sigma_yy}*$lx*$lz/$N1*0.6242E-05]
  set f1z [expr -${sigma_yz}*$lx*$lz/$N1*0.6242E-05]

  set f2x [expr  ${sigma_xy}*$lx*$lz/$N1*0.6242E-05]
  set f2y [expr  ${sigma_yy}*$lx*$lz/$N1*0.6242E-05]
  set f2z [expr  ${sigma_yz}*$lx*$lz/$N1*0.6242E-05]

  MD++ extforce = \[ 2   $f1x $f1y $f1z  $f2x $f2y $f2z \]

  #MD++ stress = \[ 0 0 ${sigma_xz} 0 0 0 ${sigma_xz} 0 ${sigma_zz} \]  stressmul = 1.0
  MD++ eval
  for { set i 1 } { $i <= 20 } { incr i } {
      MD++ conj_ftol = 1e-4 conj_itmax = 1000 conj_fevalmax = 20
      # fix box
      MD++ conj_fixbox = 1  
      MD++ conj_fixboxvec = \[ 0 0 0   1 0 1   1 0 0 \]
      MD++ relax
      MD++ eval
      MD++ finalcnfile = "ni-cross-slip-relax-$i.cn" writecn
      MD++ finalcnfile = "ni-cross-slip-relax-$i.cfg" writeatomeyecfg
  }

  MD++ finalcnfile = "ni-cross-slip-relax-${n2}_${n3}_${n4}_${n5}.cn"  writecn
  MD++ finalcnfile = "ni-cross-slip-relax-${n2}_${n3}_${n4}_${n5}.cfg" writeatomeyecfg

  exitmd

} elseif { $status == 7 } {
  # visualize dislocation structure
  #  select atoms on one plane to be analysed for cross slip (status == 7)
  MD++ setnolog
  initmd $status
  readpot

  # set atoms on cross slip plane to group 950129

  MD++ incnfile = "../ni-Rao-screw-30x30x20-6/stress_${n2}_${n3}_${n4}_${n5}/ni-cross-slip-relax-${n2}_${n3}_${n4}_${n5}.cn" readcn
  
  MD++ latticeconst = 3.52 #(A) for Ni
  set a0 3.52
  set nx 0
  set ny [expr 1.0/3.0]
  set nz [expr -sqrt(8.0)/3.0]
  MD++ SHtoR
  set NP [MD++_Get "NP"]

  puts "find atoms on cross slip plane..."
  for { set i 0 } { $i < $NP } { incr i } { 
     set ind [expr $i*3]
     set rx [MD++_Get R($ind)]
     set ind [expr $i*3+1]
     set ry [MD++_Get R($ind)]
     set ind [expr $i*3+2]
     set rz [MD++_Get R($ind)]
     set rdotn [expr $rx*$nx + $ry*$ny + $rz*$nz]
     set temp [expr $rdotn/($a0*sqrt(3.0)/4.0)]
     if { $temp > -0.5 && $temp < 0.5 } {
         MD++_Set fixed($i)  1
     }
  }
  
  openwindow
  MD++ eval

  MD++ sleep
  exitmd
  
} elseif { $status == 7.5 } {
  # visualize dislocation structure
  #  select atoms on one plane to be analysed for cross slip (status == 3)
  MD++ setnolog
  initmd $status
  readpot

  # set atoms on cross slip plane to group 950129

  #MD++ incnfile = "../ni-Rao-screw-30x30x20-6.5/stress_${n2}_${n3}_${n4}_${n5}/ni-cross-slip-relax-${n2}_${n3}_${n4}_${n5}.cn" readcn
  MD++ incnfile = "../ni-Rao-screw-30x30x20-5.5/stress_${n2}_${n3}_${n4}_${n5}/ni-screw-gp-relax-1.cn" readcn
  MD++ incnfile = "../ni-Rao-screw-30x30x20-6.5/stress_${n2}_${n3}_${n4}_${n5}/ni-cross-slip-relax-1.cn" readcn
  #MD++ incnfile = "../ni-Rao-cs-112.5_${n}_${n2}_${n3}_${n4}_${n5}/NEBfinal.cn" readcn
  
  MD++ latticeconst = 3.52 #(A) for Ni
  set a0 3.52
  set nx 0
  set ny [expr 1.0/3.0]
  set nz [expr -sqrt(8.0)/3.0]
  MD++ SHtoR
  set NP [MD++_Get "NP"]

  puts "find atoms on cross slip plane..."
if {0} {
   for { set i 0 } { $i < $NP } { incr i } { 
     set ind [expr $i*3]
     set rx [MD++_Get R($ind)]
     set ind [expr $i*3+1]
     set ry [MD++_Get R($ind)]
     set ind [expr $i*3+2]
     set rz [MD++_Get R($ind)]
     set rdotn [expr $rx*$nx + $ry*$ny + $rz*$nz]
     set temp [expr $rdotn/($a0*sqrt(3.0)/4.0)]
     if { $temp > -0.5 && $temp < 0.5 } {
         MD++_Set fixed($i)  1
     }
  }
}

  #MD++ incnfile = "../ni-Rao-screw-30x30x20-5.5/stress_${n2}_${n3}_${n4}_${n5}/ni-screw-gp-relax-1.cn" readcn
  ##MD++ incnfile = "../ni-Rao-cs-112.5_${n}_${n2}_${n3}_${n4}_${n5}/NEBinit.cn" readcn
  openwindow

#     MD++ NCS = 12 eval calcentralsymmetry 
#     MD++ eval plot 

  MD++ incnfile = "../ni-Rao-screw-30x30x20-6.5/stress_${n2}_${n3}_${n4}_${n5}/ni-cross-slip-relax-1.cn" readcn
  ##MD++ incnfile = "../ni-Rao-cs-112.5_${n}_${n2}_${n3}_${n4}_${n5}/NEBfinal.cn" readcn
  #MD++ eval
     MD++ NCS = 12 eval calcentralsymmetry
     MD++ eval plot

  MD++ sleep
  exitmd
  
} elseif { $status == 11 } {
  # Parallel NEB relaxation (for initially created state A and B)
  # to compile
  #  make eam build=R SYS=mc2_mpich SPEC+=-DMAXCONSTRAINATOMS=800001
  # to run (in PBS file)
  #  mpiexec -np $ncpu bin1/eam_mc2_mpich scripts/work/ni_disl/ni_screw_cross_slip 11 5000 0 0 0

  set NEBMAXSTEP $n
  set EQUILSTEP  [expr $NEBMAXSTEP - 100 ]
  if { $EQUILSTEP < 0 } {
     set EQUILSTEP 0
  }

  set sigma_xy   0.0
  set sigma_yz   $n2
  set sigma_zz   $n3
  set sigma_yy   $n4
  set sigma_xz   $n5

  #set sigma_e_g  $n2
  #set sigma_e_c  $n3
  #set sigma_s_c  $n4
  #set sigma_s_c  $n5

  set iselect 10
  #set iselect 30

  set ncpu [ MD++_Get numDomains ]
  if { $ncpu < 1 } { set ncpu 1 }
  puts "ncpu = $ncpu"
  
  MD++ setoverwrite 
  MD++ dirname = runs/ni-Rao-cs-${status}_${n}_${n2}_${n3}_${n4}_${n5}
  MD++ NNM = 200
  exec ln -f -s ../ni-Rao-screw-30x30x20-0/ni-screw-gp-relaxed.cn   NEBinit.cn
  #exec ln -f -s ../ni-Rao-screw-30x30x20-2/ni-cross-slip-relaxed.cn NEBfinal.cn
  exec ln -f -s "../ni-Rao-screw-30x30x20-2/ni-cross-slip-relax-${iselect}.cn" NEBfinal.cn


  MD++ writeall = 1

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

  readpot
  MD++ Broadcast_EAM_Param

  MD++ incnfile = "../ni-Rao-screw-30x30x20-0/perf.cn" readcn  writeall = 1
  set NP0 [MD++_Get "NP"]

  MD++ incnfile = NEBinit.cn  readcn setconfig1 #A
  MD++ incnfile = NEBfinal.cn readcn setconfig2 #B
  MD++ incnfile = NEBinit.cn  readcn            #A

  MD++ freeallatoms

  setup_window

  set NP [MD++_Get "NP"]
  set vacuumratio  [expr 1.0 - ($NP*1.0)/$NP0]
  puts "vacuumratio = $vacuumratio"
  MD++ vacuumratio = $vacuumratio

  MD++ fixallatoms constrain_fixedatoms freeallatoms
  MD++ chainlength = [expr $ncpu-1] 

  MD++ timestep = 0.07 printfreq = 2
  MD++ nebspec = \[ 0 1 0 1 1 0 \] totalsteps = $NEBMAXSTEP equilsteps = $EQUILSTEP

  #MD++ timestep = 0.1 printfreq = 2
  #MD++ nebspec = \[ 0 1 0 1 1 1 \] totalsteps = $NEBMAXSTEP equilsteps = $EQUILSTEP
  ## [ relax_surround redistr_freq moveleftend moverightend yesclimbimage islengthconstant ]

  MD++ Slave_chdir #Slave change to same directory as Master (dirname)
  MD++ allocchain_parallel initRchain_parallel

  if { $exist > 0 } {
    MD++ incnfile = stringrelax.chain.cn readRchain_parallel
  }

  MD++ stringrelax_parallel 
  #MD++ stringrelax_parallel_1

  MD++ finalcnfile = stringrelax.chain.cn writeRchain_parallel
  MD++ quit_all

# Visualize results
#  octave
#   data=load('stringeng.out'); figure(1); plot(data(end-10:end,3:5:end)', '*-');
#   Emax = max(data(:,3:5:end),[],2); figure(2); plot(Emax,'*-')

} elseif { $status == 11.5 } {
  # Parallel NEB relaxation (for initially created state A and B)
  # Fleischer Mechanism
  # to compile
  #  make eam build=R SYS=mc2_mpich SPEC+=-DMAXCONSTRAINATOMS=800001
  # to run (in PBS file)
  #  mpiexec -np $ncpu bin1/eam_mc2_mpich scripts/work/ni_disl/ni_screw_cross_slip 11 5000 0 0 0

  set NEBMAXSTEP $n
  set EQUILSTEP  [expr $NEBMAXSTEP - 100 ]
  if { $EQUILSTEP < 0 } {
     set EQUILSTEP 0
  }

  set sigma_xy   0.0
  set sigma_yz   $n2
  set sigma_zz   $n3
  set sigma_yy   $n4
  set sigma_xz   $n5

  #set sigma_e_g  $n2
  #set sigma_e_c  $n3
  #set sigma_s_c  $n4
  #set sigma_s_c  $n5

  set iselect 10
  #set iselect 30

  set ncpu [ MD++_Get numDomains ]
  if { $ncpu < 1 } { set ncpu 1 }
  puts "ncpu = $ncpu"
  
  MD++ setoverwrite 
  MD++ dirname = runs/ni-Rao-cs-${status}_${n}_${n2}_${n3}_${n4}_${n5}
  MD++ NNM = 200
  exec ln -f -s ../ni-Rao-screw-30x30x20-0.5/ni-screw-gp-relaxed.cn   NEBinit.cn
  #exec ln -f -s ../ni-Rao-screw-30x30x20-2.5/ni-cross-slip-relaxed.cn NEBfinal.cn
  exec ln -f -s "../ni-Rao-screw-30x30x20-2.5/ni-cross-slip-relax-${iselect}.cn" NEBfinal.cn


  MD++ writeall = 1

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

  readpot
  MD++ Broadcast_EAM_Param

  MD++ incnfile = "../ni-Rao-screw-30x30x20-0.5/perf.cn" readcn  writeall = 1
  set NP0 [MD++_Get "NP"]

  MD++ incnfile = NEBinit.cn  readcn setconfig1 #A
  MD++ incnfile = NEBfinal.cn readcn setconfig2 #B
  MD++ incnfile = NEBinit.cn  readcn            #A

  MD++ freeallatoms

  setup_window

  set NP [MD++_Get "NP"]
  set vacuumratio  [expr 1.0 - ($NP*1.0)/$NP0]
  puts "vacuumratio = $vacuumratio"
  MD++ vacuumratio = $vacuumratio

  MD++ fixallatoms constrain_fixedatoms freeallatoms
  MD++ chainlength = [expr $ncpu-1] 

  MD++ timestep = 0.07 printfreq = 2
  MD++ nebspec = \[ 0 1 0 1 1 0 \] totalsteps = $NEBMAXSTEP equilsteps = $EQUILSTEP

  #MD++ timestep = 0.1 printfreq = 2
  #MD++ nebspec = \[ 0 1 0 1 1 1 \] totalsteps = $NEBMAXSTEP equilsteps = $EQUILSTEP
  ## [ relax_surround redistr_freq moveleftend moverightend yesclimbimage islengthconstant ]

  MD++ Slave_chdir #Slave change to same directory as Master (dirname)
  MD++ allocchain_parallel initRchain_parallel

  if { $exist > 0 } {
    MD++ incnfile = stringrelax.chain.cn readRchain_parallel
  }

  MD++ stringrelax_parallel 
  #MD++ stringrelax_parallel_1

  MD++ finalcnfile = stringrelax.chain.cn writeRchain_parallel
  MD++ quit_all

# Visualize results
#  octave
#   data=load('stringeng.out'); figure(1); plot(data(end-10:end,3:5:end)', '*-');
#   Emax = max(data(:,3:5:end),[],2); figure(2); plot(Emax,'*-')

} elseif { $status == 12 } {
  # Parallel NEB relaxation (for structures relaxed under stress)
  # to compile
  #  make eam build=R SYS=mc2_mpich SPEC+=-DMAXCONSTRAINATOMS=800001
  # to run (in PBS file)
  #  mpiexec -np $ncpu bin1/eam_mc2_mpich scripts/work/ni_disl/ni_screw_cross_slip 12 5000 0 0 0

  set NEBMAXSTEP $n
  set EQUILSTEP  [expr $NEBMAXSTEP - 100 ]
  if { $EQUILSTEP < 0 } {
     set EQUILSTEP 0
  }

  set sigma_xy   0.0
  set sigma_yz   $n2
  set sigma_zz   $n3
  set sigma_yy   $n4
  set sigma_xz   $n5

 if { $sigma_yy != 0 } {
  if { $sigma_yy <= 100 } {
     set iselect 13
  } elseif { $sigma_yy <= 300 } {
     set iselect 10
  } elseif { $sigma_yy <= 600 } {
     set iselect 10
  } else {
     set iselect 6
  }
 }
  # Escaig stress
 if { $sigma_yz != 0 } {
  if { $sigma_yz <= 300 } {
     set iselect 20
  } elseif { $sigma_yz <= 500 } {
     set iselect 10
  } elseif { $sigma_yz <= 600 } {
     set iselect 8
  } else {
     set iselect 5
  }
 }
  # Schmid stress
 if { $sigma_xz != 0 } {
  if { $sigma_xz <= 100 } {
     set iselect 13
  } elseif { $sigma_xz <= 300 } {
     set iselect 11
  } elseif { $sigma_xz <= 600 } {
     set iselect 10
  } else {
     set iselect 5
  }
 }

  set ncpu [ MD++_Get numDomains ]
  if { $ncpu < 1 } { set ncpu 1 }
  puts "ncpu = $ncpu"
  
  MD++ setoverwrite 
  MD++ dirname = runs/ni-Rao-cs-${status}_${n}_${n2}_${n3}_${n4}_${n5}

  MD++ NNM = 200
  MD++ writeall = 1

  exec ln -f -s ../ni-Rao-screw-30x30x20-5/stress_${n2}_${n3}_${n4}_${n5}/ni-screw-gp-relax-${n2}_${n3}_${n4}_${n5}.cn NEBinit.cn
  exec ln -f -s ../ni-Rao-screw-30x30x20-6/stress_${n2}_${n3}_${n4}_${n5}/ni-cross-slip-relax-${iselect}.cn NEBfinal.cn

  # use new A and B'
  #exec ln -f -s chain_no_7.cn NEBinit.cn
  #exec ln -f -s chain_no_23.cn NEBfinal.cn


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

  readpot
  MD++ Broadcast_EAM_Param

  MD++ incnfile = "../ni-Rao-screw-30x30x20-0/perf.cn" readcn  writeall = 1
  set NP0 [MD++_Get "NP"]

  # a test for new state A (0_0_0_400)
  #set chain_no 6
  #MD++ incnfile = NEBinit.cn  readcn 
  #MD++ incnfile = "../ni-cs-12_${n}_${n2}_${n3}_${n4}_${n5}.sav/stringrelax.chain.cn"
  #MD++ input = $chain_no  readRchain_parallel_toCN  RHtoS  setconfig1 #A

  #MD++ incnfile = NEBfinal.cn readcn setconfig2                       #B
  #MD++ incnfile = "../ni-cs-12_${n}_${n2}_${n3}_${n4}_${n5}.sav/stringrelax.chain.cn"
  #MD++ input = $chain_no  readRchain_parallel_toCN  RHtoS             #A

  MD++ incnfile = NEBinit.cn  readcn setconfig1 #A
  MD++ incnfile = NEBfinal.cn readcn setconfig2 #B
  MD++ incnfile = NEBinit.cn  readcn            #A

  MD++ freeallatoms

  setup_window

  set NP [MD++_Get "NP"]
  set vacuumratio  [expr 1.0 - ($NP*1.0)/$NP0]
  puts "vacuumratio = $vacuumratio"
  MD++ vacuumratio = $vacuumratio

  if { 1 } { # determine number of group 1/2 atoms 
  #Top surface atoms selection
  MD++ plot_color_axis = 2 eval # 2 means CSD
  MD++ input = \[ 1 -10 10 0.2 10 -10 10 10 100 \] fixatoms_by_pos_topol
  MD++ input = 1 setfixedatomsgroup 
  set N1 [MD++_Get "NPfixed"]
  MD++ freeallatoms
  #Bottom surface atoms selection
  MD++ input = \[ 1 -10 10 -10 -0.2 -10 10 10 100 \] fixatoms_by_pos_topol
  MD++ input = 2 setfixedatomsgroup 
  set N2 [MD++_Get "NPfixed"]
  MD++ freeallatoms
  } else {
  # the number of atoms in each group should equal to the following
    set N1 4800
    set N2 4800
  }

  set lx [MD++_Get "H_11"]
  set lz [MD++_Get "H_33"]

  #Compute Fext in terms of eV/A
  set f1x [expr -${sigma_xy}*$lx*$lz/$N1*0.6242E-05]
  set f1y [expr -${sigma_yy}*$lx*$lz/$N1*0.6242E-05]
  set f1z [expr -${sigma_yz}*$lx*$lz/$N1*0.6242E-05]

  set f2x [expr  ${sigma_xy}*$lx*$lz/$N1*0.6242E-05]
  set f2y [expr  ${sigma_yy}*$lx*$lz/$N1*0.6242E-05]
  set f2z [expr  ${sigma_yz}*$lx*$lz/$N1*0.6242E-05]

  MD++ extforce = \[ 2   $f1x $f1y $f1z  $f2x $f2y $f2z \]

  #MD++ freeallatoms clearR0 refreshnnlist eval  quit_all

  # do we need to broadcast extforce ?
  MD++ Broadcast_Atoms

  MD++ fixallatoms constrain_fixedatoms freeallatoms
  MD++ chainlength = [expr $ncpu-1] 

  MD++ timestep = 0.07 printfreq = 2
  #MD++ timestep = 0.1 printfreq = 2

  # make left end not movable
  #MD++ nebspec = \[ 0 1 0 1 1 0 \] totalsteps = $NEBMAXSTEP equilsteps = $EQUILSTEP

  set NEBSTEP1 1000

  # set path length constant
  #MD++ nebspec = \[ 0 1 0 1 1 2 \] totalsteps = $NEBMAXSTEP equilsteps = $EQUILSTEP
  MD++ nebspec = \[ 0 1 0 1 1 1 \] totalsteps = $NEBMAXSTEP equilsteps = $EQUILSTEP
  #MD++ nebspec = \[ 0 1 0 1 1 1 \] totalsteps = $NEBSTEP1 equilsteps = $EQUILSTEP
  ## [ relax_surround redistr_freq moveleftend moverightend yesclimbimage islengthconstant ]

  MD++ Slave_chdir #Slave change to same directory as Master (dirname)
  MD++ allocchain_parallel initRchain_parallel

  if { $exist > 0 } {
    MD++ incnfile = stringrelax.chain.cn readRchain_parallel
  }

  MD++ eval
  MD++ stringrelax_parallel_2
  MD++ finalcnfile = stringrelax.chain.cn writeRchain_parallel
  set check [ glob -nocomplain "stringeng_step*" ]
  set nfile [ llength $check]
  exec cp stringeng.out stringeng_step[expr $nfile+1].out

  #MD++ nebspec = \[ 0 1 0 1 1 0 \] totalsteps = $NEBMAXSTEP equilsteps = $EQUILSTEP
  MD++ nebspec = \[ 0 1 0 1 1 1 \] totalsteps = $NEBSTEP1 equilsteps = $EQUILSTEP
  MD++ Broadcast_nebsetting
  MD++ stringrelax_parallel
  #MD++ stringrelax_parallel_1
  MD++ finalcnfile = stringrelax.chain.cn writeRchain_parallel
  exec cp stringeng.out stringeng_step[expr $nfile+2].out

  MD++ quit_all

} elseif { $status == 12.5 } {
  # Parallel NEB relaxation (for structures relaxed under stress)
  # Fleischer Mechanism
  # to compile
  #  make eam build=R SYS=mc2_mpich SPEC+=-DMAXCONSTRAINATOMS=800001
  # to run (in PBS file)
  #  mpiexec -np $ncpu bin1/eam_mc2_mpich scripts/work/ni_disl/ni_screw_cross_slip 12.5 5000 0 0 0

  set NEBMAXSTEP $n
  set EQUILSTEP  [expr $NEBMAXSTEP - 100 ]
  if { $EQUILSTEP < 0 } {
     set EQUILSTEP 0
  }

  set sigma_xy   0.0
  set sigma_yz   $n2
  set sigma_zz   $n3
  set sigma_yy   $n4
  set sigma_xz   $n5

 if { $sigma_yy != 0 } {
  if { $sigma_yy <= 100 } {
     set iselect 13
  } elseif { $sigma_yy <= 300 } {
     set iselect 10
  } elseif { $sigma_yy <= 600 } {
     set iselect 10
  } else {
     set iselect 6
  }
 }
  # Escaig stress
 if { $sigma_yz != 0 } {
  if { $sigma_yz <= 300 } {
     set iselect 20
  } elseif { $sigma_yz <= 500 } {
     set iselect 10
  } elseif { $sigma_yz <= 600 } {
     set iselect 8
  } else {
     set iselect 5
  }
 }
  # Schmid stress
 if { $sigma_xz != 0 } {
  if { $sigma_xz <= 100 } {
     set iselect 13
  } elseif { $sigma_xz <= 300 } {
     set iselect 11
  } elseif { $sigma_xz <= 600 } {
     set iselect 10
  } else {
     set iselect 5
  }
 }

  set ncpu [ MD++_Get numDomains ]
  if { $ncpu < 1 } { set ncpu 1 }
  puts "ncpu = $ncpu"
  
  MD++ setoverwrite 
  MD++ dirname = runs/ni-Rao-cs-${status}_${n}_${n2}_${n3}_${n4}_${n5}

  MD++ NNM = 200
  MD++ writeall = 1

  exec ln -f -s ../ni-Rao-screw-30x30x20-5.5/stress_${n2}_${n3}_${n4}_${n5}/ni-screw-gp-relax-${n2}_${n3}_${n4}_${n5}.cn NEBinit.cn
  exec ln -f -s ../ni-Rao-screw-30x30x20-6.5/stress_${n2}_${n3}_${n4}_${n5}/ni-cross-slip-relax-${iselect}.cn NEBfinal.cn

  # use new A and B'
  #exec ln -f -s chain_no_7.cn NEBinit.cn
  #exec ln -f -s chain_no_23.cn NEBfinal.cn


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

  readpot
  MD++ Broadcast_EAM_Param

  MD++ incnfile = "../ni-Rao-screw-30x30x20-0.5/perf.cn" readcn  writeall = 1
  set NP0 [MD++_Get "NP"]

  # a test for new state A (0_0_0_400)
  #set chain_no 6
  #MD++ incnfile = NEBinit.cn  readcn 
  #MD++ incnfile = "../ni-cs-12_${n}_${n2}_${n3}_${n4}_${n5}.sav/stringrelax.chain.cn"
  #MD++ input = $chain_no  readRchain_parallel_toCN  RHtoS  setconfig1 #A

  #MD++ incnfile = NEBfinal.cn readcn setconfig2                       #B
  #MD++ incnfile = "../ni-cs-12_${n}_${n2}_${n3}_${n4}_${n5}.sav/stringrelax.chain.cn"
  #MD++ input = $chain_no  readRchain_parallel_toCN  RHtoS             #A

  MD++ incnfile = NEBinit.cn  readcn setconfig1 #A
  MD++ incnfile = NEBfinal.cn readcn setconfig2 #B
  MD++ incnfile = NEBinit.cn  readcn            #A

  MD++ freeallatoms

  setup_window

  set NP [MD++_Get "NP"]
  set vacuumratio  [expr 1.0 - ($NP*1.0)/$NP0]
  puts "vacuumratio = $vacuumratio"
  MD++ vacuumratio = $vacuumratio

  if { 1 } { # determine number of group 1/2 atoms 
  #Top surface atoms selection
  MD++ plot_color_axis = 2 eval # 2 means CSD
  MD++ input = \[ 1 -10 10 0.2 10 -10 10 10 100 \] fixatoms_by_pos_topol
  MD++ input = 1 setfixedatomsgroup 
  set N1 [MD++_Get "NPfixed"]
  MD++ freeallatoms
  #Bottom surface atoms selection
  MD++ input = \[ 1 -10 10 -10 -0.2 -10 10 10 100 \] fixatoms_by_pos_topol
  MD++ input = 2 setfixedatomsgroup 
  set N2 [MD++_Get "NPfixed"]
  MD++ freeallatoms
  } else {
  # the number of atoms in each group should equal to the following
    set N1 4800
    set N2 4800
  }

  set lx [MD++_Get "H_11"]
  set lz [MD++_Get "H_33"]

  #Compute Fext in terms of eV/A
  set f1x [expr -${sigma_xy}*$lx*$lz/$N1*0.6242E-05]
  set f1y [expr -${sigma_yy}*$lx*$lz/$N1*0.6242E-05]
  set f1z [expr -${sigma_yz}*$lx*$lz/$N1*0.6242E-05]

  set f2x [expr  ${sigma_xy}*$lx*$lz/$N1*0.6242E-05]
  set f2y [expr  ${sigma_yy}*$lx*$lz/$N1*0.6242E-05]
  set f2z [expr  ${sigma_yz}*$lx*$lz/$N1*0.6242E-05]

  MD++ extforce = \[ 2   $f1x $f1y $f1z  $f2x $f2y $f2z \]

  #MD++ freeallatoms clearR0 refreshnnlist eval  quit_all

  # do we need to broadcast extforce ?
  MD++ Broadcast_Atoms

  MD++ fixallatoms constrain_fixedatoms freeallatoms
  MD++ chainlength = [expr $ncpu-1] 

  MD++ timestep = 0.07 printfreq = 2
  #MD++ timestep = 0.1 printfreq = 2

  # make left end not movable
  #MD++ nebspec = \[ 0 1 0 1 1 0 \] totalsteps = $NEBMAXSTEP equilsteps = $EQUILSTEP

  set NEBSTEP1 1000

  # set path length constant
  #MD++ nebspec = \[ 0 1 0 1 1 2 \] totalsteps = $NEBMAXSTEP equilsteps = $EQUILSTEP
  MD++ nebspec = \[ 0 1 0 1 1 1 \] totalsteps = $NEBMAXSTEP equilsteps = $EQUILSTEP
  #MD++ nebspec = \[ 0 1 0 1 1 1 \] totalsteps = $NEBSTEP1 equilsteps = $EQUILSTEP
  ## [ relax_surround redistr_freq moveleftend moverightend yesclimbimage islengthconstant ]

  MD++ Slave_chdir #Slave change to same directory as Master (dirname)
  MD++ allocchain_parallel initRchain_parallel

  if { $exist > 0 } {
    MD++ incnfile = stringrelax.chain.cn readRchain_parallel
  }

  MD++ eval
  MD++ stringrelax_parallel_1
  MD++ finalcnfile = stringrelax.chain.cn writeRchain_parallel
  set check [ glob -nocomplain "stringeng_step*" ]
  set nfile [ llength $check]
  exec cp stringeng.out stringeng_step[expr $nfile+1].out

  #MD++ nebspec = \[ 0 1 0 1 1 0 \] totalsteps = $NEBMAXSTEP equilsteps = $EQUILSTEP
  MD++ nebspec = \[ 0 1 0 1 1 1 \] totalsteps = $NEBSTEP1 equilsteps = $EQUILSTEP
  MD++ Broadcast_nebsetting
  MD++ stringrelax_parallel
  #MD++ stringrelax_parallel_1
  MD++ finalcnfile = stringrelax.chain.cn writeRchain_parallel
  exec cp stringeng.out stringeng_step[expr $nfile+2].out

  MD++ quit_all

} elseif { $status == 13 } {
  # Parallel NEB relaxation (for structures relaxed under stress)
  # using relaxed path at a different stress as initial path
  # to compile
  #  make eam build=R SYS=mc2_mpich SPEC+=-DMAXCONSTRAINATOMS=800001
  # to run (in PBS file)
  #  mpiexec -np $ncpu bin1/eam_mc2_mpich scripts/work/ni_disl/ni_screw_cross_slip 13 5000 0 0 0

  set NEBMAXSTEP $n
  set EQUILSTEP  [expr $NEBMAXSTEP - 100 ]
  if { $EQUILSTEP < 0 } {
     set EQUILSTEP 0
  }

  set sigma_xy   0.0
  set sigma_yz   $n2
  set sigma_zz   $n3
  set sigma_yy   $n4
  set sigma_xz   $n5

  set n5prev   400

  set iselect 6

  set rerun 1

  set ncpu [ MD++_Get numDomains ]
  if { $ncpu < 1 } { set ncpu 1 }
  puts "ncpu = $ncpu"
  
  MD++ setoverwrite 
  MD++ dirname = runs/ni-Rao-cs-${status}_${n}_${n2}_${n3}_${n4}_${n5}

  MD++ NNM = 200
  MD++ writeall = 1

  exec ln -f -s ../ni-Rao-screw-30x30x20-5/stress_${n2}_${n3}_${n4}_${n5}/ni-screw-gp-relax-${n2}_${n3}_${n4}_${n5}.cn NEBinit.cn
  exec ln -f -s ../ni-Rao-screw-30x30x20-6/stress_${n2}_${n3}_${n4}_${n5}/ni-cross-slip-relax-${iselect}.cn NEBfinal.cn

 if { $rerun } {
  # rerun simulation using previously relaxed path at same stress
  set check [ glob -nocomplain "stringrelax.chain.cn.cpu00" ]
  set exist [ llength $check ]
  if { $exist > 0 } {
     for { set i 0 } { $i < $ncpu } { incr i 1 } {
         set ival [ format %02d $i ]
         exec cp stringrelax.chain.cn.cpu$ival Prev_stringrelax.chain.cn.cpu$ival
     }
     exec cp stringeng.out Prev_stringeng.out
  }
 } else {
  # relax path using a previously relaxed path at lower stress
  set check [ glob -nocomplain "../ni-cs-12_${n}_${n2}_${n3}_${n4}_${n5prev}/stringrelax.chain.cn.cpu00" ]
  set exist [ llength $check ]
  if { $exist > 0 } {
     for { set i 0 } { $i < $ncpu } { incr i 1 } {
         set ival [ format %02d $i ]
         exec cp ../ni-cs-12_${n}_${n2}_${n3}_${n4}_${n5prev}/stringrelax.chain.cn.cpu$ival Prev_stringrelax.chain.cn.cpu$ival
     }
     #exec cp stringeng.out Prev_stringeng.out
     exec cp ../ni-cs-12_${n}_${n2}_${n3}_${n4}_${n5prev}/NEBinit.cn Prev_NEBinit.cn
  }
 }

  readpot
  MD++ Broadcast_EAM_Param

  MD++ incnfile = "../ni-Rao-screw-30x30x20-0/perf.cn" readcn  writeall = 1
  set NP0 [MD++_Get "NP"]

  MD++ incnfile = NEBinit.cn  readcn setconfig1 #A
  MD++ incnfile = NEBfinal.cn readcn setconfig2 #B
  MD++ incnfile = NEBinit.cn  readcn saveH      #A

  MD++ freeallatoms

  setup_window

  set NP [MD++_Get "NP"]
  set vacuumratio  [expr 1.0 - ($NP*1.0)/$NP0]
  puts "vacuumratio = $vacuumratio"
  MD++ vacuumratio = $vacuumratio

  if { 1 } { # determine number of group 1/2 atoms 
  #Top surface atoms selection
  MD++ plot_color_axis = 2 eval # 2 means CSD
  MD++ input = \[ 1 -10 10 0.2 10 -10 10 10 100 \] fixatoms_by_pos_topol
  MD++ input = 1 setfixedatomsgroup 
  set N1 [MD++_Get "NPfixed"]
  MD++ freeallatoms
  #Bottom surface atoms selection
  MD++ input = \[ 1 -10 10 -10 -0.2 -10 10 10 100 \] fixatoms_by_pos_topol
  MD++ input = 2 setfixedatomsgroup 
  set N2 [MD++_Get "NPfixed"]
  MD++ freeallatoms
  } else {
  # the number of atoms in each group should equal to the following
    set N1 4800
    set N2 4800
  }

  set lx [MD++_Get "H_11"]
  set lz [MD++_Get "H_33"]

  #Compute Fext in terms of eV/A
  set f1x [expr -${sigma_xy}*$lx*$lz/$N1*0.6242E-05]
  set f1y [expr -${sigma_yy}*$lx*$lz/$N1*0.6242E-05]
  set f1z [expr -${sigma_yz}*$lx*$lz/$N1*0.6242E-05]

  set f2x [expr  ${sigma_xy}*$lx*$lz/$N1*0.6242E-05]
  set f2y [expr  ${sigma_yy}*$lx*$lz/$N1*0.6242E-05]
  set f2z [expr  ${sigma_yz}*$lx*$lz/$N1*0.6242E-05]

  MD++ extforce = \[ 2   $f1x $f1y $f1z  $f2x $f2y $f2z \]

  # do we need to broadcast extforce ?
  MD++ Broadcast_Atoms

  MD++ fixallatoms constrain_fixedatoms freeallatoms
  MD++ chainlength = [expr $ncpu-1] 

  #MD++ timestep = 0.01 printfreq = 2
  MD++ timestep = 0.07 printfreq = 2
  #MD++ timestep = 0.1 printfreq = 2

  # make left end not movable
  #MD++ nebspec = \[ 0 1 0 1 1 0 \] totalsteps = $NEBMAXSTEP equilsteps = $EQUILSTEP

  # set path length constant
  MD++ nebspec = \[ 0 1 0 1 1 2 \] totalsteps = $NEBMAXSTEP equilsteps = $EQUILSTEP
  #MD++ nebspec = \[ 0 1 0 1 1 1 \] totalsteps = $NEBMAXSTEP equilsteps = $EQUILSTEP
  ## [ relax_surround redistr_freq moveleftend moverightend yesclimbimage islengthconstant ]

  MD++ Slave_chdir #Slave change to same directory as Master (dirname)
  MD++ allocchain_parallel initRchain_parallel

 if { $rerun } {
  # rerun simulation using previously relaxed path at same stress
  if { $exist > 0 } {
    MD++ incnfile = stringrelax.chain.cn readRchain_parallel
  }
 } else {
  # relax path using a previously relaxed path at lower stress
  if { $exist > 0 } {
    MD++ incnfile = Prev_NEBinit.cn readcn  
    MD++ incnfile = Prev_stringrelax.chain.cn readRchain_parallel
    # shift the box while keeping scaled coordinates the same
    MD++ Broadcast_H copyRchaintoCN_parallel
    MD++ restoreH Broadcast_H copyCNtoRchain_parallel
    MD++ incnfile = NEBinit.cn  readcn copyCNtoRchain_masteronly
  }
 }
  MD++ eval
  MD++ stringrelax_parallel 

  MD++ finalcnfile = stringrelax.chain.cn writeRchain_parallel
  MD++ quit_all

} elseif { $status == 13.5 } {
  # Parallel NEB relaxation (for structures relaxed under stress)
  # using relaxed path at a different stress as initial path
  # Fleischer Mechanism
  # to compile
  #  make eam build=R SYS=mc2_mpich SPEC+=-DMAXCONSTRAINATOMS=800001
  # to run (in PBS file)
  #  mpiexec -np $ncpu bin1/eam_mc2_mpich scripts/work/ni_disl/ni_screw_cross_slip 13 5000 0 0 0

  set NEBMAXSTEP $n
  set EQUILSTEP  [expr $NEBMAXSTEP - 100 ]
  if { $EQUILSTEP < 0 } {
     set EQUILSTEP 0
  }

  set sigma_xy   0.0
  set sigma_yz   $n2
  set sigma_zz   $n3
  set sigma_yy   $n4
  set sigma_xz   $n5

  set n5prev   400

  set iselect 6

  set rerun 1

  set ncpu [ MD++_Get numDomains ]
  if { $ncpu < 1 } { set ncpu 1 }
  puts "ncpu = $ncpu"
  
  MD++ setoverwrite 
  MD++ dirname = runs/ni-Rao-cs-${status}_${n}_${n2}_${n3}_${n4}_${n5}

  MD++ NNM = 200
  MD++ writeall = 1

  exec ln -f -s ../ni-Rao-screw-30x30x20-5.5/stress_${n2}_${n3}_${n4}_${n5}/ni-screw-gp-relax-${n2}_${n3}_${n4}_${n5}.cn NEBinit.cn
  exec ln -f -s ../ni-Rao-screw-30x30x20-6.5/stress_${n2}_${n3}_${n4}_${n5}/ni-cross-slip-relax-${iselect}.cn NEBfinal.cn

 if { $rerun } {
  # rerun simulation using previously relaxed path at same stress
  set check [ glob -nocomplain "stringrelax.chain.cn.cpu00" ]
  set exist [ llength $check ]
  if { $exist > 0 } {
     for { set i 0 } { $i < $ncpu } { incr i 1 } {
         set ival [ format %02d $i ]
         exec cp stringrelax.chain.cn.cpu$ival Prev_stringrelax.chain.cn.cpu$ival
     }
     exec cp stringeng.out Prev_stringeng.out
  }
 } else {
  # relax path using a previously relaxed path at lower stress
  set check [ glob -nocomplain "../ni-cs-12_${n}_${n2}_${n3}_${n4}_${n5prev}/stringrelax.chain.cn.cpu00" ]
  set exist [ llength $check ]
  if { $exist > 0 } {
     for { set i 0 } { $i < $ncpu } { incr i 1 } {
         set ival [ format %02d $i ]
         exec cp ../ni-cs-12_${n}_${n2}_${n3}_${n4}_${n5prev}/stringrelax.chain.cn.cpu$ival Prev_stringrelax.chain.cn.cpu$ival
     }
     #exec cp stringeng.out Prev_stringeng.out
     exec cp ../ni-cs-12_${n}_${n2}_${n3}_${n4}_${n5prev}/NEBinit.cn Prev_NEBinit.cn
  }
 }

  readpot
  MD++ Broadcast_EAM_Param

  MD++ incnfile = "../ni-Rao-screw-30x30x20-0.5/perf.cn" readcn  writeall = 1
  set NP0 [MD++_Get "NP"]

  MD++ incnfile = NEBinit.cn  readcn setconfig1 #A
  MD++ incnfile = NEBfinal.cn readcn setconfig2 #B
  MD++ incnfile = NEBinit.cn  readcn saveH      #A

  MD++ freeallatoms

  setup_window

  set NP [MD++_Get "NP"]
  set vacuumratio  [expr 1.0 - ($NP*1.0)/$NP0]
  puts "vacuumratio = $vacuumratio"
  MD++ vacuumratio = $vacuumratio

  if { 1 } { # determine number of group 1/2 atoms 
  #Top surface atoms selection
  MD++ plot_color_axis = 2 eval # 2 means CSD
  MD++ input = \[ 1 -10 10 0.2 10 -10 10 10 100 \] fixatoms_by_pos_topol
  MD++ input = 1 setfixedatomsgroup 
  set N1 [MD++_Get "NPfixed"]
  MD++ freeallatoms
  #Bottom surface atoms selection
  MD++ input = \[ 1 -10 10 -10 -0.2 -10 10 10 100 \] fixatoms_by_pos_topol
  MD++ input = 2 setfixedatomsgroup 
  set N2 [MD++_Get "NPfixed"]
  MD++ freeallatoms
  } else {
  # the number of atoms in each group should equal to the following
    set N1 4800
    set N2 4800
  }

  set lx [MD++_Get "H_11"]
  set lz [MD++_Get "H_33"]

  #Compute Fext in terms of eV/A
  set f1x [expr -${sigma_xy}*$lx*$lz/$N1*0.6242E-05]
  set f1y [expr -${sigma_yy}*$lx*$lz/$N1*0.6242E-05]
  set f1z [expr -${sigma_yz}*$lx*$lz/$N1*0.6242E-05]

  set f2x [expr  ${sigma_xy}*$lx*$lz/$N1*0.6242E-05]
  set f2y [expr  ${sigma_yy}*$lx*$lz/$N1*0.6242E-05]
  set f2z [expr  ${sigma_yz}*$lx*$lz/$N1*0.6242E-05]

  MD++ extforce = \[ 2   $f1x $f1y $f1z  $f2x $f2y $f2z \]

  # do we need to broadcast extforce ?
  MD++ Broadcast_Atoms

  MD++ fixallatoms constrain_fixedatoms freeallatoms
  MD++ chainlength = [expr $ncpu-1] 

  #MD++ timestep = 0.01 printfreq = 2
  MD++ timestep = 0.07 printfreq = 2
  #MD++ timestep = 0.1 printfreq = 2

  # make left end not movable
  #MD++ nebspec = \[ 0 1 0 1 1 0 \] totalsteps = $NEBMAXSTEP equilsteps = $EQUILSTEP

  # set path length constant
  MD++ nebspec = \[ 0 1 0 1 1 2 \] totalsteps = $NEBMAXSTEP equilsteps = $EQUILSTEP
  #MD++ nebspec = \[ 0 1 0 1 1 1 \] totalsteps = $NEBMAXSTEP equilsteps = $EQUILSTEP
  ## [ relax_surround redistr_freq moveleftend moverightend yesclimbimage islengthconstant ]

  MD++ Slave_chdir #Slave change to same directory as Master (dirname)
  MD++ allocchain_parallel initRchain_parallel

 if { $rerun } {
  # rerun simulation using previously relaxed path at same stress
  if { $exist > 0 } {
    MD++ incnfile = stringrelax.chain.cn readRchain_parallel
  }
 } else {
  # relax path using a previously relaxed path at lower stress
  if { $exist > 0 } {
    MD++ incnfile = Prev_NEBinit.cn readcn  
    MD++ incnfile = Prev_stringrelax.chain.cn readRchain_parallel
    # shift the box while keeping scaled coordinates the same
    MD++ Broadcast_H copyRchaintoCN_parallel
    MD++ restoreH Broadcast_H copyCNtoRchain_parallel
    MD++ incnfile = NEBinit.cn  readcn copyCNtoRchain_masteronly
  }
 }
  MD++ eval
  MD++ stringrelax_parallel 

  MD++ finalcnfile = stringrelax.chain.cn writeRchain_parallel
  MD++ quit_all

} elseif { $status == 14 } {
  # Parallel NEB relaxation (for structures relaxed under stress)
  # to compile
  #  make eam build=R SYS=mc2_mpich SPEC+=-DMAXCONSTRAINATOMS=800001
  # to run (in PBS file)
  #  mpiexec -np $ncpu bin1/eam_mc2_mpich scripts/work/ni_disl/ni_screw_cross_slip 12 5000 0 0 0

  set NEBMAXSTEP $n
  set EQUILSTEP  [expr $NEBMAXSTEP - 100 ]
  if { $EQUILSTEP < 0 } {
     set EQUILSTEP 0
  }

  set sigma_xy   0.0
  set sigma_yz   $n2
  set sigma_zz   $n3
  set sigma_yy   $n4
  set sigma_xz   $n5

 if { $sigma_yy != 0 } {
  if { $sigma_yy <= 100 } {
     set iselect 13
  } elseif { $sigma_yy <= 300 } {
     set iselect 10
  } elseif { $sigma_yy <= 600 } {
     set iselect 10
  } else {
     set iselect 6
  }
 }
  # Escaig stress
 if { $sigma_yz != 0 } {
  if { $sigma_yz <= 300 } {
     set iselect 20
  } elseif { $sigma_yz <= 500 } {
     set iselect 10
  } elseif { $sigma_yz <= 600 } {
     set iselect 8
  } else {
     set iselect 5
  }
 }
  # Schmid stress
 if { $sigma_xz != 0 } {
  if { $sigma_xz <= 100 } {
     set iselect 13
  } elseif { $sigma_xz <= 300 } {
     set iselect 11
  } elseif { $sigma_xz <= 600 } {
     set iselect 10
  } else {
     set iselect 5
  }
 }

  set ncpu [ MD++_Get numDomains ]
  if { $ncpu < 1 } { set ncpu 1 }
  puts "ncpu = $ncpu"
  
  MD++ setoverwrite 
  MD++ dirname = runs/ni-Rao-cs-${status}_${n}_${n2}_${n3}_${n4}_${n5}

  MD++ NNM = 200
  MD++ writeall = 1

  exec ln -f -s ../ni-Rao-screw-30x30x20-5/stress_${n2}_${n3}_${n4}_${n5}/ni-screw-gp-relax-${n2}_${n3}_${n4}_${n5}.cn NEBinit.cn
  exec ln -f -s ../ni-Rao-screw-30x30x20-6/stress_${n2}_${n3}_${n4}_${n5}/ni-cross-slip-relax-${iselect}.cn NEBfinal.cn

  # use new A and B'
  #exec ln -f -s chain_no_7.cn NEBinit.cn
  #exec ln -f -s chain_no_23.cn NEBfinal.cn


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

  readpot
  MD++ Broadcast_EAM_Param

  MD++ incnfile = "../ni-Rao-screw-30x30x20-0/perf.cn" readcn  writeall = 1
  set NP0 [MD++_Get "NP"]

  # a test for new state A (0_0_0_400)
  #set chain_no 6
  #MD++ incnfile = NEBinit.cn  readcn 
  #MD++ incnfile = "../ni-cs-12_${n}_${n2}_${n3}_${n4}_${n5}.sav/stringrelax.chain.cn"
  #MD++ input = $chain_no  readRchain_parallel_toCN  RHtoS  setconfig1 #A

  #MD++ incnfile = NEBfinal.cn readcn setconfig2                       #B
  #MD++ incnfile = "../ni-cs-12_${n}_${n2}_${n3}_${n4}_${n5}.sav/stringrelax.chain.cn"
  #MD++ input = $chain_no  readRchain_parallel_toCN  RHtoS             #A

  MD++ incnfile = NEBinit.cn  readcn setconfig1 #A
  MD++ incnfile = NEBfinal.cn readcn setconfig2 #B
  MD++ incnfile = NEBinit.cn  readcn            #A

  MD++ freeallatoms

  setup_window

  set NP [MD++_Get "NP"]
  set vacuumratio  [expr 1.0 - ($NP*1.0)/$NP0]
  puts "vacuumratio = $vacuumratio"
  MD++ vacuumratio = $vacuumratio

  if { 1 } { # determine number of group 1/2 atoms 
  #Top surface atoms selection
  MD++ plot_color_axis = 2 eval # 2 means CSD
  MD++ input = \[ 1 -10 10 0.2 10 -10 10 10 100 \] fixatoms_by_pos_topol
  MD++ input = 1 setfixedatomsgroup 
  set N1 [MD++_Get "NPfixed"]
  MD++ freeallatoms
  #Bottom surface atoms selection
  MD++ input = \[ 1 -10 10 -10 -0.2 -10 10 10 100 \] fixatoms_by_pos_topol
  MD++ input = 2 setfixedatomsgroup 
  set N2 [MD++_Get "NPfixed"]
  MD++ freeallatoms
  } else {
  # the number of atoms in each group should equal to the following
    set N1 4800
    set N2 4800
  }

  set lx [MD++_Get "H_11"]
  set lz [MD++_Get "H_33"]

  #Compute Fext in terms of eV/A
  set f1x [expr -${sigma_xy}*$lx*$lz/$N1*0.6242E-05]
  set f1y [expr -${sigma_yy}*$lx*$lz/$N1*0.6242E-05]
  set f1z [expr -${sigma_yz}*$lx*$lz/$N1*0.6242E-05]

  set f2x [expr  ${sigma_xy}*$lx*$lz/$N1*0.6242E-05]
  set f2y [expr  ${sigma_yy}*$lx*$lz/$N1*0.6242E-05]
  set f2z [expr  ${sigma_yz}*$lx*$lz/$N1*0.6242E-05]

  MD++ extforce = \[ 2   $f1x $f1y $f1z  $f2x $f2y $f2z \]

  #MD++ freeallatoms clearR0 refreshnnlist eval  quit_all

  # do we need to broadcast extforce ?
  MD++ Broadcast_Atoms

  MD++ fixallatoms constrain_fixedatoms freeallatoms
  MD++ chainlength = [expr $ncpu-1] 

  MD++ timestep = 0.07 printfreq = 2
  #MD++ timestep = 0.1 printfreq = 2

  # make left end not movable
  #MD++ nebspec = \[ 0 1 0 1 1 0 \] totalsteps = $NEBMAXSTEP equilsteps = $EQUILSTEP

  set NEBSTEP1 1000

  # set path length constant
  #MD++ nebspec = \[ 0 1 0 1 1 2 \] totalsteps = $NEBMAXSTEP equilsteps = $EQUILSTEP
  MD++ nebspec = \[ 0 1 0 1 1 1 \] totalsteps = $NEBMAXSTEP equilsteps = $EQUILSTEP
  #MD++ nebspec = \[ 0 1 0 1 1 1 \] totalsteps = $NEBSTEP1 equilsteps = $EQUILSTEP
  ## [ relax_surround redistr_freq moveleftend moverightend yesclimbimage islengthconstant ]

  MD++ Slave_chdir #Slave change to same directory as Master (dirname)
  MD++ allocchain_parallel initRchain_parallel

  if { $exist > 0 } {
    MD++ incnfile = stringrelax.chain.cn readRchain_parallel
  }

  MD++ eval
  MD++ reparametrization_with_trimming
  #MD++ stringrelax_parallel_2
  MD++ finalcnfile = stringrelax_trim.chain.cn writeRchain_parallel
  set check [ glob -nocomplain "stringeng_step*" ]
  set nfile [ llength $check]

  MD++ quit_all


} elseif { $status == 50 } {
  # visualize path from previous Parallel NEB relaxation 
  MD++ setnolog
  initmd "view"
  readpot

  set sigma_xy   0.0
  set sigma_yz   $n2
  set sigma_zz   $n3
  set sigma_yy   $n4
  set sigma_xz   $n5

  set total_no 23
  MD++ chainlength = $total_no

  #MD++ incnfile = "../ni-Rao-cs-12_${n}_${n2}_${n3}_${n4}_${n5}.sav/NEBinit.cn" readcn saveH
  #MD++ incnfile = "../ni-Rao-cs-12_${n}_${n2}_${n3}_${n4}_${n5}.sav/stringrelax.chain.cn"
  #MD++ incnfile = "../ni-Rao-cs-12_${n}_${n2}_${n3}_${n4}_${n5}/NEBinit.cn" readcn saveH
  #MD++ incnfile = "../ni-Rao-cs-12_${n}_${n2}_${n3}_${n4}_${n5}/stringrelax.chain.cn"
  MD++ incnfile = "../ni-Rao-cs-112_${n}_${n2}_${n3}_${n4}_${n5}/NEBinit.cn" readcn saveH
  MD++ incnfile = "../ni-Rao-cs-112_${n}_${n2}_${n3}_${n4}_${n5}/stringrelax.chain.cn"

  set N1 4800
  set N2 4800
  set lx [MD++_Get "H_11"]
  set lz [MD++_Get "H_33"]

  #Compute Fext in terms of eV/A
  set f1x [expr -${sigma_xy}*$lx*$lz/$N1*0.6242E-05]
  set f1y [expr -${sigma_yy}*$lx*$lz/$N1*0.6242E-05]
  set f1z [expr -${sigma_yz}*$lx*$lz/$N1*0.6242E-05]

  set f2x [expr  ${sigma_xy}*$lx*$lz/$N1*0.6242E-05]
  set f2y [expr  ${sigma_yy}*$lx*$lz/$N1*0.6242E-05]
  set f2z [expr  ${sigma_yz}*$lx*$lz/$N1*0.6242E-05]

  MD++ extforce = \[ 2   $f1x $f1y $f1z  $f2x $f2y $f2z \]

  # write one file into cn format
  if { 0 } {
    #set chain_no 0
    #set chain_no 6
    set chain_no 23
    #set chain_no 23
    MD++ input = $chain_no  readRchain_parallel_toCN  RHtoS
    # must use writeall = 0 to reuse cn file as initial condition, not sure why
    MD++ finalcnfile = "chain_no_${chain_no}.cn"  writeall = 0 writecn
    MD++ eval
    exitmd 
  }

  setup_window

  MD++ zipfiles = 0
  for { set iter 0 } { $iter <= $total_no } { incr iter 1 } {
     # set chain number
     set chain_no $iter
     MD++ input = $chain_no  readRchain_parallel_toCN  RHtoS
     MD++ NCS = 12 eval calcentralsymmetry 

     if { $iter == 0 } {
         openwindow
     }
     MD++ eval plot 
     # write cn or cfg file
     #MD++ finalcnfile = "chain_no_${chain_no}.cn"  writecn
     MD++ finalcnfile = "chain_no_${chain_no}.cfg" writeatomeyecfg
     MD++ sleepseconds = 1  sleep

  }
  MD++ sleepseconds = 10  sleep

  exitmd
} elseif { $status == 50.1 } {
  # visualize path from previous Parallel NEB relaxation
  # Keep rename .cfg files to maintain the stress state in the file name
  MD++ setnolog
  initmd "view"
  readpot

  set sigma_xy   0.0
  set sigma_yz   $n2
  set sigma_zz   $n3
  set sigma_yy   $n4
  set sigma_xz   $n5

  set total_no 23
  MD++ chainlength = $total_no

  #MD++ incnfile = "../ni-Rao-cs-12_${n}_${n2}_${n3}_${n4}_${n5}.sav/NEBinit.cn" readcn saveH
  #MD++ incnfile = "../ni-Rao-cs-12_${n}_${n2}_${n3}_${n4}_${n5}.sav/stringrelax.chain.cn"
  #MD++ incnfile = "../ni-Rao-cs-12_${n}_${n2}_${n3}_${n4}_${n5}/NEBinit.cn" readcn saveH
  #MD++ incnfile = "../ni-Rao-cs-12_${n}_${n2}_${n3}_${n4}_${n5}/stringrelax.chain.cn"
  MD++ incnfile = "../ni-Rao-cs-112_${n}_${n2}_${n3}_${n4}_${n5}/NEBinit.cn" readcn saveH
  MD++ incnfile = "../ni-Rao-cs-112_${n}_${n2}_${n3}_${n4}_${n5}/stringrelax.chain.cn"

  set N1 4800
  set N2 4800
  set lx [MD++_Get "H_11"]
  set lz [MD++_Get "H_33"]

  #Compute Fext in terms of eV/A
  set f1x [expr -${sigma_xy}*$lx*$lz/$N1*0.6242E-05]
  set f1y [expr -${sigma_yy}*$lx*$lz/$N1*0.6242E-05]
  set f1z [expr -${sigma_yz}*$lx*$lz/$N1*0.6242E-05]

  set f2x [expr  ${sigma_xy}*$lx*$lz/$N1*0.6242E-05]
  set f2y [expr  ${sigma_yy}*$lx*$lz/$N1*0.6242E-05]
  set f2z [expr  ${sigma_yz}*$lx*$lz/$N1*0.6242E-05]

  MD++ extforce = \[ 2   $f1x $f1y $f1z  $f2x $f2y $f2z \]

  # write one file into cn format
  if { 0 } {
    #set chain_no 0
    #set chain_no 6
    set chain_no 23
    #set chain_no 23
    MD++ input = $chain_no  readRchain_parallel_toCN  RHtoS
    # must use writeall = 0 to reuse cn file as initial condition, not sure why
    MD++ finalcnfile = "chain_no_${chain_no}.cn"  writeall = 0 writecn
    MD++ eval
    exitmd 
  }

  setup_window

  MD++ zipfiles = 0
  for { set iter 0 } { $iter <= $total_no } { incr iter 1 } {
     # set chain number
     set chain_no $iter
     MD++ input = $chain_no  readRchain_parallel_toCN  RHtoS
     MD++ NCS = 12 eval calcentralsymmetry 

     if { $iter == 0 } {
         openwindow
     }
     MD++ eval plot 
     # write cn or cfg file
     #MD++ finalcnfile = "chain_no_${chain_no}.cn"  writecn
     #MD++ finalcnfile = "chain_no_${chain_no}_112_${n2}_${n3}_${n4}_${n5}.cfg" writeatomeyecfg
     MD++ finalcnfile = "status_112_${n2}_${n3}_${n4}_${n5}_chain_no_${chain_no}.cfg" writeatomeyecfg
     MD++ sleepseconds = 1  sleep

  }
  MD++ sleepseconds = 10  sleep

  exitmd
} elseif { $status == 50.5 } {
  # visualize path from previous Parallel NEB relaxation
  # Keep rename .cfg files to maintain the stress state in the file name
  MD++ setnolog
  initmd "view"
  readpot

  set sigma_xy   0.0
  set sigma_yz   $n2
  set sigma_zz   $n3
  set sigma_yy   $n4
  set sigma_xz   $n5

  set total_no 23
  MD++ chainlength = $total_no

  #MD++ incnfile = "../ni-Rao-cs-12_${n}_${n2}_${n3}_${n4}_${n5}.sav/NEBinit.cn" readcn saveH
  #MD++ incnfile = "../ni-Rao-cs-12_${n}_${n2}_${n3}_${n4}_${n5}.sav/stringrelax.chain.cn"
  #MD++ incnfile = "../ni-Rao-cs-12_${n}_${n2}_${n3}_${n4}_${n5}/NEBinit.cn" readcn saveH
  #MD++ incnfile = "../ni-Rao-cs-12_${n}_${n2}_${n3}_${n4}_${n5}/stringrelax.chain.cn"
  MD++ incnfile = "../ni-Rao-cs-112.5_${n}_${n2}_${n3}_${n4}_${n5}/NEBinit.cn" readcn saveH
  MD++ incnfile = "../ni-Rao-cs-112.5_${n}_${n2}_${n3}_${n4}_${n5}/stringrelax.chain.cn"

  set N1 4800
  set N2 4800
  set lx [MD++_Get "H_11"]
  set lz [MD++_Get "H_33"]

  #Compute Fext in terms of eV/A
  set f1x [expr -${sigma_xy}*$lx*$lz/$N1*0.6242E-05]
  set f1y [expr -${sigma_yy}*$lx*$lz/$N1*0.6242E-05]
  set f1z [expr -${sigma_yz}*$lx*$lz/$N1*0.6242E-05]

  set f2x [expr  ${sigma_xy}*$lx*$lz/$N1*0.6242E-05]
  set f2y [expr  ${sigma_yy}*$lx*$lz/$N1*0.6242E-05]
  set f2z [expr  ${sigma_yz}*$lx*$lz/$N1*0.6242E-05]

  MD++ extforce = \[ 2   $f1x $f1y $f1z  $f2x $f2y $f2z \]

  # write one file into cn format
  if { 0 } {
    #set chain_no 0
    #set chain_no 6
    set chain_no 23
    #set chain_no 23
    MD++ input = $chain_no  readRchain_parallel_toCN  RHtoS
    # must use writeall = 0 to reuse cn file as initial condition, not sure why
    MD++ finalcnfile = "chain_no_${chain_no}.cn"  writeall = 0 writecn
    MD++ eval
    exitmd 
  }

  setup_window

  MD++ zipfiles = 0
  for { set iter 0 } { $iter <= $total_no } { incr iter 1 } {
     # set chain number
     set chain_no $iter
     MD++ input = $chain_no  readRchain_parallel_toCN  RHtoS
     MD++ NCS = 12 eval calcentralsymmetry 

     if { $iter == 0 } {
         openwindow
     }
     MD++ eval plot 
     # write cn or cfg file
     #MD++ finalcnfile = "chain_no_${chain_no}.cn"  writecn
     MD++ finalcnfile = "status_112.5_${n2}_${n3}_${n4}_${n5}_chain_no_${chain_no}.cfg" writeatomeyecfg
     MD++ sleepseconds = 1  sleep

  }
  MD++ sleepseconds = 10  sleep

  exitmd
} elseif { $status == 51 } {
  # relax a given state in relaxed path from previous Parallel NEB relaxation 
  # Generally this is state A (0), but can be set to anything using set chain_no
  MD++ setnolog
  initmd "view"
  readpot

  set sigma_xy   0.0
  set sigma_yz   $n2
  set sigma_zz   $n3
  set sigma_yy   $n4
  set sigma_xz   $n5

  set total_no 23
  MD++ chainlength = $total_no

  MD++ incnfile = "../ni-Rao-cs-112_${n}_${n2}_${n3}_${n4}_${n5}/NEBinit.cn" readcn saveH
  MD++ incnfile = "../ni-Rao-cs-112_${n}_${n2}_${n3}_${n4}_${n5}/stringrelax.chain.cn"

  set N1 4800
  set N2 4800
  set lx [MD++_Get "H_11"]
  set lz [MD++_Get "H_33"]

  #Compute Fext in terms of eV/A
  set f1x [expr -${sigma_xy}*$lx*$lz/$N1*0.6242E-05]
  set f1y [expr -${sigma_yy}*$lx*$lz/$N1*0.6242E-05]
  set f1z [expr -${sigma_yz}*$lx*$lz/$N1*0.6242E-05]

  set f2x [expr  ${sigma_xy}*$lx*$lz/$N1*0.6242E-05]
  set f2y [expr  ${sigma_yy}*$lx*$lz/$N1*0.6242E-05]
  set f2z [expr  ${sigma_yz}*$lx*$lz/$N1*0.6242E-05]

  MD++ extforce = \[ 2   $f1x $f1y $f1z  $f2x $f2y $f2z \]

  set chain_no 0
  #set chain_no 6
  #set chain_no 7
  #set chain_no 10
  #set chain_no 23

  setup_window
  #openwindow

  MD++ input = $chain_no  readRchain_parallel_toCN  RHtoS eval #quit

##
  set EPOT_A1 [MD++_Get "EPOT"]
##
  relax_fixbox

  MD++ finalcnfile = "chain_no_${chain_no}_relaxed.cn"  writeall = 0 writecn
  MD++ NCS = 12 eval calcentralsymmetry 

##
  set EPOT_A2 [MD++_Get "EPOT"]
  set En_Diff [expr $EPOT_A2 - $EPOT_A1 ]
  puts "Relaxation Difference = $En_Diff"

  set fid [open "../ni-Rao-cs-112_${n}_${n2}_${n3}_${n4}_${n5}/EArelax.out" "w"]
  puts $fid " $En_Diff "
  close $fid

##

  MD++ finalcnfile = "chain_no_${chain_no}_relaxed.cfg" writeatomeyecfg
#  MD++ sleepseconds = 100  sleep

  exitmd
  
} elseif { $status == 51.5 } {
  # relax a given state in relaxed path from previous Parallel NEB relaxation 
  # Generally this is state A (0), but can be set to anything using set chain_no
  # Set up to do relaxation for Fleischer runs
  MD++ setnolog
  initmd "view"
  readpot

  set sigma_xy   0.0
  set sigma_yz   $n2
  set sigma_zz   $n3
  set sigma_yy   $n4
  set sigma_xz   $n5

  set total_no 23
  MD++ chainlength = $total_no

  MD++ incnfile = "../ni-Rao-cs-112.5_${n}_${n2}_${n3}_${n4}_${n5}/NEBinit.cn" readcn saveH
  MD++ incnfile = "../ni-Rao-cs-112.5_${n}_${n2}_${n3}_${n4}_${n5}/stringrelax.chain.cn"

  set N1 4800
  set N2 4800
  set lx [MD++_Get "H_11"]
  set lz [MD++_Get "H_33"]

  #Compute Fext in terms of eV/A
  set f1x [expr -${sigma_xy}*$lx*$lz/$N1*0.6242E-05]
  set f1y [expr -${sigma_yy}*$lx*$lz/$N1*0.6242E-05]
  set f1z [expr -${sigma_yz}*$lx*$lz/$N1*0.6242E-05]

  set f2x [expr  ${sigma_xy}*$lx*$lz/$N1*0.6242E-05]
  set f2y [expr  ${sigma_yy}*$lx*$lz/$N1*0.6242E-05]
  set f2z [expr  ${sigma_yz}*$lx*$lz/$N1*0.6242E-05]

  MD++ extforce = \[ 2   $f1x $f1y $f1z  $f2x $f2y $f2z \]

  set chain_no 0
  #set chain_no 6
  #set chain_no 7
  #set chain_no 10
  #set chain_no 23

  setup_window
  #openwindow

  MD++ input = $chain_no  readRchain_parallel_toCN  RHtoS eval #quit

##
  set EPOT_A1 [MD++_Get "EPOT"]
##
  relax_fixbox

  MD++ finalcnfile = "chain_no_${chain_no}_relaxed.cn"  writeall = 0 writecn
  MD++ NCS = 12 eval calcentralsymmetry 

##
  set EPOT_A2 [MD++_Get "EPOT"]
  set En_Diff [expr $EPOT_A2 - $EPOT_A1 ]
  puts "Relaxation Difference = $En_Diff"

  set fid [open "../ni-Rao-cs-112.5_${n}_${n2}_${n3}_${n4}_${n5}/EArelax.out" "w"]
  puts $fid " $En_Diff "
  close $fid

##

  MD++ finalcnfile = "chain_no_${chain_no}_relaxed.cfg" writeatomeyecfg
#  MD++ sleepseconds = 100  sleep

  exitmd

} elseif { $status == 52 } {
  # visualize state A relaxed under different stress
  MD++ setnolog
  initmd "view"
  readpot

  set sigma_xy   0.0
  set sigma_yz   $n2
  set sigma_zz   $n3
  set sigma_yy   $n4
  set sigma_xz   $n5

  set total_no 23
  MD++ chainlength = $total_no

  MD++ incnfile = "../ni-Rao-screw-30x30x20-5/stress_${n2}_${n3}_${n4}_${n5}/ni-screw-gp-relax-${n2}_${n3}_${n4}_${n5}.cn" readcn
  setup_window
  openwindow

  MD++ zipfiles = 0
  MD++ finalcnfile = "ni-screw-gp-relax-${n2}_${n3}_${n4}_${n5}.cfg" writeatomeyecfg

  # write CSD data into file
  set fid [open "csd.out" "w"]
  set NP [MD++_Get "NP"]

  for { set i 0 } { $i < $NP } { incr i } {
      set x [MD++_Get R([expr 3*$i])]
      set y [MD++_Get R([expr 3*$i+1])]
      set z [MD++_Get R([expr 3*$i+2])]
      set topol [MD++_Get TOPOL($i)]
      if { ($y > -1) && ($y < 3) } {
         puts $fid " $z   $topol "
      }
  }
  close $fid

  MD++ sleepseconds = 100  sleep

  exitmd

} elseif { $status == 111 } {
  # Parallel NEB relaxation (for structures relaxed under stress)
  # to compile
  #  make eam build=R SYS=mc2_mpich SPEC+=-DMAXCONSTRAINATOMS=800001
  # to run (in PBS file)
  #  mpiexec -np $ncpu bin1/eam_mc2_mpich scripts/work/ni_disl/ni_screw_cross_slip 12 5000 0 0 0

  set NEBMAXSTEP $n
  set EQUILSTEP  [expr $NEBMAXSTEP - 100 ]
  if { $EQUILSTEP < 0 } {
     set EQUILSTEP 0
  }

  set sigma_xy   0.0
  set sigma_yz   $n2
  set sigma_zz   $n3
  set sigma_yy   $n4
  set sigma_xz   $n5

 if { $sigma_yy != 0 } {
  if { $sigma_yy <= 100 } {
     set iselect 13
  } elseif { $sigma_yy <= 300 } {
     set iselect 10
  } elseif { $sigma_yy <= 600 } {
     set iselect 10
  } else {
     set iselect 6
  }
 }
  # Escaig stress
 if { $sigma_yz != 0 } {
  if { $sigma_yz <= 300 } {
     set iselect 20
  } elseif { $sigma_yz <= 500 } {
     set iselect 10
  } elseif { $sigma_yz <= 600 } {
     set iselect 8
  } else {
     set iselect 5
  }
 }
  # Schmid stress
 if { $sigma_xz != 0 } {
  if { $sigma_xz <= 100 } {
     set iselect 13
  } elseif { $sigma_xz <= 300 } {
     set iselect 11
  } elseif { $sigma_xz <= 600 } {
     set iselect 10
  } else {
     set iselect 5
  }
 }

  set ncpu [ MD++_Get numDomains ]
  if { $ncpu < 1 } { set ncpu 1 }
  puts "ncpu = $ncpu"
  
  MD++ setoverwrite 
  MD++ dirname = runs/ni-Rao-cs-${status}_${n}_${n2}_${n3}_${n4}_${n5}

  MD++ NNM = 200
  MD++ writeall = 1

  exec ln -f -s ../ni-Rao-screw-30x30x20-5/stress_${n2}_${n3}_${n4}_${n5}/ni-screw-gp-relax-${n2}_${n3}_${n4}_${n5}.cn NEBinit.cn
  exec ln -f -s ../ni-Rao-screw-30x30x20-6/stress_${n2}_${n3}_${n4}_${n5}/ni-cross-slip-relax-${iselect}.cn NEBfinal.cn

  # use new A and B'
  #exec ln -f -s chain_no_7.cn NEBinit.cn
  #exec ln -f -s chain_no_23.cn NEBfinal.cn


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

  readpot
  MD++ Broadcast_EAM_Param

  MD++ incnfile = "../ni-Rao-screw-30x30x20-0/perf.cn" readcn  writeall = 1
  set NP0 [MD++_Get "NP"]

  # a test for new state A (0_0_0_400)
  #set chain_no 6
  #MD++ incnfile = NEBinit.cn  readcn 
  #MD++ incnfile = "../ni-cs-12_${n}_${n2}_${n3}_${n4}_${n5}.sav/stringrelax.chain.cn"
  #MD++ input = $chain_no  readRchain_parallel_toCN  RHtoS  setconfig1 #A

  #MD++ incnfile = NEBfinal.cn readcn setconfig2                       #B
  #MD++ incnfile = "../ni-cs-12_${n}_${n2}_${n3}_${n4}_${n5}.sav/stringrelax.chain.cn"
  #MD++ input = $chain_no  readRchain_parallel_toCN  RHtoS             #A

  MD++ incnfile = NEBinit.cn  readcn setconfig1 #A
  MD++ incnfile = NEBfinal.cn readcn setconfig2 #B
  MD++ incnfile = NEBinit.cn  readcn            #A

  MD++ freeallatoms

  setup_window

  set NP [MD++_Get "NP"]
  set vacuumratio  [expr 1.0 - ($NP*1.0)/$NP0]
  puts "vacuumratio = $vacuumratio"
  MD++ vacuumratio = $vacuumratio

  if { 1 } { # determine number of group 1/2 atoms 
  #Top surface atoms selection
  MD++ plot_color_axis = 2 eval # 2 means CSD
  MD++ input = \[ 1 -10 10 0.2 10 -10 10 10 100 \] fixatoms_by_pos_topol
  MD++ input = 1 setfixedatomsgroup 
  set N1 [MD++_Get "NPfixed"]
  MD++ freeallatoms
  #Bottom surface atoms selection
  MD++ input = \[ 1 -10 10 -10 -0.2 -10 10 10 100 \] fixatoms_by_pos_topol
  MD++ input = 2 setfixedatomsgroup 
  set N2 [MD++_Get "NPfixed"]
  MD++ freeallatoms
  } else {
  # the number of atoms in each group should equal to the following
    set N1 4800
    set N2 4800
  }

  set lx [MD++_Get "H_11"]
  set lz [MD++_Get "H_33"]

  #Compute Fext in terms of eV/A
  set f1x [expr -${sigma_xy}*$lx*$lz/$N1*0.6242E-05]
  set f1y [expr -${sigma_yy}*$lx*$lz/$N1*0.6242E-05]
  set f1z [expr -${sigma_yz}*$lx*$lz/$N1*0.6242E-05]

  set f2x [expr  ${sigma_xy}*$lx*$lz/$N1*0.6242E-05]
  set f2y [expr  ${sigma_yy}*$lx*$lz/$N1*0.6242E-05]
  set f2z [expr  ${sigma_yz}*$lx*$lz/$N1*0.6242E-05]

  MD++ extforce = \[ 2   $f1x $f1y $f1z  $f2x $f2y $f2z \]

  #MD++ freeallatoms clearR0 refreshnnlist eval  quit_all

  # do we need to broadcast extforce ?
  MD++ Broadcast_Atoms

  MD++ fixallatoms constrain_fixedatoms freeallatoms
  MD++ chainlength = [expr $ncpu-1] 

  MD++ timestep = 0.07 printfreq = 2
  #MD++ timestep = 0.1 printfreq = 2

  # make left end not movable
  #MD++ nebspec = \[ 0 1 0 1 1 0 \] totalsteps = $NEBMAXSTEP equilsteps = $EQUILSTEP

  set NEBSTEP1 1000
  set NEBSTEP1 0

  # set path length constant
  #MD++ nebspec = \[ 0 1 0 1 1 2 \] totalsteps = $NEBMAXSTEP equilsteps = $EQUILSTEP
  MD++ nebspec = \[ 0 1 0 1 1 1 \] totalsteps = $NEBMAXSTEP equilsteps = $EQUILSTEP
  #MD++ nebspec = \[ 0 1 0 1 1 1 \] totalsteps = $NEBSTEP1 equilsteps = $EQUILSTEP
  ## [ relax_surround redistr_freq moveleftend moverightend yesclimbimage islengthconstant ]

  MD++ Slave_chdir #Slave change to same directory as Master (dirname)
  MD++ allocchain_parallel initRchain_parallel

  if { $exist > 0 } {
    MD++ incnfile = stringrelax.chain.cn readRchain_parallel
  }

  MD++ eval
  MD++ stringrelax_parallel_1
  MD++ finalcnfile = stringrelax.chain.cn writeRchain_parallel
  set check [ glob -nocomplain "stringeng_step*" ]
  set nfile [ llength $check]
  exec cp stringeng.out stringeng_step[expr $nfile+1].out

  #MD++ nebspec = \[ 0 1 0 1 1 0 \] totalsteps = $NEBMAXSTEP equilsteps = $EQUILSTEP
  MD++ nebspec = \[ 0 1 0 1 1 1 \] totalsteps = $NEBSTEP1 equilsteps = $EQUILSTEP
  MD++ Broadcast_nebsetting
  MD++ stringrelax_parallel
  #MD++ stringrelax_parallel_1
  MD++ finalcnfile = stringrelax.chain.cn writeRchain_parallel
  exec cp stringeng.out stringeng_step[expr $nfile+2].out

  MD++ quit_all

} elseif { $status == 112 } {
  # Parallel NEB relaxation (for structures relaxed under stress)
  # to compile
  #  make eam build=R SYS=mc2_mpich SPEC+=-DMAXCONSTRAINATOMS=800001
  # to run (in PBS file)
  #  mpiexec -np $ncpu bin1/eam_mc2_mpich scripts/work/ni_disl/ni_screw_cross_slip 112 5000 0 0 0

  set NEBMAXSTEP $n
  #set NEBMAXSTEP 270 
  set EQUILSTEP  [expr $NEBMAXSTEP - 100 ]
  if { $EQUILSTEP < 0 } {
     set EQUILSTEP 0
  }

  set sigma_xy   0.0
  set sigma_yz   $n2
  set sigma_zz   $n3
  set sigma_yy   $n4
  set sigma_xz   $n5
  
  # Chooses the relaxed configuration to use as an ending point
  set iselect 20

 if { $sigma_yy != 0 } {
  if { $sigma_yy <= 100 } {
     set iselect 13
  } elseif { $sigma_yy <= 300 } {
     set iselect 10
  } elseif { $sigma_yy <= 600 } {
     set iselect 10
  } else {
     set iselect 10
     #set iselect 6
  }
 }
  # Escaig stress
 if { $sigma_yz != 0 } {
  if { $sigma_yz <= 300 } {
     if { $iselect > 20 } {
       set iselect 20
     }
  } elseif { $sigma_yz <= 500 } {
     if { $iselect > 10 } {
       set iselect 10
     }
  } elseif { $sigma_yz <= 600 } {
     if { $iselect > 8 } {
       set iselect 8
     }
  } else {
     if { $iselect > 40 } {
       set iselect 8
       #set iselect 5
     }
  }
 }
  # Schmid stress
 if { $sigma_xz != 0 } {
  if { $sigma_xz <= 100 } {
     if { $iselect > 40 } {
       set iselect 13
     }
  } elseif { $sigma_xz <= 300 } {
     if { $iselect > 40 } {
       set iselect 11
     }
  } elseif { $sigma_xz <= 600 } {
     if { $iselect > 40 } {
       set iselect 10
     }
  } else {
     if { $iselect > 40 } {
       #set iselect 5
       set iselect 10
     }
  }
 }

  set iselect 5
  #### For 4 20 15 Test
  #set iselect 3 
  puts "iselect = $iselect"

  set ncpu [ MD++_Get numDomains ]
  if { $ncpu < 1 } { set ncpu 1 }
  puts "ncpu = $ncpu"
  
  MD++ setoverwrite 
  MD++ dirname = runs/ni-Rao-cs-${status}_${n}_${n2}_${n3}_${n4}_${n5}

  MD++ NNM = 200
  MD++ writeall = 1

  exec ln -f -s ../ni-Rao-screw-30x30x20-5/stress_${n2}_${n3}_${n4}_${n5}/ni-screw-gp-relax-${n2}_${n3}_${n4}_${n5}.cn NEBinit.cn
  exec ln -f -s ../ni-Rao-screw-30x30x20-6/stress_${n2}_${n3}_${n4}_${n5}/ni-cross-slip-relax-${iselect}.cn NEBfinal.cn

  # use new A and B'
  #exec ln -f -s chain_no_7.cn NEBinit.cn
  #exec ln -f -s chain_no_23.cn NEBfinal.cn


  # rename previous result files in order to not overwrite them
  set check [ glob -nocomplain "stringrelax.chain.cn.cpu00" ]
  set exist [ llength $check ]
  if { $exist > 0 } {
     for { set i 0 } { $i < $ncpu } { incr i 1 } {
         set ival [ format %02d $i ]
         exec cp stringrelax.chain.cn.cpu$ival Prev_stringrelax.chain.cn.cpu$ival
     }
     if { [ llength [ glob -nocomplain "stringeng.out" ] ] > 0 } { 
       exec cp stringeng.out Prev_stringeng.out
     }
  }

  readpot
  MD++ Broadcast_EAM_Param

  MD++ incnfile = "../ni-Rao-screw-30x30x20-0/perf.cn" readcn  writeall = 1
  set NP0 [MD++_Get "NP"]

  # a test for new state A (0_0_0_400)
  #set chain_no 6
  #MD++ incnfile = NEBinit.cn  readcn 
  #MD++ incnfile = "../ni-cs-12_${n}_${n2}_${n3}_${n4}_${n5}.sav/stringrelax.chain.cn"
  #MD++ input = $chain_no  readRchain_parallel_toCN  RHtoS  setconfig1 #A

  #MD++ incnfile = NEBfinal.cn readcn setconfig2                       #B
  #MD++ incnfile = "../ni-cs-12_${n}_${n2}_${n3}_${n4}_${n5}.sav/stringrelax.chain.cn"
  #MD++ input = $chain_no  readRchain_parallel_toCN  RHtoS             #A

  MD++ incnfile = NEBinit.cn  readcn setconfig1 #A
  MD++ incnfile = NEBfinal.cn readcn setconfig2 #B
  MD++ incnfile = NEBinit.cn  readcn            #A

  MD++ freeallatoms

  setup_window

  set NP [MD++_Get "NP"]
  set vacuumratio  [expr 1.0 - ($NP*1.0)/$NP0]
  puts "vacuumratio = $vacuumratio"
  MD++ vacuumratio = $vacuumratio

  if { 1 } { # determine number of group 1/2 atoms 
  #Top surface atoms selection
  MD++ plot_color_axis = 2 eval # 2 means CSD
  MD++ input = \[ 1 -10 10 0.2 10 -10 10 10 100 \] fixatoms_by_pos_topol
  MD++ input = 1 setfixedatomsgroup 
  set N1 [MD++_Get "NPfixed"]
  MD++ freeallatoms
  #Bottom surface atoms selection
  MD++ input = \[ 1 -10 10 -10 -0.2 -10 10 10 100 \] fixatoms_by_pos_topol
  MD++ input = 2 setfixedatomsgroup 
  set N2 [MD++_Get "NPfixed"]
  MD++ freeallatoms
  } else {
  # the number of atoms in each group should equal to the following
    set N1 4800
    set N2 4800
  }

  set lx [MD++_Get "H_11"]
  set lz [MD++_Get "H_33"]

  #Compute Fext in terms of eV/A
  set f1x [expr -${sigma_xy}*$lx*$lz/$N1*0.6242E-05]
  set f1y [expr -${sigma_yy}*$lx*$lz/$N1*0.6242E-05]
  set f1z [expr -${sigma_yz}*$lx*$lz/$N1*0.6242E-05]

  set f2x [expr  ${sigma_xy}*$lx*$lz/$N1*0.6242E-05]
  set f2y [expr  ${sigma_yy}*$lx*$lz/$N1*0.6242E-05]
  set f2z [expr  ${sigma_yz}*$lx*$lz/$N1*0.6242E-05]

  MD++ extforce = \[ 2   $f1x $f1y $f1z  $f2x $f2y $f2z \]

  #MD++ freeallatoms clearR0 refreshnnlist eval  quit_all

  # do we need to broadcast extforce ?
  MD++ Broadcast_Atoms

  MD++ fixallatoms constrain_fixedatoms freeallatoms
  MD++ chainlength = [expr $ncpu-1] 

  MD++ timestep = 0.07 printfreq = 2
  #MD++ timestep = 0.1 printfreq = 2

  # make left end not movable
  #MD++ nebspec = \[ 0 1 0 1 1 0 \] totalsteps = $NEBMAXSTEP equilsteps = $EQUILSTEP

  #set NEBSTEP1 1000
  set NEBSTEP1 10
  set NEBSTEP1 0

  # set path length constant
  #MD++ nebspec = \[ 0 1 0 1 1 2 \] totalsteps = $NEBMAXSTEP equilsteps = $EQUILSTEP
  #MD++ nebspec = \[ 0 1 0 1 1 1 10 \] totalsteps = $NEBMAXSTEP equilsteps = $EQUILSTEP
  #MD++ nebspec = \[ 0 1 0 1 1 1 8 \] totalsteps = $NEBMAXSTEP equilsteps = $EQUILSTEP
  #MD++ nebspec = \[ 0 1 0 1 1 1 6 \] totalsteps = $NEBMAXSTEP equilsteps = $EQUILSTEP
  MD++ nebspec = \[ 0 1 0 1 1 1 7 \] totalsteps = $NEBMAXSTEP equilsteps = $EQUILSTEP
  #MD++ nebspec = \[ 0 1 0 1 1 1 \] totalsteps = $NEBSTEP1 equilsteps = $EQUILSTEP
  ## [ relax_surround redistr_freq moveleftend moverightend yesclimbimage islengthconstant ]

  MD++ Slave_chdir #Slave change to same directory as Master (dirname)
  MD++ allocchain_parallel initRchain_parallel

  if { $exist > 0 } {
    MD++ incnfile = stringrelax.chain.cn readRchain_parallel
  }

  MD++ eval
  if { $n6 < 3 } {
    MD++ stringrelax_parallel_2
  } else {
    MD++ stringrelax_parallel_3
  }
  MD++ finalcnfile = stringrelax.chain.cn writeRchain_parallel
  set check [ glob -nocomplain "stringeng_step*" ]
  set nfile [ llength $check]
  exec cp stringeng.out stringeng_step[expr $nfile+1].out

#
#  #MD++ nebspec = \[ 0 1 0 1 1 0 \] totalsteps = $NEBMAXSTEP equilsteps = $EQUILSTEP
#  MD++ nebspec = \[ 0 1 0 1 1 1 \] totalsteps = $NEBSTEP1 equilsteps = $EQUILSTEP
#  MD++ Broadcast_nebsetting
#  MD++ stringrelax_parallel
#  #MD++ stringrelax_parallel_1
#  MD++ finalcnfile = stringrelax.chain.cn writeRchain_parallel
#  exec cp stringeng.out stringeng_step[expr $nfile+2].out

  MD++ quit_all

} elseif { $status == 112.5 } {
  # Parallel NEB relaxation (for structures relaxed under stress)
  # to compile
  #  make eam build=R SYS=mc2_mpich SPEC+=-DMAXCONSTRAINATOMS=800001
  # to run (in PBS file)
  #  mpiexec -np $ncpu bin1/eam_mc2_mpich scripts/work/ni_disl/ni_screw_cross_slip 112.5 5000 0 0 0

  set NEBMAXSTEP $n
  #set NEBMAXSTEP 270 
  set EQUILSTEP  [expr $NEBMAXSTEP - 100 ]
  if { $EQUILSTEP < 0 } {
     set EQUILSTEP 0
  }

  set sigma_xy   0.0
  set sigma_yz   $n2
  set sigma_zz   $n3
  set sigma_yy   $n4
  set sigma_xz   $n5

 if { $sigma_yy != 0 } {
  if { $sigma_yy <= 100 } {
     set iselect 13
  } elseif { $sigma_yy <= 300 } {
     set iselect 10
  } elseif { $sigma_yy <= 600 } {
     set iselect 10
  } else {
     set iselect 10
     #set iselect 6
  }
 }
  # Escaig stress
 if { $sigma_yz != 0 } {
  if { $sigma_yz <= 300 } {
     set iselect 20
  } elseif { $sigma_yz <= 500 } {
     set iselect 10
  } elseif { $sigma_yz <= 600 } {
     set iselect 8
  } else {
     set iselect 8
     #set iselect 5
  }
 }
  # Schmid stress
 if { $sigma_xz != 0 } {
  if { $sigma_xz <= 100 } {
     set iselect 13
  } elseif { $sigma_xz <= 300 } {
     set iselect 11
  } elseif { $sigma_xz <= 600 } {
     set iselect 10
  } else {
     #set iselect 5
     set iselect 10 
  }
 }

  # For current Fleischer stuff:
  #set iselect 3 # Original Fleischer Work
  set iselect 2
  #set iselect 1

  puts "iselect = $iselect"

  set ncpu [ MD++_Get numDomains ]
  if { $ncpu < 1 } { set ncpu 1 }
  puts "ncpu = $ncpu"
  
  MD++ setoverwrite 
  MD++ dirname = runs/ni-Rao-cs-${status}_${n}_${n2}_${n3}_${n4}_${n5}

  MD++ NNM = 200
  MD++ writeall = 1

  exec ln -f -s ../ni-Rao-screw-30x30x20-5.5/stress_${n2}_${n3}_${n4}_${n5}/ni-screw-gp-relax-${n2}_${n3}_${n4}_${n5}.cn NEBinit.cn
  exec ln -f -s ../ni-Rao-screw-30x30x20-6.5/stress_${n2}_${n3}_${n4}_${n5}/ni-cross-slip-relax-${iselect}.cn NEBfinal.cn

  # use new A and B'
  #exec ln -f -s chain_no_7.cn NEBinit.cn
  #exec ln -f -s chain_no_23.cn NEBfinal.cn


  # rename previous result files in order to not overwrite them
  set check [ glob -nocomplain "stringrelax.chain.cn.cpu00" ]
  set exist [ llength $check ]
  if { $exist > 0 } {
     for { set i 0 } { $i < $ncpu } { incr i 1 } {
         set ival [ format %02d $i ]
         exec cp stringrelax.chain.cn.cpu$ival Prev_stringrelax.chain.cn.cpu$ival
     }
     if { [ llength [ glob -nocomplain "stringeng.out" ] ] > 0 } { 
       exec cp stringeng.out Prev_stringeng.out
     }
  }

  readpot
  MD++ Broadcast_EAM_Param

  MD++ incnfile = "../ni-Rao-screw-30x30x20-0.5/perf.cn" readcn  writeall = 1
  set NP0 [MD++_Get "NP"]

  # a test for new state A (0_0_0_400)
  #set chain_no 6
  #MD++ incnfile = NEBinit.cn  readcn 
  #MD++ incnfile = "../ni-cs-12_${n}_${n2}_${n3}_${n4}_${n5}.sav/stringrelax.chain.cn"
  #MD++ input = $chain_no  readRchain_parallel_toCN  RHtoS  setconfig1 #A

  #MD++ incnfile = NEBfinal.cn readcn setconfig2                       #B
  #MD++ incnfile = "../ni-cs-12_${n}_${n2}_${n3}_${n4}_${n5}.sav/stringrelax.chain.cn"
  #MD++ input = $chain_no  readRchain_parallel_toCN  RHtoS             #A

  MD++ incnfile = NEBinit.cn  readcn setconfig1 #A
  MD++ incnfile = NEBfinal.cn readcn setconfig2 #B
  MD++ incnfile = NEBinit.cn  readcn            #A

  MD++ freeallatoms

  setup_window

  set NP [MD++_Get "NP"]
  set vacuumratio  [expr 1.0 - ($NP*1.0)/$NP0]
  puts "vacuumratio = $vacuumratio"
  MD++ vacuumratio = $vacuumratio

  if { 1 } { # determine number of group 1/2 atoms 
  #Top surface atoms selection
  MD++ plot_color_axis = 2 eval # 2 means CSD
  MD++ input = \[ 1 -10 10 0.2 10 -10 10 10 100 \] fixatoms_by_pos_topol
  MD++ input = 1 setfixedatomsgroup 
  set N1 [MD++_Get "NPfixed"]
  MD++ freeallatoms
  #Bottom surface atoms selection
  MD++ input = \[ 1 -10 10 -10 -0.2 -10 10 10 100 \] fixatoms_by_pos_topol
  MD++ input = 2 setfixedatomsgroup 
  set N2 [MD++_Get "NPfixed"]
  MD++ freeallatoms
  } else {
  # the number of atoms in each group should equal to the following
    set N1 4800
    set N2 4800
  }

  set lx [MD++_Get "H_11"]
  set lz [MD++_Get "H_33"]

  #Compute Fext in terms of eV/A
  set f1x [expr -${sigma_xy}*$lx*$lz/$N1*0.6242E-05]
  set f1y [expr -${sigma_yy}*$lx*$lz/$N1*0.6242E-05]
  set f1z [expr -${sigma_yz}*$lx*$lz/$N1*0.6242E-05]

  set f2x [expr  ${sigma_xy}*$lx*$lz/$N1*0.6242E-05]
  set f2y [expr  ${sigma_yy}*$lx*$lz/$N1*0.6242E-05]
  set f2z [expr  ${sigma_yz}*$lx*$lz/$N1*0.6242E-05]

  MD++ extforce = \[ 2   $f1x $f1y $f1z  $f2x $f2y $f2z \]

  #MD++ freeallatoms clearR0 refreshnnlist eval  quit_all

  # do we need to broadcast extforce ?
  MD++ Broadcast_Atoms

  MD++ fixallatoms constrain_fixedatoms freeallatoms
  MD++ chainlength = [expr $ncpu-1] 

  MD++ timestep = 0.07 printfreq = 2
  #MD++ timestep = 0.1 printfreq = 2

  # make left end not movable
  #MD++ nebspec = \[ 0 1 0 1 1 0 \] totalsteps = $NEBMAXSTEP equilsteps = $EQUILSTEP

  #set NEBSTEP1 1000
  set NEBSTEP1 10
  set NEBSTEP1 0

  # set path length constant
  #MD++ nebspec = \[ 0 1 0 1 1 2 \] totalsteps = $NEBMAXSTEP equilsteps = $EQUILSTEP
  #MD++ nebspec = \[ 0 1 0 1 1 1 10 \] totalsteps = $NEBMAXSTEP equilsteps = $EQUILSTEP
  MD++ nebspec = \[ 0 1 0 1 1 1 10 \] totalsteps = $NEBMAXSTEP equilsteps = $EQUILSTEP
  #MD++ nebspec = \[ 0 1 0 1 1 1 6 \] totalsteps = $NEBMAXSTEP equilsteps = $EQUILSTEP
  #MD++ nebspec = \[ 0 1 0 1 1 1 4 \] totalsteps = $NEBMAXSTEP equilsteps = $EQUILSTEP
  #MD++ nebspec = \[ 0 1 0 1 1 1 \] totalsteps = $NEBSTEP1 equilsteps = $EQUILSTEP
  ## [ relax_surround redistr_freq moveleftend moverightend yesclimbimage islengthconstant ]

  MD++ Slave_chdir #Slave change to same directory as Master (dirname)
  MD++ allocchain_parallel initRchain_parallel

  if { $exist > 0 } {
    MD++ incnfile = stringrelax.chain.cn readRchain_parallel
  }

  MD++ eval
  if { $n6 < 3 } {
    MD++ stringrelax_parallel_2
  } else {
    MD++ stringrelax_parallel_3
  }
  MD++ finalcnfile = stringrelax.chain.cn writeRchain_parallel
  set check [ glob -nocomplain "stringeng_step*" ]
  set nfile [ llength $check]
  exec cp stringeng.out stringeng_step[expr $nfile+1].out

#
#  #MD++ nebspec = \[ 0 1 0 1 1 0 \] totalsteps = $NEBMAXSTEP equilsteps = $EQUILSTEP
#  MD++ nebspec = \[ 0 1 0 1 1 1 \] totalsteps = $NEBSTEP1 equilsteps = $EQUILSTEP
#  MD++ Broadcast_nebsetting
#  MD++ stringrelax_parallel
#  #MD++ stringrelax_parallel_1
#  MD++ finalcnfile = stringrelax.chain.cn writeRchain_parallel
#  exec cp stringeng.out stringeng_step[expr $nfile+2].out

  MD++ quit_all

} elseif { $status == 215 } {
  # MD++ run of dislocation under applied yz shear to equilibrate at 300K
  #  make eam build=R SYS=mc2_mpich SPEC+=-DMAXCONSTRAINATOMS=800001
  # to run (in PBS file)
  #  mpiexec -np $ncpu bin1/eam_mc2_mpich scripts/work/ni_disl/ni_screw_cross_slip 112 5000 0 0 0

  set totalsteps $n
  puts "totalsteps = $totalsteps"
  
  set sigma_xy   0.0
  set sigma_yz   $n2
  set sigma_zz   0.0
  set sigma_yy   0.0
  set sigma_xz   0.0
  set sigma_xx   0.0
  
  set modulus_xx 450e3 ; #MPa
  set modulus_zz 450e3 ; #MPa
  set modulus_xz 100e3 ; #MPa
  
  set seed [clock clicks]
  MD++ setoverwrite 
  MD++ dirname = runs/ni-cs/ni-Rao-cs-${status}_${n}_${n2}_T${n6}

  MD++ NNM = 200
  MD++ writeall = 1
  
  
  if { $n3 == 1 } {
  	# If previous iterations did not relax the stress waves out, continue from the previously partially relaxed state.  Otherwise start over.
  	exec ln -f -s ../ni-Rao-cs-${status}_${n}_${n2}_T300/final.cn initial.cn
  } else {
#  	exec ln -f -s ../ni-Rao-screw-30x30x20-5/stress_0_0_0_212.132/ni-screw-gp-relax-0_0_0_212.132.cn initial.cn
  	exec ln -f -s ../ni-Rao-screw-30x30x20-5/stress_100_0_0_0/ni-screw-gp-relax-100_0_0_0.cn initial.cn
  }

  readpot

  MD++ incnfile = "../ni-Rao-screw-30x30x20-0/perf.cn" readcn  writeall = 1
  set NP0 [MD++_Get "NP"]

  MD++ incnfile = initial.cn  readcn setconfig1 #A

  MD++ freeallatoms

  #setup_window
  #openwindow

  set NP [MD++_Get "NP"]
  set vacuumratio  [expr 1.0 - ($NP*1.0)/$NP0]
  puts "vacuumratio = $vacuumratio"
  MD++ vacuumratio = $vacuumratio

  if { 1 } { # determine number of group 1/2 atoms 
  #Top surface atoms selection
  MD++ plot_color_axis = 2 eval # 2 means CSD
  MD++ input = \[ 1 -10 10 0.2 10 -10 10 10 100 \] fixatoms_by_pos_topol
  MD++ input = 1 setfixedatomsgroup 
  set N1 [MD++_Get "NPfixed"]
  MD++ freeallatoms
  #Bottom surface atoms selection
  MD++ input = \[ 1 -10 10 -10 -0.2 -10 10 10 100 \] fixatoms_by_pos_topol
  MD++ input = 2 setfixedatomsgroup 
  set N2 [MD++_Get "NPfixed"]
  MD++ freeallatoms
  } else {
  # the number of atoms in each group should equal to the following
    set N1 4800
    set N2 4800
  }

  set lx [MD++_Get "H_11"]
  set lz [MD++_Get "H_33"]

  #Compute Fext in terms of eV/A
  set f1x [expr -${sigma_xy}*$lx*$lz/$N1*0.6242E-05]
  set f1y [expr -${sigma_yy}*$lx*$lz/$N1*0.6242E-05]
  set f1z [expr -${sigma_yz}*$lx*$lz/$N1*0.6242E-05]

  set f2x [expr  ${sigma_xy}*$lx*$lz/$N2*0.6242E-05]
  set f2y [expr  ${sigma_yy}*$lx*$lz/$N2*0.6242E-05]
  set f2z [expr  ${sigma_yz}*$lx*$lz/$N2*0.6242E-05]

  MD++ extforce = \[ 2   $f1x $f1y $f1z  $f2x $f2y $f2z \]
  
  MD++ eval
  
  append seed [clock clicks]
  append seed [clock seconds] [pid]
  
  setup_md 
  MD++ totalsteps = $totalsteps
  MD++ T_OBJ = 300
  MD++ ensemble_type = "NVTC" integrator_type = "VVerlet"
  MD++ NHChainLen = 4 NHMass = \[ 2e-3 2e-6 2e-6 2e-6 \]
  MD++ DOUBLE_T = 0 
  #MD++ randseed = 12345 srand48 initvelocity
  if { $n3 == 1 } {
  	# Do nothing
  } else {
  	MD++ srand48bytime initvelocity
  }
  
  MD++ savecnfreq = 10000
  
  MD++ eval
  
  # Each subiteration has 100 steps, and there are 20 subiterations
  set maxiter [expr $totalsteps / 200.0 ] 
  for { set iter 1 } { $iter <= $maxiter } { incr iter } {

    set sig_avg_xx 0
    set sig_avg_yy 0
    set sig_avg_zz 0
    set sig_avg_xy 0
    set sig_avg_yz 0
    set sig_avg_xz 0
    
    if { $iter <= 50 } {
          set factor 400.0e3
    } elseif { $iter <= 100 } {
    	  set factor 800.0e3
    } elseif { $iter <= 150 } {
          set factor 1200.0e3
    } elseif { $iter <= 200 } {
    	  set factor 1600.0e3
    } else {
    	  set factor 2000.0e3
    }
    
    

    for { set subiter 1 } { $subiter <= 20 } { incr subiter} {
       MD++ continue_curstep = 1
       MD++ totalsteps = 10 run
       set sig_xx [MD++_Get TSTRESSinMPa_xx]
       set sig_yy [MD++_Get TSTRESSinMPa_yy]
       set sig_zz [MD++_Get TSTRESSinMPa_zz]
       set sig_xy [MD++_Get TSTRESSinMPa_xy]
       set sig_yz [MD++_Get TSTRESSinMPa_yz]
       set sig_xz [MD++_Get TSTRESSinMPa_xz]
       puts "sig_xx = $sig_xx sig_yy = $sig_yy sig_zz = $sig_zz" 
       puts "sig_yz = $sig_yz sig_xz = $sig_xz sig_xy = $sig_xy"

       if { $subiter > 10 } {
          set sig_avg_xx [expr $sig_avg_xx + $sig_xx]
          set sig_avg_yy [expr $sig_avg_yy + $sig_yy]
          set sig_avg_zz [expr $sig_avg_zz + $sig_zz]
          set sig_avg_xy [expr $sig_avg_xy + $sig_xy]
          set sig_avg_yz [expr $sig_avg_yz + $sig_yz]
          set sig_avg_xz [expr $sig_avg_xz + $sig_xz]
       }
    }
    set sig_avg_xx [expr $sig_avg_xx / 10]
    set sig_avg_yy [expr $sig_avg_yy / 10]
    set sig_avg_zz [expr $sig_avg_zz / 10]
    set sig_avg_xy [expr $sig_avg_xy / 10]
    set sig_avg_yz [expr $sig_avg_yz / 10]
    set sig_avg_xz [expr $sig_avg_xz / 10]

    puts "sig_avg_xx = $sig_avg_xx sig_avg_yy = $sig_avg_yy sig_avg_zz = $sig_avg_zz" 
    puts "sig_avg_yz = $sig_avg_yz sig_avg_xz = $sig_avg_xz sig_avg_xy = $sig_avg_xy"
    if { $iter < $maxiter } {
       set xx_err [ expr ($sig_avg_xx - $sigma_xx) ]
       set yy_err [ expr ($sig_avg_yy - $sigma_yy) ]
       set zz_err [ expr ($sig_avg_zz - $sigma_zz) ]
       set xz_err [ expr ($sig_avg_xz - $sigma_xz) ]
       # Newly added to adjust applied force
       set yz_err [ expr ($sig_avg_yz - $sigma_yz) ]
       set delta_xx [ expr $xx_err / 400.0e3 ]
       set delta_yy [ expr $yy_err / 400.0e3 ]
       set delta_zz [ expr $zz_err / 400.0e3 ]
       set delta_xz [ expr $xz_err / 400.0e3 ]
       MD++ input = \[ 1 1 $delta_xx \] changeH_keepS
       MD++ input = \[ 2 2 $delta_yy \] changeH_keepS
       MD++ input = \[ 3 3 $delta_zz \] changeH_keepS
       MD++ input = \[ 3 1 $delta_xz \] changeH_keepS
       swap_velocities $NP 20
    }
  }
  
  MD++ finalcnfile = final.cn writeall = 1 writecn

  MD++ quit_all

} elseif { $status == 216 } {
  # MD++ run to equilibrate temperature and stresses
  # to compile
  #  make eam build=R SYS=mc2_mpich SPEC+=-DMAXCONSTRAINATOMS=800001
  # to run (in PBS file)
  #  mpiexec -np $ncpu bin1/eam_mc2_mpich scripts/work/ni_disl/ni_screw_cross_slip 112 5000 0 0 0

  set totalsteps $n
  puts "totalsteps = $totalsteps"
  
  set sigma_xy   0.0
  set sigma_yz   $n2
  set sigma_zz   $n3
  set sigma_yy   $n4
  set sigma_xz   $n5
  set sigma_xx   0.0
  
  set modulus_xx 450e3 ; #MPa
  set modulus_zz 450e3 ; #MPa
  set modulus_xz 100e3 ; #MPa
  
  set seed [clock clicks]
  MD++ setoverwrite 
#  MD++ dirname = /data/Cai-group/Codes/MD++.svn/runs/ni_cs_Xiaohan/ni_cs_rao99_md_215_216/ni-Rao-cs-${status}_${n}_${n2}_${n3}_${n4}_${n5}_T${n6}
  MD++ dirname = /home/xzhang11/MD++/runs/ni-cs/ni_cs_rao99_md_215_216/ni-Rao-cs-${status}_${n}_${n2}_${n3}_${n4}_${n5}_T${n6}

  MD++ NNM = 200
  MD++ writeall = 1

#  exec ln -f -s /data/Cai-group/Codes/MD++.svn/runs/ni_cs_Xiaohan/ni_cs_rao99_md_215_216/ni-Rao-cs-215_${n}_${n2}_T300/final.cn initial.cn
  exec ln -f -s /home/xzhang11/MD++/runs/ni-cs/ni_cs_rao99_md_215_216/ni-Rao-cs-215_${n}_${n2}_T300/final.cn initial.cn

  readpot

#  MD++ incnfile = "/data/Cai-group/Codes/MD++.svn/runs/ni_cs_Xiaohan/ni_cs_rao99_md_215_216/ni-Rao-screw-30x30x20-0/perf.cn" readcn  writeall = 1
  MD++ incnfile = "/home/xzhang11/MD++/runs/ni-cs/ni_cs_rao99_md_215_216/ni-Rao-screw-30x30x20-0/perf.cn" readcn  writeall = 1

  set NP0 [MD++_Get "NP"]

  MD++ incnfile = initial.cn  readcn setconfig1 #A

  if { $USEGPU == 1 } {
    MD++ cuda_memcpy_all
  }

  MD++ freeallatoms

  #setup_window
  #openwindow

  set NP [MD++_Get "NP"]
  set vacuumratio  [expr 1.0 - ($NP*1.0)/$NP0]
  puts "vacuumratio = $vacuumratio"
  MD++ vacuumratio = $vacuumratio

  if { 1 } { # determine number of group 1/2 atoms 
  #Top surface atoms selection
  puts "I am here 2.1"
  MD++ plot_color_axis = 2 
  puts "I am here 2.2"

  MD++ eval # 2 means CSD
  puts "I am here 2"
  MD++ input = \[ 1 -10 10 0.2 10 -10 10 10 100 \] fixatoms_by_pos_topol
  MD++ input = 1 setfixedatomsgroup 
  set N1 [MD++_Get "NPfixed"]
  puts "I am here 3"
  MD++ freeallatoms
  #Bottom surface atoms selection
  MD++ input = \[ 1 -10 10 -10 -0.2 -10 10 10 100 \] fixatoms_by_pos_topol
  puts "I am here 3.1"
  MD++ input = 2 setfixedatomsgroup 
  set N2 [MD++_Get "NPfixed"]
  MD++ freeallatoms
  } else {
  # the number of atoms in each group should equal to the following
    set N1 4800
    set N2 4800
  }

  puts "I am here 4"

  set lx [MD++_Get "H_11"]
  set lz [MD++_Get "H_33"]

  #Compute Fext in terms of eV/A
  set f1x [expr -${sigma_xy}*$lx*$lz/$N1*0.6242E-05]
  set f1y [expr -${sigma_yy}*$lx*$lz/$N1*0.6242E-05]
  set f1z [expr -${sigma_yz}*$lx*$lz/$N1*0.6242E-05]

  set f2x [expr  ${sigma_xy}*$lx*$lz/$N2*0.6242E-05]
  set f2y [expr  ${sigma_yy}*$lx*$lz/$N2*0.6242E-05]
  set f2z [expr  ${sigma_yz}*$lx*$lz/$N2*0.6242E-05]

  MD++ extforce = \[ 2   $f1x $f1y $f1z  $f2x $f2y $f2z \]
  puts "I am here 5"
  
  MD++ eval
  
  setup_md
  #set totalsteps [expr $totalsteps / 5.0 ];
  MD++ totalsteps = $totalsteps
  MD++ T_OBJ = $n6
  MD++ ensemble_type = "NVTC" integrator_type = "VVerlet"
  MD++ NHChainLen = 4 NHMass = \[ 2e-3 2e-6 2e-6 2e-6 \]
  MD++ DOUBLE_T = 0 
  #MD++ randseed = 12345 srand48 initvelocity
  #MD++ srand48bytime initvelocity
  
  MD++ savecnfreq = 5000
  
  MD++ eval
  
  append seed [clock clicks]
  append seed [clock seconds] [pid]
  
  # Lower totalsteps because everything is converging faster
  set totalsteps [expr $totalsteps / 4]
 #????????????????????????????????????????????????????????????????????????????????????????? 
  set totalsteps 400 
 #????????????????????????????????????????????????????????????????????????????????????????? 
  # Each subiteration has 10 steps, and there are 20 subiterations
  set maxiter [expr $totalsteps / 200.0 ]
  puts "maxiter = $maxiter"
  set delta_T [expr ($n6-300)*2 / $maxiter ]
  puts "delta_T = $delta_T"
  set newT 300

  set runC 0
  if { $runC == 1 } { 
     MD++ wrapper_run
  } else { 
  #swap_velocities $NP 1000
  for { set iter 1 } { $iter <= $maxiter } { incr iter } {

    set sig_avg_xx 0
    set sig_avg_yy 0
    set sig_avg_zz 0
    set sig_avg_xy 0
    set sig_avg_yz 0
    set sig_avg_xz 0
    
    if { $iter <= 50 } {
          set factor 400.0e3
    } elseif { $iter <= 100 } {
    	  set factor 800.0e3
    } elseif { $iter <= 150 } {
          set factor 1200.0e3
    } elseif { $iter <= 200 } {
    	  set factor 1600.0e3
    } else {
    	  set factor 2000.0e3
    }
    
    if { $newT < $n6 } {
    	set newT [expr $newT + $delta_T]
    	MD++ T_OBJ = $newT
    }
    
    for { set subiter 1 } { $subiter <= 20 } { incr subiter} {
       MD++ continue_curstep = 1
       MD++ totalsteps = 10 run
       set sig_xx [MD++_Get TSTRESSinMPa_xx]
       set sig_yy [MD++_Get TSTRESSinMPa_yy]
       set sig_zz [MD++_Get TSTRESSinMPa_zz]
       set sig_xy [MD++_Get TSTRESSinMPa_xy]
       set sig_yz [MD++_Get TSTRESSinMPa_yz]
       set sig_xz [MD++_Get TSTRESSinMPa_xz]
       puts "sig_xx = $sig_xx sig_yy = $sig_yy sig_zz = $sig_zz" 
       puts "sig_yz = $sig_yz sig_xz = $sig_xz sig_xy = $sig_xy"

       if { $subiter > 10 } {
          set sig_avg_xx [expr $sig_avg_xx + $sig_xx]
          set sig_avg_yy [expr $sig_avg_yy + $sig_yy]
          set sig_avg_zz [expr $sig_avg_zz + $sig_zz]
          set sig_avg_xy [expr $sig_avg_xy + $sig_xy]
          set sig_avg_yz [expr $sig_avg_yz + $sig_yz]
          set sig_avg_xz [expr $sig_avg_xz + $sig_xz]
       }
    }
    set sig_avg_xx [expr $sig_avg_xx / 10]
    set sig_avg_yy [expr $sig_avg_yy / 10]
    set sig_avg_zz [expr $sig_avg_zz / 10]
    set sig_avg_xy [expr $sig_avg_xy / 10]
    set sig_avg_yz [expr $sig_avg_yz / 10]
    set sig_avg_xz [expr $sig_avg_xz / 10]

    puts "sig_avg_xx = $sig_avg_xx sig_avg_yy = $sig_avg_yy sig_avg_zz = $sig_avg_zz" 
    puts "sig_avg_yz = $sig_avg_yz sig_avg_xz = $sig_avg_xz sig_avg_xy = $sig_avg_xy"
    if { $iter < $maxiter } {
       set xx_err [ expr ($sig_avg_xx - $sigma_xx) ]
       set yy_err [ expr ($sig_avg_yy - $sigma_yy) ]
       set zz_err [ expr ($sig_avg_zz - $sigma_zz) ]
       set xz_err [ expr ($sig_avg_xz - $sigma_xz) ]
       # Newly added to adjust applied force
       set yz_err [ expr ($sig_avg_yz - $sigma_yz) ]
       set delta_xx [ expr $xx_err / $factor ]
       set delta_yy [ expr $yy_err / $factor ]
       set delta_zz [ expr $zz_err / $factor ]
       set delta_xz [ expr $xz_err / $factor ]
       MD++ input = \[ 1 1 $delta_xx \] changeH_keepS
       MD++ input = \[ 2 2 $delta_yy \] changeH_keepS
       MD++ input = \[ 3 3 $delta_zz \] changeH_keepS
       MD++ input = \[ 3 1 $delta_xz \] changeH_keepS
       swap_velocities $NP 20
    }
  }
  }
  
  MD++ finalcnfile = final.cn writeall = 1 writecn
  MD++ finalcnfile = final.cfg writeatomeyecfg
  MD++ quit_all

} elseif { $status == 217 } {
  # MD++ run to see if cross-slip occurs
  # to compile
  #  make eam build=R SYS=mc2_mpich SPEC+=-DMAXCONSTRAINATOMS=800001
  # to run (in PBS file)
  #  mpiexec -np $ncpu bin1/eam_mc2_mpich scripts/work/ni_disl/ni_screw_cross_slip 112 5000 0 0 0

  set totalsteps $n
  puts "totalsteps = $totalsteps"
  set totalsteps [ expr $n * 40 ] 
  set sigma_xy   0.0
  set sigma_yz   $n2
  set sigma_zz   $n3
  set sigma_yy   $n4
  set sigma_xz   $n5
  set sigma_xx   0.0
  
  set modulus_xx 450e3 ; #MPa
  set modulus_zz 450e3 ; #MPa
  set modulus_xz 100e3 ; #MPa
  
  set seed [clock clicks]
  MD++ setoverwrite 
# MD++ dirname = /data/Cai-group/Codes/MD++.svn/runs/ni_cs_Xiaohan/ni_cs_rao99_md_217/ni-Rao-cs-${status}_${n}_${n2}_${n3}_${n4}_${n5}_T${n6}_${n7}
  MD++ dirname = /scratch/xzhang11/ni_cs_rao99_md_217/ni-Rao-cs-${status}_${n}_${n2}_${n3}_${n4}_${n5}_T${n6}_${n7}

  MD++ NNM = 200
  MD++ writeall = 1

#  exec ln -f -s /data/Cai-group/Codes/MD++.svn/runs/ni_cs_Xiaohan/ni_cs_rao99_md_215_216/ni-Rao-cs-216_${n}_${n2}_${n3}_${n4}_${n5}_T${n6}/final.cn initial.cn
  exec ln -f -s /scratch/xzhang11/ni_cs_rao99_md_215_216/ni-Rao-cs-216_${n}_${n2}_${n3}_${n4}_${n5}_T${n6}/final.cn initial.cn
  readpot

#  MD++ incnfile = "/data/Cai-group/Codes/MD++.svn/runs/ni_cs_Xiaohan/ni_cs_rao99_md_215_216/ni-Rao-screw-30x30x20-0/perf.cn" readcn  writeall = 1
  MD++ incnfile = "/scratch/xzhang11/ni_cs_rao99_md_215_216/ni-Rao-screw-30x30x20-0/perf.cn" readcn  writeall = 1
  set NP0 [MD++_Get "NP"]

  MD++ incnfile = initial.cn  readcn setconfig1 #A

  MD++ freeallatoms

  # setup_window
  #openwindow

  set NP [MD++_Get "NP"]
  set vacuumratio  [expr 1.0 - ($NP*1.0)/$NP0]
  puts "vacuumratio = $vacuumratio"
  MD++ vacuumratio = $vacuumratio

  if { 1 } { # determine number of group 1/2 atoms 
  #Top surface atoms selection
  MD++ plot_color_axis = 2 eval # 2 means CSD
  MD++ input = \[ 1 -10 10 0.2 10 -10 10 10 100 \] fixatoms_by_pos_topol
  MD++ input = 1 setfixedatomsgroup 
  set N1 [MD++_Get "NPfixed"]
  MD++ freeallatoms
  #Bottom surface atoms selection
  MD++ input = \[ 1 -10 10 -10 -0.2 -10 10 10 100 \] fixatoms_by_pos_topol
  MD++ input = 2 setfixedatomsgroup 
  set N2 [MD++_Get "NPfixed"]
  MD++ freeallatoms
  } else {
  # the number of atoms in each group should equal to the following
    set N1 4800
    set N2 4800
  }

  set lx [MD++_Get "H_11"]
  set lz [MD++_Get "H_33"]

  #Compute Fext in terms of eV/A
  set f1x [expr -${sigma_xy}*$lx*$lz/$N1*0.6242E-05]
  set f1y [expr -${sigma_yy}*$lx*$lz/$N1*0.6242E-05]
  set f1z [expr -${sigma_yz}*$lx*$lz/$N1*0.6242E-05]

  set f2x [expr  ${sigma_xy}*$lx*$lz/$N2*0.6242E-05]
  set f2y [expr  ${sigma_yy}*$lx*$lz/$N2*0.6242E-05]
  set f2z [expr  ${sigma_yz}*$lx*$lz/$N2*0.6242E-05]
  
  set sigma_yz_adj $sigma_yz

  MD++ extforce = \[ 2   $f1x $f1y $f1z  $f2x $f2y $f2z \]
  
  MD++ eval
  
  append seed [clock clicks]
  append seed [clock seconds] [pid]
  #append seed [winfo id .] [winfo geometry .] [winfo pointerxy .]
  #set hashedseed [binary format H* [sha1::sha1 $seed]]
  #srand($seed)
  set atom [expr { int(srand($seed)*$NP) } ]
  puts "atom = $atom"
  #for { set i 0 } { $i < 100 } {incr i } {
  #   set atom [expr { int(rand()*$NP) } ]
  #   puts "atom = $atom"
  #}
  
  puts "Perturb velocities ..."
  set maxAtomsPerturbed  [expr ($NP*1.0)/1000]
  set vx [MD++_Get VSR(0)]
  set vy [MD++_Get VSR(1)]
  set vz [MD++_Get VSR(2)]
  for { set i 0 } { $i < $maxAtomsPerturbed } { incr i } {
     set atom [expr { int(rand()*$NP) } ]
     set ind [expr $atom*3]
     set vx0 [MD++_Get VSR($ind)]
     MD++_Set VSR($ind) $vx
     set ind [expr $atom*3+1]
     set vy0 [MD++_Get VSR($ind)]
     MD++_Set VSR($ind) $vy
     set ind [expr $atom*3+2]
     set vz0 [MD++_Get VSR($ind)]
     MD++_Set VSR($ind) $vz
     puts "vx = $vx, vx0 = $vx0"
     puts "vy = $vy, vy0 = $vy0"
     puts "vz = $vz, vz0 = $vz0"
     set vx $vx0
     set vy $vy0
     set vz $vz0
  }
  MD++_Set VSR(0) $vx
  MD++_Set VSR(1) $vy
  MD++_Set VSR(2) $vz
  
  setup_md
  #set totalsteps [expr $totalsteps / 5.0 ];
  MD++ totalsteps = $totalsteps
  MD++ T_OBJ = $n6
  MD++ ensemble_type = "NVTC" integrator_type = "VVerlet"
  MD++ NHChainLen = 4 NHMass = \[ 2e-3 2e-6 2e-6 2e-6 \]
  MD++ DOUBLE_T = 0 
  #MD++ randseed = 12345 srand48 initvelocity
  #MD++ srand48bytime initvelocity
  
  MD++ savecnfreq = 10000
  
  MD++ eval
  
  # Each subiteration has 20 steps, and there are (totalsteps/200) subiterations
  set maxiter [expr $totalsteps / 200.0 ]
  puts "maxiter = $maxiter"
  for { set iter 1 } { $iter <= $maxiter } { incr iter } {

    set sig_avg_xx 0
    set sig_avg_yy 0
    set sig_avg_zz 0
    set sig_avg_xy 0
    set sig_avg_yz 0
    set sig_avg_xz 0
    
    if { $iter <= 50 } {
          set factor 400.0e3
    } elseif { $iter <= 100 } {
    	  set factor 800.0e3
    } elseif { $iter <= 150 } {
          set factor 1200.0e3
    } elseif { $iter <= 200 } {
    	  set factor 1600.0e3
    } else {
    	  set factor 2000.0e3
    }
    

    for { set subiter 1 } { $subiter <= 20 } { incr subiter} {
       MD++ continue_curstep = 1
       MD++ totalsteps = 10 run
       set sig_xx [MD++_Get TSTRESSinMPa_xx]
       set sig_yy [MD++_Get TSTRESSinMPa_yy]
       set sig_zz [MD++_Get TSTRESSinMPa_zz]
       set sig_xy [MD++_Get TSTRESSinMPa_xy]
       set sig_yz [MD++_Get TSTRESSinMPa_yz]
       set sig_xz [MD++_Get TSTRESSinMPa_xz]
    #   puts "sig_xx = $sig_xx sig_yy = $sig_yy sig_zz = $sig_zz" 
    #   puts "sig_yz = $sig_yz sig_xz = $sig_xz sig_xy = $sig_xy"

       if { $subiter > 10 } {
          set sig_avg_xx [expr $sig_avg_xx + $sig_xx]
          set sig_avg_yy [expr $sig_avg_yy + $sig_yy]
          set sig_avg_zz [expr $sig_avg_zz + $sig_zz]
          set sig_avg_xy [expr $sig_avg_xy + $sig_xy]
          set sig_avg_yz [expr $sig_avg_yz + $sig_yz]
          set sig_avg_xz [expr $sig_avg_xz + $sig_xz]
       }
    }
    set sig_avg_xx [expr $sig_avg_xx / 10]
    set sig_avg_yy [expr $sig_avg_yy / 10]
    set sig_avg_zz [expr $sig_avg_zz / 10]
    set sig_avg_xy [expr $sig_avg_xy / 10]
    set sig_avg_yz [expr $sig_avg_yz / 10]
    set sig_avg_xz [expr $sig_avg_xz / 10]

    #puts "sig_avg_xx = $sig_avg_xx sig_avg_yy = $sig_avg_yy sig_avg_zz = $sig_avg_zz" 
    #puts "sig_avg_yz = $sig_avg_yz sig_avg_xz = $sig_avg_xz sig_avg_xy = $sig_avg_xy"
    if { $iter < $maxiter } {
       set xx_err [ expr ($sig_avg_xx - $sigma_xx) ]
       set yy_err [ expr ($sig_avg_yy - $sigma_yy) ]
       set zz_err [ expr ($sig_avg_zz - $sigma_zz) ]
       set xz_err [ expr ($sig_avg_xz - $sigma_xz) ]
       # Newly added to adjust applied force
       set yz_err [ expr ($sig_avg_yz - $sigma_yz) ]
       set delta_xx [ expr $xx_err / 400.0e3 ]
       set delta_yy [ expr $yy_err / 400.0e3 ]
       set delta_zz [ expr $zz_err / 400.0e3 ]
       set delta_xz [ expr $xz_err / 400.0e3 ]
       MD++ input = \[ 1 1 $delta_xx \] changeH_keepS
       MD++ input = \[ 2 2 $delta_yy \] changeH_keepS
       MD++ input = \[ 3 3 $delta_zz \] changeH_keepS
       MD++ input = \[ 3 1 $delta_xz \] changeH_keepS
    }
  }
  
  MD++ finalcnfile = final.cn writeall = 1 writecn

  MD++ quit_all


} elseif { $status == 2017 } {
  # MD++ run to see if cross-slip occurs
  # to compile
  #  make eam build=R SYS=mc2_mpich SPEC+=-DMAXCONSTRAINATOMS=800001
  # to run (in PBS file)
  #  mpiexec -np $ncpu bin1/eam_mc2_mpich scripts/work/ni_disl/ni_screw_cross_slip 112 5000 0 0 0
  
  readpot

#MD++ incnfile = "runs/ni-cs/ni-Rao-cs-216_100000_400_-1449.5689_1449.5689_848.5281_T450/final.cn" readcn
#MD++ incnfile = "runs/ni-cs/ni-Rao-cs-216_100000_400_-1449.5689_1449.5689_848.5281_T500/final.cn" readcn
#MD++ incnfile = "runs/ni-cs/ni-Rao-cs-216_100000_400_-1449.5689_1449.5689_848.5281_T550/final.cn" readcn
#MD++ incnfile = "runs/ni-cs/ni-Rao-cs-216_100000_400_-1449.5689_1449.5689_848.5281_T600/final.cn" readcn
#MD++ incnfile = "runs/ni-cs/ni-Rao-cs-216_100000_200_-1520.2796_1520.2796_1060.6602_T450/final.cn" readcn
#MD++ incnfile = "runs/ni-cs/ni-Rao-cs-216_100000_200_-1520.2796_1520.2796_1060.6602_T500/final.cn" readcn
#MD++ incnfile = "runs/ni-cs/ni-Rao-cs-216_100000_200_-1520.2796_1520.2796_1060.6602_T550/final.cn" readcn
MD++ incnfile = "runs/ni-cs/ni-Rao-cs-216_100000_200_-1520.2796_1520.2796_1060.6602_T600/final.cn" readcn

  #MD++ incnfile = "runs/ni-cs/ni-Rao-cs-216_100000_400_-1449.5689_1449.5689_848.5281_T600/inter0003.cn" readcn  
  MD++ finalcnfile = "final.cfg" writeatomeyecfg

} else {
 	puts "unknown status = $status"
 	exitmd 

} 

