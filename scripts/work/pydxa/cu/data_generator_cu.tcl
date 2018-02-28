# -*-shell-script-*-
# sw_mc2 scripts/disl_nuc_hetero.tcl 0
source "$::env(MDPLUS_DIR)/scripts/Examples/Tcl/startup.tcl"

#*******************************************
# Definition of procedures
#*******************************************
proc initmd { status {T 0} {epVAL 0} {opt 0} {alt 0}  } {
MD++ setnolog
MD++ setoverwrite
MD++ dirname = runs/pydxa/cu/
MD++ atommass = 63.546 # (g/mol)
}

#------------------------------------------------------------
proc readpot { } { MD++ {
#--------------------------------------------
#Read in potential file
#
#potfile = $::env(MDPLUS_DIR)/potentials/w_pot readpot $$$$$$$$$$$$$$$$$$$$$$$
potfile = ~/Planet/Libs/MD++.svn/potentials/EAMDATA/eamdata.CuMishin
eamgrid = 5000 readeam NNM = 600
} }

# make sure the coordinate is right hand sided.
proc make_perfect_crystal { status nx ny nz } {
    MD++ crystalstructure = face-centered-cubic latticeconst = 3.615 #(A) for Cu
    if { $status == 0 } { 
      MD++ latticesize = \[  1 -2 1  $nx  1 1 1  $ny  1 0 -1  $nz \]
    } else { 
      puts "Undefined for status = $status"
      MD++ quit
    }
    MD++ makecrystal #finalcnfile = perf.cn writecn #eval
}

#0: Si. 1 : Ge.
proc set_all_atoms_species { id } {
    MD++ fixedallatoms
    MD++ input = $id setfixedatomsspecies
    MD++ freeallatoms
}

#--------------------------------------------
proc relax_fixbox { } { MD++ {
# Conjugate-Gradient relaxation
conj_ftol = 2e-6 conj_itmax = 1000 conj_fevalmax = 30
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
                   1 0 1
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
fixatomcolor = yellow backgroundcolor = gray70
#atomcolor = lightgrey highlightcolor = purple  bondcolor = darkgrey
plot_color_axis = 2  NCS = 4
plot_color_windows = [ 1
                       0.6 9.4   1  
                       9.4  10   5
                       10 20   6
                       20 50   8
                       0  0.6  4
                     ]

# plot_limits = [ 1 -10 10 -0.2 0.2 -10 10 ]

#
#xiaohan
plot_atom_info = 1 # reduced coordinates of atoms
#plot_atom_info = 2 # real coordinates of atoms
#plot_atom_info = 3 # energy of atoms
#plot_highlight = [ 0 0 1 2 3 4 5 6 7 8 9 ]
plotfreq = 10
rotateangles = [ -0 90 0 1.2 ]
#openwin alloccolors rotate saverot plot
#plot_color_axis = 0 input = [ -8 -3 10] GnuPlotHistogram
#plot_color_axis = 2 input = [ 0.6 50 50 ] GnuPlotHistogram
} }

proc openwindow { } { 
#setup_window
MD++ openwin alloccolors rotate saverot eval plot
}

#--------------------------------------------
proc exitmd { } { MD++ quit }
#end of proc exitmd
#--------------------------------------------

proc setup_md { } { MD++ {     
T_OBJ = 300 #Kelvin #add by xiaohan

equilsteps = 0  totalsteps = 5000 timestep = 0.0001 # (ps)
atommass = 28.0855 # (g/mol)
DOUBLE_T = 1
saveprop = 1 savepropfreq = 100 # openpropfile #run
savecn = 1 savecnfreq = 10000 openintercnfile
plotfreq = 10 printfreq = 100
ensemble_type = "NPTC" integrator_type = "Gear6" implementation_type = 0
#ensemble_type = "NVE" integrator_type = "VVerlet" implementation_type = 0
NHChainLen = 4 NHMass=[0.5e-1 1e-2 1e-2 1e-2 ]
vt2 = 1e28  #1e28 2e28 5e28
wallmass = 2e3     # atommass * NP = 14380
boxdamp = 1e-3     # optimal damping for 216 atoms and wallmass 1e-3
saveH # Use current H as reference (H0), needed for specifying stress

fixboxvec = [ 0 1 1 1 0 1 1 1 0 ]

output_fmt="curstep EPOT KATOM Tinst TSTRESSinMPa_xx TSTRESSinMPa_yy TSTRESSinMPa_zz TSTRESSinMPa_xy TSTRESSinMPa_yz TSTRESSinMPa_zx H_11 H_12 H_13 H_21 H_22 H_23 H_31 H_32 H_33 N_lgst_cluster"
} }
#end of proc setup_md

proc datafile_process { filename index frac fracend operation } {
   set fp [ open $filename r ]
   set data [ split [read $fp] \n]
   close $fp
   set NUM [ llength $data ]
   set NUM [ expr $NUM-1 ]
   set Nini [ expr round($NUM*$frac) ] 
   set Nend [ expr round($NUM*$fracend) ]
#   puts "total line $NUM \n"
   set Sum 0
   set k 0
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
#       puts "$datum"
       set Sum [expr $Sum+$datum]
   }
   set Ave [ expr $Sum/$k]
#   puts "$Sum $Ave"
#   puts "$MAX $MIN"

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
#   puts "$Std $Var" 
#   puts "$LIST"
   return $LIST
}


proc make_ellipse_dislocation_loop { x0 y0 z0 epsilon } {
    set Lx [MD++_Get H_11]
    set Ly [MD++_Get H_22]
    set Lz [MD++_Get H_33]
    set store 1
    set a 3.615
    set bx 0.3333
    set by 0
    set bz 0

    set bx 0
    set by 0
    set bz [expr -sqrt(2)/2 ]
    
# 0.5[0 -1 1] = 1/sqrt(6) * [1 -2 1] +  1/sqrt(3) * [1 1 1]  + 1/sqrt(2)*[1 0 -1] : Full
    set bx [expr 0.5*sqrt(6)/2]
    set by 0
    set bz [expr -0.5*sqrt(2)/2 ]

# 0.5[0 -1 1] = 1/6[-1 -1 2] + 1/6[1 -2 1]
# 1/6[ -1 -1 2 ] = 1/sqrt(6) * [1 -2 1] * (0.5)  +  1/sqrt(3) * [1 1 1] *(0)  + 1/sqrt(2)*[1 0 -1] *(-1.5)
    set bx [expr sqrt(6)*0.5/6.0]
    set by 0
    set bz [expr sqrt(2)*(-1.5)/6.0]

    set bx [expr 1/sqrt(6) ]
    set by 0
    set bz 0

    set lx 1
    set ly 0
    set lz 0
    set nx 0
    set ny 1
    set nz 0
    #set y 1.0435
    set x [ expr $x0*$Lx ]
    set y [ expr $y0*$Ly ]
    set z [ expr $z0*$Lz ] 
#0.0300 x
    #set epsilon 0.145

    set Ra [expr $Lx * $epsilon]
    set Rb [expr 0.8*$Ra]
    MD++ input= \[ 1 $Ra $Rb $a $bx $by $bz $lx $ly $lz $nx $ny $nz $x $y $z $store \] 
    MD++ makedislellipse  
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
 set n 300
} elseif { $argc > 1 } {
 set n [lindex $argv 1]
}
puts "n = $n"

if { $argc <= 2 } {
 set flag 0.110
} elseif { $argc > 2 } {
 set flag [lindex $argv 2]
}
puts "flag = $flag"

if { $argc <= 3 } {
 set opt 0
} elseif { $argc > 3 } { set opt [lindex $argv 3] }
puts "opt = $opt"

if { $status == 0 } {
# \[  1 -2 1  $nx  1 1 1  $ny  1 0 -1  $nz \]

  initmd $status
  readpot
  make_perfect_crystal $status 8 6 13
  MD++ finalcnfile = "0K_0.0_relaxed.cn" writecn
 # setup_window
 # openwindow
  set H11_0 [ MD++_Get H_11 ] ; 
  set H22_0 [ MD++_Get H_22 ] ; 
  set H33_0 [ MD++_Get H_33 ] ;
  set getcompressedconfig 1
  set strain 0.07

  if { $getcompressedconfig == 1 } { 
    set maxitereps 200
    set maxiter    100
    set C11 520000  
    set C44 160000
    set factor 0.7
    set H12_fix [ expr 1.0*$H22_0*$strain ]
    MD++ H_12 = $H12_fix
    for { set iter 0 } { $iter <= $maxiter } { incr iter 1 } {

      set sig_xx [ MD++_Get TSTRESSinMPa_xx ] ;
      set sig_yy [ MD++_Get TSTRESSinMPa_yy ] ;
      set sig_zz [ MD++_Get TSTRESSinMPa_zz ] ;
      set sig_xy [ MD++_Get TSTRESSinMPa_xy ] ;
      set sig_xz [ MD++_Get TSTRESSinMPa_xz ] ;
      set sig_yz [ MD++_Get TSTRESSinMPa_yz ] ;
      set e_xx [ expr $sig_xx / $C11 ] ;
      set e_yy [ expr $sig_yy / $C11 ] ; 
      set e_zz [ expr $sig_zz / $C11 ] ;
      set e_xy [ expr $sig_xy / $C44 ] ;
      set e_xz [ expr $sig_xz / $C44 ] ; 
      set e_yz [ expr $sig_yz / $C44 ] ;

      set H11_cur [ MD++_Get H_11 ] ; 
      set H22_cur [ MD++_Get H_22 ] ;
      set H33_cur [ MD++_Get H_33 ] ;

      set H11_new [ expr ${H11_cur}*(1.0+$e_xx*$factor) ] ; 
      MD++ H_11 = ${H11_new}
      set H22_new [ expr ${H22_cur}*(1.0+$e_yy*$factor) ] ; 
      MD++ H_22 = ${H22_new}
      set H33_new [ expr ${H33_cur}*(1.0+$e_zz*$factor) ] ;    
      MD++ H_33 = ${H33_new}
      set H12_new [ MD++_Get H_12]

      if { $iter == [expr $maxiter + 1] } {
            MD++ conj_ftol = 2e-6 conj_fixbox = 1 relax
      }
      MD++ eval
   }
   relax_fixbox
   MD++ finalcnfile = "0K_${strain}_relaxed.cn" writecn
   MD++ finalcnfile = "0K_${strain}_relaxed.lammps" writeLAMMPS
  }


  set maxitereps 10
  set maxiters 10
  set maxloaditers 30
  set index 0
  set strain0 $strain
  for { set itereps 0 } { $itereps <= $maxitereps } { incr itereps 1 } {
    set x0 0
    set y0 [ expr -0.4 + $itereps * (0.5-(-0.4))/$maxitereps ] 
    puts "y0 = $y0"
    set z0 0
    for { set iter 0 } { $iter <= $maxiters } {incr iter 1 } { 
      set epsilon [ expr 0.05 + $iter * (0.4-0.1)/$maxiters ]
      MD++ incnfile = "0K_0.0_relaxed.cn" readcn
      make_ellipse_dislocation_loop $x0 $y0 $z0 $epsilon
      MD++ incnfile = "0K_${strain0}_relaxed.cn" readcn
      MD++ commit_storedr
      MD++ finalcnfile = "init_cu_${status}_${index}.lammps" writeLAMMPS
      set strain $strain0
      for { set loaditer 0 } { $loaditer <= $maxloaditers } { incr loaditer 1 } { 
        set localmaxiter 100
        set C11 520000  
        set C44 160000
        set factor 0.7
        set H12_fix [ expr 1.0*$H22_0*$strain ]
        MD++ H_12 = $H12_fix
        for { set localiter 0 } { $localiter <= $localmaxiter } { incr localiter 1 } {

          set sig_xx [ MD++_Get TSTRESSinMPa_xx ] ;
          set sig_yy [ MD++_Get TSTRESSinMPa_yy ] ;
          set sig_zz [ MD++_Get TSTRESSinMPa_zz ] ;
          set sig_xy [ MD++_Get TSTRESSinMPa_xy ] ;
          set sig_xz [ MD++_Get TSTRESSinMPa_xz ] ;
          set sig_yz [ MD++_Get TSTRESSinMPa_yz ] ;
          set e_xx [ expr $sig_xx / $C11 ] ;
          set e_yy [ expr $sig_yy / $C11 ] ; 
          set e_zz [ expr $sig_zz / $C11 ] ;
          set e_xy [ expr $sig_xy / $C44 ] ;
          set e_xz [ expr $sig_xz / $C44 ] ; 
          set e_yz [ expr $sig_yz / $C44 ] ;

          set H11_cur [ MD++_Get H_11 ] ; 
          set H22_cur [ MD++_Get H_22 ] ;
          set H33_cur [ MD++_Get H_33 ] ;

          set H11_new [ expr ${H11_cur}*(1.0+$e_xx*$factor) ] ; 
          MD++ H_11 = ${H11_new}
          set H22_new [ expr ${H22_cur}*(1.0+$e_yy*$factor) ] ; 
          MD++ H_22 = ${H22_new}
          set H33_new [ expr ${H33_cur}*(1.0+$e_zz*$factor) ] ;    
          MD++ H_33 = ${H33_new}
          set H12_new [ MD++_Get H_12]

          if { $iter == [expr $maxiter + 1] } {
                MD++ conj_ftol = 2e-6 conj_fixbox = 1 relax
          }
          MD++ eval
        }
        MD++ conj_ftol = 2e-6 
        MD++ conj_fevalmax = 10
        MD++ conj_fixbox = 1
        MD++ relax
        set strain [ expr $strain+0.001 ]
        set strain [format "%.4f" $strain]
        MD++ finalcnfile = "cu_${status}_${index}.cn" writecn
        MD++ finalcnfile = "cu_${status}_${index}.lammps" writeLAMMPS
        set index [expr $index + 1]
      }
    }
  }
  exitmd

} elseif { $status == 20 } {
  # visualization (from stringrelax_parallel runs)
  MD++ setnolog
  initmd "view-$n" 
#  readpot
  setup_window
  openwindow
  exitmd

} else {
        
 puts "unknown status = $status"
 exitmd 

} 


