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
proc initmd { status {T 0} {epVAL 0} {opt 0} {alt 0}  } {
#MD++ setnolog
MD++ setoverwrite
#relative to the  directory of lauching 
#MD++ dirname = runs/test-stress-strain/stress-strain-$n
if { $status == 0 } {
   MD++ dirname = runs/ho-sige-DN/0K_stressstrain
}  elseif { $status == 1 } {
   MD++ dirname = runs/ho-sige-DN/${T}K_${epVAL}eps_112
} elseif { $status == 2 || $status == 22 } {
   MD++ dirname = runs/ho-sige-DN/md_${T}K_${epVAL}eps_112-clean-$opt-$alt
} elseif { $status == 3 } {
   MD++ dirname = runs/ho-sige-DN/mc_${T}K_${epVAL}eps_112
} elseif { $status == 44 } {
   MD++ dirname = runs/ho-sige-DN/${T}K_${epVAL}eps_112/NEB
} elseif { $status == 4 } {
   MD++ dirname = runs/ho-sige-DN/${T}K_${epVAL}eps_112/NEB
} elseif { $status == 5 || $status == 21 } { 
   MD++ dirname =  runs/ho-sige-DN/${T}K_${epVAL}eps_112/UMB_$opt
} elseif { $status == 10 } { 
   MD++ dirname =  runs/ho-sige-DN/view-${T}K_${epVAL}eps_112
}
exec cp -r /home/xzhang11/Planet/Libs/MD++.svn/scripts/work/SiGeHomo-sworig/ss_sige.tcl .  
set myname [ MD++_Get "myname" ]
readmeam-according-to-myname $myname
#max number neigbors
MD++ NNM = 300

}

#------------------------------------------------------------
proc readmeam-lammps { } { 
#Read in MEAM potential (Baskes format)
MD++ meamfile = "~/Planet/Libs/MD++UMB.svn3/potentials/MEAMDATA/meamf"
MD++ meafile = "~/Planet/Libs/MD++UMB.svn3/potentials/MEAMDATA/SiGe.meam" 
MD++ nspecies = 2 element0 = "Siz" element1 = "Ge" 
MD++ {
rcut = 4.5  readMEAM
NNM = 300
atommass = [28.0855 72.64]
} }

# make sure the coordinate is right hand sided.

proc make_perfect_crystal { nx ny nz } {
    MD++ crystalstructure = diamond-cubic latticeconst = 5.4309529817532409 #(A) for Si
    MD++ latticesize = \[  1 -2 1  $nx  1 1 1  $ny  1 0 -1  $nz \]
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
plot_color_windows = [ 5
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
plot_color_windows = [ 1
                       0.6 9.4   1  
                       9.4  10   5
                       10 20   6
                       20 50   8
                       0  0.6  4
                     ]

 plot_limits = [ 1 -10 10 -0.2 0.2 -10 10 ]
# plot_limits = [ 1 -10 10 -10 10 -10 -0.10 ]
# plot_limits = [ 1 -10 10 -10 10 -0.28 -0.22 ]

#
#xiaohan

#plot_color_windows = 5 plot_limits = [ 1 -10 10 -0.001 0.028 -10 10 ]
#plot_color_windows = 0 
plot_atom_info = 1 # reduced coordinates of atoms
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


proc setup_umb_window { } {
MD++ {
#Plot Configuration
atomradius = 0.9 #bondradius = 0.3 bondlength = 2.8285 #for Si
atomcolor = orange highlightcolor = purple  bondcolor = red backgroundcolor = gray
fixatomcolor = yellow
color00 = "orange"  color01 = "purple" color02 = "green"
color03 = "magenta" color04 = "cyan"   color05 = "red"
color06 = "gray80"  color07 = "white"
plot_color_axis = 9 
#plot_limits = [ 1 -10 10 -0.05 0.05 -10 10 ]
NCS = 4
plot_atom_info = 2
plot_color_windows = [ 1
#                      1.2 9.4   1
                      0.9 9.4  1  
                       9.4  10   5
                       10 20   6
                       20 50   8
                       0  0.6  4
                     ]
plot_map_pbc = 1
plotfreq = 10
rotateangles = [ 0 0 0 1.0 ]
}
}


proc openwindow { } { 
#setup_window
MD++ openwin alloccolors rotate saverot eval plot
}

#--------------------------------------------
proc exitmd { } { MD++ quit }
#end of proc exitmd
#--------------------------------------------


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
proc setup_md { } { MD++ {     
T_OBJ = 300 #Kelvin #add by xiaohan

equilsteps = 0  totalsteps = 5000 timestep = 0.0001 # (ps)
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
}

MD++ NNM = 300
MD++ atommass = \[28.0855 72.64\]
}
#end of proc setup_md

proc setup_mc { T } {
MD++ {
output_fmt="curstep EPOT MC_accept_ratio N_lgst_temp MC_UMB_accept_ratio N_lgst_cluster N_solid_P"
saveprop = 1 openpropfile 
savepropfreq = 100000 
#printfreq = 1000000
}
MD++ T_OBJ = $T  
MD++ {
#mcatom = -1 timestep = 0.0085 totalsteps = 40001 #Angstrom
mcatom = 1  
timestep = 0.15   
totalsteps = 10000000001
}
}

proc setup_UMB { UMB } {
MD++ YES_UMB = 1
#MD++ UMB_K_in_Kelvin = 90.0
MD++ UMB_K_in_Kelvin = $UMB

MD++ MC_UMB_cal_period=2500
MD++ MC_UMB_log_period=1000000
MD++ UMB_equilstep = 60000000
MD++ savecontinuecnfreq = 1000000
}

proc alloc_QLM { } {
MD++ UMB_order_param = 45
MD++ Rc_for_QLM = 4.0
#glide
#MD++ QLM_cutoff = 1.90
#shuffle

#MD++ QLM_cutoff = 0.35
MD++ QLM_cutoff = 1.8

#MD++ N_lgst_index = -2

MD++ RLIST = 4.7
#               RUMB_R0  UMB_nvec   UMB_thick    UMB_slipvec    UMB_Hmin     UMB_Hmax   UMB_NNBR
#glide
#MD++ input = \[ 0 0 0   0 -1 0  2.352  1 0 0   0.16667 0.5 3 \]
#shuffle
MD++ input = \[ 0 0.7802 0    0 -1 0     2.352        0 0 1     0.83333       1.166667  1 \]

MD++ allocUMBorder
}

proc read_lam { } {
   set fp [ open "../NEB/lam_array.dat" r ]
   set data [ read $fp ]
   close $fp
   return $data 
}

proc get_N { } {
   set Nmax 0
   for { set k 0 } { $k <= 200 } { incr k } {
	set filename ../NEB/total${k}.txt
	set check [ glob -nocomplain $filename ]
	set NUM [ llength $check ]
	if { $NUM > 0 } {
	   set Nmax $k 
        }
   }
#   puts $Nmax
   return $Nmax
}

proc distribute { ncpu } {
   set lambda [ read_lam ]
#   puts $lambda
   set Nmax [ get_N ]
   set index 0
   set lam_real 0
   set N_real 0
   set Nmax [ expr $Nmax + 2 ]
   set lam_max [ lindex $lambda $Nmax ]

   puts "lam_max = $lam_max"
   puts "Nmax = $Nmax"

   set delta_n [ expr 1.0*${lam_max}/$ncpu/2.0 ]
   set delta_n_int [ expr int(3.8*$delta_n) ]

   for { set i 0 } {$i < $ncpu } { incr i 1 } {	
        set current 10000
        if { $i == 0 } { 
            set coe 1.0 
        } else {
            set coe 2.0 
        }

	set lam_real [ expr $lam_real+$coe*$delta_n ]
#        puts "lam_real=$lam_real" 
	for { set j 2 } { $j <= $Nmax } { incr j 1 } {
           set lam_cur_n [ lindex $lambda $j ] 
	   set current_n [ expr {abs($lam_cur_n-$lam_real)} ]
           puts "i=$i ;j=$j ;lam_cur_n=$lam_cur_n ; current_n = $current_n ; Nmax = $Nmax ; lam_real = $lam_real "
           if { $current_n > $current } {
           break;
           } 
           set lam_cur $lam_cur_n
           set current $current_n
        }
        set nind [ expr $j-3 ]
        set lam_int [ expr {int($lam_real)} ]
        set fp [ open "../UMB_arg${i}.txt" w ]
#        puts ".."
#        puts $fp $i
        puts $fp $lam_int
#	set delta_n 15
        set delta_nx2 [ expr $delta_n_int*2 ]
        puts $fp $delta_nx2
        puts $fp $nind
        close $fp
   }
}

proc setup_umbrella { n } {

    set lambda [ read_lam ]
    set lam_min [ lindex $lambda 2 ]
 
    while { 1 } {
    set check [ glob -nocomplain "../UMB_arg${n}.txt" ]
    set NUM [ llength $check ]
    if { $NUM > 0 } {
	break;
    }
    exec sleep 10
    }
    
    set fp [ open "../UMB_arg${n}.txt" r ]
    set data [ split [ read $fp ]]
    set para1 [ lindex $data 0 ]
    set para2 [ lindex $data 1 ]
    set para3 [ lindex $data 2 ]
    close $fp

    set criterion [ expr 0.7*$lam_min ]
    if { $para1 < $criterion } {
#        MD++ UMB_K_in_Kelvin = -1
    }

    MD++ n_center = $para1
    MD++ delta_n = $para2

    puts "-----------------------------------------------------------"
    puts $para1
    puts $para2
    puts $para3
    set tempcenter [ MD++_Get n_center ]
    puts "n_center = $tempcenter"
    puts "-----------------------------------------------------------"
    return $para3
}


proc read_initcn { n o_T o_epVAL o_epsilon } {

MD++ incnfile="../../${o_T}K_${o_epVAL}eps_112/inter45.cn" readcn saveH
MD++ incnfile="../../NEB/disl-nuc-homo-0/0K_${o_epsilon}_relaxed_surf001.cn" readcn restoreH

relax_fixbox
MD++ finalcnfile = "final.cn" writecn
MD++ incnfile = "final.cn" readcn setconfig2
MD++ assignUMBindexref
MD++ incnfile = "../NEB/run_${n}/FFS_${n}_0001.cn" readcn restoreH

}



## data average proc
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
} elseif { $argc > 3 } {
 set opt [lindex $argv 3]
}
puts "opt = $opt"

if { $argc <= 4 } {
 set alt -1
} elseif { $argc > 4 } {
 set alt [lindex $argv 4]
}

set UMB $alt
puts "UMB = $UMB"


if { $status == 0 } {
  initmd $status
  make_perfect_crystal 6 5 10
  setup_window

  MD++ {
   conj_fixboxvec = [ 0 0 1
                    1 0 1 
                    0 0 0 ]
  }

  MD++ conj_fixbox = 1
  MD++ conj_ftol = 2e-2 conj_itmax = 1000 conj_fevalmax = 2000
  MD++ relax
  MD++ writeall = 1 finalcnfile = "0K_perfect.cn" writecn
  MD++ eval
  set H11_0 [ MD++_Get H_11 ]
  set H12_0 [ MD++_Get H_12 ]
  set H13_0 [ MD++_Get H_13 ]
  set H22_0 [ MD++_Get H_22 ]
  set H23_0 [ MD++_Get H_23 ]
  set H33_0 [ MD++_Get H_33 ]
  set H32_0 [ MD++_Get H_32 ]

  set epsilon 0.0
  set C11 170000
  set factor 0.64
  set maxiter 100
  set epsIncrement 0.002

  for { set itereps 0 } { $itereps <= 200 } { incr itereps 1 } {

      puts "\n\n epsilon=$epsilon"
#exec sleep 1

      set H32_fix [ expr 1.0*$H22_0*$epsilon ]
      MD++ H_32 = $H32_fix
      MD++ eval
      MD++ plot

# adjust other strains in small steps
      for { set iter 0 } { $iter <= $maxiter } { incr iter 1 } {

          set sig_xx [ MD++_Get TSTRESSinMPa_xx ]
	  set sig_yy [ MD++_Get TSTRESSinMPa_yy ]
	  set sig_zz [ MD++_Get TSTRESSinMPa_zz ]
	  set sig_xy [ MD++_Get TSTRESSinMPa_xy ]
	  set sig_xz [ MD++_Get TSTRESSinMPa_xz ]
	  set sig_zy [ MD++_Get TSTRESSinMPa_zy ]

	  set e_xx [ expr $sig_xx / $C11 ]
	  set e_yy [ expr $sig_yy / $C11 ]
	  set e_zz [ expr $sig_zz / $C11 ]
	  set e_zx [ expr $sig_xz / $C11 ]
	  set e_xy [ expr $sig_xy / $C11 ]

	  puts "iter = $iter"
	  puts "sig_xx = $sig_xx sig_yy = $sig_yy sig_zz = $sig_zz"
	  puts "sig_xy = $sig_xy sig_xz = $sig_xz sig_zy = $sig_zy"
	  puts "e_xx = $e_xx e_yy = $e_yy e_zz = $e_zz"

	  set H11_cur [ MD++_Get H_11 ]
	  set H22_cur [ MD++_Get H_22 ]
	  set H33_cur [ MD++_Get H_33 ]
	  set H31_cur [ MD++_Get H_31 ]
	  set H12_cur [ MD++_Get H_12 ]
	  set H32_new [ MD++_Get H_32 ]

	  set H11_new [ expr ${H11_cur}*(1.0+$e_xx*$factor) ]
	  MD++ H_11 = ${H11_new}
	  set H22_new [ expr ${H22_cur}*(1.0+$e_yy*$factor) ]
	  MD++ H_22 = ${H22_new}
	  set H33_new [ expr ${H33_cur}*(1.0+$e_zz*$factor) ]
	  MD++ H_33 = ${H33_new}
	  set H31_new [ expr ${H31_cur} +  ${H11_cur} * $e_zx*$factor ]
	  MD++ H_31 = ${H31_new}
	  set H12_new [ expr ${H12_cur} +  ${H22_cur} * $e_xy*$factor ]
	  MD++ H_12 = ${H12_new}

	  set iter_ [ expr $iter-1 ]
	  if { $iter_ == $maxiter } {
		  MD++ conj_ftol = 2e-6 conj_fixbox = 1 relax
	      }
	  MD++ eval
      }
      set SIGZY [ expr int($sig_zy) ]
      MD++ conj_ftol = 2e-6 conj_fixbox = 1 relax
      MD++ finalcnfile = "0K_${SIGZY}MPa_relaxed.cn" writecn
      MD++ finalcnfile = "NEBinit.cn" writecn
      set e_xx [ expr ($H11_new-$H11_0)/$H11_0 ]
      set e_yy [ expr ($H22_new-$H22_0)/$H22_0 ]
      set e_zz [ expr ($H33_new-$H33_0)/$H33_0 ]
      set e_zy [ expr ($H32_new)/$H22_0 ]
      set e_zx [ expr ($H31_new)/$H11_0 ]

      set fp [ open "strain.dat" a+ ]
      puts $fp "$e_xx $e_zy $e_yy $e_zz"
      close $fp
      set fp [ open "strain_info.dat" w ]
      puts $fp "e_xx e_zy e_yy e_zz"
      close $fp

      set fp [ open "stress.dat" a+ ]
      puts $fp "$sig_xx $sig_yy $sig_zz $sig_xy $sig_xz $sig_zy"
      close $fp
      set fp [ open "stress_info.dat" w ]
      puts $fp "sig_xx sig_yy sig_zz sig_xy sig_xz sig_zy"
      close $fp
      
      set fp [ open "EPOT.dat" a+ ]
      puts $fp [ MD++_Get EPOT ]
      close $fp

      set epsilon [ expr $epsilon+$epsIncrement ]
  }

} elseif { $status == 1 } {

  set T $n
  set epsilon $flag

  set epVAL [ format %03d [expr int($epsilon*1000) ] ]
  
  initmd $status $T $epVAL 
 
  MD++ srand48bytime
  MD++ initvelocity_type="Gaussian"

  make_perfect_crystal 6 5 10
  setup_md

  setup_window

  set H11_0 [ MD++_Get H_11 ]
  set H12_0 [ MD++_Get H_12 ]
  set H13_0 [ MD++_Get H_13 ]
  set H22_0 [ MD++_Get H_22 ]
  set H23_0 [ MD++_Get H_23 ]
  set H33_0 [ MD++_Get H_33 ]
  set H32_0 [ MD++_Get H_32 ]

  set H32_fix [ expr $epsilon*$H22_0 ]
  MD++_Set H_32 $H32_fix

  MD++ eval 

  relax_freebox

  MD++ fixboxvec = \[ 0 1 1 1 0 1 1 1 0 \]
  #MD++ relax
  MD++ eval saveH

  MD++ T_OBJ = $T DOUBLE_T =1 initvelocity

  MD++ totalsteps=27500
  MD++ openpropfile run
  MD++ closepropfile

  set sigxx [datafile_process "prop.out" 4 0.5 1 "AVE" ]
  set sigyy [datafile_process "prop.out" 5 0.5 1 "AVE" ]
  set sigzz [datafile_process "prop.out" 6 0.5 1 "AVE" ]
  set sigxy [datafile_process "prop.out" 7 0.5 1 "AVE" ]
  set sigyz [datafile_process "prop.out" 8 0.5 1 "AVE" ]
  set sigzx [datafile_process "prop.out" 9 0.5 1 "AVE" ]
  set H11 [datafile_process "prop.out" 10 0.5 1 "AVE" ]
  set H12 [datafile_process "prop.out" 11 0.5 1 "AVE" ]
  set H13 [datafile_process "prop.out" 12 0.5 1 "AVE" ]
  set H21 [datafile_process "prop.out" 13 0.5 1 "AVE" ]
  set H22 [datafile_process "prop.out" 14 0.5 1 "AVE" ]
  set H23 [datafile_process "prop.out" 15 0.5 1 "AVE" ]
  set H31 [datafile_process "prop.out" 16 0.5 1 "AVE" ]
  set H32 [datafile_process "prop.out" 17 0.5 1 "AVE" ]
  set H33 [datafile_process "prop.out" 18 0.5 1 "AVE" ]
  puts "$H11 $H12 $H13" 
  puts "$H21 $H22 $H23"
  puts "$H31 $H32 $H33"

  MD++ H_11 = $H11 H_22 = $H22 H_33 = $H33
  MD++ H_13 = $H13 H_31 = $H31 H_23 = $H23 H_32 = $H32_fix
  MD++ H_12 = $H12 H_21 = $H21
  MD++ ensemble_type="NVTC" totalsteps = 15000
  MD++ writeall = 1 
  MD++ initvelocity

### start steepest decent
  set C11 170000
  set factor 0.64
# increase three times bigger
  set maxiter 45

  set datafile "stress_trace.dat"
  set fileID [open $datafile a+ ]
  puts $fileID "% iter  XX    YY    ZZ    XY    YZ    ZX   (in MPa)"
  close $fileID

  for { set iter 0 } { $iter <= $maxiter } { incr iter 1 } {
#exec rm prop.out 
      set fp [ open "void.out" w ]
      puts $fp ""
      close $fp
      exec cp void.out prop.out
      MD++ openpropfile
      set e_xx [ expr $sigxx/$C11 ]
      set e_yy [ expr $sigyy/$C11 ]
      set e_zz [ expr $sigzz/$C11 ]
      set e_zx [ expr $sigzx / $C11 ]
      set e_xy [ expr $sigxy / $C11 ]

      puts "iter = $iter"
      puts "sigxx = $sigxx sigyy = $sigyy sigzz = $sigzz"
      puts "sigxy = $sigxy sigzx = $sigzx sigyz = $sigyz"
      puts "e_xx = $e_xx e_yy = $e_yy e_zz = $e_zz"

      set H11_cur [ MD++_Get H_11 ]
      set H22_cur [ MD++_Get H_22 ]
      set H33_cur [ MD++_Get H_33 ]
      set H31_cur [ MD++_Get H_31 ]
      set H12_cur [ MD++_Get H_12 ]
      set H32_new [ MD++_Get H_32 ]
      MD++ finalcnfile=inter${iter}.cn writecn

      set H11_new [ expr ${H11_cur}*(1.0+$e_xx*$factor) ]
      MD++ H_11 = ${H11_new}
      set H22_new [ expr ${H22_cur}*(1.0+$e_yy*$factor) ]
      MD++ H_22 = ${H22_new}
      set H33_new [ expr ${H33_cur}*(1.0+$e_zz*$factor) ]
      MD++ H_33 = ${H33_new}
      set H31_new [ expr ${H31_cur} +  ${H11_cur} * $e_zx*$factor ]
      MD++ H_31 = ${H31_new}
      set H12_new [ expr ${H12_cur} +  ${H22_cur} * $e_xy*$factor ]
      MD++ H_12 = ${H12_new}

      MD++ DOUBLE_T = 0 initvelocity run

      set sigxx [datafile_process "prop.out" 4 0.2 1 "AVE" ]
      set sigyy [datafile_process "prop.out" 5 0.2 1 "AVE" ]
      set sigzz [datafile_process "prop.out" 6 0.2 1 "AVE" ]
      set sigxy [datafile_process "prop.out" 7 0.2 1 "AVE" ]
      set sigyz [datafile_process "prop.out" 8 0.2 1 "AVE" ]
      set sigzx [datafile_process "prop.out" 9 0.2 1 "AVE" ]

      set datafile "stress_trace.dat"
      set fileID [open $datafile a+ ]
#puts $fileID "$iter"
#puts $fileID "XX=$sigxx YY=$sigyy ZZ=$sigzz"
#puts $fileID "XY=$sigxy YZ=$sigyz ZX=$sigzx"
      puts $fileID "$iter $sigxx $sigyy $sigzz $sigxy $sigyz $sigzx"
      close $fileID
      MD++ closepropfile      
  }

  set SIGXY [ expr int($sigxy) ]
  MD++ conj_fixbox=1 conj_ftol = 1e-12
  MD++ conj_itmax = 9000 conj_fevalmax = 9000
#MD++ relax relax
  MD++ writeall = 1
  MD++ finalcnfile = "pure_shear_stress.cn" writecn
  set e_xx [ expr ($H11_new-$H11_0)/$H11_0 ]
  set e_yy [ expr ($H22_new-$H22_0)/$H22_0 ]
  set e_zz [ expr ($H33_new-$H33_0)/$H33_0 ]
  set e_zy [ expr ($H32_new)/$H22_0 ]

  set fp [ open "strain.dat" w ]
  puts $fp "$e_xx $e_zy $e_yy $e_zz"
  close $fp

  set fp [ open "stress.dat" w ]
  puts $fp "$sigxx $sigyy $sigzz $sigxy $sigzx $sigyz"
  close $fp

  set fp [ open "strain_info.dat" w ]
  puts $fp "e_xx e_zy e_yy e_zz"
  close $fp

  set fp [ open "stress_info.dat" w ]
  puts $fp "sig_xx sig_yy sig_zz sig_xy sig_zx sig_yz"
  close $fp

#  MD++ sleep
  exitmd
  
} else {
        
 puts "unknown status = $status"
 exitmd 

} 

