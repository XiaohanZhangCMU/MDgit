# -*-shell-script-*-
# Monte Carlo simulation of Stinger-Weber Silicon

# if MODE == 1 change cell size
set MODE 1

set graphic "Off"
source "scripts/Examples/Tcl/startup.tcl"

# default = 0 
# change to 1 only for continueing run
# when we need a week to equilibrate the UMB histogram.
#MD++ UMB_continue = 1
#MD++ DLNdiv = 9 DLxyz = 1
#MD++ USE_Total_Number = 1
MD++ savecontinuecnfreq = 10000000

if { $argc==0} {
	set n 20
	set nc 100
        set P 110 
        set T 300 
} elseif { $argc > 0 } {
	set n [ lindex $argv 0 ]
	set nc [ lindex $argv 1 ]
        set P [ lindex $argv 2 ]
	set T [ lindex $argv 3 ]
}

#################

#####333
proc init_MC { P T } {
#setnolog
MD++ setoverwrite
#set PVAL [ format %03d $P ]
MD++ dirname = runs/ho-cu-DN/${T}K_${P}eps_112/Kinetic
MD++ srand48bytime
MD++ srandbytime
#------------------------------------------------------------
#Create Perfect Lattice Configuration
#
MD++ {
potfile=~/Planet/Libs/MD++.svn/potentials/EAMDATA/eamdata.CuMishin
eamgrid = 5000 readeam NNM = 600
}
}

proc setup_UMB { T nc } {

MD++ YES_UMB = 1
MD++ UMB_K = 60.0
MD++ SKIN = 1.1
MD++ n_center = $nc
MD++ delta_n = 10
MD++ MC_UMB_cal_period=2500
MD++ MC_UMB_log_period=1000000
MD++ UMB_equilstep = 100000
#MD++ savecontinuecnfreq = 10000000

}

proc alloc_QLM { } {
MD++ UMB_order_param = 41
MD++ Rc_for_QLM = 3.09
MD++ QLM_cutoff = 0.35
#MD++ N_lgst_index = -2
#MD++ HETERO_DN = 1
MD++ allocUMBorder
}


proc openwindow { } {
MD++ {
#Plot Configuration
#
atomradius = 0.9 #bondradius = 0.3 bondlength = 2.8285 #for Si
atomcolor = orange highlightcolor = purple  bondcolor = red backgroundcolor = gray
fixatomcolor = yellow
color00 = "orange"  color01 = "purple" color02 = "green"
color03 = "magenta" color04 = "cyan"   color05 = "purple"
color06 = "gray80"  color07 = "white"
#plot_color_axis = 9
plot_color_windows = [ 2
                          -2.01 -1.99 0  #color00 = orange
                          0.33 5.0  2  #color01 = purple
                     ]
plot_atom_info = 2
plotfreq = 2500
rotateangles = [ 0 0 0 1.2 ]
openwin alloccolors rotate saverot refreshnnlist eval plot
}
}



proc run_setup { T } {
MD++ {
output_fmt="curstep EPOT MC_accept_ratio N_lgst_temp MC_UMB_accept_ratio N_lgst_cluster"
#saveprop = 1 openpropfile 
#saveprop = 1
savepropfreq = 100000
plotfreq = 2500
#printfreq = 1000000
}
MD++ T_OBJ = $T  
MD++ {
#mcatom = -1 timestep = 0.0085 totalsteps = 40001 #Angstrom
mcatom = 1  
timestep = 0.2   
totalsteps = 1000000000
}
}

proc run_setup_MD { T } {
MD++ {
#	plot_color_axis = 9
	srand48bytime
	atommass = 63.556
#        totalsteps = 1000000
	timestep = 0.0005
	saveprop = 1 openpropfile
	savepropfreq = 20 plotfreq = 20
	ensemble_type="NVE" 
	integrator_type="VVerlet"
        initvelocity_type="Gaussian"
        DOUBLE_T = 0
        initvelocity
   	output_fmt="curstep EPOT KATOM Tinst TSTRESS_zz N_lgst_cluster"
}
}

proc read_lam { } {
   set fp [ open "../lam_array.dat" r ]
   set data [ read $fp ]
   close $fp
   return $data 
}

proc get_N { } {
   set Nmax 0
   for { set k 0 } { $k <= 200 } { incr k } {
	set filename total${k}.txt
	set check [ glob -nocomplain ../$filename ]
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
	set lam_max [ lindex $lambda $Nmax ]
#        puts $lam_max
 	set delta_n [ expr 1.0*${lam_max}/$ncpu/2.0 ]
        set delta_n_int [ expr int(1.4*$delta_n) ]
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
#           puts "$j $lam_cur_n $current_n"
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

proc setup_window { n } {

    set lambda [ read_lam ]
    set lam_min [ lindex $lambda 2 ]
 
    while { 1 } {
    set check [ glob -nocomplain ../UMB_arg${n}.txt ]
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
#        MD++ UMB_K = -1
    }

    MD++ n_center = $para1
    MD++ delta_n = $para2
    puts $para3
    return $para3
}


proc read_initcn2 { n } {

MD++ incnfile = "../NEB/NEBinit.cn" readcn
MD++ setconfig2

MD++ incnfile = "../UMB_${n}/continue.cn" readcn

}



#*******************************************
# Main program starts here

#*******************************************

  if { $graphic=="On" } {
  MD++ setnolog 
  MD++ printfreq = 1000
  } else {
  MD++ printfreq = 1000000000000
  } 

  init_MC $P $T
puts " i am here -2"
  setup_UMB $T $nc

  set temp1 [ MD++_Get n_center ]
  set temp2 [ MD++_Get delta_n ]
  set umbk [ MD++_Get UMB_K ] 
  puts "$temp1 $temp2 $umbk"

  read_initcn2 $n 

  puts "i am here 1"

  alloc_QLM  
 
  MD++ Kinetic_Time = 300
  MD++ Kinetic = 1
  MD++ MCequilstep = 1000000
  # must do this after read initcn

puts "i am here 2"
  if { $graphic=="On" } {
  MD++ N_lgst_index = -2
  openwindow
  }

puts "i am here 3"
for { set i 1 } { $i <= 1000000 } { incr i 1 } {
  run_setup $T
  MD++ saveprop = 0
  MD++ runMC
 
  MD++ conj_step = 0
  run_setup_MD $T
  MD++ run
}

