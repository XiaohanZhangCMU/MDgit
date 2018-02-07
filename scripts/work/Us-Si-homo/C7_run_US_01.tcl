# -*-shell-script-*-
# for different setup, change following
# 1) read_potential
# 2) alloc_QLM, UMB_order_param, Rc and cutoff
# 3) read_initcn 
# 4) definitely, the directory name

set graphic "Off"

if { $argc==0} {
	set n 0
	set ncpu 0
        set P 110
        set T 300 
} elseif { $argc > 0 } {
	set n [ lindex $argv 0 ]
	set ncpu [ lindex $argv 1 ]
        set P [ lindex $argv 2 ]
	set T [ lindex $argv 3 ]
}

if { $graphic!="On" } {
    exec sleep $n
}


proc read_potential { } {
    MD++ {
    potfile=~/Codes/MD++/potentials/EAMDATA/eamdata.CuMishin
    eamgrid = 5000 readeam NNM = 600 }
}

proc init_MC { n P T } {
#MD++ setnolog
MD++ setoverwrite
MD++ dirname = runs/ho-cu-DN/${T}K_${P}eps_112/UMBnobias_${n}
MD++ srand48bytime
MD++ srandbytime
read_potential
}

proc setup_UMB { } {
MD++ YES_UMB = 1
MD++ UMB_K = 60.0
MD++ SKIN = 1.1
MD++ MC_UMB_cal_period=2500
MD++ MC_UMB_log_period=1000000
MD++ UMB_equilstep = 20000000
MD++ savecontinuecnfreq = 10000000
}

proc alloc_QLM { } {
MD++ UMB_order_param = 41
MD++ Rc_for_QLM = 3.09
MD++ QLM_cutoff = 0.35
#MD++ N_lgst_index = -2
MD++ allocUMBorder
}

proc openwindow { } {
MD++ {
#Plot Configuration
atomradius = 0.9 #bondradius = 0.3 bondlength = 2.8285 #for Si
atomcolor = orange highlightcolor = purple  bondcolor = red backgroundcolor = gray
fixatomcolor = yellow
color00 = "orange"  color01 = "purple" color02 = "green"
color03 = "magenta" color04 = "cyan"   color05 = "red"
color06 = "gray80"  color07 = "white"
plot_color_windows = [ 2
                       -10   -1.99  0
#                       -1.0 -0.01   5
#                      -0.01 0.01   4
#                       0.01 1.0    3
                       0.33 10.0   2
                     ]
plot_atom_info = 2
plot_map_pbc = 1
plotfreq = 10000
rotateangles = [ 0 0 0 1.1 ]
openwin alloccolors rotate saverot refreshnnlist eval plot
}
}

proc run_setup { T } {
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
timestep = 0.2   
totalsteps = 10000000001
}
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
	set filename total${k}.txt
	set check [ glob -nocomplain ../NEB/$filename ]
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


proc read_initcn { n T P } {
MD++ incnfile = "../pure_shear_stress.cn" readcn
set H11 [ MD++_Get H_11 ]
set H12 [ MD++_Get H_12 ]
set H13 [ MD++_Get H_13 ]
set H21 [ MD++_Get H_21 ]
set H22 [ MD++_Get H_22 ]
set H23 [ MD++_Get H_23 ]
set H31 [ MD++_Get H_31 ]
set H32 [ MD++_Get H_32 ]
set H33 [ MD++_Get H_33 ]
MD++ setconfig2

set check [ glob -nocomplain "continue.cn" ]
set exist [ llength $check ]
if { $exist < 1 } {
  MD++ UMB_continue = 0
} else {
  MD++ UMB_continue = 1
}

set UMBcontinue [ MD++_Get UMB_continue ]
if { $UMBcontinue > 0 } {
MD++ readcontinuecn
} else {
MD++ incnfile = "../NEB/run_${n}/FFS_${n}_0001.cn" readcn
MD++ H_11 = $H11 H_12 = $H12 H_13 = $H13 
MD++ H_21 = $H21 H_22 = $H22 H_23 = $H23
MD++ H_31 = $H31 H_32 = $H32 H_33 = $H33
#MD++ finalcnfile = "../run_${n}/FFS_${n}_0001.cn" writecn readcn
}

}



#*******************************************
# Main program starts here

#*******************************************
if { $graphic=="On" } {
  MD++ setnolog 
  MD++ printfreq = 10000 
} else {
  MD++ printfreq = 1000000000000
} 

init_MC $n $P $T

# distribute windows run after original stuff
#if { $n == 0 } {
#     distribute $ncpu 
#}

# assign directory number, according to window
set direc_no [setup_window $n]

if { 0 } {
if { ${direc_no} >= 4 } {
set direc_no 4
}
}

# setup UMB parameters
setup_UMB

# print out n_center, width, UMB_K for the current job
set temp1 [ MD++_Get n_center ]
set temp2 [ MD++_Get delta_n ]
set umbk [ MD++_Get UMB_K ] 
puts "$temp1 $temp2 $umbk"
 
# read reference cn and initcn
read_initcn ${direc_no} $T $P

if { $n <= 1 } {
    set direc_no 0
	if { 1 } {
	    MD++ UMB_K = -1
	    MD++ MC_UMB_cal_period = 100000
	}
}

# you have to alloc_QLM after read_cn !  ### never change order!
alloc_QLM  

# openwindow if graphic is On
if { $graphic=="On" } {
MD++ N_lgst_index = -2
openwindow
}

#MD++ sleep

# run UMB MC
run_setup $T
MD++ runMC

