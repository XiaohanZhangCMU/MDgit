# -*-shell-script-*-
# Ising 3D Monte Carlo simulation
#-------------------------------------------

set graphic "Off" 
set lam_auto "On"
source "scripts/Examples/Tcl/startup.tcl"

if { $argc==0} {
	set n 0
	set ncpu 0
        set h 0.55
        set t 2.71111 
} elseif { $argc > 0 } {
	set n [ lindex $argv 0 ] 
        set ncenter [ lindex $argv 1 ] 
	set deltan [ lindex $argv 2 ]
        set h [ lindex $argv 3 ] 
	set t [ lindex $argv 4 ]
}

set check [ glob -nocomplain runs/2D_ISING ]
set key [ llength $check ]
if { $key < 1 } {
exec mkdir runs/2D_ISING
}

set check [ glob -nocomplain runs/2D_ISING/FFS_${h}_${t} ]
set key [ llength $check ]
if { $key < 1 } {
exec mkdir runs/2D_ISING/FFS_${h}_${t}
}


#*******************************************
# Definition of procedures
#*******************************************
proc init_ising { n h t } {
MD++ setnolog
MD++ setoverwrite
MD++ dirname = "runs/2D_ISING/FFS_${h}_${t}/UMB_${n}"
MD++ srand48bytime 
MD++ srandbytime
#--------------------------------------------
# Initializtion
#
MD++ NX = 100  NY = 100 NZ = 1 input = -1  initspin
#incnfile = ../ising2d/s100.cn readcn
#--------------------------------------------
# Physical parameters
MD++ J = 1     # interaction strength
MD++ H = ${h}  # applied field
MD++ kBT = ${t} # temperature
set kumb [ expr 0.16*$t ]
MD++ UMB_K = ${kumb}
MD++ zipfiles = 1
} 

#--------------------------------------------
proc openwindow { } { MD++ {
# Plot Configuration
#
atomradius = 1.0 win_width = 480 win_height = 480 rotateangles = [ 0 0 0 1.7]
color00 = white color01 = blue backgroundcolor = gray80
openwin alloccolors rotate saverot plot
} }
#end of proc openwindow

proc run_setup { } {
MD++ {
	saveprop = 1 savepropfreq=1 plotfreq = 100
        savecn   = 0 calcryfreq=1 calcryfreq1 = 250
	openpropfile
}
}	


proc read_lam { } {
   set fp [ open "../lam_array.txt" r ]
   set data [ read $fp ]
   close $fp
   return $data 
}


proc setup_window1 { ncenter } {

    set lambda [ read_lam ]
    set lam_min [ lindex $lambda 2 ]

    set criterion [ expr 0.7*$lam_min ]
    if { $ncenter < $criterion } {
        MD++ UMB_K = -1
    }

    for { set j 2 } { $j <= 10 } { incr j 1 } {
	set lam_now [ lindex $lambda $j ]
        set cur_diff [ expr {abs($lam_now-$ncenter)}]
        if { $j > 2 } {
	    if { $cur_diff > $pre_diff } {
                break;
            }
        }
        set pre_diff $cur_diff
    }
   
    set nindex [ expr $j-3 ]
    return $nindex
}

proc get_init_cn { m } {
	set mm [ expr $m - 1 ]
	set fp [ open ../total${mm}.txt r ]
	set data [ split [ read $fp ] ]
	close $fp
	
	set fp [ open ../weight${mm}.txt r ]
	set weights [ split [ read $fp] \n]
	close $fp
	set NUM [ llength $weights ]
	set NUM [ expr $NUM-1 ]

	set ran_max [ lindex $data 2 ]
	set n_rand [ expr { rand()*$ran_max } ]


	set sum 0
	for { set i 0 } { $i < $NUM } { incr i 1 } {
		set weight [ lindex $weights $i ]
		set sum [ expr $sum+$weight ]
		if { $sum >= $n_rand } {
			set cn_no [ expr $i+1 ]
#			puts "yes"
			break;
		}
#		puts "$i $weight $sum"
	} 
#	puts $i
#	puts $sum
#	puts ${cn_no}
	set cn_no [ expr {int($cn_no)} ]
	return $cn_no
}



#*******************************************
# Main program starts here

#*******************************************

  # Monte Carlo simulation (brute force)

 
  if { $graphic=="On" } {
  MD++ setnolog 
  MD++ printfreq = 100 
  } else {
  MD++ printfreq = 1000000000000
  } 

  init_ising $n $h $t

  set direc_no [setup_window1 $ncenter]

  MD++ YES_UMB = 1
  MD++ n_center =$ncenter
  MD++ delta_n =$deltan
  set umbk [ MD++_Get UMB_K ] 
  puts "n_center=$ncenter delta_n=$deltan US_K=$umbk"

  MD++ incnfile = "../run_${direc_no}/FFS_${direc_no}_0001.cn"
  MD++ readcn

  if { $graphic=="On" } {
  openwindow
  }

  MD++ allocDFS 
  run_setup   
  MD++ totalsteps = 20000 MCrun

  MD++ quit

