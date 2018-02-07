
set graphic "Off"

if { $argc==0} {
        set eps 110
        set T 0
} elseif { $argc > 0 } {
        set eps [ lindex $argv 0 ]
        set T [ lindex $argv 1 ]
}

puts "eps = $eps"

proc read_potential { } {
    MD++ {
    potfile=~/Codes/MD++/potentials/EAMDATA/eamdata.CuMishin
    eamgrid = 5000 readeam NNM = 600 }
}

proc init_MC { eps T } {
#MD++ setnolog
MD++ setoverwrite
MD++ dirname = runs/ho-cu-DN/${T}K_${eps}eps_112/NEB
MD++ srand48bytime
MD++ srandbytime
read_potential
}

proc alloc_QLM { } {
#MD++ YES_UMB = 1
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
rotateangles = [ 0 0 0 1.0 ]
openwin alloccolors rotate saverot refreshnnlist eval plot
}
}

init_MC ${eps} $T
MD++ incnfile=NEBinit.cn readcn setconfig2
alloc_QLM

MD++ incnfile=NEBinit.cn readcn
openwindow

set EPOT0 [MD++_Get EPOT ]

set fp [ open "lam_array.dat" w ]
puts $fp "\["
puts $fp "200"

set index 0
set index1 0

for { set iter 0 } { $iter <=25 } { incr iter 1 } {

MD++ incnfile="chain_no_${iter}.cn" readcn 
alloc_QLM
openwindow

set Nc [MD++_Get N_lgst_cluster ]

if { $Nc > 0 } {
set index [expr $index + 1 ]
set direc_no [ expr $index - 1]
exec mkdir run_${direc_no}
MD++ finalcnfile="run_${direc_no}/FFS_${direc_no}_0001.cn" writecn
puts $fp "$Nc"
set fp1 [ open "total${direc_no}.txt" w ]
close $fp1
}

set EPOT [ MD++_Get EPOT ]
set dE [ expr $EPOT - $EPOT0 ]
puts "iter=$iter Nc = $Nc dE = $dE index = $index"

if { $iter > 0 } {
if { $dE < $dE1 } {
set index1 [ expr $index1 + 1 ]
}
if { $index1 > 2 } {
break
}
}

set dE1 $dE

}

for { set iter 0 } { $iter <=200 } { incr iter 1 } {
set Nc [ expr $Nc + 10 ]
puts $fp "$Nc"
}
puts $fp "\]"
close $fp

MD++ sleep
