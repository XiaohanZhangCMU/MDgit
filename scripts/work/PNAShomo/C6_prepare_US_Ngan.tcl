
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
    potfile=~/Planet/Libs/MD++.svn/potentials/EAMDATA/eamdata.CuMishin
    eamgrid = 5000 readeam NNM = 600 }
}

proc init_MC { eps T } {
#MD++ setnolog
MD++ setoverwrite
MD++ dirname = runs/ho-cu-DN/${T}K_${eps}eps_112/NEB_old
MD++ srand48bytime
MD++ srandbytime
read_potential
}

proc alloc_QLM { T  } {
#MD++ YES_UMB = 1
MD++ UMB_order_param = 10
MD++ Rc_for_QLM = 3.4
if { $T < 450 } {
MD++ QLM_cutoff = 0.33
} 
if { $T < 550 && $T >= 450 } {
MD++ QLM_cutoff = 0.38
}
if { $T >= 550 } {
MD++ QLM_cutoff = 0.43
}
#MD++ N_lgst_index = -2

set umbordergroup 102506
MD++ UMB_order_group = $umbordergroup

if { 1 } {
  set nx 0.0
  set ny [expr 1.0]
  set nz 0.0
  MD++ SHtoR
  set NP [MD++_Get "NP"]

  set Height [ MD++_Get H_22 ]
  set dH [ expr 1.0*$Height/18.0 ]
  set r0x 0.0
  set r0y 0.0
  set r0z 0.0
}

if { 1 } {
#  puts "set atoms on the slip plane to group 102506 ..."
  set atoms 0
  for { set i 0 } { $i < $NP } { incr i } {
     set ind [expr $i*3]
     set rx [MD++_Get R($ind)]
     set ind [expr $i*3+1]
     set ry [MD++_Get R($ind)]
     set ind [expr $i*3+2]
     set rz [MD++_Get R($ind)]
     set rdotn [expr ($rx-$r0x)*$nx + ($ry-$r0y)*$ny + ($rz-$r0z)*$nz]
     set temp [expr $rdotn/($dH)]
     if { $temp > -0.5 && $temp < 1.5 } {
         MD++_Set group($i) $umbordergroup
         set atoms [ expr $atoms + 1 ]
#         puts "$atoms $i"
     }
  }
}
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
MD++ incnfile=../NEB/NEBinit.cn readcn setconfig2
alloc_QLM $T

MD++ incnfile=../NEB/NEBinit.cn readcn
openwindow

set EPOT0 [MD++_Get EPOT ]

set fp [ open "lam_array.dat" w ]
puts $fp "\["
puts $fp "200"

set index 0
set index1 0

for { set iter 0 } { $iter <=25 } { incr iter 1 } {

MD++ incnfile="../NEB/chain_no_${iter}.cn" readcn 
alloc_QLM $T
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
