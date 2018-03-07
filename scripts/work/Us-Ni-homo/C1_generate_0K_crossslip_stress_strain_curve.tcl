# -*-shell-script-*-
# MD code for making FCC Cu pure shear stress 
#------------------------------------------------------------

set maxiter 100

MD++ {
setnolog
setoverwrite
}

set check [ glob -nocomplain runs/shear_copper ]
set key [ llength $check ]
if { $key < 1 } {
exec mkdir runs/shear_copper
}


MD++ dirname = runs/ho-Ni-DN/0K_cross_slip_stressstrain
MD++ {
#--------------------------------------------
#Read in potential file
#
potfile = ~/Planet/Libs/MD++.svn/potentials/EAMDATA/eamdata.Ni.Rao99

eamgrid = 5000 readeam NNM = 600
#------------------------------------------------------------
#make a crystal structure 
#
latticestructure = face-centered-cubic latticeconst = 3.52
latticesize = [  1 -2  1  8
                 1  1  1  6
                 1  0 -1  13 ]
makecn finalcnfile = au-DN-0.cn writecn
#-------------------------------------------------------------
#Plot Configuration
atomradius = 1.0 bondradius = 0.3 bondlength = 0
atomcolor = cyan highlightcolor = purple backgroundcolor = gray
bondcolor = red fixatomcolor = yellow
color00 = "orange"  color01 = "red"    color02 = "green"
color03 = "magenta" color04 = "cyan"   color05 = "purple"
color06 = "gray80"  color07 = "white"  color08 = "blue"
#------------------------
#Use local energy to plot
plot_color_axis = 2  #2: use CSD (default 0: use local energy)
#NCS = 12 #number of neighboring atoms in CSD analysis (default: 12 - fcc)
plot_color_windows = [ 1
                       0    20   1 
                     ]
plot_atom_info = 1  plotfreq = 10 
fixatomcolor = red
rotateangles = [ 0 0 0 1.5 ]
win_width = 600 win_height = 600
openwin alloccolors rotate saverot refreshnnlist eval plot

#-------------------------------------------------------------
#Molecular Dynamics setting (equilibration run, NVE)
timestep = 0.005 # (ps)
}
MD++ eval saveH

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


set epsilon 0.0
set C11 200000
set modulus_xx 450e3 ; #MPa
set modulus_zz 450e3 ; #MPa
set C44 100e3 ; #MPa

#set C11 170000
set factor 0.7

for { set itereps 0 } { $itereps <= 200 } { incr itereps 1 } {

puts "\n\n epsilon=$epsilon"
#exec sleep 1

set scaleStress [expr 2262.741700*$epsilon]
set targetStressMpa_xx [expr 0 ];
set targetStressMpa_yy [expr  $scaleStress];
set targetStressMpa_zz [expr -1 * $scaleStress];
set targetStressMpa_xy [expr 0 ];
set targetStressMpa_zx [expr 424.264069 /2262.741700 * $scaleStress];
set targetStressMpa_yz [expr 800/2262.7414700 * $scaleStress];

set dsig_xx [expr 0]
set dsig_yy [expr 0]
set dsig_zz [expr 0]
set dsig_xy [expr 0]
set dsig_zx [expr 0]
set dsig_yz [expr 0]

MD++ eval
MD++ plot

# adjust other strains in small steps
for { set iter 0 } { $iter <= $maxiter } { incr iter 1 } {

set sig_xx [ MD++_Get TSTRESSinMPa_xx ]
set sig_yy [ MD++_Get TSTRESSinMPa_yy ]
set sig_zz [ MD++_Get TSTRESSinMPa_zz ]
set sig_xy [ MD++_Get TSTRESSinMPa_xy ]
set sig_zx [ MD++_Get TSTRESSinMPa_zx ]
set sig_yz [ MD++_Get TSTRESSinMPa_yz ]

set dsig_xx [expr -1.0*( $targetStressMpa_xx -  $sig_xx )] ; 
set dsig_yy [expr -1.0*( $targetStressMpa_yy -  $sig_yy )]
set dsig_zz [expr -1.0*( $targetStressMpa_zz -  $sig_zz )] ;
set dsig_xy [expr -1.0*( $targetStressMpa_xy -  $sig_xy )]
set dsig_zx [expr -1.0*( $targetStressMpa_zx -  $sig_zx )] ; 
set dsig_yz [expr -1.0*( $targetStressMpa_yz -  $sig_yz )]

set e_xx [ expr $dsig_xx / $C11 ] ; set e_yy [ expr $dsig_yy / $C11 ] ; 
set e_zz [ expr $dsig_zz / $C11 ] ;
set e_xy [ expr $dsig_xy / $C44 ] ; set e_zx [ expr $dsig_zx / $C44 ] ; 
set e_yz [ expr $dsig_yz / $C44 ]

puts "iter = $iter"
puts "sig_xx = $sig_xx sig_yy = $sig_yy sig_zz = $sig_zz"
puts "sig_xy = $sig_xy sig_zx = $sig_zx sig_yz = $sig_yz"
puts "e_xx = $e_xx e_yy = $e_yy e_zz = $e_zz"

set H11_cur [ MD++_Get H_11 ]
set H22_cur [ MD++_Get H_22 ]
set H33_cur [ MD++_Get H_33 ]
set H13_cur [ MD++_Get H_13 ]
set H12_cur [ MD++_Get H_12 ]
set H23_cur [ MD++_Get H_23 ]

set H11_new [ expr ${H11_cur}*(1.0+$e_xx*$factor) ]; MD++ H_11 = ${H11_new}
set H22_new [ expr ${H22_cur}*(1.0+$e_yy*$factor) ]; MD++ H_22 = ${H22_new}
set H33_new [ expr ${H33_cur}*(1.0+$e_zz*$factor) ]; MD++ H_33 = ${H33_new}
set H13_new [ expr ${H13_cur} + ${H33_cur}*$e_zx*$factor ] ; MD++ H_13 = ${H13_new}
set H12_new [ expr ${H12_cur} + ${H22_cur}*$e_xy*$factor ] ; MD++ H_12 = ${H12_new}
set H23_new [ expr ${H23_cur} + ${H33_cur}*$e_yz*$factor ] ; MD++ H_23 = ${H23_new}

set iter_ [ expr $iter-1 ]
if { $iter_ == $maxiter } {
MD++ conj_ftol = 2e-12 conj_fixbox = 1 relax
}
MD++ eval
}

set SIGXY [ expr int($sig_xy) ]
MD++ conj_ftol = 2e-10 conj_fixbox = 1 relax
MD++ finalcnfile = "0K_${SIGXY}MPa_relaxed.cn" writecn
MD++ finalcnfile = "NEBinit.cn" writecn
set e_xx [ expr ($H11_new-$H11_0)/$H11_0 ]
set e_yy [ expr ($H22_new-$H22_0)/$H22_0 ]
set e_zz [ expr ($H33_new-$H33_0)/$H33_0 ]
set e_xy [ expr ($H12_new)/$H22_0 ]
set e_zx [ expr ($H13_new)/$H33_0 ]
set e_yz [ expr ($H23_new)/$H33_0 ]
set fp [ open "strain.dat" a+ ]
puts $fp "$e_xx $e_yy $e_zz $e_xy $e_zx $e_yz"
close $fp
set fp [ open "strain_info.dat" w ]
puts $fp "e_xx e_yy e_zz e_xy e_zx e_yz"
close $fp

set fp [ open "stress.dat" a+ ]
puts $fp "$sig_xx $sig_yy $sig_zz $sig_xy $sig_zx $sig_yz"
close $fp
set fp [ open "stress_info.dat" w ]
puts $fp "sig_xx sig_yy sig_zz sig_xy sig_zx sig_yz"
close $fp

set epsilon [ expr $epsilon+0.01 ]

}

MD++ quit
#------------------------------------------------------------
