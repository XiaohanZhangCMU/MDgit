# -*-shell-script-*-
# MD code for making FCC Cu pure shear stress 
#------------------------------------------------------------

set maxiter 100

MD++ {
setnolog
setoverwrite
}

set check [ glob -nocomplain runs/shear_silicon ]
set key [ llength $check ]
if { $key < 1 } {
exec mkdir runs/shear_silicon
}


MD++ dirname = runs/ho-si-DN/0K_stressstrain
MD++ {
#--------------------------------------------
#Read in potential file
#
potfile = $::env(MDPLUS_DIR)/potentials/w_pot readpot
#------------------------------------------------------------
#make a crystal structure 
#
crystalstructure = diamond-cubic latticeconst = 5.4309529817532409 #(A) for Si
latticesize = [  1  1  0  12
                -1  1  0  12
                 0  0  1  16 ]
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
set C11 170000
set factor 0.7

for { set itereps 0 } { $itereps <= 20 } { incr itereps 1 } {

puts "\n\n epsilon=$epsilon"
#exec sleep 1

set H12_fix [ expr 1.0*$H22_0*$epsilon ]
MD++ H_12 = $H12_fix
MD++ eval
MD++ plot

# adjust other strains in small steps
for { set iter 0 } { $iter <= $maxiter } { incr iter 1 } {

set sig_xx [ MD++_Get TSTRESSinMPa_xx ]
set sig_yy [ MD++_Get TSTRESSinMPa_yy ]
set sig_zz [ MD++_Get TSTRESSinMPa_zz ]
set sig_xy [ MD++_Get TSTRESSinMPa_xy ]
set sig_xz [ MD++_Get TSTRESSinMPa_xz ]
set sig_yz [ MD++_Get TSTRESSinMPa_yz ]
set e_xx [ expr $sig_xx / $C11 ]
set e_yy [ expr $sig_yy / $C11 ]
set e_zz [ expr $sig_zz / $C11 ]
puts "iter = $iter"
puts "sig_xx = $sig_xx sig_yy = $sig_yy sig_zz = $sig_zz"
puts "sig_xy = $sig_xy sig_xz = $sig_xz sig_yz = $sig_yz"
puts "e_xx = $e_xx e_yy = $e_yy e_zz = $e_zz"

set H11_cur [ MD++_Get H_11 ]
set H22_cur [ MD++_Get H_22 ]
set H33_cur [ MD++_Get H_33 ]
set H12_new [ MD++_Get H_12 ]

set H11_new [ expr ${H11_cur}*(1.0+$e_xx*$factor) ]
MD++ H_11 = ${H11_new}
set H22_new [ expr ${H22_cur}*(1.0+$e_yy*$factor) ]
MD++ H_22 = ${H22_new}
set H33_new [ expr ${H33_cur}*(1.0+$e_zz*$factor) ]
MD++ H_33 = ${H33_new}
set iter_ [ expr $iter-1 ]
if { $iter_ == $maxiter } {
MD++ conj_ftol = 2e-6 conj_fixbox = 1 relax
}
MD++ eval
}

set SIGXY [ expr int($sig_xy) ]
MD++ conj_ftol = 2e-6 conj_fixbox = 1 relax
MD++ finalcnfile = "0K_${SIGXY}MPa_relaxed.cn" writecn
MD++ finalcnfile = "NEBinit.cn" writecn
set e_xx [ expr ($H11_new-$H11_0)/$H11_0 ]
set e_yy [ expr ($H22_new-$H22_0)/$H22_0 ]
set e_zz [ expr ($H33_new-$H33_0)/$H33_0 ]
set e_xy [ expr ($H12_new)/$H22_0 ]
set fp [ open "strain.dat" a+ ]
puts $fp "$e_xx $e_xy $e_yy $e_zz"
close $fp
set fp [ open "strain_info.dat" w ]
puts $fp "e_xx e_xy e_yy e_zz"
close $fp

set fp [ open "stress.dat" a+ ]
puts $fp "$sig_xx $sig_yy $sig_zz $sig_xy $sig_xz $sig_yz"
close $fp
set fp [ open "stress_info.dat" w ]
puts $fp "sig_xx sig_yy sig_zz sig_xy sig_xz sig_yz"
close $fp

set epsilon [ expr $epsilon+0.01 ]

}

MD++ quit
#------------------------------------------------------------
