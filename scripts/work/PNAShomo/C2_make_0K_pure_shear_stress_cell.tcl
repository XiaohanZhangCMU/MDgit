# MD++ code for making 0K pure shear stress cell
#------------------------------------------------------------

if { $argc == 0 } {
   set epsilon 0.110
} elseif { $argc > 0 } {
   set epsilon [ lindex $argv 0 ]
} 
   set epVAL [ format %03d [ expr int($epsilon*1000) ] ]


MD++ srand48bytime
MD++ initvelocity_type="Gaussian"

# initialize MD++
MD++ {
##
setnolog
##
setoverwrite
}
MD++ dirname = runs/ho-cu-DN/0K_${epVAL}eps_112
MD++ {
#--------------------------------------------
#Read in potential file
#
potfile = ~/Planet/Libs/MD++.svn/potentials/EAMDATA/eamdata.CuMishin

 eamgrid = 5000 readeam NNM = 600

#------------------------------------------------------------
#Read in structure (created by mo-edge-make.script)
#
#incnfile = au-DN-0.cn readcn
latticestructure = face-centered-cubic latticeconst = 3.615
latticesize = [  1 -2  1  8
                 1  1  1  6
                 1  0 -1  13 ]
makecn finalcnfile = au-DN-0.cn writecn
#input = [ 1 -10 10 0.4 10 -10 10 ] fixatoms_by_position
#removefixedatoms
#input = [ 1 -10 10 -10 -0.4 -10 10 ] fixatoms_by_position
#removefixedatoms
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
#plot_color_axis = 0  #0: use local energy (default)
#plot_color_windows = [ 2
#                       -6.725  -6.0 8 #color08 = blue
#                       -6.0     0   6 #color06 = gray80
#                      ]
#------------------------
#Use CSD to plot
plot_color_axis = 2  #2: use CSD (default 0: use local energy)
#NCS = 12 #number of neighboring atoms in CSD analysis (default: 12 - fcc)
plot_color_windows = [ 3
                       0    20   1 
                       5    10  4  #color04 = cyan
                       10   20  6  #color06 = gray80
                     ]
plot_atom_info = 1  plotfreq = 10 
fixatomcolor = red
rotateangles = [ 0 0 0 1.5 ]
#
win_width = 600 win_height = 600

openwin alloccolors rotate saverot refreshnnlist eval plot

}

MD++ eval saveH

set C11 170000
set factor 0.7
set maxiter 100

MD++ conj_ftol = 2e-10 conj_fixbox = 1 relax
set H11_0 [ MD++_Get H_11 ]
set H22_0 [ MD++_Get H_22 ]
set H33_0 [ MD++_Get H_33 ]

set H12_fix [ expr 1.0*$H22_0*$epsilon ]
MD++ H_12 = $H12_fix

MD++ plot
#MD++ sleep
#MD++ incnfile = "0K_-2500MPa_1strelax.cn" readcn
#MD++ conj_fixbox = 1
#MD++ eval

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
MD++ conj_ftol = 2e-10 conj_fixbox = 1 relax
MD++ finalcnfile = "0K_${SIGXY}MPa_relaxed.cn" writecn
MD++ finalcnfile = "pure_shear_stress.cn" writecn
set e_xx [ expr ($H11_new-$H11_0)/$H11_0 ]
set e_yy [ expr ($H22_new-$H22_0)/$H22_0 ]
set e_zz [ expr ($H33_new-$H33_0)/$H33_0 ]
set e_xy [ expr ($H12_new)/$H22_0 ]
set fp [ open "strain.dat" w ]
puts $fp "$e_xx $e_xy $e_yy $e_zz"
close $fp
set fp [ open "strain_info.dat" w ]
puts $fp "e_xx e_xy e_yy e_zz"
close $fp

set fp [ open "stress.dat" w ]
puts $fp "$sig_xx $sig_yy $sig_zz $sig_xy $sig_xz $sig_yz"
close $fp
set fp [ open "stress_info.dat" w ]
puts $fp "sig_xx sig_yy sig_zz sig_xy sig_xz sig_yz"
close $fp
#exec mkdir NEBconfig_free
#exec cp NEBinit.cn NEBconfig_free/
#puts ${TSTRESS_xx_cur}
#MD++ quit 
#MD++ run
#MD++ quit
#------------------------------------------------------------
