
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
plot_color_windows = [ 0
                       0.6 9.4   1  
                       0.3  10   5
                       10 20   6
                       20 50   8
                       0  0.6  4
                     ]
# plot_limits = [ 1 -10 10 -10 10 0.38 10 ]

#  plot_limits = [ 1 -10 -0.47 -10 10 -10 10 ]
# plot_limits = [ 1 -10 10 -10 10 -10 -0.38 ]

#
#xiaohan

#plot_color_windows = 5 plot_limits = [ 1 -10 10 -0.001 0.028 -10 10 ]
#plot_color_windows = 0 
#plot_atom_info = 1 # reduced coordinates of atoms
plot_atom_info = 2 # real coordinates of atoms
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

proc openwindow { } { 
setup_window
MD++ openwin alloccolors rotate saverot eval plot
}


if { $argc > 0 } {
set NO [ lindex $argv 0 ]
set eps [ lindex $argv 1 ]
set T [ lindex $argv 2 ]
} else {
set NO 2000
#set NO 30
set eps 110
set T 0
}


MD++ setnolog
MD++ setoverwrite
MD++ dirname = runs/ho-cu-DN/${T}K_${eps}eps_112/NEB
exec cp ../../../../scripts/work/PNAShomo/NEB/NEBinit.cn .

if { $eps <= 0.12 } {
  exec cp ../../../../scripts/work/PNAShomo/NEB/NEBfinal.cn .
} else {
  exec cp ../../../../scripts/work/PNAShomo/NEB/NEBfinal1.cn NEBfinal.cn
}

#--------------------------------------------
#Read in potential file
#MD++ potfile = "~/Planet/Libs/MD++.svn/trunk/potentials/EAMDATA/eamdata.CuMishin" 
#MD++ eamgrid=5000 readeam NNM = 600

MD++ potfile=~/Planet/Libs/MD++.svn/potentials/EAMDATA/eamdata.CuMishin
MD++ eamgrid = 5000 readeam NNM = 600

#------------------------------------------------------------
#Read in structure
#
MD++ incnfile = ../pure_shear_stress.cn readcn
set H11 [ MD++_Get H_11 ]
set H12 [ MD++_Get H_12 ]
set H13 [ MD++_Get H_13 ] 
set H21 [ MD++_Get H_21 ]
set H22 [ MD++_Get H_22 ]
set H23 [ MD++_Get H_23 ]
set H31 [ MD++_Get H_31 ]
set H32 [ MD++_Get H_32 ]
set H33 [ MD++_Get H_33 ]
MD++ conj_ftol = 2e-10 conj_fixbox = 1 relax relax eval
MD++ finalcnfile = NEBinit.cn writecn
MD++ incnfile = NEBinit.cn readcn setconfig1 #A
MD++ incnfile = NEBfinal.cn readcn 


MD++ H_11 = $H11 H_12 = $H12 H_13 = $H13
MD++ H_21 = $H21 H_22 = $H22 H_23 = $H23
MD++ H_31 = $H31 H_32 = $H32 H_33 = $H33
MD++ finalcnfile = NEBfinal.cn writecn 
MD++ incnfile = NEBfinal.cn readcn


#setup_window
#openwindow
MD++ conj_ftol = 1e-4 conj_itmax = 3800 conj_fevalmax =  120
MD++ conj_fixbox = 1  relax
MD++ finalcnfile = "debug.cfg" writeatomeyecfg
MD++ sleep



MD++ eval setconfig2 #B

MD++ incnfile = NEBinit.cn readcn #A


MD++ {
#------------------------------------------------------------
#
#Plot Configuration
#
atomradius = 0.9 bondradius = 0.3 bondlength = 2.8285 #for Si
atomcolor = orange highlightcolor = purple  bondcolor = red
backgroundcolor = gray fixatomcolor = yellow
color00 = "orange"  color01 = "purple" color02 = "green"
color03 = "magenta" color04 = "cyan"   color05 = "purple"
color06 = "gray80"  color07 = "white"
plot_limits = [ 1 -10 10 -10 10  -0.05  0.05 ]
plot_color_windows = [ 2
                          -10 -4.5 0  #color00 = orange
                          -4.5 10  1  #color01 = purple
                     ]
plot_atom_info = 1
plotfreq = 10
rotateangles = [ 0 0 0 1.5 ]
win_width = 600 win_height = 600
#openwin alloccolors rotate saverot refreshnnlist eval plot
}

MD++ fixallatoms constrain_fixedatoms freeallatoms
MD++ chainlength = 25 totalsteps = $NO equilsteps = 20

MD++ timestep = 0.04 printfreq = 2
MD++ allocchain initRchain

MD++ {
nebspec = [ 0 1 0 1 1 ]
#nebspec = [ 0 1 0 0 1 ]
stringrelax
}

MD++ finalcnfile = stringrelax.chain.cn writeRchain
MD++ quit

