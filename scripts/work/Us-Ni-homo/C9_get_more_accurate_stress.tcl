# -*-shell-script-*-
# MD code of FCC Cu dislocation
#  equlibration run
#------------------------------------------------------------

if { $argc == 0 } {
   set T 300
   set epsilon 0.110
} elseif { $argc > 0 } {
   set T [ lindex $argv 0 ]
   set epsilon [ lindex $argv 1 ]
} 

set epVAL [ format %03d [expr int($epsilon*1000) ] ]

proc relax_freebox { } { MD++ {
# Conjugate-Gradient relaxation
conj_ftol = 1e-10 conj_itmax = 10000 conj_fevalmax = 10000
conj_fixbox = 0
conj_fixboxvec = [ 0 1 1
                   1 0 1
                   1 1 0 ]
relax
} }

MD++ srand48bytime
MD++ initvelocity_type="Gaussian"


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

# initialize proc
MD++ {
#setnolog
setoverwrite
}
MD++ dirname = runs/ho-cu-DN/${T}K_${epVAL}eps_112

# plotting proc
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
#-------------------------------------------------------------
#Plot Configuration
atomradius = 1.0 bondradius = 0.3 bondlength = 0
atomcolor = cyan highlightcolor = purple backgroundcolor = gray
bondcolor = red fixatomcolor = yellow
color00 = "orange"  color01 = "red"    color02 = "green"
color03 = "magenta" color04 = "cyan"   color05 = "purple"
color06 = "gray80"  color07 = "white"  color08 = "blue"
#------------------------
#Use CSD to plot
plot_color_axis = 2  #2: use CSD (default 0: use local energy)
#NCS = 12 #number of neighboring atoms in CSD analysis (default: 12 - fcc)
plot_color_windows = [ 3
                       0    20   1 
                       5    10  4  #color04 = cyan
                       10   20  6  #color06 = gray80
                     ]
plot_atom_info = 1  plotfreq = 100
fixatomcolor = red
rotateangles = [ 0 0 0 1.5 ]
win_width = 600 win_height = 600
#openwin alloccolors rotate saverot refreshnnlist eval plot
}

# MD setting
MD++ {
#-------------------------------------------------------------
#Molecular Dynamics setting (equilibration run, NVE)
timestep = 0.0005 # (ps)
atommass = 63.55  #Mo (g/mol)
wallmass = 8e2
#### vt2=1e26
boxdamp=9e-3
NHChainLen = 4 NHMass=[0.5e-1 1e-2 1e-2 1e-2 ]
# for equilibration
#srand48bytime             #randomize random number generator
DOUBLE_T = 1 initvelocity
integrator_type = 1        #Velocity Verlet Algorithm
saveprop = 1 savepropfreq = 100 openpropfile 
}

relax_freebox
MD++ finalcnfile = struc0K.cn writecn saveH
MD++ T_OBJ = $T DOUBLE_T = 0 initvelocity
MD++ ensemble_type="NPTC" integrator_type="Gear6"
MD++ totalsteps = 100000 run
set H11_0 [datafile_process "prop.out" 10 0.3 1 "AVE" ]
set H22_0 [datafile_process "prop.out" 14 0.3 1 "AVE" ]
set H33_0 [datafile_process "prop.out" 18 0.3 1 "AVE" ]

MD++ {
output_fmt="curstep EPOT KATOM Tinst TSTRESSinMPa_xx TSTRESSinMPa_yy TSTRESSinMPa_zz TSTRESSinMPa_xy TSTRESSinMPa_yz TSTRESSinMPa_zx H_11 H_12 H_13 H_21 H_22 H_23 H_31 H_32 H_33"
}

MD++ incnfile = "pure_shear_stress.cn" readcn
set H11_new [ MD++_Get H_11 ]
set H12_new [ MD++_Get H_12 ]
set H13_new [ MD++_Get H_13 ]
set H22_new [ MD++_Get H_22 ]
set H23_new [ MD++_Get H_23 ]
set H33_new [ MD++_Get H_33 ]

MD++ T_OBJ = $T DOUBLE_T =0 initvelocity
MD++ ensemble_type="NVTC" integrator_type="Gear6"

MD++ totalsteps=100000
MD++ run

set sigxx [datafile_process "prop.out" 4 0.3 1 "AVE" ]
set sigyy [datafile_process "prop.out" 5 0.3 1 "AVE" ]
set sigzz [datafile_process "prop.out" 6 0.3 1 "AVE" ]
set sigxy [datafile_process "prop.out" 7 0.3 1 "AVE" ]
set sigyz [datafile_process "prop.out" 8 0.3 1 "AVE" ]
set sigzx [datafile_process "prop.out" 9 0.3 1 "AVE" ]

set e_xx [ expr ($H11_new-$H11_0)/$H11_0 ]
set e_yy [ expr ($H22_new-$H22_0)/$H22_0 ]
set e_zz [ expr ($H33_new-$H33_0)/$H33_0 ]
set e_xy [ expr ($H12_new)/$H22_0 ]

set fp [ open "strain_long.dat" w ]
puts $fp "$e_xx $e_xy $e_yy $e_zz"
close $fp

set fp [ open "stress_long.dat" w ]
puts $fp "$sigxx $sigyy $sigzz $sigxy $sigzx $sigyz"
close $fp

MD++ quit 
#------------------------------------------------------------
