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
MD++ dirname = runs/ho-Ni-DN/${T}K_${epVAL}eps_112

# plotting proc
MD++ {
#--------------------------------------------
#Read in potential file
#
potfile = ~/Planet/Libs/MD++.svn/potentials/EAMDATA/eamdata.Ni.Rao99

eamgrid = 5000 readeam NNM = 600

#------------------------------------------------------------
#Read in structure (created by mo-edge-make.script)
#
#incnfile = au-DN-0.cn readcn
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

set H11_0 [ MD++_Get H_11 ]
set H12_0 [ MD++_Get H_12 ]
set H13_0 [ MD++_Get H_13 ]
set H22_0 [ MD++_Get H_22 ]
set H23_0 [ MD++_Get H_23 ]
set H33_0 [ MD++_Get H_33 ]

set H12_fix [ expr $epsilon*$H22_0 ]
MD++_Set H_12 $H12_fix

MD++ conj_fixboxvec = \[ 0 1 1  1 0 1  1 1 0 \]
MD++ conj_fixbox=0 relax
MD++ fixboxvec = \[ 0 1 1 1 0 1 1 1 0 \]
MD++ relax
MD++ eval saveH
MD++ {
output_fmt="curstep EPOT KATOM Tinst TSTRESSinMPa_xx TSTRESSinMPa_yy TSTRESSinMPa_zz TSTRESSinMPa_xy TSTRESSinMPa_yz TSTRESSinMPa_zx H_11 H_12 H_13 H_21 H_22 H_23 H_31 H_32 H_33"
}

MD++ T_OBJ = $T DOUBLE_T =1 initvelocity
MD++ ensemble_type="NPTC" integrator_type="Gear6"

MD++ totalsteps=12500
MD++ run

set sigxx [datafile_process "prop.out" 4 0.5 1 "AVE" ]
set sigyy [datafile_process "prop.out" 5 0.5 1 "AVE" ]
set sigzz [datafile_process "prop.out" 6 0.5 1 "AVE" ]
set sigxy [datafile_process "prop.out" 7 0.5 1 "AVE" ]
set sigyz [datafile_process "prop.out" 8 0.5 1 "AVE" ]
set sigzx [datafile_process "prop.out" 9 0.5 1 "AVE" ]
set H11 [datafile_process "prop.out" 10 0.5 1 "AVE" ]
set H12 [datafile_process "prop.out" 11 0.5 1 "AVE" ]
set H13 [datafile_process "prop.out" 12 0.5 1 "AVE" ]
set H21 [datafile_process "prop.out" 13 0.5 1 "AVE" ]
set H22 [datafile_process "prop.out" 14 0.5 1 "AVE" ]
set H23 [datafile_process "prop.out" 15 0.5 1 "AVE" ]
set H31 [datafile_process "prop.out" 16 0.5 1 "AVE" ]
set H32 [datafile_process "prop.out" 17 0.5 1 "AVE" ]
set H33 [datafile_process "prop.out" 18 0.5 1 "AVE" ]
puts "$H11 $H12 $H13" 
puts "$H21 $H22 $H23"
puts "$H31 $H32 $H33"

MD++ H_11 = $H11 H_22 = $H22 H_33 = $H33
MD++ H_13 = $H13 H_31 = $H31 H_23 = $H23 H_32 = $H32
MD++ H_12 = $H12_fix H_21 = $H21
MD++ ensemble_type="NVTC" totalsteps = 15000
MD++ writeall = 1 
MD++ initvelocity

### start steepest decent
set C11 200000
set factor 0.64
# increase three times bigger
set maxiter 60

set datafile "stress_trace.dat"
set fileID [open $datafile a+ ]
puts $fileID "% iter  XX    YY    ZZ    XY    YZ    ZX   (in MPa)"
close $fileID

for { set iter 0 } { $iter <= $maxiter } { incr iter 1 } {
#exec rm prop.out
set fp [ open "void.out" w ]
puts $fp ""
close $fp
exec cp void.out prop.out
MD++ openpropfile
set e_xx [ expr $sigxx/$C11 ]
set e_yy [ expr $sigyy/$C11 ]
set e_zz [ expr $sigzz/$C11 ]
puts "iter = $iter"
puts "sigxx = $sigxx sigyy = $sigyy sigzz = $sigzz"
puts "sigxy = $sigxy sigzx = $sigzx sigyz = $sigyz"
puts "e_xx = $e_xx e_yy = $e_yy e_zz = $e_zz"

set H11_cur [ MD++_Get H_11 ]
set H22_cur [ MD++_Get H_22 ]
set H33_cur [ MD++_Get H_33 ]
set H12_new [ MD++_Get H_12 ]

MD++ finalcnfile=inter${iter}.cn writecn

set H11_new [ expr ${H11_cur}*(1.0+$e_xx*$factor) ]
MD++ H_11 = ${H11_new}
set H22_new [ expr ${H22_cur}*(1.0+$e_yy*$factor) ]
MD++ H_22 = ${H22_new}
set H33_new [ expr ${H33_cur}*(1.0+$e_zz*$factor) ]
MD++ H_33 = ${H33_new}

MD++ DOUBLE_T = 0 initvelocity run

set sigxx [datafile_process "prop.out" 4 0.2 1 "AVE" ]
set sigyy [datafile_process "prop.out" 5 0.2 1 "AVE" ]
set sigzz [datafile_process "prop.out" 6 0.2 1 "AVE" ]
set sigxy [datafile_process "prop.out" 7 0.2 1 "AVE" ]
set sigyz [datafile_process "prop.out" 8 0.2 1 "AVE" ]
set sigzx [datafile_process "prop.out" 9 0.2 1 "AVE" ]

set datafile "stress_trace.dat"
set fileID [open $datafile a+ ]
#puts $fileID "$iter"
#puts $fileID "XX=$sigxx YY=$sigyy ZZ=$sigzz"
#puts $fileID "XY=$sigxy YZ=$sigyz ZX=$sigzx"
puts $fileID "$iter $sigxx $sigyy $sigzz $sigxy $sigyz $sigzx"
close $fileID


}

set SIGXY [ expr int($sigxy) ]
MD++ conj_fixbox=1 conj_ftol = 1e-12
MD++ conj_itmax = 9000 conj_fevalmax = 9000
MD++ relax relax
MD++ writeall = 1
MD++ finalcnfile = "pure_shear_stress.cn" writecn
set e_xx [ expr ($H11_new-$H11_0)/$H11_0 ]
set e_yy [ expr ($H22_new-$H22_0)/$H22_0 ]
set e_zz [ expr ($H33_new-$H33_0)/$H33_0 ]
set e_xy [ expr ($H12_new)/$H22_0 ]

set fp [ open "strain.dat" w ]
puts $fp "$e_xx $e_xy $e_yy $e_zz"
close $fp

set fp [ open "stress.dat" w ]
puts $fp "$sigxx $sigyy $sigzz $sigxy $sigzx $sigyz"
close $fp

set fp [ open "strain_info.dat" w ]
puts $fp "e_xx e_xy e_yy e_zz"
close $fp

set fp [ open "stress_info.dat" w ]
puts $fp "sig_xx sig_yy sig_zz sig_xy sig_zx sig_yz"
close $fp

#exec cp ${T}K_${P}MPa.cn ..
#exec cp ${T}K_${P}MPa_R2.cn ..
MD++ quit 
#MD++ run
#MD++ quit
#------------------------------------------------------------
