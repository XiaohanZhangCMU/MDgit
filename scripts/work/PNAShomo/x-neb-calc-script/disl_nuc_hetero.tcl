# -*-shell-script-*-
# compute energy barrier of homogeneous dislocation nucleation in BCC W
# relax all other stress components (except yz) to zero
#
# make fs build=R SYS=mc2
# sw_mc2_mpich scripts/disl_nuc_hetero.tcl 0 0 001 1     -generates 0.01 0.02 0.03 strain perfect states A
# sw_mc2_mpich scripts/disl_nuc_hetero.tcl 1 0.030 001 1 -generates 0.03 state B.
# mpirun -np $ncpu sw_mc2_mpich scripts/disl_nuc_hetero.tcl 4 0.030 001 1 ----string_relax_parallel
source "$::env(MDPLUS_DIR)/scripts/Examples/Tcl/startup.tcl"

#*******************************************
# Definition of procedures
#*******************************************
proc initmd { n } {
#MD++ setnolog
MD++ setoverwrite
#relative to the  directory of lauching 
MD++ dirname = runs/ho-cu-DN/x-neb/cu-disl-homo-$n
#max number neigbors
MD++ NNM = 200
}

proc readpot { } { MD++ {
#--------------------------------------------
#Read in potential file
#
#potfile = $::env(MDPLUS_DIR)/potentials/w_pot readpot $$$$$$$$$$$$$$$$$$$$$$$
potfile = ~/Planet/Libs/MD++.svn/potentials/EAMDATA/eamdata.CuMishin
eamgrid = 5000 readeam NNM = 600
} }

# make sure the coordinate is right hand sided.

proc make_perfect_crystal { nx ny nz } {
#    MD++ crystalstructure = diamond-cubic latticeconst = 5.4309529817532409 #(A) for Si$$$$$$$$$$$$$$$$$$$$
#    MD++ latticesize = \[  1 1 0  $nx  -1 1 0  $ny  0 0 1  $nz \] $$$$$$$$$$$$$$$$

    MD++ crystalstructure = face-centered-cubic latticeconst = 3.615
    MD++ latticesize = \[  1 -2  1  $nx  1  1  1  $ny   1  0 -1  $nz ]

    MD++ makecrystal #finalcnfile = perf.cn writecn #eval
}


proc make_dislocation_loop { flag epsilon } {
    set Lx [MD++_Get H_11]
    set Ly [MD++_Get H_22]
    set Lz [MD++_Get H_33]
    set store 1
    set a 3.615
    set bx 0.3333
    set by 0
    set bz 0


    set bx 0
    set by 0
    set bz [expr -sqrt(2)/2 ]
    
# 0.5[0 -1 1] = 1/sqrt(6) * [1 -2 1] +  1/sqrt(3) * [1 1 1]  + 1/sqrt(2)*[1 0 -1] : Full
    set bx [expr 0.5*sqrt(6)/2]
    set by 0
    set bz [expr -0.5*sqrt(2)/2 ]

# 0.5[0 -1 1] = 1/6[-1 -1 2] + 1/6[1 -2 1]
# 1/6[ -1 -1 2 ] = 1/sqrt(6) * [1 -2 1] * (0.5)  +  1/sqrt(3) * [1 1 1] *(0)  + 1/sqrt(2)*[1 0 -1] *(-1.5)
    set bx [expr sqrt(6)*0.5/6.0]
    set by 0
    set bz [expr sqrt(2)*(-1.5)/6.0]

    set bx [expr 1/sqrt(6) ]
    set by 0
    set bz 0


    set lx 1
    set ly 0
    set lz 0
    set nx 0
    set ny 1
    set nz 0
    set x0 0 
    set y0 1.0435
    set z0 0
#0.0300 x

    if { $epsilon <= 0.0300 } {  
        set Ra [expr $Lx * 0.42]

    } elseif { $epsilon <= 0.0900 } {
	set Ra [expr $Lx * 0.182]
    } elseif { $epsilon <= 0.0960 } {
	set Ra [expr $Lx * 0.172]
    } elseif { $epsilon <= 0.1020 } {
        set Ra [expr $Lx * 0.162]
    } elseif { $epsilon <= 0.1080 } {
        set Ra [expr $Lx * 0.152]
    } elseif { $epsilon <= 0.1140 } {
        set Ra [expr $Lx * 0.142]
    } elseif { $epsilon <= 0.1200 } {
        set Ra [expr $Lx * 0.140]
    } elseif { $epsilon <= 0.1260 } {
        set Ra [expr $Lx * 0.138]
    } elseif { $epsilon <= 0.1360 } {
        set Ra [expr $Lx * 0.132]
    } elseif { $epsilon <= 0.1460 } {
        set Ra [expr $Lx * 0.130]

    } else {
	puts 'unexpected load'
	MD++ quit
    }
    set Rb [expr 0.8*$Ra]
    MD++ input= \[ 1 $Ra $Rb $a $bx $by $bz $lx $ly $lz $nx $ny $nz $x0 $y0 $z0 $store \] 
    MD++ makedislellipse  
}

#--------------------------------------------
proc relax_fixbox { } { MD++ {
# Conjugate-Gradient relaxation
conj_ftol = 1e-4 conj_itmax = 1000 conj_fevalmax = 1000
conj_fixbox = 1
relax
} }
#end of proc relax_fixbox

#--------------------------------------------
proc relax_freebox { } { MD++ {
# Conjugate-Gradient relaxation
conj_ftol = 1e-4 conj_itmax = 1000 conj_fevalmax = 1000
conj_fixbox = 0
conj_fixboxvec = [ 0 0 1
                   1 0 1
                   0 0 0 ]
relax
} }
#end of proc relax_fixbox

proc setup_window { } { MD++ {
#------------------------------------------------------------
#colors for Central symmetry view
color00 = "red" color01 = "blue" color02 = "green"
color03 = "magenta" color04 = "cyan" color05 = "purple"
color06 = "gray80" color07 = "white" color08 = "orange"
#--------------------------------------------
# Plot Configuration
#
atomradius = [1.0 0.78] bondradius = 0.3 bondlength = 0 #for Si
#atomradius = 0.9 bondradius = 0.3 bondlength = 0 #2.8285 #for Si
atomcolor = orange highlightcolor = purple  bondcolor = red

plot_color_axis = 2 
# NCS = 4
plot_color_windows = [ 1
                       0.6 9.4   1  
                       3  5   5
                       10 20   6
                       20 50   8
                       0  0.6  4
                     ]
plot_atom_info = 2 # real coordinates of atoms
#plot_atom_info = 3 # energy of atoms
#plot_highlight = [ 0 0 1 2 3 4 5 6 7 8 9 ]
plotfreq = 10
#

rotateangles = [ -0 90 0 1.2 ]

#rotateangles = [ 0 0 0 1.7 ]
#rotateangles = [ 0 -90 0 1.7 ]
win_width = 600 win_height = 600
openwin alloccolors rotate saverot plot
#plot_color_axis = 0 input = [ -8 -3 10] GnuPlotHistogram
#plot_color_axis = 2 input = [ 0.6 50 50 ] GnuPlotHistogram
#plot_limits = [ 1 -10 10 -0.25 0.25 -10 10 ]
} }

proc openwindow { } { 
setup_window
MD++ openwin alloccolors rotate saverot eval plot
}

#--------------------------------------------
proc exitmd { } { MD++ quit }
#end of proc exitmd
#--------------------------------------------

#--------------------------------------------
proc setup_md { } { MD++ {     
T_OBJ = 300 #Kelvin #add by xiaohan

equilsteps = 0  totalsteps = 5000 timestep = 0.0001 # (ps)
atommass = 28.0855 # (g/mol)
DOUBLE_T = 1
saveprop = 1 savepropfreq = 100 openpropfile #run
savecn = 1 savecnfreq = 10000 openintercnfile
plotfreq = 100 printfreq = 100
#ensemble_type = "NPH" integrator_type = "Gear6" implementation_type = 0
ensemble_type = "NVE" integrator_type = "VVerlet" implementation_type = 0
vt2 = 1e28  #1e28 2e28 5e28
wallmass = 2e3     # atommass * NP = 14380
boxdamp = 1e-3     # optimal damping for 216 atoms and wallmass 1e-3
saveH # Use current H as reference (H0), needed for specifying stress
fixboxvec = [ 0 0 1
              1 0 1
              0 0 0 ]
output_fmt = "curstep EPOT KATOM Tinst HELM HELMP TSTRESS_xx TSTRESS_yy TSTRESS_zz H_11 H_22 H_33" 
} }
#end of proc setup_md


#*******************************************
# Main program starts here
#*******************************************
# status 0:
#        1:
#        2:
#
# read in status from command line argument
if { $argc == 0 } {
 set status 0
} elseif { $argc > 0 } {
 set status [lindex $argv 0]
}
puts "status = $status"

if { $argc <= 1 } {
 set n 0
} elseif { $argc > 1 } {
 set n [lindex $argv 1]
}
puts "n = $n"

if { $argc <= 2 } {
 set flag 0
} elseif { $argc > 2 } {
 set flag [lindex $argv 2]
}
puts "flag = $flag"

if { $argc <= 3 } {
 set opt 0
} elseif { $argc > 3 } {
 set opt [lindex $argv 3]
}
puts "opt = $opt"


if { $status == 0 } {
  # generate zero temperature normal stress-strain curve
  MD++ setnolog
  initmd $n
  readpot

  make_perfect_crystal 8 6 13

  MD++ finalcnfile = "w-perf.cn" writecn
  set strain_data_file "strain_surf_112.dat"
  set stress_data_file "stress_surf_112.dat"

  MD++ eval saveH conj_fixboxvec =  \[ 0 1 0  0 0 0   1 1 0 \]
  MD++ conj_fixbox = 1 conj_ftol = 2e-2 conj_itmax = 1000 conj_fevalmax = 2000
  MD++ relax finalcnfile = "0K_perfect.cn" writecn
  MD++ eval

  set H11_0 [ MD++_Get H_11 ] ; set H22_0 [ MD++_Get H_22 ] ; set H33_0 [ MD++_Get H_33 ]

  set epsilon 0.0
  set C11 520000  
  set C44 160000
  set factor 0.7

  set maxitereps 200
  set maxiter    100
  for { set itereps 0 } { $itereps <= $maxitereps } { incr itereps 1 } {

      puts "\n\n epsilon=$epsilon"

      MD++ eval plot
      set H12_fix [expr $epsilon*$H22_0]
      MD++ H_12 = $H12_fix

      if { 1 } {

      # adjust other strains in small steps
      for { set iter 0 } { $iter <= $maxiter } { incr iter 1 } {

          set sig_xx [ MD++_Get TSTRESSinMPa_xx ] ; set sig_yy [ MD++_Get TSTRESSinMPa_yy ]
          set sig_zz [ MD++_Get TSTRESSinMPa_zz ] ; set sig_xy [ MD++_Get TSTRESSinMPa_xy ]
          set sig_xz [ MD++_Get TSTRESSinMPa_xz ] ; set sig_yz [ MD++_Get TSTRESSinMPa_yz ]
          set e_xx [ expr $sig_xx / $C11 ] ; set e_yy [ expr $sig_yy / $C11 ] ; set e_zz [ expr $sig_zz / $C11 ]
          set e_xy [ expr $sig_xy / $C44 ] ; set e_xz [ expr $sig_xz / $C44 ] ; set e_yz [ expr $sig_yz / $C44 ]

          set H11_cur [ MD++_Get H_11 ] ; set H22_cur [ MD++_Get H_22 ] ; set H33_cur [ MD++_Get H_33 ] 
          set H31_cur [ MD++_Get H_31 ] ; set H12_cur [ MD++_Get H_12 ] ; set H23_cur [ MD++_Get H_23 ]

          set H11_new [ expr ${H11_cur}*(1.0+$e_xx*$factor) ] ;        MD++ H_11 = ${H11_new}
          set H22_new [ expr ${H22_cur}*(1.0+$e_yy*$factor) ] ;        MD++ H_22 = ${H22_new}
          set H33_new [ expr ${H33_cur}*(1.0+$e_zz*$factor) ] ;        MD++ H_33 = ${H33_new}
          set H31_new [ expr ${H31_cur} + ${H33_cur}*$e_xz*$factor ] ; MD++ H_31 = ${H31_new}
          #set H12_new [ expr ${H12_cur} + ${H22_cur}*$e_xy*$factor ] ; MD++ H_12 = ${H12_new}
          set H23_new [ expr ${H23_cur} + ${H33_cur}*$e_yz*$factor ] ; MD++ H_23 = ${H23_new}

          set H12_new ${H12_cur}

          if { $iter == [expr $maxiter + 1] } {
              MD++ conj_ftol = 2e-6 conj_fixbox = 1 relax
          }
          MD++ eval
     }
     set SIGYY [ expr int($sig_yy) ]
     #MD++ conj_ftol = 2e-6 conj_fixbox = 1 relax
     }
     relax_fixbox

     MD++ finalcnfile = "0K_${epsilon}_relaxed_surf${flag}.cn" writecn
     MD++ finalcnfile = "0K_${epsilon}_relaxed_surf${flag}.cfg" writeatomeyecfg
     set e_xx [ expr ($H11_new-$H11_0)/$H11_0 ] ; set e_yy [ expr ($H22_new-$H22_0)/$H22_0 ]
     set e_zz [ expr ($H33_new-$H33_0)/$H33_0 ] ; set e_yz [ expr ($H23_new)/$H33_0 ]
     set e_xy [ expr ($H12_new)/$H22_0 ] ;        set e_xz [ expr ($H31_new)/$H33_0 ]
 
     set fp [ open $strain_data_file a+ ] ; puts $fp "$e_xx $e_yy $e_zz $e_xy $e_xz $e_yz" ; close $fp
     set fp [ open "strain_info.dat" w ] ; puts $fp "e_xx e_yy e_zz e_xy e_xz e_yz" ; close $fp

     set fp [ open $stress_data_file a+ ] ; puts $fp "$sig_xx $sig_yy $sig_zz $sig_xy $sig_xz $sig_yz" ; close $fp
     set fp [ open "stress_info.dat" w ] ; puts $fp "sig_xx sig_yy sig_zz sig_xy sig_xz sig_yz" ;   close $fp

#     setup_window
#     openwindow
#     MD++ sleep 

     set epsilon [ expr $epsilon+0.002 ]
     set epsilon [format "%.4f" $epsilon]
  }

  MD++ sleep
  exitmd

# To plot shear stress strain curve
#    cd runs/w-disl-nuc-hetero-0
#    octave
#         load stress.dat ; load strain.dat
#         plot(-strain(:,2),stress(:,2)/1000,'o-');
#         xlabel('\epsilon_{yy}');  ylabel('\sigma_{yy}  (GPa)');
  
} elseif { $status == 1 } {
  # prepare and relax dislocation structure
  MD++ setnolog
  initmd 0
  readpot

  set epsilon $n
  set epVAL [ format %03d [expr int($epsilon*1000) ] ]

  if { $flag == 112 } {
    #set epsilon 0.04 
    #set epsilon 0.074
    set FEVALMAX 20
  } elseif { $flag == 111 } {
    #set epsilon 0.04
    #set epsilon 0.074
    set FEVALMAX 20
  } else {
     puts "unknown flag = $flag, must be 112 or 110"
     exitmd
  }


#make dislocation after reconstru but before the straining
  MD++ incnfile = "0K_0.0_relaxed_surf${flag}.cn" readcn

  make_dislocation_loop $flag $epsilon

  MD++ incnfile = "../../0K_${epVAL}eps_112/pure_shear_stress.cn" readcn
  MD++ commit_storedr
  
  MD++ finalcnfile = "w-loop-init.cn" writecn
  MD++ finalcnfile = "w-loop-init.cfg" writeatomeyecfg

#  setup_window
#  openwindow


  MD++ eval plot 
  # this is to confirm the structure is beyond critical
  MD++ conj_ftol = 1e-4 conj_itmax = 3800 conj_fevalmax =  20

  MD++ conj_fixbox = 1  relax

#  MD++ conj_ftol = 1e-5 conj_itmax = 3800 conj_fevalmax = 30
#  MD++ conj_fixbox = 1 relax

  MD++ finalcnfile = "w-loop-compression-112-eps${epVAL}.cn"  writecn
  MD++ finalcnfile = "w-loop-compression-112-eps${epVAL}.cfg" writeatomeyecfg

  MD++ eval plot

#  MD++ sleep
  exitmd
  
} elseif { $status == 2 } {
  # find minimum energy path for dislocation loop nucleation
  #MD++ setnolog
  initmd $flag-$n
  readpot

  set epsilon $n
  exec cp /home/xzhang11/Planet/Libs/MD++.svn/scripts/work/PNAShomo/NEB/NEBfinal.cn .
  MD++ incnfile ="NEBfinal.cn" readcn 
  MD++ finalcnfile = "NEBfinal.cfg" writeatomeyecfg
  
} elseif { $status == 3 } {

  # using parallel run to find minimum energy path for dislocation loop nucleation
  #MD++ setnolog
  set ncpu [ MD++_Get numDomains ]
  if { $ncpu < 1 } { set ncpu 1 }
  puts "ncpu = $ncpu"

  initmd $flag-$n.$ncpu

  set epsilon $n

  MD++ Slave_chdir #Slave change to same directory as Master (dirname)
#  readpot

  MD++ incnfile = "../w-disl-nuc-hetero-0/0K_${epsilon}_relaxed_surf${flag}.cn" readcn  saveH
  MD++ incnfile = "../w-disl-nuc-hetero-0/w-loop-compression-surf${flag}.cn" readcn   restoreH SHtoR  setconfig2
  MD++ eval
  MD++ incnfile = "../w-disl-nuc-hetero-0/0K_${epsilon}_relaxed_surf${flag}.cn" readcn  setconfig1
  MD++ eval 

  #setup_window
  #openwindow

  set chainlength [expr $ncpu - 1]
  MD++ chainlength = $chainlength  
  MD++ timestep = 0.01 printfreq = 10

  if { $opt == 0 } {
   if { $n >= 0.09 } {
    MD++ { nebspec = [ 0 1 0 1 0 0 ] totalsteps = 200 equilsteps = 2000 }
   } else {
    MD++ { nebspec = [ 0 1 0 1 0 0 ] totalsteps = 200 equilsteps = 2000 }
   }
  } else {
    MD++ { nebspec = [ 0 1 0 1 0 1 ] totalsteps = 5000 equilsteps = 2000 }
    ## [ relax_surround redistr_freq moveleftend moverightend yesclimbimage islengthconstant ]
  }
  MD++ Broadcast_EAM_Param

  MD++ fixallatoms constrain_fixedatoms freeallatoms
  MD++ allocchain_parallel initRchain_parallel

  if { $opt == 1 } {
    MD++ incnfile = "neb.chain.1000" readRchain_parallel
  }

#  MD++ incnfile = "../w-disl-nuc-hetero-$flag-$n.$ncpu/neb.chain.500" readRchain_parallel

#  MD++ { nebspec = [ 0 1 0 1 0 0 ] totalsteps = 20000 }
#  MD++ stringrelax_parallel
#  MD++ finalcnfile = neb.chain.500 writeRchain_parallel

#  MD++ incnfile = "neb.chain.500" readRchain_parallel
#  MD++ { nebspec = [ 0 1 0 1 0 1 ] totalsteps = 20000 }

  MD++ stringrelax_parallel

  if { $opt == 0 } {
    MD++ finalcnfile = neb.chain.1000 writeRchain_parallel
  } else {
    MD++ finalcnfile = neb.chain.2000 writeRchain_parallel
  }

  if { $ncpu == 1 }   {
        MD++ quit
  } else {
        MD++ quit_all
  }
  
} elseif { $status == 4 } {
  # using parallel run to find minimum energy path for dislocation loop nucleation
  # two step relaxation (first with length constant, then with length free)
  #MD++ setnolog
  set ncpu [ MD++_Get numDomains ]
  if { $ncpu < 1 } { set ncpu 1 }
  puts "ncpu = $ncpu"
  #set ncpu 24

  initmd $flag-$n.$ncpu

  set epsilon $n

  MD++ Slave_chdir #Slave change to same directory as Master (dirname)
  readpot

  set epVAL [ format %03d [expr int($epsilon*1000) ] ]
  MD++ incnfile = "../../0K_${epVAL}eps_112/pure_shear_stress.cn" readcn saveH
  MD++ incnfile = "../cu-disl-homo-0/w-loop-compression-112-eps${epVAL}.cn" readcn   restoreH SHtoR  setconfig2 
  MD++ eval
  MD++ incnfile = "../../0K_${epVAL}eps_112/pure_shear_stress.cn" readcn  setconfig1
  MD++ eval 

  # rename previous result files in order to not overwrite them
  set check [ glob -nocomplain "stringrelax.chain.cn.cpu00" ]
  set exist [ llength $check ]
  if { $exist > 0 } {
     for { set i 0 } { $i < $ncpu } { incr i 1 } {
         set ival [ format %02d $i ]
         exec cp stringrelax.chain.cn.cpu$ival Prev_stringrelax.chain.cn.cpu$ival
     }
     exec cp stringeng.out Prev_stringeng.out
  }

#  setup_window
  #openwindow

  set chainlength [expr $ncpu - 1]
  MD++ chainlength = $chainlength  
  MD++ timestep = 0.01 printfreq = 10
  #MD++ timestep = 0.002 printfreq = 10

  # First relaxation
  if { 0 } {
   if { $n >= 0.2 } {
    MD++ { nebspec = [ 0 1 0 1 0 1 ] totalsteps = 200 equilsteps = 2000 }
  ## [ relax_surround redistr_freq moveleftend moverightend yesclimbimage islengthconstant ]
   } else {
    MD++ { nebspec = [ 0 1 0 1 0 1 ] totalsteps = 1000 equilsteps = 2000 }
   }
  } else {
    MD++ nebspec = \[ 0 1 0 1 0 0 1 \] totalsteps = 10000 equilsteps = 10000
  }   

  MD++ Broadcast_EAM_Param

  MD++ fixallatoms constrain_fixedatoms freeallatoms
  MD++ allocchain_parallel initRchain_parallel

  if { $exist > 0 } {
    MD++ incnfile = stringrelax.chain.cn readRchain_parallel
  }

  MD++ eval
  #MD++ stringrelax_parallel
  MD++ stringrelax_parallel_2
  MD++ finalcnfile = stringrelax.chain.cn writeRchain_parallel

  if { 0 } {
  exec cp stringeng.out stringeng_step1.out

  MD++ nebspec = \[ 0 1 0 1 0 1 \] totalsteps = 200000 equilsteps = 200000
  MD++ Broadcast_nebsetting
  MD++ stringrelax_parallel_1
  #we adopt William's revised formulation
  MD++ finalcnfile = stringrelax.chain.cn writeRchain_parallel
  exec cp stringeng.out stringeng_step2.out
  }

  if { $ncpu == 1 }   {
        MD++ quit
  } else {
        MD++ quit_all
  }
  
} elseif { $status == 5 } {
  # anneal surface with MD
  MD++ setnolog
  set ncpu [ MD++_Get numDomains ]
  if { $ncpu < 1 } { set ncpu 1 }
  puts "ncpu = $ncpu"

  initmd $flag-$n.$ncpu

  set epsilon $n

  readpot

  MD++ incnfile = "../w-disl-nuc-hetero-0/0K_${epsilon}_relaxed_surf${flag}.cn" readcn  saveH
  MD++ eval

  setup_window
  MD++ plot_color_windows = 4
#  openwindow

  setup_md
  MD++  T_OBJ = 100   initvelocity  totalsteps = 1000 fixbox = 1  run

  MD++ eval
  relax_fixbox 

  MD++ sleep

  exitmd 
} elseif { $status == 10 } {
  # visualization (from stringrelax serial runs)
  MD++ setnolog
  initmd "view" 
  readpot

#  set epsilon 0.09
  set epsilon $n

  MD++ incnfile = "../w-disl-nuc-hetero-0/0K_${epsilon}_relaxed_surf${flag}.cn" readcn  saveH
  MD++ incnfile = "../w-disl-nuc-hetero-0/w-loop-compression-surf${flag}.cn" readcn restoreH SHtoR setconfig2
  MD++ eval
  MD++ incnfile = "../w-disl-nuc-hetero-0/0K_${epsilon}_relaxed_surf${flag}.cn" readcn  setconfig1
  MD++ eval 

  if { 0 } {
    setup_window
    MD++ NCS = 8 eval calcentralsymmetry 
    MD++ incnfile = "../w-disl-nuc-hetero0/w-loop-compression-surf${flag}.cn" readcn   restoreH SHtoR 
    MD++ finalcnfile = "w-loop-compression-surf${flag}.cfg" writeatomeyecfg
    openwindow
    MD++ sleep
    exitmd
  }

  MD++ zipfiles = 0
  for { set iter 0 } { $iter <= 25 } { incr iter 1 } {
     MD++ incnfile ="../w-disl-nuc-hetero-${flag}-${epsilon}/stringrelax.chain.cn.iter0" readRchain
     MD++ finalcnfile = stringrelax.chain.cn.iter0 writeRchain

     # set chain number
     set chain_no $iter
     set total_no 25

     MD++ input = $iter  copyRchaintoCN eval

     if { $iter == 0 } {
         openwindow
     }
     MD++ eval plot 
     # write cn or cfg file
     MD++ finalcnfile = "chain_no_${chain_no}.cn"  writecn
     MD++ finalcnfile = "chain_no_${chain_no}.cfg" writeatomeyecfg
     MD++ sleepseconds = 1  sleep

  }
  MD++ sleepseconds = 100  sleep
  exitmd

}  elseif { $status == 20 } {
  # visualization (from stringrelax_parallel runs)
  MD++ setnolog
  initmd "view" 
  readpot

  #set epsilon 0.084
  set epsilon $n
  set chain_no 0
  #set chain_no 8
  set total_no 23
  MD++ chainlength = $total_no

  MD++ incnfile = "../w-disl-nuc-hetero-0/0K_${epsilon}_relaxed_surf${flag}.cn" readcn  saveH
  MD++ incnfile = "../w-disl-nuc-hetero-${flag}-${epsilon}.24/stringrelax.chain.cn" 
  MD++ input = $chain_no  readRchain_parallel_toCN  RHtoS

  #MD++ clearR0 refreshnnlist eval quit

#  setup_window

  MD++ zipfiles = 0
  for { set iter 0 } { $iter <= $total_no } { incr iter 1 } {
     # set chain number
     set chain_no $iter
     MD++ input = $chain_no  readRchain_parallel_toCN  RHtoS
     MD++ NCS = 8 eval calcentralsymmetry 

     if { $iter == 0 } {
#         openwindow
     }
     MD++ eval plot 

     if { $iter == 1} {
     #    relax_fixbox
     #    MD++ eval plot sleepseconds = 100 sleep
     }
     # write cn or cfg file
     MD++ finalcnfile = "chain_no_${chain_no}.cn"  writecn
     MD++ finalcnfile = "chain_no_${chain_no}.cfg" writeatomeyecfg
     MD++ sleepseconds = 1  sleep

  }
  MD++ sleepseconds = 100  sleep
  exitmd

} elseif { $status == 21 } {

  MD++ setnolog
  initmd "relax" 

  if { $flag == 112 } {
    #set epsilon 0.04 
    #set epsilon 0.074
    set FEVALMAX 20
  } elseif { $flag == 110 } {
    #set epsilon 0.04
    #set epsilon 0.074
    set FEVALMAX 20
  } else {
     puts "unknown flag = $flag, must be 112 or 110"
     exitmd
  }

  MD++ incnfile = "../w-disl-nuc-hetero-view/chain_no_11.cn" readcn

  setup_window
  openwindow

  MD++ eval
  relax_fixbox

} else {
        
 puts "unknown status = $status"
 exitmd 

} 

