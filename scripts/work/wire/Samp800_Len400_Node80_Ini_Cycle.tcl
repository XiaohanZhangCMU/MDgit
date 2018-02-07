# -*-shell-script-*-
# MD code of fibers (coarse grained model of carbon nanotubes)
# computing Young's modulus
#
# Findings:
#  1. lowest energy state for a free-end bendingless chain 
#      is a self-folded structures, folding initiates from the end
#
#
# Ideas:
#  1. introduce random (permanent) kinkness to the fibers
#  2. introduce multiple minima in bending potential 
#  3. finite temperature MD simulations in tension
#

source "scripts/Examples/Tcl/startup.tcl"

### definition of parameters

#*******************************************
# Definition of procedures
#*******************************************
proc initmd { status nfiber orient ncase } { 
MD++ setnolog
MD++ setoverwrite
MD++ dirname = runs/wire/Samp800_Len400_Node80_Ini_Cycle_${nfiber}_${ncase}
# density pref
MD++ initvelocity_type = "Gaussian"
}

proc initialize_potential { } { 
#--------------------------------------------
global bondlen

# in GPa (100-5000)
set Youngs_Modulus 5000 
# in Angstrom
# set CNT_radius     5-20    
set CNT_radius     10    
# CNT_radius out=CNT_radius
# CNT_radius in=CNT_radius out=0.5

set cross_section  [expr 3.1415926 * $CNT_radius * $CNT_radius *0.14] 
set moment_inertia [expr 3.1415926 * $CNT_radius * $CNT_radius * $CNT_radius * $CNT_radius *0.2803 / 4 ] 

# C-C interaction  http://www.csb.yale.edu/userguides/datamanip/autodock/html/Using_AutoDock_305.a.html
# in kcal mol^-1 A^12
set C12_per_atom 2516582.400
# in kcal mol^-1 A^6
set C6_per_atom  1228.800000 
# conversion factor http://mukamel.ps.uci.edu/~zyli/science/units.html
set kcal_per_mol_to_eV 0.04336 

# in Angstrom^2

# 5.2388 for 2 atom
set atom_area     2.6194
set natoms_per_bond [expr 3.1415924 * 2.0 * $CNT_radius * $bondlen] 
set stiffness [expr $Youngs_Modulus / 160.2 * $cross_section / $bondlen]

MD++ C12_00 = [expr $C12_per_atom * $kcal_per_mol_to_eV * $natoms_per_bond ]
MD++ C6_00 =  [expr $C6_per_atom  * $kcal_per_mol_to_eV * $natoms_per_bond ]

# put FEM as a new species and turn off the lj interactions
MD++ C12_03 = 0 C6_03 = 0 
MD++ C12_33 = 0 C6_33 = 0

MD++ BOND_R0 = $bondlen   # Angstrom
MD++ BOND_K  = $stiffness
MD++ BOND_B  = [expr $Youngs_Modulus / 160.2 * $moment_inertia / $bondlen ]     # eV/A^2
MD++ Rcut    = 50   # Angstrom

MD++ substrate_Y = 0
MD++ substrate_REP = $stiffness
MD++ substrate_ATR = $stiffness
MD++ substrate_CON = $stiffness

MD++ initLJ
MD++ NNM = 500
MD++ atommass = \[ 36 36 36 36 \] # mass in g/mol
# MD++ atommass = \[ 12 12 12 12 \] # mass in g/mol

puts "BOND_R0 = [MD++_Get BOND_R0]"
puts "BOND_K  = [MD++_Get BOND_K ]"
puts "BOND_B  = [MD++_Get BOND_B ]"
puts "C12     = [MD++_Get C12_00 ]"
puts "C6      = [MD++_Get C6_00  ]"
puts "Rcut    = [MD++_Get Rcut   ]"

after 20
} 

proc initialize_config { ntype nfiber chainlen bondlen Lx Ly Lz orient_pref } { 
#--------------------------------------------
#Create a fiber structure
#

set flatten 1
MD++ input = \[ $ntype $nfiber $chainlen $bondlen $Lx $Ly $Lz $orient_pref 0 0 $flatten \]
MD++ makefibers 
#finalcnfile = fibers-init.cn writeall = 1 writecn
} 
#end of proc initialize_config

#--------------------------------------------
proc openwindow { } { MD++ {
# Plot Configuration
#
atomradius = 25 bondradius = 5 bondlength = 0
#atomradius = 7.5 bondradius = 7.5 bondlength = 0
atomcolor0 = orange atomcolor1 = green atomcolor2 = red atomcolor3 = black
bondcolor = blue backgroundcolor = white fixatomcolor = yellow
plotfreq = 50  rotateangles = [ 0 90 0 1.2 ] 
plot_map_pbc = 1
win_width = 800 win_height = 800
openwin  alloccolors rotate saverot eval plot
} }
#end of proc openwindow

#--------------------------------------------
proc exitmd { } { MD++ quit }
#end of proc exitmd

#--------------------------------------------
proc relax_freebox { } { MD++ {
# Conjugate-Gradient relaxation
conj_ftol = 1e-4 conj_itmax = 1 conj_fevalmax = 1
# originally 1000000
conj_fixbox = 0  conj_fixboxvec = [ 0 1 1
                                    1 0 1
                                    1 1 0 ]
relax
} }
#end of proc relax_freebox

#--------------------------------------------
proc relax_fixbox { } { MD++ {
# Conjugate-Gradient relaxation
conj_ftol = 1e-4 conj_itmax = 100 conj_fevalmax = 100
# originally 1000000
conj_fixbox = 1  conj_fixboxvec = [ 0 1 1
                                    1 0 1
                                    1 1 0 ]
relax
} }
#end of proc relax_freebox
#--------------------------------------------
proc setup_md { } { MD++ {     
equilsteps = 0  timestep = 0.001
DOUBLE_T = 0
saveprop = 1 savepropfreq = 100    openpropfile 
savecn   = 1   savecnfreq = 10000  openintercnfile
savecfg  = 1  savecfgfreq = 10000  intercfgfile = "inter.cfg"
vt2 = 2e28  
wallmass = 2e3     # atommass * NP = 14380
boxdamp = 1e-3     # optimal damping for 216 atoms and wallmass 1e-3
saveH # Use current H as reference (H0), needed for specifying stress
fixboxvec  = [ 0 1 1 
               1 0 1
               1 1 0 ]
stress = [ 0 0 0
           0 0 0
           0 0 0 ]
output_fmt = "curstep EPOT KATOM Tinst HELM HELMP TSTRESSinMPa_xx TSTRESSinMPa_yy TSTRESSinMPa_zz"
} }
#end of proc setup_md
#--------------------------------------------
proc read_stress {fid} { 
set stress1 [MD++_Get "TSTRESSinMPa_xx"]
set stress2 [MD++_Get "TSTRESSinMPa_yy"]
set stress3 [MD++_Get "TSTRESSinMPa_zz"]
puts $fid "$stress1 $stress2 $stress3"
flush $fid
}
#end of proc read_stress
#--------------------------------------------
proc read_stress_pot {fid} { 
set stress1 [MD++_Get "TSTRESSinMPa_xx"]
set stress2 [MD++_Get "TSTRESSinMPa_yy"]
set stress3 [MD++_Get "TSTRESSinMPa_zz"]
set epot [MD++_Get "EPOT"]
set helm [MD++_Get "HELM"]
set helmp [MD++_Get "HELMP"]
puts $fid "$stress1 $stress2 $stress3 $epot $helm $helmp"
flush $fid
}
#end of proc read_stress_pot
#--------------------------------------------
proc read_stress_H {fid} { 
set stress1 [MD++_Get "TSTRESSinMPa_xx"]
set stress2 [MD++_Get "TSTRESSinMPa_yy"]
set stress3 [MD++_Get "TSTRESSinMPa_zz"]
set H1 [MD++_Get "H_11"]
set H3 [MD++_Get "H_33"]
puts $fid "$stress1 $stress2 $stress3 $H1 $H3"
flush $fid
}
#end of proc read_stress_H
#--------------------------------------------
proc uniaxial_stretch { stretch Lx0 Lz0 } { 
MD++ H_11=[expr $Lx0*$stretch]
MD++ H_33=[expr $Lz0/sqrt($stretch)]

#MD++ input= \[1 1 $strain_rate \] changeH_keepS
#MD++ input= \[3 3 [expr -$strain_rate/2] \] changeH_keepS
}
#end of proc uniaxial_stretch
#--------------------------------------------
proc comp {stretch  Ly0 } {
MD++ H_22=[expr $Ly0*$stretch]
#MD++ H_11=[expr $Ly0/sqrt($stretch)]
#MD++ H_33=[expr $Ly0/sqrt($stretch)]
 
#MD++ input= \[2 2 -$strain_rate \] changeH_keepS
#MD++ input= \[1 1 [expr $strain_rate/2] \] changeH_keepS
#MD++ input= \[3 3 [expr $strain_rate/2] \] changeH_keepS
}
#end of proc comp
#--------------------------------------------
proc reshape {strainx strainz Lx Lz} {
MD++ H_11=[expr $Lx*(1+$strainx)]
MD++ H_33=[expr $Lz*(1+$strainz)]
}
#end of proc reshape




# This process is added by Xiaohan, to perform fem on substrate
# proc initialize_membrane
proc initialize_substrate { } {

# Read in FEM infos and put on contraints
# finalcnfile needs to be moved to c code
#MD++ finalcnfile = "./config1.cn"
set NP0 [MD++_Get "NP"] 
#MD++ incnfile = "./config12.cn" read2cn

if {1} {

MD++ incnfile = "../../../scripts/work/wire/membrane2d-8x5element.cn" read2cn
MD++ map23 = \[ 0 2 \]
MD++ fem_coeff_file = "../../../scripts/work/wire/membrane2d-8x5element-CPE4-Q4.dat"  read_fem_coeff
MD++ elements_file = "../../../scripts/work/wire/membrane2d-8x5element.ele"  read_elements
MD++ fem_bdy_nodes_file = "../../../scripts/work/wire/membrane2d-8x5element.bdy" read_bdy_nodes
MD++ input = $NP0 shift_fem_node_id 
} else {
MD++ incnfile = "../../../scripts/work/wire/membrane2d-5element.cn" read2cn
MD++ map23 = \[ 0 2 \]
MD++ fem_coeff_file = "../../../scripts/work/wire/membrane2d-5element-CPE4-Q4.dat"  read_fem_coeff
MD++ elements_file = "../../../scripts/work/wire/membrane2d-5element.ele"  read_elements
MD++ fem_bdy_nodes_file = "../../../scripts/work/wire/membrane2d-5element.bdy" 

MD++ read_bdy_nodes
MD++ input = $NP0 shift_fem_node_id 
}

MD++ RtoRref    

set n_fem_nodes [MD++_Get n_fem_nodes]
set n_bdy_nodes [MD++_Get n_bdy_nodes]
set n_bdy_tags [MD++_Get n_bdy_tags]
puts "n_bdy_nodes = $n_bdy_nodes"
puts "n_bdy_tags = $n_bdy_tags"
puts "n_fem_nodes = $n_fem_nodes"

for { set niter 0 } { $niter < $n_fem_nodes } { incr niter 1 } {
    set bd [expr $niter + $NP0 ]
    puts "niter = $niter, bd=$bd is set to be species(3)"
#    MD++ fixed($bd) = -1	    
    MD++ species($bd) = 3
}

if { 0 } {
for { set niter 0 } { $niter < $n_bdy_nodes } { incr niter 1 } {
    set bnd [ MD++_Get bNds($niter) ]
    set bid [ MD++_Get bTags($niter) ]
    set xb  [ MD++_Get bXs($niter) ]
    set yb  [ MD++_Get bYs($niter) ]
    set zb  [ MD++_Get bZs($niter) ]
#    puts "gets $bnd is fixed"
    if { $bid == 0 || $bid == 2 || $bid == 3 || $bid == 1 } {
       puts "$bnd is fixed"
       MD++ fixed($bnd) = 1	    
    }
}
}

#MD++ dockbeads

}




#*******************************************
# Main program starts here
#*******************************************
# status 0:Equilibrate, 1:Expand box, 10:Apply attrac,
#        20:Apply wall, 30:Turn on alkanes, 40:Remove attrac
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
 set n2 0
} elseif { $argc > 2 } {
 set n2 [lindex $argv 2]
}
puts "n2 = $n2"

if { $argc <= 3 } {
 set n3 0
} elseif { $argc > 3 } {
 set n3 [lindex $argv 3]
}
puts "n3 = $n3"

# fiber geometry
set ntype      1

# we want to simulate nfiber CNTs
# fibers4_d30p10_2
#set nfiber    [expr 102*1.4 ]
#set nfiber    [expr 128*3.0 ]
set nfiber $n

# set the degree of orientation preference along y axis
#set orient_pref 1.1
#set orient_pref 1.0
#set orient_pref 0.6
set orient_pref $n2

set ncase $n3

#puts "nfiber=$nfiber"
#puts "orient_pref=$orient_pref"
#puts "ncase=$ncase"

# we want to discretize each CNT by 100 segments
set chainlen  20
set L	8000
set Lx $L
set Ly 1000
set Lz $L

# each CNT has length of 1 micron
set CNT_Length 2000
#800  

# find how long is each chain
set bondlen [expr $CNT_Length / $chainlen ]



puts "bondlen = $bondlen"
puts "Lx = $Lx"
puts "Ly = $Ly"
puts "Lz = $Lz"

set density [expr 1.0 * $nfiber * $CNT_Length / ($Lx * $Ly * $Lz) ]
puts "number density = $density per Angstorm^2"

after 20

#*******************************************
if { $status == 0 } {
  # Step 0: creates random structure and relax -- time consuming
  initmd $status $nfiber $orient_pref $ncase
  initialize_potential
  #only leaving the substrate_CON to be non-zero
  MD++ substrate_REP = 0 substrate_ATR = 0
  MD++ randseed = 33333 srand48

# These values are added by Xiaohan, to perform fem on substrate
# all atoms on thread are anchord. cell is anchored as well which has two beads
  MD++ NANCHOR =  [expr $chainlen + 2]
  MD++ K_ANCHOR = 1e-5
  MD++ nnode_per_element = 4


  initialize_config  $ntype $nfiber $chainlen $bondlen $Lx $Ly $Lz $orient_pref
  

# This process is added by Xiaohan, to perform fem on substrate
  initialize_substrate
  MD++ input = \[ $ntype $nfiber $chainlen $bondlen $Lx $Ly $Lz \] linkfibers writeall = 100
  openwindow
  # FIXBOX
  MD++ sleep
  for {set x 1} {$x<=5} {incr x} {
    relax_fixbox
  }
  MD++ eval
  MD++ finalcnfile = "fibers-relaxed-fixbox.cn" writecn

  MD++ freeallatoms
  MD++ input = \[ 1 0 \] fixatoms_by_ID
  MD++ input = 1 setfixedatomsgroup
  MD++ freeallatoms
  
  set cellid1 $chainlen
  set cellid2 [expr $chainlen+1]
  MD++ input = \[ 2 $cellid1 $cellid2 \] fixatoms_by_ID


if { 1 } {
set n_bdy_nodes [MD++_Get n_bdy_nodes]
for { set niter 0 } { $niter < $n_bdy_nodes } { incr niter 1 } {
    set bnd [ MD++_Get bNds($niter) ]
    set bid [ MD++_Get bTags($niter) ]
#    puts "gets $bnd is fixed"
    if { $bid == 0 || $bid == 2 || $bid == 3 || $bid == 1 } {
       puts "$bnd is fixed"
       MD++ fixed($bnd) = 1	    
    }
}
}

  set maxiter 100
  set strainRate 0.001
  for { set iter 1 } { $iter <= $maxiter } { incr iter 1} {
      set f1x 0
      set f1y 0
#      set f1z -0.021
      set f1z [ expr -1*$strainRate*$iter ]

      MD++ extforce = \[ 1 $f1x $f1y $f1z \]
      MD++ conj_itmax = 1000 conj_fevalmax = 1000
      MD++ conj_fixbox = 1 
      MD++ relax
      MD++ finalcnfile = "fibers-zload-$iter.cn" writecn
  }
  MD++ sleep
  exitmd

#*******************************************
} elseif { $status == 100 } {
  # Visualization
  initmd $status $nfiber $orient_pref $ncase

  initialize_potential
  MD++ substrate_REP = 0 substrate_ATR = 0
  MD++ randseed = 33333 srand48

# These values are added by Xiaohan, to perform fem on substrate
# all atoms on thread are anchord. cell is anchored as well which has two beads
  MD++ NANCHOR = [expr $chainlen+2]
  MD++ K_ANCHOR = 1e-1
  MD++ nnode_per_element = 4


  initialize_config  $ntype $nfiber $chainlen $bondlen $Lx $Ly $Lz $orient_pref
  

# This process is added by Xiaohan, to perform fem on substrate
#  initialize_substrate
#  MD++ input = \[ $ntype $nfiber $chainlen $bondlen $Lx $Ly $Lz \] linkfibers writeall = 100

  MD++ eval
  openwindow
 
  set maxiter 100
  set strainRate 0.001
  for { set iter 1 } { $iter <= $maxiter } { incr iter 1} {
      MD++ incnfile = "fibers-zload-$iter.cn" readcn
      MD++ plot
      MD++ eval
      if { $iter == 19 } {
#	 MD++ sleep 
      }
  }
  MD++ sleep 
  exitmd

#*******************************************
} else {

  puts "unknown status = $status"
  exitmd 

} 
