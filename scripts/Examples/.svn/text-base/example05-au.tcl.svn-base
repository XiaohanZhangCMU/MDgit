# -*-shell-script-*-
# test the MEAM potential of Au (TCL)
# run by meam-baskes
#------------------------------------------------------------

source "scripts/Examples/Tcl/startup.tcl"

#*******************************************
# Definition of procedures
#*******************************************
proc initmd { } { MD++ {
setnolog
setoverwrite
#dirname = runs/au-test-s0-md-300K
dirname = runs/au-test
zipfiles = 0 writeall = 1
} }
#end of proc initmd

#------------------------------------------------------------
proc readmeam-baskes { } { MD++ {
#Read in MEAM potential (Baskes format)
meamfile = "~/Codes/MD++/Fortran/MEAM-Baskes/meamf"
meafile  = "~/Codes/MD++/Fortran/MEAM-Baskes/meafile_Au_orig"
ntypes = 1  ename0 = "Au1  " #leave spaces to be compatible with fortran

rcut = 6.0  kode0 = "library"
readMEAM
NNM = 300
} }

#------------------------------------------------------------
proc readmeam-lammps { } { MD++ {
#Read in MEAM potential (Baskes format)
meamfile = "~/Codes/MD++/Fortran/MEAM-Lammps/meamf"
nspecies = 1  element0 = "Au"  rcut = 6.0  readMEAM
NNM = 300
} }

#-------------------------------------------------------------
proc readmeam { } { MD++ {
#Read in MEAM potential (MD++ format)
potfile = "~/Codes/MD++/potentials/MEAMDATA/meamdata.Au2" readMEAM
#potfile = "~/Codes/MD++/potentials/MEAMDATA/meamdata.Siz" readMEAM
#potfile = "~/Codes/MD++/potentials/MEAMDATA/meamdata.Si" readMEAM
#xzbl = -10 xzbl0 = -9 #turn off correction
} }

#-------------------------------------------------------------
proc readeam { } { MD++ {
#Read in EAM potential (MD++ format)
potfile = "~/Codes/MD++/potentials/EAMDATA/eamdata.AuFoiles" eamgrid = 500 readeam
NNM = 300
} }

#--------------------------------------------
proc create_FCC_Au { } { MD++ {
# Create Perfect Lattice Configuration
crystalstructure = face-centered-cubic  latticeconst = 4.07 #(A) FCC  structure
element0 = "Au"   element1 = "Si"  
atommass = [ 196.96655 28.0855 ] # Au, Si (g/mol)
latticesize = [ 1 0 0  4
                0 1 0  4
                0 0 1  4 ]
makecrystal finalcnfile = fcc-au-4x4x4.cn writecn
incnfile = fcc-au-4x4x4.cn writecn
} }        

#--------------------------------------------
proc create_BCC_Au { } { MD++ {
# Create Perfect Lattice Configuration
crystalstructure = body-centered-cubic  latticeconst = 3.0 #(A) BCC  structure
element0 = "Au"   element1 = "Si"  
atommass = [ 196.96655 28.0855 ] # Au, Si (g/mol)
latticesize = [ 1 0 0  5
                0 1 0  5
                0 0 1  5 ]
makecrystal finalcnfile = bcc-au-4x4x4.cn writecn
incnfile = bcc-au-4x4x4.cn writecn
} }        

#--------------------------------------------
proc create_HCP_Au { } { MD++ {
# Create Perfect Lattice Configuration
crystalstructure = hexagonal-ortho  latticeconst = [ 2.87792 2.87792 4.6996 ] #(A) HCP  structure
element0 = "Au"   element1 = "Si"  
atommass = [ 196.96655 28.0855 ] # Au, Si (g/mol)
latticesize = [ 1 0 0  5
                0 1 0  5
                0 0 1  5 ]
makecrystal finalcnfile = hcp-au-5x5x5.cn writecn
incnfile = hcp-au-5x5x5.cn writecn
} }        

#-------------------------------------------------------------
proc openwindow { } { MD++ {
# Plot Configuration
#
atomradius = [ 0.67 0.34 ] bondradius = 0.3 bondlength = 1.0
atomcolor = orange bondcolor = red backgroundcolor = gray70 fixatomcolor = yellow
plotfreq = 100  rotateangles = [ 0 0 0 1.25 ]
color00 = "orange"  color01 = "red"    color02 = "green"
color03 = "magenta" color04 = "cyan"   color05 = "purple"
color06 = "gray80"  color07 = "white"  color08 = "blue"
#
#plot_color_axis = 2  #2: use CSD (default 0: use local energy)
#
plot_color_windows = [ 0
                      1.5  3  8 #color08 = blue
                       3  20  6 #color06 = gray80
		       20 40  7
                       40 100 1
                      ]
#plot_limits = [ 1 -10 10 -0.05 10 -10 10 ]
plot_atom_info = 1  plotfreq = 10 
#
openwin  alloccolors rotate saverot eval plot
} }
#end of proc openwindow

#--------------------------------------------
proc exitmd { } { MD++ quit }
#end of proc exitmd

#--------------------------------------------
proc relax_freebox { } { MD++ {
# Conjugate-Gradient relaxation
conj_ftol = 2e-6 conj_itmax = 1000 conj_fevalmax = 10000
conj_fixbox = 0  conj_fixboxvec = [ 0 1 1
                                    1 0 1
                                    1 1 0 ]
relax
} }
#end of proc relax_freebox
#--------------------------------------------
proc relax_fixbox { } { MD++ {
# Conjugate-Gradient relaxation
conj_ftol = 2e-6 conj_itmax = 1000 conj_fevalmax = 10000
conj_fixbox = 1  conj_fixboxvec = [ 0 1 1
                                    1 0 1
                                    1 1 0 ]
relax
} }
#end of proc relax_freebox

#--------------------------------------------
proc setup_md { } { MD++ {     
equilsteps = 0  timestep = 0.001 # (ps)
atommass = [ 196.96655 28.0855 ] # Au, Si (g/mol)
DOUBLE_T = 0
saveprop = 1 savepropfreq = 100 openpropfile #run
savecn = 1 savecnfreq = 1000 openintercnfile
plotfreq = 100
vt2 = 2e28  #1e28 2e28 5e28
wallmass = 2e3     # atommass * NP = 14380
boxdamp = 1e-3     # optimal damping for 216 atoms and wallmass 1e-3
saveH # Use current H as reference (H0), needed for specifying stress
fixboxvec  = [ 0 1 1 
               1 0 1
               1 1 0 ]
stress = [ 0 0 0
           0 0 0
           0 0 0 ]
output_fmt = "curstep EPOT KATOM Tinst HELM HELMP TSTRESS_xx TSTRESS_yy TSTRESS_zz H_11 H_22 H_33" 
} }
#end of proc setup_md



#---------------------------------------------
# Compute elastic constants for cubic crystal
proc compute_elast_const { datafile structname argv0 } {    
relax_freebox
MD++ finalcnfile = relaxed.cn writecn

set nx [MD++_Get latticesize  3]
set ny [MD++_Get latticesize  7]
set nz [MD++_Get latticesize 11]

set Lx [MD++_Get H 0]
set Ly [MD++_Get H 4]
set Lz [MD++_Get H 8]

set NP [MD++_Get NP]

set a0 [expr $Lx/$nx]

set nx0 [expr $Lx/$a0]
set ny0 [expr $Ly/$a0]
set nz0 [expr $Lz/$a0]

#-----------------------
# Compute bulk modulus
#-----------------------
MD++ saveH conj_fixbox = 1

set fileID [open $datafile w ]

puts $fileID "%This file is generated by MD++ ($argv0)\n\n"
puts $fileID "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
puts $fileID "% $structname "
puts $fileID "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
puts $fileID "$structname = struct();"
puts $fileID "$structname.N = $NP;"
puts $fileID "$structname.a0 = $a0;"
puts $fileID "$structname.C1 = \[ $nx0 0 0 \];"
puts $fileID "$structname.C2 = \[ 0 $ny0 0 \];"
puts $fileID "$structname.C3 = \[ 0 0 $nz0 \];"

puts $fileID "$structname.data = \[ %lattice_const.(A)  Epot(eV)"

for {set x 0.997} {$x<=1.003} {set x [expr $x+0.001]} {

    MD++ restoreH input = $x scaleH relax eval

    puts $fileID "[format %.3f%3d%24.13e $x 1 [MD++_Get EPOT]]"
}

puts $fileID "];"
puts $fileID "$structname.data(:,1) = $structname.data(:,1)*$structname.a0;\n%"
puts $fileID "$structname = fit_aEB($structname);\n%"
puts $fileID "%disp(sprintf('$structname: a0 = %.6f A, Ecoh = %.6f eV, B = %.6f GPa', ..."
puts $fileID "%             $structname.a0, $structname.Ecoh, $structname.B));\n"

#-----------------------
# Compute C11 C12
#-----------------------
puts $fileID "$structname.data = \[ %lattice_const.(A)  Epot(eV)"

for {set x -0.003} {$x<=0.003} {set x [expr $x+0.001]} {

    MD++ restoreH input = \[ 1 1 $x \] shiftbox relax eval

    puts $fileID "[format %.3f%3d%24.13e [expr $x+1] 1 [MD++_Get EPOT]]"
}

puts $fileID "];"
puts $fileID "$structname = fit_C11($structname);\n%"
puts $fileID "%disp(sprintf('$structname: B = %.6f GPa, C11 = %.6f GPa, C12 = %.6f GPa', ..."
puts $fileID "%             $structname.B, $structname.C11, $structname.C12));\n"



#-----------------------
# Compute C44
#-----------------------
# rotate crystal by 45 degrees around z axis
MD++  {
    restoreH relax
    finalcnfile = perf.cn writeall = 1 writecn
    incnfile = perf.cn input = [ 1 0 ] splicecn
    input = [ 2 1 0.5 ] redefinepbc
    input = [ 1 2 -1  ] redefinepbc  reorientH 
    relax eval saveH
    finalcnfile = perf2.cn writecn
    incnfile = perf2.cn readcn
}

puts $fileID "$structname.N = [expr $NP*2];"
puts $fileID "$structname.C1 = \[  $nx0 -$ny0 0 \];"
puts $fileID "$structname.C2 = \[  $nx0 $ny0 0 \];"
puts $fileID "$structname.C3 = \[  0 0 $nz0 \];"

puts $fileID "$structname.data = \[ %lattice_const.(A)  Epot(eV)"

MD++ fixallatoms  #Note: only for SC and BCC crystals!

for {set x -0.003} {$x<=0.003} {set x [expr $x+0.001]} {

    MD++ restoreH input = \[ 1 1 $x \] shiftbox relax eval

    puts $fileID "[format %.3f%3d%24.13e [expr $x+1] 1 [MD++_Get EPOT]]"
}

puts $fileID "];"
puts $fileID "$structname = fit_C44($structname);\n%"
puts $fileID "disp(sprintf('$structname: a0 = %.6f A, Ecoh = %.6f eV, B = %.6f GPa', ..."
puts $fileID "             $structname.a0, $structname.Ecoh, $structname.B));"
puts $fileID "disp(sprintf('        C11 = %.6f GPa, C12 = %.6f GPa, C44 = %.6f GPa', ..."
puts $fileID "             $structname.C11, $structname.C12, $structname.C44));\n"


close $fileID
}     



#*******************************************
# Main program starts here
#*******************************************
# read in status from command line argument
if { $argc == 0 } {
 set status 0
} elseif { $argc > 0 } {
 set status [lindex $argv 0]
}
puts "status = $status"

set myname  [ MD++_Get "myname"]
puts "myname = $myname"

initmd

if { [ string match "*baskes*" $myname ] } {
  puts "readmeam-baskes"
  readmeam-baskes
} elseif { [ string match "*lammps*" $myname ] } {
  puts "readmeam-lammps"
  readmeam-lammps
} elseif { [ string match "*meam*" $myname ] } {
  puts "readmeam"
  readmeam
} elseif { [ string match "*eam*" $myname ] } {
  puts "readeam"
  readeam
}

if { $status == 0 } {
  #computing Ecoh
  create_FCC_Au
  #openwindow
  MD++ eval

  set N  [MD++_Get "NP"]
  set E1 [MD++_Get "EPOT"]

  puts "N = $N  E1 = $E1"

  MD++ fixallatoms  
  relax_freebox  

  MD++ eval
  set E0 [MD++_Get "EPOT"]
  set Ecoh [expr $E0/$N]
  set Lx [MD++_Get H_11]
  set Nx [MD++_Get latticesize(3)]
  set a0 [expr $Lx/$Nx]

  puts "Lx = $Lx Nx = $Nx a0 = $a0"
  puts "E0 = $E0  Ecoh = $Ecoh"


  #MD++ sleep

  
} elseif { $status == 1 } {
  #computing unrelaxed vacancy energy
  create_FCC_Au
  MD++ eval
  set E0 [MD++_Get "EPOT"]
  set N  [MD++_Get "NP"]
  set Ecoh [expr $E0/$N]

  MD++ input = \[ 1 0 \] fixatoms_by_ID removefixedatoms
  #openwindow
  MD++ eval
  set N1 [MD++_Get "NP"]
  set E1 [MD++_Get "EPOT"]
  set Ev [expr $E1 - $E0 + $Ecoh]

  puts "N = $N  E0 = $E0 E1 = $E1  Ecoh = $Ecoh  Ev = $Ev"

  relax_fixbox  
  set E2 [MD++_Get "EPOT"]
  set Ev_relaxed [expr $E2 - $E0 + $Ecoh]
  puts "Ev_relaxed = $Ev_relaxed"

  #MD++ sleep

} elseif { $status == 2 } {
  #BCC structure
  create_BCC_Au
  #openwindow
  MD++ eval

  set N  [MD++_Get "NP"]
  set E1 [MD++_Get "EPOT"]

  puts "N = $N  E1 = $E1"

  MD++ fixallatoms  
  relax_freebox  

  MD++ eval
  set E0 [MD++_Get "EPOT"]
  set Ecoh [expr $E0/$N]
  set Lx [MD++_Get H_11]
  set Nx [MD++_Get latticesize(3)]
  set a0 [expr $Lx/$Nx]

  puts "Lx = $Lx Nx = $Nx a0 = $a0"
  puts "E0 = $E0  Ecoh = $Ecoh"


  #MD++ sleep

} elseif { $status == 3 } {
  #HCP structure
  create_HCP_Au
  #openwindow
  #MD++ sleep
  MD++ eval

  set N  [MD++_Get "NP"]
  set E1 [MD++_Get "EPOT"]

  puts "N = $N  E1 = $E1"
  
  MD++ fixallatoms  
  relax_freebox  

  MD++ eval
  set E0 [MD++_Get "EPOT"]
  set Ecoh [expr $E0/$N]
  set Lx [MD++_Get H_11]
  set Nx [MD++_Get latticesize(3)]
  set a0 [expr $Lx/$Nx]
  set Lz [MD++_Get H_33]
  set Nz [MD++_Get latticesize(11)]
  set c0 [expr $Lz/$Nz]
  set covera [expr $c0/$a0]
  puts "Lx = $Lx Nx = $Nx a0 = $a0"
  puts "Lz = $Lz Nz = $Nz c0 = $c0"
  puts "E0 = $E0  Ecoh = $Ecoh  covera = $covera"


  #MD++ sleep

} elseif { $status == 4 } {
  #FCC structure with stacking fault

 MD++ {
  crystalstructure = face-centered-cubic  
  latticeconst = 4.07 #(A) FCC  structure
  element0 = "Au"   element1 = "Si"  
  atommass = [ 196.96655 28.0855 ] # Au, Si (g/mol)
  latticesize = [ 1 -1 0   4
                  1  1 1   6
                 -1 -1 2   4 ]
  makecrystal finalcnfile = fcc-au-4x4x4.cn writecn
  incnfile = fcc-au-4x4x4.cn writecn

  input = [ 1 -10 10 0.412 10   -10 10 ] fixatoms_by_position
  input = [ 1 -10 10 -10 -0.454 -10 10 ] fixatoms_by_position
  removefixedatoms
 }

 MD++ eval
 set E1 [MD++_Get "EPOT"]
 #MD++ plot_color_axis = 2 eval input = \[ 0.1 10 100 \] GnuPlotHistogram

 #creating stacking fault
 MD++ {
  input = [ 1 -10 10 0.112 10   -10 10 ] fixatoms_by_position
  input = 1 setfixedatomsgroup freeallatoms
 }

 set Lx [MD++_Get H_11]
 set Lz [MD++_Get H_33]
 set A  [expr $Lx*$Lz]

 set dz [expr -$Lz/24]
 MD++ input = \[ 1  0 0 $dz  1 \] movegroup
 
  openwindow

  MD++ plot_color_axis = 2 plot_color_windows = 4 

  MD++ eval plot

  set E2 [MD++_Get "EPOT"]
  set Esf [expr ($E2-$E1)/$A]
  set Esf_in_Jm2 [expr $Esf * 1.602e4]
  puts "Esf = $Esf (eV/A) = $Esf_in_Jm2 (mJ/m^2)"

  #MD++ plot_color_axis = 2 eval input = \[ 0.1 10 100 \] GnuPlotHistogram

  MD++ sleep

}  elseif { $status == 5 } {
  #FCC structure with edge dislocation

 MD++ {
  crystalstructure = face-centered-cubic  
  latticeconst = 4.07 #(A) FCC  structure
  element0 = "Au"   element1 = "Si"  
  atommass = [ 196.96655 28.0855 ] # Au, Si (g/mol)
  latticesize = [ 1 -1 0   20
                  1  1 1   12
                 -1 -1 2   2 ]
  makecrystal finalcnfile = fcc-au-4x4x4.cn writecn
  incnfile = fcc-au-4x4x4.cn writecn

  input = [ 1 -10 10 0.412 10   -10 10 ] fixatoms_by_position
  input = [ 1 -10 10 -10 -0.454 -10 10 ] fixatoms_by_position
  removefixedatoms
 }

 openwindow
 MD++ plot_color_axis = 2 plot_color_windows = 4 plot
 #MD++ sleep

 MD++ {
  input = [ 3 2 #z(dislocation line), y(dipole direction)
            0.025 0 0 #(burgers vector relative to box)
            0.00625 -0.5139  -0.0139  
            0.305 #\nu
            -10 10 -10 10 1  #number of images
            0 0 0 0 ]
  makedipole finalcnfile = makedp.cn writecn
 } 
  
 
 #MD++ plot_color_axis = 2 eval input = \[ 0.1 10 100 \] GnuPlotHistogram
 MD++ eval plot
 #MD++ sleep

  relax_fixbox
  MD++ finalcnfile = "relaxed.cn" writecn
  MD++ sleep
}  elseif { $status == 6 } {
  #FCC structure with edge dislocation (dissociated)

 MD++ {
  crystalstructure = face-centered-cubic  
  latticeconst = 4.07 #(A) FCC  structure
  element0 = "Au"   element1 = "Si"  
  atommass = [ 196.96655 28.0855 ] # Au, Si (g/mol)
  latticesize = [ 1 -1 0   20
                  1  1 1   12
                 -1 -1 2   2 ]
  makecrystal finalcnfile = fcc-au-4x4x4.cn writecn
  incnfile = fcc-au-4x4x4.cn writecn

  input = [ 1 -10 10 0.412 10   -10 10 ] fixatoms_by_position
  input = [ 1 -10 10 -10 -0.454 -10 10 ] fixatoms_by_position
  removefixedatoms
 }

 #openwindow
 #MD++ plot_color_axis = 2 plot_color_windows = 4 plot

 #MD++ sleep

 MD++ {
  input = [ 3 2 #z(dislocation line), y(dipole direction)
            0.025 0 0 #(burgers vector relative to box)
            0.00625 -0.5139  -0.0139  
            0.305 #\nu
            -10 10 -10 10 1  #number of images
            0 0 0 1 ]
  makedipole finalcnfile = makedp.cn writecn
 } 

 #bp = ([1 -1 0]/4 + [-1 -1 2]/12 ) = [1 -2 1]/6 
 #bp (in scaled coordinate) = [1/40 0 1/24 ] 
 MD++ {
  input = [ 3 1 #z(dislocation line), x(dipole direction)
            -0.0125 0 -0.0417 #(burgers vector relative to box)
            -0.0139 0.00625 0.35
            0.305 #\nu
            -10 10 -10 10 0  #number of images
            0 0 0 0 ]
  makedipole finalcnfile = makedp-2.cn writecn
 } 
  
 
 #MD++ plot_color_axis = 2 eval input = \[ 0.1 10 100 \] GnuPlotHistogram
 MD++ plot_color_axis = 2 eval plot
 #MD++ sleep

  relax_fixbox
  MD++ finalcnfile = "relaxed.cn" writecn
  #MD++ sleep


  #E (undissociated) = -9.66955548666e3   dE =  0      eV  sep = 0
  #E (dissociated)   = -9.6709018352e3    dE = -1.3463 eV  sep = 0.1
  #E (dissociated)   = -9.6725472859e3    dE = -2.9918 eV  sep = 0.2
  #E (dissociated)   = -9.6725819649e3    dE = -3.0265 eV  sep = 0.3
  #E (dissociated)   = -9.6724053235e3    dE = -2.8498 eV  sep = 0.4
  #E (dissociated)   = -9.6719741529e3    dE = -2.4187 eV  sep = 0.5

}  elseif { $status == 7 } {
  #FCC structure with edge dislocation (dissociated) : MD simulation

  MD++ incnfile = "../au-test/relaxed-2x2.cn" readcn
  #openwindow
  MD++ plot_color_axis = 2 plot_color_windows = 4 eval plot

  setup_md
  MD++ T_OBJ = 300 initvelocity ensemble_type = "NVT"
  MD++ integrator_type = "VVerlet" implementation_type = 2
  MD++ totalsteps = 10000 run
  MD++ finalcnfile = "au_perfedge_equil.cn" writecn
  
  relax_fixbox
  MD++ finalcnfile = "relaxed.cn" writecn
  #MD++ sleep

}  elseif { $status == -1 } {
  #load structure to view

 set inputdir ../au-test-s0-md-300K

 MD++ incnfile = "relaxed-2x2.cn" readcn
 #MD++ incnfile = "relaxed-1x1-partials.cn" readcn
 #MD++ incnfile = "$inputdir/relaxed.cn" readcn
 #MD++ incnfile = "$inputdir/au_perfedge_equil.cn" readcn
 openwindow
  
 #MD++ sleep

 #MD++ plot_color_axis = 2 eval input = \[ 0.1 10 100 \] GnuPlotHistogram
 MD++ plot_color_axis = 2 eval plot

 for { set i 1 } { $i <= 11 } { incr i 1 } {
     set num [format %04d $i]
     MD++ incnfile = "$inputdir/inter$num.cn" readcn
     MD++ eval plot
     MD++ finalcnfile = "$inputdir/inter$num.cfg" writeatomeyecfg
 }
 
 MD++ sleep
}

exitmd

