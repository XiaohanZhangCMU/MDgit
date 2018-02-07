# -*-shell-script-*-
# TCL script file for Fe simulation using FS potential
#
# Compile code by, e.g.
#  make fs build=R SYS=intel
#
# Run the code by,
#  bin/fs_intel scripts/ME346/fe-md.tcl
#
# The output is stored in the prop.out file whose format is given by
#  output_fmt = "curstep EPOT KATOM Tinst HELM HELMP TSTRESS_xx TSTRESS_yy TSTRESS_zz H_11 H_22 H_33"
#
# The best way to extract the thermal expansion coefficient is to write
#  a separate matlab file to read prop.out and compute the averaged volume at each temperature
#  after the simulation has reached thermal equilibrium

#*******************************************
# Definition of procedures
#*******************************************
proc initmd { n } { 
MD++ setnolog
MD++ setoverwrite
MD++ dirname = "runs/fe-md-$n"
}
#end of proc initmd

proc read_fs { } { 
#-----------------------------------------
# Read the potential file
MD++ potfile = ~/Codes/MD++/potentials/fe_pot readpot
}
#end of proc read_fs

proc create_crystal { } { MD++ {
#--------------------------------------------
# Create Perfect Lattice Configuration
#
latticestructure = body-centered-cubic
latticeconst = 2.8665 #(A) for Fe
latticesize  = [ 1 0 0  5
                 0 1 0  5
                 0 0 1  5 ]
makecrystal finalcnfile = perf.cn writecn
} }
#end of proc create_crystal

#--------------------------------------------
proc setup_md { } { MD++ {     
equilsteps = 0  timestep = 0.0001 # (ps)
atommass = 55.845 # (g/mol)
DOUBLE_T = 0  equilsteps = 0
saveprop = 1 savepropfreq = 100 openpropfile #run
savecn = 1 savecnfreq = 10000 openintercnfile
plotfreq = 20
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
writeall = 1
} }
#end of proc setup_md

#-------------------------------------------------------------
proc openwindow { } { MD++ {
# Plot Configuration
#
atomradius = 0.67 bondradius = 0.3 bondlength = 0
atomcolor = cyan bondcolor = red backgroundcolor = gray70
plotfreq = 100  rotateangles = [ 0 0 0 1.25 ]
openwin  alloccolors rotate saverot eval plot
} }
#end of proc openwindow

#--------------------------------------------
proc exitmd { } { MD++ quit }
#end of proc exitmd

#--------------------------------------------
proc relax_freebox { } { MD++ {
# Conjugate-Gradient relaxation
conj_ftol = 1e-6 conj_itmax = 1000 conj_fevalmax = 1000
conj_fixbox = 0  relax
} }
#end of proc relax_freebox

#--------------------------------------------
# for diagnosis purposes
proc printfile {fname} { puts [read [open $fname r]] }

proc index3 { ID coord } {
    set ind.x 1; set ind.y 2; set ind.z 3
    expr { $ID * 3 + [set ind.$coord] - 1 }
}

proc MD++_GetVector { name ID coord } {
    MD++_Get $name [index3 $ID $coord]
}

proc MD++_PrintVar { name {unit \?} {fmt %11.4e} } {
    puts "$name\t= [format $fmt [MD++_Get  $name]] (in $unit)"
}

proc MD++_PrintMatrix { name {unit \?} {fmt %11.4e} } {
    puts "$name\t= [format $fmt [MD++_Get $name 0]] [format $fmt [MD++_Get $name 1]] [format $fmt [MD++_Get $name 2]]  (in $unit)"
    puts      "\t  [format $fmt [MD++_Get $name 3]] [format $fmt [MD++_Get $name 4]] [format $fmt [MD++_Get $name 5]]"
    puts      "\t  [format $fmt [MD++_Get $name 6]] [format $fmt [MD++_Get $name 7]] [format $fmt [MD++_Get $name 8]]"
}

proc MD++_PrintArray { name {unit \?} {i0 0} {i1 5} {fmt %11.4e} } {
    for {set ID $i0} {$ID < $i1} {incr ID 1} {
        puts "$name\($ID\)= [format $fmt [MD++_Get  $name $ID]] [expr {$ID==$i0?"(in $unit)":" "}]"
    }
}

proc MD++_PrintVectorArray { name {unit \}?} {i0 0} {i1 0} {fmt %11.4e} } {
    for {set ID $i0} {$ID < $i1} {incr ID 1} {
        puts "$name\($ID\)= ([format $fmt [MD++_GetVector $name $ID x]]\
                           [format $fmt [MD++_GetVector $name $ID y]]\
                           [format $fmt [MD++_GetVector $name $ID z]])\
                           [expr {$ID==$i0?"(in $unit)":" "}]"
    }
}
# end of diagnosis procs

    



#*******************************************
# Main program starts here
#*******************************************
if { $argc == 0 } {
 set status 0
} elseif { $argc > 0 } {
 set status [lindex $argv 0]
}
puts "status = $status"



initmd $status
read_fs
create_crystal
setup_md

openwindow

MD++ randseed = 12345 srand48
#MD++  srand48bytime 

MD++ totalsteps = 50000

MD++ T_OBJ = 100 initvelocity  ensemble_type = "NVE" integrator_type = "VVerlet"
    
MD++ run

exitmd

