# -*-shell-script-*-
#*******************************************
# Definition of procedures
#*******************************************
proc initmd { } { MD++ {
setnolog
setoverwrite
dirname = runs/si-example
zipfiles = 1   # zip output files
#
NIC = 200 NNM = 200
#--------------------------------------------
# Create Perfect Lattice Configuration
#
element0 = Si
crystalstructure = diamond-cubic
latticeconst = 5.430949821 #(A) for Si
latticesize = [ 1 0 0 4
                0 1 0 4
                0 0 1 4 ]
} }

proc MD++_PrintVar { name {unit \?} {fmt %11.4e} } {
    puts "$name\t= [format $fmt [MD++_Get  $name]] (in $unit)"
}

#*******************************************
# Main program starts here
#*******************************************
initmd
MD++ makecrystal writecn

MD++ saveH

for {set x 997} {$x <= 1003} {incr x} {
  set y [expr $x/1000.0]
  MD++ restoreH input = \[ $y \] scaleH
  MD++ eval
  MD++_PrintVar OMEGA "A^3" "%21.14e"
  MD++_PrintVar EPOT "eV" "%21.14e"
}

MD++ quit

