# -*-shell-script-*-
#*******************************************
# Definition of procedures
#*******************************************
proc initmd { } { MD++ {
setnolog
setoverwrite
dirname = runs/si-example
zipfiles = 1   # zip output files
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

#--------------------------------------------
proc relax_fixbox { } { MD++ {
# Conjugate-Gradient relaxation
conj_ftol = 1e-6 conj_itmax = 1000 conj_fevalmax = 1000
conj_fixbox = 1  relax
} }

#*******************************************
# Main program starts here
#*******************************************
initmd
MD++ makecrystal
if {$argc > 0} {
  set do_relax [ lindex $argv 0 ]
  puts "do_relax = $do_relax"
}

MD++ saveH
set fileID [open "EvsVol.dat" w]
for {set x 997} {$x <= 1003} {incr x} {
  set y [expr $x/1000.0]
  MD++ restoreH input = \[ $y \] scaleH
  if { $do_relax == 1 } { relax_fixbox }
  MD++ eval
  puts $fileID "[format %21.14e [MD++_Get OMEGA]]\t \
                [format %21.14e [MD++_Get EPOT]]"
}
close $fileID
MD++ quit

