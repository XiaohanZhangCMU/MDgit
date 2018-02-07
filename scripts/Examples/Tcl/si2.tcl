# -*-shell-script-*-

source "scripts/Examples/Tcl/startup.tcl"

#*******************************************
# Definition of procedures
#*******************************************
proc initmd { n } { MD++ {
setnolog
setoverwrite
dirname = runs/si-example-$n
zipfiles = 1               # zip output files
NIC = 200 NNM = 200
#--------------------------------------------
# Create Perfect Lattice Configuration
#
element0 = Si
crystalstructure = diamond-cubic
latticeconst = 5.430949821 #(A) for Si
latticesize = [ 1 0 0 2
                0 1 0 2
                0 0 1 3 ]
} }

#--------------------------------------------
proc openwindow { } { MD++ {
# Plot Configuration
#
atomradius = 0.67 bondradius = 0.3 bondlength = 2.8285 #for Si
atomcolor = orange highlightcolor = purple
bondcolor = red backgroundcolor = gray70
plotfreq = 10  rotateangles = [ 0 0 0 1.25 ]
openwin alloccolors rotate saverot eval plot
} }

#*******************************************
# Main program starts here
#*******************************************
initmd 1
MD++ makecrystal writecn
openwindow
MD++ sleep quit

