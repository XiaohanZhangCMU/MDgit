# -*-shell-script-*-

source "scripts/Examples/Tcl/startup.tcl"

MD++ setnolog
MD++ setoverwrite
MD++ dirname = runs/si-example

#------------------------------------------------------------
#Create Perfect Lattice Configuration
#
MD++  element0 = Si
MD++  crystalstructure = diamond-cubic
MD++  latticeconst = 5.4309529817532409 #(A) for Si
MD++  {
  latticesize = [ 1 0 0 2
                  0 1 0 2
                  0 0 1 3 ]
}
MD++  makecrystal  writecn

#------------------------------------------------------------
#Plot Configuration
#
MD++  atomradius = 0.67 bondradius = 0.3 bondlength = 2.8285
MD++  atomcolor = orange highlightcolor = purple
MD++  bondcolor = red backgroundcolor = white
MD++  plotfreq = 10  rotateangles = \[ 0 0 0 1.25 \]
MD++  openwin alloccolors rotate saverot eval plot
MD++  sleep quit

