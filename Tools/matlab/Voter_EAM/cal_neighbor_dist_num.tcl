# -*-shell-script-*-
# MD++ implementation of Mishin Cu

source "scripts/Examples/Tcl/startup.tcl"

#*******************************************
# Definition of procedures
#*******************************************
proc initmd { } { MD++ {
#  setnolog
  setoverwrite
  dirname = runs/cu-bulk
  zipfiles = 1
  NNM = 300 # number of neighbors per atom (default: 60)
}}

#-----------------------------------------------------------
proc readeam { } { MD++ {
#Read in EAM Mishin Cu potential, PRB (2001) 63 224106        
potfile = "~/Codes/MD++/potentials/EAMDATA/eamdata.CuMishin" eamgrid = 5000  readeam
}}

#-----------------------------------------------------------
proc makecrystal { } { MD++ {
#Create Perfect Lattice Configuration
#
  latticestructure = face-centered-cubic
  element0 = Cu latticeconst = 3.6149 #(A) for Cu
  latticesize = [  1 0 0   10    #(x)
                   0 1 0   10    #(y)
                   0 0 1   10  ] #(z)            
  makecn 
}}

#-----------------------------------------------------------
proc relax_fixbox_setup { } { MD++ {
# Conjugate-Gradient relaxation
conj_ftol = 1e-10 conj_itmax = 1000 conj_fevalmax = 1000
conj_fixbox = 1             
} }
#end of proc relax_fixbox_setup

#-----------------------------------------------------------
proc relax_freebox_setup { } { MD++ {
# Conjugate-Gradient relaxation
conj_ftol = 1e-10 conj_itmax = 1000 conj_fevalmax = 1000
conj_fixbox = 0             #allow box to change shape
conj_fixboxvec = [ 0 0 1    #fix three components of the 
                   1 0 1    # H matrix to prevent rotation
                   0 0 0 ]
} }
#end of proc relax_freebox_setup


#-----------------------------------------------------------
proc exitmd { } { MD++ quit }
#end of proc exitmd


#*******************************************
# Main program starts here
#*******************************************
# status 0: 
#        1: 
#        2: Visualize
#
# read in status from command line argument
if { $argc == 0 } {
 set status 0
} elseif { $argc > 0 } {
 set status [lindex $argv 0]
}
puts "status = $status"

if { $status == 0 } {
# 
  initmd
  readeam

  makecrystal  
  MD++ eval

  set NP [MD++_Get NP]
  set Lx [MD++_Get H_11];
  set Ly [MD++_Get H_22];
  set Lz [MD++_Get H_33];
  set sx0 [MD++_GetVector SR 0 x]
  set sy0 [MD++_GetVector SR 0 y]
  set sz0 [MD++_GetVector SR 0 z]
  set fid [open "interatomic_dist.fcc" w]
  for { set i 1 } { $i < $NP } { incr i } {

    set sx [MD++_GetVector SR $i x]
    set sy [MD++_GetVector SR $i y]
    set sz [MD++_GetVector SR $i z]

    set dsx [expr $sx - $sx0]; set dsx [expr $dsx - round($dsx)]
    set dsy [expr $sy - $sy0]; set dsy [expr $dsy - round($dsy)]
    set dsz [expr $sz - $sz0]; set dsz [expr $dsz - round($dsz)]
 
    set dx [expr $Lx * $dsx]
    set dy [expr $Ly * $dsy]
    set dz [expr $Lz * $dsz]

    set dr [expr sqrt($dx*$dx + $dy*$dy + $dz*$dz)]
    set dr [expr $dr/3.6149]
    puts $fid "$dr"
  }

  MD++ quit

} elseif { $status == 1 } {

} else {
        
  puts "unknown status = $status"
  exitmd 

} 
