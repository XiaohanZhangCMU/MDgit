#!/usr/bin/tclsh
set x 1
set y [expr 1.0/(99+$x)]
set z [expr 1/(99+$x)]
puts "y is $y and z is $z when x=$x."
puts "y is [format %20.13e $y]\
      and z is [format %5d $z] when x=$x."
