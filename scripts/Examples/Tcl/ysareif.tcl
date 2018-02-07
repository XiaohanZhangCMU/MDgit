#!/usr/bin/tclsh
puts "The name of this script is $argv0"
if {$argc > 0} {
  puts "There are $argc arguments to this script."
  puts "The argument is : $argv"
  set flag $argv
} else {
  puts "There are no argument to this script."
  puts "Default argument = \"G(radual)\" is used."
  set flag "G"
}

switch $flag {
  G {
    puts "In case of argument = \"G(radual)\","
    for {set x 1} {$x < 100} {incr x} {
      set y [expr 1.0/(99+$x)]
      puts "y = [format %20.13e $y] at x = [format %5d $x]"
    }
  }
  I {
    puts "In case of argument = \"I(nstantaneous)\","
    for {set x 1} {$x < 100} {incr x} {
      set y [expr $x/100.0]
      puts "y = [format %20.13e $y] at x = [format %5d $x]"
    }
  }
  default {
    puts "No action specified for the given argument"
  }
}

