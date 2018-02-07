#!/usr/bin/wish

if { [info commands button] == "button" } {

button .hello   -text Hello -command { puts stdout "Hello, World!" }
button .goodbye -text Goodbye -command { exit }
pack   .hello   -padx 10 -pady 10
pack   .goodbye -padx 10 -pady 30

} else {

      puts "\a"
      puts "***********************************************************"
      puts "create_buttons command can only be used when Tk is enabled."
      puts "compile with 'make ... Tk=yes'"
      puts "***********************************************************"

}

