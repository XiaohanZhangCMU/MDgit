puts -nonewline "Hello! "

#*******************************************
# Definition of procedures
#*******************************************

#--------------------------------------------
proc exitmd { } { MD++ quit }
#end of proc exitmd

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

#--------------------------------------------
# Run the program and arrange to read its input

proc Run {} {
	global command input log but
#	if [catch {open "|$command |& cat"} input] {
#		$log insert end $input\n
#	} else {
#		fileevent $input readable Log
#		$log insert end $command\n
#		$but config -text Stop -command Stop
#	}
        if [catch { set output [eval "$command"]} input] {
		$log insert end $input\n
	} else {
                $log insert end "% $command\n"
		#$but config -text Stop -command Stop
                $log insert end "$output\n"
                $log see end
	}
}

# Read and log output from the program
# (not working yet...  requires thread to enable stopping a running command)
#proc Log {} {
#	global input log
#	if [eof $input] {
#		Stop
#	} else {
#		gets $input line
#		$log insert end $line\n
#		$log see end
#	}
#}
#
# Stop the program and fix up the button
# (not working yet...  requires thread to enable stopping a running command)
#proc Stop {} {
#	global input but
#	#catch {close $input}
#	$but config -text "Run it" -command Run
#}

#--------------------------------------------
# Tk: prepare menu window

if { [info commands button] == "button" } {

#wm title . "MD++"

# Create a frame for buttons and entry.

frame .top -borderwidth 10
pack .top -side top -fill x

# Create the command buttons.

button .top.quit -text Quit -command exit
set but [button .top.run -text "Run it" -command Run]
pack .top.quit .top.run -side right

# Create a labeled entry for the command

label .top.l -text Command: -padx 0
entry .top.cmd -width 20 -relief sunken \
	-textvariable command
pack .top.l -side left
pack .top.cmd -side left -fill x -expand true

# Set up key binding equivalents to the buttons

bind .top.cmd <Return> Run
bind .top.cmd <Control-c> Stop
focus .top.cmd

# Create a text widget to log the output

frame .t
set log [text .t.log -width 60 -height 10 \
	-borderwidth 2 -relief raised -setgrid true \
	-yscrollcommand {.t.scroll set}]
scrollbar .t.scroll -command {.t.log yview}
pack .t.scroll -side right -fill y
pack .t.log -side left -fill both -expand true
pack .t -side top -fill both -expand true

}
#end of Tk: prepare menu window


#---------------------------------------------
proc create_menu_window { } {
# create buttons in Tk widget
   if { [info commands button] == "button" } {
      button .hello   -text Hello -command { puts stdout "Hello, World!" }
      button .goodbye -text Goodbye -command { exit }
      pack   .hello   -padx 10 -pady 10
      pack   .goodbye -padx 10 -pady 10
   } else {
      puts "\a"
      puts "***********************************************************"
      puts "create_buttons command can only be used when Tk is enabled."
      puts "compile with 'make ... Tk=yes'"
      puts "***********************************************************"
   }
}


#----------------------------------------------
puts "startup.tcl loaded."
