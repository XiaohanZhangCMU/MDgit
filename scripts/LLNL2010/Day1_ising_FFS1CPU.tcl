# -*-shell-script-*-
# Ising 2D Monte Carlo simulation
#-------------------------------------------

set graphic "On" 
set lam_auto "On"
# set number of configurations saved in a interface
set Ncut 50

# stop FFS at interface NMAX
set NMAX 100
source "scripts/Examples/Tcl/startup.tcl"

if { $argc==0} {
# default
	set n 1
	set m 0
        set m 100
        set h 0.06
        set t 1.5 
} elseif { $argc > 0 } {
	set n 1
	set m [ lindex $argv 0 ]
        set m1 [ lindex $argv 1 ]
        set h [ lindex $argv 2 ]
	set t [ lindex $argv 3 ]

set NMAX [ expr $m1-$m]

}

set check [ glob -nocomplain runs/2D_ISING ]
set key [ llength $check ]
if { $key < 1 } {
exec mkdir runs/2D_ISING
}

set check [ glob -nocomplain runs/2D_ISING/FFS_${h}_${t} ]
set key [ llength $check ]
if { $key < 1 } {
exec mkdir runs/2D_ISING/FFS_${h}_${t}
}

#### FFS setup for default lambda_A and lambda_B ####
MD++ lambda_A = 15
MD++ lambda_B = 5000
#### assign interface number
MD++ FFSoption = ${m}

### read when pruning technique is used
MD++ FFSpruning=0
MD++ FFS_Pp=1.0

# for initial brute force MC step totalsteps
set step0totstep 30000 

# for second & further stages
set singlestepslimit 30000
set totalstepslimit 800000

### default interface of lambdas 
set lambdasequence "\[ 44 64 84 104 124 144 164 184 204 224 \]" 

#*******************************************
# Definition of procedures
#*******************************************
proc init_ising { n h t } {
#MD++ setnolog
MD++ setoverwrite
MD++ dirname = "runs/2D_ISING/FFS_${h}_${t}/CPUno_${n}"
MD++ srand48bytime 
MD++ srandbytime
#--------------------------------------------
# Initializtion
#
MD++ NX = 100  NY = 100 input = -1  initspin
#incnfile = ../ising2d/s100.cn readcn
#--------------------------------------------
# Physical parameters
MD++ J = 1     # interaction strength
MD++ H = ${h}  # applied field
MD++ kBT = ${t} # temperature
MD++ zipfiles = 1
} 

#--------------------------------------------
proc openwindow { } { MD++ {
# Plot Configuration
#
atomradius = 1.0 win_width = 480 win_height = 480 rotateangles = [ 0 0 0 1.7]
color00 = white color01 = blue backgroundcolor = gray80
openwin alloccolors rotate saverot plot
} }
#end of proc openwindow

proc run_setup { } {
MD++ {
	saveFFScn= 1
	saveprop = 1 savepropfreq=1 plotfreq = 100
        savecn   = 0 calcryfreq=1 calcryfreq1=250
	openpropfile
	openFFScnfile
}
}	

proc equil_run { } {
   puts "run MC to obtain mean largest cluster size"
   MD++ {
        saveprop = 0 savecn = 0
	plotfreq = 100
        calcryfreq1 = 250
	calcryfreq = 1 
        totalsteps = 10000
        MCrun
   }
        set Nc_avg [ MD++_Get N_lgst_avg ]
   puts "mean largest cluster size = $Nc_avg"
}

proc lam_set { t h n } {
   set Nc_avg [ MD++_Get N_lgst_avg ]
   if { $Nc_avg < 1 } {
   set Nc_avg 1
   }

   set l_A [ expr $Nc_avg*2.5]
   set l_A [ expr int($l_A) ]
   if { $l_A < 2 } {
      set l_A 2
   }

### Assign lambda_A
   MD++ lambda_A = $l_A
  
   set quo [ expr $l_A/4 ]
   if { $quo < 2 } {
   set quo 1
   }

### compute interface from lambda_0 to lambda_200
   set fp [ open "../lam_array.txt" w ]
   puts $fp "\["
   set a(0) $l_A
   puts $fp 200
   for { set i 1 } { $i <= 200 } { incr i 1 } {
      set quo1 [ expr $i/4 ]
      set i_ [ expr $i-1 ]
      set a($i) [ expr $a($i_)+$l_A+$quo*$quo1 ]
#      puts $a($i)
      puts $fp $a($i)
   }
   puts $fp "\]"
   close $fp
}

proc read_lam { } {
   while { 1 } {
	set check [ glob -nocomplain ../lam_array.txt ]
	set NUM [ llength $check ]
	if { $NUM > 0 } {
	    break;
	}  
        exec sleep 10
   }
   set fp [ open "../lam_array.txt" r ]
   set data [ read $fp ]
   close $fp
   return $data
}

proc get_init_cn { m } {

### select configuration file randomly from previous interface
	set mm [ expr $m - 1 ]
	set fp [ open ../total${mm}.txt r ]
	set data [ split [ read $fp ] ]
	close $fp
	
	set fp [ open ../weight${mm}.txt r ]
	set weights [ split [ read $fp] \n]
	close $fp
	set NUM [ llength $weights ]
	set NUM [ expr $NUM-1 ]

	set ran_max [ lindex $data 2 ]
	set n_rand [ expr { rand()*$ran_max } ]

	set sum 0
	for { set i 0 } { $i < $NUM } { incr i 1 } {
		set weight [ lindex $weights $i ]
		set sum [ expr $sum+$weight ]
		if { $sum >= $n_rand } {
			set cn_no [ expr $i+1 ]
			break;
		}
	} 
	set cn_no [ expr {int($cn_no)} ]
	return $cn_no
}

proc move_savedcn { m } {
   set check [ glob -nocomplain "../run_${m}" ]
   set existcheck [ llength $check ]
   if { $existcheck < 1 } {
      exec mkdir ../run_${m}
   }

   set count 0
   set check [ glob -nocomplain FFS*.cn.gz ]
   set NUM [ llength $check ]
   for { set i 0 } { $i < $NUM } { incr i } {
        incr count
	set y [ format %04d $count ]
	exec cp [ lindex $check $i ] ../run_${m}/FFS_${m}_${y}.cn.gz
	exec rm [ lindex $check $i ]
        if { $m == 0 } {
	    set fp [ open "../weight0.txt" a+ ]
            puts $fp "1"
            close $fp
        }
   }
   return $NUM
}


#*******************************************
# Main program starts here
#*******************************************

  # Monte Carlo simulation (brute force)
   
  if { $graphic=="On" } {
#  MD++ setnolog 
  MD++ printfreq = 100
  } else {
  MD++ printfreq = 1000000000000
  } 

  init_ising $n $h $t

  if { $graphic=="On" } {
  openwindow
  }   
   
  set mm [ expr $m-1 ]
  MD++ input = $lambdasequence assign_Lam

  if { $m == 0 } {
  MD++ allocDFS
  equil_run 
  }

### if lam_auto == "On" compute interface automatically using MCrun
  if { $lam_auto=="On" } {
     if { $m==0 } {
     lam_set $t $h $n
     }
     set lambdasequence [ read_lam ] 
     MD++ input = $lambdasequence assign_Lam
   }

###########################
## initial flux rate run
###########################

if { $m < 1 } {
  MD++ allocDFS
  run_setup 
  MD++ totalsteps = $step0totstep

  set lambda_init [ lindex $lambdasequence 2 ]
  puts "run step 0 toward lambda0 = $lambda_init"

  MD++ MCrun
  MD++ continue_curstep = 1

## ignor this part for single CPU jobs
  set fp [ open "../weight0.txt" w]
  close $fp
  set fp [ open "../history0.txt" w ]
  close $fp

  set count [ move_savedcn $m ]

## save initial flux steps ##
  set init_step [ MD++_Get FFScurstep ]
  set fp [ open "../total0.txt" w ]
  puts $fp "$init_step $count $count"
  close $fp
  set m 1
}

####################################
## 2nd and higher interfaces probability
#################################### 

for { set N 1 } { $N <= $NMAX } { incr N } {
set count 0
set initstep [ MD++_Get "FFScurstep" ]
MD++ FFSoption ${m}
set mmm [ expr $m + 1 ]
set mmmm [ expr $m + 2 ]
set mm  [ expr $m - 1 ]

set lambda0 [ lindex $lambdasequence ${mmm} ]
set lambdanext [ lindex $lambdasequence ${mmmm} ]
# to check current lambda that we're shooting from

# wait for processing previous step
## only used for multi-CPU runs, ignore if using single CPU
while { 1 } {
set check [ glob -nocomplain ../total${mm}.txt ]
set key [ llength $check ]
set check_ [ glob -nocomplain ../weight${mm}.txt ]
set key_ [ llength $check_ ]
set crosskey [ expr $key * $key_ ]
if { $crosskey > 0 } { break }
exec sleep 20
}

### if sussess probability of previous step
### is bigger than 90%, directly proceed to state B
set to_stateB 0
if { $m > 1 } {
   set fp [ open "../total${mm}.txt" r ]
   set data [ split [ read $fp ] ] 
   set suc [ lindex $data 1 ]
   set tot [ lindex $data 0 ]
   set P [ expr 1.0*$suc/$tot ]
   if { $P > 0.90 } {
   MD++ FFSautoend = 1
   set totalstepslimit [ expr $totalstepslimit*3 ]
   set singlestepslimit [ expr $singlestepslimit*3 ]
   set to_stateB 1
   }
}

if { ${to_stateB} == 0 } {
    puts "shoot from interface ${mm} from lambda${mm}=$lambda0 to lambda${m}=$lambdanext"
} else {
    puts "now shooting to state B"
}

### Shooting trial moves from previous configurations
for { set i 1 } { $i <= 100000 } { incr i } {

### obtain configuration number randomly
   set n_rand [ get_init_cn $m ]
   set n_rand_format [ format %04d $n_rand ]

### read the configuration file saved from previous step
   set filename FFS_${mm}_${n_rand_format}.cn
   MD++ incnfile="../run_${mm}/$filename"
   MD++ readcn
   MD++ allocDFS
   run_setup
   
### set interface number and lambda value each trial ####
   MD++ N_lgst_cluster = $lambda0
   MD++ FFScn_weight = 1
   MD++ FFS0_check = 0
   MD++ FFSoption = ${m}

### stop running if totalstep limit is reached ##
   set grosstotalstep [ MD++_Get "FFScurstep" ]

   set lim [ expr $totalstepslimit + $initstep ]
   if { $grosstotalstep > $lim || $count > $Ncut } {
	break;
   }
  
### after all setting, shoot trial move!
   MD++ totalsteps = $singlestepslimit
   MD++ MCrun

### processing trial shooting results
### we will get FFS0_check = 1 from md++ 
### if this trial reaches to next interface
   set check [ MD++_Get "FFS0_check" ]
   set cn_weight [ MD++_Get "FFScn_weight" ]

   # only for success case
   if { $check > 0 } {
   incr count
   
   set fp [ open "weight.txt" a+ ]
   puts $fp "$cn_weight"
   close $fp

   set fp [ open "success.txt" a+ ]
   puts $fp "${n_rand}"
   close $fp
   }

### log if this trial is success or not
   set Nfinal [ MD++_Get N_lgst_cluster ]
   set fp [ open "history.txt" a+ ]
   puts $fp "${n_rand} $check $lambda0 $Nfinal"
   close $fp

### log total # of trial and # of success
   set fp [ open "status.txt" w ]
   puts $fp "$i $count"
   close $fp

### keep curstep for further simulation
   MD++ continue_curstep = 1
}
## end of for loop of i, end of shooting trials ##

## after all trials movie change datafile names
exec cp weight.txt ../weight${m}.txt
exec rm weight.txt
exec cp success.txt ../success${m}.txt
exec rm success.txt
exec cp history.txt ../history${m}.txt
exec rm history.txt
exec rm status.txt

set count [ move_savedcn $m ]
set fp [ open "../total${m}.txt" w ]
puts $fp "$i $count $count"
close $fp

## increase m to do next step simulation
set m [ expr $m + 1 ]

## if current trial is toward state_B, stop simulation! it's done!
   set autoend [ MD++_Get FFSautoend ]
   if { $autoend== 1} {
        if { $n < 2 } {
        set fp [ open "../autoend.txt" w ]
        close $fp
        }
	break;
   }

}

MD++ quit

