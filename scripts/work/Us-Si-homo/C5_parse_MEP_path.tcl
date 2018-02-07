
if { $argc > 0 } {
set eps [ lindex $argv 0 ]
set T [ lindex $argv 1 ]
} else {
set eps 110
set T 0
}



MD++ setnolog
MD++ setoverwrite
MD++ dirname = runs/ho-cu-DN/${T}K_${eps}eps_112/NEB

#--------------------------------------------
#Read in potential file
MD++ potfile = "~/Planet/Libs/MD++/potentials/EAMDATA/eamdata.CuMishin" eamgrid=5000 readeam
MD++ NNM = 600
#MD++ potfile = ~/Codes/MD++/potentials/ta_pot readpot

#------------------------------------------------------------
#Read in structure
#
MD++ incnfile = NEBinit.cn readcn setconfig1 #A
MD++ incnfile = NEBfinal.cn readcn setconfig2 #B
MD++ incnfile = NEBinit.cn readcn #A

MD++ {
#------------------------------------------------------------
#
#Plot Configuration
#
atomradius = 0.9 bondradius = 0.3 bondlength = 2.8285 #for Si
atomcolor = orange highlightcolor = purple  bondcolor = red
backgroundcolor = gray fixatomcolor = yellow
color00 = "orange"  color01 = "purple" color02 = "green"
color03 = "magenta" color04 = "cyan"   color05 = "purple"
color06 = "gray80"  color07 = "white"
plot_atom_info = 1
plotfreq = 10
rotateangles = [ 0 0 0 1.0 ]
win_width = 600 win_height = 600
#openwin alloccolors rotate saverot refreshnnlist eval plot
}

MD++ zipfiles = 0
for { set iter 0 } { $iter <= 25 } { incr iter 1 } {
MD++ incnfile ="stringrelax.chain.cn" readRchain

# set chain number
set chain_no $iter
set total_no 25

MD++ input = $iter  copyRchaintoCN eval

MD++ openwin alloccolors rotate saverot plot 

# write cn or cfg file
MD++ finalcnfile = "chain_no_${chain_no}.cn" writecn

}
MD++ sleep



MD++ quit
