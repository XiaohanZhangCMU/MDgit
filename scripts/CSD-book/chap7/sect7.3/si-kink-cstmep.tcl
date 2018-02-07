# -*-shell-script-*-
# MD code of Stillinger-Weber Si
#  constrained minimization to find saddle points
#  kink migration on 30 partial dislocations
#
# Compile by: make sw build=R
# run by:     bin/sw_gpp scripts/CSD-book/chap7/sect7.3/si-kink-cstmep.script
#

MD++ {
#setnolog
setoverwrite
dirname = runs/si-kink-cstmep
}

MD++ {
#------------------------------------------------------------
#Read in structure
#
incnfile = ~/Codes/MD++.01-07-2007/structures/CSD-book/chap7/partial30dkp-6x6x4relax.cn readcn setconfig1 #A
incnfile = ~/Codes/MD++.01-07-2007/structures/CSD-book/chap7/partial30rkp-6x6x4relax.cn readcn setconfig2 #B
incnfile = ~/Codes/MD++.01-07-2007/structures/CSD-book/chap7/partial30dkp-6x6x4relax.cn readcn #A (current configuration)
NIC = 200
}

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
plot_limits = [ 1 -10 10 -10 10  -0.05  0.05 ]
plot_color_windows = [ 2
                          -10 -4.5 0  #color00 = orange
                          -4.5 10  1  #color01 = purple
                     ]
plot_atom_info = 1
plotfreq = 10
rotateangles = [ 0 0 0 1.5 ]
win_width = 600 win_height = 600
openwin alloccolors rotate saverot refreshnnlist eval plot
#sleep quit
}

MD++ {
#-------------------------------------------------------------
# Conjugate-Gradient relaxation with constraints
#
conj_ftol = 1e-7 conj_itmax = 1000 conj_fevalmax = 3000 conj_fixbox = 1

constrainatoms = [ 109 284 285 1406 1407 1409 1410 1411 1412 1413 1414
1415 1417 1421 1424 1425 1427 1428 1429 1430 1431 1434 1435 1436 1437
1438 1439 1442 1443 1444 1450 1451 1452 1453 1454 1455 1456 1457 1458
1459 1460 1461 1463 1472 1473 1474 1482 1483 1484 1486 1610 1611 1612
1614 1615 1632 1633 1638 1640 1642 1643 1656 2548 2558 2565 2566 2573
2576 2577 2578 2579 2580 2581 2582 2583 2584 2585 2589 2590 2591 2594
2595 2596 2597 2604 2605 2606 2607 2608 2609 2610 2612 2613 2614 2624
2625 2626 2627 2628 2629 2745 2762 2766 2784 2785 2787 2789 2790 2791
2792 ] #R=8A
}

set fileID [open "EpotvsS.dat" w]
for { set x 0 } { $x <= 100 } { incr x 5 } {
    set y [expr $x/100.0]
    puts "constrained minimization s = $y"
    MD++ constrainS = $y constrainedrelax finalcnfile = relaxed$x.cn writecn
    puts $fileID "$y [format %21.14e [MD++_Get EPOT]]"
    puts         "E = [format %21.14e [MD++_Get EPOT]]"
}
close $fileID

MD++ quit
