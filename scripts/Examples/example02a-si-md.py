import mdpp

# Run MD++ script
mdpp.cmd("""
# MD code of Stinger-Weber Silicon
setnolog
setoverwrite
dirname = runs/si-example

#Create Perfect Lattice Configuration
crystalstructure = diamond-cubic
latticeconst = 5.4309529817532409 #(A) for Si
latticesize  = [ 1 0 0 4
                 0 1 0 4
                 0 0 1 4 ]
makecrystal  writecn

#Plot Configuration
atomradius = 0.67 bondradius = 0.3 bondlength = 2.8285 #for Si
atomcolor = orange highlightcolor = purple  bondcolor = red backgroundcolor = gray70
#plot_color_bar = [ 1 -4.85 -4.50 ]  highlightcolor = red
plotfreq = 10  rotateangles = [ 0 0 0 1.25 ] #[ 0 -90 -90 1.5 ]
openwin  alloccolors rotate saverot eval plot

#Conjugate-Gradient relaxation
#
conj_ftol = 1e-7 conj_itmax = 1000 conj_fevalmax = 10000
conj_fixbox = 1 #conj_monitor = 1 conj_summary = 1
relax finalcnfile = relaxed.cn writecn

#MD settings           
T_OBJ = 300 #Kelvin
equilsteps = 0  totalsteps = 100 timestep = 0.0001 # (ps)
atommass = 28.0855 # (g/mol)
#atomTcpl = 200.0 boxTcpl = 20.0
DOUBLE_T = 0
srand48bytime
initvelocity totalsteps = 10000 saveprop = 0
#integrator_type = 1 #0: Gear6, 1: Velocity Verlet
saveprop = 1 savepropfreq = 10 openpropfile #run
#usenosehoover = 1  vt2 = 2e28  #1e28 2e28 5e28
ensemble_type = NVE integrator_type = VVerlet

#Test MD run
totalsteps = 10000
saveprop = 1 savepropfreq = 10 openpropfile
plotfreq = 1 #autowritegiffreq = 10
run finalcnfile = si100.cn writecn

""")

# Plot simulation data
print "Plot simulation data using pylab"
from pylab import *
prop = loadtxt('prop.out')
plot( prop[:,0], prop[:,9] )
xlabel("steps")
ylabel("T (K)")
ion()
show()

# Sleep 
import time
sleep_seconds = 60
print "Python is going to sleep for " + str(sleep_seconds) + " seconds."
time.sleep(sleep_seconds)
