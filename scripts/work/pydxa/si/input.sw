# Input file for uniaxial tensile loading of Ag PentaTwin
# Amin Aghaei March 2014

# ------------------------ INITIALIZATION ----------------------------
units 		metal
dimension	3
boundary	p	p	p
atom_style	atomic
# processors  2 2 *

# ----------------------- ATOM DEFINITION ----------------------------
read_data    makeHole/SiGeHole_sw.lammps
# read_restart  geom.restart.20000

# ------------------------ FORCE FIELDS ------------------------------
pair_style sw
pair_coeff * * Si.sw Si 

# ------------------------- SETTINGS ---------------------------------
compute mytemp  all temp
compute atpoten all pe/atom
compute atstres all stress/atom mytemp ke pair
compute csym    all centro/atom 4
# compute atvoron all voronoi/atom
velocity        all create 1000 87287 loop geom

# ---------------------- BOUNDARY CONDITIONS -------------------------
region fixBC1 block INF INF INF INF INF 5.0 units box
region fixBC2 block INF INF -2.0 2.0 INF 5.0 units box
group botFace region fixBC1
group botLine region fixBC2

######################################

variable tmp equal "lx"
variable Lx0 equal ${tmp}
variable tmp equal "ly"
variable Ly0 equal ${tmp}
variable tmp equal "lz"
variable Lz0 equal ${tmp}

variable srate equal 2.0e8
variable srate1 equal "v_srate / 1.0e12"

print "Initial X Length, Lx0: ${Lx0}"
print "Initial Y Length, Ly0: ${Ly0}"
print "Initial Z Length, Lz0: ${Lz0}"
print "Strain rate 1/ps: ${srate1}"

variable Xstrain equal "(lx - v_Lx0) / v_Lx0 * 100.0"
variable Ystrain equal "(ly - v_Ly0) / v_Ly0 * 100.0"
variable Zstrain equal "(lz - v_Lz0) / v_Lz0 * 100.0"
variable XYstrain equal "2.0 + step * v_srate1 * 0.001 * 100.0"
variable stepN equal "step"
variable tm equal "temp"
variable ti equal "time / 1000.0"  ## Time in nano-seconds
variable l1 equal "lx"
variable l2 equal "ly"
variable l3 equal "lz"
variable px equal "-pxx/10000"    ### Stress in GPa (pxx is in bar)
variable py equal "-pyy/10000"
variable pz equal "-pzz/10000"
variable pxy equal "-pxy/10000"
variable pyz equal "-pyz/10000"
variable pxz equal "-pxz/10000"

variable aStresX atom "c_atstres[1] / 10000 / (20.9)"   ### Converting to GPa (dividing by the volume per atom = 20.9)
variable aStresY atom "c_atstres[2] / 10000 / (20.9)"   ### Instead of 20.9 you can use (vol/360045)
variable aStresZ atom "c_atstres[3] / 10000 / (20.9)" 

######################################

reset_timestep  0
timestep 0.001

fix  centCy all ave/atom 1 1000 10000 c_csym
fix  astrsX all ave/atom 1 1000 10000 v_aStresX
fix  astrsY all ave/atom 1 1000 10000 v_aStresY
fix  astrsZ all ave/atom 1 1000 10000 v_aStresZ

######################################
# RUN THE SIMULATION

variable srate2 equal "-v_srate1" 
velocity botFace set NULL NULL 0.0
# velocity botLine set NULL 0.0  0.0
fix 100 botFace setforce NULL NULL 0.0
# fix 101 botLine setforce NULL 0.0  0.0

# (FIRST STEP) FOR HAVING THERMAL RELAXATION ALONG Z UNCOMMNET THE FOLLOWING LINE
#fix 1 all npt temp 800.0 800.0 5.0 z 0.0 0.0 5.0
fix 1 all nvt temp 800.0 800.0 5.0  

# (SECOND STEP) FOR STRETCHING ALONG Z WITH CONSTANT TEMPERATURE UNCOMMENT THE FOLLOWING TWO LINES
# NOTE: strain rate of 1e8 1/s means 0.1% strain after 10000 steps
# fix 1 all nvt temp 800.0 800.0 5.0
fix 2 all deform 1 x erate ${srate2} units box remap x

# (THIRD STEP) FOR HAVING THERMAL RELAXATION AT CONSTANT TEMPERATURE AND LENGTH UNCOMMNET THE FOLLOWING LINE
# fix 1 all nvt temp 1000.0 1000.0 5.0

# dump          1 all custom/gz 10000 DXA/dump/x*.gz id type x y z c_csym v_aStresX v_aStresY v_aStresZ
dump 1 all cfg 2000 DXA/dump_sw/dump.snap.*.cfg mass type xs ys zs
#dump_modify   1 format "%d %d %15.8e %15.8e %15.8e %14.7e %14.7e %14.7e %14.7e"

restart  10000 geom.restart
thermo   1000
thermo_style	custom step v_Zstrain v_pz temp press

fix  def1 all print 1 "${stepN} ${ti} ${tm} ${Xstrain} ${Ystrain} ${Zstrain} ${px} ${py} ${pz}" append load_sw.out screen no

run    255000
unfix  2
run    4000000
undump 1

######################################
# SIMULATION DONE
print "All done"
