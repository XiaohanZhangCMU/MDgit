import os, sys
import numpy as np
import random
import pickle
from subprocess import call
from numpy import linalg as LA

from bitarray import bitarray

# exec(open("./global_search.py").read())

'''
 About this script:
 1) Require disl_nuc_hetero.tcl in the same folder.
    all configurations are dumped to MD++ dirname folder.
 2) status = 0  create free surface, surface reconstructed. applied strain.
    status = 1  read trenched configuration, relax with a little strain help
 3) Select two points (x,z) in the reference configuration.
    One on top surface one in the body. Both on (111) glide-set.
 4) Check slice.cfg(cn) and init_trench.cn for correctness

 Some notations:
   cfg: atom configuration for all atoms without any trench
   atomIdx: [0,1,2...NP], NP = # of atoms of cfg
   totIdx: indices for all atoms on the slice
   nucleus: indices for atoms that are removed from the slice
   nbrlist: nbrlist for all atoms on the slice (totIdx)
'''

# icd: initial cluster dimension
# make Nmax tries
# dataset: a dictionary { bits : potential energy }
global dirname, icdx, icdy, icdz, Nmax, dataset
dirname = 'runs/frankMEP/'
icdx = 0.42 #0.32
icdy = 0.42 #0.32
icdz = 0.42 #0.32
Nmax = 100
dataset = { }

# Return atom Id within an ecllipse given _SR
# pt1: center coordinate in _SR coordinate
# a, b: width and length
# normal: of the plane ecllipse lies on

def select_atoms_within_ecllipse(pt1, a, b, normal, cfg, atomIdx, totIdx):
  print("Select nucleus within an elipse.")
  e0 = [1,0,0]; e1 = [0,1,0];                e2 = [0,0,1];
  ep0 = e1;     ep1 = np.cross(normal, ep0); ep2 = normal; 
  
  # coordinate transform
  Q = np.zeros((3,3))

  Q[0][0] = np.dot(ep0,e0)
  Q[1][1] = np.dot(ep1,e1)
  Q[2][2] = np.dot(ep2,e2)

  Q[0][1] = np.dot(ep0, e1)
  Q[0][2] = np.dot(ep0, e2)

  Q[1][0] = np.dot(ep1, e0)
  Q[1][2] = np.dot(ep1, e2)

  Q[2][0] = np.dot(ep2, e0)
  Q[2][1] = np.dot(ep2, e1)

  nucleus = []; 
  for i in range(cfg.shape[0]):
    _x0 = cfg[i,:] - pt1; # origin shifts to the picked point
    _x = np.dot(Q,_x0);
    x = _x[0]; y = _x[1]; z = _x[2]; 
    if (x)*(x)/a/a + (y)*(y)/b/b <= 1 :
      nucleus.append(i);

  nucleus = np.array(nucleus)
  return np.intersect1d(nucleus, totIdx)
  

def select_initial_nucleus(pt1, cfg, atomIdx, totIdx):
  print("Select initial nucleus.")
  x0 = pt1[0]; y0 = pt1[1]; z0 = pt1[2];
  print("P0 = [ {0}, {1}, {2} ]".format(x0, y0, z0))
  cond = np.abs(cfg[:,2]-z0)<icdz;idx=np.extract(cond,np.arange(cfg.shape[0]));
  zatomIdx=np.extract(cond, atomIdx); zcfg = cfg[idx, :]
  print("I am here 1 ")
  cond=np.abs(zcfg[:,0]-x0)<icdx;idx=np.extract(cond,np.arange(zcfg.shape[0]));
  xatomIdx=np.extract(cond, zatomIdx); xcfg = zcfg[idx, :]
  print("I am here 2 ")
  cond=np.abs(xcfg[:,1]-y0)<icdy;idx=np.extract(cond,np.arange(xcfg.shape[0]));
  print("I am here 3 ")
  yatomIdx=np.extract(cond, xatomIdx); ycfg = xcfg[idx,:]
  # Check if yatomIdx and ycfg are consistent with each other
  print(yatomIdx)
  print(ycfg)
  assert(np.amax(ycfg-cfg[yatomIdx])==0)
  # The intersection  set of yatomIdx and totIdx
  return np.intersect1d(yatomIdx, totIdx)

def select_totalatoms_on_slice(eps, pt1, cfg, normal, d, atomIdx):
  print(normal)
  D = np.ones(cfg[:,0].shape) * d ; 
  cond = np.abs(normal[0] * cfg[:,0] + normal[1] * cfg[:,1]  + normal[2] * cfg[:,2] -D) < eps ;
  idx = np.extract(cond,atomIdx)
  return idx


def write_SR(data, finalcnfile):
  with open(dirname+finalcnfile+".dat", 'w') as the_file:
    for i in range(data.shape[0]):
      for j in range(3):
        the_file.write(str(data[i,j]) + " ")
      the_file.write('\n') 

def writecfg_fromdata(data, H, finalcnfile):
  with open(dirname+finalcnfile+".cn", 'w') as the_file:
    the_file.write(str(data.shape[0])+'\n')
    for i in range(data.shape[0]):
      for j in range(3):
        the_file.write(str(data[i,j]) + "\t")
      the_file.write('\n') 
    for i in range(3):
      for j in range(3):
        the_file.write(str(H[i,j]) + "\t")
      the_file.write('\n') 
    the_file.write('1 Mo\n')
    the_file.write(' 0 0')
  call([ "sw_mc2_mpich",  "scripts/work/MEP_frank_partial/disl_nuc_hetero.tcl", "100", finalcnfile ])
  
def removeNucleusAtoms(cfg, nucleus, atomIdx):
  idx = np.setdiff1d(atomIdx, nucleus)
  cfg = cfg[idx,:]
  return cfg
  
def find_bdy_atoms(nbrlist, nucleus, cfg):
  print("find boundary atoms of current nucleus.")
  bdy_list = []
  for atomId in nucleus:
    if (np.in1d(nbrlist[atomId,:], nucleus)).all():
      continue;
    else:
      bdy_list.append(atomId)
  return np.array(bdy_list)

# arg1: list of atoms to build nbrlist
# opt: False -- no build just load. True -- build nbrlist
def build_nbrlist(totIdx, cfg, maxnbrs, cutoff, opt):
  if not opt:
    nbrlist = np.load(dirname+"nbrlist.npy")
  else:
    print("Buidling neighbor list for atoms on slice.")
    print("Patience......")
    nbrlist = -1* np.ones((cfg.shape[0], maxnbrs)) 
    nbrlist = nbrlist.astype(int)
    for atom_I in totIdx:
      cnt = 0
      for atom_J in totIdx:
        if atom_I != atom_J:
          a = LA.norm(cfg[atom_I, :]-cfg[atom_J,:], 2) 
          if a < cutoff:
            nbrlist[atom_I, cnt] = atom_J
            cnt += 1
            assert (cnt < maxnbrs), "Need to increase maxnbrs."
            assert (cnt > 0), "Need to increase cutoff."
    np.save(dirname+"nbrlist.npy", nbrlist)
  return nbrlist

def save_obj(obj, name ):
  with open(dirname + name + '.pkl', 'wb') as f:
    pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
			
def load_obj(name ):
  with open(dirname + name + '.pkl', 'rb') as f:
    return pickle.load(f) 

def set_bits_for_nucleus(nucleus, totIdx):
  bits = np.in1d(totIdx, nucleus)
  bits = bits.astype(int)
  return bits
  
'''                        
   Main Program Starts Here

'''                       
status = 1
if status == 0: 
  call([ "sw_mc2_mpich",  "scripts/work/MEP_frank_partial/disl_nuc_hetero.tcl", "0", "n", "001", "1" ])
  exit(0)
 
elif status == 1:
  data=np.genfromtxt(dirname+"0K_0.0_relaxed_surf001.cn",skip_header=1,skip_footer=2);
  cfg  = data[:-3,:]
  H    = data[-3::,:]
  atomIdx = np.arange(cfg.shape[0])

  pt1d = np.array([-0.0208,0, 0.3765]); pt1u = np.array([-0.0208,0, 0.3908])
  pt2d = np.array([0.1250, 0, 0.1953]);  pt2u = np.array([0.1250, 0, 0.2083])
  pt1 = 0.5*(pt1d + pt1u);  pt2 = 0.5*(pt2d + pt2u)
  
  # Create a slice beneath free surface and dump to slice.cn(cfg)
  eps = 0.02
  v = pt2- pt1
  v = v/LA.norm(v)
  assert(np.abs(np.sqrt(v[0]*v[0]+v[1]*v[1] + v[2]*v[2])- LA.norm(v))<1e-10)
  normal = np.cross(v, [0,1,0])
  x0 = pt1[0]; y0 = pt1[1]; z0 = pt1[2];
  D = normal[0] * x0 + normal[1] * y0 + normal[2] * z0;
  totIdx = select_totalatoms_on_slice(eps, pt1, cfg, normal, D, atomIdx)
  writecfg_fromdata(cfg[totIdx,:], H, "slice")
  write_SR(cfg[totIdx,:], "slice")

  # Choose initial nucleus at the top free surface along x=0  
  # The safer way is to pick atom (x,z) from X-window
  #nucleus = select_initial_nucleus(pt1, cfg, atomIdx, totIdx)
  a0 = 0.1; b0 = 0.1;
  nucleus = select_atoms_within_ecllipse(pt1, a0, b0, normal, cfg, atomIdx, totIdx)   

  print("I am here")
  
  # Build neighborhood list for all atoms on the slice
  cutoff = 0.07
  maxnbrs = 40
  nbrlist = build_nbrlist(totIdx, cfg, maxnbrs, cutoff, False);
  
  # Check the following configurations to make sure the script is working correctly
  init_trenched = removeNucleusAtoms(cfg, nucleus, atomIdx)
  totl_trenched = removeNucleusAtoms(cfg, totIdx,  atomIdx)
  writecfg_fromdata(init_trenched, H, "trench_init")
  writecfg_fromdata(totl_trenched, H, "trench_tot")
  
  # Randomly choose atoms beyond (but connecting to) initial nucleus. 
  # Generating a graph of atomistic states each denoted by missing atoms list
  
  for step in range(Nmax):

    # Write out trench and removed atoms (nucleus) positions
    writecfg_fromdata(cfg[np.setdiff1d(totIdx,nucleus),:],H,"slice"+str(step))
    #np.setdiff1d(nucleus, atom_j)
    cfg_step = removeNucleusAtoms(cfg, nucleus, atomIdx)
    writecfg_fromdata(cfg_step, H, "trench_"+str(step))
    write_SR(cfg[nucleus,:], "nucleus-"+str(step))

    # Read trench_$step.cn then relax the state with help of applied strain
    # Calculate potential energy of the final state write to file "trench_energy.dat" on line step
    call([ "sw_mc2_mpich","scripts/work/MEP_frank_partial/disl_nuc_hetero.tcl","1",str(step),str(normal[0]), str(normal[1]), str(normal[2]), str(D), str(pt1[0]), str(pt1[1]), str(pt1[2]), str(icdx), str(icdy), str(icdz) ]);

    bits = set_bits_for_nucleus(nucleus, totIdx)
    data = np.loadtxt(dirname+'EPOT_2.dat')
    mybits=""
    for k in bits:
      mybits+=str(k)
    dataset[mybits] = data[step]
    print("----------------------------------------------------------size of dataset is {0}".format( len(dataset)) )

    a0 += (icdy-a0)/Nmax; b0 = a0;
    print('a0 = {0}, b0 = {1}'.format(a0,b0));
    nucleus = select_atoms_within_ecllipse(pt1, a0, b0, normal, cfg, atomIdx, totIdx)   

  save_obj(dataset, 'dataset')
#    # Find an atom i randomly in current nucleus
#    bdy_list = find_bdy_atoms(nbrlist,nucleus,cfg)
#    atom_i = random.sample(set(bdy_list),1)
#
#    # Find an atom j randomly attached to atom i
#    atom_j = -1
#    tries = 0
#    while(True):
#      inbrlist = nbrlist[atom_i,:].reshape(nbrlist[atom_i,:].size)
#      atom_j = random.sample(set(inbrlist),1)
#      tries += 1
#      if (not np.in1d(atom_j, nucleus)) or tries ==100:
#        break;
#    assert(tries<100 and atom_j !=-1), "Cannot find new place to remove atom."
#    # Remove atom j. update nucleus. 
#    nucleus = np.append(nucleus, atom_j)
 
  
  # end of for step in range(Nmax)



elif status == 2:
  # Read all energy states and connectivities and form graph
  exit(0)



# Find lowest energy barrier in the graph



