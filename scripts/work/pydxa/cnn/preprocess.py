''' 

Read in .cfg or .lammps or .npy file in folder ``base_folder'' 
Extract a chunk of the system with a fixed number of atoms
Put the extracted file in a destination folder with index. 
such as silicon/silicon_1.dat, silicon_2.dat, or nickel, etc.

'''

# from lammps.lib import lammps 
import numpy as np
import os 
import sys

global PY3, N, base_folder, destination_folder

if (sys.version_info[0] < 3):
    PY3 = False
else:
    PY3 = True
            
N = 1000 # Total number of extracted subsystem.
base_folder = "base_Folder"
destination_folder = "lattice/"

def extract(md_data):
    Natoms = md_data[0]
    cfg = md_data[1]
    H = md_data[2] 


    #convert coordinate system to scaled


    #adjust window isotropically, starting from the middle
    while ( not converged ):
        #adjust 3d window size Wx,Wy,Wz


        #check if atoms bounded by [W] is > or < N


        #set ratio of adjustment of new step 

    #adjust window in the direction that W_i * W_j is smaller
    while ( not converged ):
        #adjust 3d window size Wx,Wy,Wz


        #check if atoms bounded by [W] is > or < N


        #set ratio of adjustment of new step 
    

    #remove atoms from the box face to make box of N atoms


    return (Natoms, cfg, H)

def read_cfg_file(filename):
    lines = [line.rstrip('\n') for line in open(filename)]
    Natoms = int(lines[0].split()[1])
    H = np.zeros((3,3))
    H[0][0] = float(lines[2].split()[1])
    H[1][1] = float(lines[6].split()[1])
    H[2][2] = float(lines[10].split()[1])
    cfg = np.array([i.split(' ')[0:3] for i in lines[16::3]]).astype('float')
    assert(cfg.shape[0] == Natoms), "Filename = {0} is not read in properly".format(filename)
    return (Natoms, cfg, H)

def read_lammps_file(filename):
    if PY3:
        lines = list(filter(None, [line.rstrip('\n') for line in open(filename)]))
    else:
        lines = filter(None, [line.rstrip('\n') for line in open(filename)])

    Natoms = int(lines[1].split(' ')[0])
    Nspecies = int(lines[2].split(' ')[0])
    H = np.zeros((3,3))
    H[0][0] = float(lines[3].split()[1]) - float(lines[3].split()[0])
    H[1][1] = float(lines[4].split()[1]) - float(lines[4].split()[0])
    H[2][2] = float(lines[5].split()[1]) - float(lines[5].split()[0])
    cfg = np.array([i.split(' ')[2:] for i in lines[8+Nspecies::]]).astype('float')
    assert(cfg.shape[0] == Natoms), "Filename = {0} is not read in properly".format(filename)
    return (Natoms, cfg, H)

def read_npy_file(filename):
    cfg = load(filename)
    Natoms = len(cfg)
    cfg = np.array([i[3] for i in cfg)]).astype('float')
    H[0][0] = np.max(cfg[:,0])-np.min(cfg[:,0])
    H[1][1] = np.max(cfg[:,1])-np.min(cfg[:,1])
    H[2][2] = np.max(cfg[:,2])-np.min(cfg[:,2])
    return (Natoms, cfg, H)

def main(argv):
    for fentry in os.listdir(os.path.expanduser(base_folder)):
        filename, file_ext = os.path.splitext(fentry)

        if (file_ext == ".cfg"):
            cfg = read_cfg_file(filename)
        elif (file_ext == ".lammps"):
            cfg = read_lammps_file(filename)
        elif (file_ext == ".npy"):
            cfg = read_npy_file(filename)
        assert(cfg.shape[0] > N), "filename ={0} has {1} atoms which is less than prescribed {1}".format(filename, cfg.shape[0], N)

        subdata = extract(cfg)

        save(subdata, destination+"/lattice_"+index.str()+".npy")

if __name == '__main__':
    main(sys.argv[1:])
