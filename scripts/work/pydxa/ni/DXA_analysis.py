# this script is only verified with ovito-2.7.1
from ovito.io import import_file, export_file
import ovito.data as data
from ovito.modifiers import DislocationAnalysisModifier
#from scipy.optimize import leastsq
import numpy as np
import os
import sys

def f(x, A, B): # this is your 'straight line' y=f(x)
    return A*x + B

def f_min(X,p):
    plane_xyz = p[0:3]
    distance = (plane_xyz*X.T).sum(axis=1) + p[3]
    return distance / np.linalg.norm(plane_xyz)
	
def residuals(params, signal, X):
    return f_min(X, params)

dirname = "./"
print(sys.argv)

observationstep = int(sys.argv[1])

if (observationstep < 1):
    print("Data points fewer than one!")
    exit()

loop_length_recorder = [] 

label_crystal = 1 #Diamond=0. FCC=1. BCC=2.
label_dislocation = 1 #perfect=0. partial=1.
bigvector= np.empty((0,))
datafilename = dirname + '/data'

for snap in range(0,1+observationstep):
    print("processing snap {0}".format(snap))
    dumpfile=dirname + "/ni_0_"+str(snap)+".lammps.gz"
    print(dumpfile)
    snap_str = "{0:0>4}".format(snap)    

    npyfilename = dirname+'/segpts'+snap_str

    if not os.path.isfile(dumpfile):
        print("There is no lammps dumpfile in "+dirname)
        continue
#    if os.path.isfile(npyfilename+".npy"):
#        print(npyfilename+" already exists in "+dirname)
#        continue

    node = import_file(dumpfile)
# Extract dislocation lines from a crystal with diamond structure:
    modifier = DislocationAnalysisModifier()
    modifier.input_crystal_structure = DislocationAnalysisModifier.Lattice.FCC
    node.modifiers.append(modifier)
    data = node.compute()

    total_line_length = node.output.attributes['DislocationAnalysis.total_line_length']
    dislocation_density = total_line_length / node.output.cell.volume
    print("Dislocation density: %f" % (dislocation_density))
#    print("Dislocation total length: %f" % total_line_length)
#    print(node.output.attributes)

    # Print list of dislocation lines:
    network = node.output.dislocations
    print("Found %i dislocation segments" % len(network.segments))
    coord = data.particle_properties.position.array
    Na = coord.shape[0]

    bigvector = np.concatenate((bigvector, [label_crystal],[label_dislocation], [total_line_length], [Na*3], coord[:,0], coord[:,1], coord[:,2]), axis=0) 

    if ( len(network.segments) != 0):
        dislpts = []
        for segment in network.segments:
#            print("Segment %i: length=%f, Burgers vector=%s" % (segment.id, segment.length, segment.true_burgers_vector))
#            print(segment.points)  
            if (segment.id == 0):
                dislpts = segment.points
            else:
                dislpts = np.concatenate((dislpts, segment.points), axis=0)
        np.save(dirname+'segpts'+snap_str, dislpts) 

np.save(datafilename, bigvector)

