import sys
sys.path.append('/home/arnaud/sparse_demo/build')
import numpy as np
import pyAssembly as ass

# Read connectivity matrix
connectivity = np.loadtxt("fixed2.txt", dtype='int', delimiter=',')

stiffness = ass.assembleblas(connectivity, nnodes=np.amax(connectivity[:,1:]), nelts=connectivity.shape[0], nnodeperelement=connectivity.shape[1]-1)
