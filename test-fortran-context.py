import sys
sys.path.append('build')
sys.path.append('.')
sys.path.append("/home/arnaud/petsc/arch-linux2-c-debug/lib")
import numpy as np
import pyAssembly
from mpi4py import MPI

filePath = sys.argv[1]

connectivity = np.loadtxt(filePath, dtype='int')


# Preallocation parameter for sparse matrices
NNZ = 500

nNodePerElement = connectivity.shape[1]-1
nNodes = np.amax(connectivity[:,1:])
nElts = connectivity.shape[0]

#s = pyAssembly.buildmat.assemble(connectivity, nNodes, nElts, nNodePerElement)
#print(s)


# Assemble matrix in a Fortran context
pyAssembly.buildmat.startpetsc()
pyAssembly.buildmat.assemblepetscblock(connectivity, nelts=nElts, nnodes=nNodes, nnodeperelement=nNodePerElement, nz=NNZ)
pyAssembly.buildmat.stoppetsc()


