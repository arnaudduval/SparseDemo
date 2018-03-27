import sys
sys.path.append('build')
sys.path.append('.')
sys.path.append("/home/arnaud/petsc/arch-linux2-c-debug/lib")
import numpy as np
import pyAssembly
from mpi4py import MPI
import petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc


filePath = sys.argv[1]


connectivity = np.loadtxt(filePath, dtype='int')


# Preallocation parameter for sparse matrices
NNZ = 500

nNodePerElement = connectivity.shape[1]-1
nNodes = np.amax(connectivity[:,1:])
nElts = connectivity.shape[0]

#s = pyAssembly.buildmat.assemble(connectivity, nNodes, nElts, nNodePerElement)
#print(s)


# Assemble matrix in a Python context

A = PETSc.Mat()
A.create(PETSc.COMM_SELF)
A.setSizes([nNodes, nNodes])
A.setType("seqaij")
A.setPreallocationNNZ(NNZ)

for iElement in range(nElts):
    m,n,v = pyAssembly.buildmat.element(connectivity, nelts=nElts, nnodes=nNodes, nnodeperelement=nNodePerElement, ielement=iElement)
    A.setValues(m,n,v,addv=True)

A.assemblyBegin()
A.assemblyEnd()
