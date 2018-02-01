import petsc4py, sys

petsc4py.init(sys.argv)

from petsc4py import PETSc

# tailles
m,n = 32, 32
hx = 1.0/(m-1)
hy = 1.0/(n-1)

# creation de matrice
A = PETSc.Mat()
A.create(PETSc.COMM_WORLD)
A.setSizes([m*n,m*n])
A.setType('aij') # sparse
A.setUp() # A ajouter pour ne pas avoir d'erreur MPI

# precompute values
diagv = 2.0/hx**2 + 2.0/hy**2
offdx = -1.0/hx**2
offdy = -1.0/hy**2

# loop
Istart, Iend = A.getOwnershipRange()
print(Istart, Iend)
for I in range(Istart, Iend):
    A[I,I] = diagv
    i = I//n
    j = I-i*n
    if i>0:
        J = I-n
        A[I,J] = offdx
    if i < m-1:
        J = I+n
        A[I,J] = offdx
    if j > 0:
        J = I-1
        A[I,J] = offdy
    if j < n-1:
        J = I+1
        A[I,J] = offdy
        
A.assemblyBegin()
A.assemblyEnd()


# Create Linear Solver
ksp = PETSc.KSP()
ksp.create(PETSc.COMM_WORLD)

# Use conjugate gradients
ksp.setType('cg')
# And incomplete Cholesky
ksp.getPC().setType('icc')

# obtain sol & rhs vector
x, b = A.getVecs()
x.set(0)
b.set(1)

# and solve
ksp.setOperators(A)
ksp.setFromOptions()
ksp.solve(b, x)

try:
    from matplotlib import pylab
except:
    raise SystemExit("matplotlib not available")
    
from numpy import mgrid

X, Y = mgrid[0:1:1j*m, 0:1:1j*m]
Z = x[...].reshape(m,n)
pylab.figure()
pylab.contourf(X,Y,Z)
pylab.plot(X.ravel(), Y.ravel(), '.k')
pylab.axis('equal')
pylab.colorbar()
pylab.show()

