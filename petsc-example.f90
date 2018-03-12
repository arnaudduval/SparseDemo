!! Basic Petsc example
!! compile with gfortran petsc-example.f90 -I/usr/lib/petscdir/3.7.6/x86_64-linux-gnu-real/include/petsc/ -I/usr/lib/petscdir/3.7.6/x86_64-linux-gnu-real/include -I/usr/lib/petscdir/3.7.6/x86_64-linux-gnu-real/include/petsc/mpiuni -cpp -lpetsc -o petsc-example

program main

    implicit none
    
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"

    PetscErrorCode ierr
    Vec x
    
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    call VecCreateSeq(PETSC_COMM_SELF, 100, x, ierr)
    call VecSet(x, 1., ierr)
    call PetscFinalize(ierr)
    

end program main
