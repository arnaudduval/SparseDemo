module UserModule
    !! User defined contex that contains all the data structures used in the linear solution process
#include <petsc/finclude/petscksp.h>
    use petscksp
    type User
        Vec x
        Vec b
        Mat A
        KSP ksp
        PetscInt m
        PetscInt n
    end type User     
end module UserModule
    
program main
    use UserModule
    implicit none
    
    external UserInitializeLinearSolver
    external UserFinalizeLinearSolver
    external UserDoLinearSolver
    
    PetscScalar hx, hy, x, y
    type(User) userctx
    PetscErrorCode ierr
    PetscInt m, n, t, tmax, i, j
    PetscBool flg
    PetscMPIInt size, rank
    PetscReal enorm
    PetscScalar cnorm
    PetscScalar, allocatable :: userx(:,:)
    PetscScalar, allocatable :: userb(:,:)
    PetscScalar, allocatable :: solution(:,:)
    PetscScalar, allocatable :: rho(:,:)
    
    PetscReal hx2, hy2
    common /param/ hx2, hy2
    
    tmax = 2
    m = 6
    n = 7
    
    
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    if(ierr .ne. 0) then
        print*,"Unable to initialize PETSc"
        stop
    endif
    
   
    


end program main
