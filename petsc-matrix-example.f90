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
    
    call MPI_Comm_size(PETSC_COMM_WORLD, size, ierr)
    if (size .ne. 1) then
        call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
        if(rank .eq. 0) then
            write(6,*) "This is a uniprocessor example only!"
        endif
        SETERRA(PETSC_COMM_WORLD, 1, " ")
    endif
    
    !! For testing only, allow the user to chosse grid size at runtime
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-m', m, flg, ierr);CHKERRA(ierr)
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-n', n, flg, ierr);CHKERRA(ierr)
    
    !! Create empty sparse matrix and linear solver data structures
    call UserInitializeLinearSolver(m, n, userctx, ierr);CHKERRA(ierr)
    
    !! Allocate arrays
    allocate(userx(m,n), userb(m,n), solution(m,n))
    allocate(rho(m,n))
    
    !! Fill the array rho with the function rho(x,y) = x
    !! Fill the RHS b and the solution with a known problem for testing
    hx = 1.0/real(m+1)
    hy = 1.0/real(n+1)
    y = hy
    do 20 j=1, n
        x = hx
        do 10 i=1, m
            rho(i,j) = x
            solution(i,j) = sin(2.*PETSC_PI*x)*sin(2.*PETSC_PI*y)
            userb(i,j) = -2.*PETSC_PI*cos(2.*PETSC_PI*x)*       &
                         sin(2.*PETSC_PI*y) +                   &
                         8*PETSC_PI*PETSC_PI*x*                 &
                         sin(2.*PETSC_PI*x)*sin(2.*PETSC_PI*y)
            x = x + hx
 10     continue
        y = y + hy
 20 continue
 
    !! Loop over time steps
    do 100 t=1, tmax
        call UserDoLinearSolver(rho, userctx, userb, userx, ierr);CHKERRA(ierr)    

        !! Compute error
        cnorm = 0.0
        do 90 j=1,n
            do 80 i=1,n
                cnorm = cnorm +                                                             &
                        PetscConj(solution(i,j)-userx(i,j))*(solution(i,j)-userx(i,j))
 80         continue
 90     continue
        enorm = PetscRealPart(cnorm*hx*hy)
        write(6,115) m, n, enorm
 115    format('m = ',I2,' n = ',I2,' error norm = ',1PE11.4)
 
 100 continue

    !! Clean up data structures
    deallocate(userx, userb, solution, rho)
    
    call UserFinalizeLinearSolver(userctx, ierr);CHKERRA(ierr)
    call PetscFinalize(ierr)    
    
end program main


subroutine UserInitializeLinearSolver(m, n, userctx, ierr)
    use UserModule
    implicit none
    
    PetscInt m, n
    PetscErrorCode ierr
    type(User) userctx
    
    common /param/ hx2, hy2
    PetscReal hx2, hy2
    
    Mat A
    Vec b, x
    KSP ksp
    PetscInt Ntot, five, one
    
    hx2 = (m+1)*(m+1)
    hy2 = (n+1)*(n+1)
    Ntot = m*n
    
    five = 5
    one = 1

    !! Create sparse matrix, preallocate 5 nonzeros per row
    call MatCreateSeqAIJ(PETSC_COMM_SELF, Ntot, Ntot, five, PETSC_NULL_INTEGER, A, ierr);CHKERRQ(ierr)
    
    !! Create vectors with no memory allocated
    call VecCreateSeqWithArray(PETSC_COMM_SELF, one, Ntot, PETSC_NULL_SCALAR, b, ierr);CHKERRQ(ierr)
    call VecDuplicate(b, x, ierr);CHKERRQ(ierr)
    
    !! Create linear solver context
    call KSPCreate(PETSC_COMM_SELF, ksp, ierr);CHKERRQ(ierr)
    
    userctx%x = x
    userctx%b = b
    userctx%A = A
    userctx%ksp = ksp
    userctx%m = m
    userctx%n = n
    
    return
end subroutine UserInitializeLinearSolver

!! Solves -div (rho grad psi) = F using finite differences
subroutine UserDoLinearSolver(rho, userctx, userb, userx, ierr)
    use UserModule
    implicit none
    
    PetscErrorCode ierr
    type(User) userctx
    PetscScalar rho(*), userb(*), userx(*)
    
    common /param/ hx2, hy2
    PetscReal hx2, hy2
    
    PC pc
    KSP ksp
    Vec b, x
    Mat A
    PetscInt m, n, one    
    PetscInt i, j, II, JJ
    PetscScalar v
    
    one = 1
    x = userctx%x
    b = userctx%b
    A = userctx%A
    ksp = userctx%ksp
    m = userctx%m
    n = userctx%n
    
    !! Compute operator -div rho grad
    II = 0
    do 110 j=1, n
        do 100 i=1, m
            if (j .gt. 1) then
                JJ = II - m
                v = -0.5*(rho(II+1) + rho(JJ+1))*hy2
                call MatSetValues(A, one, II, one, JJ, v, INSERT_VALUES, ierr);CHKERRQ(ierr)
            endif
            if (j .lt. n) then
                JJ = II + m
                v = -0.5*(rho(II+1) + rho(JJ+1))*hy2
                call MatSetValues(A, one, II, one, JJ, v, INSERT_VALUES, ierr);CHKERRQ(ierr)
            endif
            if(i .gt. 1) then
                JJ = II - 1
                v = -0.5*(rho(II+1) + rho(JJ+1))*hx2
                call MatSetValues(A, one, II, one, JJ, v, INSERT_VALUES, ierr);CHKERRQ(ierr)      
            endif
            if(i .lt. m) then
                JJ = II + 1
                v = -0.5*(rho(II+1) + rho(JJ+1))*hx2
                call MatSetValues(A, one, II, one, JJ, v, INSERT_VALUES, ierr);CHKERRQ(ierr)
            endif
            v = 2.*rho(II+1)*(hx2*hy2)
            call MatSetValues(A, one, II, one, II, v, INSERT_VALUES, ierr);CHKERRQ(ierr)                
            II = II + 1
 100    continue
 110 continue
 
    !! Assemble matrix
    call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr);CHKERRQ(ierr)
    
    !! Set operators
    call KSPSetOperators(ksp, A, A, ierr);CHKERRQ(ierr)
    
    !! Set linear solver defaults
    call KSPGetPC(ksp, pc, ierr);CHKERRQ(ierr)
    call PCSetType(pc, PCLU, ierr);CHKERRQ(ierr)
    
    !! Set Runtime options
    !! ex : -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    call KSPSetFromOptions(ksp, ierr);CHKERRQ(ierr)
    
    !! Allow PETSc linear solver to directly compute in user's array rather than PETSc vectors
    !! it is not recommanded for beginners
    call VecPlaceArray(x, userx, ierr);CHKERRQ(ierr)
    call VecPlaceArray(b, userb, ierr);CHKERRQ(ierr)
    
    !! Solve the linear system
    call KSPSolve(ksp, b, x, ierr);CHKERRQ(ierr)
    
    call VecResetArray(x, ierr);CHKERRQ(ierr)
    call VecResetArray(b, ierr);CHKERRQ(ierr)
    
    return
end subroutine UserDoLinearSolver


subroutine UserFinalizeLinearSolver(userctx, ierr)
    use UserModule
    implicit none
    
    PetscErrorCode ierr
    type(User) userctx
    
    call VecDestroy(userctx%x, ierr);CHKERRQ(ierr)
    call VecDestroy(userctx%b, ierr);CHKERRQ(ierr)
    call MatDestroy(userctx%A, ierr);CHKERRQ(ierr)
    call KSPDestroy(userctx%ksp, ierr);CHKERRQ(ierr)
    
    return

end subroutine UserFinalizeLinearSolver

    
