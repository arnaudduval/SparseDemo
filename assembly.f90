!! Programme qui permet de faire l'assemblage d'une matrice type raideur

module buildMat
    
    real(8), pointer :: pA(:)
    integer, pointer :: pia(:), pja(:)
contains

subroutine Assemble(connectivity, nElts, nNodes, nNodePerElement, stiffness)
    implicit none
    integer, intent(in)                                         :: nElts, nNodePerElement, nNodes
    integer,dimension(nElts, nNodePerElement+1), intent(in)     :: connectivity
    double precision, dimension(nNodes, nNodes), intent(out)      :: stiffness
    
    
    integer :: i, j, iElement
    
    stiffness = 0.0
    
    do iElement = 1, nElts
        do i = 1, nNodePerElement
            do j = 1, nNodePerElement
                stiffness(connectivity(iElement, i+1), connectivity(iElement, j+1)) &
                         & = stiffness(connectivity(iElement, i+1), connectivity(iElement, j+1)) + 1.0
            enddo
        enddo
    enddo
    
end subroutine Assemble

subroutine AssemblePETSc(connectivity, nElts, nNodes, nNodePerElement)
#include <petsc/finclude/petscmat.h>
    use petscmat
    use iso_c_binding
    
    implicit none
    integer, intent(in) :: nElts, nNodePerElement, nNodes
    integer, dimension(nElts, nNodePerElement+1), intent(in)    :: connectivity

    


    !! Variable locale pour le moment
    Mat stiffness, B
    
    integer :: i, j, iElement, cpt
    PetscErrorCode :: ierr
    PetscScalar :: val(1)
    PetscInt :: II(1), JJ(1)
    PetscScalar :: one
    PetscInt :: n
    PetscBool :: done
    
    one = 1.0
    
    !! Initialisation la matrice
    !!call MatCreateSeqAIJ(PETSC_COMM_SELF, nNodes, nNodes, PETSC_DEFAULT_INTEGER, PETSC_NULL_INTEGER, stiffness, ierr);CHKERRA(ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, nNodes, nNodes, 500, PETSC_NULL_INTEGER, stiffness, ierr);CHKERRA(ierr)
    
    do iElement = 1, nElts
        write(*,*) iElement, '/', nElts
        do i = 1, nNodePerElement
            do j = 1, nNodePerElement
                !! -1 sur les indices car on fait du stockage qui commence à l'indice 0
                !! En principe, on devrait pouvoir faire du stockage qui commence à l'indice 1
                call MatSetValue(stiffness, connectivity(iElement, i+1)-1, connectivity(iElement, j+1)-1, 1.0_PETSC_REAL_KIND, ADD_VALUES, ierr);CHKERRA(ierr)
            enddo
        enddo
    enddo 
    
    !! Assemblage
    call MatAssemblyBegin(stiffness, MAT_FINAL_ASSEMBLY, ierr);CHKERRA(ierr)
    call MatAssemblyEnd(stiffness, MAT_FINAL_ASSEMBLY, ierr);CHKERRA(ierr)
    
   
    call MatSeqAIJGetArrayF90(stiffness, pA, ierr)
    !!write(*,*) "array"
    !!write(*,100) A
    call MatGetRowIJF90(stiffness, 0, PETSC_FALSE, PETSC_TRUE, n, pIA, pJA, done, ierr)
    
    write(*,*) done
    !!write(*,*) 'n', n
    !!write(*,*) 'IA', ia
    !!write(*,*) 'JA', ja
    
    
    
    call MatCreateSeqAIJWithArrays(PETSC_COMM_SELF, nNodes, nNodes, pIA, pJA, pA, B, ierr)
    write(*,*) lbound(pA), ubound(pA)

    
    !!A = pA
    !!IA = pIA
    !!JA = pJA
  
    call MatEqual(stiffness, B, done, ierr)
    write(*,*) "matrices egales", done
    
 100 format( 2000(F3.1,'   '))

end subroutine AssemblePETSc

subroutine StartPETSc
#include <petsc/finclude/petscsys.h>
    use petscsys
    PetscErrorCode ierr
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    if(ierr .ne. 0) then
        print *, "Unable to initialize PETSc"
        stop
    endif    

end subroutine StartPETSc



end module


program test
#include <petsc/finclude/petscmat.h>
    use petscmat    
    use BuildMat
    implicit none
    
    integer, allocatable, dimension(:,:) :: connectivity
    integer :: nNodes, nNodePerElement, nElts
    integer :: i, j, dummy
    double precision, allocatable :: stiffness(:,:)
    
    real(PETSC_REAL_KIND), allocatable :: A(:)
    integer, allocatable :: ia(:), ja(:)
    
    Mat K
    PetscErrorCode ierr
     
     
!!    nNodes = 36926
!!    nNodePerElement = 64
!!    nElts = 11712
    nNodes = 9
    nNodePerElement = 4
    nElts = 4
    allocate(connectivity(nElts, nNodePerElement+1))
    !!allocate(stiffness(nNodes, nNodes))

    
    !!open(11, file="/home/arnaud/sparse_demo/fixed.txt", form="formatted")
    open(11, file="/home/arnaud/sparse_demo/fixed3.txt", form="formatted")
    do i = 1, nElts
        read(11,*) (connectivity(i,j), j=1, nNodePerElement+1)
    enddo
    
    !!call Assemble(connectivity, nElts, nNodes, nNodePerElement, stiffness)
    
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    if(ierr .ne. 0) then
        print *, "Unable to initialize PETSc"
        stop
    endif
    

    call AssemblePETSc(connectivity, nElts, nNodes, nNodePerElement)
    
    !!deallocate(stiffness, connectivity)
    
    call PetscFinalize(ierr)

end program test

