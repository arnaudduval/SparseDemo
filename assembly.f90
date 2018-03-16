!! Programme qui permet de faire l'assemblage d'une matrice type raideur

module buildMat
#include <petsc/finclude/petscmat.h>
    use petscmat    
    real(8), pointer :: pA(:)
    integer, pointer :: pia(:), pja(:)
    Mat :: stiffness, K
    
contains

subroutine Element(connectivity, nElts, nNodes, nNodePerElement, iElement, m, n, v)
    integer, intent(in) :: nElts, nNodes, nNodePerElement, iElement
    integer, dimension(nElts, nNodePerElement+1), intent(in) :: connectivity
    
    integer, dimension(nNodePerElement), intent(out) :: m, n
    real(8), dimension(nNodePerElement*nNodePerElement), intent(out) :: v
    
    integer i, j, idx
      
    do i = 1, nNodePerElement
        do j = 1, nNodePerElement
            idx = j + (i-1)*nNodePerElement !! Vérifier si on est en Row major ou pas
            v(idx) = 1.0
            m(i) = connectivity(iElement, i+1) - 1
            n(j) = connectivity(iElement, j+1) - 1
        enddo
    enddo

end subroutine Element

subroutine Assemble(connectivity, nElts, nNodes, nNodePerElement, Dstiffness)
    implicit none
    integer, intent(in)                                         :: nElts, nNodePerElement, nNodes
    integer,dimension(nElts, nNodePerElement+1), intent(in)     :: connectivity
    double precision, dimension(nNodes, nNodes), intent(out)      :: Dstiffness
    
    
    integer :: i, j, iElement
    
    Dstiffness = 0.0
    
    do iElement = 1, nElts
        do i = 1, nNodePerElement
            do j = 1, nNodePerElement
                Dstiffness(connectivity(iElement, i+1), connectivity(iElement, j+1)) &
                         & = Dstiffness(connectivity(iElement, i+1), connectivity(iElement, j+1)) + 1.0
            enddo
        enddo
    enddo
    
end subroutine Assemble

subroutine AssemblePETSc(connectivity, nElts, nNodes, nNodePerElement, nz)
#include <petsc/finclude/petscmat.h>
    use petscmat
    
    implicit none
    integer, intent(in) :: nElts, nNodePerElement, nNodes, nz
    integer, dimension(nElts, nNodePerElement+1), intent(in)    :: connectivity
    
    integer :: i, j, iElement, cpt
    PetscErrorCode :: ierr
    PetscScalar :: val(1)
    PetscInt :: II(1), JJ(1)
    PetscScalar :: one
    PetscInt :: n
    PetscBool :: done
    
    one = 1.0
    
    !! Initialisation la matrice
    !! nz = nombre de valeur non nulles sur chaque ligne, utilisé pour la préallocation
    call MatCreateSeqAIJ(PETSC_COMM_SELF, nNodes, nNodes, nz, PETSC_NULL_INTEGER, stiffness, ierr);CHKERRA(ierr)
    
    do iElement = 1, nElts
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

end subroutine AssemblePETSc


subroutine AssemblePETScBlock(connectivity, nElts, nNodes, nNodePerElement, nz)
#include <petsc/finclude/petscmat.h>
    use petscmat
    
    implicit none
    integer, intent(in) :: nElts, nNodePerElement, nNodes, nz
    integer, dimension(nElts, nNodePerElement+1), intent(in) :: connectivity
    
    integer :: iElement, cpt
    PetscErrorCode :: ierr
    PetscInt :: n
   
    real(8), dimension(nNodePerElement*nNodePerElement) :: val
    integer, dimension(nNodePerElement) :: II, JJ
    
        !! Initialisation la matrice
    !! nz = nombre de valeur non nulles sur chaque ligne, utilisé pour la préallocation
    call MatCreateSeqAIJ(PETSC_COMM_SELF, nNodes, nNodes, nz, PETSC_NULL_INTEGER, stiffness, ierr);CHKERRA(ierr)
    
    do iElement = 1, nElts
        call Element(connectivity, nElts, nNodes, nNodePerElement, iElement, II, JJ, val)    
        call MatSetValues(stiffness, nNodePerElement, II, nNodePerElement, JJ, val, ADD_VALUES, ierr);CHKERRA(ierr)
    enddo 
    
    !! Assemblage
    call MatAssemblyBegin(stiffness, MAT_FINAL_ASSEMBLY, ierr);CHKERRA(ierr)
    call MatAssemblyEnd(stiffness, MAT_FINAL_ASSEMBLY, ierr);CHKERRA(ierr)
    
end subroutine AssemblePETScBlock

subroutine ExtractCSRStructure
    PetscErrorCode  :: ierr
    PetscBool       :: done
    PetscInt        :: n
    
    call MatSeqAIJGetArrayF90(stiffness, pA, ierr)
    call MatGetRowIJF90(stiffness, 0, PETSC_FALSE, PETSC_TRUE, n, pIA, pJA, done, ierr)

end subroutine ExtractCSRStructure

subroutine GetValue(i,j,v)
    integer, intent(in) :: i, j
    real(8), intent(out) :: v
    
    integer :: ii(1), jj(1)
    real(8) :: vv(1)
    
    integer::  ierr
    
    ii(1) = i
    jj(1) = j
    
    
    call MatGetValues(stiffness, 1, ii, 1, jj, vv, ierr)

    v = vv(1)
end subroutine GetValue

subroutine GetNRows(nrows)
    implicit none
    integer, intent(out) :: nrows
    integer              :: ncols
    PetscErrorCode       :: ierr
    
    call MatGetSize(stiffness, nrows, ncols, ierr)

end subroutine GetNRows

subroutine GetNz(nz)
    integer, intent(out) :: nz
    
    nz = ubound(pA,1)

end subroutine GetNz

subroutine GetCSR(nz, nrows, A, IA, JA)
    integer, intent(in) :: nz
    integer, intent(in) :: nrows
    real(8), dimension(nz), intent(out) :: A
    integer, dimension(nrows+1), intent(out) :: IA
    integer, dimension(nz), intent(out) :: JA
    
    A = pA
    IA = pIA
    JA = pJA


end subroutine GetCSR


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
    double precision, allocatable :: Dstiffness(:,:)
    
    real(PETSC_REAL_KIND), allocatable :: A(:)
    integer, allocatable :: ia(:), ja(:)
    
    PetscErrorCode ierr
    PetscBool flg
     
     
    nNodes = 36926
    nNodePerElement = 64
    nElts = 11712
    !!nNodes = 9
    !!nNodePerElement = 4
    !!nElts = 4
    allocate(connectivity(nElts, nNodePerElement+1))
    !!allocate(stiffness(nNodes, nNodes))

    
    open(11, file="/home/arnaud/sparse_demo/fixed.txt", form="formatted")
    !!open(11, file="/home/arnaud/sparse_demo/fixed3.txt", form="formatted")
    do i = 1, nElts
        read(11,*) (connectivity(i,j), j=1, nNodePerElement+1)
    enddo
    
    !!call Assemble(connectivity, nElts, nNodes, nNodePerElement, stiffness)
    
    call StartPETSc
    

    call AssemblePETSc(connectivity, nElts, nNodes, nNodePerElement, 500)
    !! Copie de stiffness
    call MatConvert(stiffness, MATSAME, MAT_INITIAL_MATRIX, K, ierr);CHKERRA(ierr)
    call MatDestroy(stiffness, ierr);CHKERRA(ierr)
    
    !! Assemblage avec blocs
    call AssemblePETScBlock(connectivity, nElts, nNodes, nNodePerElement, 500)

    call MatEqual(stiffness, K, flg, ierr);CHKERRA(ierr)

    
    write(*,*) 'MatEqual', flg
    
    !!deallocate(stiffness, connectivity)
    
    call PetscFinalize(ierr);CHKERRA(ierr)

end program test

