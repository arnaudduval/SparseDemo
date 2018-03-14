!! Programme qui permet de faire l'assemblage d'une matrice type raideur

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





program test
    
    implicit none
    
    integer, allocatable, dimension(:,:) :: connectivity
    integer :: nNodes, nNodePerElement, nElts
    integer :: i, j, dummy
    double precision, allocatable :: stiffness(:,:)
     
     
    !!nNodes = 36926
    !!nNodePerElement = 64
    !!nElts = 11712
    nNodes = 16
    nNodePerElement = 9
    nElts = 4
    allocate(connectivity(nElts, nNodePerElement+1))
    allocate(stiffness(nNodes, nNodes))
    
    
    !!open(11, file="/home/arnaud/sparse_demo/fixed.txt", form="formatted")
    open(11, file="/home/arnaud/sparse_demo/fixed2.txt", form="formatted")
    do i = 1, nElts
        read(11,*) (connectivity(i,j), j=1, nNodePerElement+1)
    enddo
    
    
    
    call Assemble(connectivity, nElts, nNodes, nNodePerElement, stiffness)
    
    deallocate(stiffness, connectivity)

end program test

