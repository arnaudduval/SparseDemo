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

    use blas_sparse

    implicit none
    
    integer nmax, nnz
    parameter(nmax=4, nnz=6)
    integer :: i, n, a, istat
    integer, dimension(:), allocatable :: indx, jndx
    double precision, dimension(:), allocatable :: val, x, y
    
    allocate(val(nnz), x(nmax), y(nmax), indx(nnz), jndx(nnz))
    
    indx=(/1,2,2,3,4,4/)
    jndx=(/1,2,4,3,1,4/)
    val=(/1.1,2.2,2.4,3.3,4.1,4.4/)
    
    
    write(*,*) "Hello world !"
    
    n = nmax
    
    !! 1) create sparse blas handle
    call duscr_begin(n, n, a, istat)
    
    !! 2) insert entries one by one
    
    do i =1, nnz
        call uscr_insert_entry(A, val(i), indx(i), jndx(i), istat)
    enddo
    
    !! 3) complete construction of sparse matrix
    
    call uscr_end(a, istat)
    
    !! 4) compute matrix vector product y = A*x
    
    call usmv(a,x,y,istat)
    
    !! 5) Release matrix handle
    
    call usds(a, istat)
    
    

end program test
