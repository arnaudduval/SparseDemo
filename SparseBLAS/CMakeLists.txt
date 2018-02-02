#file(DOWNLOAD "http://www.netlib.org/toms/818.gz" "${CMAKE_BINARY_DIR}/SparseBLAS.gz")

#add_custom_target(SparseBLASExpand)
#add_custom_command(TARGET SparseBLASExpand COMMAND "gunzip ${CMAKE_BINARY_DIR}/SparseBLAS.gz" DEPENDS "${CMAKE_BINARY_DIR}/SparseBLAS.gz" WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

#execute_process(COMMAND ${CMAKE_COMMAND} -E gunzip ${CMAKE_BINARY_DIR}/SparseBLAS.gz WORKING_DIRECTORY ${CMAKE_BINARY_DIR})


file( GLOB SPARSE_BLAS_SRCS src/*.f90)

add_library(SparseBLAS SHARED ${SPARSE_BLAS_SRCS})



set_target_properties(SparseBLAS PROPERTIES LINKER_LANGUAGE Fortran)
