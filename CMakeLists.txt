cmake_minimum_required (VERSION 3.3)
project(sparse_demo)
enable_language (Fortran)


# Repertoire d'installation des modules Fortran
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/fortran)

# Modules de rechrech CMake supplémentaires
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")


set(PETSC_EXECUTABLE_RUNS ON)

find_package(PETSc REQUIRED)

message(STATUS "PETSc_DIR : ${PETSC_DIR}")
message(STATUS "PETSc_ARCH : ${PETSC_ARCH}")
message(STATUS "PETSc_FOUND : ${PETSc_FOUND}")
message(STATUS "PETSC_LIBRARIES : ${PETSC_LIBRARIES}")
message(STATUS "PETSC_INCLUDES : ${PETSC_INCLUDES}")
message(STATUS "PETSC_DEFINITIONS : ${PETSC_DEFINITIONS}")



include_directories(${PETSC_INCLUDES})
add_definitions(${PETSC_DEFINITIONS})



set(CMAKE_Fortran_FLAGS "-cpp -ffree-line-length-none")





# Génération d'un module Python avec f2py
#add_custom_target( pyAssembly ALL DEPENDS plop )

#message(${DEMO_SRCS})
#add_custom_command( OUTPUT plop COMMAND f2py -m pyAssembly -c ${CMAKE_SOURCE_DIR}/${DEMO_SRCS} -L${CMAKE_BINARY_DIR}/SparseBLAS -lSparseBLAS -I${CMAKE_BINARY_DIR}/fortran)


set( PETSC_VECTOR_EXAMPLE_SRC petsc-vector-example.f90)
add_executable( petsc-vector-example ${PETSC_VECTOR_EXAMPLE_SRC} )
target_link_libraries( petsc-vector-example ${PETSC_LIBRARIES})

set( PETSC_MATRIX_EXAMPLE_SRC petsc-matrix-example.f90)
add_executable( petsc-matrix-example ${PETSC_MATRIX_EXAMPLE_SRC} )
target_link_libraries( petsc-matrix-example ${PETSC_LIBRARIES})

set( ASSEMBLY_SRCS assembly.f90)
add_executable( assembly ${ASSEMBLY_SRCS} )


