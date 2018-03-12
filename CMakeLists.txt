cmake_minimum_required (VERSION 3.3)
project(sparse_demo)
enable_language (Fortran)


# Repertoire d'installation des modules Fortran
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/fortran)



#find_package(PETSc)

find_package(PkgConfig)
pkg_search_module(PETSC REQUIRED PETSc)
message(STATUS "TEST : ${PETSC_FOUND}")
message(STATUS "TEST2 : ${PETSC_INCLUDE_DIRS}")

set(CMAKE_Fortran_FLAGS "-cpp")



# Génération d'un module Python avec f2py
#add_custom_target( pyAssembly ALL DEPENDS plop )

#message(${DEMO_SRCS})
#add_custom_command( OUTPUT plop COMMAND f2py -m pyAssembly -c ${CMAKE_SOURCE_DIR}/${DEMO_SRCS} -L${CMAKE_BINARY_DIR}/SparseBLAS -lSparseBLAS -I${CMAKE_BINARY_DIR}/fortran)


set( PETSC_EXAMPLE_SRC petsc-example.f90)
add_executable( petsc-example ${PETSC_EXAMPLE_SRC} )
target_include_directories( petsc-example PUBLIC ${PETSC_INCLUDE_DIRS})
target_link_libraries( petsc-example petsc)


