cmake_minimum_required (VERSION 3.3)
project(sparse_demo)
enable_language (Fortran)


# Repertoire d'installation des modules Fortran
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/fortran)

# Modules de rechrech CMake supplémentaires
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")




find_package(PETSc REQUIRED)

#find_package(PkgConfig)
#pkg_search_module(PETSC REQUIRED PETSc)
message(STATUS "PETSc_DIR : ${PETSC_DIR}")
message(STATUS "PETSc_ARCH : ${PETSC_ARCH}")
message(STATUS "PETSc_FOUND : ${PETSc_FOUND}")
message(STATUS "PETSC_LIBRARIES : ${PETSC_LIBRARIES}")
message(STATUS "PETSC_INCLUDES : ${PETSC_INCLUDES}")
message(STATUS "PETSC_DEFINITIONS : ${PETSC_DEFINITIONS}")



include_directories(${PETSC_INCLUDES})
add_definitions(${PETSC_DEFINITIONS})



set(CMAKE_Fortran_FLAGS "-cpp")





# Génération d'un module Python avec f2py
#add_custom_target( pyAssembly ALL DEPENDS plop )

#message(${DEMO_SRCS})
#add_custom_command( OUTPUT plop COMMAND f2py -m pyAssembly -c ${CMAKE_SOURCE_DIR}/${DEMO_SRCS} -L${CMAKE_BINARY_DIR}/SparseBLAS -lSparseBLAS -I${CMAKE_BINARY_DIR}/fortran)


set( PETSC_EXAMPLE_SRC petsc-example.f90)
add_executable( petsc-example ${PETSC_EXAMPLE_SRC} )
#target_include_directories( petsc-example PUBLIC ${PETSC_INCLUDE_DIRS})
target_link_libraries( petsc-example ${PETSC_LIBRARIES})


