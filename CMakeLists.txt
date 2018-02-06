cmake_minimum_required (VERSION 2.8)
project(sparse_demo)
enable_language (Fortran)


# Repertoire d'installation des modules Fortran
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/fortran)



add_subdirectory(SparseBLAS)




set( DEMO_SRCS assembly.f90)

# Executable compilé à partor du Fortran
add_executable( sparsedemo ${DEMO_SRCS} )
target_link_libraries( sparsedemo SparseBLAS )

# Génération d'un module Python avec f2py
add_custom_target( pyAssembly ALL DEPENDS plop )

message(${DEMO_SRCS})
add_custom_command( OUTPUT plop COMMAND f2py -m pyAssembly -c ${CMAKE_SOURCE_DIR}/${DEMO_SRCS} -L${CMAKE_BINARY_DIR}/SparseBLAS -lSparseBLAS -I${CMAKE_BINARY_DIR}/fortran)


