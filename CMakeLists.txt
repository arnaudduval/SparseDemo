cmake_minimum_required (VERSION 2.8)
project(sparse_demo)
enable_language (Fortran)


# Repertoire d'installation des modules Fortran
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/fortran)



add_subdirectory(blas_interface)
