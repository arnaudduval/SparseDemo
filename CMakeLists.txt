cmake_minimum_required (VERSION 3.3)
project(sparse_demo)
enable_language (Fortran)


# Repertoire d'installation des modules Fortran
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/fortran)

# Modules de rechrech CMake supplémentaires
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")


set(PETSC_EXECUTABLE_RUNS ON)

find_package(PETSc REQUIRED)
find_package(F2PY REQUIRED)

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
target_link_libraries( assembly ${PETSC_LIBRARIES})

#add_custom_target( pyAssembly ALL DEPENDS forwrap )
#add_custom_command( OUTPUT forwrap COMMAND f2py -m pyAssembly -c ${CMAKE_SOURCE_DIR}/${ASSEMBLY_SRCS} -I${PETSC_INCLUDES})

set(f2py_module_name "pyAssembly")
set(fortran_src_file "${CMAKE_CURRENT_SOURCE_DIR}/${ASSEMBLY_SRCS}")
set(generated_module_file ${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}${PYTHON_EXTENSION_MODULE_SUFFIX})
add_custom_target(${f2py_module_name} ALL DEPENDS ${generated_module_file})

#generate Include string for petsc
set(PETSC_INCLUDE_STRING )
foreach(dir ${PETSC_INCLUDES})
    message(STATUS ${dir})
    set(PETSC_INCLUDE_STRING ${PETSC_INCLUDE_STRING} -I${dir})
endforeach()
message(${PETSC_INCLUDE_STRING})

#generate Link string for petsc
set(PETSC_LINK_STRING)
foreach(lib ${PETSC_LIBRARIES})
    message(STATUS ${lib})
    get_filename_component(dir ${lib} DIRECTORY)
    message(STATUS ${dir})
    get_filename_component(libshort ${lib} NAME_WE)
    string(REPLACE "lib" "" libshort ${libshort})
    message(STATUS ${libshort})
    set(PETSC_LINK_STRING ${PETSC_LINK_STRING} -L${dir} -l${libshort})
endforeach()
message(${PETSC_LINK_STRING})

add_custom_command(OUTPUT ${generated_module_file} 
                   COMMAND ${F2PY_EXECUTABLE}
                   -m ${f2py_module_name}
                   -c
                   ${fortran_src_file}
                   ${PETSC_INCLUDE_STRING}
                   ${PETSC_LINK_STRING}
                   --f90flags="-cpp -ffree-line-length-none"
                   WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                   )



