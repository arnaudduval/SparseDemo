find_program(F2PY_EXECUTABLE NAMES f2py f2py${PYTHON_VERSION_MAJOR})

if(F2PY_EXECUTABLE)
    #extract the version string
    execute_process(COMMAND "${F2PY_EXECUTABLE}" -v
                    OUTPUT_VARIABLE F2PY_VERSION_STRING
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    if("${F2PY_VERSION_STRING}" MATCHES "^([0-9]+)(.([0-9+]))?(.([0-9+]))?$")
        set(F2PY_VERSION_MAJOR ${CMAKE_MATCH_1})
        set(F2PY_VERSION_MINOR ${CMAKE_MATCH_3})
        set(F2PY_VERSION_PATCH ${CMAKE_MATCH_5})
    endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(F2PY
    REQUIRED_VARS F2PY_EXECUTABLE
    VERSION_VAR F2PY_VERSION_STRING
    )
    
mark_as_advanced(F2PY_EXECUTABLE)
