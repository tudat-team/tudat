# FindCSPICE.cmake - Attempt to find CSPICE libraries and include paths

# This module defines:
#   CSPICE_FOUND - True if headers and requested libraries were found
#   CSPICE_INCLUDE_DIRS - Where to find the headers
#   CSPICE_LIBRARIES - List of libraries when using CSPICE
#   CSPICE_VERSION - The version of the found CSPICE

if(NOT TARGET CSPICE::cspice)
    # Define the version of CSPICE you are looking for
    set(CSPICE_VERSION "67.0.0")

    # Search for the include directory containing cspice/SpiceUsr.h
    find_path(CSPICE_INCLUDE_DIR_PARENT
            NAMES cspice/SpiceUsr.h
            PATHS ${CMAKE_PREFIX_PATH} ${CMAKE_INSTALL_PREFIX}
            PATH_SUFFIXES include
            DOC "Parent directory where cspice/SpiceUsr.h can be found"
    )

    # Adjust CSPICE_INCLUDE_DIRS to the parent of the found include directory
    if(CSPICE_INCLUDE_DIR_PARENT)
        set(CSPICE_INCLUDE_DIRS "${CSPICE_INCLUDE_DIR_PARENT}")
    endif()

    # Search for the main CSPICE library
    find_library(CSPICE_LIBRARY
            NAMES cspice
            PATHS ${CMAKE_PREFIX_PATH} ${CMAKE_INSTALL_PREFIX}
            PATH_SUFFIXES lib
            DOC "Main CSPICE library"
    )

    # Search for the CSPICE support library
    find_library(CSPICE_SUPPORT_LIBRARY
            NAMES csupport.66 csupport.67 csupport csupport.a
            PATHS ${CMAKE_PREFIX_PATH} ${CMAKE_INSTALL_PREFIX}
            PATH_SUFFIXES lib
            DOC "CSPICE support library"
    )

    # Aggregate found components
    include(FindPackageHandleStandardArgs)
    message(STATUS "CSPICE_LIBRARY: ${CSPICE_LIBRARY}")
    message(STATUS "CSPICE_SUPPORT_LIBRARY: ${CSPICE_SUPPORT_LIBRARY}")
    message(STATUS "CSPICE_INCLUDE_DIRS: ${CSPICE_INCLUDE_DIRS}")
    message(STATUS "CSPICE_VERSION: ${CSPICE_VERSION}")
    # Print general CMake search paths
    message(STATUS "CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}")
    message(STATUS "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
    message(STATUS "CMAKE_SYSTEM_PREFIX_PATH: ${CMAKE_SYSTEM_PREFIX_PATH}")
    message(STATUS "CMAKE_SYSTEM_FRAMEWORK_PATH: ${CMAKE_SYSTEM_FRAMEWORK_PATH}")
    message(STATUS "CMAKE_SYSTEM_APPBUNDLE_PATH: ${CMAKE_SYSTEM_APPBUNDLE_PATH}")
    message(STATUS "CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")
    message(STATUS "CMAKE_LIBRARY_PATH: ${CMAKE_LIBRARY_PATH}")
    message(STATUS "CMAKE_INCLUDE_PATH: ${CMAKE_INCLUDE_PATH}")
    message(STATUS "CMAKE_PROGRAM_PATH: ${CMAKE_PROGRAM_PATH}")
    message(STATUS "PKG_CONFIG_PATH: $ENV{PKG_CONFIG_PATH}")
    message(STATUS "CMAKE_FIND_ROOT_PATH: ${CMAKE_FIND_ROOT_PATH}")
    find_package_handle_standard_args(CSPICE
            REQUIRED_VARS CSPICE_LIBRARY CSPICE_INCLUDE_DIRS
            VERSION_VAR CSPICE_VERSION
    )

    if(CSPICE_FOUND)
        # Create imported targets for CSPICE
        add_library(CSPICE::cspice UNKNOWN IMPORTED)
        set_target_properties(CSPICE::cspice PROPERTIES
                IMPORTED_LOCATION "${CSPICE_LIBRARY}"
                INTERFACE_INCLUDE_DIRECTORIES "${CSPICE_INCLUDE_DIRS}"
        )

        add_library(CSPICE::csupport UNKNOWN IMPORTED)
        set_target_properties(CSPICE::csupport PROPERTIES
                IMPORTED_LOCATION "${CSPICE_SUPPORT_LIBRARY}"
                INTERFACE_INCLUDE_DIRECTORIES "${CSPICE_INCLUDE_DIRS}"
        )

        # Aggregate the libraries to a variable for external usage
        set(CSPICE_LIBRARIES CSPICE::cspice CSPICE::csupport)

        message(STATUS "Found CSPICE: ${CSPICE_INCLUDE_DIRS} (found version ${CSPICE_VERSION})")
    else()
        message(STATUS "Could not find CSPICE")
    endif()

    mark_as_advanced(CSPICE_INCLUDE_DIRS CSPICE_LIBRARY CSPICE_SUPPORT_LIBRARY CSPICE_LIBRARIES CSPICE_VERSION)
endif()
