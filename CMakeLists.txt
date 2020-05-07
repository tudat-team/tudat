#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.


# Specify minimum CMake version required to compile Tudat.
cmake_minimum_required( VERSION 2.8.10 )

# Clean cache

if(EXISTS "${CMAKE_BINARY_DIR}/CMakeCache.txt")
    file(REMOVE "${CMAKE_BINARY_DIR}/CMakeCache.txt")
endif()

# Specific Tudat project name.
project( TudatBundle )

# ://stackoverflow.com/questions/38296756/what-is-the-idiomatic-way-in-cmake-to-add-the-fpic-compiler-option
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

# Package options
option(TUDAT_WITH_CSPICE "Build CSPICE libary, tudat with SPICE support and example applications." ON)
#option(USE_JSON "Build JSON library and example applications." ON)
option(TUDAT_WITH_NRLMSISE00 "Build NRLMISE-00 library and tudat with NRLMSISE support." ON)
option(TUDAT_WITH_SOFA "Build Sofa library." ON)
option(TUDAT_WITH_PAGMO "Build PaGMO2 library. Forced to ON if optimization package is enabled" OFF)
option(TUDAT_WITH_PYGMO "Build PaGMO library with Python bindings." OFF)
option(TUDAT_WITH_NLOPT "Build PaGMO library with Python bindings." OFF)
option(TUDAT_WITH_TUDAT_EXAMPLE_APPLICATIONS "Build Tudat with example applications." ON)

# Set root-directory for code to current source directory.
set(PROJECTROOT "${CMAKE_CURRENT_SOURCE_DIR}/")
message(STATUS ${PROJECT_SOURCE_DIR})
enable_testing( )

# Set the flag to force to use include_directories macro on Apple too
option(APPLE_INCLUDE_FORCE "Force include_directories on Apple too." ON)

# Let CMake find our custom cmake modules/scripts
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/external/CMake" "${CMAKE_CURRENT_SOURCE_DIR}/tudat/Tudat/External/CMake")

#
# Boost
#
# Set the required boost components (for Tudat)
set(BoostComponents filesystem system regex date_time unit_test_framework thread)
# Add additional required boost componets for other libraries
# if these are used.
if(USE_PAGMO)
    list(APPEND BoostComponents serialization chrono atomic)
endif()
if(USE_PYGMO)
    list(APPEND BoostComponents python python3)
endif()
# Pick a boost version
set(BoostVersion 1.64.0)
set(BoostSHA1 51421ef259a4530edea0fbfc448460fcc5c64edb)
#set(BoostVersion 1.60.0)
#set(BoostSHA1 7f56ab507d3258610391b47fef6b11635861175a)
# set(BoostVersion 1.57.0)
# set(BoostSHA1 e151557ae47afd1b43dc3fac46f8b04a8fe51c12)
# set(BoostVersion 1.53.0)
# set(BoostSHA1 e6dd1b62ceed0a51add3dda6f3fc3ce0f636a7f3)
# set(BoostVersion 1.48.0)
# set(BoostSHA1 27aced5086e96c6f7b2b684bda2bd515e115da35)

#
# TudatPy
#
#if(USE_TUDATPY)
#
#    set(pybind11_INCLUDE_DIR "${PROJECTROOT}pybind11/include")
#    add_subdirectory("${PROJECTROOT}/pybind11/")
#    list(APPEND BoostComponents python3)
#    find_package(pybind11 2.0 REQUIRED)
#
#endif()

# Minimum boost version required for PagMO2 is 1.55
if(TUDAT_WITH_PAGMO)
    # Explode the BoostVersion string ##.##.## into boost_version;boost_release;##
    string(REPLACE "." ";" boost_v_list ${BoostVersion})
    list(GET boost_v_list 0 boost_version)
    list(GET boost_v_list 1 boost_release)
    math(EXPR boost_v_number "${boost_version} * 100000 + ${boost_release} * 100")
    # If calculated version < 105500 change to BoostVersion to 1.57.0
    if( boost_v_number LESS 105500 )
        set(BoostVersion 1.57.0)
        set(BoostSHA1 e151557ae47afd1b43dc3fac46f8b04a8fe51c12)
        message(STATUS "Setting boost version to ${BoostVersion}")
    endif()
    unset(boost_v_list)
    unset(boost_version)
    unset(boost_release)
endif()

# Build boost
include(add_boost)

#
# SPICE
#
if(TUDAT_WITH_CSPICE)
    add_subdirectory( "${PROJECTROOT}/cspice/")
    set_property(TARGET cspice APPEND PROPERTY COMPILE_OPTIONS "-Wno-implicit-int")
    set_property(TARGET cspice APPEND PROPERTY COMPILE_OPTIONS "-Wno-incompatible-pointer-types")
    set_property(TARGET cspice APPEND PROPERTY COMPILE_OPTIONS "-Wno-pointer-to-int-cast")
    set_property(TARGET cspice APPEND PROPERTY COMPILE_OPTIONS "-Wno-shift-op-parentheses")
    set_property(TARGET cspice APPEND PROPERTY COMPILE_OPTIONS "-Wno-logical-op-parentheses")
    set_property(TARGET cspice APPEND PROPERTY COMPILE_OPTIONS "-Wno-parentheses")
    set_property(TARGET cspice APPEND PROPERTY COMPILE_OPTIONS "-Wno-tautological-compare")
    set_property(TARGET cspice APPEND PROPERTY COMPILE_OPTIONS "-Wno-format")
    set_property(TARGET cspice APPEND PROPERTY COMPILE_OPTIONS "-Wno-unsequenced")
endif()

#
# JSON
#
#if(USE_JSON)
#    add_subdirectory( "${PROJECTROOT}/json/")
#endif()

#
# NRLMSISE
#
# Add the NRLMSISE project and build
if(TUDAT_WITH_NRLMSISE00)
    set(NRLMSISE00_WITH_TESTS OFF CACHE BOOL "Building of NRLMSISE00 tests.")
    add_subdirectory( "${PROJECTROOT}/nrlmsise-00/")
endif()

#
# SOFA
#
if(TUDAT_WITH_SOFA)
    add_subdirectory( "${PROJECTROOT}/sofa/src")
endif()

#
# PaGMO2
#
if(TUDAT_WITH_PAGMO)
    option(PAGMO_BUILD_TESTS   "Build PaGMO2 with test set." ON)
    option(PAGMO_BUILD_EXAMPLES "Build PaGMO2 with examples." ON)
    option(PAGMO_BUILD_TUTORIALS "Build PaGMO2 with tutorials." ON)
    option(PAGMO_WITH_EIGEN3 "Build PaGMO2 with Eigen 3." ON)
    option(PAGMO_USE_TBB "Build with TBB." OFF)
    option(PAGMO_AS_STATIC_LIB "Compile pagmo as a static instead of a dynamic library." ON)

    if(PAGMO_AS_STATIC_LIB)
        add_definitions(-DPAGMO_STATIC_LIB=1)
    else()
        add_definitions(-DPAGMO_STATIC_LIB=0)
    endif()


    if(TUDAT_WITH_NLOPT)
        option(PAGMO_WITH_NLOPT "Build PaGMO2 with NLOPT." ON)
    endif()

    if(TUDAT_WITH_PYGMO)
        option(PAGMO_BUILD_PAGMO "Do not build PaGMO2" OFF)
        option(PAGMO_BUILD_PYGMO "Build PaGMO2 with Python bindings" ON)
        message(STATUS "this")
    else(TUDAT_WITH_PYGMO)
        option(PAGMO_BUILD_PAGMO "Build PaGMO2" ON)
        option(PAGMO_BUILD_PYGMO "Do not build PaGMO2 with Python bindings" OFF)
        message(STATUS "not this ${PAGMO_BUILD_PAGMO}")
    endif(TUDAT_WITH_PYGMO)

    message(STATUS "The things are: ${PAGMO_BUILD_PAGMO} ${PAGMO_BUILD_PYGMO} ${PAGMO_BUILD_TESTS} ${PAGMO_BUILD_EXAMPLES}")
    unset(pagmo_LIB_DEPENDS CACHE)

    if(PAGMO_USE_TBB)
        set(TBB_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/tbb")
        add_subdirectory( "${PROJECTROOT}tbb/")
    endif()

    add_subdirectory( "${PROJECTROOT}pagmo2/" )

elseif(TUDAT_WITH_PYGMO)
    option(PAGMO_BUILD_PAGMO "Do not build PaGMO2" OFF)
    option(PAGMO_BUILD_PYGMO "Do not build PaGMO2 with Python bindings" ON)
    unset(pagmo_LIB_DEPENDS CACHE)
    add_subdirectory( "${PROJECTROOT}pagmo2/" )
endif()

#
# NLOPT
#
if(TUDAT_WITH_NLOPT)
    if(PAGMO_AS_STATIC_LIB)
        option(BUILD_SHARED_LIBS "Building NLOPT as static libraries" OFF)
    else()
        option(BUILD_SHARED_LIBS "Building NLOPT as shared libraries" ON)
    endif()

    option(NLOPT_TESTS "Enabling NLOPT tests" ON)

    add_subdirectory( "${PROJECTROOT}/nlopt/")
endif()


#
# Tudat and applications
#
option(BUILD_WITH_ESTIMATION_TOOLS "Build Tudat with state estimation functionality." ON)
if(NOT BUILD_WITH_ESTIMATION_TOOLS)
    add_definitions(-DBUILD_WITH_ESTIMATION_TOOLS=0)
else()
    add_definitions(-DBUILD_WITH_ESTIMATION_TOOLS=1)
endif()
set(TUDAT_BUILD_DIR "${PROJECTROOT}/tudat")
add_subdirectory( "${PROJECTROOT}/Tudat/")

if(USE_TUDAT_EXAMPLE_APPLICATIONS)
    add_subdirectory( "${PROJECTROOT}/tudatExampleApplications/satellitePropagatorExamples/SatellitePropagatorExamples/")

    if(TUDAT_WITH_PAGMO)
        add_subdirectory( "${PROJECTROOT}/tudatExampleApplications/libraryExamples/PaGMOEx")
    endif()
endif()

#
##
## TudatPy
##
#if(USE_TUDATPY)
#    add_subdirectory("${PROJECTROOT}/tudatpy")
#endif()
#

add_subdirectory( "${PROJECTROOT}/tudatExampleApplications/templateApplication/TemplateApplication")