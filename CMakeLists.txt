#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.

# Specify minimum cmake version required to compile Tudat.
cmake_minimum_required(VERSION 3.3)

# Trick for recompiling changes for BUILD-TREE.
if (EXISTS "${CMAKE_BINARY_DIR}/CMakeCache.txt")
    file(REMOVE "${CMAKE_BINARY_DIR}/CMakeCache.txt")
endif ()

# Set cmake policy.
cmake_policy(SET CMP0069 NEW)

# Get version from file.
file(READ "version" Tudat_VERSION_LABEL)
string(REGEX MATCH "^([0-9]+)\\.([0-9]+)\\.([0-9]+)" Tudat_VERSION "${Tudat_VERSION_LABEL}")

# Specific project name.
project(Tudat VERSION ${Tudat_VERSION} LANGUAGES C CXX)

# Get lower string name of PROJECT_NAME for directory configuration.
string(TOLOWER ${PROJECT_NAME} PROJECT_NAME_LOWER)
string(TOUPPER ${PROJECT_NAME} PROJECT_NAME_UPPER)
set(CMAKE_POSITION_INDEPENDENT_CODE True)

# Add local module path.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/yolo")

# General system information.
message(STATUS "System name: ${CMAKE_SYSTEM_NAME}")

# +============================================================================
# BUILD OPTIONS
#  Offer the user the choice of defining the build variation.
# +============================================================================
# Build option: enable the test suite.
option(TUDAT_BUILD_TESTS "Build the test suite." ON)

# Build with propagation tests.
option(TUDAT_BUILD_WITH_PROPAGATION_TESTS "Build tudat with propagation tests. (>30 s propagations - Total test time > 10 minutes.)" OFF)

# Build with estimation tools.
option(TUDAT_BUILD_WITH_ESTIMATION_TOOLS "Build tudat with estimation tools." ON)

# Build option: enable tutorials.
option(TUDAT_BUILD_TUDAT_TUTORIALS "Build tutorials." ON)

# Build static library instead of dynamic.
option(TUDAT_BUILD_STATIC_LIBRARY "Build tudat as a static library, instead of dynamic." ON)

# Build static library instead of dynamic.
option(TUDAT_BUILD_WITH_FILTERS "Build Tudat with preliminary Kalman Filter functionality." OFF)

# Build sofa interface.
option(TUDAT_BUILD_WITH_SOFA_INTERFACE "Build Tudat with sofa interface." ON)

# Build json interface.
option(TUDAT_BUILD_WITH_JSON_INTERFACE "Build Tudat with json interface." OFF)

# Build with nrlmsise-00 atmosphere model.
option(TUDAT_BUILD_WITH_NRLMSISE00 "Build with nrlmsise-00 atmosphere model." OFF)

# Build pagmo-dependent code.
option(TUDAT_BUILD_WITH_PAGMO "Build Tudat with pagmo." OFF)

# Build extended precision propagation tools.
option(TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS "Build tudat with extended precision propagation tools." OFF)

# Add all build settings as definitions for downstream projects.
set(Tudat_DEFINITIONS "${Tudat_DEFINITIONS} -DTUDAT_BUILD_WITH_FILTERS=${TUDAT_BUILD_WITH_FILTERS}")
set(Tudat_DEFINITIONS "${Tudat_DEFINITIONS} -DTUDAT_BUILD_WITH_NRLMSISE00=${TUDAT_BUILD_WITH_NRLMSISE00}")
set(Tudat_DEFINITIONS "${Tudat_DEFINITIONS} -DTUDAT_BUILD_WITH_SOFA_INTERFACE=${TUDAT_BUILD_WITH_SOFA_INTERFACE}")
set(Tudat_DEFINITIONS "${Tudat_DEFINITIONS} -DTUDAT_BUILD_WITH_JSON_INTERFACE=${TUDAT_BUILD_WITH_JSON_INTERFACE}")
set(Tudat_DEFINITIONS "${Tudat_DEFINITIONS} -DTUDAT_BUILD_WITH_ESTIMATION_TOOLS=${TUDAT_BUILD_WITH_ESTIMATION_TOOLS}")
set(Tudat_DEFINITIONS "${Tudat_DEFINITIONS} -DTUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS=${TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS}")
add_definitions(${Tudat_DEFINITIONS})

# +============================================================================
# INSTALL TREE CONFIGURATION (Project name independent)
#  Offer the user the choice of overriding the installation directories.
# +============================================================================
# TODO: Consider if this can work with <tudat/paths.hpp> going forward.
set(INSTALL_LIB_DIR "${CMAKE_INSTALL_PREFIX}/lib/" CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR "${CMAKE_INSTALL_PREFIX}/bin/" CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}/include/" CACHE PATH "Installation directory for header files")
set(INSTALL_DATA_DIR "${CMAKE_INSTALL_PREFIX}/data/" CACHE PATH "Installation directory for data files")
set(INSTALL_MAN_DIR "${CMAKE_INSTALL_PREFIX}/" CACHE PATH "Installation directory for man documentation")
#set(INSTALL_TESTS_DIR "${CMAKE_INSTALL_PREFIX}/tests/" CACHE PATH "Installation directory for tests (default=OFF)")
set(INSTALL_CMAKE_DIR "${INSTALL_LIB_DIR}/cmake/${PROJECT_NAME_LOWER}" CACHE PATH "Installation directory for cmake config files")

# Make relative paths absolute (needed later on)
foreach (p LIB BIN INCLUDE DATA CMAKE)
    set(var INSTALL_${p}_DIR)
    set(RELATIVE_INSTALL_${p}_DIR ${INSTALL_${p}_DIR})
    if (NOT IS_ABSOLUTE "${${var}}")
        set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
    endif ()
endforeach ()

# Find relative paths for cmake-config file. (for find_package)
file(RELATIVE_PATH rel_include_dir "${INSTALL_CMAKE_DIR}" "${INSTALL_INCLUDE_DIR}")
file(RELATIVE_PATH rel_lib_dir "${INSTALL_CMAKE_DIR}" "${INSTALL_LIB_DIR}")
file(RELATIVE_PATH rel_data_dir "${INSTALL_CMAKE_DIR}" "${INSTALL_DATA_DIR}")
file(RELATIVE_PATH rel_bin_dir "${INSTALL_CMAKE_DIR}" "${INSTALL_BIN_DIR}")
file(RELATIVE_PATH rel_man_dir "${INSTALL_CMAKE_DIR}" "${INSTALL_MAN_DIR}")

# Set relative paths for config.cmake.
foreach (p include lib data bin man)
    string(TOUPPER ${p} P)
    set(RELATIVE_INSTALL_${P}_DIR ${rel_${p}_dir})
endforeach ()

# +============================================================================
# DEPENDENCIES SPECIFIC TO BUILD OPTIONS
# +============================================================================
if (NOT TUDAT_BUILD_STATIC_LIBRARY)
    # Required settings for shared library build type. See provided link for
    # more info.
    # https://stackoverflow.com/questions/30398238/cmake-rpath-not-working-could-not-find-shared-object-file
    set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif ()

# Set required boost components. This variable is used to indicate required
# boost libraries inside the `tudat/cmake_modules/TudatFindBoost.cmake` file.
set(_TUDAT_REQUIRED_BOOST_LIBS filesystem system regex date_time thread chrono atomic)

if (TUDAT_BUILD_TESTS)
    # If tests are built, then the boost unit_test_framework is needed.
    list(APPEND _TUDAT_REQUIRED_BOOST_LIBS unit_test_framework)
endif ()

# Set cmake build-type. If it not supplied by the user via a build argument,
# the default built type is "Release".
if ((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE STREQUAL "Release"))
    set(CMAKE_BUILD_TYPE Release)
elseif (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_BUILD_TYPE Debug)
endif ()

# Summarize the build configuration related section for debugging.
message(STATUS "=============================================================")
message(STATUS "******************** BUILD CONFIGURATION ********************")
message(STATUS "=============================================================")
message(STATUS "PROJECT_NAME ........................................ ${PROJECT_NAME}")
message(STATUS "CMAKE_BUILD_TYPE .................................... ${CMAKE_BUILD_TYPE}")
message(STATUS "BUILD_STYLE ......................................... ${BUILD_STYLE}")
message(STATUS "TUDAT_BUILD_TESTS ................................... ${TUDAT_BUILD_TESTS}")
message(STATUS "TUDAT_BUILD_WITH_PROPAGATION_TESTS .................. ${TUDAT_BUILD_WITH_PROPAGATION_TESTS}")
message(STATUS "TUDAT_BUILD_WITH_ESTIMATION_TOOLS ................... ${TUDAT_BUILD_WITH_ESTIMATION_TOOLS}")
message(STATUS "TUDAT_BUILD_TUDAT_TUTORIALS ......................... ${TUDAT_BUILD_TUDAT_TUTORIALS}")
message(STATUS "TUDAT_BUILD_STATIC_LIBRARY .......................... ${TUDAT_BUILD_STATIC_LIBRARY}")
message(STATUS "TUDAT_BUILD_WITH_FILTERS ............................ ${TUDAT_BUILD_WITH_FILTERS}")
message(STATUS "TUDAT_BUILD_WITH_SOFA_INTERFACE ..................... ${TUDAT_BUILD_WITH_SOFA_INTERFACE}")
message(STATUS "TUDAT_BUILD_WITH_JSON_INTERFACE ..................... ${TUDAT_BUILD_WITH_JSON_INTERFACE}")
message(STATUS "TUDAT_BUILD_WITH_NRLMSISE00 ......................... ${TUDAT_BUILD_WITH_NRLMSISE00}")
message(STATUS "TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS ${TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS}")

# +============================================================================
# COMPILER + SYSTEM SPECIFIC
# +============================================================================

# Include compiler.cmake which takes care of ALL compiler related CMake
# configurations, including known warning suppression for specific compilers
# and systems.
include(compiler)

# Summarize the build configuration related section for debugging.
message(STATUS "=============================================================")
message(STATUS "************* COMPILER/SYSTEM CONFIGURATION *****************")
message(STATUS "=============================================================")
message(STATUS "CMAKE_C_COMPILER .................................... ${CMAKE_C_COMPILER}")
message(STATUS "CMAKE_CXX_COMPILER .................................. ${CMAKE_CXX_COMPILER}")
message(STATUS "CMAKE_CXX_COMPILER_ID ............................... ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CMAKE_CXX_SIMULATE_ID ............................... ${CMAKE_CXX_SIMULATE_ID}")
message(STATUS "CMAKE_SYSTEM_PROCESSOR .............................. ${CMAKE_SYSTEM_PROCESSOR}")
message(STATUS "CMAKE_CXX_COMPILER_VERSION .......................... ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "CMAKE_SYSTEM_NAME ................................... ${CMAKE_SYSTEM_NAME}")
message(STATUS "WIN32 ............................................... ${WIN32}")
message(STATUS "UNIX ................................................ ${UNIX}")
message(STATUS "APPLE ............................................... ${APPLE}")

# +============================================================================
# MANDATORY DEPENDENCIES
# +============================================================================

# Include TudatFindBoost.cmake which is adapted from Pagmo directly.
include(TudatFindBoost)
if (UNIX OR (WIN32 AND CLANG))
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${Boost_INCLUDE_DIRS}\"")
endif ()

# Get the Tudat resources dependency.
find_package(TudatResources REQUIRED)
include_directories(TudatResources_INCLUDE_DIRS)

# Find Eigen3 library on local system and include directories.
find_package(Eigen3 REQUIRED)
include_directories(SYSTEM AFTER "${EIGEN3_INCLUDE_DIR}")
if (APPLE OR (WIN32 AND CLANG))
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${EIGEN3_INCLUDE_DIR}\"")
endif ()

# CSpice dependency
find_package(CSpice REQUIRED 1.0.0)

# Summarize the build configuration related section for debugging.
message(STATUS "=============================================================")
message(STATUS "***************** MANDATORY DEPENDENCIES ********************")
message(STATUS "=============================================================")
message(STATUS "CSpice_LIBRARIES .................................... ${CSpice_LIBRARIES}")
message(STATUS "CSpice_VERSION ...................................... ${CSpice_VERSION}")
message(STATUS "TudatResources_INCLUDE_DIRS ......................... ${TudatResources_INCLUDE_DIRS}")
message(STATUS "TudatResources_VERSION .............................. ${TudatResources_VERSION}")
message(STATUS "EIGEN3_INCLUDE_DIR: ................................. ${EIGEN3_INCLUDE_DIR}")
message(STATUS "EIGEN3_VERSION ...................................... ${EIGEN3_VERSION}")
message(STATUS "Boost_LIBRARIES  .................................... ${Boost_LIBRARIES}")
message(STATUS "Boost_VERSION_STRING  ............................... ${Boost_VERSION_STRING}")


# +============================================================================
# OPTIONAL DEPENDENCIES
# +============================================================================

# Sofa dependency if in build settings.
if (TUDAT_BUILD_WITH_PAGMO)
    # If Tudat builds with Pagmo, then the Sims Flannagan method is available.
    find_package(Pagmo REQUIRED)
endif ()

if (TUDAT_BUILD_WITH_SOFA_INTERFACE)
    find_package(Sofa REQUIRED 1.0.0)
endif ()

if (TUDAT_BUILD_WITH_NRLMSISE00)
    find_package(NRLMSISE00 REQUIRED 0.1)
endif ()

if (TUDAT_BUILD_WITH_JSON_INTERFACE)
    find_package(nlohmann_json REQUIRED 3.7.3)
endif ()

# +============================================================================
# SETUP ALL LIBRARIES AND CMAKE RELATED TOOLS
# +============================================================================

# Set compiler based on preferences (e.g. USE_CLANG) and system.
include(TudatLinkLibraries)

# Include all the custom scripts which provide the following functions for
# making life easier:
# - TUDAT_ADD_EXECUTABLE
# - TUDAT_ADD_LIBRARY
# - TUDAT_ADD_TEST_CASE
include(YOLOProjectAddTestCase)
include(YOLOProjectAddLibrary)
include(YOLOProjectAddExecutable)

# Add to parent scope if project is added via add_subdirectory().
if (NOT ${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
    set(Tudat_PROPAGATION_LIBRARIES ${Tudat_PROPAGATION_LIBRARIES} PARENT_SCOPE)
    set(Tudat_ESTIMATION_LIBRARIES ${Tudat_ESTIMATION_LIBRARIES} PARENT_SCOPE)
endif ()

add_subdirectory(src/astro)
add_subdirectory(src/basics)
add_subdirectory(src/interface)
add_subdirectory(src/math)
add_subdirectory(src/simulation)
add_subdirectory(src/io)
add_subdirectory(examples)


if (TUDAT_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif ()

file(RELATIVE_PATH TUDAT_DATA_DIR_RELATIVE_TO_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" "${INSTALL_DATA_DIR}/tudat")

configure_file(
        "${PROJECT_SOURCE_DIR}/config.hpp.in"
        "${PROJECT_BINARY_DIR}/include/tudat/config.hpp" @ONLY
)

# Installation of the header files.
install(DIRECTORY include/ DESTINATION "${INSTALL_INCLUDE_DIR}")
install(FILES
        "${PROJECT_BINARY_DIR}/include/tudat/config.hpp"
        DESTINATION ${INSTALL_INCLUDE_DIR}/tudat)

# Add all targets to the build-tree export set
export(TARGETS ${${PROJECT_NAME}_LIBRARIES} NAMESPACE ${PROJECT_NAME}:: FILE "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME_LOWER}_export.cmake")

# Installation of the export modules.
install(EXPORT ${PROJECT_NAME_LOWER}_export NAMESPACE ${PROJECT_NAME}:: DESTINATION "${RELATIVE_INSTALL_CMAKE_DIR}")

# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
export(PACKAGE ${PROJECT_NAME})

# +============================================================================
# PROJECT CONFIG
#  For external use of project.
# +============================================================================
include(CMakePackageConfigHelpers)

configure_package_config_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME_LOWER}-config.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME_LOWER}-config.cmake"
        INSTALL_DESTINATION "${INSTALL_CMAKE_DIR}")

configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME_LOWER}-config.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME_LOWER}-config.cmake" @ONLY
)

file(
        COPY "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/TudatFindBoost.cmake"
        DESTINATION "${CMAKE_CURRENT_BINARY_DIR}"
)

install(
        FILES
        "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME_LOWER}-config.cmake"
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/TudatFindBoost.cmake"
        DESTINATION "${INSTALL_CMAKE_DIR}"
)

# Take care of versioning.

# NOTE: SameMinorVersion available only
# since CMake 3.11.
if (${CMAKE_VERSION} VERSION_LESS "3.11.0")
    write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME_LOWER}-config-version.cmake" VERSION ${${PROJECT_NAME}_VERSION}
            COMPATIBILITY SameMajorVersion)
else ()
    write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME_LOWER}-config-version.cmake" VERSION ${${PROJECT_NAME}_VERSION}
            COMPATIBILITY SameMinorVersion)
endif ()
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME_LOWER}-config-version.cmake" DESTINATION "${INSTALL_CMAKE_DIR}")

# ========
message(STATUS "TUDAT_DATA_DIR_RELATIVE_TO_INSTALL_PREFIX: ${TUDAT_DATA_DIR_RELATIVE_TO_INSTALL_PREFIX}")

# +============================================================================
# CLEAN UP (Project name independent)
#  Cleanup in case project is not top level.
# +============================================================================
unset(INSTALL_LIB_DIR)
unset(INSTALL_BIN_DIR)
unset(INSTALL_INCLUDE_DIR)
unset(INSTALL_DATA_DIR)
unset(INSTALL_MAN_DIR)
unset(INSTALL_TESTS_DIR)
unset(INSTALL_CMAKE_DIR)

unset(RELATIVE_INSTALL_LIB_DIR)
unset(RELATIVE_INSTALL_BIN_DIR)
unset(RELATIVE_INSTALL_INCLUDE_DIR)
unset(RELATIVE_INSTALL_DATA_DIR)
unset(RELATIVE_INSTALL_MAN_DIR)
unset(RELATIVE_INSTALL_TESTS_DIR)
unset(RELATIVE_INSTALL_CMAKE_DIR)

