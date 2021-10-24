#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#
# Add source files.
set(numerical_simulation_HEADERS
        helpers.h
        BaseSimulator.h
        SingleArcSimulator.h
        BaseObserver.h
        ComputationExpense.h
        )

# Add header files.
set(numerical_simulation_SOURCES
        BaseSimulator.cpp
        SingleArcSimulator.cpp
        BaseObserver.cpp
        ComputationExpense.cpp
        )

# Add library.
TUDAT_ADD_LIBRARY("numerical_simulation"
        "${numerical_simulation_SOURCES}"
        "${numerical_simulation_HEADERS}"
        )

add_subdirectory(environment_setup)
add_subdirectory(propagation_setup)
if (TUDAT_BUILD_WITH_ESTIMATION_TOOLS)
    add_subdirectory(estimation_setup)
endif ()


