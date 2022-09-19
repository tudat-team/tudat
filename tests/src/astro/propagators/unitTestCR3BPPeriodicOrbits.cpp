/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

//#define BOOST_TEST_DYN_LINK
//#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>
#include <boost/bind/bind.hpp>
using namespace std::placeholders;

#include <boost/make_shared.hpp>
#include <memory>

#include "tudat/io/basicInputOutput.h"
#include "tudat/simulation/propagation_setup/createCR3BPPeriodicOrbits.h"
#include "tudat/math/integrators/createNumericalIntegrator.h"
#include "tudat/astro/gravitation/librationPoint.h"

//namespace tudat
//{

//namespace unit_tests
//{

using namespace tudat;
using namespace propagators;
using namespace numerical_integrators;
using namespace circular_restricted_three_body_problem;

//BOOST_AUTO_TEST_SUITE( test_cr3bp_periodic_orbits )

//BOOST_AUTO_TEST_CASE( testCr3bpPeriodicOrbits )
int main( )
{
    double minimumStepSize   = std::numeric_limits<double>::epsilon( ); // 2.22044604925031e-16
    double maximumStepSize   = 100.0;//std::numeric_limits<double>::infinity( ); // 2.22044604925031e-16

    const double relativeErrorTolerance = 1.0E-10;
    const double absoluteErrorTolerance = 1.0E-14;

    std::shared_ptr< IntegratorSettings< double > > integratorSettings =
                        std::make_shared< RungeKuttaVariableStepSizeSettings< > >
                        ( 0.0, 1.0E-5,
                          rungeKutta87DormandPrince, minimumStepSize, maximumStepSize,
                          relativeErrorTolerance, absoluteErrorTolerance );


    double massParameter = computeMassParameter( 3.986004418E14, 1.32712440018e20 / ( 328900.56 * ( 1.0 + 81.30059 ) ) );

    int maximumDifferentialCorrections = 25;
    double maximumPositionDeviation = 1.0E-12;
    double maximumVelocityDeviation = 1.0E-12;

    int maximumNumberOfOrbits = 10000;
    double maximumEigenvalueDeviation = 1.0E-3;

    int librationPointNumber = 1;
    CR3BPPeriodicOrbitTypes orbitType = halo_orbit;
    CR3BPPeriodicOrbitGenerationSettings orbitSettings(
                massParameter, orbitType, librationPointNumber,maximumDifferentialCorrections,
                maximumPositionDeviation, maximumVelocityDeviation, maximumNumberOfOrbits, maximumEigenvalueDeviation );

    double amplitudeFirstGuess;
    std::pair< Eigen::Vector6d, double > periodicOrbitInitialGuess;
    std::vector< CR3BPPeriodicOrbitConditions > periodicOrbits;

    for( int i = 0; i < 2; i ++ )
    {
        amplitudeFirstGuess = initializeEarthMoonPeriodicOrbitAmplitude(
                    librationPointNumber, orbitType, i );
        periodicOrbitInitialGuess = richardsonApproximationLibrationPointPeriodicOrbit(
                    massParameter, orbitType, librationPointNumber, amplitudeFirstGuess );
        periodicOrbits.push_back( createCR3BPPeriodicOrbit(
                                      periodicOrbitInitialGuess.first, periodicOrbitInitialGuess.second,
                                      orbitSettings, integratorSettings ) );
    }

    createCR3BPPeriodicOrbitsThroughNumericalContinuation(
            periodicOrbits, integratorSettings, orbitSettings );
    std::cout<<periodicOrbits.size( )<<std::endl;

    std::string outputFolder = "/home/dominic/Software/periodicOrbitResults/";
    for( unsigned int i = 0; i < periodicOrbits.size( ); i++ )
    {

        std::map< double, Eigen::Vector6d > currentOrbit = propagatePeriodicOrbit(
                periodicOrbits.at( i ), integratorSettings );
        input_output::writeDataMapToTextFile(
                   currentOrbit, "periodic_orbit_" + std::to_string( i ) + ".dat", outputFolder );

    }

}

//BOOST_AUTO_TEST_SUITE_END( )

//}

//}
