/* git    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <string>
#include <thread>

#include <boost/test/unit_test.hpp>


#include "tudat/basics/testMacros.h"
#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/unitConversions.h"

#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/math/integrators/rungeKuttaCoefficients.h"
#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/astro/ephemerides/keplerEphemeris.h"

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/estimation_setup/variationalEquationsSolver.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createBodies.h"
#include "tudat/simulation/environment_setup/createSystemModel.h"
#include "tudat/simulation/estimation_setup/createNumericalSimulator.h"
#include "tudat/simulation/estimation_setup/createEstimatableParameters.h"

namespace tudat
{

namespace unit_tests
{

//Using declarations.
using namespace tudat;
using namespace tudat::estimatable_parameters;
using namespace tudat::orbit_determination;
using namespace tudat::interpolators;
using namespace tudat::numerical_integrators;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::basic_astrodynamics;
using namespace tudat::orbital_element_conversions;
using namespace tudat::ephemerides;
using namespace tudat::propagators;

BOOST_AUTO_TEST_SUITE( test_variational_equation_calculation )

BOOST_AUTO_TEST_CASE( testStateParameterBodiesDifference )
{

    spice_interface::loadStandardSpiceKernels( );

    std::map< double, Eigen::MatrixXd > fullStateTransitionMatrixHistory;
    std::map< double, Eigen::VectorXd > fullStateHistory;

    for( int test = 0; test < 3; test++ )
    {
        // Define bodies in simulation
        std::vector< std::string > bodyNames;
        bodyNames.push_back( "Earth" );
        bodyNames.push_back( "Sun" );
        bodyNames.push_back( "Moon" );
        bodyNames.push_back( "Mars" );

        // Specify initial time
        double initialEphemerisTime = double( 1.0E7 );
        double finalEphemerisTime = initialEphemerisTime + 4.0 * 3600.0;

        // Create bodies needed in simulation
        BodyListSettings bodySettings =
                getDefaultBodySettings( bodyNames );
        SystemOfBodies bodies = createSystemOfBodies( bodySettings );
        bodies.createEmptyBody( "Vehicle1" );
        bodies.createEmptyBody( "Vehicle2" );
        bodies.at( "Vehicle1" )->setConstantBodyMass( 400.0 );
        bodies.at( "Vehicle2" )->setConstantBodyMass( 100.0 );


        // Create radiation pressure settings
        double referenceAreaRadiation = 4.0;
        double radiationPressureCoefficient = 1.2;
        std::vector< std::string > occultingBodies;
        occultingBodies.push_back( "Earth" );
        std::shared_ptr< RadiationPressureInterfaceSettings > asterixRadiationPressureSettings =
                std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                    "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

        // Create and set radiation pressure settings
        bodies.at( "Vehicle1" )->setRadiationPressureInterface(
                    "Sun", createRadiationPressureInterface(
                        asterixRadiationPressureSettings, "Vehicle1", bodies ) );
        bodies.at( "Vehicle2" )->setRadiationPressureInterface(
            "Sun", createRadiationPressureInterface(
                asterixRadiationPressureSettings, "Vehicle2", bodies ) );

        // Set accelerations on Vehicle that are to be taken into account.
        SelectedAccelerationMap accelerationMap;
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
        accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 2, 2 ) );

        accelerationsOfVehicle[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                       basic_astrodynamics::point_mass_gravity ) );
        accelerationsOfVehicle[ "Moon" ].push_back( std::make_shared< AccelerationSettings >(
                                                        basic_astrodynamics::point_mass_gravity ) );
        accelerationsOfVehicle[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                       basic_astrodynamics::cannon_ball_radiation_pressure ) );
        accelerationMap[ "Vehicle1" ] = accelerationsOfVehicle;
        accelerationMap[ "Vehicle2" ] = accelerationsOfVehicle;


        // Set bodies for which initial state is to be estimated and integrated.
        std::vector< std::string > bodiesToIntegrate;
        std::vector< std::string > centralBodies;
        bodiesToIntegrate.push_back( "Vehicle1" );
        bodiesToIntegrate.push_back( "Vehicle2" );

        centralBodies.push_back( "Earth" );
        centralBodies.push_back( "Earth" );

        // Create acceleration models
        AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodies, accelerationMap, bodiesToIntegrate, centralBodies );

        // Create integrator settings
        std::shared_ptr< IntegratorSettings< double > > integratorSettings =
                std::make_shared< IntegratorSettings< double > >
                ( rungeKutta4, initialEphemerisTime, 15.0 );

        // Set Keplerian elements for Asterix.
        Eigen::Vector6d asterixInitialStateInKeplerianElements;
        asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7500.0E3;
        asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
        asterixInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
        asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
                = unit_conversions::convertDegreesToRadians( 235.7 );
        asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
                = unit_conversions::convertDegreesToRadians( 23.4 );
        asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

        double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );

        // Set (perturbed) initial state.
        Eigen::Vector6d initialTranslationalState1 = convertKeplerianToCartesianElements(
                    asterixInitialStateInKeplerianElements, earthGravitationalParameter );
        asterixInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 25.3 );
        Eigen::Vector6d initialTranslationalState2 = convertKeplerianToCartesianElements(
            asterixInitialStateInKeplerianElements, earthGravitationalParameter );
        Eigen::VectorXd initialTranslationalState = Eigen::VectorXd::Zero( 12 );
        initialTranslationalState << initialTranslationalState1, initialTranslationalState2;

        // Create propagator settings
        std::shared_ptr< TranslationalStatePropagatorSettings< double, double > > propagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double, double > >
                ( centralBodies, accelerationModelMap, bodiesToIntegrate, initialTranslationalState, initialEphemerisTime, integratorSettings,
                  propagationTimeTerminationSettings( finalEphemerisTime ) );
        propagatorSettings->getPrintSettings( )->setPrintInitialAndFinalConditions( true );

        // Define parameters.
        std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
            getInitialStateParameterSettings< double >( propagatorSettings, bodies );
        if( test > 0 )
        {
            parameterNames = { parameterNames.at( test - 1 ) };
        }
        // Create parameters
        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
                createParametersToEstimate( parameterNames, bodies );


        SingleArcVariationalEquationsSolver< double, double > dynamicsSimulator =
                SingleArcVariationalEquationsSolver< double, double >(
                    bodies, propagatorSettings, parametersToEstimate );

        if( test == 0 )
        {
            fullStateTransitionMatrixHistory = dynamicsSimulator.getSingleArcVariationalPropagationResults( )->getStateTransitionSolution( );
            fullStateHistory = dynamicsSimulator.getSingleArcVariationalPropagationResults( )->getDynamicsResults( )->getEquationsOfMotionNumericalSolution( );
        }
        else
        {
            std::map< double, Eigen::MatrixXd > partialStateTransitionMatrixHistory = dynamicsSimulator.getSingleArcVariationalPropagationResults( )->getStateTransitionSolution( );
            std::map< double, Eigen::VectorXd > partialStateHistory = dynamicsSimulator.getSingleArcVariationalPropagationResults( )->getDynamicsResults( )->getEquationsOfMotionNumericalSolution( );
            for( auto it : fullStateTransitionMatrixHistory )
            {
                BOOST_CHECK( partialStateTransitionMatrixHistory.count( it.first ) > 0 );
                BOOST_CHECK( partialStateHistory.count( it.first ) > 0 );
                BOOST_CHECK( fullStateHistory.count( it.first ) > 0 );

                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ( fullStateHistory.at( it.first ) ), ( partialStateHistory.at( it.first ) ), std::numeric_limits<double >::epsilon( ) );
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ( fullStateTransitionMatrixHistory.at( it.first ).block( 6 * ( test - 1 ), 6 * ( test - 1 ), 6, 6 ) ), ( partialStateTransitionMatrixHistory.at( it.first ) ), std::numeric_limits<double >::epsilon( ) );

            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}

}

