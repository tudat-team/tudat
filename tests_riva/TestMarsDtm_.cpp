/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *
 */

//#define BOOST_TEST_DYN_LINK
//#define BOOST_TEST_MAIN

#include <limits>
#include "fstream"
#include "iostream"

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/astro/aerodynamics/customAerodynamicCoefficientInterface.h"
#include "tudat/astro/aerodynamics/aerodynamicAcceleration.h"
#include "tudat/astro/reference_frames/aerodynamicAngleCalculator.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"
#include "tudat/interface/spice/spiceEphemeris.h"
#include "tudat/interface/spice/spiceRotationalEphemeris.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/estimation_setup/createNumericalSimulator.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/astro/basic_astro/celestialBodyConstants.h"

using namespace tudat::aerodynamics;

int main( )
{
    using namespace tudat;
    using namespace aerodynamics;
    using namespace simulation_setup;
    using namespace numerical_integrators;
    using namespace simulation_setup;
    using namespace basic_astrodynamics;
    using namespace propagators;
    using namespace basic_mathematics;
    using namespace basic_astrodynamics;
    using namespace orbital_element_conversions;
    using namespace unit_conversions;

    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Create Earth object
    BodyListSettings defaultBodySettings =
            getDefaultBodySettings( { "Mars" } );
    std::string filename ="/Users/ralkahal/OneDrive - Delft University of Technology/PhD/Programs/atmodensitydtm/dtm_mars";;
    defaultBodySettings.at( "Mars" )->atmosphereSettings = marsDtmAtmosphereSettings( filename, 3378.0E3);
    SystemOfBodies bodies = createSystemOfBodies( defaultBodySettings );

    // Create vehicle object.
    double vehicleMass = 400;
    bodies.createEmptyBody( "Vehicle" );
    bodies.at( "Vehicle" )->setConstantBodyMass( vehicleMass );

    // Set aerodynamic coefficients.
    std::shared_ptr<AerodynamicCoefficientSettings> aerodynamicCoefficientSettings =
            std::make_shared<ConstantAerodynamicCoefficientSettings>(
                    2.0, 4.0, Eigen::Vector3d::Zero( ), Eigen::Vector3d::UnitX( ), Eigen::Vector3d::Zero( ),
                    negative_aerodynamic_frame_coefficients, negative_aerodynamic_frame_coefficients );
    bodies.at( "Vehicle" )->setAerodynamicCoefficientInterface(
            createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Vehicle", bodies ));

    // Define acceleration model settings.
    SelectedAccelerationMap accelerationMap;
    std::vector<std::string> bodiesToPropagate;
    std::vector<std::string> centralBodies;
    std::map<std::string, std::vector<std::shared_ptr<AccelerationSettings> > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Mars" ].push_back( std::make_shared<AccelerationSettings>( point_mass_gravity ));
    accelerationsOfVehicle[ "Mars" ].push_back( std::make_shared<AccelerationSettings>( aerodynamic ));
    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;
    bodiesToPropagate.push_back( "Vehicle" );
    centralBodies.push_back( "Mars" );

    // Set initial state
    //Eigen::Vector6d systemInitialState;
    //systemInitialState << 3378.0E3, 3000E3, 4200E3, 0.0, 0.0, 0.0;

    Eigen::Vector6d vehicleInitialStateInKeplerianElements;
    vehicleInitialStateInKeplerianElements( semiMajorAxisIndex ) = 3389.0E3+400.0E3;
    vehicleInitialStateInKeplerianElements( eccentricityIndex ) = 0;
    vehicleInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 85.3 );
    vehicleInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) =
            convertDegreesToRadians( 235.7 );
    vehicleInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) =
            convertDegreesToRadians( 23.4 );
    vehicleInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 139.87 );

    double marsGravitationalParameter = celestial_body_constants::MARS_GRAVITATIONAL_PARAMETER;

    Eigen::VectorXd systemInitialState = convertKeplerianToCartesianElements(
            vehicleInitialStateInKeplerianElements,
            marsGravitationalParameter );

    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
            bodies, accelerationMap, bodiesToPropagate, centralBodies );

    // Set variables to save
    std::vector<std::shared_ptr<SingleDependentVariableSaveSettings> > dependentVariables;
    dependentVariables.push_back(
            std::make_shared<SingleDependentVariableSaveSettings>(
                    altitude_dependent_variable, "Vehicle", "Mars" ));
    dependentVariables.push_back(
            std::make_shared<BodyAerodynamicAngleVariableSaveSettings>(
                    "Vehicle", reference_frames::longitude_angle ));
    dependentVariables.push_back(
            std::make_shared<BodyAerodynamicAngleVariableSaveSettings>(
                    "Vehicle", reference_frames::latitude_angle ));
    dependentVariables.push_back(
            std::make_shared<SingleDependentVariableSaveSettings>(
                    local_density_dependent_variable, "Vehicle", "Mars" ));
    dependentVariables.push_back(
            std::make_shared<SingleDependentVariableSaveSettings>(
                    solar_longitude,"Mars"));
    dependentVariables.push_back(
            std::make_shared<SingleAccelerationDependentVariableSaveSettings>(
                    aerodynamic, "Vehicle", "Mars", 0 ));


    // Set propagation/integration settings
    std::shared_ptr<PropagationTimeTerminationSettings> terminationSettings =
            std::make_shared<propagators::PropagationTimeTerminationSettings>( 365*24*3600.0 );
    std::shared_ptr<IntegratorSettings<> > integratorSettings =
            std::make_shared<IntegratorSettings<> >
                    ( rungeKutta4, 0.0, 60.0 );
    std::shared_ptr<tudat::propagators::TranslationalStatePropagatorSettings<double> >
            translationalPropagatorSettings =
            std::make_shared<TranslationalStatePropagatorSettings < double> >
                    ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, 0.0,
                      integratorSettings, terminationSettings,
                      cowell, dependentVariables );

    std::cout << "Propagation started" << std::endl;
    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator<> dynamicsSimulator(
            bodies, translationalPropagatorSettings );

    std::map< double, Eigen::VectorXd > dependentVariableHistory = dynamicsSimulator.getDependentVariableHistory( );
    //std::map< double, Eigen::VectorXd > stateHistory = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > dependentVariableOutput =
            dynamicsSimulator.getDependentVariableHistory( );

    std::vector<double> times;
    std::vector<double> densityValues;
    std::vector<double> solarLongValues;
    std::cout << "Saving densities started" << std::endl;
    // Iterate over the data
    std::for_each(dependentVariableOutput.begin(), dependentVariableOutput.end(),
                  [&](const auto& it) {
                      // Extract the density value
                      double density = it.second(3); // Assuming density is at index 4
                      double solarLong = it.second(4); // Assuming density is at index 4

                      // Store time and density value
                      //convert time to days
                      times.push_back(it.first/86400.0);

                      densityValues.push_back(density);
                      solarLongValues.push_back(solarLong);

                  });
    std::ofstream outputFile(
            "/Users/ralkahal/OneDrive - Delft University of Technology/PhD/Programs/atmodensitydtm/density_output_test.txt");
    if (!outputFile.is_open()) {
        std::cerr << "Unable to open the file!" << std::endl;
        return 1; // return an error code
    }
    for (int i = 0; i < times.size(); ++i) {
        outputFile << times[i] << " " << densityValues[i] << " " << solarLongValues[i] << std::endl;
    }
    outputFile.close();
    std::cout << "Density computation and output written to file successfully." << std::endl;


}
//
//namespace tudat
//{
//namespace unit_tests
//{
//
//BOOST_AUTO_TEST_SUITE( test_mars_dtm_atmosphere )
//
//BOOST_AUTO_TEST_CASE( testMarsDtmAtmosphere )
//{
//
//}
//
//BOOST_AUTO_TEST_SUITE_END( )
//
//} // namespace unit_tests
//} // namespace tudat
