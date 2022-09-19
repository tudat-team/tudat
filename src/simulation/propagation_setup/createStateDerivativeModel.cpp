/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/propagation_setup/createStateDerivativeModel.h"

#include <tudat/simulation/estimation.h>

#include "tudat/astro/gravitation/librationPoint.h"
#include "tudat/astro/gravitation/unitConversionsCircularRestrictedThreeBodyProblem.h"
#include "tudat/math/interpolators/createInterpolator.h"
#include "tudat/astro/basic_astro/celestialBodyConstants.h"

namespace tudat
{

namespace propagators
{

using namespace tudat::simulation_setup;



//! Function to create an integrator to propagate the dynamics (in normalized units) in CR3BP
std::shared_ptr< numerical_integrators::NumericalIntegrator< double, Eigen::Vector6d > > createCR3BPIntegrator(
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const double massParameter,
        const Eigen::Vector6d& initialState )
{
    // Create state derivative model
    std::function< Eigen::Vector6d( const double, const Eigen::Vector6d& ) > stateDerivativeFunction =
            std::bind( &computeCr3bpStateDerivative,
                       std::placeholders::_1, std::placeholders::_2, massParameter );

    // Create integrator object
    return numerical_integrators::createIntegrator< double, Eigen::Vector6d >(
                stateDerivativeFunction, initialState, integratorSettings );
}


//! Function to create an integrator to propagate the dynamics (in normalized units) in CR3BP
std::shared_ptr< numerical_integrators::NumericalIntegrator< double, Eigen::MatrixXd > > createCR3BPWithStmIntegrator(
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const double massParameter,
        const Eigen::Vector6d& initialState )
{
    // Create state derivative model
    std::function< Eigen::MatrixXd( const double, const Eigen::MatrixXd& ) > stateDerivativeFunction =
            std::bind( &computeStateDerivativeWithStateTransitionMatrix,
                       std::placeholders::_1, std::placeholders::_2, massParameter );

    Eigen::MatrixXd fullInitialState = Eigen::MatrixXd::Zero( 6, 7 );
    fullInitialState.block( 0, 0, 6, 1 ) = initialState;
    fullInitialState.block( 0, 1, 6, 6 ).setIdentity( );

    // Create integrator object
    return numerical_integrators::createIntegrator< double, Eigen::MatrixXd >(
                stateDerivativeFunction, fullInitialState, integratorSettings );
}

//! Function to propagate the dynamics (in normalized units) in CR3BP
std::map< double, Eigen::Vector6d > performCR3BPIntegration(
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const double massParameter,
        const Eigen::Vector6d& initialState,
        const double finalTime,
        const bool propagateToExactFinalTime,
        const bool saveFullStateHistory )
{
    // Create integrator object
    std::shared_ptr< numerical_integrators::NumericalIntegrator< double, Eigen::Vector6d > > integrator =
            createCR3BPIntegrator( integratorSettings, massParameter, initialState );

    return propagateCustomDynamics(
            integrator, integratorSettings->initialTimeStep_, finalTime, propagateToExactFinalTime, saveFullStateHistory );

}

std::map< double, Eigen::MatrixXd > performCR3BPWithStmIntegration(
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const double massParameter,
        const Eigen::Vector6d& initialState,
        const double finalTime,
        const bool propagateToExactFinalTime,
        const bool saveFullStateHistory )
{
    // Create integrator object
    std::shared_ptr< numerical_integrators::NumericalIntegrator< double, Eigen::MatrixXd > > integrator =
            createCR3BPWithStmIntegrator( integratorSettings, massParameter, initialState );

    return propagateCustomDynamics(
            integrator, integratorSettings->initialTimeStep_, finalTime, propagateToExactFinalTime, saveFullStateHistory );

}



//template std::vector< std::shared_ptr< SingleStateTypeDerivative< double, double > > > createStateDerivativeModels< double, double >(
//        const std::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings,
//        const simulation_setup::SystemOfBodies& bodies,
//        const double propagationStartTime );
//template std::shared_ptr< SingleStateTypeDerivative< double, double > > createStateDerivativeModel< double, double >(
//        const std::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings,
//        const simulation_setup::SystemOfBodies& bodies,
//        const double propagationStartTime );

//#if( TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )
//template std::vector< std::shared_ptr< SingleStateTypeDerivative< long double, double > > > createStateDerivativeModels< long double, double >(
//        const std::shared_ptr< SingleArcPropagatorSettings< long double > > propagatorSettings,
//        const simulation_setup::SystemOfBodies& bodies,
//        const double propagationStartTime );
//template std::vector< std::shared_ptr< SingleStateTypeDerivative< double, Time > > > createStateDerivativeModels< double, Time >(
//        const std::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings,
//        const simulation_setup::SystemOfBodies& bodies,
//        const Time propagationStartTime );
//template std::vector< std::shared_ptr< SingleStateTypeDerivative< long double, Time > > > createStateDerivativeModels< long double, Time >(
//        const std::shared_ptr< SingleArcPropagatorSettings< long double > > propagatorSettings,
//        const simulation_setup::SystemOfBodies& bodies,
//        const Time propagationStartTime );
//template std::shared_ptr< SingleStateTypeDerivative< long double, double > > createStateDerivativeModel< long double, double >(
//        const std::shared_ptr< SingleArcPropagatorSettings< long double > > propagatorSettings,
//        const simulation_setup::SystemOfBodies& bodies,
//        const double propagationStartTime );
//template std::shared_ptr< SingleStateTypeDerivative< double, Time > > createStateDerivativeModel< double, Time >(
//        const std::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings,
//        const simulation_setup::SystemOfBodies& bodies,
//        const Time propagationStartTime );
//template std::shared_ptr< SingleStateTypeDerivative< long double, Time > > createStateDerivativeModel< long double, Time >(
//        const std::shared_ptr< SingleArcPropagatorSettings< long double > > propagatorSettings,
//        const simulation_setup::SystemOfBodies& bodies,
//        const Time propagationStartTime );
//#endif

}

}
