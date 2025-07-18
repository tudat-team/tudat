/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/acceleration_partials/rtgAccelerationPartial.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/rtgForceVector.h"

namespace tudat
{

namespace acceleration_partials
{
/*
using namespace gravitation;


//! Function determine the numerical partial derivative of the true anomaly wrt the elements of the Cartesian state
Eigen::Matrix< double, 1, 6 > calculateNumericalPartialOfTrueAnomalyWrtState( const Eigen::Vector6d& cartesianElements,
                                                                              const double gravitationalParameter,
                                                                              const Eigen::Vector6d& cartesianStateElementPerturbations )
{
    using namespace tudat::orbital_element_conversions;

    // Initialize partial to zero
    Eigen::Matrix< double, 1, 6 > partial = Eigen::Matrix< double, 1, 6 >::Zero( );

    // Iterate over all six elements and calculate partials.
    Eigen::Vector6d perturbedCartesianElements;
    double upPerturbedTrueAnomaly, downPerturbedTrueAnomaly;
    for( int i = 0; i < 6; i++ )
    {
        // Calculate true anomaly at up-perturbed cartesian element i
        perturbedCartesianElements = cartesianElements;
        perturbedCartesianElements( i ) += cartesianStateElementPerturbations( i );
        upPerturbedTrueAnomaly =
                convertCartesianToKeplerianElements( perturbedCartesianElements, gravitationalParameter )( trueAnomalyIndex );

        // Check validity of result
        if( upPerturbedTrueAnomaly != upPerturbedTrueAnomaly )
        {
            std::cerr << "Error 1 in partial of true anomaly wrt cartesian state, element" << i << ", keplerian state: " << std::endl
                      << convertCartesianToKeplerianElements( perturbedCartesianElements, gravitationalParameter ).transpose( ) << std::endl
                      << perturbedCartesianElements.transpose( ) << std::endl;
        }

        // Calculate true anomaly at down-perturbed cartesian element i
        perturbedCartesianElements = cartesianElements;
        perturbedCartesianElements( i ) -= cartesianStateElementPerturbations( i );
        downPerturbedTrueAnomaly =
                convertCartesianToKeplerianElements( perturbedCartesianElements, gravitationalParameter )( trueAnomalyIndex );

        // Check validity of result
        if( downPerturbedTrueAnomaly != downPerturbedTrueAnomaly )
        {
            std::cerr << "Error 2 in partial of true anomaly wrt cartesian state, element" << i << ", keplerian state: " << std::endl
                      << convertCartesianToKeplerianElements( perturbedCartesianElements, gravitationalParameter ).transpose( ) << std::endl
                      << perturbedCartesianElements.transpose( ) << std::endl;
        }

        // Calculate central difference of term i
        partial( 0, i ) = ( upPerturbedTrueAnomaly - downPerturbedTrueAnomaly ) / ( 2.0 * cartesianStateElementPerturbations( i ) );

        // Check validity of result
        if( partial( 0, i ) != partial( 0, i ) )
        {
            std::cerr << "Error 2 in partial of true anomaly wrt cartesian state, element" << i << std::endl;
        }
    }

    return partial;
}

 */

//! Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
std::pair< std::function< void( Eigen::MatrixXd& ) >, int > RTGAccelerationPartial::getParameterPartialFunction(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
{
    using namespace tudat::estimatable_parameters;

    std::function< void( Eigen::MatrixXd& ) > partialFunction;
    int numberOfRows = 0;

    if( parameter->getParameterName( ).second.first == acceleratedBody_ )
    {
        switch( parameter->getParameterName( ).first )
        {
            case rtg_force_vector: {
                if( parameter->getParameterName( ).second.second == acceleratingBody_ )
                {
                    partialFunction =
                            std::bind( &RTGAccelerationPartial::wrtRTGForceVector,
                                       this,
                                       std::placeholders::_1 );
                    numberOfRows = parameter->getParameterSize( );
                }
                break;
            }
            default:
                break;
        }
    }

    return std::make_pair( partialFunction, numberOfRows );
}


//! Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
std::pair< std::function< void( Eigen::MatrixXd& ) >, int > RTGAccelerationPartial::getParameterPartialFunction(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
{
    using namespace tudat::estimatable_parameters;

    std::function< void( Eigen::MatrixXd& ) > partialFunction;
    int numberOfRows = 0;

    if( parameter->getParameterName( ).second.first == acceleratedBody_ )
    {
        switch( parameter->getParameterName( ).first )
        {
            case rtg_force_vector_magnitude: {
                if( parameter->getParameterName( ).second.second == acceleratingBody_ )
                {
                    partialFunction =
                            std::bind( &RTGAccelerationPartial::wrtRTGForceVectorMagnitude,
                                       this,
                                       std::placeholders::_1 );
                    numberOfRows = parameter->getParameterSize( );
                }
                break;
            }
            default:
                break;
        }
    }

    return std::make_pair( partialFunction, numberOfRows );
}

//! Function for updating common blocks of partial to current state.
void RTGAccelerationPartial::update( const double currentTime )
{
    if( !( currentTime_ == currentTime ) )
    {
        using namespace tudat::basic_mathematics;
        using namespace tudat::linear_algebra;

        // reduced to doing nothing, since partials are const 0

        /*
        // Get current state and associated data.
        Eigen::Vector6d currentState = rtgAcceleration_->getCurrentState( );
        Eigen::Vector3d angularMomentumVector =
                Eigen::Vector3d( currentState.segment( 0, 3 ) ).cross( Eigen::Vector3d( currentState.segment( 3, 3 ) ) );
        Eigen::Vector3d crossVector = angularMomentumVector.cross( Eigen::Vector3d( currentState.segment( 0, 3 ) ) );

        // Compute constituent geometric partials.
        Eigen::Matrix3d normPositionWrtPosition =
                calculatePartialOfNormalizedVector( Eigen::Matrix3d::Identity( ), currentState.segment( 0, 3 ) );
        Eigen::Matrix3d angularMomentumWrtPosition = -getCrossProductMatrix( currentState.segment( 3, 3 ) );
        Eigen::Matrix3d angularMomentumWrtVelocity = getCrossProductMatrix( currentState.segment( 0, 3 ) );
        Eigen::Matrix3d normAngularMomentumWrtPosition =
                calculatePartialOfNormalizedVector( angularMomentumWrtPosition, angularMomentumVector );
        Eigen::Matrix3d normAngularMomentumWrtVelocity =
                calculatePartialOfNormalizedVector( angularMomentumWrtVelocity, angularMomentumVector );
        Eigen::Matrix3d crossVectorWrtPosition =
                getCrossProductMatrix( angularMomentumVector ) - angularMomentumWrtVelocity * angularMomentumWrtPosition;
        Eigen::Matrix3d crossVectorWrtVelocity =
                -getCrossProductMatrix( currentState.segment( 0, 3 ) ) * getCrossProductMatrix( currentState.segment( 0, 3 ) );
        Eigen::Matrix3d normCrossVectorWrtPosition = calculatePartialOfNormalizedVector( crossVectorWrtPosition, crossVector );
        Eigen::Matrix3d normCrossVectorWrtVelocity = calculatePartialOfNormalizedVector( crossVectorWrtVelocity, crossVector );

        // Retrieve current empirical acceleration
        Eigen::Vector3d localAcceleration = empiricalAcceleration_->getCurrentLocalAcceleration( );

        // Compute partial derivative contribution of derivative of rotation matrix from RSW to inertial frame
        currentPositionPartial_ = localAcceleration.x( ) * normPositionWrtPosition + localAcceleration.y( ) * normCrossVectorWrtPosition +
                localAcceleration.z( ) * normAngularMomentumWrtPosition;
        currentVelocityPartial_ =
                localAcceleration.y( ) * normCrossVectorWrtVelocity + localAcceleration.z( ) * normAngularMomentumWrtVelocity;


        // Compute partial derivative contribution of derivative of true anomaly
        Eigen::Matrix< double, 1, 6 > localTrueAnomalyPartial =
                calculateNumericalPartialOfTrueAnomalyWrtState( empiricalAcceleration_->getCurrentState( ),
                                                                empiricalAcceleration_->getCurrentGravitationalParameter( ),
                                                                cartesianStateElementPerturbations );
        Eigen::Matrix< double, 1, 6 > trueAnomalyPartial = localTrueAnomalyPartial;
        currentPositionPartial_ += empiricalAcceleration_->getCurrentToInertialFrame( ) *
                ( empiricalAcceleration_->getCurrentAccelerationComponent( basic_astrodynamics::sine_empirical ) *
                          std::cos( empiricalAcceleration_->getCurrentTrueAnomaly( ) ) -
                  empiricalAcceleration_->getCurrentAccelerationComponent( basic_astrodynamics::cosine_empirical ) *
                          std::sin( empiricalAcceleration_->getCurrentTrueAnomaly( ) ) ) *
                trueAnomalyPartial.block( 0, 0, 1, 3 );
        currentVelocityPartial_ += ( empiricalAcceleration_->getCurrentToInertialFrame( ) *
                                     ( empiricalAcceleration_->getCurrentAccelerationComponent( basic_astrodynamics::sine_empirical ) *
                                               std::cos( empiricalAcceleration_->getCurrentTrueAnomaly( ) ) -
                                       empiricalAcceleration_->getCurrentAccelerationComponent( basic_astrodynamics::cosine_empirical ) *
                                               std::sin( empiricalAcceleration_->getCurrentTrueAnomaly( ) ) ) ) *
                trueAnomalyPartial.block( 0, 3, 1, 3 );
         */
        currentTime_ = currentTime;

        // Check output.
        if( currentPositionPartial_ != currentPositionPartial_ )
        {
            throw std::runtime_error( "Error, empirical acceleration position partials are NaN" );
        }
        if( currentVelocityPartial_ != currentVelocityPartial_ )
        {
            throw std::runtime_error( "Error, empirical acceleration velocity partials are NaN" );
        }
    }
}

//! Function to compute the partial w.r.t. reference force vector
void RTGAccelerationPartial::wrtRTGForceVector( Eigen::MatrixXd& partialDerivativeMatrix )
{
    // Compute partial derivative w.r.t. reference force vector
    double partialWrtReferenceForceVector;
    partialWrtReferenceForceVector = rtgAcceleration_->getCurrentDecayTerm( );
    partialDerivativeMatrix = partialWrtReferenceForceVector * Eigen::MatrixXd::Identity( 3, 3 );

}

//! Function to compute the partial w.r.t. magnitude of reference force vector
void RTGAccelerationPartial::wrtRTGForceVectorMagnitude(Eigen::MatrixXd& partialDerivativeMatrix )
{
    // Compute partial derivative w.r.t. magnitude of reference force vector
    Eigen::Vector3d partialWrtReferenceForceMagnitude;
    partialWrtReferenceForceMagnitude = rtgAcceleration_->getCurrentDecayTerm( ) * rtgAcceleration_->getBodyFixedForceUnitVectorAtReferenceEpoch( ) ;
    partialDerivativeMatrix = partialWrtReferenceForceMagnitude;

}

}  // namespace acceleration_partials

}  // namespace tudat