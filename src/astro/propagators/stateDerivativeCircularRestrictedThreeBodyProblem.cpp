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
 *        Wakker, K.F., "astro I, AE4-874", Delft University of Technology, 2007.
 *
 */

#include <cmath>
#include <iostream>

#include "tudat/basics/utilityMacros.h"

#include "tudat/astro/propagators/stateDerivativeCircularRestrictedThreeBodyProblem.h"
#include "tudat/astro/basic_astro/stateVectorIndices.h"

namespace tudat
{
namespace propagators
{

//! Compute state derivative.
Eigen::Vector6d computeCr3bpStateDerivative(
        const double time, const Eigen::Vector6d& cartesianState, const double massParameter )
{
    using namespace orbital_element_conversions;

    TUDAT_UNUSED_PARAMETER( time );

    // Compute distance to primary body.
    const double xCoordinateToPrimaryBodySquared =
            ( cartesianState( xCartesianPositionIndex ) + massParameter )
            * ( cartesianState( xCartesianPositionIndex ) + massParameter );

    const double yCoordinateSquared = cartesianState( yCartesianPositionIndex )
            * cartesianState( yCartesianPositionIndex );

    const double zCoordinateSquared = cartesianState( zCartesianPositionIndex )
            * cartesianState( zCartesianPositionIndex );

    const double normDistanceToPrimaryBodyCubed = pow(
                xCoordinateToPrimaryBodySquared + yCoordinateSquared + zCoordinateSquared, 1.5 );

    // Compute distance to secondary body.
    const double xCoordinateSecondaryBodySquared =
            ( cartesianState( xCartesianPositionIndex ) - ( 1.0 - massParameter ) )
            * ( cartesianState( xCartesianPositionIndex ) - ( 1.0 - massParameter ) );

    double normDistanceToSecondaryBodyCubed = pow(
                xCoordinateSecondaryBodySquared + yCoordinateSquared + zCoordinateSquared, 1.5 );

    // Compute derivative of state.
    Eigen::VectorXd stateDerivative( 6 );

    stateDerivative.segment( xCartesianPositionIndex, 3 ) = cartesianState.segment( xCartesianVelocityIndex, 3 );

    stateDerivative( 3 ) = cartesianState( xCartesianPositionIndex )
            - ( ( 1.0 - massParameter ) / normDistanceToPrimaryBodyCubed )
            * ( cartesianState( xCartesianPositionIndex ) + massParameter )
            - ( massParameter / normDistanceToSecondaryBodyCubed )
            * ( cartesianState( xCartesianPositionIndex ) - ( 1.0 - massParameter ) )
            + 2.0 * cartesianState( yCartesianVelocityIndex );
    stateDerivative( 4 ) = cartesianState( yCartesianPositionIndex )
            * ( 1.0 - ( ( 1.0 - massParameter ) / normDistanceToPrimaryBodyCubed )
                - ( massParameter / normDistanceToSecondaryBodyCubed ) )
            - 2.0 * cartesianState( xCartesianVelocityIndex );
    stateDerivative( 5 ) = -cartesianState( zCartesianPositionIndex )
            * ( ( ( 1.0 - massParameter ) / normDistanceToPrimaryBodyCubed )
                + ( massParameter / normDistanceToSecondaryBodyCubed ) );

    // Return computed state derivative.
    return stateDerivative;
}


Eigen::MatrixXd computeStateDerivativeWithStateTransitionMatrix(
        const double time, const Eigen::MatrixXd& cartesianState, const double massParameter )
{
//    std::cout<<"STM and state derivative : "<<time<<" "<<massParameter<<std::endl<<
//               cartesianState<<" "<<std::endl<<std::endl;
    // Declare state derivative vector with same length as the state.
    Eigen::MatrixXd stateDerivative = Eigen::MatrixXd::Zero( 6, 7 );

    // Set the derivative of the position equal to the velocities.
    stateDerivative.block< 3, 1 >( 0, 0 ) = cartesianState.block( 3, 0, 3, 1 );

    double xPositionScaledSquared = ( cartesianState( 0 ) + massParameter ) * ( cartesianState( 0 ) + massParameter );
    double xPositionScaledSquared2 = ( 1.0 - massParameter - cartesianState( 0 ) ) *
            ( 1.0 - massParameter - cartesianState( 0 ) );
    double yPositionScaledSquared = ( cartesianState( 1 ) * cartesianState( 1 ) );
    double zPositionScaledSquared = ( cartesianState( 2 ) * cartesianState( 2 ) );

    // Compute distances to primaries.
    double distanceToPrimaryBody = sqrt(xPositionScaledSquared + yPositionScaledSquared + zPositionScaledSquared);
    double distanceToSecondaryBody = sqrt(xPositionScaledSquared2 + yPositionScaledSquared + zPositionScaledSquared);

    double distanceToPrimaryCubed = distanceToPrimaryBody * distanceToPrimaryBody * distanceToPrimaryBody;
    double distanceToSecondaryCubed = distanceToSecondaryBody * distanceToSecondaryBody * distanceToSecondaryBody;

    double distanceToPrimaryToFifthPower = distanceToPrimaryCubed * distanceToPrimaryBody * distanceToPrimaryBody;
    double distanceToSecondaryToFifthPower = distanceToSecondaryCubed * distanceToSecondaryBody * distanceToSecondaryBody;

    // Set the derivative of the velocities to the accelerations.
    double termRelatedToPrimaryBody   = ( 1.0 - massParameter ) / distanceToPrimaryCubed;
    double termRelatedToSecondaryBody = massParameter / distanceToSecondaryCubed;
    stateDerivative( 3, 0 ) = -termRelatedToPrimaryBody * ( massParameter+cartesianState( 0 ) ) +
            termRelatedToSecondaryBody * ( 1.0 - massParameter - cartesianState( 0 ) ) +
            cartesianState( 0 ) + 2.0 *cartesianState( 4 );
    stateDerivative( 4, 0 ) = -termRelatedToPrimaryBody *cartesianState( 1 ) -
            termRelatedToSecondaryBody *cartesianState( 1 ) + cartesianState( 1 ) - 2.0 *cartesianState( 3 );
    stateDerivative( 5, 0 ) = -termRelatedToPrimaryBody *cartesianState( 2 ) -
            termRelatedToSecondaryBody *cartesianState( 2 );

    // Compute partial derivatives of the potential.
    double Uxx = ( 3.0 * ( 1.0 - massParameter ) *xPositionScaledSquared ) / distanceToPrimaryToFifthPower+
            ( 3.0  *massParameter * xPositionScaledSquared2 ) / distanceToSecondaryToFifthPower -
            ( 1.0 - massParameter ) / distanceToPrimaryCubed - massParameter/distanceToSecondaryCubed + 1.0;
    double Uxy = ( 3.0 * ( 1.0 - massParameter ) * ( cartesianState( 0 ) + massParameter ) *cartesianState( 1 ) ) /
            distanceToPrimaryToFifthPower - ( 3.0 *massParameter * ( 1.0 - massParameter - cartesianState( 0 ) ) *
                                              cartesianState( 1 ) ) / distanceToSecondaryToFifthPower;
    double Uxz = ( 3.0 * ( 1.0 - massParameter ) * ( cartesianState( 0 ) + massParameter ) *cartesianState( 2 ) ) /
            distanceToPrimaryToFifthPower - ( 3.0 *massParameter * ( 1.0 - massParameter - cartesianState( 0 ) ) *
                                              cartesianState( 2 ) ) / distanceToSecondaryToFifthPower;
    double Uyx = Uxy;
    double Uyy = ( 3.0 * ( 1.0 - massParameter ) * yPositionScaledSquared ) / distanceToPrimaryToFifthPower +
            ( 3.0 *massParameter * yPositionScaledSquared ) / distanceToSecondaryToFifthPower -
            ( 1.0 - massParameter ) / distanceToPrimaryCubed - massParameter/distanceToSecondaryCubed + 1.0 ;
    double Uyz = ( 3.0 * ( 1.0 - massParameter ) * cartesianState( 1 ) *cartesianState( 2 ) ) / distanceToPrimaryToFifthPower +
            ( 3.0 *massParameter * cartesianState( 1 ) *cartesianState( 2 ) ) / distanceToSecondaryToFifthPower;
    double Uzx = Uxz;
    double Uzy = Uyz;
    double Uzz = ( 3.0 * ( 1.0 - massParameter ) *zPositionScaledSquared ) / distanceToPrimaryToFifthPower +
            ( 3.0 *massParameter * zPositionScaledSquared ) / distanceToSecondaryToFifthPower -
            ( 1.0 - massParameter ) / distanceToPrimaryCubed - massParameter/distanceToSecondaryCubed ;

    Eigen::MatrixXd stateTransitionDerivative_ = Eigen::MatrixXd::Zero( 6, 6 );
    stateTransitionDerivative_.block( 0, 3, 3, 3 ) = Eigen::Matrix3d::Identity( );
    stateTransitionDerivative_( 3, 4 ) = 2.0;
    stateTransitionDerivative_( 4, 3 ) = -2.0;
    stateTransitionDerivative_.block< 3, 3 >( 3, 0 ) << Uxx, Uxy, Uxz,
                                                      Uyx, Uyy, Uyz,
                                                      Uzx, Uzy, Uzz;
    // Differentiate the STM.
    stateDerivative.block( 0, 1, 6, 6 ) = stateTransitionDerivative_ * cartesianState.block( 0, 1, 6, 6 );

//    std::cout<<stateDerivative<<std::endl<<std::endl<<std::endl;

    return stateDerivative;

}



} // namespace propagators

} // namespace tudat
