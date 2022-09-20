/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/gravitation/librationPointFunctions.h"
#include "tudat/simulation/propagation_setup/createCR3BPPeriodicOrbits.h"

namespace tudat
{

namespace propagators
{

std::string getPeriodicOrbitName( enum CR3BPPeriodicOrbitTypes orbitType )
{
    std::string orbitName;
    switch( orbitType )
    {
    case horizontal_lyapunov_orbit:
        orbitName = "HL";
        break;
    case vertical_lyapunov_orbit:
        orbitName = "VL";
        break;
    case halo_orbit:
        orbitName = "Halo";
        break;
    case axial_orbit:
        orbitName = "Axials";
        break;
    default:
        throw std::runtime_error( "Error when getting periodic orbit type name, id not recognized" );
    }
    return orbitName;
}

std::map< double, Eigen::Vector6d > propagatePeriodicOrbit(
        const CR3BPPeriodicOrbitConditions& orbitDefinition,
        const std::shared_ptr< tudat::numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const double numberOfPeriods )
{
    return performCR3BPIntegration(
                integratorSettings, orbitDefinition.massParameter( ), orbitDefinition.initialState_,
                orbitDefinition.orbitPeriod_ * numberOfPeriods, true, true );
}


std::pair< Eigen::Vector6d, double >  richardsonApproximationLibrationPointPeriodicOrbit(
        const double massParameter,
        const CR3BPPeriodicOrbitTypes orbitType,
        int librationPointNr, double amplitude, double n )
{
    using namespace circular_restricted_three_body_problem;

    Eigen::Vector6d initialStateVector = Eigen::Vector6d::Zero( );

    double gammaL;
    double c2;
    double c3;
    double c4;
    double Ax = 0.0;
    double Az = 0.0;

    if (librationPointNr == 1)
    {
        // Create object containing the functions.
        std::shared_ptr< LibrationPointLocationFunction1 > LibrationPointLocationFunction =
                std::make_shared< LibrationPointLocationFunction1 >( 1, massParameter );

        // The termination condition.
        tudat::root_finders::NewtonRaphson< >::TerminationFunction terminationConditionFunction =
                std::bind( &tudat::root_finders::RootAbsoluteToleranceTerminationCondition< double >::checkTerminationCondition,
                           std::make_shared< tudat::root_finders::RootAbsoluteToleranceTerminationCondition< double > >(
                               LibrationPointLocationFunction->getTrueRootAccuracy( ) ),
                           std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5 );

        // Test Newton-Raphson object.
        tudat::root_finders::NewtonRaphson< > newtonRaphson( terminationConditionFunction );

        // Let Newton-Raphson search for the root.
        gammaL = newtonRaphson.execute( LibrationPointLocationFunction, LibrationPointLocationFunction->getInitialGuess( ) );

        c2 = 1.0 / pow(gammaL, 3.0) * (pow(1.0,2.0) * massParameter + pow(-1.0,2.0) * (1.0 - massParameter) * pow(gammaL, 2.0+1.0) / pow((1.0 - gammaL), (2.0+1.0)));
        c3 = 1.0 / pow(gammaL, 3.0) * (pow(1.0,3.0) * massParameter + pow(-1.0,3.0) * (1.0 - massParameter) * pow(gammaL, 3.0+1.0) / pow((1.0 - gammaL), (3.0+1.0)));
        c4 = 1.0 / pow(gammaL, 3.0) * (pow(1.0,4.0) * massParameter + pow(-1.0,4.0) * (1.0 - massParameter) * pow(gammaL, 4.0+1.0) / pow((1.0 - gammaL), (4.0+1.0)));
    }
    else
    {
        // Create object containing the functions.
        std::shared_ptr< LibrationPointLocationFunction2 > LibrationPointLocationFunction =
                std::make_shared< LibrationPointLocationFunction2 >( 1, massParameter );

        // The termination condition.
        tudat::root_finders::NewtonRaphson< >::TerminationFunction terminationConditionFunction =
                std::bind( &tudat::root_finders::RootAbsoluteToleranceTerminationCondition< double >::checkTerminationCondition,
                           std::make_shared< tudat::root_finders::RootAbsoluteToleranceTerminationCondition< double > >(
                               LibrationPointLocationFunction->getTrueRootAccuracy( ) ),
                           std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5 );
        // Test Newton-Raphson object.
        tudat::root_finders::NewtonRaphson< > newtonRaphson( terminationConditionFunction );

        // Let Newton-Raphson search for the root.
        gammaL = newtonRaphson.execute( LibrationPointLocationFunction, LibrationPointLocationFunction->getInitialGuess( ) );

        c2 = 1 / pow(gammaL, 3.0) * (pow(-1.0, 2.0) * massParameter + pow(-1.0, 2.0) * (1.0 - massParameter) * pow(gammaL, 2.0+1.0) / pow((1.0 + gammaL), (2.0 + 1.0)));
        c3 = 1 / pow(gammaL, 3.0) * (pow(-1.0, 3.0) * massParameter + pow(-1.0, 3.0) * (1.0 - massParameter) * pow(gammaL, 3.0+1.0) / pow((1.0 + gammaL), (3.0 + 1.0)));
        c4 = 1 / pow(gammaL, 3.0) * (pow(-1.0, 4.0) * massParameter + pow(-1.0, 4.0) * (1.0 - massParameter) * pow(gammaL, 4.0+1.0) / pow((1.0 + gammaL), (4.0 + 1.0)));
    }

    double lambda = pow(1.0 - c2/2.0 + 1.0/2.0*pow((pow((c2-2.0), 2.0) + 4.0*(c2-1.0)*(1.0+2.0*c2)), 0.5), 0.5);
    double k      = 2.0 * lambda / (pow(lambda,2.0) + 1.0 - c2);
    double delta  = pow(lambda,2.0) - c2;

    double d1  = 3.0 * pow(lambda, 2.0) / k * (k * (6.0 * pow(lambda, 2.0) - 1.0) - 2.0 * lambda);
    double d2  = 8.0 * pow(lambda, 2.0) / k * (k * (11.0 * pow(lambda, 2.0) - 1.0) - 2.0 * lambda);

    double a21 = 3.0 * c3 * (pow(k,2.0) - 2.0) / (4.0 * (1.0 + 2.0 * c2));
    double a22 = 3.0 * c3 / (4.0 * (1.0 + 2.0 * c2));
    double a23 = -3.0 * c3 * lambda / (4.0 * k * d1) * (3.0 * pow(k,3.0) * lambda - 6.0 * k * (k - lambda) + 4.0);
    double a24 = -3.0 * c3 * lambda / (4.0 * k * d1) * (2.0 + 3.0 * k * lambda);

    double b21 = -3.0 * c3 * lambda / (2.0 * d1) * (3.0 * k * lambda - 4.0);
    double b22 = 3.0 * c3 * lambda / d1;

    double d21 = -c3 / (2.0 * pow(lambda,2.0));

    double a31 = -9.0 * lambda / (4.0 * d2) * (4.0 * c3 * (k * a23 - b21) + k * c4 * (4.0 + pow(k, 2.0)))
            + (9.0 * pow(lambda, 2.0) + 1.0 - c2) / (2.0 * d2) * (3.0 * c3 * (2.0 * a23 - k * b21) + c4 * (2.0 + 3.0 * pow(k,2.0)));
    double a32 = -1.0 / d2 * (9.0 * lambda / 4.0 * (4.0 * c3 * (k * a24 - b22) + k * c4)
                              + 3.0 / 2.0 * (9.0 * pow(lambda, 2.0) + 1.0 - c2) * (c3 * (k * b22 + d21 - 2.0 * a24) - c4));

    double b31 = 3.0 / (8.0 * d2) * (8.0 * lambda * (3.0 * c3 * (k * b21 - 2.0 * a23) - c4 * (2.0 + 3.0 * pow(k, 2.0)))
                                     + (9.0 * pow(lambda, 2.0) + 1.0 + 2.0 * c2) * (4.0 * c3 * (k * a23 - b21) + k * c4 * (4.0 + pow(k, 2.0))));
    double b32 = 1.0 / d2 * (9.0 * lambda * (c3 * (k * b22 + d21 - 2 * a24) - c4)
                             + 3.0/8.0 * (9.0 * pow(lambda, 2.0) + 1.0 + 2.0 * c2) * (4.0 * c3 * (k * a24 - b22) + k * c4));

    double d31 = 3.0 / (64.0 * pow(lambda, 2.0)) * (4.0 * c3 * a24 + c4);
    double d32 = 3.0 / (64.0 * pow(lambda, 2.0)) * (4.0 * c3 * (a23 - d21) + c4 * (4.0 + pow(k,2.0)));

    double a1  = -3.0/2.0 * c3 * (2.0 * a21 + a23 + 5.0 * d21) - 3.0/8.0 * c4 * (12.0 - pow(k,2.0));
    double a2  = 3.0/2.0 * c3 * (a24 - 2.0 * a22) + 9.0/8.0 * c4;

    double s1  = 1.0 / (2.0 * lambda * (lambda * (1.0 + pow(k, 2.0)) - 2.0 * k)) * (
                3.0 / 2.0 * c3 * (2.0 * a21 * (pow(k, 2.0) - 2.0) - a23 * (pow(k, 2.0) + 2.0) - 2.0 * k * b21) - 3.0/8.0 * c4 * (
                    3.0 * pow(k, 4.0) - 8.0 * pow(k, 2.0) + 8.0));
    double s2  = 1.0 / (2.0 * lambda * (lambda * (1.0 + pow(k, 2.0)) - 2.0 * k)) * (
                3.0/2.0 * c3 * (2.0 * a22 * (pow(k, 2.0) - 2.0) + a24 * (pow(k, 2.0) + 2.0) + 2.0 * k * b22 + 5.0 * d21) + 3.0/8.0 * c4 * (
                    12.0 - pow(k, 2.0)));

    double l1  = a1 + 2.0 * pow(lambda, 2.0) * s1;
    double l2  = a2 + 2.0 * pow(lambda, 2.0) * s2;

    if (orbitType == horizontal_lyapunov_orbit )
    {
        Ax = amplitude;
        Az = 0.0;
    }
    if (orbitType == vertical_lyapunov_orbit )
    {
        Ax = 0.0;
        Az = amplitude;
    }
    if (orbitType == halo_orbit )
    {
        Az = amplitude;
        Ax = pow(((-delta - l2 * pow(Az, 2.0)) / l1), 0.5);
    }

    //std::cout << "Ax = " << Ax << ", Az =  = " << Az << std::endl;

    double omega1 = 0.0;
    double omega2 = s1 * pow(Ax, 2.0) + s2 * pow(Az, 2.0);
    double omega  = 1.0 + omega1 + omega2;

    double tau1   = 0.0;
    double deltan = 2.0 - n;

    double x             = a21 * pow(Ax, 2.0) + a22 * pow(Az, 2.0) - Ax * std::cos(tau1) + (a23 * pow(Ax, 2.0) - a24 * pow(Az, 2.0)) * std::cos(2.0 * tau1) + (a31 * pow(Ax, 3.0) - a32 * Ax * pow(Az, 2.0)) * std::cos(3.0 * tau1);
    double y             = k * Ax * std::sin(tau1) + (b21 * pow(Ax, 2.0) - b22 * pow(Az, 2.0)) * std::sin(2 * tau1) + (b31 * pow(Ax, 3.0) - b32 * Ax * pow(Az, 2.0)) * std::sin(3.0 * tau1);
    double z             = deltan * Az * std::cos(tau1) + deltan * d21 * Ax * Az * (std::cos(2.0 * tau1) - 3.0) + deltan * (d32 * Az * pow(Ax, 2.0) - d31 * pow(Az, 3.0)) * std::cos(3.0 * tau1);
    double xdot          = omega * lambda * Ax * std::sin(tau1) - 2.0 * omega * lambda * (a23 * pow(Ax, 2.0) - a24 * pow(Az, 2.0)) * std::sin(2.0 * tau1) - 3.0 * omega * lambda * (a31 * pow(Ax, 3.0) - a32 * Ax * pow(Az, 2.0)) * std::sin(3.0 * tau1);
    double ydot          = omega * lambda * (k * Ax * std::cos(tau1) + 2.0 * (b21 * pow(Ax, 2.0) - b22 * pow(Az, 2.0)) * std::cos(2.0 * tau1) + 3.0 * (b31 * pow(Ax, 3.0) - b32 * Ax * pow(Az, 2.0)) * std::cos(3.0 * tau1));
    double zdot          = -1.0 * omega * lambda * deltan * Az * std::sin(tau1) - 2.0 * omega * lambda * deltan * d21 * Ax * Az * std::sin(2.0 * tau1) - 3.0 * omega * lambda * deltan * (d32 * Az * pow(Ax, 2.0) - d31 * pow(Az, 3.0)) * std::sin(3.0 * tau1);
    double orbitalPeriod = 2.0 * tudat::mathematical_constants::PI / (lambda * omega);

    if (librationPointNr == 1)
    {
        initialStateVector(0) = (x - 1.0) * gammaL + 1.0 - massParameter;
    }
    if (librationPointNr == 2){
        initialStateVector(0) = (x + 1.0) * gammaL + 1.0 - massParameter;
    }

    initialStateVector(1) = y * gammaL;
    initialStateVector(2) = z * gammaL;
    initialStateVector(3) = xdot * gammaL;
    initialStateVector(4) = ydot * gammaL;
    initialStateVector(5) = zdot * gammaL;

    return std::make_pair( initialStateVector, orbitalPeriod );
}

double initializeEarthMoonPeriodicOrbitAmplitude(
        const int librationPointNumber, const CR3BPPeriodicOrbitTypes orbitType, const int guessIteration )
{
    if( librationPointNumber < 1 ||librationPointNumber > 2 )
    {
        throw std::runtime_error( "Error when getting Earth-Moon periodic orbit amplitude, libration point " +
                                  std::to_string( librationPointNumber ) + " is not supported." );
    }
    double amplitude = TUDAT_NAN;
    if( guessIteration == 0 )
    {
        if (orbitType == horizontal_lyapunov_orbit )
        {
            if (librationPointNumber == 1 )
            {
                amplitude = 1.0e-3;
            }
            else if (librationPointNumber == 2 )
            {
                amplitude = 1.0e-4;
            }
        }
        else if (orbitType == vertical_lyapunov_orbit )
        {
            if (librationPointNumber == 1)
            {
                amplitude = 1.0e-1;
            }
            else if (librationPointNumber == 2)
            {
                amplitude = 1.0e-1;
            }
        }
        else if (orbitType == halo_orbit )
        {
            if (librationPointNumber == 1)
            {
                amplitude = -1.1e-1;
            }
            else if (librationPointNumber == 2)
            {
                amplitude = 1.5e-1;
            }
        }
    }
    else if( guessIteration == 1 )
    {

        if (orbitType == horizontal_lyapunov_orbit )
        {
            if (librationPointNumber == 1)
            {
                amplitude = 1.0e-4;
            }
            else if (librationPointNumber == 2)
            {
                amplitude = 1.0e-3;
            }
        }
        else if (orbitType == vertical_lyapunov_orbit )
        {
            if (librationPointNumber == 1)
            {
                amplitude = 2.0e-1;
            }
            else if (librationPointNumber == 2)
            {
                amplitude = 2.0e-1;
            }
        }
        else if (orbitType == halo_orbit )
        {
            if (librationPointNumber == 1)
            {
                amplitude = -1.2e-1;
            }
            else if (librationPointNumber == 2)
            {
                amplitude = 1.6e-1;
            }
        }
    }
    else
    {
        throw std::runtime_error( "Error when getting Earth-Moon periodic orbit amplitude, iteration " +
                                  std::to_string( guessIteration ) + " is not supported." );
    }

    return amplitude;
}



bool continueCR3BPDifferentialCorrection(
        const Eigen::Vector6d stateVector, const int numberOfIterations,
        const CR3BPPeriodicOrbitGenerationSettings& periodicOrbitSettings )
{
    double maximumPositionDeviationToUse, maximumVelocityDeviationToUse;

    // If the maximum number of iterations has been reached, return a zero vector to stop the numerical continuation
    if ( numberOfIterations > periodicOrbitSettings.maximumDifferentialCorrections_ )
    {
        // Relax the periodicity constraints after exceeding the maximum number of iterations instead of termination
        maximumPositionDeviationToUse = 10.0 * periodicOrbitSettings.maximumPositionDeviation_;
        maximumVelocityDeviationToUse = 10.0 * periodicOrbitSettings.maximumVelocityDeviation_;
    }
    // Relax the maximum deviation requirements to compute the horizontal Lyapunov family in L2
    else if ( numberOfIterations > 10 &&
              periodicOrbitSettings.orbitType_ == horizontal_lyapunov_orbit &&
              periodicOrbitSettings.librationPointNumber_ == 2)
    {
        maximumPositionDeviationToUse = 10.0 * periodicOrbitSettings.maximumPositionDeviation_;
        maximumVelocityDeviationToUse = 10.0 * periodicOrbitSettings.maximumVelocityDeviation_;
    }
    else
    {
        maximumPositionDeviationToUse = periodicOrbitSettings.maximumPositionDeviation_;
        maximumVelocityDeviationToUse = periodicOrbitSettings.maximumVelocityDeviation_;
    }

    double positionDeviationFromPeriodicOrbit, velocityDeviationFromPeriodicOrbit;
    if ( periodicOrbitSettings.orbitType_ == axial_orbit )
    {
        // Initial condition for axial family should be [x, 0, 0, 0, ydot, zdot]
        positionDeviationFromPeriodicOrbit = sqrt(pow(stateVector(1), 2) + pow(stateVector(2), 2));
        velocityDeviationFromPeriodicOrbit = sqrt(pow(stateVector(3), 2));
    }
    else
    {
        // Initial condition for other families should be [x, 0, y, 0, ydot, 0]
        positionDeviationFromPeriodicOrbit = sqrt(pow(stateVector(1), 2));
        velocityDeviationFromPeriodicOrbit = sqrt(pow(stateVector(3), 2) + pow(stateVector(5), 2));
    }
    return( ( positionDeviationFromPeriodicOrbit > maximumPositionDeviationToUse ) ||
            ( velocityDeviationFromPeriodicOrbit > maximumVelocityDeviationToUse ) );
}


Eigen::Vector7d computeCR3BPPeriodicOrbitsDifferentialCorrection(
        const Eigen::Matrix6d& stateTransitionMatrix,
        const Eigen::Vector6d& cartesianState,
        const double currentTime,
        const int currentIteration,
        const CR3BPPeriodicOrbitGenerationSettings& periodicOrbitSettings )
{
    Eigen::Vector7d differentialCorrection;
    Eigen::Vector3d corrections;
    Eigen::Matrix3d updateMatrix;
    Eigen::Vector3d multiplicationMatrix;

    // Compute the accelerations and velocities on the spacecraft
    Eigen::Vector6d cartesianAccelerations =
            computeCr3bpStateDerivative( currentTime, cartesianState, periodicOrbitSettings.massParameter_ );

    bool xPositionFixed = false;
    if ( currentIteration > 10 &&
         periodicOrbitSettings.orbitType_ == axial_orbit &&
         periodicOrbitSettings.librationPointNumber_ == 2 )
    {
        xPositionFixed = true;
    }

    // If type is axial, the desired state vector has the form [x, 0, 0, 0, ydot, zdot] and
    // requires a differential correction for {x, ydot, T/2}
    if ( periodicOrbitSettings.orbitType_ == axial_orbit )
    {
        // Check which deviation is larger: x-velocity or z-position.
        if ( std::abs(cartesianState(2)) < std::abs(cartesianState(3)) && !xPositionFixed )
        {
            // Correction on {x, ydot, T/2} for constant {zdot}

            // Set the correct multiplication matrix (state at T/2)
            multiplicationMatrix << cartesianState(1), cartesianState(2), cartesianState(3);

            // Compute the update matrix.
            updateMatrix << stateTransitionMatrix(1, 0), stateTransitionMatrix(1, 4), cartesianState(4),
                    stateTransitionMatrix(2, 0), stateTransitionMatrix(2, 4), cartesianState(5),
                    stateTransitionMatrix(3, 0), stateTransitionMatrix(3, 4), cartesianAccelerations(3);

            // Compute the necessary differential correction.
            corrections = updateMatrix.inverse() * multiplicationMatrix;

            // Put corrections in correct format.
            differentialCorrection.setZero( );
            differentialCorrection(0) = -corrections(0);
            differentialCorrection(4) = -corrections(1);
            differentialCorrection(6) = -corrections(2);
        }
        else
        {
            // Correction on {ydot, zdot T/2} for constant {x}

            // Set the correct multiplication matrix (state at T/2)
            multiplicationMatrix << cartesianState(1), cartesianState(2), cartesianState(3);

            // Compute the update matrix.
            updateMatrix << stateTransitionMatrix(1, 4), stateTransitionMatrix(1, 5), cartesianState(4),
                    stateTransitionMatrix(2, 4), stateTransitionMatrix(2, 5), cartesianState(5),
                    stateTransitionMatrix(3, 4), stateTransitionMatrix(3, 5), cartesianAccelerations(3);

            // Compute the necessary differential correction.
            corrections = updateMatrix.inverse() * multiplicationMatrix;

            // Put corrections in correct format.
            differentialCorrection.setZero( );
            differentialCorrection(4) = -corrections(0);
            differentialCorrection(5) = -corrections(1);
            differentialCorrection(6) = -corrections(2);
        }
    }

    // If type is not axial, the desired state vector has the form [x, 0, z, 0, ydot, 0]
    // and requires a differential correction for either {z, ydot, T/2} or {x, ydot, T/2}
    else
    {

        // Check which deviation is larger: x-velocity or z-velocity.
        if ( std::abs(cartesianState(3)) < std::abs(cartesianState(5)) ||
             periodicOrbitSettings.orbitType_ == horizontal_lyapunov_orbit )
        {
            // Correction on {z, ydot, T/2} for constant {x}

            // Set the correct multiplication matrix (state at T/2)
            multiplicationMatrix << cartesianState(1), cartesianState(3), cartesianState(5);

            // Compute the update matrix.
            updateMatrix         << stateTransitionMatrix(1, 2), stateTransitionMatrix(1, 4), cartesianState(4),
                    stateTransitionMatrix(3, 2), stateTransitionMatrix(3, 4), cartesianAccelerations(3),
                    stateTransitionMatrix(5, 2), stateTransitionMatrix(5, 4), cartesianAccelerations(5);

            // Compute the necessary differential correction.
            corrections = updateMatrix.inverse() * multiplicationMatrix;

            // Put corrections in correct format.
            differentialCorrection.setZero();
            differentialCorrection(2) = -corrections(0);
            differentialCorrection(4) = -corrections(1);
            differentialCorrection(6) = -corrections(2);
        }
        else
        {
            // Correction on {x, ydot, T/2} for constant {z}
            // Set the correct multiplication matrix (state at T/2)
            multiplicationMatrix << cartesianState(1), cartesianState(3), cartesianState(5);

            // Compute the update matrix.
            updateMatrix << stateTransitionMatrix(1, 0), stateTransitionMatrix(1, 4), cartesianState(4),
                    stateTransitionMatrix(3, 0), stateTransitionMatrix(3, 4), cartesianAccelerations(3),
                    stateTransitionMatrix(5, 0), stateTransitionMatrix(5, 4), cartesianAccelerations(5);

            // Compute the necessary differential correction.
            corrections = updateMatrix.inverse() * multiplicationMatrix;

            // Put corrections in correct format.
            differentialCorrection.setZero();
            differentialCorrection(0) = -corrections(0);
            differentialCorrection(4) = -corrections(1);
            differentialCorrection(6) = -corrections(2);
        }
    }

    // Return differential correction.
    return differentialCorrection;
}



CR3BPPeriodicOrbitConditions createCR3BPPeriodicOrbit(
        const Eigen::Vector6d& initialStateVector,
        double orbitalPeriod,
        const CR3BPPeriodicOrbitGenerationSettings& periodicOrbitSettings,
        const std::shared_ptr< tudat::numerical_integrators::IntegratorSettings< double > > integratorSettings )
{
    // Propagate CR3BP with state transition matrix
    std::map< double, Eigen::MatrixXd > fullStateHistoryUpToHalfPeriod =
            performCR3BPWithStmIntegration(
                integratorSettings, periodicOrbitSettings.massParameter_, initialStateVector, orbitalPeriod / 2.0, true, false );

    // Retrieve numerical results
    double currentTime = fullStateHistoryUpToHalfPeriod.begin( )->first;
    Eigen::Vector6d finalStateVector = fullStateHistoryUpToHalfPeriod.begin( )->second.block( 0, 0, 6, 1 );
    Eigen::Matrix6d finalStateStateTransitionMatix = fullStateHistoryUpToHalfPeriod.begin( )->second.block( 0, 1, 6, 6 );

    // Apply differential correction and propagate to half-period point until converged.
    int numberOfIterations = 0;
    Eigen::Vector7d differentialCorrection = Eigen::Vector7d::Zero( );
    Eigen::Vector6d correctedInitialState = initialStateVector;
    do
    {
        // Perform differential correction
        differentialCorrection = computeCR3BPPeriodicOrbitsDifferentialCorrection(
                    finalStateStateTransitionMatix, finalStateVector, currentTime, numberOfIterations, periodicOrbitSettings );
        correctedInitialState += differentialCorrection.segment( 0, 6 );
        orbitalPeriod  = orbitalPeriod + 2.0 * differentialCorrection( 6 );

        // Repropagate with corrected initial state
        fullStateHistoryUpToHalfPeriod = performCR3BPWithStmIntegration(
                    integratorSettings, periodicOrbitSettings.massParameter_,
                    correctedInitialState, orbitalPeriod / 2.0, true, false );

        // Retrieve numerical results
        currentTime = fullStateHistoryUpToHalfPeriod.begin( )->first;
        finalStateVector = fullStateHistoryUpToHalfPeriod.begin( )->second.block( 0, 0, 6, 1 );
        finalStateStateTransitionMatix = fullStateHistoryUpToHalfPeriod.begin( )->second.block( 0, 1, 6, 6 );

        numberOfIterations += 1;
    }
    while ( continueCR3BPDifferentialCorrection( finalStateVector, numberOfIterations, periodicOrbitSettings ) );

    std::map< double, Eigen::MatrixXd > fullOrbitNumericalResults = performCR3BPWithStmIntegration(
                integratorSettings, periodicOrbitSettings.massParameter_,
                correctedInitialState, orbitalPeriod, true, false );
    Eigen::Matrix6d monodromyMatrix = fullOrbitNumericalResults.rbegin( )->second.block( 0, 1, 6, 6 );

    return CR3BPPeriodicOrbitConditions(
                correctedInitialState, finalStateVector, monodromyMatrix, orbitalPeriod, numberOfIterations,
                periodicOrbitSettings );
}

bool checkForMonodromyUnitEigenvalues(
        const Eigen::Matrix6d &monodromyMatrix,
        const CR3BPPeriodicOrbitGenerationSettings periodicOrbitSettings )
{
    bool moduleOneInsteadOfRealOne;

    // Exception for the horizontal Lyapunov family in Earth-Moon L2: eigenvalue may be of module
    // one instead of a real one to compute a more extensive family
    if ( ( periodicOrbitSettings.librationPointNumber_ == 2 ) && ( periodicOrbitSettings.orbitType_ == horizontal_lyapunov_orbit ) )
    {
        moduleOneInsteadOfRealOne = true;
    }
    else
    {
        moduleOneInsteadOfRealOne = false;
    }

    // Initialize variables
    bool eigenvalueRealOne = false;

    // Reshape the STM for one period to matrix form and compute the eigenvalues
    Eigen::EigenSolver<Eigen::Matrix6d> eig( monodromyMatrix );

    // Determine whether the monodromy matrix contains at least one eigenvalue of real
    // one within the maximumAllowedEigenvalueDeviation_
    for (int i = 0; i <= 5; i++)
    {
        if (std::abs(eig.eigenvalues().imag()(i)) < periodicOrbitSettings.maximumEigenvalueDeviation_ )
        {
            if (std::abs(eig.eigenvalues().real()(i) - 1.0) < periodicOrbitSettings.maximumEigenvalueDeviation_ )
            {
                eigenvalueRealOne = true;
            }
        }
    }

    // Optional argument to generalize the test from real one to module one
    if (moduleOneInsteadOfRealOne == true)
    {
        for (int i = 0; i <= 5; i++)
        {
            if (std::abs(std::abs(eig.eigenvalues()(i)) - 1.0 ) < periodicOrbitSettings.maximumEigenvalueDeviation_ )
            {
                eigenvalueRealOne = true;
            }
        }
    }

    return eigenvalueRealOne;
}


bool continueCR3BPNumericalContinuation(
        const Eigen::Matrix6d& monodromyMatrix,
        const int numberOfOrbitsGenerated,
        const CR3BPPeriodicOrbitGenerationSettings periodicOrbitSettings )
{
    // Check termination conditions
    bool continueNumericalContinuation = true;
    if( monodromyMatrix != monodromyMatrix )
    {
        continueNumericalContinuation = false;
        std::cout << "\n\nNUMERICAL CONTINUATION STOPPED DUE TO NaN STATE TRANSITION MATROX \n\n" << std::endl;
    }
    else if ( numberOfOrbitsGenerated >= periodicOrbitSettings.maximumNumericalContinuations_ )
    {
        std::cout<<numberOfOrbitsGenerated<<" "<<periodicOrbitSettings.maximumNumericalContinuations_<<std::endl;
        continueNumericalContinuation = false;
        std::cout << "\n\nNUMERICAL CONTINUATION STOPPED DUE TO EXCEEDING MAXIMUM NUMBER OF ITERATIONS\n\n" << std::endl;
    }
    else
    {        // Check eigenvalue condition (at least one pair equalling a real one)
        continueNumericalContinuation = checkForMonodromyUnitEigenvalues( monodromyMatrix, periodicOrbitSettings );
        if( !continueNumericalContinuation )
        {
            std::cout << "\n\nNUMERICAL CONTINUATION STOPPED DUE TO EIGENVALUES OF MONODROMY MATRIX\n\n" << std::endl;
        }
    }
    return continueNumericalContinuation;
}


double getDefaultPseudoArcLength(
        const double distanceIncrement,
        const Eigen::Vector6d& currentState )
{
    return std::max( distanceIncrement / currentState.segment( 0, 3 ).norm( ), 0.01 );
}

void createCR3BPPeriodicOrbitsThroughNumericalContinuation(
        std::vector< CR3BPPeriodicOrbitConditions >& periodicOrbits,
        const std::shared_ptr< tudat::numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const CR3BPPeriodicOrbitGenerationSettings periodicOrbitSettings,
        const std::function< double( const Eigen::Vector6d& ) > pseudoArcLengthFunction )

{
    if( periodicOrbits.size( ) < 2 )
    {
        throw std::runtime_error( "Error when creating periodic orbits through numerical continuation, requires at least two orbits as input, but only " +
                                  std::to_string( periodicOrbits.size( ) ) + " provided." );
    }
    // Set exit parameters of continuation procedure
    int numberOfInitialConditions = periodicOrbits.size( );
    int maximumNumberOfInitialConditions = 10000;

    // Initialize state vectors and orbital periods
    Eigen::Vector6d initialStateVector = Eigen::Vector6d::Zero( );

    // Generate periodic orbits until termination
    double orbitalPeriod  = 0.0, periodIncrement = 0.0, pseudoArcLengthCorrection = 0.0;
    bool continueContinuation = true;
    Eigen::Vector6d stateIncrement;

    while( ( numberOfInitialConditions < maximumNumberOfInitialConditions ) && continueContinuation )
    {
        // Determine increments to state and time
        stateIncrement = periodicOrbits[ numberOfInitialConditions - 1 ].initialState_ -
                periodicOrbits[ numberOfInitialConditions - 2 ].initialState_;
        periodIncrement = periodicOrbits[ numberOfInitialConditions - 1 ].orbitPeriod_ -
                periodicOrbits[ numberOfInitialConditions - 2 ].orbitPeriod_;
        pseudoArcLengthCorrection =
                pseudoArcLengthFunction( stateIncrement );

        // Apply numerical continuation
        initialStateVector = periodicOrbits[ numberOfInitialConditions - 1 ].initialState_ +
                stateIncrement * pseudoArcLengthCorrection;
        orbitalPeriod = periodicOrbits[ numberOfInitialConditions - 1 ].orbitPeriod_ +
                periodIncrement * pseudoArcLengthCorrection;

//        std::cout<<numberOfInitialConditions<<" "<<periodicOrbits[ numberOfInitialConditions - 1 ].orbitPeriod_<<std::endl;
        try
        {
            periodicOrbits.push_back( createCR3BPPeriodicOrbit(
                                          initialStateVector, orbitalPeriod, periodicOrbitSettings,
                                          integratorSettings ) );
            numberOfInitialConditions += 1;

            continueContinuation =
                    continueCR3BPNumericalContinuation(
                        periodicOrbits[ numberOfInitialConditions - 1 ].monodromyMatrix_,
                    numberOfInitialConditions,
                    periodicOrbitSettings );
        }
        catch( ... )
        {
            std::cout << "\n\nNUMERICAL CONTINUATION STOPPED DUE TO EXCEPTION IN INTEGRATION\n\n" << std::endl;
            continueContinuation = false;

        }

    }

    //    writeFinalResultsToFiles( periodicOrbitModel->getLibrationPointNumber( ), periodicOrbitModel->getOrbitType( ), initialConditions, differentialCorrections );
}

//#endif

}

}
