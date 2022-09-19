/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/gravitation/jacobiEnergy.h"
#include "tudat/math/interpolators/lookupScheme.h"
#include "tudat/simulation/propagation_setup/createCR3BPManifolds.h"

namespace tudat
{

namespace propagators
{

double determineEigenvectorSign( const Eigen::Vector6d& eigenvector )
{
    double eigenvectorSign;

    if ( eigenvector( 0 ) > 0.0 )
    {
        eigenvectorSign = 1.0;
    }
    else
    {
        eigenvectorSign = -1.0;
    }

    return eigenvectorSign;
}


void determineStableUnstableEigenvectors(
        Eigen::Vector6d& stableEigenvector,
        Eigen::Vector6d& unstableEigenvector,
        const Eigen::MatrixXd& monodromyMatrix,
        double maxEigenvalueDeviation )
{
    int indexMaximumEigenvalue;
    int indexMinimumEigenvalue;
    double maximumEigenvalue = 0.0;
    double minimumEigenvalue = 1000.0;

    // Compute eigenvectors of the monodromy matrix
    Eigen::EigenSolver< Eigen::MatrixXd > eig(monodromyMatrix);

    // Find the minimum (maximum) eigenvalue, corresponding to the stable (unstable) subspace
    auto monodromyEigenvalues =  eig.eigenvalues( );
    for ( int i = 0; i <= 5; i++ )
    {
        if ( monodromyEigenvalues.real( )( i ) > maximumEigenvalue &&
             std::abs( monodromyEigenvalues.imag( )( i ) ) < maxEigenvalueDeviation )
        {
            maximumEigenvalue = monodromyEigenvalues.real( )( i );
            indexMaximumEigenvalue = i;
        }
        if ( std::abs( monodromyEigenvalues.real( )( i ) ) < minimumEigenvalue &&
             std::abs( monodromyEigenvalues.imag( )( i )) < maxEigenvalueDeviation )
        {
            minimumEigenvalue = std::abs( monodromyEigenvalues.real( )( i ) );
            indexMinimumEigenvalue = i;
        }
    }

    unstableEigenvector = eig.eigenvectors( ).real( ).col( indexMaximumEigenvalue );
    stableEigenvector = eig.eigenvectors( ).real( ).col( indexMinimumEigenvalue );

    // Check whether the two selected eigenvalues belong to the same reciprocal pair
    if ( ( 1.0 / minimumEigenvalue - maximumEigenvalue) > maxEigenvalueDeviation )
    {
        std::cout << "\n\n\nERROR - EIGENVALUES MIGHT NOT BELONG TO SAME RECIPROCAL PAIR" << std::endl;
        throw std::exception( );
    }
}

void computeManifoldSetFromSinglePoint(
        std::vector< std::map< double, Eigen::Vector6d > >& manifoldStateHistories,
        const Eigen::MatrixXd& stateIncludingStm,
        const CR3BPPeriodicOrbitConditions periodicOrbitConditions,
        const CR3BPManifoldSettings manifoldSettings,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings )
{
    double initialTime = integratorSettings->initialTime_;
    double propagationTime = 2.0 * periodicOrbitConditions.orbitPeriod_;
    double currentInitialTimeStep = integratorSettings->initialTimeStep_;

    Eigen::Matrix6d localStateTransitionMatrix = stateIncludingStm.block( 0, 1, 6, 6 );
    Eigen::Vector6d localState = stateIncludingStm.block( 0, 0, 6, 1 );

    for( int i = 0; i < 4; i++ )
    {
        std::cout<<"Manifold index "<<std::endl;
        // Apply displacement epsilon from the periodic orbit at <numberOfTrajectoriesPerManifold> locations on the final orbit.
        Eigen::Vector6d monodromyMatrixEigenvector = manifoldSettings.getEigenVectorForInitialDisplacement( i );
        Eigen::Vector6d localNormalizedEigenvector = ( localStateTransitionMatrix * monodromyMatrixEigenvector ).normalized( );
        Eigen::Vector6d manifoldInitialState = localState + manifoldSettings.getEigenVectorDisplacementScaling( i ) *
                localNormalizedEigenvector;
        std::cout<<"Initial state deviation: "<<( manifoldSettings.getEigenVectorDisplacementScaling( i ) *
                   localNormalizedEigenvector ).transpose( )<<std::endl;
        std::cout<<"Norm. eigenvector: "<<localNormalizedEigenvector.transpose( )<<std::endl;
        double direction = manifoldSettings.getIntegrationDirection( i );
        integratorSettings->initialTime_ = initialTime;
        integratorSettings->initialTimeStep_ = direction * currentInitialTimeStep;
        manifoldStateHistories.push_back(
                    performCR3BPIntegration(
                        integratorSettings, periodicOrbitConditions.massParameter( ),
                        manifoldInitialState, initialTime + direction * propagationTime, false, true ) );
    }
}

std::vector< double > createManifoldDeparturePoints(
        const double orbitalPeriod,
        const int numberOfDeparturePoints,
        const std::map< double, Eigen::MatrixXd >& stateTransitionMatrixHistory )
{
    double idealDeparturePointSpacing = orbitalPeriod /
            ( static_cast< double >( numberOfDeparturePoints ) - 1.0 );

    std::vector< double > departurePoints;
    interpolators::HuntingAlgorithmLookupScheme< double > lookupScheme(
                utilities::createVectorFromMapKeys( stateTransitionMatrixHistory ) );

    double currentIdealDeparturePoint = 0.0;
    bool ignoreCurrentDataPoint = false;
    for( int i = 0; i < numberOfDeparturePoints; i++ )
    {
        double nearestTime = lookupScheme.getNearestIndependentVariableValue( currentIdealDeparturePoint );

        ignoreCurrentDataPoint = false;
        if( departurePoints.size( ) > 0 )
        {
            if( nearestTime == departurePoints.at( departurePoints.size( ) - 1 ) )
            {
                ignoreCurrentDataPoint = true;
            }
        }
        if( !ignoreCurrentDataPoint )
        {
            departurePoints.push_back( nearestTime );
        }
        currentIdealDeparturePoint += idealDeparturePointSpacing;
    }
    return departurePoints;
}
void computeManifolds( std::vector< std::vector< std::map< double, Eigen::Vector6d > > >& fullManifoldStateHistories,
                       const CR3BPPeriodicOrbitConditions periodicOrbitConditions,
                       const double eigenvectorDisplacementFromOrbit,
                       const int numberOfDeparturePoints,
                       const double maxEigenvalueDeviation,
                       const std::shared_ptr< tudat::numerical_integrators::IntegratorSettings< double > > integratorSettings )
{
    double jacobiEnergyOnOrbit = tudat::gravitation::computeJacobiEnergy(
                periodicOrbitConditions.massParameter( ), periodicOrbitConditions.initialState_ );
    std::cout << "\nInitial state vector:" << std::endl << periodicOrbitConditions.initialState_.transpose( ) <<std::endl
              << "\nwith C: " << jacobiEnergyOnOrbit    << ", T: " << periodicOrbitConditions.massParameter( ) << std::endl;;


    // Determine the eigenvector directions of the (un)stable subspace of the monodromy matrix
    Eigen::MatrixXd monodromyMatrix = periodicOrbitConditions.monodromyMatrix_;
    Eigen::Vector6d stableEigenvector;
    Eigen::Vector6d unstableEigenvector;
    determineStableUnstableEigenvectors( stableEigenvector, unstableEigenvector, monodromyMatrix, maxEigenvalueDeviation );

    // Propagate the periodicOrbitConditions.initialState_ for a full period and write output to file.
    std::map< double, Eigen::MatrixXd > stateTransitionMatrixHistory = performCR3BPWithStmIntegration(
                integratorSettings,
                periodicOrbitConditions.massParameter( ),
                periodicOrbitConditions.initialState_,
                periodicOrbitConditions.orbitPeriod_,
                true, true );

    std::vector< double > departurePoints = createManifoldDeparturePoints(
                periodicOrbitConditions.orbitPeriod_,
                numberOfDeparturePoints,
                stateTransitionMatrixHistory );

    CR3BPManifoldSettings manifoldSettings( stableEigenvector, stableEigenvector, eigenvectorDisplacementFromOrbit );
    for( unsigned int i = 0; i < departurePoints.size( ); i++ )
    {
        std::cout<<"Departure point "<<i<<" "<<departurePoints.at( i )<<std::endl;
        std::vector< std::map< double, Eigen::Vector6d > > currentManifoldStateHistories;
        computeManifoldSetFromSinglePoint( currentManifoldStateHistories,
                                           stateTransitionMatrixHistory.at( departurePoints.at( i ) ),
                                           periodicOrbitConditions, manifoldSettings, integratorSettings );
        fullManifoldStateHistories.push_back( currentManifoldStateHistories );
    }

}

}

}
