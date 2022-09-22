/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CR3BPPERIODICORBITS_H
#define TUDAT_CR3BPPERIODICORBITS_H

#include <string>
#include <boost/bind/bind.hpp>

#include "tudat/simulation/propagation_setup/createStateDerivativeModel.h"

namespace tudat
{

namespace propagators
{


enum CR3BPPeriodicOrbitTypes
{
    horizontal_lyapunov_orbit = 0,
    vertical_lyapunov_orbit = 1,
    halo_orbit = 2,
    axial_orbit = 3
};


std::string getPeriodicOrbitName( enum CR3BPPeriodicOrbitTypes orbitType );

struct CR3BPPeriodicOrbitGenerationSettings
{
    CR3BPPeriodicOrbitGenerationSettings(
            const double massParameter,
            const CR3BPPeriodicOrbitTypes orbitType,
            const int librationPointNumber,
            const int maximumDifferentialCorrections,
            const double maximumPositionDeviation,
            const double maximumVelocityDeviation ,
            const int maximumNumericalContinuations,
            const double maximumEigenvalueDeviation ):
        massParameter_( massParameter ), orbitType_( orbitType ), librationPointNumber_( librationPointNumber ),
        maximumDifferentialCorrections_( maximumDifferentialCorrections ),
        maximumPositionDeviation_( maximumPositionDeviation ),
        maximumVelocityDeviation_( maximumVelocityDeviation ),
        maximumNumericalContinuations_( maximumNumericalContinuations ),
        maximumEigenvalueDeviation_( maximumEigenvalueDeviation ){ }

    double massParameter_;
    CR3BPPeriodicOrbitTypes orbitType_;
    int librationPointNumber_;

    int maximumDifferentialCorrections_;
    double maximumPositionDeviation_;
    double maximumVelocityDeviation_;

    int maximumNumericalContinuations_;
    double maximumEigenvalueDeviation_;

};


class CR3BPPeriodicOrbitConditions
{
public:
    CR3BPPeriodicOrbitConditions(
            const Eigen::Vector6d& initialState,
            const double orbitPeriod,
            const double massParameter ):
        initialState_( initialState ),
        orbitPeriod_( orbitPeriod ),
        massParameter_( massParameter ){ }

    Eigen::Vector6d getInitialState( )
    {
        return initialState_;
    }

    double getOrbitPeriod( )
    {
        return orbitPeriod_;
    }

    double getMassParameter( )
    {
        return massParameter_;
    }

protected:

    const Eigen::Vector6d initialState_;
    const double orbitPeriod_;
    const double massParameter_;
};

class PropagatedCR3BPPeriodicOrbitConditions: public CR3BPPeriodicOrbitConditions
{
public:
    PropagatedCR3BPPeriodicOrbitConditions(
            const Eigen::Vector6d& initialState,
            const double orbitPeriod,
            const double massParameter,
            const Eigen::Matrix6d monodromyMatrix,
            const CR3BPPeriodicOrbitTypes orbitType,
            const int librationPoint ):
        CR3BPPeriodicOrbitConditions( initialState, orbitPeriod, massParameter ),
        monodromyMatrix_( monodromyMatrix ),
        orbitType_( orbitType ),
        librationPoint_( librationPoint ){ }

    Eigen::Matrix6d getMonodromyMatrix( )
    {
        return monodromyMatrix_;
    }

    CR3BPPeriodicOrbitTypes getOrbitType( )
    {
        return orbitType_;
    }

    int getLibrationPointNumber( )
    {
        return librationPoint_;
    }

protected:

    Eigen::Matrix6d monodromyMatrix_;
    CR3BPPeriodicOrbitTypes orbitType_;
    int librationPoint_;

};

class GeneratedCR3BPPeriodicOrbitConditions: public PropagatedCR3BPPeriodicOrbitConditions
{
public:
    GeneratedCR3BPPeriodicOrbitConditions(
            const Eigen::Vector6d initialState,
            const Eigen::Vector6d halfPeriodStateVector,
            const Eigen::Matrix6d monodromyMatrix,
            const double orbitPeriod,
            const int numberOfDifferentialCorrections,
            const CR3BPPeriodicOrbitGenerationSettings& generationSettings ):
        PropagatedCR3BPPeriodicOrbitConditions(
            initialState, orbitPeriod, generationSettings.massParameter_, monodromyMatrix,
            generationSettings.orbitType_, generationSettings.librationPointNumber_ ),
        halfPeriodStateVector_( halfPeriodStateVector),
        numberOfDifferentialCorrections_( numberOfDifferentialCorrections ),
        generationSettings_( generationSettings ){ }

    Eigen::Vector6d getHalfPeriodStateVector( )
    {
        return halfPeriodStateVector_;
    }

    int getNumberOfDifferentialCorrections( )
    {
        return numberOfDifferentialCorrections_;
    }

    CR3BPPeriodicOrbitGenerationSettings getGenerationSettings( )
    {
        return generationSettings_;
    }

private:

    Eigen::Vector6d halfPeriodStateVector_;

    int numberOfDifferentialCorrections_;

    CR3BPPeriodicOrbitGenerationSettings generationSettings_;
};

std::map< double, std::shared_ptr< CR3BPPeriodicOrbitConditions > > sortPeriodicOrbitsByInitialJacobiEnergy(
        const std::vector< std::shared_ptr< CR3BPPeriodicOrbitConditions > >& periodicOrbits,
        const double massParameter );

std::map< double, double > getOrbitJacobiEnergyHistory(
        const std::map< double, Eigen::Vector6d >& orbitStateHistory,
        const double massParameter );

std::map< double, Eigen::Vector6d > propagatePeriodicOrbit(
        const std::shared_ptr< CR3BPPeriodicOrbitConditions > orbitDefinition,
        const std::shared_ptr< tudat::numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const double numberOfPeriods = 1.0 );


std::pair< Eigen::Vector6d, double >  richardsonApproximationLibrationPointPeriodicOrbit(
        const double massParameter,
        const CR3BPPeriodicOrbitTypes orbitType,
        int librationPointNr, double amplitude, double n = 1.0 );

std::pair< Eigen::Vector6d, double >  richardsonApproximationEarthMoonLibrationPointPeriodicOrbit(
        const double massParameter,
        const CR3BPPeriodicOrbitTypes orbitType,
        int librationPointNr, const int guessIteration, double n = 1.0 );

double initializeEarthMoonPeriodicOrbitAmplitude(
        const int librationPointNumber, const CR3BPPeriodicOrbitTypes orbitType, const int guessIteration );


bool continueCR3BPDifferentialCorrection(
        const Eigen::Vector6d stateVector, const int numberOfIterations,
        const CR3BPPeriodicOrbitGenerationSettings& periodicOrbitSettings );


Eigen::Vector7d computeCR3BPPeriodicOrbitsDifferentialCorrection(
        const Eigen::Matrix6d& stateTransitionMatrix,
        const Eigen::Vector6d& cartesianState,
        const double currentTime,
        const int currentIteration,
        const CR3BPPeriodicOrbitGenerationSettings& periodicOrbitSettings );

std::shared_ptr< GeneratedCR3BPPeriodicOrbitConditions > createCR3BPPeriodicOrbit(
        const Eigen::Vector6d& initialStateVector,
        double orbitalPeriod,
        const CR3BPPeriodicOrbitGenerationSettings& periodicOrbitSettings,
        const std::shared_ptr< tudat::numerical_integrators::IntegratorSettings< double > > integratorSettings );

bool checkForMonodromyUnitEigenvalues( const Eigen::Matrix6d &monodromyMatrix );

bool continueCR3BPNumericalContinuation(
        const Eigen::Matrix6d& monodromyMatrix,
        const Eigen::Vector6d& cartesianState,
        const int numberOfOrbitsGenerated  );

double getDefaultPseudoArcLength(
        const double distanceIncrement,
        const Eigen::Vector6d& currentState );

//void createCR3BPPeriodicOrbitsThroughNumericalContinuation(
//        std::vector< std::shared_ptr< CR3BPPeriodicOrbitConditions > >& periodicOrbits,
//        const std::shared_ptr< tudat::numerical_integrators::IntegratorSettings< double > > integratorSettings,
//        const CR3BPPeriodicOrbitGenerationSettings periodicOrbitSettings,
//        const std::function< double( const Eigen::Vector6d& ) > pseudoArcLengthFunction =
//        std::bind( &getDefaultPseudoArcLength, 1.0E-4, std::placeholders::_1 ) );

void createCR3BPPeriodicOrbitsThroughNumericalContinuation(
        std::vector< std::shared_ptr< PropagatedCR3BPPeriodicOrbitConditions > >& periodicOrbits,
        const std::shared_ptr< tudat::numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const CR3BPPeriodicOrbitGenerationSettings periodicOrbitSettings,
        const std::function< double( const Eigen::Vector6d& ) > pseudoArcLengthFunction =
        std::bind( &getDefaultPseudoArcLength, 1.0E-4, std::placeholders::_1 ));

} // namespace propagators

} // namespace tudat

#endif // TUDAT_CR3BPPERIODICORBITS_H
