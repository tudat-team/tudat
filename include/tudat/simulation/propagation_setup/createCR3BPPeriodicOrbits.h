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
    horizontal_lyapunov_orbit,
    vertical_lyapunov_orbit,
    halo_orbit,
    axial_orbit
};

struct CR3BPPeriodicOrbitGenerationSettings
{
    CR3BPPeriodicOrbitGenerationSettings(
            const double massParameter,
            const CR3BPPeriodicOrbitTypes orbitType,
            const int librationPointNumber,
            const int maximumDifferentialCorrections,
            const double maximumPositionDeviation,
            const double maximumVelocityDeviation ,
            const double maximumEigenvalueDeviation ):
        massParameter_( massParameter ), orbitType_( orbitType ), librationPointNumber_( librationPointNumber ),
        maximumDifferentialCorrections_( maximumDifferentialCorrections ),
        maximumPositionDeviation_( maximumPositionDeviation ),
        maximumVelocityDeviation_( maximumVelocityDeviation ),
        maximumEigenvalueDeviation_( maximumEigenvalueDeviation ){ }

    double massParameter_;
    CR3BPPeriodicOrbitTypes orbitType_;
    int librationPointNumber_;

    int maximumDifferentialCorrections_;
    double maximumPositionDeviation_;
    double maximumVelocityDeviation_;
    double maximumEigenvalueDeviation_;

};


struct CR3BPPeriodicOrbitConditions
{
    CR3BPPeriodicOrbitConditions(
            const Eigen::Vector6d initialState,
            const Eigen::Vector6d halfPeriodStateVector,
            const double orbitPeriod,
            const int numberOfDifferentialCorrections ):
        initialState_( initialState ),
        halfPeriodStateVector_( halfPeriodStateVector ),
        orbitPeriod_( orbitPeriod ),
        numberOfDifferentialCorrections_( numberOfDifferentialCorrections ){ }

    Eigen::Vector6d initialState_;
    Eigen::Vector6d halfPeriodStateVector_;
    double orbitPeriod_;
    int numberOfDifferentialCorrections_;
};

std::pair< Eigen::Vector6d, double >  richardsonApproximationLibrationPointPeriodicOrbit(
        const double massParameter,
        const CR3BPPeriodicOrbitTypes orbitType,
        int librationPointNr, double amplitude, double n = 1.0 );

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

CR3BPPeriodicOrbitConditions createCR3BPPeriodicOrbit(
        const Eigen::Vector6d& initialStateVector,
        double orbitalPeriod,
        double massParameter,
        const CR3BPPeriodicOrbitGenerationSettings& periodicOrbitSettings,
        const std::shared_ptr< tudat::numerical_integrators::IntegratorSettings< double > > integratorSettings );

//bool checkForMonodromyUnitEigenvalues( const Eigen::Matrix6d &monodromyMatrix );

//bool terminateCR3BPNumericalContinuation(
//        const Eigen::Matrix6d& monodromyMatrix,
//        const Eigen::Vector6d& cartesianState,
//        const int numberOfOrbitsGenerated  );

//std::vector< CR3BPPeriodicOrbitConditions > createCR3BPPeriodicOrbitsThroughNumericalContinuation(
//        const std::shared_ptr< tudat::numerical_integrators::IntegratorSettings< double > > integratorSettings,
//        const CR3BPPeriodicOrbitConditions& firstPeriodicOrbit,
//        const CR3BPPeriodicOrbitConditions& secondPeriodicOrbit,
//        const CR3BPPeriodicOrbitGenerationSettings periodicOrbitSettings,
//        const std::function< double( const Eigen::Vector6d& ) > pseudoArcLengthFunction );

} // namespace propagators

} // namespace tudat

#endif // TUDAT_CR3BPPERIODICORBITS_H
