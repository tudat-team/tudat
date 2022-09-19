/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CR3BPMANIFOLDS_H
#define TUDAT_CR3BPMANIFOLDS_H

#include <string>
#include <boost/bind/bind.hpp>

#include "tudat/simulation/propagation_setup/createCR3BPPeriodicOrbits.h"

namespace tudat
{

namespace propagators
{


double determineEigenvectorSign( const Eigen::Vector6d& eigenvector );

void determineStableUnstableEigenvectors(
        Eigen::Vector6d& stableEigenvector,
        Eigen::Vector6d& unstableEigenvector,
        const Eigen::MatrixXd& monodromyMatrix,
        double maxEigenvalueDeviation );


struct CR3BPManifoldSettings
{
    CR3BPManifoldSettings(
            const Eigen::Vector6d& stableEigenvector,
            const Eigen::Vector6d& unstableEigenvector,
            const double eigenvectorDisplacementFromOrbit ):
        stableEigenvector_( stableEigenvector ),
        unstableEigenvector_( unstableEigenvector ),
        eigenvectorDisplacementFromOrbit_( eigenvectorDisplacementFromOrbit )
    {
        stableEigenvectorSign = determineEigenvectorSign( stableEigenvector );
        unstableEigenvectorSign = determineEigenvectorSign( unstableEigenvector );
        offsetSigns_  = { 1.0 * stableEigenvectorSign, -1.0 * stableEigenvectorSign,
                          1.0 * unstableEigenvectorSign, -1.0 * unstableEigenvectorSign };

    }

    Eigen::Vector6d getEigenVectorForInitialDisplacement( const int manifoldIndex ) const
    {
        return ( manifoldIndex < 2 ) ? stableEigenvector_ : unstableEigenvector_;
    }

    double getEigenVectorDisplacementScaling( const int manifoldIndex ) const
    {
        return eigenvectorDisplacementFromOrbit_ * offsetSigns_.at( manifoldIndex );
    }

    double getIntegrationDirection( const int manifoldIndex ) const
    {
        return integrationDirections_.at( manifoldIndex );
    }

    Eigen::Vector6d stableEigenvector_;
    Eigen::Vector6d unstableEigenvector_;
    double eigenvectorDisplacementFromOrbit_;

    double stableEigenvectorSign;
    double unstableEigenvectorSign;
    std::vector< double > offsetSigns_;
    std::vector< double > integrationDirections_ = { -1.0, -1.0, 1.0, 1.0 };

};


void computeManifoldSetFromSinglePoint(
        std::vector< std::map< double, Eigen::Vector6d > >& manifoldStateHistories,
        const Eigen::MatrixXd& stateIncludingStm,
        const CR3BPPeriodicOrbitConditions periodicOrbitConditions,
        const CR3BPManifoldSettings manifoldSettings,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings );

std::vector< double > createManifoldDeparturePoints(
        const double orbitalPeriod,
        const int numberOfDeparturePoints,
        const std::map< double, Eigen::MatrixXd >& stateTransitionMatrixHistory );

void computeManifolds( std::vector< std::vector< std::map< double, Eigen::Vector6d > > >& fullManifoldStateHistories,
                       const CR3BPPeriodicOrbitConditions periodicOrbitConditions,
                       const double eigenvectorDisplacementFromOrbit,
                       const int numberOfDeparturePoints,
                       const double maxEigenvalueDeviation,
                       const std::shared_ptr< tudat::numerical_integrators::IntegratorSettings< double > > integratorSettings );



//class ManifoldTerminationConditionU1U4
//{
//public:
//    ManifoldTerminationConditionU1U4( )
//    {
//        ySignSet_ = false;
//        ySign = TUDAT_NAN;
//    }

//    bool isPropagationTerminated(
//            const Eigen::MatrixXd& currentState, const double currentTime )
//    {
//        // Determine sign of y when crossing x = 0  (U1, U4)
//        if ( (currentState( 0, 0 ) < 0 ) && !ySignSet_ )
//        {
//            if ( currentState( 1, 0 ) < 0 )
//            {
//                ySign = -1.0;
//            }
//            else if ( currentState( 1, 0 ) > 0 )
//            {
//                ySign = 1.0;
//            }
//            else
//            {
//                throw std::runtime_error( "Error when determining sign of y when crossing x=0 for U1 and U4 manifolds" );
//            }
//            ySignSet_ = true;
//        }

//        if ( (currentState( 1, 0 ) * ySign < 0 ) && ySignSet_ )
//        {
//            return true;
//        }
//        else
//        {
//            return false;
//        }
//    }

//private:

//    bool ySignSet_;

//    double ySign;
//};

//class ManifoldTerminationConditionU2U3
//{
//public:
//    ManifoldTerminationConditionU2U3(
//            const double massParameter ):
//        massParameter_( massParameter )
//    {
//        xDiffSignSet_ = false;
//        xDiffSign_ = TUDAT_NAN;
//    }

//    bool isPropagationTerminated(
//            const Eigen::MatrixXd& currentState, const double currentTime )
//    {
//        // Determine whether the trajectory approaches U2, U3 from the right or left (U2, U3)
//        if ( !xDiffSignSet_ )
//        {
//            if ( ( currentState(0, 0) - (1.0 - massParameter_ ) ) < 0 )
//            {
//                xDiffSign_ = -1.0;
//            }
//            else if ( ( currentState(0, 0) - (1.0 - massParameter_ ) ) > 0 )
//            {
//                xDiffSign_ = 1.0;
//            }
//            else
//            {
//                throw std::runtime_error( "Error when determining whether the trajectory approaches U2, U3 from the right or left" );
//            }
//        }


//        if ( (currentState( 1, 0 ) * ySign < 0 ) && ySignSet_ )
//        {
//            return true;
//        }
//        else
//        {
//            return false;
//        }
//    }

//private:

//    double massParameter_;

//    bool xDiffSignSet_;

//    double xDiffSign_;
//};


} // namespace propagators

} // namespace tudat

#endif // TUDAT_CR3BPMANIFOLDS_H
