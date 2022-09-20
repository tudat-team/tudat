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


template< typename StateType = Eigen::Vector6d >
class ManifoldTerminationConditionU1U4: public CustomPropagationTerminationCondition< StateType, double >
{
public:
    ManifoldTerminationConditionU1U4(
            const std::shared_ptr< root_finders::RootFinderSettings > terminationRootFinderSettings = nullptr ):
        CustomPropagationTerminationCondition< StateType, double >( terminationRootFinderSettings )
    {
        ySignSet_ = false;
        ySign = TUDAT_NAN;
    }

    bool terminatePropagationDerived(
            const Eigen::Vector6d& currentState, const double currentTime )
    {
        // Determine sign of y when crossing x = 0  (U1, U4)
        if ( (currentState( 0, 0 ) < 0 ) && !ySignSet_ )
        {
            if( currentState( 0, 0 ) == 0.0 )
            {  throw std::runtime_error( "Error when determining sign of y when crossing x=0 for U1 and U4 manifolds" ); }

            currentState( 1, 0 ) < 0 ? ySign = -1.0 : ySign = 1.0;
            ySignSet_ = true;
        }

        // Check if termination is reached
        if ( ( computeStopConditionError( currentState,currentTime  ) < 0 ) && ySignSet_ )
        {
            std::cout<<"U1U4 termination "<<currentTime<<std::endl;
            return true;
        }
        else
        {
            return false;
        }
    }

    double computeStopConditionError(
            const StateType& currentState, const double currentTime )
    {
        return currentState( 1, 0 ) * ySign;
    }

private:

    bool ySignSet_;

    double ySign;
};

template< typename StateType = Eigen::Vector6d >
class ManifoldTerminationConditionU2U3: public CustomPropagationTerminationCondition< StateType, double >
{
public:
    ManifoldTerminationConditionU2U3(
            const double massParameter,
            int librationPoint,
            int manifoldNumber,
            const std::shared_ptr< root_finders::RootFinderSettings > terminationRootFinderSettings = nullptr ):
        CustomPropagationTerminationCondition< StateType, double >( terminationRootFinderSettings ),
        massParameter_( massParameter ),
        librationPoint_( librationPoint ),
        manifoldNumber_( manifoldNumber )
    {
        xDiffSignSet_ = false;
        xDiffSign_ = TUDAT_NAN;
    }

    bool terminatePropagationDerived(
            const Eigen::Vector6d& currentState, const double currentTime )
    {
        // Determine whether the trajectory approaches U2, U3 from the right or left (U2, U3)
        if ( !xDiffSignSet_ )
        {
            ( ( currentState(0, 0) - (1.0 - massParameter_ ) ) < 0 ) ? xDiffSign_ = -1.0 : xDiffSign_ = 1.0;

            if ( ( currentState(0, 0) - (1.0 - massParameter_ ) ) == 0 )
            {  throw std::runtime_error( "Error when determining whether the trajectory approaches U2, U3 from the right or left" ); }
            xDiffSignSet_ = true;
        }


        if ( ( computeStopConditionError( currentState, currentTime ) < 0 ) &&
             ( ( librationPoint_ == 1 && ( manifoldNumber_ == 0 || manifoldNumber_ == 2) ) ||
               ( librationPoint_ == 2 && ( manifoldNumber_ == 1 || manifoldNumber_ == 3) ) ) )
        {
            std::cout<<"U2U3 termination "<<currentTime<<std::endl;
            return true;
        }
        else
        {
            return false;
        }
    }

    double computeStopConditionError(
            const StateType& currentState, const double currentTime )
    {
        return ( currentState(0, 0) - ( 1.0 - massParameter_ ) ) * xDiffSign_;
    }

private:

    double massParameter_;

    int librationPoint_;

    int manifoldNumber_;

    bool xDiffSignSet_;

    double xDiffSign_;
};


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


template< typename StateType = Eigen::Vector6d >
std::shared_ptr< MultipleCustomPropagationTerminationCondition< StateType, double > > getPoincareSectionsTerminationConditions(
        const std::shared_ptr< PropagatedCR3BPPeriodicOrbitConditions > periodicOrbitConditions,
        const int manifoldNumber )
{
    std::shared_ptr< root_finders::RootFinderSettings > terminationRootFinderSettings =
            root_finders::bisectionRootFinderSettings( 1.0E-14, 1.0E-14, 1.0E-14, 1000, root_finders::accept_result_with_warning );
    std::shared_ptr< ManifoldTerminationConditionU1U4< StateType > > u1u4Termination =
            std::make_shared< ManifoldTerminationConditionU1U4< StateType > >( terminationRootFinderSettings );
    std::shared_ptr< ManifoldTerminationConditionU2U3< StateType > > u2u3Termination =
            std::make_shared< ManifoldTerminationConditionU2U3< StateType > >(
                periodicOrbitConditions->massParameter_,
                periodicOrbitConditions->librationPoint_,
                manifoldNumber, terminationRootFinderSettings );
    return std::make_shared<  MultipleCustomPropagationTerminationCondition< StateType, double > >(
                std::vector< std::shared_ptr< CustomPropagationTerminationCondition< StateType, double > > >( { u1u4Termination, u2u3Termination } ) );
}

void computeManifoldSetFromSinglePoint(
        std::vector< std::map< double, Eigen::Vector6d > >& manifoldStateHistories,
        const Eigen::MatrixXd& stateIncludingStm,
        const std::shared_ptr< PropagatedCR3BPPeriodicOrbitConditions > periodicOrbitConditions,
        const CR3BPManifoldSettings manifoldSettings,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings );

std::vector< double > createManifoldDeparturePoints(
        const double orbitalPeriod,
        const int numberOfDeparturePoints,
        const std::map< double, Eigen::MatrixXd >& stateTransitionMatrixHistory );

void computeManifolds(std::vector< std::vector< std::map< double, Eigen::Vector6d > > >& fullManifoldStateHistories,
                      const std::shared_ptr< PropagatedCR3BPPeriodicOrbitConditions > periodicOrbitConditions,
                      const double eigenvectorDisplacementFromOrbit,
                      const int numberOfDeparturePoints,
                      const double maxEigenvalueDeviation,
                      const std::shared_ptr< tudat::numerical_integrators::IntegratorSettings< double > > integratorSettings ,
                      const std::shared_ptr<tudat::numerical_integrators::IntegratorSettings<double> > manifoldIntegratorSettings );




} // namespace propagators

} // namespace tudat

#endif // TUDAT_CR3BPMANIFOLDS_H
