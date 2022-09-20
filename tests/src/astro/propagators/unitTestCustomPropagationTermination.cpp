/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "tudat/simulation/simulation.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat;
using namespace propagators;
using namespace mathematical_constants;
using namespace root_finders;

BOOST_AUTO_TEST_SUITE( test_custom_propagation_termination )

class PositionTerminationCondition: public CustomPropagationTerminationCondition< Eigen::Vector6d, double >
{
public:
    PositionTerminationCondition(
            const std::shared_ptr< root_finders::RootFinderSettings > terminationRootFinderSettings ):
        CustomPropagationTerminationCondition< Eigen::Vector6d, double >( terminationRootFinderSettings ){ }

    bool terminatePropagationDerived(
            const Eigen::Vector6d& currentState, const double currentTime )
    {
        if( computeStopConditionError( currentState, currentTime ) > 0.0 )
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    double computeStopConditionError(
            const Eigen::Vector6d& currentState, const double currentTime )
    {
        return currentState( 0, 0 ) - 1.5;
    }
};

class TimeTerminationCondition: public CustomPropagationTerminationCondition< Eigen::Vector6d, double >
{
public:
    TimeTerminationCondition(
            const std::shared_ptr< root_finders::RootFinderSettings > terminationRootFinderSettings ):
        CustomPropagationTerminationCondition< Eigen::Vector6d, double >( terminationRootFinderSettings ){ }

    bool terminatePropagationDerived(
            const Eigen::Vector6d& currentState, const double currentTime )
    {
        if( computeStopConditionError( currentState, currentTime ) > 0.0 )
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    double computeStopConditionError(
            const Eigen::Vector6d& currentState, const double currentTime )
    {
        return currentTime - 0.500102;
    }
};

BOOST_AUTO_TEST_CASE( testCR3BPPropagationExactTermination )
{
    using namespace tudat::propagators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::basic_mathematics;
    using namespace tudat::numerical_integrators;
    using namespace tudat::orbital_element_conversions;

    // Define initial normalized state for each case
    Eigen::Vector6d massLessBodyInitialCartesianState;
    massLessBodyInitialCartesianState[ xCartesianPositionIndex ] = 0.994;
    massLessBodyInitialCartesianState[ yCartesianPositionIndex ] = 0.853;
    massLessBodyInitialCartesianState[ zCartesianPositionIndex ] = 0.312;
    massLessBodyInitialCartesianState[ xCartesianVelocityIndex ] = 0.195;
    massLessBodyInitialCartesianState[ yCartesianVelocityIndex ] = -0.211;
    massLessBodyInitialCartesianState[ zCartesianVelocityIndex ] = 0.15;

    // Set normalized-mass parameter and propagation intervals
    double massParameter = 2.528e-5;
    double simulationStartEpoch = 0.0;
    double  timeStep = 0.0001;

    // Set integrator settings
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, simulationStartEpoch, timeStep );

    for( int terminateExactly = 0; terminateExactly < 2; terminateExactly++ )
    {
        bool terminateOnFinalCondition = static_cast< bool >( terminateExactly );
        for( int testCase = 0; testCase < 2; testCase++ )
        {
            std::cout<<"Case: "<<terminateExactly<<" "<<testCase<<std::endl;
            std::shared_ptr< CustomPropagationTerminationCondition< Eigen::Vector6d, double > > customTermination;

            if( testCase == 0 )
            {
                customTermination = std::make_shared< PositionTerminationCondition >(
                            bisectionRootFinderSettings(
                                1.0E-9, 1.0E-9, 1.0E-12 ) );
            }
            else
            {
                customTermination = std::make_shared< TimeTerminationCondition >(
                            bisectionRootFinderSettings(
                                1.0E-15, 1.0E-15, 1.0E-15 ) );
            }

            std::map< double, Eigen::Vector6d > stateHistory = performCR3BPIntegration(
                        integratorSettings, massParameter,
                        massLessBodyInitialCartesianState,
                        customTermination,
                        terminateOnFinalCondition, true );

            auto it = stateHistory.rbegin( );
            double lastTime = it->first;
            Eigen::Vector6d lastState = it->second;

            it++;
            double secondToLastTime = it->first;
            Eigen::Vector6d secondToLastState = it->second;

            if( testCase == 0 )
            {
                if( terminateOnFinalCondition )
                {
                    BOOST_CHECK_EQUAL( secondToLastState( 0, 0 ) < 1.5, true );
                    BOOST_CHECK_SMALL( std::fabs( lastState( 0, 0 ) - 1.5 ), 1.0E-9 );
                }
                else
                {
                    BOOST_CHECK_EQUAL( secondToLastState( 0, 0 ) < 1.5, true );
                    BOOST_CHECK_EQUAL( lastState( 0, 0 ) > 1.5, true );
                }
            }
            else
            {
                double terminationTime = 0.500102;
                if( terminateOnFinalCondition )
                {
                    BOOST_CHECK_EQUAL( secondToLastTime < terminationTime, true );
                    BOOST_CHECK_SMALL( std::fabs( lastTime - terminationTime ), 1.0E-14 );
                }
                else
                {
                    BOOST_CHECK_EQUAL( secondToLastTime < terminationTime, true );
                    BOOST_CHECK_EQUAL( lastTime > terminationTime, true );
                }
            }
        }
    }
}


BOOST_AUTO_TEST_SUITE_END( )

}

}
