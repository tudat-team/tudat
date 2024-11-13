/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_FREQUENCYOFARRIVALDIFFERENCEOBSERVATIONMODEL_H
#define TUDAT_FREQUENCYOFARRIVALDIFFERENCEOBSERVATIONMODEL_H

#include <map>
#include <Eigen/Core>

#include "tudat/math/basic/coordinateConversions.h"
#include "tudat/astro/observation_models/lightTimeSolution.h"
#include "tudat/astro/observation_models/oneWayDopplerObservationModel.h

namespace tudat
{

namespace observation_models
{

//! Class for simulating time of arrival difference observables.
template< typename ObservationScalarType = double, typename TimeType = double >
class FrequencyOfArrivalDifferenceObservationModel: public ObservationModel< 1, ObservationScalarType, TimeType >
{
public:

    typedef Eigen::Matrix< ObservationScalarType, 6, 1 > StateType;
    typedef Eigen::Matrix< ObservationScalarType, 6, 1 > PositionType;

    //! Constructor.
    /*!
     *  Constructor,
     *  \param dopplerModelFirstReceiver Object to compute the light-time (including any corrections w.r.t. Euclidean case)
     *  between source and first receiver
     *  \param dopplerModelSecondReceiver Object to compute the light-time (including any corrections w.r.t. Euclidean case)
     *  between source and second receiver
     *  \param observationBiasCalculator Object for calculating system-dependent errors in the
     *  observable, i.e. deviations from the physically ideal observable between reference points (default none).
     */
    FrequencyOfArrivalDifferenceObservationModel(
            const LinkEnds linkEnds,
            const std::shared_ptr< observation_models::OneWayDopplerObservationModel< ObservationScalarType, TimeType > > dopplerModelFirstReceiver,
            const std::shared_ptr< observation_models::OneWayDopplerObservationModel< ObservationScalarType, TimeType > > dopplerModelSecondReceiver,
            const std::shared_ptr< ObservationBias< 1 > > observationBiasCalculator = nullptr ):
        ObservationModel< 1, ObservationScalarType, TimeType >( time_difference_of_arrival, linkEnds, observationBiasCalculator ),
        dopplerModelFirstReceiver_( dopplerModelFirstReceiver ), dopplerModelSecondReceiver_( dopplerModelSecondReceiver ) { }

    //! Destructor
    ~FrequencyOfArrivalDifferenceObservationModel( ){ }

    //! Function to compute ideal TOA difference observation at given time, between two transmitters.
    /*!
     *  This function compute ideal TOA difference observation at a given time, between two transmitters. The time argument can be either the
     *  reception or transmission time (defined by linkEndAssociatedWithTime input).
     *  Note that this observable does include e.g. light-time corrections, which represent physically true corrections.
     *  It does not include e.g. system-dependent measurement.
     *  The times and states of the link ends are also returned in full precision (determined by class template
     *  arguments). These states and times are returned by reference.
     *  \param time Time at which observation is to be simulated
     *  \param linkEndAssociatedWithTime Link end at which given time is valid, i.e. link end for which associated time
     *  is kept constant (to input value)
     *  \param linkEndTimes List of times at each link end during observation (returned by reference).
     *  \param linkEndStates List of states at each link end during observation (returned by reference).
     *  \return Calculated TOA difference observable values.
     */
    Eigen::Matrix< ObservationScalarType, 1, 1 > computeIdealObservationsWithLinkEndData(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime,
            std::vector< double >& linkEndTimes,
            std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const std::shared_ptr< ObservationAncilliarySimulationSettings > ancilliarySetings = nullptr  )

    {
//        Eigen::Matrix< ObservationScalarType, 6, 1 > firstReceiverState;
//        Eigen::Matrix< ObservationScalarType, 6, 1 > secondReceiverState;
//        Eigen::Matrix< ObservationScalarType, 6, 1 > transmitterState;
//
//        // Compute light-times and receiver/transmitters states.
//        ObservationScalarType lightTimeFirstReceiver;
//        ObservationScalarType lightTimeSecondReceiver;
//
//        switch( linkEndAssociatedWithTime )
//        {
//        case transmitter:
//            lightTimeFirstReceiver = dopplerModelFirstReceiver_->calculateLightTimeWithLinkEndsStates(
//                firstReceiverState, transmitterState, time, false, ancilliarySetings );
//            lightTimeSecondReceiver = dopplerModelSecondReceiver_->calculateLightTimeWithLinkEndsStates(
//                secondReceiverState, transmitterState, time, false, ancilliarySetings );
//
//            linkEndTimes.push_back( static_cast< double >( time ) );
//            linkEndTimes.push_back( static_cast< double >( time + lightTimeFirstReceiver ) );
//            linkEndTimes.push_back( static_cast< double >( time + lightTimeSecondReceiver ) );
//
//            break;
//        case receiver1:
//            lightTimeFirstReceiver = dopplerModelFirstReceiver_->calculateLightTimeWithLinkEndsStates(
//                firstReceiverState, transmitterState, time, true, ancilliarySetings );
//            lightTimeSecondReceiver = dopplerModelSecondReceiver_->calculateLightTimeWithLinkEndsStates(
//                secondReceiverState, transmitterState, time - lightTimeFirstReceiver, false, ancilliarySetings );
//
//            linkEndTimes.push_back( static_cast< double >( time - lightTimeFirstReceiver ) );
//            linkEndTimes.push_back( static_cast< double >( time ) );
//            linkEndTimes.push_back( static_cast< double >( time - lightTimeFirstReceiver + lightTimeSecondReceiver ) );
//
//            break;
//        case receiver2:
//            lightTimeSecondReceiver = dopplerModelSecondReceiver_->calculateLightTimeWithLinkEndsStates(
//                secondReceiverState, transmitterState, time, true, ancilliarySetings );
//            lightTimeFirstReceiver = dopplerModelFirstReceiver_->calculateLightTimeWithLinkEndsStates(
//                firstReceiverState, transmitterState, time - lightTimeSecondReceiver, false, ancilliarySetings );
//            linkEndTimes.push_back( static_cast< double >( time - lightTimeSecondReceiver ) );
//            linkEndTimes.push_back( static_cast< double >( time - lightTimeSecondReceiver + lightTimeFirstReceiver ) );
//            linkEndTimes.push_back( static_cast< double >( time ) );
//
//            break;
//        default:
//            throw std::runtime_error( "Error when computing TOA difference observable, did not recognize reference link end " + getLinkEndTypeString( linkEndAssociatedWithTime ) );
//        }
//
//        // Set link end times and states.
//        linkEndTimes.clear( );
//        linkEndStates.clear( );
//
//        linkEndStates.push_back( transmitterState.template cast< double >( ) );
//        linkEndStates.push_back( firstReceiverState.template cast< double >( ) );
//        linkEndStates.push_back( secondReceiverState.template cast< double >( ) );
//
//        // Return observable
//        return ( Eigen::Matrix< ObservationScalarType, 1, 1 >( ) << lightTimeSecondReceiver - lightTimeFirstReceiver ).finished( );
    }

    //! Function to get the object to calculate light time between first transmitter and receiver.
    /*!
     * Function to get the object to calculate light time between first transmitter and receiver.
     * \return Object to calculate light time between first transmitter and receiver.
     */
    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > getDopplerModelFirstReceiver( )
    {
        return dopplerModelFirstReceiver_;
    }

    //! Function to get the object to calculate light time between second transmitter and receiver.
    /*!
     * Function to get the object to calculate light time between second transmitter and receiver.
     * \return Object to calculate light time between second transmitter and receiver.
     */
    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > getDopplerModelSecondReceiver( )
    {
        return dopplerModelSecondReceiver_;
    }

    LinkEnds getFirstLinkEnds( )
    {
        LinkEnds firstLinkEnds;
        firstLinkEnds[ transmitter ] = this->linkEnds_[ transmitter ];
        firstLinkEnds[ receiver ] = this->linkEnds_[ receiver1 ];
        return firstLinkEnds;
    }

    LinkEnds getSecondLinkEnds( )
    {
        LinkEnds secondLinkEnds;
        secondLinkEnds[ transmitter ] = this->linkEnds_[ transmitter ];
        secondLinkEnds[ receiver ] = this->linkEnds_[ receiver2 ];
        return secondLinkEnds;
    }
private:

    //! Object to calculate light time between the first transmitter and receiver.
    /*!
     *  Object to calculate light time between the first transmitter and receiver, including possible corrections from troposphere, relativistic corrections, etc.
     */
    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > dopplerModelFirstReceiver_;

    //! Object to calculate light time between the second transmitter and receiver.
    /*!
     *  Object to calculate light time between the second transmitter and receiver, including possible corrections from troposphere, relativistic corrections, etc.
     */
    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > dopplerModelSecondReceiver_;
};

} // namespace observation_models

} // namespace tudat

#endif // TUDAT_FREQUENCYOFARRIVALDIFFERENCEOBSERVATIONMODEL_H
