/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References: 820-013, TRK-2-18 Tracking System Interfaces Orbit Data File Interface, Revision E, 2008, JPL/DSN
 */

#ifndef TUDAT_PROCESSMPCFILE_H
#define TUDAT_PROCESSMPCFILE_H

#include "tudat/io/readMpcFile.h"
#include "tudat/simulation/estimation_setup/observations.h"
#include "tudat/simulation/estimation_setup/observationSimulationSettings.h"

namespace tudat
{

namespace observation_models
{

template< typename ObservationScalarType = double, typename TimeType = double >
std::vector< std::shared_ptr< observation_models::SingleObservationSet< ObservationScalarType, TimeType > > > createMpcSingleObservationSets(
    const input_output::MpcFileContents& mpcFileContents )
{
    std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > timeScaleConverter = earth_orientation::createDefaultTimeConverter( );
    std::vector< std::shared_ptr< observation_models::SingleObservationSet< ObservationScalarType, TimeType > > > singleObservationSet;

    std::vector< std::string > bodies = mpcFileContents.getListOfBodies( );
    for( unsigned int i = 0; i < bodies.size( ); i++ )
    {
        std::map< std::string, std::vector< std::pair< Eigen::Vector2d, double > > > currentData;
        mpcFileContents.getObservationsForBody(
            bodies.at( i ), currentData );
        for( auto it : currentData)
        {
            std::string stationName = it.first;
            std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observations;
            std::vector< TimeType > observationTimes;
            for( unsigned int j = 0; j < it.second.size( ); j++ )
            {
                observations.push_back( it.second.at( j ).first.template cast< ObservationScalarType >( ) );
                observationTimes.push_back( static_cast< TimeType >(
                    timeScaleConverter->getCurrentTime< double >( basic_astrodynamics::utc_scale, basic_astrodynamics::tdb_scale, it.second.at( j ).second ) ) );
            }
            LinkEnds currentLinkEnds;
            currentLinkEnds[ transmitter ] = LinkEndId( bodies.at( i ), "" );
            currentLinkEnds[ receiver ] = LinkEndId( "Earth", stationName );

            LinkDefinition currentLinkDefinition( currentLinkEnds );
            singleObservationSet.push_back( std::make_shared< SingleObservationSet< ObservationScalarType, TimeType > >(
                angular_position, currentLinkDefinition,
                observations, observationTimes, receiver ) );
        }
    }

    return singleObservationSet;
}

} // namespace observation_models

} // namespace tudat

#endif // TUDAT_PROCESSMPCFILE_H
