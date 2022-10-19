/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_MPC_FILE_READERS_H
#define TUDAT_MPC_FILE_READERS_H

#include <string>
#include <map>
#include <fstream>
#include <sstream>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/trim_all.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/lexical_cast.hpp>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/astro/basic_astro/celestialBodyConstants.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/basic_astro/stateRepresentationConversions.h"
#include "tudat/io/matrixTextFileReader.h"
#include "tudat/io/streamFilters.h"
#include "tudat/basics/timeType.h"

namespace tudat
{
namespace input_output
{

struct MpcObservationRecords
{
    MpcObservationRecords( const std::string& observationLine )
    {
        std::string currentSubstring;
        note1_ = observationLine.at( 13 );
        note2_ = observationLine.at( 14 );
        int year = boost::lexical_cast< int >( observationLine.substr( 15, 4 ) );
        int month = boost::lexical_cast< int >( observationLine.substr( 20, 2 ) );
        int day = boost::lexical_cast< int >( observationLine.substr( 23, 2 ) );

        currentSubstring = observationLine.substr( 25, 7 );

        boost::trim( currentSubstring );

        double dayFraction = boost::lexical_cast< double >( currentSubstring );

        observationTime_ = basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch< long double >(
                    year, month, day, dayFraction, basic_astrodynamics::JULIAN_DAY_ON_J2000 ) *
                physical_constants::getJulianDay< long double >( );

        currentSubstring = observationLine.substr( 38, 6 );
        boost::trim( currentSubstring );
        rightAscenscion_ = unit_conversions::convertDegreesToRadians(
                ( boost::lexical_cast< int >( observationLine.substr( 32, 2 ) ) * 15 +
                boost::lexical_cast< int >( observationLine.substr( 35, 2 ) ) / 60 +
                boost::lexical_cast< double >( currentSubstring ) / 3600 ) );

        std::string declinationSign = observationLine.substr( 44, 1 );

        currentSubstring = observationLine.substr( 51, 5 );
        boost::trim( currentSubstring );
        declination_ = unit_conversions::convertDegreesToRadians(
                ( boost::lexical_cast< int >( observationLine.substr( 45, 2 ) ) +
                boost::lexical_cast< int >( observationLine.substr( 48, 2 ) ) / 60 +
                boost::lexical_cast< double >( currentSubstring ) / 3600 ) );

        if( declinationSign == "-" )
        {
            declination_ *= -1.0;
        }
        else if( declinationSign != "+" )
        {
            throw std::runtime_error( "Error, did not recognize declination sign in MPC data line: " + observationLine );
        }

        magnitudeBand_ = observationLine.substr( 65, 5 );
        stationId_ = boost::lexical_cast< std::string >( observationLine.substr( 77, 3 ) );
        std::cout<<observationTime_/86400.0/365.25<<" "<<rightAscenscion_<<" "<<declination_<<" "<<stationId_<<std::endl;

    }

    long double  observationTime_;
    double rightAscenscion_;
    double declination_;
    std::string stationId_;

    std::string magnitudeBand_;
    std::string reductionCatalog_;
    std::string note1_;
    std::string note2_;
};

struct MpcStationInformation
{
public:
    MpcStationInformation(
            const int stationIndex,
            const double longitude,
            const double radiusCosineLatitude,
            const double radiusSineLatitude,
            const std::string stationName ):
        stationIndex_( stationIndex ), stationName_( stationName )
    {
        double stationLatitude = std::atan2( radiusSineLatitude, radiusCosineLatitude );
        sphericalPosition_ =
                ( Eigen::Vector3d( ) << celestial_body_constants::EARTH_EQUATORIAL_RADIUS,
                                   stationLatitude, longitude ).finished( ) ;
        cartesianPosition_ =
                coordinate_conversions::convertPositionElements(
                    sphericalPosition_, coordinate_conversions::spherical_position, coordinate_conversions::cartesian_position );
    }

    const int stationIndex_;
    const std::string stationName_;

    Eigen::Vector3d sphericalPosition_;
    Eigen::Vector3d cartesianPosition_;

};

std::vector< MpcStationInformation > readMpcStationFile(
        const std::string& fileName )
{
    // Open file
    std::fstream stream( fileName.c_str( ), std::ios::in );
    if ( stream.fail( ) )
    {
        throw std::runtime_error( "MPC station file could not be opened: " + fileName );
    }

    // Read the filtered stream into lines.
    std::vector< std::string > lines;
    while ( !stream.fail( ) && !stream.eof( ) )
    {
        std::string line;
        getline( stream, line );
        if ( !line.empty( ) )
        {
            boost::trim_all( line );
            lines.push_back( line );
        }
    }

    std::vector< MpcStationInformation > stationDefinitions;

    int stationIndex;
    std::string stationName;
    std::string currentSubstring;
    for ( unsigned int rowIndex = 0; rowIndex < lines.size( ); rowIndex++ )
    {

        std::string currentLine = lines.at( rowIndex );
        try
        {
            currentSubstring = currentLine.substr( 0, 6 ) ;
            boost::trim( currentSubstring );
            stationIndex = boost::lexical_cast< int >( currentSubstring );

            currentSubstring = currentLine.substr( 6, 8 ) ;
            boost::trim( currentSubstring );
            double longitude = boost::lexical_cast< double >( currentSubstring );

            currentSubstring = currentLine.substr( 14, 9 ) ;
            boost::trim( currentSubstring );
            double radiusCosineLatitude = boost::lexical_cast< double >( currentSubstring );

            currentSubstring = currentLine.substr( 23, 10 ) ;
            boost::trim( currentSubstring );
            double radiusSineLatitude = boost::lexical_cast< double >( currentSubstring );

            currentSubstring = currentLine.substr( 33, currentLine.size( ) - 33 ) ;
            boost::trim( currentSubstring );
            stationName = boost::lexical_cast< int >( currentSubstring );

            stationDefinitions.push_back(
                        MpcStationInformation(
                            stationIndex, longitude, radiusCosineLatitude, radiusSineLatitude, stationName ) );
        }
        catch( ... ){}
    }
    return stationDefinitions;
}


std::map< std::string, std::vector< MpcObservationRecords > > readMpcObservations(
        const std::string& fileName )
{
    // Open file
    std::fstream stream( fileName.c_str( ), std::ios::in );
    if ( stream.fail( ) )
    {
        throw std::runtime_error( "MPC data file could not be opened: " + fileName );
    }

    // Read the filtered stream into lines.
    std::vector< std::string > lines;
    while ( !stream.fail( ) && !stream.eof( ) )
    {
        std::string line;
        getline( stream, line );
        if ( !line.empty( ) )
        {
            lines.push_back( line );
        }
    }

    std::map< std::string, std::vector< MpcObservationRecords > > stationDefinitions;

    for ( unsigned int rowIndex = 0; rowIndex < lines.size( ); rowIndex++ )
    {

        std::string currentLine = lines.at( rowIndex );

        MpcObservationRecords currentRecord( currentLine );
        stationDefinitions[ currentRecord.stationId_ ].push_back( currentRecord );
    }
    return stationDefinitions;
}

} // namespace input_output
} // namespace tudat

#endif // TUDAT_MPC_FILE_READERS_H
