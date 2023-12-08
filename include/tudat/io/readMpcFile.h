/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References: https://www.minorplanetcenter.net/iau/info/OpticalObs.html
 *
 */

#ifndef TUDAT_READ_MPC_FILE_H
#define TUDAT_READ_MPC_FILE_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <map>

#include <boost/algorithm/string/trim.hpp>

#include "tudat/astro/basic_astro/dateTime.h"
#include "tudat/io/basicInputOutput.h"

namespace tudat
{
namespace input_output
{

enum MPCObjectType
{
    mpc_comet,
    mpc_natural_satellite,
    mpc_minor_planet
};

struct MPCFileSingleLine
{
public:
    MPCFileSingleLine(
        const double secondsSinceJ2000Utc,
        const double rightAscension,
        const double declination,
        const std::string& observatoryCode,
        const std::string& bodyIdentifier = "",
        const std::string& bodyProvisionalIdentifier = "",
        const std::string& catalogId = "",
        const std::string& magnitudeAndBand = "",
        const std::string& note1 = "",
        const std::string& note2 = "",
        const bool isDiscovery = false ):
        secondsSinceJ2000Utc_( secondsSinceJ2000Utc ),
        rightAscension_( rightAscension ),
        declination_( declination ),
        observatoryCode_( observatoryCode ),
        bodyIdentifier_( bodyIdentifier ),
        bodyProvisionalIdentifier_( bodyProvisionalIdentifier ),
        catalogId_( catalogId ),
        magnitudeAndBand_( magnitudeAndBand ),
        note1_( note1 ),
        note2_( note2 ),
        isDiscovery_( isDiscovery )
        {
            if( bodyIdentifier_ == "" && bodyProvisionalIdentifier_ == "" )
            {
                throw std::runtime_error( "Error when creating MPC file line contents; no body ID or provisional ID provided" );
            }
        }

    const double secondsSinceJ2000Utc_;
    const double rightAscension_;
    const double declination_;
    const std::string observatoryCode_;
    const std::string bodyIdentifier_;
    const std::string bodyProvisionalIdentifier_;
    const std::string catalogId_;
    const std::string magnitudeAndBand_;
    const std::string note1_;
    const std::string note2_;
    const bool isDiscovery_;

};

class MpcFileContents
{
public:

    MpcFileContents( const std::string fileName, const MPCObjectType objectType )
    {
        readMpcFile( fileName, objectType );
    }

    std::map< std::string, std::string > getProvisionalToProperNames( )
    {
        return provisionalToProperNames_;
    }

    std::map< std::string, std::vector< int > > getBodyNameEntries( )
    {
        return bodyNameEntries_;
    }

    std::vector< std::string > getProvisionalOnlyNames( )
    {
        return provisionalOnlyNames_;
    }

    void filterForStations(
        const std::vector< std::string >& stationsToRemove = std::vector< std::string >( ),
        const std::vector< std::string >& stationsToRetain = std::vector< std::string >( ) )
    {

    }

    void filterWithTime( const double startEpoch = TUDAT_NAN, const double endEpoch = TUDAT_NAN )
    {

    }

    void getObservationsForBody(
        const std::string& bodyName,
        std::map< std::string, std::vector< std::pair< Eigen::Vector2d, double > > > angularPositionsPerStation ) const
    {
        if( bodyNameEntries_.count( bodyName ) == 0 )
        {
            throw std::runtime_error( "Error when retrieving MPC data for " + bodyName + " no such body found." );
        }
        std::vector< int > requestedBodyEntries = bodyNameEntries_.at( bodyName );

        for( unsigned int i = 0; i < requestedBodyEntries.size( ); i++ )
        {
            const MPCFileSingleLine& currentLine = mpcFileLines_.at( requestedBodyEntries.at( i ) );
            angularPositionsPerStation[ currentLine.observatoryCode_ ].push_back( std::make_pair(
                ( Eigen::Vector2d( ) << currentLine.rightAscension_, currentLine.declination_ ).finished( ),
                currentLine.secondsSinceJ2000Utc_ ) );
        }
    }

    std::vector< std::string > getListOfBodies( ) const
    {
        return utilities::createVectorFromMapKeys( bodyNameEntries_ );
    }
private:

    double getAngleInRadiansFromString( const std::string& angleString, const bool isNegative, const bool isRightAscension );

    void readMpcFile( const std::string fileName, const MPCObjectType objectType );

    void postProcessProvisionalBodyEntries( );

    std::map< std::string, std::vector< int > > bodyNameEntries_;

    std::map< std::string, std::vector< int > > bodyProvisionalNameEntries_;

    std::vector< std::string > provisionalOnlyNames_;

    std::map< std::string, std::string > provisionalToProperNames_;

    std::vector< MPCFileSingleLine > mpcFileLines_;

};

} // namespace input_output

} // namespace tudat

#endif // TUDAT_READ_MPC_FILE_H
