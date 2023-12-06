/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/io/readMpcFile.h"

namespace tudat
{
namespace input_output
{


double MpcFileContents::getAngleInRadiansFromString( const std::string& angleString, const bool isNegative, const bool isRightAscension )
{
    double minute = std::stod( angleString.substr( 3, 2 ) );
    double seconds = std::stod( angleString.substr( 6, 6 ) );
    if( isRightAscension )
    {
        double hour = std::stod( angleString.substr( 0, 2 ) );
        return 2.0 * mathematical_constants::PI * ( hour + ( minute + seconds / 60.0 ) / 60.0 ) / 24.0;
    }
    else
    {
        double hour = std::stod( angleString.substr( 0, 2 ) );
        return 2.0 * mathematical_constants::PI * ( hour + ( minute + seconds / 60.0 ) / 60.0 ) / 360.0 * ( isNegative ? -1.0 : 1.0 );
    }
}

void MpcFileContents::readMpcFile( const std::string fileName, const MPCObjectType objectType )
{
    std::ifstream dataFile_;
    dataFile_.open( fileName.c_str( ) );

    // Check if file could be opened. Throw exception with error message if file could not be
    // opened.
    if ( !dataFile_.is_open( ) )
    {
        throw std::runtime_error( "MPC Data file could not be opened. Provided file path is:" + fileName );
    }

    // Read file line-by-line
    std::string currentLine;
    int lineCounter = 0;

    std::string bodyIdentifier = "", bodyProvisionalIdentifier = "", note1 = "", note2 = "";
    std::string yearString, monthString, dayString, rightAscensionString, declinationSignString, declinationString, magnitudeAndBand, catalogId, observatoryCode;
    int year, month, fullDay;
    double day, fractionalDay;
    double secondsSinceJ2000Utc;
    bool declinationIsNegative;
    double rightAscension, declination;
    while( !dataFile_.eof( ) )
    {
        // Get next line of data from data file and store in a string
        std::getline( dataFile_, currentLine );
        try
        {
            if( currentLine.length( ) != 80 )
            {
                std::cerr<<"Warning, MPC line " + std::to_string( lineCounter ) + ", which reads:" + currentLine +": has wrong number of characters: " + std::to_string( currentLine.length( ) ) <<std::endl;
            }

            bodyIdentifier = currentLine.substr( 0, 5 );
            bodyProvisionalIdentifier = currentLine.substr( 5, 7 );
            note1 = currentLine.at( 13 );
            note2 = currentLine.at( 14 );

            boost::trim( bodyIdentifier );
            boost::trim( bodyProvisionalIdentifier );
            boost::trim( note1 );
            boost::trim( note2 );

            if( !( ( note2 == "V" ) || ( note2 == "v" ) || ( note2 == "R" ) || ( note2 == "r" ) || ( note2 == "S" ) || ( note2 == "s" ) ) )
            {
                if( objectType == mpc_natural_satellite )
                {
                    if( bodyIdentifier.size( ) == 0 )
                    {
                        throw std::runtime_error(
                            "Error when reading MPC file, satellite name of line " + std::to_string( lineCounter ) +
                            " is malformed (size 0), line is:" + currentLine );
                    }

                    if( bodyIdentifier.back( ) != 'S' )
                    {
                        throw std::runtime_error(
                            "Error when reading MPC file, satellite name of line " + std::to_string( lineCounter ) +
                            " is malformed (does not end in S), line is:" + currentLine );
                    }
                    bodyIdentifier.pop_back();
                }

                if( bodyIdentifier.size( ) == 0 && bodyProvisionalIdentifier.size( ) == 0 )
                {
                    throw std::runtime_error(
                        "Error when reading MPC file, no body name defined in line " + std::to_string( lineCounter ) +
                        ", line is:" + currentLine );
                }

                yearString = currentLine.substr( 15, 4 );
                monthString = currentLine.substr( 20, 2 );
                dayString = currentLine.substr( 23, 9 );

                boost::trim( yearString );
                boost::trim( monthString );
                boost::trim( dayString );

                year = std::stoi( yearString );
                month = std::stoi( monthString );
                day = std::stod( dayString );
                fullDay = std::floor( day );
                fractionalDay = day - static_cast< double >( fullDay );

                secondsSinceJ2000Utc = basic_astrodynamics::timeFromDecomposedDateTime<double>( year, month, day, 0, 0,
                                                                                                physical_constants::JULIAN_DAY *
                                                                                                fractionalDay );

                rightAscensionString = currentLine.substr( 32, 12 );
                boost::trim( rightAscensionString );


                declinationSignString = currentLine.at( 44 );
                if ( declinationSignString == "+" )
                {
                    declinationIsNegative = false;
                }
                else if ( declinationSignString == "-" )
                {
                    declinationIsNegative = true;
                }
                else
                {
                    throw std::runtime_error(
                        "Error when reading MPC file, declination sign of line " + std::to_string( lineCounter ) +
                        " is " + declinationSignString );
                }
                declinationString = currentLine.substr( 45, 11 );
                boost::trim( declinationString );

                rightAscension = getAngleInRadiansFromString( rightAscensionString, false, true );
                declination = getAngleInRadiansFromString( declinationString, declinationIsNegative, false );

                magnitudeAndBand = currentLine.substr( 65, 6 );
                catalogId = currentLine.at( 71 );
                observatoryCode = currentLine.substr( 77, 3 );

                if( bodyIdentifier.size( ) > 0 )
                {
                    bodyNameEntries_[ bodyIdentifier ].push_back( mpcFileLines_.size( ) );
                }
                else
                {
                    bodyProvisionalNameEntries_[ bodyProvisionalIdentifier ].push_back( mpcFileLines_.size( ) );
                }

                if( bodyProvisionalIdentifier.size( ) > 0 && bodyIdentifier.size( ) > 0 )
                {
                    if( provisionalToProperNames_.count( bodyProvisionalIdentifier ) == 0 )
                    {
                        provisionalToProperNames_[ bodyProvisionalIdentifier ] = bodyIdentifier;
                    }
                }
                mpcFileLines_.push_back(
                    MPCFileSingleLine( secondsSinceJ2000Utc, rightAscension, declination, observatoryCode,
                                       bodyIdentifier, bodyProvisionalIdentifier, catalogId,
                                       magnitudeAndBand, note1, note2 ));
            }
        }
        catch( ... )
        {
            std::cerr<<"Could not read line number "<<lineCounter<<" with contents:"<<currentLine<<std::endl;
        }
        lineCounter++;
    }

    postProcessProvisionalBodyEntries( );
}

void MpcFileContents::postProcessProvisionalBodyEntries( )
{
    std::set< std::string > provisionalNameList = utilities::vectorToSet(
        utilities::createVectorFromMapKeys( bodyProvisionalNameEntries_ ) );

    for( auto it : provisionalToProperNames_ )
    {
        if( bodyNameEntries_.count( it.second ) == 0 )
        {
            bodyNameEntries_[ it.second ] = bodyProvisionalNameEntries_.at( it.first );
        }
        else if( bodyProvisionalNameEntries_.count( it.first ) > 0 )
        {
            std::vector< int > indexList = bodyNameEntries_.at( it.second );
            std::vector< int > indexListProvisional = bodyProvisionalNameEntries_.at( it.first );

            indexList.insert( indexList.end( ), indexListProvisional.begin( ), indexListProvisional.end( ) );
            std::sort( indexList.begin( ), indexList.end( ) );
            bodyNameEntries_[ it.second ] = indexList;
        }

        if( bodyProvisionalNameEntries_.count( it.first ) > 0 )
        {
            provisionalNameList.erase( it.first );
        }
    }

    provisionalOnlyNames_ = utilities::setToVector( provisionalNameList );
}


} // namespace input_output

} // namespace tudat
