/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

//#define BOOST_TEST_DYN_LINK
//#define BOOST_TEST_MAIN

#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <vector>

#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "tudat/basics/testMacros.h"
#include "tudat/io/streamFilters.h"

#include "tudat/io/mpcFileReading.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/io/aerodynamicCoefficientReader.h"

using namespace tudat::input_output;
// Test if multi-array file reader is working correctly
//BOOST_AUTO_TEST_CASE( testMultiArrayReader )
int main( )
{
    std::string directory = "/home/dominic/Software/tudat-bundle/tudat-bundle/tudat/src/io/";
    std::string fileName = "karelwakker.txt";


    std::map< std::string, std::vector< MpcObservationRecords > > mpcData = readMpcObservations(
                directory + fileName );
    std::cout<<mpcData.size( )<<std::endl;
    int counter = 0;
    for( auto it : mpcData )
    {
        std::map< double, Eigen::Vector2d > perStationData;
        for( int i = 0; i < it.second.size( ); i++ )
        {
            MpcObservationRecords currentRecord = it.second.at( i );
            perStationData[ currentRecord.observationTime_ ] =
                    ( Eigen::Vector2d( ) << currentRecord.rightAscenscion_, currentRecord.declination_  ).finished( );
        }
        std::string outputFile = "testData" + std::to_string( counter ) + ".dat" ;

        writeDataMapToTextFile( perStationData, outputFile, "/home/dominic/Software/tudat-bundle/tudat-bundle/tudat/src/io/" );
        counter++;

    }
}

//BOOST_AUTO_TEST_SUITE_END( )

//} // namespace unit_tests
//} // namespace tudat


