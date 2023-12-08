/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

//#define BOOST_TEST_DYN_LINK
//#define BOOST_TEST_MAIN

#include <vector>

#include "tudat/basics/testMacros.h"
#include "tudat/io/readMpcFile.h"
#include "tudat/simulation/estimation_setup/processMpcFile.h"

//
//namespace tudat
//{
//namespace unit_tests
//{

using namespace tudat;
using namespace tudat::input_output;
using namespace tudat::observation_models;

//
//BOOST_AUTO_TEST_SUITE( test_mpc_file_reader )
//
//BOOST_AUTO_TEST_CASE( testSingleMpcFileReader )
int main( )
{
    std::string fileName = "/home/dominic/Downloads/SatObs.txt";
    MpcFileContents mpcFileContents( fileName, mpc_natural_satellite );

    std::map< std::string, std::string > provisionalToProperNames = mpcFileContents.getProvisionalToProperNames( );
    for( auto it : provisionalToProperNames )
    {
        std::cout<<it.first<<", "<<it.second<<std::endl;
    }

    std::map< std::string, std::vector< int > > bodyNameEntries = mpcFileContents.getBodyNameEntries( );
    for( auto it : bodyNameEntries )
    {
        std::cout<<it.first<<", "<<it.second.size( )<<std::endl;
    }

    std::vector< std::string > provisionalOnlyNames = mpcFileContents.getProvisionalOnlyNames( );
    for( unsigned int i = 0; i < provisionalOnlyNames.size( ); i++ )
    {
        std::cout<<provisionalOnlyNames.at( i )<<std::endl;
    }

    std::vector< std::shared_ptr< SingleObservationSet< > > > observationSets = createMpcSingleObservationSets( mpcFileContents );
    ObservationCollection< > observationCollection( observationSets );




}

//BOOST_AUTO_TEST_SUITE_END( )
//
//}
//
//}