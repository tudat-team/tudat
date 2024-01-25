/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *
 */

//#define BOOST_TEST_DYN_LINK
//#define BOOST_TEST_MAIN

#include <limits>
#include "fstream"
#include "iostream"

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/astro/aerodynamics/customAerodynamicCoefficientInterface.h"
#include "tudat/astro/aerodynamics/aerodynamicAcceleration.h"
#include "tudat/astro/reference_frames/aerodynamicAngleCalculator.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"
#include "tudat/interface/spice/spiceEphemeris.h"
#include "tudat/interface/spice/spiceRotationalEphemeris.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/estimation_setup/createNumericalSimulator.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"

using namespace tudat::aerodynamics;


#include "tudat/astro/basic_astro/timeConversions.h"

#include "tudat/astro/aerodynamics/marsDtmAtmosphereModel.h"

int main( ) {
    //std::string filename = "/Users/ralkahal/Documents/PhD/atmodensitypds/dtm_mars";
    std::string filename = "/Users/ralkahal/OneDrive - Delft University of Technology/PhD/Programs/atmodensitydtm/dtm_mars";
    MarsDtmAtmosphereModel atmosphereModel = MarsDtmAtmosphereModel(3376.78E3, filename);
    //atmosphereModel.computelocalsolartime(0.0,16 ,12,2000, 0.0, 1.0697333);
    //std::cout<<atmosphereModel.computeGl(0.0,0.0,1.0697333,0.0,16,12,2000,1)<<std::endl;
    //std::cout<<std::abs((tudat::basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch(
    //       2021, 2, 7, 0.0, 0.0, 0.0, tudat::basic_astrodynamics::JULIAN_DAY_ON_J2000) -
    //tudat::basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch(2022,2,18,0.0,0.0,0.0,tudat::basic_astrodynamics::JULIAN_DAY_ON_J2000))*24/24.63)<<std::endl;

    //marsDate date2 = marsDate(2022, 2, 18, 0.0, 0.0, 0.0);
    //std::cout << date2.marsDayofYear(date2) << std::endl;
    //std::cout<<"Geopotential height "<< atmosphereModel.computeGeopotentialAltitude( 255.0E3 )<<std::endl;
    //std::cout<<atmosphereModel.computeCurrentTemperature( 0.0, 0.0, 1.0697333, 0.0, 16 ,12, 2000, 1)<<std::endl;
    //std::cout<<atmosphereModel.computeGamma( 0.0, 0.0, 1.0697333, 0.0, 16 ,12, 2000, 1)<<std::endl;
    //std::cout<< atmosphereModel.heightDistributionFunction(255.0E3, 0.0, 0.0, 1.0697333, 0.0, 16 ,12, 2000, 1)<<std::endl;
    std::ofstream outputFile(
            "/Users/ralkahal/OneDrive - Delft University of Technology/PhD/Programs/atmodensitydtm/density_output.txt");
    // Check if the file is opened successfully

    if (!outputFile.is_open()) {
        std::cerr << "Unable to open the file!" << std::endl;
        return 1; // return an error code
    }
/*
    for (int altitude = 138E3; altitude <= 1000E3; altitude += 10E3) {
        double alt_km = static_cast<double>(altitude);

        double rho = atmosphereModel.getTotalDensity(alt_km, 0.0, 0.0, 1.0697333, 0.0, 16, 12, 2000);

        // Write altitude and corresponding density to the file
        outputFile << alt_km << " " << rho << "\n";
    }
  */
    //for (int altitude = 138E3; altitude <= 1000E3; altitude += 10E3) {
      //  double alt_km = static_cast<double>(altitude);
      int alt_km = 138E3;
        for (int time = 0; time <= 5 * 24 * 60 * 60; time += 60) {

            int minutes = time % (60 * 60);
            int hours = (time / (60 * 60)) % 24;
            int days = ((time / (60 * 60 * 24)) % 31) + 1; // Adding 1 to ensure days are between 1 and 31
            int months = ((time / (60 * 60 * 24 * 31)) % 12) + 1; // Adding 1 to ensure months are between 1 and 12
            int years = (time / (60 * 60 * 24 * 31 * 12)) + 2000;

            double rho = atmosphereModel.getTotalDensity(alt_km, 0.0, 0.0, minutes, hours, days, months, years);
            double Ls = atmosphereModel.getSolarLongitude();
            double currentF107 = atmosphereModel.getSolarFluxIndex();
            // Write altitude and corresponding density to the file
            std::cout << alt_km << "altitude" << rho << "density" << Ls << "Ls" << currentF107 << "F107" << std::endl;
            outputFile << alt_km << " " << rho << "\n";

        //}
    }// Close the file when you are done
    outputFile.close();
    std::cout << "Density computation and output written to file successfully." << std::endl;
        //std::cout << atmosphereModel.getTotalDensity( 1000.0E3, 0.0, 0.0, 1.0697333, 0.0, 16 ,12, 2000) << std::endl;

//    return 0;
        //std::cout << atmosphereModel.getTotalDensity( 138.0E3, 0.0, 0.0, 1.0697333, 0.0, 16 ,12, 2000) <<std::endl;

//
//namespace tudat
//{
//namespace unit_tests
//{
//
//BOOST_AUTO_TEST_SUITE( test_mars_dtm_atmosphere )
//
//BOOST_AUTO_TEST_CASE( testMarsDtmAtmosphere )
//{
//
//}
//
//BOOST_AUTO_TEST_SUITE_END( )
//
//} // namespace unit_tests
//} // namespace tudat
}
