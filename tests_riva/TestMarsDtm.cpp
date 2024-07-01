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
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/astro/aerodynamics/marsDtmAtmosphereModel.h"
#include "tudat/io/solarActivityData.h"
#include "tudat/astro/basic_astro/unitConversions.h"


int main( ) {
    using namespace tudat;
    using namespace aerodynamics;
    using namespace simulation_setup;
    using namespace numerical_integrators;
    using namespace simulation_setup;
    using namespace basic_astrodynamics;
    using namespace propagators;
    using namespace basic_mathematics;
    using namespace basic_astrodynamics;


    //std::string filename = "/Users/ralkahal/Documents/PhD/atmodensitypds/dtm_mars";
    std::string filename = "/Users/ralkahal/OneDrive - Delft University of Technology/PhD/Programs/atmodensitydtm/dtm_mars";

    //marsDtmAtmosphereSettings =
    //std::shared_ptr< AtmosphereSettings > atmosphereSettings;
    std::shared_ptr< AtmosphereSettings > marsDtmAtmosphereSettings;
    //std::shared_ptr< AtmosphereModel > atmosphereModel =
    //       std::dynamic_pointer_cast< MarsDtmAtmosphereModel >(
    //              createAtmosphereModel( marsDtmAtmosphereSettings, "Mars" ) );
    std::string spaceWeatherFilePath = paths::getSpaceWeatherDataPath( ) + "/sw19571001.txt";
    tudat::input_output::solar_activity::SolarActivityDataMap solarActivityData =
            tudat::input_output::solar_activity::readSolarActivityData( spaceWeatherFilePath );
    std::shared_ptr< input_output::solar_activity::SolarActivityContainer > solarActivityContainer =
            std::make_shared< input_output::solar_activity::SolarActivityContainer >( solarActivityData );
    std::function< double( const double ) > f107Function = [=](const double time)
    {
        return solarActivityContainer->getSolarActivityData( time )->solarRadioFlux107Observed / 2.25;
    };

    marsDtmAtmosphereSettings = std::make_shared< MarsDtmAtmosphereSettings >(
            filename, 3378.0E3);

    std::shared_ptr< MarsDtmAtmosphereModel > atmosphereModel =
            std::dynamic_pointer_cast< MarsDtmAtmosphereModel >(
                    createAtmosphereModel( marsDtmAtmosphereSettings, "Mars" ) );

    /*if( atmosphereModel == nullptr )
    {
        std::cerr<< "Atmosphere model is null" << std::endl;
    }
     */
    //MarsDtmAtmosphereModel(3376.78E3, filename);
    //atmosphereModel->computelocalsolartime(0.0,16 ,12,2000, 0.0, 1.0697333);
    //std::cout<<atmosphereModel->computeGl(0.0,0.0,1.0697333,0.0,16,12,2000,1)<<std::endl;
    //std::cout<<std::abs((tudat::basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch(
    //       2021, 2, 7, 0.0, 0.0, 0.0, tudat::basic_astrodynamics::JULIAN_DAY_ON_J2000) -
    //tudat::basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch(2022,2,18,0.0,0.0,0.0,tudat::basic_astrodynamics::JULIAN_DAY_ON_J2000))*24/24.63)<<std::endl;

    //marsDate date2 = marsDate(2022, 2, 18, 0.0, 0.0, 0.0);
    //std::cout << date2.marsDayofYear(date2) << std::endl;
    //std::cout<<"Geopotential height "<< atmosphereModel->computeGeopotentialAltitude( 255.0E3 )<<std::endl;
    //std::cout<<atmosphereModel->computeCurrentTemperature( 0.0, 0.0, 1.0697333, 0.0, 16 ,12, 2000, 1)<<std::endl;
    //std::cout<<atmosphereModel->computeGamma( 0.0, 0.0, 1.0697333, 0.0, 16 ,12, 2000, 1)<<std::endl;
    //std::cout<< atmosphereModel->heightDistributionFunction(255.0E3, 0.0, 0.0, 1.0697333, 0.0, 16 ,12, 2000, 1)<<std::endl;
  /* std::ofstream outputFile(
           "/Users/ralkahal/OneDrive - Delft University of Technology/PhD/Programs/atmodensitydtm/density_output_pure.txt");
    // Check if the file is opened successfully

    if (!outputFile.is_open()) {
        std::cerr << "Unable to open the file!" << std::endl;
        return 1; // return an error code
    }
*/
 /*
    for (int altitude = 138E3; altitude <= 1000E3; altitude += 10E3) {
        double alt_km = static_cast<double>(altitude);

        double rho = atmosphereModel->getTotalDensity(alt_km, 0.0, 0.0, 1.0697333, 0.0, 16, 12, 2000);

        // Write altitude and corresponding density to the file
        outputFile << alt_km << " " << rho << "\n";
    }
  */
    //for (int altitude = 138E3; altitude <= 1000E3; altitude += 10E3) {
    //  double alt_km = static_cast<double>(altitude);
    int alt_km = 400E3;
    //668 * 24 * 60 * 60
    for (int time = 0.0; time <= 345600 ; time += 86400) {
        //for (int time = 0; time <= 365 * 24 * 60 * 60; time += 60) {

        basic_astrodynamics::DateTime currentDateTime_ = basic_astrodynamics::getCalendarDateFromTime(time);
        int minutes = currentDateTime_.getMinute();
        int hours = currentDateTime_.getHour();
        int days = currentDateTime_.getDay();
        int months = currentDateTime_.getMonth();
        int years = currentDateTime_.getYear();

        std::string filename =
                "/Users/ralkahal/OneDrive - Delft University of Technology/PhD/atmodensitypds/densityFiles/" +
                std::to_string(time) + ".txt";
        std::ofstream outfile(filename);

        //std::cout << "Time: " << time << " minutes: " << minutes << " hours: " << hours << " days: " << days << " months: " << months << " years: " << years << std::endl;
        //double latitude = 0.0;
        //double longitude = 0.0;
        if (outfile.is_open()) {
            for (double latitude = -90.0; latitude <= 90.0; latitude += 1.0) {
                for (double longitude = 0.0; longitude <= 360.0; longitude += 1.0) {
                    double rho = atmosphereModel->getTotalDensity(alt_km, unit_conversions::convertDegreesToRadians(latitude), unit_conversions::convertDegreesToRadians(longitude), minutes, hours, days,
                                                                  months, years);
                    auto results = computeSolarLongitude(longitude, days, months, years, hours, minutes);
                    double Ls = std::get<0>(results);
                    double currentF107 = 0.0; // Assuming this is defined somewhere
                    outfile << latitude << " " << longitude << " " << rho << " " << Ls << "\n";
                }
            }

        } else {
            std::cerr << "Error creating file " << filename << std::endl;
        }
        outfile.close();




        //std::cout << "Density computation and output written to file successfully." << std::endl;
        //std::cout << atmosphereModel->getTotalDensity( 1000.0E3, 0.0, 0.0, 1.0697333, 0.0, 16 ,12, 2000) << std::endl;
    }
    return 0;
    //std::cout << atmosphereModel->getTotalDensity( 138.0E3, 0.0, 0.0, 1.0697333, 0.0, 16 ,12, 2000) <<std::endl;

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
// namespace tudat
//
}

