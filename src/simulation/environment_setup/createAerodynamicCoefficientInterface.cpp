/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <boost/lambda/lambda.hpp>

#include "tudat/simulation/environment_setup/createAerodynamicCoefficientInterface.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to create aerodynamic coefficient settings from coefficients stored in data files
std::shared_ptr< AerodynamicCoefficientSettings > readTabulatedAerodynamicCoefficientsFromFiles(
        const std::map< int, std::string > forceCoefficientFiles,
        const std::map< int, std::string > momentCoefficientFiles,
        const double referenceLength,
        const double referenceArea,
        const std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > independentVariableNames,
        const aerodynamics::AerodynamicCoefficientFrames forceCoefficientFrame,
        const aerodynamics::AerodynamicCoefficientFrames momentCoefficientFrame,
        const Eigen::Vector3d& momentReferencePoint,
        const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings )
{
    // Retrieve number of independent variables from file.
    int numberOfIndependentVariables =
            input_output::getNumberOfIndependentVariablesInCoefficientFile( forceCoefficientFiles.begin( )->second );

    // Call approriate file reading function for N independent variables
    std::shared_ptr< AerodynamicCoefficientSettings > coefficientSettings;
    if( numberOfIndependentVariables == 1 )
    {
        coefficientSettings = readGivenSizeTabulatedAerodynamicCoefficientsFromFiles< 1 >(
                    forceCoefficientFiles, momentCoefficientFiles, referenceLength, referenceArea,
                    momentReferencePoint, independentVariableNames,
                    forceCoefficientFrame,
                    momentCoefficientFrame, interpolatorSettings );
    }
    else if( numberOfIndependentVariables == 2 )
    {
        coefficientSettings = readGivenSizeTabulatedAerodynamicCoefficientsFromFiles< 2 >(
                    forceCoefficientFiles, momentCoefficientFiles, referenceLength, referenceArea,
                    momentReferencePoint, independentVariableNames,
                    forceCoefficientFrame,
                    momentCoefficientFrame, interpolatorSettings );
    }
    else if( numberOfIndependentVariables == 3 )
    {
        coefficientSettings = readGivenSizeTabulatedAerodynamicCoefficientsFromFiles< 3 >(
                    forceCoefficientFiles, momentCoefficientFiles, referenceLength, referenceArea,
                    momentReferencePoint, independentVariableNames,
                    forceCoefficientFrame,
                    momentCoefficientFrame, interpolatorSettings );
    }
    else
    {
        throw std::runtime_error( "Error when reading aerodynamic coefficient settings from file, found " +
                                  std::to_string( numberOfIndependentVariables ) +
                                  " independent variables, up to 3 currently supported" );
    }

    if( !momentReferencePoint.hasNaN( ) )
    {
        coefficientSettings->setAddForceContributionToMoments( true );
    }
    return coefficientSettings;
}

std::shared_ptr< AerodynamicCoefficientSettings > readTabulatedAerodynamicCoefficientsFromFilesDeprecated(
        const std::map< int, std::string > forceCoefficientFiles,
        const std::map< int, std::string > momentCoefficientFiles,
        const double referenceLength,
        const double referenceArea,
        const double lateralReferenceLength,
        const Eigen::Vector3d& momentReferencePoint,
        const std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > independentVariableNames,
        const bool areCoefficientsInAerodynamicFrame,
        const bool areCoefficientsInNegativeAxisDirection,
        const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings )
{
    if( referenceLength != lateralReferenceLength && lateralReferenceLength == lateralReferenceLength )
    {
        throw std::runtime_error( "Error when using deprecated reading of aerodynamic coefficients from files, lateral reference length must be equal to regular reference length" );
    }
    return readTabulatedAerodynamicCoefficientsFromFiles(
            forceCoefficientFiles, momentCoefficientFiles, referenceLength, referenceArea,
            independentVariableNames,
            aerodynamics::getAerodynamicCoefficientFrame( areCoefficientsInAerodynamicFrame, areCoefficientsInNegativeAxisDirection ),
            aerodynamics::getAerodynamicCoefficientFrame( areCoefficientsInAerodynamicFrame, areCoefficientsInNegativeAxisDirection ),
            momentReferencePoint,
            interpolatorSettings );
}


//! Function to create aerodynamic coefficient settings from coefficients stored in data files
std::shared_ptr< AerodynamicCoefficientSettings >
readTabulatedAerodynamicCoefficientsFromFiles(
        const std::map< int, std::string > forceCoefficientFiles,
        const double referenceArea,
        const std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > independentVariableNames,
        const aerodynamics::AerodynamicCoefficientFrames forceCoefficientFrame,
        const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings )
{
    // Retrieve number of independent variables from file.
    int numberOfIndependentVariables =
            input_output::getNumberOfIndependentVariablesInCoefficientFile( forceCoefficientFiles.begin( )->second );

    // Call approriate file reading function for N independent variables
    std::shared_ptr< AerodynamicCoefficientSettings > coefficientSettings;
    if( numberOfIndependentVariables == 1 )
    {
        coefficientSettings = readGivenSizeTabulatedAerodynamicCoefficientsFromFiles< 1 >(
                forceCoefficientFiles, referenceArea, independentVariableNames,
                forceCoefficientFrame, interpolatorSettings );
    }
    else if( numberOfIndependentVariables == 2 )
    {
        coefficientSettings = readGivenSizeTabulatedAerodynamicCoefficientsFromFiles< 2 >(
                forceCoefficientFiles, referenceArea, independentVariableNames,
                forceCoefficientFrame, interpolatorSettings );
    }
    else if( numberOfIndependentVariables == 3 )
    {
        coefficientSettings = readGivenSizeTabulatedAerodynamicCoefficientsFromFiles< 3 >(
                forceCoefficientFiles, referenceArea, independentVariableNames,
                forceCoefficientFrame, interpolatorSettings );
    }
    else
    {
        throw std::runtime_error( "Error when reading aerodynamic coefficient settings from file, found " +
                                  std::to_string( numberOfIndependentVariables ) +
                                  " independent variables, up to 3 currently supported" );
    }
    return coefficientSettings;
}

//! Function to create aerodynamic coefficient settings from coefficients stored in data files
std::shared_ptr< AerodynamicCoefficientSettings >
readTabulatedAerodynamicCoefficientsFromFilesDeprecated(
        const std::map< int, std::string > forceCoefficientFiles,
        const double referenceArea,
        const std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > independentVariableNames,
        const bool areCoefficientsInAerodynamicFrame,
        const bool areCoefficientsInNegativeAxisDirection,
        const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings )
{
    return readTabulatedAerodynamicCoefficientsFromFiles(
            forceCoefficientFiles,  referenceArea, independentVariableNames,  aerodynamics::getAerodynamicCoefficientFrame(
                    areCoefficientsInAerodynamicFrame, areCoefficientsInNegativeAxisDirection ), interpolatorSettings );
}

//! Function to create an aerodynamic coefficient interface containing constant coefficients.
std::shared_ptr< aerodynamics::AerodynamicCoefficientInterface >
createConstantCoefficientAerodynamicCoefficientInterface(
        const Eigen::Vector3d constantForceCoefficient,
        const Eigen::Vector3d constantMomentCoefficient,
        const double referenceLength,
        const double referenceArea,
        const Eigen::Vector3d& momentReferencePoint,
        const aerodynamics::AerodynamicCoefficientFrames forceCoefficientsFrame,
        const aerodynamics::AerodynamicCoefficientFrames momentCoefficientsFrame  )
{
    // Create coefficient interface
    std::shared_ptr< aerodynamics::AerodynamicCoefficientInterface > coefficientInterface =
            std::make_shared< aerodynamics::CustomAerodynamicCoefficientInterface >(
                [ = ]( const std::vector< double >& ){ return constantForceCoefficient; },
                [ = ]( const std::vector< double >& ){ return constantMomentCoefficient; },
                referenceLength, referenceArea, momentReferencePoint,
                std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables >( ),
                forceCoefficientsFrame, momentCoefficientsFrame );
    coefficientInterface->updateFullCurrentCoefficients( std::vector< double >( ) );
    return coefficientInterface;
}

//! Function to create an aerodynamic coefficient interface containing constant coefficients.
std::shared_ptr< aerodynamics::AerodynamicCoefficientInterface >
createZeroParameterAerodynamicCoefficientInterface(
        const std::function< Eigen::Vector3d( ) > constantForceCoefficientFunction,
        const std::function< Eigen::Vector3d( ) > constantMomentCoefficientFunction,
        const double referenceLength,
        const double referenceArea,
        const Eigen::Vector3d& momentReferencePoint,
        const aerodynamics::AerodynamicCoefficientFrames forceCoefficientsFrame,
        const aerodynamics::AerodynamicCoefficientFrames momentCoefficientsFrame )
{
    // Create coefficient interface
    std::shared_ptr< aerodynamics::AerodynamicCoefficientInterface > coefficientInterface =
            std::make_shared< aerodynamics::CustomAerodynamicCoefficientInterface >(
                [ = ]( const std::vector< double >& ){ return constantForceCoefficientFunction( ); },
                [ = ]( const std::vector< double >& ){ return constantMomentCoefficientFunction( ); },
                referenceLength, referenceArea, momentReferencePoint,
                std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables >( ),
                forceCoefficientsFrame, momentCoefficientsFrame );
    coefficientInterface->updateFullCurrentCoefficients( std::vector< double >( ) );
    return coefficientInterface;
}

//! Factory function for tabulated (1-D independent variables) aerodynamic coefficient interface from coefficient settings.
std::shared_ptr< aerodynamics::AerodynamicCoefficientInterface >
createUnivariateTabulatedCoefficientAerodynamicCoefficientInterface(
        const std::shared_ptr< AerodynamicCoefficientSettings > coefficientSettings,
        const std::string& body )
{
    using namespace tudat::interpolators;

    // Check consistency of type.
    std::shared_ptr< TabulatedAerodynamicCoefficientSettings< 1 > > tabulatedCoefficientSettings =
            std::dynamic_pointer_cast< TabulatedAerodynamicCoefficientSettings< 1 > >( coefficientSettings );
    if( tabulatedCoefficientSettings == nullptr )
    {
        throw std::runtime_error(
                    "Error, expected tabulated aerodynamic coefficients of size " +
                    std::to_string( 1 ) + "for body " + body );
    }
    else
    {

        // Retrieve or generate interpolation settings
        std::shared_ptr< OneDimensionalInterpolator< double, Eigen::Vector3d > > forceInterpolator;
        std::shared_ptr< OneDimensionalInterpolator< double, Eigen::Vector3d > > momentInterpolator;
        if ( tabulatedCoefficientSettings->getInterpolatorSettings( ) == nullptr )
        {
            forceInterpolator = createOneDimensionalInterpolator( tabulatedCoefficientSettings->getForceCoefficients( ),
                                                                  std::make_shared< InterpolatorSettings >( linear_interpolator ) );
            momentInterpolator = createOneDimensionalInterpolator( tabulatedCoefficientSettings->getForceCoefficients( ),
                                                                   std::make_shared< InterpolatorSettings >( linear_interpolator ) );
        }
        else
        {
            forceInterpolator = createOneDimensionalInterpolator( tabulatedCoefficientSettings->getForceCoefficients( ),
                                                                  std::dynamic_pointer_cast< InterpolatorSettings >(
                                                                      tabulatedCoefficientSettings->getInterpolatorSettings( ) ) );
            momentInterpolator = createOneDimensionalInterpolator( tabulatedCoefficientSettings->getForceCoefficients( ),
                                                                   std::dynamic_pointer_cast< InterpolatorSettings >(
                                                                       tabulatedCoefficientSettings->getInterpolatorSettings( ) ) );
        }

        // Create aerodynamic coefficient interface.
        return  std::make_shared< aerodynamics::CustomAerodynamicCoefficientInterface >(
                    std::bind( &Interpolator< double, Eigen::Vector3d >::interpolate, forceInterpolator, std::placeholders::_1 ),
                    std::bind( &Interpolator< double, Eigen::Vector3d >::interpolate, momentInterpolator, std::placeholders::_1 ),
                    tabulatedCoefficientSettings->getReferenceLength( ),
                    tabulatedCoefficientSettings->getReferenceArea( ),
                    tabulatedCoefficientSettings->getMomentReferencePoint( ),
                    tabulatedCoefficientSettings->getIndependentVariableNames( ),
                    tabulatedCoefficientSettings->getForceCoefficientsFrame( ),
                    tabulatedCoefficientSettings->getMomentCoefficientsFrame( ) );
    }
}


// std::shared_ptr< aerodynamics::AerodynamicCoefficientInterface >
// createBridgedAerodynamicCoefficientInterface(
//     const std::shared_ptr< AerodynamicCoefficientSettings > coefficientSettings,
//     const std::string& body,
//     const SystemOfBodies& bodies
//     )
// {

//     using namespace tudat::aerodynamics;

//     // Check consistency of type.
//     std::shared_ptr< BridgedModelsAerodynamicCoefficientSettings > bridgedCoefficientSettings =
//             std::dynamic_pointer_cast< BridgedModelsAerodynamicCoefficientSettings >(
//                 coefficientSettings );
//     if( bridgedCoefficientSettings == nullptr )
//     {
//         throw std::runtime_error(
//                     "Error, expected bridged aerodynamic coefficients for body " + body );
//     }
//     else
//     {

//         // Create coefficient interface.
//         // First, create both coefficient interfaces to be bridged

//         std::shared_ptr< AerodynamicCoefficientInterface > coefficientInterface_1 =
//             createAerodynamicCoefficientInterface(
//                 bridgedCoefficientSettings->getCoefficientSettings1( ), body, bodies );

//         std::shared_ptr< AerodynamicCoefficientInterface > coefficientInterface_2 =
//             createAerodynamicCoefficientInterface(
//                 bridgedCoefficientSettings->getCoefficientSettings2( ), body, bodies );

//         // Perform checks on the consistency of the settings
//         // Check if the models use the same reference area
//         if ( coefficientInterface_1->getReferenceArea( ) != coefficientInterface_2->getReferenceArea( ) )
//         {
//             throw std::runtime_error(
//                 "Error, expected the same reference area for both models for bridged aerodynamic coefficients for body " + body );
//         }

//         // Check if the models use the same reference length
//         if ( coefficientInterface_1->getReferenceLength( ) != coefficientInterface_2->getReferenceLength( ) )
//         {
//             throw std::runtime_error(
//                 "Error, expected the same reference length for both models for bridged aerodynamic coefficients for body " + body );
//         }

//         // Check if the models use the same moment reference point
//         if ( coefficientInterface_1->getMomentReferencePoint( ) != coefficientInterface_2->getMomentReferencePoint( ) )
//         {
//             throw std::runtime_error(
//                 "Error, expected the same moment reference point for both models for bridged aerodynamic coefficients for body " + body );
//         }

//         // Check if the models use the same force coefficients frame
//         if ( coefficientInterface_1->getForceCoefficientsFrame( ) != coefficientInterface_2->getForceCoefficientsFrame( ) )
//         {
//             throw std::runtime_error(
//                 "Error, expected the same force coefficients frame for both models for bridged aerodynamic coefficients for body " + body );
//         }

//         // Check if the models use the same moment coefficients frame
//         if ( coefficientInterface_1->getMomentCoefficientsFrame( ) != coefficientInterface_2->getMomentCoefficientsFrame( ) )
//         {
//             throw std::runtime_error(
//                 "Error, expected the same moment coefficients frame for both models for bridged aerodynamic coefficients for body " + body );
//         }

//         int numberOfBridgingVars = bridgedCoefficientSettings->getBridgingVariableNames( ).size( );
//         if ( numberOfBridgingVars != 1 )
//         {
//             throw std::runtime_error(
//                 "Error, expected one bridging variable for bridged aerodynamic coefficients for body " + body );
//         }

//         int numberOfModel1Vars = coefficientInterface_1->getIndependentVariableNames( ).size( );
//         int numberOfModel2Vars = coefficientInterface_2->getIndependentVariableNames( ).size( );

//         // Create force coefficient function
//         std::function< Eigen::Vector6d( const std::vector< double >& ) > aerodynamicCoefficientsFunction =
//             [ = ]( const std::vector< double >& independentVariables )
//             {
//                 // check if we are in the first, second or bridged region 

//                 if ( independentVariables[ 0 ] < bridgedCoefficientSettings->getBridgingFunctionLimits( ).first )
//                 {
//                     // update the coefficients
//                     coefficientInterface_1->updateCurrentCoefficients( std::vector<double>(independentVariables.begin() + 1, independentVariables.begin() + numberOfModel1Vars + 1) );

//                     return coefficientInterface_1->getCurrentAerodynamicCoefficients( );
//                 }
//                 else if ( independentVariables[ 0 ] > bridgedCoefficientSettings->getBridgingFunctionLimits( ).second )
//                 {
//                     // update the coefficients
//                     coefficientInterface_2->updateCurrentCoefficients( std::vector<double>(independentVariables.begin() + numberOfModel1Vars + 1, independentVariables.begin() + numberOfModel1Vars + numberOfModel2Vars + 1) );

//                     return coefficientInterface_2->getCurrentAerodynamicCoefficients( );
//                 }
//                 else
//                 {
//                     // bridge between the two models
//                     // update both models
//                     coefficientInterface_1->updateCurrentCoefficients( std::vector<double>(independentVariables.begin() + 1, independentVariables.begin() + numberOfModel1Vars + 1) );
//                     coefficientInterface_2->updateCurrentCoefficients( std::vector<double>(independentVariables.begin() + numberOfModel1Vars + 1, independentVariables.begin() + numberOfModel1Vars + numberOfModel2Vars + 1) );

//                     // get the bridging value
//                     double bridgingValue = bridgedCoefficientSettings->getBridgingFunction( )( independentVariables[ 0 ] );

//                     // get the force coefficients from both models
//                     Eigen::Vector6d aerodynamicCoefficients_1 = coefficientInterface_1->getCurrentAerodynamicCoefficients( );
//                     Eigen::Vector6d aerodynamicCoefficients_2 = coefficientInterface_2->getCurrentAerodynamicCoefficients( );

//                     // interpolate between the two
//                     Eigen::Vector6d bridgedCoeffs = (aerodynamicCoefficients_1.array() * bridgingValue + aerodynamicCoefficients_2.array() * (1.0 - bridgingValue)).matrix();
//                     return bridgedCoeffs;
//                 }
//             };
               
//         std::function< Eigen::Vector3d( const std::vector< double >& ) > forceCoefficientFunction =
//             [ = ]( const std::vector< double >& independentVariables )
//             {
//                 return aerodynamicCoefficientsFunction( independentVariables ).segment( 0, 3 );
//             };

//         std::function< Eigen::Vector3d( const std::vector< double >& ) > momentCoefficientFunction =
//             [ = ]( const std::vector< double >& independentVariables )
//             {
//                 return aerodynamicCoefficientsFunction( independentVariables ).segment( 3, 3 );
//             };


//         std::vector< AerodynamicCoefficientsIndependentVariables > vec1 = bridgedCoefficientSettings->getBridgingVariableNames( );
//         std::vector< AerodynamicCoefficientsIndependentVariables > vec2 = coefficientInterface_1->getIndependentVariableNames( );
//         std::vector< AerodynamicCoefficientsIndependentVariables > vec3 = coefficientInterface_2->getIndependentVariableNames( );

//         std::vector< AerodynamicCoefficientsIndependentVariables > allIndependentVariableNames;
//         // Combine vec1 into allIndependentVariableNames
//         allIndependentVariableNames.insert(allIndependentVariableNames.end(), vec1.begin(), vec1.end());

//         // Combine vec2 into allIndependentVariableNames
//         allIndependentVariableNames.insert(allIndependentVariableNames.end(), vec2.begin(), vec2.end());

//         // Combine vec3 into allIndependentVariableNames
//         allIndependentVariableNames.insert(allIndependentVariableNames.end(), vec3.begin(), vec3.end());
        

//         // Create aerodynamic coefficient interface.
//         return  std::make_shared< CustomAerodynamicCoefficientInterface >(
//                 forceCoefficientFunction, 
//                 momentCoefficientFunction,
//                 coefficientInterface_1->getReferenceLength( ),
//                 coefficientInterface_1->getReferenceArea( ),
//                 coefficientInterface_1->getMomentReferencePoint( ),
//                 allIndependentVariableNames, // independent variables
//                 coefficientInterface_1->getForceCoefficientsFrame( ),
//                 coefficientInterface_1->getMomentCoefficientsFrame( ) );
//     }
// }

std::shared_ptr< aerodynamics::AerodynamicMomentContributionInterface > createMomentContributionInterface(
            const aerodynamics::AerodynamicCoefficientFrames forceCoefficientFrame,
            const aerodynamics::AerodynamicCoefficientFrames momentCoefficientFrame,
            const std::shared_ptr< Body > body )
{
    if( body->getFlightConditions( ) == nullptr )
    {
        throw std::runtime_error( "Error when creating moment contribution interface, body " + body->getBodyName( ) + " has no flight conditions." );
    }

    std::pair< reference_frames::AerodynamicsReferenceFrames, double > forceCoefficientFrameId =
            convertCoefficientFrameToGeneralAerodynamicFrame( forceCoefficientFrame );
    std::pair< reference_frames::AerodynamicsReferenceFrames, double > momentCoefficientFrameId =
            convertCoefficientFrameToGeneralAerodynamicFrame( momentCoefficientFrame );
    std::function< Eigen::Matrix3d( ) > coefficientRotationFunction;
    std::function< Eigen::Matrix3d( ) > armRotationFunction;

    if( forceCoefficientFrameId.first == momentCoefficientFrameId.first )
    {
        coefficientRotationFunction = [=](){ return Eigen::Matrix3d::Identity( ); };
    }
    else
    {
        coefficientRotationFunction = [=](){ return body->getFlightConditions( )->getAerodynamicAngleCalculator( )->getRotationMatrixBetweenFrames(
                forceCoefficientFrameId.first, momentCoefficientFrameId.first ); };
    }

    if( momentCoefficientFrameId.first == reference_frames::body_frame )
    {
        armRotationFunction = [=](){ return Eigen::Matrix3d::Identity( ); };
    }
    else
    {
        armRotationFunction = [=]( )
            {
                return body->getFlightConditions( )->getAerodynamicAngleCalculator( )->getRotationMatrixBetweenFrames(
                    reference_frames::body_frame, momentCoefficientFrameId.first );
            };
    }

    return std::make_shared< aerodynamics::AerodynamicMomentContributionInterface >(
            coefficientRotationFunction, armRotationFunction, std::bind( &Body::getBodyFixedCenterOfMass, body ), forceCoefficientFrameId.second );
}

std::shared_ptr< aerodynamics::AerodynamicMomentContributionInterface > createMomentContributionInterface(
    const std::shared_ptr< AerodynamicCoefficientSettings > coefficientSettings,
    const std::shared_ptr< Body > body )
{
    return createMomentContributionInterface(
        coefficientSettings->getForceCoefficientsFrame( ),
        coefficientSettings->getMomentCoefficientsFrame( ),
        body );
}

//! Function to create and aerodynamic coefficient interface.
std::shared_ptr< aerodynamics::AerodynamicCoefficientInterface >
createAerodynamicCoefficientInterface(
        const std::shared_ptr< AerodynamicCoefficientSettings > coefficientSettings,
        const std::string& body,
        const SystemOfBodies& bodies )
{

    using namespace tudat::aerodynamics;

    std::shared_ptr< AerodynamicCoefficientInterface > coefficientInterface;

    // Check type of interface that is to be created.
    switch( coefficientSettings->getAerodynamicCoefficientType( ) )
    {
    case constant_aerodynamic_coefficients:
    {
        // Check consistency of type.
        std::shared_ptr< ConstantAerodynamicCoefficientSettings > constantCoefficientSettings =
                std::dynamic_pointer_cast< ConstantAerodynamicCoefficientSettings >(
                    coefficientSettings );
        if( constantCoefficientSettings == nullptr )
        {
            throw std::runtime_error(
                        "Error, expected constant aerodynamic coefficients for body " + body );
        }
        else
        {
            // create constant interface.
            coefficientInterface = createConstantCoefficientAerodynamicCoefficientInterface(
                        constantCoefficientSettings->getConstantForceCoefficient( ),
                        constantCoefficientSettings->getConstantMomentCoefficient( ),
                        constantCoefficientSettings->getReferenceLength( ),
                        constantCoefficientSettings->getReferenceArea( ),
                        constantCoefficientSettings->getMomentReferencePoint( ),
                        constantCoefficientSettings->getForceCoefficientsFrame( ),
                        constantCoefficientSettings->getMomentCoefficientsFrame( ) );
        }
        break;
    }
    case custom_aerodynamic_coefficients:
    {
        // Check consistency of type.
        std::shared_ptr< CustomAerodynamicCoefficientSettings > customCoefficientSettings =
                std::dynamic_pointer_cast< CustomAerodynamicCoefficientSettings >(
                    coefficientSettings );
        if( customCoefficientSettings == nullptr )
        {
            throw std::runtime_error(
                        "Error, expected custom aerodynamic coefficients for body " + body );
        }
        else
        {
            // create constant interface.
            coefficientInterface = std::make_shared< CustomAerodynamicCoefficientInterface >(
                        customCoefficientSettings->getForceCoefficientFunction( ),
                        customCoefficientSettings->getMomentCoefficientFunction( ),
                        customCoefficientSettings->getReferenceLength( ),
                        customCoefficientSettings->getReferenceArea( ),
                        customCoefficientSettings->getMomentReferencePoint( ),
                        customCoefficientSettings->getIndependentVariableNames( ),
                        customCoefficientSettings->getForceCoefficientsFrame( ),
                        customCoefficientSettings->getMomentCoefficientsFrame( ) );
        }
        break;
    }
    case tabulated_coefficients:
    {
        // Check number of dimensions of tabulated coefficients.
        int numberOfDimensions = coefficientSettings->getIndependentVariableNames( ).size( );
        switch( numberOfDimensions )
        {
        case 1:
        {
            coefficientInterface = createUnivariateTabulatedCoefficientAerodynamicCoefficientInterface(
                        coefficientSettings, body );
            break;
        }
        case 2:
        {
            coefficientInterface = createTabulatedCoefficientAerodynamicCoefficientInterface< 2 >(
                        coefficientSettings, body );
            break;
        }
        case 3:
        {
            coefficientInterface = createTabulatedCoefficientAerodynamicCoefficientInterface< 3 >(
                        coefficientSettings, body );
            break;
        }
        case 4:
        {
            coefficientInterface = createTabulatedCoefficientAerodynamicCoefficientInterface< 4 >(
                        coefficientSettings, body );
            break;
        }
        case 5:
        {
            coefficientInterface = createTabulatedCoefficientAerodynamicCoefficientInterface< 5 >(
                        coefficientSettings, body );
            break;
        }
        case 6:
        {
            coefficientInterface = createTabulatedCoefficientAerodynamicCoefficientInterface< 6 >(
                        coefficientSettings, body );
            break;
        }
        default:
            throw std::runtime_error( "Error when making tabulated aerodynamic coefficient interface, " +
                                      std::to_string( numberOfDimensions ) + " dimensions not yet implemented" );
        }
        break;
    }
    case scaled_coefficients:
    {
        // Check consistency of type and class.
        std::shared_ptr< ScaledAerodynamicCoefficientInterfaceSettings > scaledCoefficientSettings =
                std::dynamic_pointer_cast< ScaledAerodynamicCoefficientInterfaceSettings >(
                    coefficientSettings );
        if( coefficientSettings->getAddForceContributionToMoments( ) )
        {
            throw std::runtime_error( "Error when creating scaled aerodynamic coefficients, force contribution to moments not permitted." );
        }
        if( scaledCoefficientSettings == nullptr )
        {
            throw std::runtime_error(
                        "Error, expected scaled aerodynamic coefficient settings for body " + body );
        }
        else
        {
            std::shared_ptr< AerodynamicCoefficientInterface > baseInterface = createAerodynamicCoefficientInterface(
                        scaledCoefficientSettings->getBaseSettings( ), body, bodies );
            coefficientInterface = std::make_shared< ScaledAerodynamicCoefficientInterface >(
                        baseInterface, scaledCoefficientSettings->getForceScaling( ),
                        scaledCoefficientSettings->getMomentScaling( ), scaledCoefficientSettings->getIsScalingAbsolute( ) );
        }
        break;
    }
    case rarefied_flow_aerodynamic_coefficients:
    {
        coefficientInterface = std::make_shared< RarefiedFlowAerodynamicCoefficientInterface >(
                    bodies.at( body )->getVehicleSystems( ), coefficientSettings->getReferenceLength( ),
                    coefficientSettings->getReferenceArea( ),
                    coefficientSettings->getMomentReferencePoint( ),
                    coefficientSettings->getIndependentVariableNames( ),
                    coefficientSettings->getForceCoefficientsFrame( ),
                    coefficientSettings->getMomentCoefficientsFrame( ));

        // print sucess to console
        std::cout<<"Rarefied flow aerodynamic coefficients created successfully"<<std::endl;

        
        break;
    }
    // case hypersonic_flow_aerodynamic_coefficients:
    // {
    //     coefficientInterface = std::make_shared< HypersonicFlowAerodynamicCoefficientInterface >(
    //                 bodies.at( body )->getVehicleSystems( ), coefficientSettings->getReferenceLength( ),
    //                 coefficientSettings->getReferenceArea( ),
    //                 coefficientSettings->getMomentReferencePoint( ),
    //                 coefficientSettings->getIndependentVariableNames( ),
    //                 coefficientSettings->getForceCoefficientsFrame( ),
    //                 coefficientSettings->getMomentCoefficientsFrame( ) );
    //     break;
    // }
    // case bridged_models_aerodynamic_coefficients:
    // {

    //     coefficientInterface = createBridgedAerodynamicCoefficientInterface(
    //                 coefficientSettings, body, bodies );
    //     break;
    // }    
    default:
        throw std::runtime_error( "Error, do not recognize aerodynamic coefficient settings for " + body );
    }

    // Create and set control surfaces
    if( coefficientSettings->getControlSurfaceSettings( ).size( ) != 0 )
    {
        std::map< std::string, std::shared_ptr< ControlSurfaceIncrementAerodynamicInterface > >
                controlSurfaceIncrementInterfaces;
        std::map< std::string, std::shared_ptr< ControlSurfaceIncrementAerodynamicCoefficientSettings > >
                controlSurfaceSettings = coefficientSettings->getControlSurfaceSettings( );
        for( std::map< std::string, std::shared_ptr< ControlSurfaceIncrementAerodynamicCoefficientSettings > >::iterator
             settingIterator = controlSurfaceSettings.begin( ); settingIterator != controlSurfaceSettings.end( );
             settingIterator++ )
        {
            controlSurfaceIncrementInterfaces[ settingIterator->first ] =
                    createControlSurfaceIncrementAerodynamicCoefficientInterface(
                        settingIterator->second, body );
        }
        coefficientInterface->setControlSurfaceIncrements( controlSurfaceIncrementInterfaces );

    }

    std::shared_ptr< AerodynamicMomentContributionInterface > momentContributionInterface;
    if( coefficientSettings->getAddForceContributionToMoments( ) )
    {
        if( bodies.count( body ) == 0 )
        {
            throw std::runtime_error( "Error when making aerodynamic moment correction interface, no body " + body + " was found." );
        }
        momentContributionInterface = createMomentContributionInterface( coefficientSettings, bodies.at( body ) );
        coefficientInterface->setMomentContributionInterface( momentContributionInterface );
    }

    return coefficientInterface;
}

std::shared_ptr< aerodynamics::AerodynamicCoefficientInterface >
createAerodynamicCoefficientInterfaceDeprecated(
        const std::shared_ptr< AerodynamicCoefficientSettings > coefficientSettings,
        const std::string& body )
{
    std::cerr<<"Warning, you are using an outdated version of the create_aerodynamic_coefficient_interface function. Please use the updated one, which has the bodies (as SystemOfBodies object) as third input argument"<<std::endl;
    SystemOfBodies bodies = SystemOfBodies( );
    return createAerodynamicCoefficientInterface(
            coefficientSettings, body, bodies );
}
} // simulation_setup

} // tudat
