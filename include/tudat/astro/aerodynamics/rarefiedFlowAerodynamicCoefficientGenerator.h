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

#ifndef TUDAT_RAREFIEDFLOW_AERODYNAMIC_COEFFICIENT_GENERATOR_H
#define TUDAT_RAREFIEDFLOW_AERODYNAMIC_COEFFICIENT_GENERATOR_H

#include <map>
#include <string>
#include <vector>

#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <memory>

#include <Eigen/Core>

#include "tudat/astro/aerodynamics/aerodynamicCoefficientGenerator.h"
#include "tudat/basics/basicTypedefs.h"
#include "tudat/astro/system_models/vehicleExteriorPanels.h"

namespace tudat
{
namespace aerodynamics
{

template< unsigned int NumberOfIndependentVariables >
class RarefiedFlowAerodynamicCoefficientGenerator: public AerodynamicCoefficientGenerator< NumberOfIndependentVariables, 6 >
{
public:

    RarefiedFlowAerodynamicCoefficientGenerator(
            const std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > > vehiclePanels,
            const std::map< AerodynamicCoefficientsIndependentVariables, std::vector< double > >& dataPointsOfIndependentVariables,
            const double referenceArea,
            const double referenceLength,
            const Eigen::Vector3d& momentReferencePoint,
            const bool savePressureCoefficients = false ) : AerodynamicCoefficientGenerator< NumberOfIndependentVariables, 6 >(
        dataPointsOfIndependentVariables, referenceLength, referenceArea,
        momentReferencePoint, utilities::createVectorFromMapKeys( dataPointsOfIndependentVariables ),
        positive_aerodynamic_frame_coefficients, positive_aerodynamic_frame_coefficients ),
    vehiclePanels_( vehiclePanels ), dataPointsOfIndependentVariables_( utilities::createVectorFromMapValues( dataPointsOfIndependentVariables ) ),
    angleOfAttackIndex_( -1 ), sideslipAngleIndex_( -1 )
    {
        if( dataPointsOfIndependentVariables.count( angle_of_attack_dependent ) > 0 )
        {
            angleOfAttackIndex_ = std::distance( dataPointsOfIndependentVariables.begin( ),
                                     std::find(dataPointsOfIndependentVariables.begin(), dataPointsOfIndependentVariables.end(), angle_of_attack_dependent ) );
        }

        if( dataPointsOfIndependentVariables.count( angle_of_sideslip_dependent ) > 0 )
        {
            sideslipAngleIndex_ = std::distance( dataPointsOfIndependentVariables.begin( ),
                                                 std::find(dataPointsOfIndependentVariables.begin(), dataPointsOfIndependentVariables.end(), angle_of_attack_dependent ) );
        }

        // Allocate memory for panel inclinations and currentPressureCoefficients_.
        currentPanelInclinations_.resize( vehiclePanels_.size( ) );
        currentPressureCoefficients_.resize( vehiclePanels_.size( ) );
        
        for( int i = 0; i < NumberOfIndependentVariables; i++ )
        {
            numberOfPointsPerIndependentVariables_[ i ] =
                dataPointsOfIndependentVariables_[ i ].size( );
        }

        isCoefficientGenerated_.resize( numberOfPointsPerIndependentVariables_ );

        std::fill( isCoefficientGenerated_.origin( ),
                   isCoefficientGenerated_.origin( ) + isCoefficientGenerated_.num_elements( ), 0 );

        generateCoefficients( );
        this->createInterpolator( );
    }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~RarefiedFlowAerodynamicCoefficientGenerator( ) { }

    //! Get aerodynamic coefficients.
    /*!
     *  Returns aerodynamic coefficients.
     *  The physical meaning of each of the three independent variables is: 0 = mach number,
     *  1 = angle of attack, 2 = angle of sideslip.
     * \param independentVariables Array of values of independent variable
     *          indices in dataPointsOfIndependentVariables_.
     * \return vector of coefficients at specified independent variable indices.
     */
    Eigen::Vector6d getAerodynamicCoefficientsDataPoint(
            const boost::array< int, NumberOfIndependentVariables > independentVariables )
    {
        if( isCoefficientGenerated_( independentVariables ) == 0 )
        {
            determineVehicleCoefficients( independentVariables );
        }

        // Return requested coefficients.
        return aerodynamicCoefficients_( independentVariables );
    }

    //! Determine inclination angles of panels on a given part.
    /*!
     * Determines panel inclinations for all panels on all parts for given attitude.
     * Outward pointing surface-normals are assumed!
     * \param angleOfAttack Angle of attack at which to determine inclination angles.
     * \param angleOfSideslip Angle of sideslip at which to determine inclination angles.
     */
    void determineInclinations( const double angleOfAttack,
                                const double angleOfSideslip )

    {
        // Declare free-stream velocity vector.
        Eigen::Vector3d freestreamVelocityDirection;

        // Set freestream velocity vector in body frame.
        double freestreamVelocityDirectionX = std::cos( angleOfAttack )* std::cos( angleOfSideslip );
        double freestreamVelocityDirectionY = std::sin( angleOfSideslip );
        double freestreamVelocityDirectionZ = std::sin( angleOfAttack ) * std::cos( angleOfSideslip );
        freestreamVelocityDirection( 0 ) = freestreamVelocityDirectionX;
        freestreamVelocityDirection( 1 ) = freestreamVelocityDirectionY;
        freestreamVelocityDirection( 2 ) = freestreamVelocityDirectionZ;

        // Declare cosine of inclination angle.
        double cosineOfInclination;

        // Loop over all panels of given vehicle part and set inclination angles.
        for( unsigned int k = 0; k < vehiclePanels_.size( ); k++ )
        {
            cosineOfInclination = vehiclePanels_[ k ]->getFrameFixedSurfaceNormal( )( ).dot( freestreamVelocityDirection );
            // Set inclination angle.
            currentPanelInclinations_[ k ] = mathematical_constants::PI / 2.0 - std::acos( cosineOfInclination );
   
        }
    }
    
    std::vector< std::vector< std::vector< double > > > getPressureCoefficientList(
            const boost::array< int, 3 > independentVariables )
    {
        return pressureCoefficientList_.at( independentVariables );
    }

    void clearData( )
    {
        currentPanelInclinations_.clear( );
        currentPressureCoefficients_.clear( );
        previouslyComputedInclinations_.clear( );
        this->clearBaseData( );
    }


private:

    void generateCoefficients( )
    {
        // Allocate variable to pass to coefficient determination for independent
        // variable indices.
        boost::array<int, NumberOfIndependentVariables> independentVariableIndices;
        for( int i = 0; i < NumberOfIndependentVariables; i++ )
        {
            independentVariableIndices[ i ] = 0;
        }

        iterateOverAllIndependentVariables( independentVariableIndices );
    }

    void iterateOverAllIndependentVariables( boost::array< int, NumberOfIndependentVariables >& independentVariableIndices)
    {
        TODO write recursive for loop over variable dimension
        {
            determineVehicleCoefficients( independentVariableIndices );
        }
    }

    //! Generate aerodynamic coefficients at a single set of independent variables.
    /*!
     * Generates aerodynamic coefficients at a single set of independent variables.
     * Determines values and sets corresponding entry in vehicleCoefficients_ array.
     * \param independentVariableIndices Array of indices from lists of Mach number,
     *          angle of attack and angle of sideslip points at which to perform analysis.
     */
    void determineVehicleCoefficients( const boost::array< int, NumberOfIndependentVariables >& independentVariableIndices )
    {
        // Declare coefficients vector and initialize to zeros.
        Eigen::Vector6d coefficients = Eigen::Vector6d::Zero( );

        // Loop over all vehicle parts, calculate aerodynamic coefficients and add
        // to aerodynamicCoefficients_.
        double angleOfAttack = 0.0;
        if( angleOfAttackIndex_ >= 0 )
        {
            angleOfAttack =
                dataPointsOfIndependentVariables_[ angleOfAttackIndex_ ][ independentVariableIndices[ angleOfAttackIndex_ ]];
        }

        double angleOfSideslip = 0.0;
        if( sideslipAngleIndex_ >= 0 )
        {
            angleOfSideslip =
                dataPointsOfIndependentVariables_[ sideslipAngleIndex_ ][ independentVariableIndices[ sideslipAngleIndex_ ]];
        }

        // Check whether the inclinations of the vehicle part have already been computed.
        if ( previouslyComputedInclinations_.count( std::make_pair( angleOfAttack, angleOfSideslip ) ) == 0 )
        {
            // Determine panel inclinations for part.
            determineInclinations( angleOfAttack, angleOfSideslip );

            // Add panel inclinations to container
            previouslyComputedInclinations_[ std::pair< double, double >(
                angleOfAttack, angleOfSideslip ) ] = currentPanelInclinations_;
        }

        else
        {
            // Fetch inclinations from container
            currentPanelInclinations_ = previouslyComputedInclinations_[ std::make_pair(
                angleOfAttack, angleOfSideslip ) ];
        }

        // Set currentPressureCoefficients_ array for given independent variables.
        determinePressureCoefficients( independentVariableIndices );

        // Calculate force coefficients from pressure coefficients.
        coefficients.segment( 0, 3 ) = calculateForceCoefficients( );

        // Calculate moment coefficients from pressure coefficients.
        coefficients.segment( 3, 3 ) = calculateMomentCoefficients( );
        

        if( savePressureCoefficients_ )
        {
            pressureCoefficientList_[ independentVariableIndices ] = currentPressureCoefficients_;
        }

        this->aerodynamicCoefficients_( independentVariableIndices ) = coefficients;
        isCoefficientGenerated_( independentVariableIndices ) = true;
    }
    
    //! Determine pressure coefficients on a given part.
    /*!
     * Determines pressure coefficients on a single vehicle part.
     * Calls the updateExpansionPressures and updateCompressionPressures for given vehicle part.
     * \param partNumber Index from vehicleParts_ array for which to determine coefficients.
     * \param independentVariableIndices Array of indices of independent variables.
     */
    void determinePressureCoefficients( const boost::array< int, NumberOfIndependentVariables >& independentVariableIndices )
    {
        for( unsigned int i = 0; i < vehiclePanels_.size( ); i++ )
        {
            double currentPanelInclination = currentPanelInclinations_.at( i );
            currentPressureCoefficients_[ i ] = vehiclePanels_.at( i )->computePressureCoefficient( currentPanelInclination, TODO add indendepnt variables );
        }
    }

    //! Determine force coefficients of a part.
    /*!
     * Sums the pressure coefficients of given part and determines force coefficients from it by
     * non-dimensionalization with reference area.
     * \param partNumber Index from vehicleParts_ array for which determine coefficients.
     * \return Force coefficients for requested vehicle part.
     */
    Eigen::Vector3d calculateForceCoefficients( const int partNumber )
    {
        // Declare force coefficient vector and intialize to zeros.
        Eigen::Vector3d forceCoefficients = Eigen::Vector3d::Zero( );

        // Loop over all panels and add pressures, scaled by panel area, to force
        // coefficients.
        for ( int i = 0 ; i < vehiclePanels_.size( ) ; i++ )
        {
            forceCoefficients -=
                currentPressureCoefficients_.at( i ) *
                vehiclePanels_.at( i )->getPanelArea( ) *
                vehiclePanels_.at( i )->getFrameFixedSurfaceNormal( )( );
        }

        // Normalize result by reference area.
        forceCoefficients /= this->referenceArea_;

        return forceCoefficients;
    }


    //! Determine moment coefficients of a part.
    /*!
     * Determines the moment coefficients of a given part by summing the contributions of all
     * panels on the part. Moment arms are taken from panel centroid to momentReferencePoint. Non-
     * dimensionalization is performed by product of referenceLength and referenceArea.
     * \param partNumber Index from vehicleParts_ array for which to determine coefficients.
     * \return Moment coefficients for requested vehicle part.
     */
    Eigen::Vector3d calculateMomentCoefficients( const int partNumber )
    {

        // Declare moment coefficient vector and intialize to zeros.
        Eigen::Vector3d momentCoefficients = Eigen::Vector3d::Zero( );

        // Declare moment arm for panel moment determination.
        Eigen::Vector3d referenceDistance;

        // Loop over all panels and add moments due pressures.
        for ( int i = 0 ; i < vehiclePanels_.size( ) ; i++ )
        {
            // Determine moment arm for given panel centroid.
            referenceDistance = ( vehiclePanels_.at( i )->getPanelCentroid(  ) -
                                  this->momentReferencePoint_ );

            momentCoefficients -=
                currentPressureCoefficients_.at( i ) *
                vehiclePanels_.at( i )->getPanelArea( ) *
                    referenceDistance.cross( vehiclePanels_.at( i )->getFrameFixedSurfaceNormal( )( ) );
        }


        // Scale result by reference length and area.
        momentCoefficients /= ( this->referenceLength_ * this->referenceArea_ );

        return momentCoefficients;
    }


    //! Multi-array as which indicates which coefficients have been calculated already.
    /*!
     * Multi-array as which indicates which coefficients have been calculated already. Indices of
     * entries coincide with indices of aerodynamicCoefficients_.
     */
    boost::multi_array< bool, NumberOfIndependentVariables > isCoefficientGenerated_;
    
    std::vector< double > currentPanelInclinations_;

    //! Map of angle of attack and -sideslip pair and associated panel inclinations.
    /*!
     * Map of angle of attack and -sideslip pair and associated panel inclinations.
     */
    std::map< std::pair< double, double >, std::vector< double > > previouslyComputedInclinations_;
    
    std::vector< double > currentPressureCoefficients_;

    std::map< boost::array< int, NumberOfIndependentVariables >, std::vector< double >  > pressureCoefficientList_;
    

    std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > > vehiclePanels_;
    
    std::vector< std::vector< double > > dataPointsOfIndependentVariables_;

    boost::array< int, NumberOfIndependentVariables > numberOfPointsPerIndependentVariables_;

    int angleOfAttackIndex_;
    
    int sideslipAngleIndex_;
};


} // namespace aerodynamics
} // namespace tudat

#endif // TUDAT_RAREFIEDFLOW_AERODYNAMIC_COEFFICIENT_GENERATOR_H
