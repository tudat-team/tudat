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
 *      Doornbos, E. N. Thermospheric Density and Wind Determination from Satellite Dynamics, 2011. 
 *      (Page 66+)
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
#include "tudat/astro/aerodynamics/rarefiedFlowInteractionModel.h"

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
            const bool savePressureCoefficients = false,
            ) : AerodynamicCoefficientGenerator< NumberOfIndependentVariables, 6 >(
        dataPointsOfIndependentVariables, referenceLength, referenceArea,
        momentReferencePoint, utilities::createVectorFromMapKeys( dataPointsOfIndependentVariables ),
        positive_aerodynamic_frame_coefficients, positive_aerodynamic_frame_coefficients ),
    vehiclePanels_( vehiclePanels ), dataPointsOfIndependentVariables_( utilities::createVectorFromMapValues( dataPointsOfIndependentVariables ) ),
    angleOfAttackIndex_( -1 ), sideslipAngleIndex_( -1 )
    {

        calculatePartForceCoefficientVector( );
        calculatePartMomentCoefficientVector( );   

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

        // define drag unit vector
        Eigen::Vector3d u_drag = freestreamVelocityDirection.normalized();
        // Declare lift unit vector
        Eigen::Vector3d u_lift;

        // Loop over all panels of given vehicle part and set inclination angles.
        for( unsigned int k = 0; k < vehiclePanels_.size( ); k++ )
        {

            // Declare cosine of angle between panel normal vector and drag unit vector (=velocity vector).
            double cosineOfNormalDragAngle;
            // Declare cosine of angle between panel normal vector and lift unit vector.
            double cosineOfNormalLiftAngle;

            cosineOfNormalDragAngle = -(u_drag).dot(vehiclePanels_[ k ]->getFrameFixedSurfaceNormal( )( )); // gammai in Doornbos

            u_lift = -(u_drag.cross(vehiclePanels_[ k ]->getFrameFixedSurfaceNormal( )( ))).cross(u_drag);
            u_lift.normalize();

            cosineOfNormalLiftAngle = -(u_lift).dot(vehiclePanels_[ k ]->getFrameFixedSurfaceNormal( )( )); // li in Doornbos

            // Set cosine of drag and lift angles (wrt to panel normal).
            currentCosineOfNormalDragAngles_[k] = cosineOfNormalDragAngle;
            currentCosineOfNormalLiftAngles_[k] = cosineOfNormalLiftAngle;
            
            // Set inclination angle -> unused by rarefied flow model.
            // currentPanelInclinations_[ k ] = mathematical_constants::PI / 2.0 - std::acos( cosineOfInclination );

        }
    }
    
    std::vector< std::vector< std::vector< double > > > getPressureCoefficientList(
            const boost::array< int, 3 > independentVariables )
    {
        return pressureCoefficientList_.at( independentVariables );
    }

    void clearData( )
    {
        currentCosineOfNormalDragAngles_.clear( );
        currentCosineOfNormalLiftAngles_.clear( );
        currentForceCoefficientvectors_.clear( );
        previouslyComputedInclinations_.clear( );
        this->clearBaseData( );
    }

private:
    
    //! Determine force coefficient vectors on a given part.
    /*!
     * Determines force coefficient vectors on a single vehicle part.
     * Calls something to be described
     * \param partNumber Index from vehicleParts_ array for which to determine coefficients.
     * \param independentVariableIndices Array of indices of independent variables.
     */
    void determinePanelForceCoefficientVectors( const boost::array< int, NumberOfIndependentVariables >& independentVariableIndices )
    {
        for( unsigned int i = 0; i < vehiclePanels_.size( ); i++ )
        {

            double currentCosineOfNormalDragAngle = currentCosineOfNormalDragAngles_.at( i );
            double currentCosineOfNormalLiftAngle = currentCosineOfNormalLiftAngles_.at( i );
            currentPanelForceCoefficientvectors_[ i ] = vehiclePanels_.at( i )->computePanelForceCoefficientVector( 
                currentCosineOfNormalDragAngle, currentCosineOfNormalLiftAngle, vehiclePanels_.at( i ).getPanelArea( ),
                u_lift, u_drag, Vinf, T_atm, number_densities, total_number_density, Aref);

        }
    }

    //! Determine force coefficients of a part.
    /*!
     * Sums the pressure coefficients of given part and determines force coefficients from it by
     * non-dimensionalization with reference area.
     * \param partNumber Index from vehicleParts_ array for which determine coefficients.
     * \return Force coefficients for requested vehicle part.
     */
    Eigen::Vector3d calculatePartForceCoefficientVector( const int partNumber )
    {
        // Declare force coefficient vector and intialize to zeros.
        Eigen::Vector3d partForceCoefficientVector = Eigen::Vector3d::Zero( );

        // Loop over all panels and add pressures, scaled by panel area, to force
        // coefficients.
        for ( int i = 0 ; i < vehiclePanels_.size( ) ; i++ )
        {
            partForceCoefficientVector += currentPanelForceCoefficientvectors_.at( i );
        }
        return partForceCoefficientVector;
    }

    //! Determine moment coefficients of a part.
    /*!
     * Determines the moment coefficients of a given part by summing the contributions of all
     * panels on the part. Moment arms are taken from panel centroid to momentReferencePoint. Non-
     * dimensionalization is performed by product of referenceLength and referenceArea.
     * \param partNumber Index from vehicleParts_ array for which to determine coefficients.
     * \return Moment coefficients for requested vehicle part.
     */
    Eigen::Vector3d calculatePartMomentCoefficientVector( const int partNumber )
    {

        // Declare moment coefficient vector and intialize to zeros.
        Eigen::Vector3d partMomentCoefficientVector = Eigen::Vector3d::Zero( );

        // Declare moment arm for panel moment determination.
        Eigen::Vector3d panelPositionVector;

        // Loop over all panels and add moments due pressures.
        for ( int i = 0 ; i < vehiclePanels_.size( ) ; i++ )
        {
            // Determine moment arm for given panel centroid.
            panelPositionVector = ( vehiclePanels_.at( i )->getPanelCentroid(  ) -
                                  this->momentReferencePoint_ );

            partMomentCoefficientVector += rarefiedFlowInteractionModel::computePanelMomentCoefficientVector( 
                currentForceCoefficientVectors_.at( i ), panelPositionVector, this->referenceLength_ );
        }

        return partMomentCoefficientVector;
    }


    
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
}
}


} // namespace aerodynamics
} // namespace tudat

#endif // TUDAT_RAREFIEDFLOW_AERODYNAMIC_COEFFICIENT_GENERATOR_H
