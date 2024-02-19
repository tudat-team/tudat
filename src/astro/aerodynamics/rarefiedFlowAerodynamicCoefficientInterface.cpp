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

#include <map>
#include <string>
#include <vector>

#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <memory>

#include <Eigen/Core>

#include "tudat/astro/aerodynamics/aerodynamicCoefficientInterface.h"
#include "tudat/basics/basicTypedefs.h"
#include "tudat/astro/system_models/vehicleExteriorPanels.h"
#include "tudat/astro/aerodynamics/rarefiedFlowInteractionModel.h"

namespace tudat
{
namespace aerodynamics
{

template< unsigned int NumberOfIndependentVariables >
class RarefiedFlowAerodynamicCoefficientInterface: public AerodynamicCoefficientInterface
{
public:
    /*!
     * Constructor of rarefied flow aerodynamic coefficient interface.
     * \param vehicleExteriorPanels Vehicle panels
     * \param vehiclePartOrientation Vehicle part orientation
     * \param referenceLength Reference length
     * \param referenceArea Reference area
     * \param momentReferencePoint Moment reference point
     * \param independentVariableNames Independent variable names
     * \param forceCoefficientsFrame Force coefficients frame
     * \param momentCoefficientsFrame Moment coefficients frame
     * \param accountForShadedPanels Account for shaded panels
     * \param dataPointsOfInclinationsForShading Data points of inclinations for shading
     */
    RarefiedFlowAerodynamicCoefficientInterface(
        const std::map< std::string, std::vector< std::shared_ptr< VehicleExteriorPanel > > > vehicleExteriorPanels,
        const std::map< std::string, std::shared_ptr< ephemerides::RotationalEphemeris > > vehiclePartOrientation,
        const double referenceLength,
        const double referenceArea,
        const Eigen::Vector3d& momentReferencePoint,
        const std::vector< AerodynamicCoefficientsIndependentVariables > independentVariableNames,
        const AerodynamicCoefficientFrames forceCoefficientsFrame = negative_aerodynamic_frame_coefficients,
        const AerodynamicCoefficientFrames momentCoefficientsFrame = body_fixed_frame_coefficients,
        const bool accountForShadedPanels = false,
        const std::map< int, std::vector< double > > dataPointsOfInclinationsForShading = std::map< int, std::vector< double > >( ),
        ) : 
        AerodynamicCoefficientInterface(
            referenceLength, referenceLength, momentReferencePoint, independentVariableNames, forceCoefficientsFrame, momentCoefficientsFrame
            ),
        //defining the member variables
        vehicleExteriorPanels_( vehicleExteriorPanels ), vehiclePartOrientation_( vehiclePartOrientation ), referenceLength_( referenceLength ), referenceArea_( referenceArea ),
        momentReferencePoint_( momentReferencePoint ), independentVariableNames_( independentVariableNames ), forceCoefficientsFrame_( forceCoefficientsFrame ),
        momentCoefficientsFrame_( momentCoefficientsFrame ), accountForShadedPanels_( accountForShadedPanels ), dataPointsOfInclinationsForShading_( dataPointsOfInclinationsForShading )
    {

        
        

    }



    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~RarefiedFlowAerodynamicCoefficientGenerator( ) { }   

    //! Compute the aerodynamic coefficients of the body itself (without control surfaces) at current flight condition.
    /*!
     *  Computes the current force and moment coefficients of the body itself (without control surfaces) and is to be
     *  implemented in derived classes. Input is a set of independent variables
     *  (doubles) which represent the variables from which the coefficients are calculated
     *  \param independentVariables Independent variables of force and moment coefficient
     *  determination implemented by derived class
     *  \param currentTime Time to which coefficients are to be updated
     */
    virtual void updateCurrentCoefficients(
        const std::vector< double >& independentVariables,
        const double currentTime = TUDAT_NAN )
    {

        determineIncinations( currentTime, independentVariables.at(0), independentVariables.at(1) );

        std::vector< double > species_number_densities;

        for ( int i = 4; i < independentVariables.size(); i++ )
        {
            species_number_densities.push_back( independentVariables.at(i) );
        }

        determinePanelForceCoefficientVectors( independentVariables.at(2), independentVariables.at(3), species_number_densities );

        determinePanelMomentCoefficientVectors( currentTime );

        // sum all panel force and moment coefficient vectors to get the total force and moment coefficient vectors
        Eigen::Vector6d totalAerodynamicCoefficients = Eigen::Vector6d::Zero();

        for ( std::string vehiclePartName: vehicleExteriorPanels_ )
        {
            for ( int i = 0; i < vehicleExteriorPanels_.at( vehiclePartName ).size(); i++ )
            {
                totalAerodynamicCoefficients.segment( 0, 3 ) += vehiclePanelForceCoefficientVectors_[ vehiclePartName ].at( i );
                totalAerodynamicCoefficients.segment( 3, 7 ) += vehiclePanelMomentCoefficientVectors_[ vehiclePartName ].at( i );
            }
        }

    }

    private:

    void determineIncinations(
        double secondsSinceEpoch,
        double angleOfAttack,
        double angleOfSideslip,
    ){
        // Declare free-stream velocity unit vector

        Eigen::Vector3d freestreamVelocityDirection = Eigen::Vector3d::Zero();

        freestreamVelocityDirection(0) = std::cos( angleOfAttack ) * std::cos( angleOfSideslip );
        freestreamVelocityDirection(1) = std::sin( angleOfSideslip );
        freestreamVelocityDirection(2) = std::sin( angleOfAttack ) * std::cos( angleOfSideslip );

        // Drag unit vector is parallel to the freestream velocity direction

        dragUnitVector_ = freestreamVelocityDirection.normalized();

        // Loop over all vehicle part names in vehicleExteriorPanels_
        for ( std::string vehiclePartName: vehicleExteriorPanels_ )
        {
            // Loop over all vehicle panels in vehicleExteriorPanels_ for this vehicle part
            for ( std::shared_ptr< VehicleExteriorPanel > vehiclePanel: vehicleExteriorPanels_.at( vehiclePartName ) )
            {
                
                // Determine panel normal vector and rotate to base frame
                Eigen::Vector3d panelNormalVector =  vehiclePartOrientation_.at( vehiclePartName )->getRotationMatrixToBaseFrame(secondsSinceEpoch) * vehiclePanel->getFrameFixedSurfaceNormal();

                // Determine cosine of angle between panel normal vector and freestream velocity direction (drag)
                // gamma_i in Doornbos (2011)
                double panelCosineDragAngle = -dragUnitVector.dot( panelNormalVector );

                // Determine lift unit vector for this panel
                Eigen::Vector3d liftUnitVector = -( ( dragUnitVector.cross( panelNormalVector ) ).cross( dragUnitVector ) ).normalized();

                // Store lift unit vector
                vehiclePanelLiftUnitVectors_[ vehiclePartName ].push_back( liftUnitVector );

                // Determine cosine of angle between panel normal vector and lift unit vector
                // l_i in Doornbos (2011)
                double panelCosineLiftAngle = -liftUnitVector.dot( panelNormalVector );

                // Store panel cosines of lift and drag angles
                vehiclePanelCosinesOfLiftAndDragAngles_[ vehiclePartName ].push_back( std::make_pair( panelCosineLiftAngle, panelCosineDragAngle ) );

            }
        }
    }

    void determinePanelForceCoefficientVectors(
        double Vinf, double Tinf, double species_number_densities
    ){

        double total_number_density = 0.0;
        for (int j_species = 0; j_species < species_number_densities.size(); j_species++)
        {
            total_number_density += species_number_densities[j_species];
        }

        // Loop over all vehicle part names in vehicleExteriorPanels_
        for ( std::string vehiclePartName: vehicleExteriorPanels_ )
        {
            // Loop over all vehicle panels in vehicleExteriorPanels_ for this vehicle part
            for ( std::shared_ptr< VehicleExteriorPanel > vehiclePanel: vehicleExteriorPanels_.at( vehiclePartName ) )
            {
                
                vehiclePanelForceCoefficientVectors_[ vehiclePartName ].push_back( 
                    RarefiedFlowInteractionModel::computePanelForceCoefficientVector( 
                        vehiclePanelCosinesOfLiftAndDragAngles_[ vehiclePartName ].at( vehiclePanel->getPanelIndex() ).first,
                        vehiclePanelCosinesOfLiftAndDragAngles_[ vehiclePartName ].at( vehiclePanel->getPanelIndex() ).second,
                        vehiclePanel->getPanelArea(),
                        vehiclePanelLiftUnitVectors_[ vehiclePartName ].at( vehiclePanel->getPanelIndex() ),
                        dragUnitVector_,
                        Vinf,
                        Tinf,
                        species_number_densities,
                        total_number_density,
                        referenceArea_
                    )
                )

            }
        }
    }   

    void determinePanelMomentCoefficientVectors(double secondsSinceEpoch){
            
        // Loop over all vehicle part names in vehicleExteriorPanels_
        for ( std::string vehiclePartName: vehicleExteriorPanels_ )
        {
            // Loop over all vehicle panels in vehicleExteriorPanels_ for this vehicle part
            for ( std::shared_ptr< VehicleExteriorPanel > vehiclePanel: vehicleExteriorPanels_.at( vehiclePartName ) )
            {
                // Calculate panel position vector in base frame
                Eigen::Vector3d panelPositionVector = 
                    vehiclePartOrientation_.at( vehiclePartName )->getRotationMatrixToBaseFrame(secondsSinceEpoch) * vehiclePanel->getPanelPositionVector()
                    + vehiclePartOrientation_.at( vehiclePartName )->getFrameFixedPositionVector();
                
                vehiclePanelMomentCoefficientVectors_[ vehiclePartName ].push_back( 
                    RarefiedFlowInteractionModel::computePanelMomentCoefficientVector( 
                        vehiclePanelForceCoefficientVectors_[ vehiclePartName ].at( vehiclePanel->getPanelIndex() ),
                        panelPositionVector,
                        referenceLength_
                    )
                )
            }
        }
    }


    // Declaration of member variables

    //! Vehicle panels
    std::map< std::string, std::vector< std::shared_ptr< VehicleExteriorPanel > > > vehicleExteriorPanels_;

    //! Vehicle part orientation
    std::map< std::string, std::shared_ptr< ephemerides::RotationalEphemeris > > vehiclePartOrientation_;

    //! Vehicle panel cosines of lift and drag angles
    std::map< std::string, std::vector< std::pair< double, double > > > vehiclePanelCosinesOfLiftAndDragAngles_;

    //! Vehicle panel lift unit vecotrs
    std::map< std::string, std::vector< Eigen::Vector3d > > vehiclePanelLiftUnitVectors_;

    //! drag unit vector
    Eigen::Vector3d dragUnitVector_;

    //! Vehicle panel force coefficient vectors
    std::map< std::string, std::vector< Eigen::Vector3d > > vehiclePanelForceCoefficientVectors_;

    //! Reference length
    double referenceLength_;
    //! Reference area
    double referenceArea_;
    //! Moment reference point
    Eigen::Vector3d momentReferencePoint_;

    //! Independent variable names
    std::vector< AerodynamicCoefficientsIndependentVariables > independentVariableNames_;

    //! Force coefficients frame
    AerodynamicCoefficientFrames forceCoefficientsFrame_;
    //! Moment coefficients frame
    AerodynamicCoefficientFrames momentCoefficientsFrame_;

    //! Account for shaded panels
    bool accountForShadedPanels_;

    //! Data points of inclinations for shading
    std::map< int, std::vector< double > > dataPointsOfInclinationsForShading_;



}


} // namespace aerodynamics
} // namespace tudat

#endif // TUDAT_RAREFIEDFLOW_AERODYNAMIC_COEFFICIENT_GENERATOR_H
