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

#ifndef TUDAT_RAREFIEDFLOW_AERODYNAMIC_COEFFICIENT_INTERFACE_H
#define TUDAT_RAREFIEDFLOW_AERODYNAMIC_COEFFICIENT_INTERFACE_H

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
        const double currentTime = TUDAT_NAN ) = 0;
            

    private:
    //! Determines variables to define the incinations required to calculate the aerodynamic force coefficients
    
    void determineIncinations(    )

    





    // Declaration of member variables

    //! Vehicle panels
    std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > > vehicleExteriorPanels_;

    //! Vehicle part orientation
    std::map< std::string, std::shared_ptr< ephemerides::RotationalEphemeris > > vehiclePartOrientation_;

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
