/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_RADIATIONPRESSURECOEFFICIENT_H
#define TUDAT_RADIATIONPRESSURECOEFFICIENT_H

#include <cmath>

#include <Eigen/Core>

#include "tudat/basics/utilities.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"
#include "tudat/astro/electromagnetism/radiationPressureInterface.h"
#include "tudat/math/interpolators/piecewiseConstantInterpolator.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Interface class for the estimation of a radiation pressure coefficient
class RadiationPressureCoefficient : public EstimatableParameter< double >
{
public:
    //! Constructor.
    /*!
     * Constructor
     * \param radiationPressureInterface Object containing the radiation pressure coefficient to be estimated.
     * \param associatedBody Name of body containing the radiationPressureInterface object
     */
    RadiationPressureCoefficient(
            const std::shared_ptr< electromagnetism::CannonballRadiationPressureTargetModel > radiationPressureInterface,
            const std::string& associatedBody ):
        EstimatableParameter< double >( radiation_pressure_coefficient, associatedBody ),
        radiationPressureInterface_( radiationPressureInterface )
    {
        if( std::isnan( radiationPressureInterface_->getCoefficient( ) ) )
        {
            throw std::runtime_error( "Error when creating estimated constant Cr coefficient for " + associatedBody +
                                      ", current Cr not initialized" );
        }

        if( radiationPressureInterface_->getCoefficientFunction( ) != nullptr )
        {
            throw std::runtime_error( "Error when creating estimated constant Cr coefficient for " + associatedBody +
                                      ", time-variable Cr function defined" );
        }
    }

    //! Destructor.
    ~RadiationPressureCoefficient( ) { }

    //! Function to get the current value of the radiation pressure coefficient that is to be estimated.
    /*!
     * Function to get the current value of the radiation pressure coefficient that is to be estimated.
     * \return Current value of the radiation pressure coefficient that is to be estimated.
     */
    double getParameterValue( )
    {
        return radiationPressureInterface_->getCoefficient( );
    }

    //! Function to reset the value of the radiation pressure coefficient that is to be estimated.
    /*!
     * Function to reset the value of the radiation pressure coefficient that is to be estimated.
     * \param parameterValue New value of the radiation pressure coefficient that is to be estimated.
     */
    void setParameterValue( double parameterValue )
    {
        radiationPressureInterface_->resetCoefficient( parameterValue );
    }

    //! Function to retrieve the size of the parameter (always 1).
    /*!
     *  Function to retrieve the size of the parameter (always 1).
     *  \return Size of parameter value (always 1).
     */
    int getParameterSize( )
    {
        return 1;
    }

protected:
private:
    //! Object containing the radiation pressure coefficient to be estimated.
    std::shared_ptr< electromagnetism::CannonballRadiationPressureTargetModel > radiationPressureInterface_;
};

class RadiationPressureScalingFactor : public EstimatableParameter< double >
{
public:
    RadiationPressureScalingFactor( const std::shared_ptr< electromagnetism::RadiationPressureAcceleration > radiationPressureAcceleration,
                                    const EstimatebleParametersEnum parameterType,
                                    const std::string& associatedBody,
                                    const std::string& exertingBody ):
        EstimatableParameter< double >( parameterType, associatedBody, exertingBody ),
        radiationPressureAcceleration_( radiationPressureAcceleration )
    {
        if( ( parameterType != source_direction_radiation_pressure_scaling_factor ) &&
            ( parameterType != source_perpendicular_direction_radiation_pressure_scaling_factor ) )
        {
            throw std::runtime_error( "Error when creating radiation pressure scaling parameter, type is inconsistent: " +
                                      std::to_string( parameterType ) );
        }
    }

    ~RadiationPressureScalingFactor( ) { }

    double getParameterValue( )
    {
        if( parameterName_.first == source_direction_radiation_pressure_scaling_factor )
        {
            return radiationPressureAcceleration_->getSourceDirectionScaling( );
        }
        else if( parameterName_.first == source_perpendicular_direction_radiation_pressure_scaling_factor )
        {
            return radiationPressureAcceleration_->getPerpendicularSourceDirectionScaling( );
        }
        else
        {
            throw std::runtime_error( "Error when getting radiation pressure scaling parameter, type is inconsistent: " +
                                      std::to_string( parameterName_.first ) );
        }
    }

    void setParameterValue( double parameterValue )
    {
        if( parameterName_.first == source_direction_radiation_pressure_scaling_factor )
        {
            radiationPressureAcceleration_->setSourceDirectionScaling( parameterValue );
        }
        else if( parameterName_.first == source_perpendicular_direction_radiation_pressure_scaling_factor )
        {
            radiationPressureAcceleration_->setPerpendicularSourceDirectionScaling( parameterValue );
        }
        else
        {
            throw std::runtime_error( "Error when setting radiation pressure scaling parameter, type is inconsistent: " +
                                      std::to_string( parameterName_.first ) );
        }
    }

    int getParameterSize( )
    {
        return 1;
    }

protected:
private:
    std::shared_ptr< electromagnetism::RadiationPressureAcceleration > radiationPressureAcceleration_;
};

//! Interface class for the estimation of an arc-wise (piecewise constant) radiation pressure coefficient
class ArcWiseRadiationPressureCoefficient : public EstimatableParameter< Eigen::VectorXd >
{
public:
    //! Constructor.
    /*!
     * Constructor
     * \param radiationPressureInterface Object containing the radiation pressure coefficient to be estimated.
     * \param timeLimits Times at which the arcs are to start.
     * \param associatedBody Name of body containing the radiationPressureInterface object
     */
    ArcWiseRadiationPressureCoefficient(
            const std::shared_ptr< electromagnetism::CannonballRadiationPressureTargetModel > radiationPressureInterface,
            const std::vector< double > timeLimits,
            const std::string& associatedBody ):
        EstimatableParameter< Eigen::VectorXd >( arc_wise_radiation_pressure_coefficient, associatedBody ),
        radiationPressureInterface_( radiationPressureInterface ), timeLimits_( timeLimits )
    {
        if( std::isnan( radiationPressureInterface_->getCoefficient( ) ) )
        {
            throw std::runtime_error( "Error when creating estimated arcwise Cr coefficient for " + associatedBody +
                                      ", current Cr not initialized" );
        }

        if( radiationPressureInterface_->getCoefficientFunction( ) != nullptr )
        {
            throw std::runtime_error( "Error when creating estimated arcwise Cr coefficient for " + associatedBody +
                                      ", time-variable Cr function defined" );
        }

        double radiationPressureCoefficient = radiationPressureInterface->getCoefficient( );
        for( unsigned int i = 0; i < timeLimits.size( ); i++ )
        {
            radiationPressureCoefficients_.push_back( radiationPressureCoefficient );
        }

        timeLimits_.push_back( std::numeric_limits< double >::max( ) );
        fullRadiationPressureCoefficients_ = radiationPressureCoefficients_;
        fullRadiationPressureCoefficients_.push_back( radiationPressureCoefficient );

        coefficientInterpolator_ = std::make_shared< interpolators::PiecewiseConstantInterpolator< double, double > >(
                timeLimits_, fullRadiationPressureCoefficients_ );

        typedef interpolators::OneDimensionalInterpolator< double, double > LocalInterpolator;
        radiationPressureInterface->resetCoefficientFunction(
                std::bind( static_cast< double ( LocalInterpolator::* )( const double ) >( &LocalInterpolator::interpolate ),
                           coefficientInterpolator_,
                           std::placeholders::_1 ) );
    }

    //! Destructor.
    ~ArcWiseRadiationPressureCoefficient( ) { }

    //! Function to get the current value of the radiation pressure coefficient that is to be estimated.
    /*!
     * Function to get the current value of the radiation pressure coefficient that is to be estimated.
     * \return Current value of the radiation pressure coefficient that is to be estimated.
     */
    Eigen::VectorXd getParameterValue( )
    {
        return utilities::convertStlVectorToEigenVector( radiationPressureCoefficients_ );
    }

    //! Function to reset the value of the radiation pressure coefficient that is to be estimated.
    /*!
     * Function to reset the value of the radiation pressure coefficient that is to be estimated.
     * \param parameterValue New value of the radiation pressure coefficient that is to be estimated.
     */
    void setParameterValue( Eigen::VectorXd parameterValue )
    {
        if( static_cast< int >( radiationPressureCoefficients_.size( ) ) != static_cast< int >( parameterValue.rows( ) ) )
        {
            throw std::runtime_error( "Error when resetting arc-wise radiation pressure coefficients, sizes are incompatible" );
        }

        radiationPressureCoefficients_ = utilities::convertEigenVectorToStlVector( parameterValue );
        for( unsigned int i = 0; i < radiationPressureCoefficients_.size( ); i++ )
        {
            fullRadiationPressureCoefficients_[ i ] = radiationPressureCoefficients_[ i ];
        }
        fullRadiationPressureCoefficients_[ radiationPressureCoefficients_.size( ) ] =
                radiationPressureCoefficients_.at( radiationPressureCoefficients_.size( ) - 1 );
        coefficientInterpolator_->resetDependentValues( fullRadiationPressureCoefficients_ );
    }

    //! Function to retrieve the size of the parameter
    /*!
     *  Function to retrieve the size of the parameter
     *  \return Size of parameter value
     */
    int getParameterSize( )
    {
        return radiationPressureCoefficients_.size( );
    }

    std::shared_ptr< interpolators::LookUpScheme< double > > getArcTimeLookupScheme( )
    {
        return coefficientInterpolator_->getLookUpScheme( );
    }

protected:
private:
    //! Object containing the radiation pressure coefficient to be estimated.
    std::shared_ptr< electromagnetism::CannonballRadiationPressureTargetModel > radiationPressureInterface_;

    //! Times at which the arcs are to start (including end time at maximum double value).
    std::vector< double > timeLimits_;

    //! Values of radiation pressure coefficients in each arc.
    std::vector< double > radiationPressureCoefficients_;

    //! Values of radiation pressure coefficients in each arc, with additional value copied at end.
    std::vector< double > fullRadiationPressureCoefficients_;

    //! Interpolator that returns the radiation pressure coefficient as a function of time.
    std::shared_ptr< interpolators::PiecewiseConstantInterpolator< double, double > > coefficientInterpolator_;
};

class ArcWiseRadiationPressureScalingFactor : public EstimatableParameter< Eigen::VectorXd >
{
public:

    ArcWiseRadiationPressureScalingFactor(
        const std::shared_ptr< electromagnetism::RadiationPressureAcceleration >& acceleration,
        const std::vector< double >& timeLimits,
        const EstimatebleParametersEnum parameterType,
        const std::string& associatedBody,
        const std::string& exertingBody )
        : EstimatableParameter< Eigen::VectorXd >( parameterType, associatedBody, exertingBody ),
          acceleration_( acceleration ), timeLimits_( timeLimits )
    {

        if (parameterType != arcwise_source_direction_radiation_pressure_scaling_factor &&
            parameterType != arcwise_source_perpendicular_direction_radiation_pressure_scaling_factor)
        {
            throw std::runtime_error("Inconsistent parameter type for arc-wise radiation pressure scaling factor");
        }

        parameterSize_ = static_cast< int >( timeLimits.size() );

        // Ensure coverage for all time
        timeLimits_.push_back(std::numeric_limits<double>::max());

        double initialScaling = (parameterType == arcwise_source_direction_radiation_pressure_scaling_factor) ?
                                    acceleration_->getSourceDirectionScaling() :
                                    acceleration_->getPerpendicularSourceDirectionScaling();

        parameterValues_ = std::vector< double >( parameterSize_, initialScaling );
        fullParameterValues_ = parameterValues_;
        fullParameterValues_.push_back(initialScaling);  // For extrapolation past last arc

        interpolator_ = std::make_shared< interpolators::PiecewiseConstantInterpolator< double, double > >(
            timeLimits_, fullParameterValues_ );

        typedef interpolators::OneDimensionalInterpolator< double, double > LocalInterpolator;
        std::function< double( double ) > interpolatorFunction =
            std::bind( static_cast< double ( LocalInterpolator::* )( const double ) >( &LocalInterpolator::interpolate ),
                    interpolator_, std::placeholders::_1 );

        if (parameterType == arcwise_source_direction_radiation_pressure_scaling_factor)
        {
            acceleration_->setSourceDirectionScalingFunction(interpolatorFunction);
        }
        else
        {
            acceleration_->setPerpendicularSourceDirectionScalingFunction(interpolatorFunction);
        }
    }

    //! Destructor
    ~ArcWiseRadiationPressureScalingFactor() {}

    //! Get current parameter value as an Eigen vector
    Eigen::VectorXd getParameterValue( ) override
    {
        return utilities::convertStlVectorToEigenVector( parameterValues_ );
    }

    //! Get the lookup scheme of the interpolator (useful for unit tests)
    std::shared_ptr< interpolators::LookUpScheme< double > > getArcTimeLookupScheme( )
    {
        return interpolator_->getLookUpScheme( );
    }

    //! Set parameter values and update interpolator
    void setParameterValue(Eigen::VectorXd parameterValue) override
    {
        if (parameterValue.size() != parameterSize_)
        {
            throw std::runtime_error("Error: Arc-wise scaling vector size mismatch.");
        }

        parameterValues_ = utilities::convertEigenVectorToStlVector(parameterValue);
        for (int i = 0; i < parameterSize_; i++)
        {
            fullParameterValues_[i] = parameterValues_[i];
        }

        // Last value again for overflow
        fullParameterValues_[parameterSize_] = parameterValues_.back();

        interpolator_->resetDependentValues(fullParameterValues_);
    }

    int getParameterSize( ) override
    {
        return parameterSize_;
    }

private:

    std::shared_ptr< electromagnetism::RadiationPressureAcceleration > acceleration_;
    std::vector< double > timeLimits_;
    std::vector< double > parameterValues_;
    std::vector< double > fullParameterValues_;
    std::shared_ptr< interpolators::PiecewiseConstantInterpolator< double, double > > interpolator_;

    int parameterSize_;
};


}  // namespace estimatable_parameters

}  // namespace tudat

#endif  // TUDAT_RADIATIONPRESSURECOEFFICIENT_H
