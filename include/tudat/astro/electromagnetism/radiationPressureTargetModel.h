/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Oliver Montenbruck, et al. Satellite Orbits. Springer Berlin Heidelberg, 2000.
 */

#ifndef TUDAT_RADIATIONPRESSURETARGETMODEL_H
#define TUDAT_RADIATIONPRESSURETARGETMODEL_H

#include <map>
#include <memory>
#include <vector>

#include <Eigen/Core>

#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/astro/electromagnetism/reflectionLaw.h"
#include "tudat/astro/system_models/vehicleExteriorPanels.h"


namespace tudat
{
namespace electromagnetism
{

/*!
 * Class modeling a target that is accelerated by radiation pressure.
 *
 * Parameters to all functions are in local target-fixed coordinates. Proper translation and rotation is ensured by the
 * RadiationPressureAcceleration class.
 */
class RadiationPressureTargetModel
{
public:
    explicit RadiationPressureTargetModel(
        const std::map<std::string, std::vector<std::string>>& sourceToTargetOccultingBodies = {}) :
        sourceToTargetOccultingBodies_(sourceToTargetOccultingBodies), radiationPressure_( TUDAT_NAN ) {}

    virtual ~RadiationPressureTargetModel() = default;

    void updateMembers(double currentTime);

    /*!
     * Calculate radiation pressure force from incident radiation using geometrical/optical target properties.
     *
     * @param sourceIrradiance Incident irradiance magnitude [W/m²]
     * @param sourceToTargetDirectionLocalFrame Direction of incoming radiation
     * @return Radiation pressure force vector in local (i.e. target-fixed) coordinates [N]
     */
    virtual Eigen::Vector3d evaluateRadiationPressureForce(
            const double sourceIrradiance, const Eigen::Vector3d& sourceToTargetDirection) = 0;
    
    std::map<std::string, std::vector<std::string>> getSourceToTargetOccultingBodies() const
    {
        return sourceToTargetOccultingBodies_;
    }

    virtual bool forceFunctionRequiresLocalFrameInputs( ) = 0;

    double getRadiationPressure() const
    {
        return radiationPressure_;
    }

protected:
    virtual void updateMembers_(const double currentTime) {};

    double currentTime_{TUDAT_NAN};
    // Only needed to transfer occultation settings from body setup to acceleration setup
    std::map<std::string, std::vector<std::string>> sourceToTargetOccultingBodies_;

    mutable double radiationPressure_;
};

/*!
 * Class modeling a target as sphere ("cannonball"). The sphere has an isotropic radiation pressure coefficient and
 * the same cross-sectional area from any direction.
 */
class CannonballRadiationPressureTargetModel : public RadiationPressureTargetModel
{
public:
    /*!
     * Constructor.
     *
     * @param area Cross-sectional area (i.e. projected area of the sphere)
     * @param coefficient Radiation pressure coefficient (between 1 [pure absorption] and 2 [pure specular reflection])
     * @param sourceToTargetOccultingBodies Map (source name -> list of occulting body names) of bodies
     *      to occult sources as seen from this target
     */
    CannonballRadiationPressureTargetModel(
        double area, double coefficient,
        const std::map<std::string, std::vector<std::string>>& sourceToTargetOccultingBodies = {}) :
        RadiationPressureTargetModel(sourceToTargetOccultingBodies),
        area_(area), coefficientFunction_( nullptr ),
        currentCoefficient_(coefficient){ }

    CannonballRadiationPressureTargetModel(
        double area, std::function< double( const double ) > coefficientFunction,
        const std::map<std::string, std::vector<std::string>>& sourceToTargetOccultingBodies = {}) :
        RadiationPressureTargetModel(sourceToTargetOccultingBodies),
        area_(area), coefficientFunction_( coefficientFunction),
        currentCoefficient_( TUDAT_NAN ){ }

    Eigen::Vector3d evaluateRadiationPressureForce(
            double sourceIrradiance,
            const Eigen::Vector3d& sourceToTargetDirection ) override;

    double getArea() const
    {
        return area_;
    }

    double getCoefficient() const
    {
        return currentCoefficient_;
    }

    void resetCoefficient( const double coefficient )
    {
        currentCoefficient_ = coefficient;
    }

    std::function< double( const double ) > getCoefficientFunction( )
    {
        return coefficientFunction_;
    }

    void resetCoefficientFunction( const std::function< double( const double ) > coefficientFunction )
    {
        coefficientFunction_ = coefficientFunction;
    }

    bool forceFunctionRequiresLocalFrameInputs( ) override
    {
        return false;
    }

private:

    virtual void updateMembers_(const double currentTime) override
    {
        if( coefficientFunction_ != nullptr )
        {
            currentCoefficient_ = coefficientFunction_( currentTime );
        }
        radiationPressure_ = 0.0;
    };

    double area_;

    std::function< double( const double ) > coefficientFunction_;

    double currentCoefficient_;
};

/*!
 * Class modeling a target as collection of panels, e.g., representing the box body and solar panels.
 */
class PaneledRadiationPressureTargetModel : public RadiationPressureTargetModel
{
public:

    /*!
     * Constructor.
     *
     * @param panels Panels comprising this paneled target
     * @param sourceToTargetOccultingBodies Map (source name -> list of occulting body names) of bodies
     *      to occult sources as seen from this target
     */
    explicit PaneledRadiationPressureTargetModel(
        const std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > >& bodyFixedPanels,
        const std::map<std::string, std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > > >& segmentFixedPanels =
            std::map<std::string, std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > > >( ),
        const std::map<std::string, std::function< Eigen::Quaterniond( ) > >& segmentFixedToBodyFixedRotations =
            std::map<std::string, std::function< Eigen::Quaterniond( ) > >( ),
        const std::map<std::string, std::vector<std::string> >& sourceToTargetOccultingBodies = { } ) :
        RadiationPressureTargetModel( sourceToTargetOccultingBodies ),
        bodyFixedPanels_( bodyFixedPanels ),
        segmentFixedPanels_( segmentFixedPanels ),
        segmentFixedToBodyFixedRotations_( segmentFixedToBodyFixedRotations )
    {
        totalNumberOfPanels_ = bodyFixedPanels_.size( );
        fullPanels_ = bodyFixedPanels_;

        int counter = 1;
        for( auto it : segmentFixedPanels_ )
        {
            totalNumberOfPanels_ += it.second.size( );
            fullPanels_.insert( fullPanels_.end( ), it.second.begin( ), it.second.end( ) );
        }

        panelForces_.resize( totalNumberOfPanels_ );
        surfacePanelCosines_.resize( totalNumberOfPanels_ );
        surfaceNormals_.resize( totalNumberOfPanels_ );
    }

    Eigen::Vector3d evaluateRadiationPressureForce(
        double sourceIrradiance,
        const Eigen::Vector3d &sourceToTargetDirectionLocalFrame ) override;

    Eigen::Vector3d evaluateRadiationPressureForcePartialWrtDiffuseReflectivity(
        double sourceIrradiance,
        const Eigen::Vector3d &sourceToTargetDirectionLocalFrame );

    Eigen::Vector3d evaluateRadiationPressureForcePartialWrtSpecularReflectivity(
        double sourceIrradiance,
        const Eigen::Vector3d &sourceToTargetDirectionLocalFrame );

    std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > >& getBodyFixedPanels( )
    {
        return bodyFixedPanels_;
    }

    std::map< std::string, std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > > >& getSegmentFixedPanels( )
    {
        return segmentFixedPanels_;
    }

    bool forceFunctionRequiresLocalFrameInputs( ) override
    {
        return true;
    }

    std::vector< Eigen::Vector3d >& getSurfaceNormals( )
    {
        return surfaceNormals_;
    }


    std::vector< double >& getSurfacePanelCosines( )
    {
        return surfacePanelCosines_;
    }

    std::vector< Eigen::Vector3d >& getPanelForces( )
    {
        return panelForces_;
    }

    std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > >& getFullPanels( )
    {
        return fullPanels_;
    }


    int getTotalNumberOfPanels( )
    {
        return totalNumberOfPanels_;
    }

    double getCoefficient() const
    {
        return currentCoefficient_;
    }

    void resetCoefficient( const double coefficient )
    {
        currentCoefficient_ = coefficient;
    }

    std::function< double( const double ) > getCoefficientFunction( )
    {
        return coefficientFunction_;
    }

    void resetCoefficientFunction( const std::function< double( const double ) > coefficientFunction )
    {
        coefficientFunction_ = coefficientFunction;
    }

    std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > > getPanelsFromId(
            const std::string& panelTypeId)
    {
        std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > > panels;
        for( unsigned int i = 0; i < static_cast<unsigned int>(totalNumberOfPanels_); i++ )
        {
            if(fullPanels_.at(i)->getPanelTypeId() == panelTypeId)
            {
                panels.push_back(fullPanels_.at(i));
            }
        }
        return panels;
    }

    double getAverageDiffuseReflectivity(const std::string& panelTypeId)
    {
        std::vector<double> coefficients;
        std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > > panelsFromId;
        panelsFromId = this->getPanelsFromId(panelTypeId);
        for( unsigned int i = 0; i < panelsFromId.size(); i++ ){
            {
                auto reflectionLaw = std::dynamic_pointer_cast<SpecularDiffuseMixReflectionLaw>(panelsFromId.at(i)->getReflectionLaw());
                coefficients.push_back(reflectionLaw->getDiffuseReflectivity());
            }
        }

        // Check if coefficients are consistent
        if (!coefficients.empty()) {
            double firstCoefficient = coefficients[0];
            for (size_t i = 1; i < coefficients.size(); ++i) {
                if (coefficients[i] != firstCoefficient) {
                    std::cerr << "Warning: Specular reflectivity coefficients for panel group "
                              << panelTypeId << " are not consistent. Returning average value" << std::endl;
                    break;  } } }

        // Calculate average coefficient
        double averageCoefficient = 0.0;
        if (!coefficients.empty())
        {
            averageCoefficient = std::accumulate(coefficients.begin(), coefficients.end(), 0.0) / coefficients.size();
        }

        return averageCoefficient;
    };


    double getAverageSpecularReflectivity(const std::string& panelTypeId)
    {
        std::vector<double> coefficients;
        std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > > panelsFromId;
        panelsFromId = this->getPanelsFromId(panelTypeId);
        for( unsigned int i = 0; i < panelsFromId.size(); i++ ){
            {
                auto reflectionLaw = std::dynamic_pointer_cast<SpecularDiffuseMixReflectionLaw>(panelsFromId.at(i)->getReflectionLaw());
                coefficients.push_back(reflectionLaw->getSpecularReflectivity());
            }
        }

        // Check if coefficients are consistent
        if (!coefficients.empty()) {
            double firstCoefficient = coefficients[0];
            for (size_t i = 1; i < coefficients.size(); ++i) {
                if (coefficients[i] != firstCoefficient) {
                    std::cerr << "Warning: Diffuse reflectivity coefficients for for panel group "
                              << panelTypeId << " are not consistent. Returning average value" << std::endl;
                    break;  } } }

        // Calculate average coefficient
        double averageCoefficient = 0.0;
        if (!coefficients.empty())
        {
            averageCoefficient = std::accumulate(coefficients.begin(), coefficients.end(), 0.0) / coefficients.size();
        }

        return averageCoefficient;
    };

    void setGroupSpecularReflectivity(const std::string& panelTypeId, double specularReflectivity)
    {
        std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > > panelsFromId;
        panelsFromId = this->getPanelsFromId(panelTypeId);
        for( unsigned int i = 0; i < panelsFromId.size(); i++ )
        {
            if(panelsFromId.at(i)->getReflectionLaw() != nullptr )
            {
                auto reflectionLaw = std::dynamic_pointer_cast<SpecularDiffuseMixReflectionLaw>(fullPanels_.at(i)->getReflectionLaw());
                reflectionLaw->setSpecularReflectivity(specularReflectivity);
            }
        }
    };
    void setGroupDiffuseReflectivity(const std::string& panelTypeId, double diffuseReflectivity)
    {
        std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > > panelsFromId;
        panelsFromId = this->getPanelsFromId(panelTypeId);
        for( unsigned int i = 0; i < panelsFromId.size(); i++ )
        {
            if(panelsFromId.at(i)->getReflectionLaw() != nullptr )
            {
                auto reflectionLaw = std::dynamic_pointer_cast<SpecularDiffuseMixReflectionLaw>(fullPanels_.at(i)->getReflectionLaw());
                reflectionLaw->setDiffuseReflectivity(diffuseReflectivity);
            }
        }
    };

    void setGroupAbsorptivity(const std::string& panelTypeId, double absorbtivity)
    {
        std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > > panelsFromId;
        panelsFromId = this->getPanelsFromId(panelTypeId);
        for( unsigned int i = 0; i < panelsFromId.size(); i++ )
        {
            if(panelsFromId.at(i)->getReflectionLaw() != nullptr )
            {
                auto reflectionLaw = std::dynamic_pointer_cast<SpecularDiffuseMixReflectionLaw>(fullPanels_.at(i)->getReflectionLaw());
                reflectionLaw->setAbsorptivity(absorbtivity);
            }
        }
    };


private:
    //void updateMembers_( double currentTime ) override;

    virtual void updateMembers_(const double currentTime) override
    {
        if( coefficientFunction_ != nullptr )
        {
            currentCoefficient_ = coefficientFunction_( currentTime );
        }
        radiationPressure_ = 0.0;
    };
    

    std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > > bodyFixedPanels_;

    std::map< std::string, std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > > > segmentFixedPanels_;

    std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > > fullPanels_;

    std::map< std::string, std::function< Eigen::Quaterniond( ) > > segmentFixedToBodyFixedRotations_;

    int totalNumberOfPanels_;
    
    std::vector< Eigen::Vector3d > surfaceNormals_;

    std::vector< double > surfacePanelCosines_;

    std::vector< Eigen::Vector3d > panelForces_;

    std::function< double( const double ) > coefficientFunction_;

    double currentCoefficient_ = 1;



};

} // tudat
} // electromagnetism

#endif //TUDAT_RADIATIONPRESSURETARGETMODEL_H
