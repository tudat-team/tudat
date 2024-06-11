/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/acceleration_partials/fullRadiationPressureAccelerationPartial.h"

namespace tudat
{

namespace acceleration_partials
{

RadiationPressureAccelerationPartial::RadiationPressureAccelerationPartial(
    const std::string acceleratedBody,
    const std::string acceleratingBody,
    const std::shared_ptr< electromagnetism::PaneledSourceRadiationPressureAcceleration > radiationPressureAcceleration,
    const std::shared_ptr< estimatable_parameters::CustomSingleAccelerationPartialCalculatorSet > customAccelerationPartialSet ):
    AccelerationPartial( acceleratedBody, acceleratingBody, basic_astrodynamics::custom_acceleration ),
    radiationPressureAcceleration_( radiationPressureAcceleration ),
    customAccelerationPartialSet_( customAccelerationPartialSet )
{
    currentPartialWrtUndergoingState_.setZero( );
    currentPartialWrtExertingState_.setZero( );

    estimatable_parameters::EstimatebleParameterIdentifier undergoingBodyIdentifier =
        std::make_pair( estimatable_parameters::initial_body_state, std::make_pair( acceleratedBody, "" ) );
    if( customAccelerationPartialSet->customInitialStatePartials_.count( undergoingBodyIdentifier ) > 0 )
    {
        bodyUndergoingPositionPartial_ = customAccelerationPartialSet->customInitialStatePartials_.at( undergoingBodyIdentifier );
    }

    estimatable_parameters::EstimatebleParameterIdentifier exertingBodyIdentifier =
        std::make_pair( estimatable_parameters::initial_body_state, std::make_pair( acceleratingBody, "" ) );
    if( customAccelerationPartialSet->customInitialStatePartials_.count( exertingBodyIdentifier ) > 0 && ( exertingBodyIdentifier != undergoingBodyIdentifier ) )
    {
        bodyExertingPositionPartial_ = customAccelerationPartialSet->customInitialStatePartials_.at( exertingBodyIdentifier );
    }
}

void RadiationPressureAccelerationPartial::update( const double currentTime )
{
    radiationPressureAcceleration_->updateMembers( currentTime );

    if( !( currentTime_ == currentTime ) )
    {
        if( bodyUndergoingPositionPartial_ != nullptr )
        {
            currentPartialWrtUndergoingState_ = bodyUndergoingPositionPartial_->computePartial(
                currentTime, radiationPressureAcceleration_->getAcceleration( ), radiationPressureAcceleration_ );
        }

        if( bodyExertingPositionPartial_ != nullptr )
        {
            currentPartialWrtExertingState_ = bodyExertingPositionPartial_->computePartial(
                currentTime, radiationPressureAcceleration_->getAcceleration( ), radiationPressureAcceleration_ );
        }

        currentTime_ = currentTime;
    }
}

void RadiationPressureAccelerationPartial::wrtRadiationPressureCoefficient(
    Eigen::MatrixXd& partial, std::shared_ptr< electromagnetism::CannonballRadiationPressureTargetModel > targetModel )
{
    if( targetModel->getCoefficient( ) == 0.0 )
    {
        throw std::runtime_error( "Error in full radiation pressure partial w.r.t. Cr, partial is only implemented for non-zero coefficient" );
    }
    partial = radiationPressureAcceleration_->getAcceleration( ) / targetModel->getCoefficient( );

}

void RadiationPressureAccelerationPartial::wrtPanelledRadiationPressureCoefficient(
    Eigen::MatrixXd& partial, std::shared_ptr< electromagnetism::PaneledRadiationPressureTargetModel > targetModel )
{
    if( targetModel->getCoefficient( ) == 0.0 )
    {
        throw std::runtime_error( "Error in full radiation pressure partial w.r.t. Cr, partial is only implemented for non-zero coefficient" );
    }
    partial = radiationPressureAcceleration_->getAcceleration( ) / targetModel->getCoefficient( );

}

void RadiationPressureAccelerationPartial::wrtSpecularReflectivity(
        Eigen::MatrixXd& partial,
        std::shared_ptr< electromagnetism::PaneledRadiationPressureTargetModel > targetModel,
        const std::string& panelTypeId)
{
    std::function<double()> targetMassFunction = radiationPressureAcceleration_->getTargetMassFunction();
    double spacecraftMass = targetMassFunction();
    std::function<Eigen::Quaterniond()> targetRotationFromLocalToGlobalFrameFunction = radiationPressureAcceleration_->getTargetRotationFromLocalToGlobalFrameFunction();
    Eigen::Quaterniond targetRotationFromGlobalToLocalFrame = targetRotationFromLocalToGlobalFrameFunction().inverse( );
    std::function<Eigen::Vector3d()> sourceCenterPositionInGlobalFrameFunction = radiationPressureAcceleration_->getSourcePositionFunction();
    std::function<Eigen::Vector3d()> targetCenterPositionInGlobalFrameFunction = radiationPressureAcceleration_->getTargetPositionFunction();
    Eigen::Vector3d targetCenterPositionInGlobalFrame = targetCenterPositionInGlobalFrameFunction() - sourceCenterPositionInGlobalFrameFunction();
    Eigen::Vector3d sourceToTargetDirectionLocalFrame = targetRotationFromGlobalToLocalFrame * targetCenterPositionInGlobalFrame.normalized( );

    std::map< int, std::shared_ptr<system_models::VehicleExteriorPanel>> panelIndexMap = targetModel->getPanelIndexMap();
    std::vector< Eigen::Vector3d >& panelForces = targetModel->getPanelForces();
    std::vector< Eigen::Vector3d >& surfaceNormals = targetModel->getSurfaceNormals();
    for (auto it = panelIndexMap.begin(); it != panelIndexMap.end(); ++it)
    {
        // If panel is part of group, add the partial contribution
        if (it->second->getPanelTypeId() == panelTypeId)
        {
            // To get panel force without the reaction vector, divide by it. Then multiply by partial contribution. Divide by Sc mass to get acceleration
            Eigen::Vector3d panelForce = panelForces.at(it->first);
            Eigen::Vector3d surfaceNormal = surfaceNormals.at(it->first);
            Eigen::Vector3d reactionVector = it->second->getReflectionLaw()->evaluateReactionVector(surfaceNormal, sourceToTargetDirectionLocalFrame );

            Eigen::Vector3d reactionVectorPartialWrtSpecularReflectivity = it->second->getReflectionLaw()
                    ->evaluateReactionVectorPartialWrtSpecularReflectivity(surfaceNormal, sourceToTargetDirectionLocalFrame);

            if (reactionVectorPartialWrtSpecularReflectivity != Eigen::Vector3d::Zero() && reactionVector.norm() > 0)
            {
                partial += panelForce.cwiseQuotient(reactionVector).cwiseProduct(reactionVectorPartialWrtSpecularReflectivity) / spacecraftMass;
            }

        }
    }
}
void RadiationPressureAccelerationPartial::wrtDiffuseReflectivity(
        Eigen::MatrixXd& partial,
        std::shared_ptr< electromagnetism::PaneledRadiationPressureTargetModel > targetModel,
        const std::string& panelTypeId)
{
    std::function<double()> targetMassFunction = radiationPressureAcceleration_->getTargetMassFunction();
    double spacecraftMass = targetMassFunction();
    std::function<Eigen::Quaterniond()> targetRotationFromLocalToGlobalFrameFunction = radiationPressureAcceleration_->getTargetRotationFromLocalToGlobalFrameFunction();
    Eigen::Quaterniond targetRotationFromGlobalToLocalFrame = targetRotationFromLocalToGlobalFrameFunction().inverse( );
    std::function<Eigen::Vector3d()> sourceCenterPositionInGlobalFrameFunction = radiationPressureAcceleration_->getSourcePositionFunction();
    std::function<Eigen::Vector3d()> targetCenterPositionInGlobalFrameFunction = radiationPressureAcceleration_->getTargetPositionFunction();
    Eigen::Vector3d targetCenterPositionInGlobalFrame = targetCenterPositionInGlobalFrameFunction() - sourceCenterPositionInGlobalFrameFunction();
    Eigen::Vector3d sourceToTargetDirectionLocalFrame = targetRotationFromGlobalToLocalFrame * targetCenterPositionInGlobalFrame.normalized( );

    std::map< int, std::shared_ptr<system_models::VehicleExteriorPanel>> panelIndexMap = targetModel->getPanelIndexMap();
    std::vector< Eigen::Vector3d >& panelForces = targetModel->getPanelForces();
    std::vector< Eigen::Vector3d >& surfaceNormals = targetModel->getSurfaceNormals();
    for (auto it = panelIndexMap.begin(); it != panelIndexMap.end(); ++it)
    {
        // If panel is part of group, add the partial contribution
        if (it->second->getPanelTypeId() == panelTypeId)
        {
            // To get panel force without the reaction vector, divide by it. Then multiply by partial contribution. Divide by Sc mass to get acceleration
            Eigen::Vector3d panelForce = panelForces.at(it->first);
            Eigen::Vector3d surfaceNormal = surfaceNormals.at(it->first);
            Eigen::Vector3d reactionVector = it->second->getReflectionLaw()->evaluateReactionVector(surfaceNormal, sourceToTargetDirectionLocalFrame );

            Eigen::Vector3d reactionVectorPartialWrtDiffuseReflectivity = it->second->getReflectionLaw()
                    ->evaluateReactionVectorPartialWrtDiffuseReflectivity(surfaceNormal, sourceToTargetDirectionLocalFrame);

            if (reactionVectorPartialWrtDiffuseReflectivity != Eigen::Vector3d::Zero() && reactionVector.norm() > 0)
            {
                partial += panelForce.cwiseQuotient(reactionVector).cwiseProduct(reactionVectorPartialWrtDiffuseReflectivity) / spacecraftMass;

            }

        }
    }
}



}

}
