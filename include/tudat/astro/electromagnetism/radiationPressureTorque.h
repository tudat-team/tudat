/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_RADIATIONPRESSURETORQUE_H
#define TUDAT_RADIATIONPRESSURETORQUE_H

#include <functional>
#include <memory>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "tudat/astro/basic_astro/torqueModel.h"
#include "tudat/astro/electromagnetism/radiationPressureAcceleration.h"


namespace tudat
{
namespace electromagnetism
{

/*!
 * Class modeling radiation pressure torque. Radiation pressure accelerates a target due to electromagnetic
 * radiation from a source.
 */
class RadiationPressureTorque: public basic_astrodynamics::TorqueModel
{
public:

    RadiationPressureTorque(
        const std::shared_ptr< RadiationPressureAcceleration > radiationPressureAcceleration,
        const std::function< Eigen::Vector3d( ) > centerOfMassFunction ):
        radiationPressureAcceleration_( radiationPressureAcceleration ),
        centerOfMassFunction_( centerOfMassFunction )
    {
        radiationPressureAcceleration_->getTargetModel( )->enableTorqueComputation( centerOfMassFunction );
    }

    ~RadiationPressureTorque( ){ }

    /*!
     * Update class members.
     *
     * @param currentTime Current simulation time
     */
    void updateMembers(double currentTime) override;

    void resetAccelerationModel( const std::shared_ptr< RadiationPressureAcceleration > radiationPressureAcceleration )
    {
        radiationPressureAcceleration_ = radiationPressureAcceleration;
    }

    Eigen::Vector3d getTorque( ) override
    {
        return currentTorque_;
    }

protected:

    std::shared_ptr< RadiationPressureAcceleration > radiationPressureAcceleration_;

    const std::function< Eigen::Vector3d( ) > centerOfMassFunction_;

    Eigen::Vector3d currentTorque_;
};

} // tudat
} // electromagnetism

#endif //TUDAT_RADIATIONPRESSURETORQUE_H
