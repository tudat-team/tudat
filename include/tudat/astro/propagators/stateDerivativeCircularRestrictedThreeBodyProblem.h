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
 *        Wakker, K.F., "astro I, AE4-874", Delft University of Technology, 2007.
 *
 */

#ifndef TUDAT_STATE_DERIVATIVE_CIRCULAR_RESTRICTED_THREE_BODY_PROBLEM_H
#define TUDAT_STATE_DERIVATIVE_CIRCULAR_RESTRICTED_THREE_BODY_PROBLEM_H

#include <memory>

#include <Eigen/Core>

#include "tudat/basics/basicTypedefs.h"

namespace tudat
{
namespace propagators
{

Eigen::Vector6d computeCr3bpStateDerivative(
        const double time, const Eigen::Vector6d& cartesianState, const double massParameter );

Eigen::MatrixXd computeStateDerivativeWithStateTransitionMatrix(
        const double time, const Eigen::MatrixXd& cartesianState, const double massParameter );

} // namespace propagators

} // namespace tudat

#endif // TUDAT_STATE_DERIVATIVE_CIRCULAR_RESTRICTED_THREE_BODY_PROBLEM_H
