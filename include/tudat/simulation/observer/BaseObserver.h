//
// Created by ggarr on 23/10/2021.
//

#ifndef TUDAT_BASEOBSERVER_H
#define TUDAT_BASEOBSERVER_H

/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <tudat/simulation/BaseSimulator.h>
#include <tudat/simulation/helpers.h>

#include <boost/make_shared.hpp>
#include <chrono>
#include <memory>
#include <string>
#include <vector>

#include "tudat/astro/ephemerides/frameManager.h"
#include "tudat/astro/propagators/dynamicsStateDerivativeModel.h"
#include "tudat/astro/propagators/integrateEquations.h"
#include "tudat/astro/propagators/nBodyStateDerivative.h"
#include "tudat/basics/tudatTypeTraits.h"
#include "tudat/basics/utilities.h"
#include "tudat/math/interpolators/lagrangeInterpolator.h"
#include "tudat/simulation/BaseSimulator.h"
#include "tudat/simulation/propagation_setup/createEnvironmentUpdater.h"
#include "tudat/simulation/propagation_setup/createStateDerivativeModel.h"
#include "tudat/simulation/propagation_setup/propagationSettings.h"
#include "tudat/simulation/propagation_setup/propagationTermination.h"
#include "tudat/simulation/propagation_setup/setNumericallyIntegratedStates.h"

namespace tudat {

namespace numerical_simulation {


template <typename StateScalarType, typename TimeType, typename TimeStepType,
          typename ObservedType>
class BaseObserver {
 public:
  BaseObserver(){};

  virtual void update(BaseSimulator<StateScalarType, TimeType> simulation) {};
};

}  // namespace numerical_simulation

}  // namespace tudat

#endif  // TUDAT_BASEOBSERVER_H
