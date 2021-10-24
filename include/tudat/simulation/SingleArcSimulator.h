//! get_docstring

#ifndef TUDAT_SINGLEARCSIMULATOR_H
#define TUDAT_SINGLEARCSIMULATOR_H

#include <tudat/simulation/BaseSimulator.h>
#include <tudat/simulation/helpers.h>

#include <boost/make_shared.hpp>
#include <chrono>
#include <string>
#include <vector>

#include "tudat/astro/ephemerides/frameManager.h"
#include "tudat/astro/propagators/dynamicsStateDerivativeModel.h"
#include "tudat/astro/propagators/integrateEquations.h"
#include "tudat/astro/propagators/nBodyStateDerivative.h"
#include "tudat/basics/tudatTypeTraits.h"
#include "tudat/basics/utilities.h"
#include "tudat/math/interpolators/lagrangeInterpolator.h"
#include "tudat/simulation/propagation_setup/createEnvironmentUpdater.h"
#include "tudat/simulation/propagation_setup/createStateDerivativeModel.h"
#include "tudat/simulation/propagation_setup/propagationSettings.h"
#include "tudat/simulation/propagation_setup/propagationTermination.h"
#include "tudat/simulation/propagation_setup/setNumericallyIntegratedStates.h"

namespace tudat {

namespace numerical_simulation {

using namespace propagators;
using namespace numerical_integrators;
using namespace simulation_setup;
using namespace std::chrono;

//! @get_docstring(SingleArcSimulator)
template <typename StateScalarType = double, typename TimeType = double,
          typename TimeStepType = double>
class SingleArcSimulator : public BaseSimulator<StateScalarType, TimeType> {
 public:
  using BaseSimulator<StateScalarType, TimeType>::bodies_;
  using BaseSimulator<StateScalarType, TimeType>::clearNumericalSolutions_;
  using BaseSimulator<StateScalarType, TimeType>::setIntegratedResult_;

  using StateType = Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1>;
  using IntegratorType =
      std::shared_ptr<numerical_integrators::NumericalIntegrator<
          TimeType, StateType, StateType, TimeStepType>>;
  using StopPropagationFunctionType =
      std::function<bool(const double, const double)>;

  //! @get_docstring(SingleArcSimulator.ctor, 0)
  SingleArcSimulator(
      const simulation_setup::SystemOfBodies &bodies,
      const std::shared_ptr<numerical_integrators::IntegratorSettings<TimeType>>
          integratorSettings,
      const std::shared_ptr<PropagatorSettings<StateScalarType>>
          propagatorSettings,
      const std::vector<
          std::shared_ptr<SingleStateTypeDerivative<StateScalarType, TimeType>>>
          &stateDerivativeModels,
      const bool integrateToTermination = true,
      const bool clearNumericalSolutions = false,
      const bool setIntegratedResult = false,
      const bool printNumberOfFunctionEvaluations = false,
      const steady_clock::time_point initialClockTime = steady_clock::now(),
      const bool printDependentVariableData = true)
      : BaseSimulator<StateScalarType, TimeType>(
            bodies, clearNumericalSolutions, setIntegratedResult),
        integratorSettings_(integratorSettings),
        propagatorSettings_(std::dynamic_pointer_cast<
                            SingleArcPropagatorSettings<StateScalarType>>(
            propagatorSettings)),
        initialPropagationTime_(integratorSettings_->initialTime_),
        printNumberOfFunctionEvaluations_(printNumberOfFunctionEvaluations),
        initialClockTime_(initialClockTime),
        propagationTerminationReason_(
            std::make_shared<PropagationTerminationDetails>(
                propagation_never_run)),
        printDependentVariableData_(printDependentVariableData) {
    if (propagatorSettings == nullptr) {
      throw std::runtime_error(
          "Error in dynamics simulator, propagator settings not defined.");
    } else if (std::dynamic_pointer_cast<
                   SingleArcPropagatorSettings<StateScalarType>>(
                   propagatorSettings) == nullptr) {
      throw std::runtime_error(
          "Error in dynamics simulator, input must be single-arc.");
    }

    if (integratorSettings == nullptr) {
      throw std::runtime_error(
          "Error in dynamics simulator, integrator settings not defined.");
    }

    if (setIntegratedResult_) {
      frameManager_ = createFrameManager(bodies.getMap());
      integratedStateProcessors_ =
          createIntegratedStateProcessors<TimeType, StateScalarType>(
              propagatorSettings_, bodies_, frameManager_);
    }

    try {
      environmentUpdater_ =
          createEnvironmentUpdaterForDynamicalEquations<StateScalarType,
                                                        TimeType>(
              propagatorSettings_, bodies_);
    } catch (const std::runtime_error &error) {
      throw std::runtime_error("Error when creating environment updater: " +
                               std::string(error.what()));
    }

    if (stateDerivativeModels.size() == 0) {
      dynamicsStateDerivative_ = std::make_shared<
          DynamicsStateDerivativeModel<TimeType, StateScalarType>>(
          createStateDerivativeModels<StateScalarType, TimeType>(
              propagatorSettings_, bodies_, initialPropagationTime_),
          std::bind(
              &EnvironmentUpdater<StateScalarType, TimeType>::updateEnvironment,
              environmentUpdater_, std::placeholders::_1, std::placeholders::_2,
              std::placeholders::_3));
    } else {
      dynamicsStateDerivative_ = std::make_shared<
          DynamicsStateDerivativeModel<TimeType, StateScalarType>>(
          stateDerivativeModels,
          std::bind(
              &EnvironmentUpdater<StateScalarType, TimeType>::updateEnvironment,
              environmentUpdater_, std::placeholders::_1, std::placeholders::_2,
              std::placeholders::_3));
    }

    propagationTerminationCondition_ = createPropagationTerminationConditions(
        propagatorSettings_->getTerminationSettings(), bodies_,
        integratorSettings->initialTimeStep_,
        dynamicsStateDerivative_->getStateDerivativeModels());

    if (propagatorSettings_->getDependentVariablesToSave() != nullptr) {
      std::pair<std::function<Eigen::VectorXd()>, std::map<int, std::string>>
          dependentVariableData =
              createDependentVariableListFunction<TimeType, StateScalarType>(
                  propagatorSettings_->getDependentVariablesToSave(), bodies_,
                  dynamicsStateDerivative_->getStateDerivativeModels());
      dependentVariablesFunctions_ = dependentVariableData.first;
      dependentVariableIds_ = dependentVariableData.second;

      if (propagatorSettings_->getDependentVariablesToSave()
              ->printDependentVariableTypes_ &&
          printDependentVariableData_) {
        std::cout << "Dependent variables being saved, output vectors contain: "
                  << std::endl
                  << "Vector entry, Vector contents" << std::endl;
        utilities::printMapContents(dependentVariableIds_);
      }
    }

    stateDerivativeFunction_ = std::bind(
        &DynamicsStateDerivativeModel<TimeType,
                                      StateScalarType>::computeStateDerivative,
        dynamicsStateDerivative_, std::placeholders::_1, std::placeholders::_2);

    doubleStateDerivativeFunction_ = std::bind(
        &DynamicsStateDerivativeModel<
            TimeType, StateScalarType>::computeStateDoubleDerivative,
        dynamicsStateDerivative_, std::placeholders::_1, std::placeholders::_2);

    statePostProcessingFunction_ = std::bind(
        &DynamicsStateDerivativeModel<TimeType,
                                      StateScalarType>::postProcessState,
        dynamicsStateDerivative_, std::placeholders::_1);

    stopPropagationFunction_ =
        std::bind(&PropagationTerminationCondition::checkStopCondition,
                  propagationTerminationCondition_, std::placeholders::_1,
                  std::placeholders::_2);

    reset();
    // Integrate equations of motion if required.
    if (integrateToTermination) {
      this->integrateToTermination();
    }
  }

  SingleArcSimulator(
      const simulation_setup::SystemOfBodies &bodies,
      const std::shared_ptr<numerical_integrators::IntegratorSettings<TimeType>>
          integratorSettings,
      const std::shared_ptr<PropagatorSettings<StateScalarType>>
          propagatorSettings,
      const bool integrateToTermination = true,
      const bool clearNumericalSolutions = false,
      const bool setIntegratedResult = false,
      const bool printDependentVariableData = true)
      : SingleArcSimulator(
            bodies, integratorSettings, propagatorSettings,
            std::vector<std::shared_ptr<
                SingleStateTypeDerivative<StateScalarType, TimeType>>>(),
            integrateToTermination, clearNumericalSolutions,
            setIntegratedResult, printDependentVariableData) {}

  //! Destructor
  ~SingleArcSimulator() = default;

  //! @get_docstring(SingleArcSimulator.reset)
  void reset() {
    // Empty solution maps
    equationsOfMotionNumericalSolution_.clear();
    equationsOfMotionNumericalSolutionRaw_.clear();

    // Reset functions
    dynamicsStateDerivative_->setPropagationSettings(
        std::vector<IntegratedStateType>(), 1, 0);
    dynamicsStateDerivative_->resetFunctionEvaluationCounter();
    dynamicsStateDerivative_->resetCumulativeFunctionEvaluationCounter();

    // Reset initial time to ensure consistency with multi-arc propagation.
    integratorSettings_->initialTime_ = this->initialPropagationTime_;

    // Integrate equations of motion numerically.
    resetPropagationTerminationConditions();
    setAreBodiesInPropagation(bodies_, true);

    //                stateDerivativeFunction_ =
    //                dynamicsStateDerivative_->convertFromOutputSolution(
    //                        propagatorSettings_->getInitialStates(),
    //                        this->initialPropagationTime_ );

    integrator_ = createIntegrator<TimeType, StateType, TimeStepType>(
        stateDerivativeFunction_, propagatorSettings_->getInitialStates(),
        integratorSettings_);

    if (integratorSettings_->assessTerminationOnMinorSteps_) {
      integrator_->setPropagationTerminationFunction(stopPropagationFunction_);
    }

    integrationInterface_ =
        std::make_shared<IntegrationInterface<StateScalarType, StateType,
                                              TimeType, TimeStepType>>(
            integrator_, integratorSettings_->initialTimeStep_,
            propagationTerminationCondition_, dependentVariablesFunctions_,
            equationsOfMotionNumericalSolutionRaw_, dependentVariableHistory_,
            cumulativeComputationTimeHistory_, statePostProcessingFunction_,
            integratorSettings_->saveFrequency_);
  }

  bool isTerminal() { return integrationInterface_->isTerminal(); }

  void integrateToTermination() {
    reset();
    propagationTerminationReason_ =
        EquationIntegrationInterface<StateType, TimeType>::integrateEquations(
            stateDerivativeFunction_, equationsOfMotionNumericalSolutionRaw_,
            dynamicsStateDerivative_->convertFromOutputSolution(
                propagatorSettings_->getInitialStates(),
                this->initialPropagationTime_),
            integratorSettings_, propagationTerminationCondition_,
            dependentVariableHistory_, cumulativeComputationTimeHistory_,
            dependentVariablesFunctions_, statePostProcessingFunction_,
            propagatorSettings_->getPrintInterval(), initialClockTime_);
    simulation_setup::setAreBodiesInPropagation(bodies_, false);

    // Convert numerical solution to conventional state
    dynamicsStateDerivative_->convertNumericalStateSolutionsToOutputSolutions(
        equationsOfMotionNumericalSolution_,
        equationsOfMotionNumericalSolutionRaw_);

    // Retrieve number of cumulative function evaluations
    cumulativeNumberOfFunctionEvaluations_ =
        dynamicsStateDerivative_->getCumulativeNumberOfFunctionEvaluations();

    // Retrieve and print number of total function evaluations
    if (printNumberOfFunctionEvaluations_) {
      std::cout << "Total Number of Function Evaluations: "
                << dynamicsStateDerivative_->getNumberOfFunctionEvaluations()
                << std::endl;
    }

    if (this->setIntegratedResult_) {
      processNumericalEquationsOfMotionSolution();
    }
  }

  void integrateByStep(const int steps) { integrationInterface_->step(steps); }

  //! @get_docstring(SingleArcSimulator.integrateEquationsOfMotion)
  [[deprecated]] void integrateEquationsOfMotion() { integrateToTermination(); }

  //!
  const std::map<TimeType, StateType> &getEquationsOfMotionNumericalSolution() {
    // Convert numerical solution to conventional state
    dynamicsStateDerivative_->convertNumericalStateSolutionsToOutputSolutions(
        equationsOfMotionNumericalSolution_,
        equationsOfMotionNumericalSolutionRaw_);
    return equationsOfMotionNumericalSolution_;
  }

  //!
  const std::map<TimeType, StateType>
      &getEquationsOfMotionNumericalSolutionRaw() {
    return equationsOfMotionNumericalSolutionRaw_;
  }

  //! Function to return the map of dependent variable history that was saved
  //! during numerical propagation.
  /*!
   * Function to return the map of dependent variable history that was saved
   * during numerical propagation. \return Map of dependent variable history
   * that was saved during numerical propagation.
   */
  const std::map<TimeType, Eigen::VectorXd> &getDependentVariableHistory() {
    return dependentVariableHistory_;
  }

  //! Function to return the map of cumulative computation time history that was
  //! saved during numerical propagation.
  /*!
   * Function to return the map of cumulative computation time history that was
   * saved during numerical propagation. \return Map of cumulative computation
   * time history that was saved during numerical propagation.
   */
  std::map<TimeType, double> getCumulativeComputationTimeHistory() {
    return cumulativeComputationTimeHistory_;
  }

  //! Function to return the map of number of cumulative function evaluations
  //! that was saved during numerical propagation.
  /*!
   * Function to return the map of cumulative number of function evaluations
   * that was saved during numerical propagation. \return Map of cumulative
   * number of function evaluations that was saved during numerical propagation.
   */
  std::map<TimeType, unsigned int> getCumulativeNumberOfFunctionEvaluations() {
    return cumulativeNumberOfFunctionEvaluations_;
  }

  //! Function to return the map of state history of numerically integrated
  //! bodies (base class interface).
  /*!
   * Function to return the map of state history of numerically integrated
   * bodies (base class interface). \return Vector is size 1, with entry: map of
   * state history of numerically integrated bodies.
   */
  std::vector<
      std::map<TimeType, Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1>>>
  getEquationsOfMotionNumericalSolutionBase() {
    return std::vector<
        std::map<TimeType, Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1>>>(
        {getEquationsOfMotionNumericalSolution()});
  }

  //! Function to return the map of dependent variable history that was saved
  //! during numerical propagation (base class interface)
  /*!
   * Function to return the map of dependent variable history that was saved
   * during numerical propagation (base class interface) \return Vector is size
   * 1, with entry: map of dependent variable history that was saved during
   * numerical propagation.
   */
  std::vector<std::map<TimeType, Eigen::VectorXd>>
  getDependentVariableNumericalSolutionBase() {
    return std::vector<std::map<TimeType, Eigen::VectorXd>>(
        {getDependentVariableHistory()});
  }

  //! Function to return the map of cumulative computation time history that was
  //! saved during numerical propagation.
  /*!
   * Function to return the map of cumulative computation time history that was
   * saved during numerical propagation (base class interface). \return Vector
   * is size 1, with entry: map of cumulative computation time history that was
   * saved during numerical propagation.
   */
  std::vector<std::map<TimeType, double>>
  getCumulativeComputationTimeHistoryBase() {
    return std::vector<std::map<TimeType, double>>(
        {getCumulativeComputationTimeHistory()});
  }

  //! Function to reset the environment from an externally generated state
  //! history.
  /*!
   * Function to reset the environment from an externally generated state
   * history, the order of the entries in the state vectors are proscribed by
   * propagatorSettings \param equationsOfMotionNumericalSolution Externally
   * generated state history. \param processSolution True if the new solution is
   * to be immediately processed (default true). \param dependentVariableHistory
   * Externally generated dependent variable history.
   */
  void manuallySetAndProcessRawNumericalEquationsOfMotionSolution(
      const std::map<TimeType,
                     Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1>>
          &equationsOfMotionNumericalSolution,
      const std::map<TimeType, Eigen::VectorXd> &dependentVariableHistory,
      const bool processSolution = true) {
    equationsOfMotionNumericalSolution_ = equationsOfMotionNumericalSolution;
    if (processSolution) {
      processNumericalEquationsOfMotionSolution();
    }

    dependentVariableHistory_ = dependentVariableHistory;
  }

  //! Function to get the settings for the numerical integrator.
  /*!
   * Function to get the settings for the numerical integrator.
   * \return The settings for the numerical integrator.
   */
  std::shared_ptr<numerical_integrators::IntegratorSettings<TimeType>>
  getIntegratorSettings() {
    return integratorSettings_;
  }

  //! Function to get the function that performs a single state derivative
  //! function evaluation.
  /*!
   * Function to get the function that performs a single state derivative
   * function evaluation. \return Function that performs a single state
   * derivative function evaluation.
   */
  std::function<Eigen::Matrix<StateScalarType, Eigen::Dynamic, Eigen::Dynamic>(
      const TimeType,
      const Eigen::Matrix<StateScalarType, Eigen::Dynamic, Eigen::Dynamic> &)>
  getStateDerivativeFunction() {
    return stateDerivativeFunction_;
  }

  //! Function to get the function that performs a single state derivative
  //! function evaluation with double precision.
  /*!
   * Function to get the function that performs a single state derivative
   * function evaluation with double precision, regardless of template
   * arguments. \return Function that performs a single state derivative
   * function evaluation with double precision.
   */
  std::function<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(
      const double,
      const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &)>
  getDoubleStateDerivativeFunction() {
    return doubleStateDerivativeFunction_;
  }

  //! Function to get the settings for the propagator.
  /*!
   * Function to get the settings for the propagator.
   * \return The settings for the propagator.
   */
  std::shared_ptr<SingleArcPropagatorSettings<StateScalarType>>

  getPropagatorSettings() {
    return propagatorSettings_;
  }

  //! Function to get the object that updates the environment.
  /*!
   * Function to get the object responsible for updating the environment based
   * on the current state and time. \return Object responsible for updating the
   * environment based on the current state and time.
   */
  std::shared_ptr<EnvironmentUpdater<StateScalarType, TimeType>>

  getEnvironmentUpdater() {
    return environmentUpdater_;
  }

  //! Function to get the object that updates and returns state derivative
  /*!
   * Function to get the object that updates current environment and returns
   * state derivative from single function call \return Object that updates
   * current environment and returns state derivative from single function call
   */
  std::shared_ptr<DynamicsStateDerivativeModel<TimeType, StateScalarType>>

  getDynamicsStateDerivative() {
    return dynamicsStateDerivative_;
  }

  //! Function to retrieve the object defining when the propagation is to be
  //! terminated.
  /*!
   * Function to retrieve the object defining when the propagation is to be
   * terminated. \return Object defining when the propagation is to be
   * terminated.
   */
  std::shared_ptr<PropagationTerminationCondition>
  getPropagationTerminationCondition() {
    return propagationTerminationCondition_;
  }

  //! Function to retrieve the list of object that process the integrated
  //! numerical solution by updating the environment
  /*!
   * Function to retrieve the List of object (per dynamics type) that process
   * the integrated numerical solution by updating the environment \return List
   * of object (per dynamics type) that process the integrated numerical
   * solution by updating the environment
   */
  std::map<IntegratedStateType,
           std::vector<std::shared_ptr<
               IntegratedStateProcessor<TimeType, StateScalarType>>>>

  getIntegratedStateProcessors() {
    return integratedStateProcessors_;
  }

  //! Function to retrieve the event that triggered the termination of the last
  //! propagation
  /*!
   * Function to retrieve the event that triggered the termination of the last
   * propagation \return Event that triggered the termination of the last
   * propagation
   */
  std::shared_ptr<PropagationTerminationDetails>
  getPropagationTerminationReason() {
    return propagationTerminationReason_;
  }

  //! Get whether the integration was completed successfully.
  /*!
   * Get whether the integration was completed successfully.
   * \return Whether the integration was completed successfully by reaching the
   * termination condition.
   */
  virtual bool integrationCompletedSuccessfully() const {
    return (propagationTerminationReason_->getPropagationTerminationReason() ==
            termination_condition_reached);
  }

  //!
  std::map<int, std::string> getDependentVariableIds() {
    return dependentVariableIds_;
  }

  //!
  double getInitialPropagationTime() { return this->initialPropagationTime_; }

  //!
  void resetInitialPropagationTime(const double initialPropagationTime) {
    initialPropagationTime_ = initialPropagationTime;
  }

  //!
  std::function<Eigen::VectorXd()> getDependentVariablesFunctions() {
    return dependentVariablesFunctions_;
  }

  //!
  void resetPropagationTerminationConditions() {
    propagationTerminationCondition_ = createPropagationTerminationConditions(
        propagatorSettings_->getTerminationSettings(), bodies_,
        integratorSettings_->initialTimeStep_,
        dynamicsStateDerivative_->getStateDerivativeModels());
  }

  //!
  void processNumericalEquationsOfMotionSolution() {
    // Create and set interpolators for ephemerides
    resetIntegratedStates(equationsOfMotionNumericalSolution_,
                          integratedStateProcessors_);

    // Clear numerical solution if so required.
    if (clearNumericalSolutions_) {
      equationsOfMotionNumericalSolution_.clear();
      equationsOfMotionNumericalSolutionRaw_.clear();
    }

    for (auto bodyIterator : bodies_.getMap()) {
      bodyIterator.second->updateConstantEphemerisDependentMemberQuantities();
    }
  }

  void suppressDependentVariableDataPrinting() {
    printDependentVariableData_ = false;
  }

  void enableDependentVariableDataPrinting() {
    printDependentVariableData_ = true;
  }

 protected:
  StopPropagationFunctionType stopPropagationFunction_;

  std::shared_ptr<
      IntegrationInterface<StateScalarType, StateType, TimeType, TimeStepType>>
      integrationInterface_;

  IntegratorType integrator_;

  //! List of object (per dynamics type) that process the integrated numerical
  //! solution by updating the environment
  std::map<IntegratedStateType,
           std::vector<std::shared_ptr<
               IntegratedStateProcessor<TimeType, StateScalarType>>>>
      integratedStateProcessors_;

  //! Object responsible for updating the environment based on the current state
  //! and time.
  /*!
   *  Object responsible for updating the environment based on the current state
   * and time. Calling the updateEnvironment function automatically updates all
   * dependent variables that are needed to calulate the state derivative.
   */
  std::shared_ptr<EnvironmentUpdater<StateScalarType, TimeType>>
      environmentUpdater_;

  //! Interface object that updates current environment and returns state
  //! derivative from single function call.
  std::shared_ptr<DynamicsStateDerivativeModel<TimeType, StateScalarType>>
      dynamicsStateDerivative_;

  //! Function that performs a single state derivative function evaluation.
  /*!
   *  Function that performs a single state derivative function evaluation, will
   * typically be set to DynamicsStateDerivativeModel< TimeType, StateScalarType
   * >::computeStateDerivative function. Calling this function will first update
   * the environment (using environmentUpdater_) and then calculate the full
   * system state derivative.
   */
  std::function<Eigen::Matrix<StateScalarType, Eigen::Dynamic, Eigen::Dynamic>(
      const TimeType,
      const Eigen::Matrix<StateScalarType, Eigen::Dynamic, Eigen::Dynamic> &)>
      stateDerivativeFunction_;

  //! Function that performs a single state derivative function evaluation with
  //! double precision.
  /*!
   *  Function that performs a single state derivative function evaluation with
   * double precision. \sa stateDerivativeFunction_
   */
  std::function<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(
      const double,
      const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &)>
      doubleStateDerivativeFunction_;

  //! Settings for numerical integrator.
  std::shared_ptr<numerical_integrators::IntegratorSettings<TimeType>>
      integratorSettings_;

  //! Settings for propagator.
  std::shared_ptr<SingleArcPropagatorSettings<StateScalarType>>
      propagatorSettings_;

  //! Object defining when the propagation is to be terminated.
  std::shared_ptr<PropagationTerminationCondition>
      propagationTerminationCondition_;

  //! Function returning dependent variables (during numerical propagation)
  std::function<Eigen::VectorXd()> dependentVariablesFunctions_;

  //! Function to post-process state (during numerical propagation)
  std::function<void(Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1> &)>
      statePostProcessingFunction_;

  //! Map listing starting entry of dependent variables in output vector, along
  //! with associated ID.
  std::map<int, std::string> dependentVariableIds_;

  //! Object for retrieving ephemerides for transformation of reference frame
  //! (origins)
  std::shared_ptr<ephemerides::ReferenceFrameManager> frameManager_;

  //! Map of state history of numerically integrated bodies.
  /*!
   *  Map of state history of numerically integrated bodies, i.e. the result of
   * the numerical integration, transformed into the 'conventional form' (\sa
   * SingleStateTypeDerivative::convertToOutputSolution). Key of map denotes
   * time, values are concatenated vectors of integrated body states (order
   * defined by propagatorSettings_). NOTE: this map is empty if
   * clearNumericalSolutions_ is set to true.
   */
  std::map<TimeType, Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1>>
      equationsOfMotionNumericalSolution_;

  //! Map of state history of numerically integrated bodies.
  /*!
   *  Map of state history of numerically integrated bodies, i.e. the result of
   * the numerical integration, in the original propagation coordinates. Key of
   * map denotes time, values are concatenated vectors of integrated body states
   * (order defined by propagatorSettings_). NOTE: this map is empty if
   * clearNumericalSolutions_ is set to true.
   */
  std::map<TimeType, Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1>>
      equationsOfMotionNumericalSolutionRaw_;

  //! Map of dependent variable history that was saved during numerical
  //! propagation.
  std::map<TimeType, Eigen::VectorXd> dependentVariableHistory_;

  //! Map of cumulative computation time history that was saved during numerical
  //! propagation.
  std::map<TimeType, double> cumulativeComputationTimeHistory_;

  //! Map of cumulative number of function evaluations that was saved during
  //! numerical propagation.
  std::map<TimeType, unsigned int> cumulativeNumberOfFunctionEvaluations_;

  //! Initial time of propagation
  double initialPropagationTime_;

  //! Boolean denoting whether the number of function evaluations should be
  //! printed at the end of propagation.
  bool printNumberOfFunctionEvaluations_;

  //! Initial clock time
  steady_clock::time_point initialClockTime_;

  //! Event that triggered the termination of the propagation
  std::shared_ptr<PropagationTerminationDetails> propagationTerminationReason_;

  bool printDependentVariableData_;
};

}  // namespace numerical_simulation

}  // namespace tudat

#endif  // TUDAT_SINGLEARCSIMULATOR_H
