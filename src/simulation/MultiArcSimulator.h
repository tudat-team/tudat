//
// Created by ggarr on 18/10/2021.
//

#ifndef TUDATBUNDLE_MULTIARCSIMULATOR_H
#define TUDATBUNDLE_MULTIARCSIMULATOR_H

//! Class for performing full numerical integration of a dynamical system over
//! multiple arcs.
/*!
 *  Class for performing full numerical integration of a dynamical system over
 * multiple arcs, equations of motion are set up for each arc (and need not be
 * equal for each arc). In this class, the governing equations are set once, but
 * can be re-integrated for different initial conditions using the same instance
 * of the class.
 */
template <typename StateScalarType = double, typename TimeType = double>
class MultiArcDynamicsSimulator
    : public BaseSimulator<StateScalarType, TimeType> {
 public:
  using BaseSimulator<StateScalarType, TimeType>::bodies_;
  using BaseSimulator<StateScalarType, TimeType>::clearNumericalSolutions_;

  //! Constructor of multi-arc simulator for same integration settings per arc.
  /*!
   *  Constructor of multi-arc simulator for same integration settings per arc.
   *  \param bodies Map of bodies (with names) of all bodies in integration.
   *  \param integratorSettings Integrator settings for numerical integrator,
   * used for all arcs. \param propagatorSettings Propagator settings for
   * dynamics (must be of multi arc type) \param arcStartTimes Times at which
   * the separate arcs start \param areEquationsOfMotionToBeIntegrated Boolean
   * to denote whether equations of motion should be integrated at the end of
   * the contructor or not. \param clearNumericalSolutions Boolean to determine
   * whether to clear the raw numerical solution member variables after
   * propagation and resetting ephemerides (default true). \param
   * setIntegratedResult Boolean to determine whether to automatically use the
   * integrated results to set ephemerides (default true).
   */
  MultiArcDynamicsSimulator(
      const simulation_setup::SystemOfBodies &bodies,
      const std::shared_ptr<
          numerical_integrators::IntegratorSettings<TimeType> >
          integratorSettings,
      const std::shared_ptr<PropagatorSettings<StateScalarType> >
          propagatorSettings,
      const std::vector<double> arcStartTimes,
      const bool areEquationsOfMotionToBeIntegrated = true,
      const bool clearNumericalSolutions = true,
      const bool setIntegratedResult = true)
      : BaseSimulator<StateScalarType, TimeType>(
            bodies, clearNumericalSolutions, setIntegratedResult) {
    multiArcPropagatorSettings_ =
        std::dynamic_pointer_cast<MultiArcPropagatorSettings<StateScalarType> >(
            propagatorSettings);
    if (multiArcPropagatorSettings_ == nullptr) {
      throw std::runtime_error(
          "Error when creating multi-arc dynamics simulator, input is not "
          "multi arc");
    } else {
      std::vector<
          std::shared_ptr<SingleArcPropagatorSettings<StateScalarType> > >
          singleArcSettings =
              multiArcPropagatorSettings_->getSingleArcSettings();

      arcStartTimes_.resize(arcStartTimes.size());

      if (singleArcSettings.size() != arcStartTimes.size()) {
        throw std::runtime_error(
            "Error when creating multi-arc dynamics simulator, input is "
            "inconsistent");
      }
      // Create dynamics simulators
      for (unsigned int i = 0; i < singleArcSettings.size(); i++) {
        integratorSettings->initialTime_ = arcStartTimes.at(i);

        singleArcDynamicsSimulators_.push_back(
            std::make_shared<SingleArcSimulator<StateScalarType, TimeType> >(
                bodies, integratorSettings, singleArcSettings.at(i), false,
                false, true));
        singleArcDynamicsSimulators_[i]->resetSetIntegratedResult(false);
      }

      equationsOfMotionNumericalSolution_.resize(arcStartTimes.size());
      dependentVariableHistory_.resize(arcStartTimes.size());
      cumulativeComputationTimeHistory_.resize(arcStartTimes.size());
      propagationTerminationReasons_.resize(arcStartTimes.size());

      // Integrate equations of motion if required.
      if (areEquationsOfMotionToBeIntegrated) {
        integrateEquationsOfMotion(
            multiArcPropagatorSettings_->getInitialStates());
      }
    }
  }

  //! Constructor of multi-arc simulator for different integration settings per
  //! arc.
  /*!
   *  Constructor of multi-arc simulator for different integration settings per
   * arc. \param bodies Map of bodies (with names) of all bodies in integration.
   *  \param integratorSettings List of integrator settings for numerical
   * integrator, defined per arc. \param propagatorSettings Propagator settings
   * for dynamics (must be of multi arc type) \param
   * areEquationsOfMotionToBeIntegrated Boolean to denote whether equations of
   * motion should be integrated at the end of the contructor or not. \param
   * clearNumericalSolutions Boolean to determine whether to clear the raw
   * numerical solution member variables after propagation and resetting
   * ephemerides (default true). \param setIntegratedResult Boolean to determine
   * whether to automatically use the integrated results to set ephemerides
   * (default true).
   */
  MultiArcDynamicsSimulator(
      const simulation_setup::SystemOfBodies &bodies,
      const std::vector<std::shared_ptr<
          numerical_integrators::IntegratorSettings<TimeType> > >
          integratorSettings,
      const std::shared_ptr<PropagatorSettings<StateScalarType> >
          propagatorSettings,
      const bool areEquationsOfMotionToBeIntegrated = true,
      const bool clearNumericalSolutions = true,
      const bool setIntegratedResult = true)
      : BaseSimulator<StateScalarType, TimeType>(
            bodies, clearNumericalSolutions, setIntegratedResult) {
    multiArcPropagatorSettings_ =
        std::dynamic_pointer_cast<MultiArcPropagatorSettings<StateScalarType> >(
            propagatorSettings);
    if (multiArcPropagatorSettings_ == nullptr) {
      throw std::runtime_error(
          "Error when creating multi-arc dynamics simulator, input is not "
          "multi arc");
    } else {
      std::vector<
          std::shared_ptr<SingleArcPropagatorSettings<StateScalarType> > >
          singleArcSettings =
              multiArcPropagatorSettings_->getSingleArcSettings();

      if (singleArcSettings.size() != integratorSettings.size()) {
        throw std::runtime_error(
            "Error when creating multi-arc dynamics simulator, input sizes are "
            "inconsistent");
      }

      arcStartTimes_.resize(singleArcSettings.size());

      // Create dynamics simulators
      for (unsigned int i = 0; i < singleArcSettings.size(); i++) {
        singleArcDynamicsSimulators_.push_back(
            std::make_shared<SingleArcSimulator<StateScalarType, TimeType> >(
                bodies, integratorSettings.at(i), singleArcSettings.at(i),
                false, false, true));
        singleArcDynamicsSimulators_[i]->resetSetIntegratedResult(false);
      }

      equationsOfMotionNumericalSolution_.resize(singleArcSettings.size());
      dependentVariableHistory_.resize(singleArcSettings.size());
      cumulativeComputationTimeHistory_.resize(singleArcSettings.size());
      propagationTerminationReasons_.resize(singleArcSettings.size());

      // Integrate equations of motion if required.
      if (areEquationsOfMotionToBeIntegrated) {
        integrateEquationsOfMotion(
            multiArcPropagatorSettings_->getInitialStates());
      }
    }
  }

  //! Destructor
  ~MultiArcSimulator() {}

  //! This function numerically (re-)integrates the equations of motion, using
  //! concatenated states for all arcs
  /*!
   *  This function numerically (re-)integrates the equations of motion, using
   * the settings set through the constructor and a new initial state vector
   * provided here. The raw results are set in the
   * equationsOfMotionNumericalSolution_ \param concatenatedInitialStates
   * Initial state vector that is to be used for numerical integration. Note
   * that this state should be in the correct frame (i.e. corresponding to
   * centralBodies in propagatorSettings_), but not in the propagator- specific
   * form (i.e Encke, Gauss, etc. for translational dynamics). The states for
   * all arcs must be concatenated in order into a single Eigen Vector.
   */
  void integrateEquationsOfMotion(
      const Eigen::Matrix<StateScalarType, Eigen::Dynamic, Eigen::Dynamic>
          &concatenatedInitialStates) {
    std::vector<Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1> >
        splitInitialState;

    int currentIndex = 0;
    for (unsigned int i = 0; i < singleArcDynamicsSimulators_.size(); i++) {
      int currentSize = singleArcDynamicsSimulators_.at(i)
                            ->getPropagatorSettings()
                            ->getConventionalStateSize();
      splitInitialState.push_back(
          concatenatedInitialStates.block(currentIndex, 0, currentSize, 1));
      currentIndex += currentSize;
    }

    if (currentIndex != concatenatedInitialStates.rows()) {
      throw std::runtime_error(
          "Error when doing multi-arc integration, input state vector size is "
          "incompatible with settings");
    }

    integrateEquationsOfMotion(splitInitialState);
  }

  //! This function numerically (re-)integrates the equations of motion, using
  //! separate states for all arcs
  /*!
   *  This function numerically (re-)integrates the equations of motion, using
   * the settings set through the constructor and a new initial state vector
   * provided here. The raw results are set in the
   * equationsOfMotionNumericalSolution_ \param initialStatesList Initial state
   * vector that is to be used for numerical integration. Note that this state
   * should be in the correct frame (i.e. corresponding to centralBodies in
   * propagatorSettings_), but not in the propagator- specific form (i.e Encke,
   * Gauss, etc. for translational dynamics). The states for all stored, in
   * order, in the input std vector.
   */
  void integrateEquationsOfMotion(
      const std::vector<Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1> >
          &initialStatesList) {
    // Clear existing solution (if any)
    for (unsigned int i = 0; i < equationsOfMotionNumericalSolution_.size();
         i++) {
      equationsOfMotionNumericalSolution_.at(i).clear();
    }

    for (unsigned int i = 0; i < dependentVariableHistory_.size(); i++) {
      dependentVariableHistory_.at(i).clear();
    }

    for (unsigned int i = 0; i < cumulativeComputationTimeHistory_.size();
         i++) {
      cumulativeComputationTimeHistory_.at(i).clear();
    }

    Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1> currentArcInitialState;
    std::vector<Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1> >
        arcInitialStateList;
    bool updateInitialStates = false;

    // Propagate dynamics for each arc
    for (unsigned int i = 0; i < singleArcDynamicsSimulators_.size(); i++) {
      // Get arc initial state. If initial state is NaN, this signals that the
      // initial state is to be taken from previous arc
      if ((i == 0) || (!linear_algebra::doesMatrixHaveNanEntries(
                          initialStatesList.at(i)))) {
        currentArcInitialState = initialStatesList.at(i);
      } else {
        currentArcInitialState = getArcInitialStateFromPreviousArcResult(
            equationsOfMotionNumericalSolution_.at(i - 1),
            singleArcDynamicsSimulators_.at(i)->getInitialPropagationTime());

        // If arc initial state is taken from previous arc, this indicates that
        // the initial states in propagator settings need to be updated.
        updateInitialStates = true;
      }
      arcInitialStateList.push_back(currentArcInitialState);

      singleArcDynamicsSimulators_.at(i)->integrateEquationsOfMotion(
          currentArcInitialState);
      equationsOfMotionNumericalSolution_[i] =
          std::move(singleArcDynamicsSimulators_.at(i)
                        ->getEquationsOfMotionNumericalSolution());
      dependentVariableHistory_[i] = std::move(
          singleArcDynamicsSimulators_.at(i)->getDependentVariableHistory());
      cumulativeComputationTimeHistory_[i] =
          std::move(singleArcDynamicsSimulators_.at(i)
                        ->getCumulativeComputationTimeHistory());
      propagationTerminationReasons_[i] =
          singleArcDynamicsSimulators_.at(i)->getPropagationTerminationReason();
      arcStartTimes_[i] = equationsOfMotionNumericalSolution_[i].begin()->first;
    }

    if (updateInitialStates) {
      multiArcPropagatorSettings_->resetInitialStatesList(arcInitialStateList);
    }

    if (this->setIntegratedResult_) {
      processNumericalEquationsOfMotionSolution();
    }
  }

  //! Function to return the numerical solution to the equations of motion.
  /*!
   *  Function to return the numerical solution to the equations of motion for
   * last numerical integration. Each vector entry denotes one arc. Key of map
   * denotes time, values are full propagated state vectors. \return List of
   * maps of history of numerically integrated states.
   */
  std::vector<
      std::map<TimeType, Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1> > >
  getEquationsOfMotionNumericalSolution() {
    return equationsOfMotionNumericalSolution_;
  }

  //! Function to return the numerical solution of the dependent variables
  /*!
   *  Function to return the numerical solution of the dependent variables for
   * last numerical integration. Each vector entry denotes one arc. Key of map
   * denotes time, values are dependent variable vectors \return List of maps of
   * dependent variable history
   */
  std::vector<std::map<TimeType, Eigen::VectorXd> >
  getDependentVariableHistory() {
    return dependentVariableHistory_;
  }

  std::vector<std::map<TimeType, double> >
  getCumulativeComputationTimeHistory() {
    return cumulativeComputationTimeHistory_;
  }

  //! Function to return the numerical solution to the equations of motion (base
  //! class interface).
  /*!
   *  Function to return the numerical solution to the equations of motion for
   * last numerical integration. Each vector entry denotes one arc. Key of map
   * denotes time, values are full propagated state vectors. \return List of
   * maps of history of numerically integrated states.
   */
  std::vector<
      std::map<TimeType, Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1> > >
  getEquationsOfMotionNumericalSolutionBase() {
    return getEquationsOfMotionNumericalSolution();
  }

  //! Function to return the numerical solution of the dependent variables (base
  //! class interface)
  /*!
   *  Function to return the numerical solution of the dependent variables for
   * last numerical integration. Each vector entry denotes one arc. Key of map
   * denotes time, values are dependent variable vectors \return List of maps of
   * dependent variable history
   */
  std::vector<std::map<TimeType, Eigen::VectorXd> >
  getDependentVariableNumericalSolutionBase() {
    return getDependentVariableHistory();
  }

  std::vector<std::map<TimeType, double> >
  getCumulativeComputationTimeHistoryBase() {
    return getCumulativeComputationTimeHistory();
  }

  //! Function to reset the environment using an externally provided list of
  //! (numerically integrated) states
  /*!
   *  Function to reset the environment using an externally provided list of
   * (numerically integrated) states, for instance provided by a variational
   * equations solver. \param equationsOfMotionNumericalSolution Vector of state
   * histories (externally provided equationsOfMotionNumericalSolution_) \param
   * dependentVariableHistory Vector of dependent variable histories (externally
   * provided dependentVariableHistory_) \param processSolution True if the new
   * solution is to be immediately processed (default true).
   */
  void manuallySetAndProcessRawNumericalEquationsOfMotionSolution(
      std::vector<std::map<TimeType,
                           Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1> > >
          &equationsOfMotionNumericalSolution,
      std::vector<std::map<TimeType, Eigen::VectorXd> >
          &dependentVariableHistory,
      const bool processSolution = true) {
    // Set equationsOfMotionNumericalSolution_
    equationsOfMotionNumericalSolution_.resize(
        equationsOfMotionNumericalSolution.size());

    for (unsigned int i = 0; i < equationsOfMotionNumericalSolution.size();
         i++) {
      equationsOfMotionNumericalSolution_[i].clear();
      equationsOfMotionNumericalSolution_[i] =
          std::move(equationsOfMotionNumericalSolution[i]);
      arcStartTimes_[i] = equationsOfMotionNumericalSolution_[i].begin()->first;
    }

    // Reset environment with new states.
    if (processSolution) {
      processNumericalEquationsOfMotionSolution();
    }

    dependentVariableHistory_.resize(dependentVariableHistory.size());

    for (unsigned int i = 0; i < dependentVariableHistory.size(); i++) {
      dependentVariableHistory_[i].clear();
      dependentVariableHistory_[i] = std::move(dependentVariableHistory[i]);
    }
  }

  //! Function to get the list of DynamicsStateDerivativeModel objects used for
  //! each arc
  /*!
   * Function to get the list of DynamicsStateDerivativeModel objects used for
   * each arc \return List of DynamicsStateDerivativeModel objects used for each
   * arc
   */
  std::vector<std::shared_ptr<
      DynamicsStateDerivativeModel<TimeType, StateScalarType> > >
  getDynamicsStateDerivative() {
    std::vector<std::shared_ptr<
        DynamicsStateDerivativeModel<TimeType, StateScalarType> > >
        dynamicsStateDerivatives;
    for (unsigned int i = 0; i < singleArcDynamicsSimulators_.size(); i++) {
      dynamicsStateDerivatives.push_back(
          singleArcDynamicsSimulators_.at(i)->getDynamicsStateDerivative());
    }
    return dynamicsStateDerivatives;
  }

  //! Function to get the list of DynamicsSimulator objects used for each arc
  /*!
   * Function to get the list of DynamicsSimulator objects used for each arc
   * \return List of DynamicsSimulator objects used for each arc
   */
  std::vector<
      std::shared_ptr<SingleArcDynamicsSimulator<StateScalarType, TimeType> > >
  getSingleArcDynamicsSimulators() {
    return singleArcDynamicsSimulators_;
  }

  //! Function to retrieve the current state and end times of the arcs
  /*!
   * Function to retrieve the current state and end times of the arcs
   * \return The current state and end times of the arcs
   */
  std::vector<double> getArcStartTimes() { return arcStartTimes_; }

  //! Get whether the integration was completed successfully.
  /*!
   * @copybrief integrationCompletedSuccessfully
   * \return Whether the integration was completed successfully by reaching the
   * termination condition.
   */
  virtual bool integrationCompletedSuccessfully() const {
    for (const std::shared_ptr<
             SingleArcDynamicsSimulator<StateScalarType, TimeType> >
             singleArcDynamicsSimulator : singleArcDynamicsSimulators_) {
      if (!singleArcDynamicsSimulator->integrationCompletedSuccessfully()) {
        return false;
      }
    }
    return true;
  }

  //! This function updates the environment with the numerical solution of the
  //! propagation.
  /*!
   *  This function updates the environment with the numerical solution of the
   * propagation. It sets the propagated dynamics solution as the new input for
   * e.g., the ephemeris object of the boies that were propagated (for
   * translational states).
   */
  void processNumericalEquationsOfMotionSolution() {
    resetIntegratedMultiArcStatesWithEqualArcDynamics(
        equationsOfMotionNumericalSolution_,
        singleArcDynamicsSimulators_.at(0)->getIntegratedStateProcessors(),
        arcStartTimes_);

    if (clearNumericalSolutions_) {
      for (unsigned int i = 0; i < equationsOfMotionNumericalSolution_.size();
           i++) {
        equationsOfMotionNumericalSolution_.at(i).clear();
      }
      equationsOfMotionNumericalSolution_.clear();
    }
  }

 protected:
  //! List of maps of state history of numerically integrated states.
  /*!
   *  List of maps of state history of numerically integrated states. Each entry
   * in the list contains data on a single arc. Key of map denotes time, values
   * are concatenated vectors of body states in order of bodiesToIntegrate
   */
  std::vector<
      std::map<TimeType, Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1> > >
      equationsOfMotionNumericalSolution_;

  //! List of maps of dependent variable history that was saved during numerical
  //! propagation.
  std::vector<std::map<TimeType, Eigen::VectorXd> > dependentVariableHistory_;

  std::vector<std::map<TimeType, double> > cumulativeComputationTimeHistory_;

  //! Objects used to compute the dynamics of the sepatrate arcs
  std::vector<
      std::shared_ptr<SingleArcDynamicsSimulator<StateScalarType, TimeType> > >
      singleArcDynamicsSimulators_;

  //! List of start times of each arc. NOTE: This list is updated after every
  //! propagation.
  std::vector<double> arcStartTimes_;

  //! Event that triggered the termination of the propagation
  std::vector<std::shared_ptr<PropagationTerminationDetails> >
      propagationTerminationReasons_;

  //! Propagator settings used by this objec
  std::shared_ptr<MultiArcPropagatorSettings<StateScalarType> >
      multiArcPropagatorSettings_;
};

#endif  // TUDATBUNDLE_MULTIARCSIMULATOR_H
