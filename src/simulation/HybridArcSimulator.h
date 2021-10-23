//
// Created by ggarr on 18/10/2021.
//

#ifndef TUDATBUNDLE_HYBRIDARCSIMULATOR_H
#define TUDATBUNDLE_HYBRIDARCSIMULATOR_H

//! @get_docstring(HybridArcDynamicsSimulator)
template<typename StateScalarType = double, typename TimeType = double>
class HybridArcDynamicsSimulator : public DynamicsSimulator<StateScalarType, TimeType> {
public:

    //! Using statemebts
    using DynamicsSimulator<StateScalarType, TimeType>::bodies_;
    using DynamicsSimulator<StateScalarType, TimeType>::clearNumericalSolutions_;

    //! Constructor of multi-arc simulator for same integration settings per arc.
    /*!
     *  Constructor of multi-arc simulator for same integration settings per arc.
     *  \param bodies Map of bodies (with names) of all bodies in integration.
     *  \param integratorSettings Integrator settings for numerical integrator, used for all arcs.
     *  \param propagatorSettings Propagator settings for dynamics (must be of multi arc type)
     *  \param arcStartTimes Times at which the separate arcs start, for the multi-arc case
     *  \param areEquationsOfMotionToBeIntegrated Boolean to denote whether equations of motion should be integrated at
     *  the end of the contructor or not.
     *  \param clearNumericalSolutions Boolean to determine whether to clear the raw numerical solution member variables
     *  after propagation and resetting ephemerides (default true).
     *  \param setIntegratedResult Boolean to determine whether to automatically use the integrated results to set
     *  ephemerides (default true).
     *  \param addSingleArcBodiesToMultiArcDynamics Boolean denoting whether to add single arc bodies to multi-arc
     *  dynamics (default false).
     */
    HybridArcDynamicsSimulator(
            const simulation_setup::SystemOfBodies &bodies,
            const std::shared_ptr<numerical_integrators::IntegratorSettings<TimeType> > integratorSettings,
            const std::shared_ptr<PropagatorSettings<StateScalarType> > propagatorSettings,
            const std::vector<double> arcStartTimes,
            const bool areEquationsOfMotionToBeIntegrated = true,
            const bool clearNumericalSolutions = true,
            const bool setIntegratedResult = true,
            const bool addSingleArcBodiesToMultiArcDynamics = false) :
            DynamicsSimulator<StateScalarType, TimeType>(
                    bodies, clearNumericalSolutions, setIntegratedResult) {
        std::shared_ptr<HybridArcPropagatorSettings<StateScalarType> > hybridPropagatorSettings =
                std::dynamic_pointer_cast<HybridArcPropagatorSettings<StateScalarType> >(propagatorSettings);
        if (hybridPropagatorSettings == nullptr) {
            throw std::runtime_error(
                    "Error when making HybridArcDynamicsSimulator, propagator settings are incompatible");
        }

        singleArcDynamicsSize_ = hybridPropagatorSettings->getSingleArcPropagatorSettings()->getPropagatedStateSize();
        multiArcDynamicsSize_ = hybridPropagatorSettings->getMultiArcPropagatorSettings()->getPropagatedStateSize();

        if (!addSingleArcBodiesToMultiArcDynamics) {
            if (!setIntegratedResult) {
                std::cerr
                        << "Warning in hybrid dynamics simulator, setIntegratedResult is false, but single arc propagation "
                           "result will be set in environment for consistency with multi-arc " << std::endl;
            }
            singleArcDynamicsSimulator_ = std::make_shared<SingleArcDynamicsSimulator<StateScalarType, TimeType> >(
                    bodies, integratorSettings, hybridPropagatorSettings->getSingleArcPropagatorSettings(),
                            false, false, true);
            multiArcDynamicsSimulator_ = std::make_shared<MultiArcDynamicsSimulator<StateScalarType, TimeType> >(
                    bodies, integratorSettings, hybridPropagatorSettings->getMultiArcPropagatorSettings(),
                            arcStartTimes,
                            false, false, setIntegratedResult);
        } else {
            throw std::runtime_error("Cannot yet add single-arc bodies to multi-arc propagation");
        }

        if (areEquationsOfMotionToBeIntegrated) {
            integrateEquationsOfMotion(hybridPropagatorSettings->getInitialStates());
        }
    }

    HybridArcDynamicsSimulator(
            const simulation_setup::SystemOfBodies &bodies,
            const std::shared_ptr<numerical_integrators::IntegratorSettings<TimeType> > singleArcIntegratorSettings,
            const std::shared_ptr<numerical_integrators::IntegratorSettings<TimeType> > multiArcIntegratorSettings,
            const std::shared_ptr<PropagatorSettings<StateScalarType> > propagatorSettings,
            const std::vector<double> arcStartTimes,
            const bool areEquationsOfMotionToBeIntegrated = true,
            const bool clearNumericalSolutions = true,
            const bool setIntegratedResult = true,
            const bool addSingleArcBodiesToMultiArcDynamics = false) :
            DynamicsSimulator<StateScalarType, TimeType>(
                    bodies, clearNumericalSolutions, setIntegratedResult) {
        std::shared_ptr<HybridArcPropagatorSettings<StateScalarType> > hybridPropagatorSettings =
                std::dynamic_pointer_cast<HybridArcPropagatorSettings<StateScalarType> >(propagatorSettings);
        if (hybridPropagatorSettings == nullptr) {
            throw std::runtime_error(
                    "Error when making HybridArcDynamicsSimulator, propagator settings are incompatible");
        }

        singleArcDynamicsSize_ = hybridPropagatorSettings->getSingleArcPropagatorSettings()->getPropagatedStateSize();
        multiArcDynamicsSize_ = hybridPropagatorSettings->getMultiArcPropagatorSettings()->getPropagatedStateSize();

        if (!addSingleArcBodiesToMultiArcDynamics) {
            if (!setIntegratedResult) {
                std::cerr
                        << "Warning in hybrid dynamics simulator, setIntegratedResult is false, but single arc propagation "
                           "result will be set in environment for consistency with multi-arc " << std::endl;
            }
            singleArcDynamicsSimulator_ = std::make_shared<SingleArcDynamicsSimulator<StateScalarType, TimeType> >(
                    bodies, singleArcIntegratorSettings,
                            hybridPropagatorSettings->getSingleArcPropagatorSettings(),
                            false, false, true);
            multiArcDynamicsSimulator_ = std::make_shared<MultiArcDynamicsSimulator<StateScalarType, TimeType> >(
                    bodies, multiArcIntegratorSettings,
                            hybridPropagatorSettings->getMultiArcPropagatorSettings(), arcStartTimes,
                            false, false, setIntegratedResult);
        } else {
            throw std::runtime_error("Cannot yet add single-arc bodies to multi-arc propagation");
        }

        if (areEquationsOfMotionToBeIntegrated) {
            integrateEquationsOfMotion(hybridPropagatorSettings->getInitialStates());
        }
    }

    //! Destructor
    ~HybridArcDynamicsSimulator() {}

    //! This function numerically (re-)integrates the equations of motion, using concatenated states for single and multi-arcs
    /*!
     *  This function numerically (re-)integrates the equations of motion, using the settings set through the constructor
     *  and a new initial state vector provided here. The raw results are set in the equationsOfMotionNumericalSolution_
     *  \param initialGlobalStates Initial state vector that is to be used for numerical integration. Note that this state
     *  should be in the correct frame (i.e. corresponding to centralBodies in propagatorSettings_), but not in the propagator-
     *  specific form (i.e Encke, Gauss, etc. for translational dynamics). The states for all arcs must be concatenated in
     *  order into a single Eigen Vector, starting with the single-arc states, followed by the mulit-arc states
     */
    void integrateEquationsOfMotion(
            const Eigen::Matrix<StateScalarType, Eigen::Dynamic, Eigen::Dynamic> &initialGlobalStates) {
        singleArcDynamicsSimulator_->integrateEquationsOfMotion(
                initialGlobalStates.block(0, 0, singleArcDynamicsSize_, 1));
        multiArcDynamicsSimulator_->integrateEquationsOfMotion(
                initialGlobalStates.block(singleArcDynamicsSize_, 0, multiArcDynamicsSize_, 1));

    }

    //! This function updates the environment with the numerical solution of the propagation
    /*!
     *  This function updates the environment with the numerical solution of the propagation
     *  (no additional functionality in hybrid arc). Function may be used to process manually updtaed propagation results in
     *  single and/or multi-arc model
     */
    void processNumericalEquationsOfMotionSolution() {
        singleArcDynamicsSimulator_->processNumericalEquationsOfMotionSolution();
        multiArcDynamicsSimulator_->processNumericalEquationsOfMotionSolution();
    }

    //! Function to retrieve the single-arc dynamics simulator
    /*!
     * Function to retrieve the single-arc dynamics simulator
     * \return Single-arc dynamics simulator
     */
    std::shared_ptr<SingleArcDynamicsSimulator<StateScalarType, TimeType> > getSingleArcDynamicsSimulator() {
        return singleArcDynamicsSimulator_;
    }

    //! Function to retrieve the multi-arc dynamics simulator
    /*!
     * Function to retrieve the multi-arc dynamics simulator
     * \return Multi-arc dynamics simulator
     */
    std::shared_ptr<MultiArcDynamicsSimulator<StateScalarType, TimeType> > getMultiArcDynamicsSimulator() {
        return multiArcDynamicsSimulator_;
    }

    //! Get whether the integration was completed successfully.
    /*!
     * Get whether the integration was completed successfully.
     * \return Whether the integration was completed successfully by reaching the termination condition.
     */
    virtual bool integrationCompletedSuccessfully() const {
        return !(singleArcDynamicsSimulator_->integrationCompletedSuccessfully() == false ||
                 multiArcDynamicsSimulator_->integrationCompletedSuccessfully() == false);

    }

    //! Function to return the numerical solution to the equations of motion (base class interface).
    /*!
     *  Function to return the numerical solution to the equations of motion for last numerical integration. First vector entry
     *  contains single-arc results. Each subsequent vector entry contains one of the multi-arcs. Key of map denotes time,
     *  values are full propagated state vectors.
     *  \return List of maps of history of numerically integrated states.
     */
    std::vector<std::map<TimeType, Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1> > >
    getEquationsOfMotionNumericalSolutionBase() {
        std::vector<std::map<TimeType, Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1> > >
                                                                                      numericalSolution = singleArcDynamicsSimulator_->getEquationsOfMotionNumericalSolutionBase();
        std::vector<std::map<TimeType, Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1> > >
                                                                                      multiArcNumericalSolution = multiArcDynamicsSimulator_->getEquationsOfMotionNumericalSolutionBase();
        numericalSolution.insert(numericalSolution.end(), multiArcNumericalSolution.begin(),
                                 multiArcNumericalSolution.end());

        return numericalSolution;
    }

    //! Function to return the numerical solution to the dependent variables (base class interface).
    /*!
     *  Function to return the numerical solution to the dependent variables for last numerical integration. First vector entry
     *  contains single-arc results. Each subsequent vector entry contains one of the multi-arcs. Key of map denotes time,
     *  values are dependent variables vectors
     *  \return List of maps of dependent variable history
     */
    std::vector<std::map<TimeType, Eigen::VectorXd> > getDependentVariableNumericalSolutionBase() {
        std::vector<std::map<TimeType, Eigen::VectorXd> >
                numericalSolution = singleArcDynamicsSimulator_->getDependentVariableNumericalSolutionBase();
        std::vector<std::map<TimeType, Eigen::VectorXd> >
                multiArcNumericalSolution = multiArcDynamicsSimulator_->getDependentVariableNumericalSolutionBase();
        numericalSolution.insert(numericalSolution.end(), multiArcNumericalSolution.begin(),
                                 multiArcNumericalSolution.end());

        return numericalSolution;
    }

    //! Function to return the map of cumulative computation time history that was saved during numerical propagation.
    /*!
     *  Function to return the map of cumulative computation time history that was saved during numerical propagation.  First vector
     *  entry contains single-arc results. Each subsequent vector entry contains one of the multi-arcs. Key of map denotes time,
     *  values are computation times.
     *  \return Vector is size 1, with entry: map of cumulative computation time history that was saved during numerical propagation.
     */
    std::vector<std::map<TimeType, double> > getCumulativeComputationTimeHistoryBase() {
        std::vector<std::map<TimeType, double> >
                computationTime = singleArcDynamicsSimulator_->getCumulativeComputationTimeHistoryBase();
        std::vector<std::map<TimeType, double> >
                multiArcComputationTime = multiArcDynamicsSimulator_->getCumulativeComputationTimeHistoryBase();
        computationTime.insert(computationTime.end(), multiArcComputationTime.begin(),
                               multiArcComputationTime.end());

        return computationTime;
    }

protected:

    //! Object used to propagate single-arc dynamics
    std::shared_ptr<SingleArcDynamicsSimulator<StateScalarType, TimeType> > singleArcDynamicsSimulator_;

    //! Object used to propagate multi-arc dynamics
    std::shared_ptr<MultiArcDynamicsSimulator<StateScalarType, TimeType> > multiArcDynamicsSimulator_;

    //! Size of single-arc (initial) state vector
    int singleArcDynamicsSize_;

    //! Size of multi-arc concatenated initial state vector
    int multiArcDynamicsSize_;
};


extern template
class SingleArcDynamicsSimulator<double, double>;

extern template
class MultiArcDynamicsSimulator<double, double>;

extern template
class HybridArcDynamicsSimulator<double, double>;


#endif //TUDATBUNDLE_HYBRIDARCSIMULATOR_H
