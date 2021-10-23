//
// Created by ggarr on 18/10/2021.
//

#ifndef TUDAT_HELPERS_H
#define TUDAT_HELPERS_H

namespace tudat {

    namespace numerical_simulation {

        //! @get_docstring(getInitialStatesOfBodiesFromFrameManager)
        template<typename TimeType = double, typename StateScalarType = double>
        Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1> getInitialStatesOfBodiesFromFrameManager(
                const std::vector <std::string> &bodiesToIntegrate,
                const std::vector <std::string> &centralBodies,
                const simulation_setup::SystemOfBodies &bodies,
                const TimeType initialTime,
                const std::shared_ptr <ephemerides::ReferenceFrameManager> frameManager) {
            // Set initial states of bodies to integrate.
            Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1> systemInitialState =
                    Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1>::Zero(bodiesToIntegrate.size() * 6, 1);
            std::shared_ptr <ephemerides::Ephemeris> ephemerisOfCurrentBody;

            // Iterate over all bodies.
            for (unsigned int i = 0; i < bodiesToIntegrate.size(); i++) {
                ephemerisOfCurrentBody = bodies.at(bodiesToIntegrate.at(i))->getEphemeris();

                if (!ephemerisOfCurrentBody) {
                    throw std::runtime_error("Could not determine initial state for body " + bodiesToIntegrate.at(i) +
                                             " because it does not have a valid Ephemeris object.");
                }

                // Get body initial state from ephemeris
                systemInitialState.segment(i * 6, 6) = ephemerisOfCurrentBody->getTemplatedStateFromEphemeris<
                        StateScalarType, TimeType>(initialTime);

                // Correct initial state if integration origin and ephemeris origin are not equal.
                if (centralBodies.at(i) != ephemerisOfCurrentBody->getReferenceFrameOrigin()) {
                    std::shared_ptr <ephemerides::Ephemeris> correctionEphemeris =
                            frameManager->getEphemeris(ephemerisOfCurrentBody->getReferenceFrameOrigin(),
                                                       centralBodies.at(i));
                    systemInitialState.segment(i * 6, 6) -= correctionEphemeris->getTemplatedStateFromEphemeris<
                            StateScalarType, TimeType>(initialTime);
                }
            }
            return systemInitialState;
        }

        //! @get_docstring(getInitialRotationalStatesOfBodies)
        template<typename TimeType = double, typename StateScalarType = double>
        Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1> getInitialRotationalStatesOfBodies(
                const std::vector <std::string> &bodiesToIntegrate,
                const std::vector <std::string> &baseOrientations,
                const simulation_setup::SystemOfBodies &bodies,
                const TimeType initialTime) {
            // Set initial states of bodies to integrate.
            Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1> systemInitialState =
                    Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1>::Zero(bodiesToIntegrate.size() * 7, 1);
            std::shared_ptr <ephemerides::RotationalEphemeris> rotationModelOfCurrentBody;

            // Iterate over all bodies.
            for (unsigned int i = 0; i < bodiesToIntegrate.size(); i++) {
                rotationModelOfCurrentBody = bodies.at(bodiesToIntegrate.at(i))->getRotationalEphemeris();

                if (!rotationModelOfCurrentBody) {
                    throw std::runtime_error("Could not determine initial state for body " + bodiesToIntegrate.at(i) +
                                             " because it does not have a valid RotationalEphemeris object.");
                }

                // Get body initial state from ephemeris
                systemInitialState.segment(i * 7, 7) = rotationModelOfCurrentBody->getRotationStateVector(
                        initialTime).template cast<StateScalarType>();

                // Correct initial state if integration origin and rotation model origin are not equal.
                if (baseOrientations.at(i) != rotationModelOfCurrentBody->getBaseFrameOrientation()) {
                    throw std::runtime_error("Error, cannot get initial rotational state w.r.t. non-base frame");
                }
            }
            return systemInitialState;
        }

        //! @get_docstring(getInitialStatesOfBodies)
        template<typename TimeType = double, typename StateScalarType = double>
        Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1> getInitialStatesOfBodies(
                const std::vector <std::string> &bodiesToIntegrate,
                const std::vector <std::string> &centralBodies,
                const simulation_setup::SystemOfBodies &bodies,
                const TimeType initialTime) {
            // Create ReferenceFrameManager and call overloaded function.
            return getInitialStatesOfBodiesFromFrameManager<TimeType, StateScalarType>(
                    bodiesToIntegrate, centralBodies, bodies, initialTime,
                    simulation_setup::createFrameManager(bodies.getMap()));
        }

        //! @get_docstring(getInitialStateOfBody)
        template<typename TimeType = double, typename StateScalarType = double>
        Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1> getInitialStateOfBody(
                const std::string &bodyToIntegrate,
                const std::string &centralBody,
                const simulation_setup::SystemOfBodies &bodies,
                const TimeType initialTime) {
            return getInitialStatesOfBodies<TimeType, StateScalarType>(
                    {bodyToIntegrate}, {centralBody}, bodies, initialTime);
        }

        //! @get_docstring(getInitialRotationalStateOfBody)
        template<typename TimeType = double, typename StateScalarType = double>
        Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1> getInitialRotationalStateOfBody(
                const std::string &bodyToIntegrate,
                const std::string &baseOrientation,
                const simulation_setup::SystemOfBodies &bodies,
                const TimeType initialTime) {
            return getInitialRotationalStatesOfBodies<TimeType, StateScalarType>(
                    std::vector < std::string > {bodyToIntegrate}, std::vector < std::string > {baseOrientation},
                    bodies, initialTime);
        }

        //! @get_docstring(getInitialArcWiseStateOfBody)
        template<typename TimeType = double, typename StateScalarType = double>
        Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1> getInitialArcWiseStateOfBody(
                const std::string &bodyToIntegrate,
                const std::vector <std::string> &centralBodies,
                const simulation_setup::SystemOfBodies &bodies,
                const std::vector <TimeType> arcStartTimes) {
            Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1> initialStates = Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1>::Zero(
                    6 * arcStartTimes.size(), 1);
            for (unsigned int i = 0; i < arcStartTimes.size(); i++) {
                initialStates.block(6 * i, 0, 6, 1) = getInitialStateOfBody<double, StateScalarType>(
                        bodyToIntegrate, centralBodies.at(i), bodies, arcStartTimes.at(i));
            }
            return initialStates;
        }

////! Function to get a vector of initial states from a vector of propagator settings
///*!
// *  Function to get a vector of initial states from a vector of propagator settings.
// *  \param propagatorSettings List of propagator settings
// *  \return List of initial states, as retrieved from propagatorSettings list.
// */
//        template<typename StateScalarType = double>
//        std::vector<Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1> > getInitialStatesPerArc(
//                const std::vector<std::shared_ptr<PropagatorSettings<StateScalarType> > > propagatorSettings) {
//        std::vector<Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1> > initialStatesList;
//        for (unsigned int i = 0; i < propagatorSettings.size(); i++) {
//        initialStatesList.push_back(propagatorSettings.at(i)->getInitialStates());
//    }
//
//    return initialStatesList;
//}
//
////! Function to get the initial state of a translational state arc from the previous state's numerical solution
///*!
// *  Function to get the initial state of a translational state arc from the previous state's numerical solution
// *  \param previousArcDynamicsSolution Numerical solution of previous arc
// *  \param currentArcInitialTime Start time of current arc
// *  \return Interpolated initial state of current arc
// */
//template<typename StateScalarType = double, typename TimeType = double>
//Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1> getArcInitialStateFromPreviousArcResult(
//        const std::map<TimeType, Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1> > &previousArcDynamicsSolution,
//        const double currentArcInitialTime) {
//    Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1> currentArcInitialState;
//    {
//        // Check if overlap exists
//        if (previousArcDynamicsSolution.rbegin()->first < currentArcInitialTime) {
//            throw std::runtime_error(
//                    "Error when getting initial arc state from previous arc: no arc overlap");
//        } else {
//            int currentIndex = 0;
//            int initialTimeIndex = -1;
//
//            std::map<TimeType, Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1> > initialStateInterpolationMap;
//
//            // Set sub-part of previous arc to interpolate for current arc
//            for (typename std::map<TimeType, Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1> >::
//            const_reverse_iterator previousArcIterator = previousArcDynamicsSolution.rbegin();
//                    previousArcIterator != previousArcDynamicsSolution.rend(); previousArcIterator++) {
//                initialStateInterpolationMap[previousArcIterator->first] = previousArcIterator->second;
//                if (initialTimeIndex < 0) {
//                    if (previousArcIterator->first < currentArcInitialTime) {
//                        initialTimeIndex = currentIndex;
//                    }
//                } else {
//                    if (currentIndex - initialTimeIndex > 5) {
//                        break;
//                    }
//                }
//                currentIndex++;
//            }
//
//            // Interpolate to obtain initial state of current arc
//            currentArcInitialState =
//                    std::make_shared<interpolators::LagrangeInterpolator<
//                    TimeType, Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1>, long double> >(
//                    initialStateInterpolationMap, 8)->interpolate(currentArcInitialTime);
//
//        }
//    }
//    return currentArcInitialState;
//}

    } // namespace numerical_simulation

} // namespace tudat

#endif //TUDAT_HELPERS_H
