//
// Created by ggarr on 23/10/2021.
//

#ifndef TUDAT_COMPUTATIONEXPENSE_H
#define TUDAT_COMPUTATIONEXPENSE_H

#include <tudat/simulation/BaseObserver.h>

namespace tudat {

namespace numerical_simulation {

//! @get_docstring(ComputationExpense)
template <typename StateScalarType, typename TimeType, typename TimeStepType,
          typename ObservedType>
class ComputationExpense : public BaseObserver<StateScalarType, TimeType,
                                               TimeStepType, ObservedType> {
 public:
  //! @get_docstring(ComputationExpense.constructor)
  ComputationExpense();

  //! @get_docstring(ComputationExpense.startTimedInterval)
  void startTimedInterval();

  //! @get_docstring(ComputationExpense.endTimedInterval)
  void endTimedInterval();

  //! @get_docstring(ComputationExpense.getCumulativeHistory)
  std::map<TimeType, double> getCumulativeHistory();

  //! @get_docstring(ComputationExpense.getIntervalHistory)
  std::map<TimeType, double> getTimedIntervalHistory();

  //! @get_docstring(ComputationExpense.getIntervalHistory)
  //  void update() override;

 private:
  //! @get_docstring(ComputationExpense.timedIntervalHistory_)
  std::map<TimeType, double> timedIntervalHistory_;
};

}  // namespace numerical_simulation

}  // namespace tudat

#endif  // TUDAT_COMPUTATIONEXPENSE_H
