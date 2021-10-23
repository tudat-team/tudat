//
// Created by ggarr on 18/10/2021.
//

#include <tudat/simulation/SingleArcSimulator.h>

namespace tudat {

    namespace numerical_simulation {

        extern template
        class SingleArcSimulator<double, double, double>;

#if(TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS)
        extern template class SingleArcSimulator< long double, double , long double>;
        extern template class SingleArcSimulator< double, Time , long double>;
        extern template class SingleArcSimulator< long double, Time , long double>;
#endif

    } // namespace tudat

} // namespace numerical_simulation
