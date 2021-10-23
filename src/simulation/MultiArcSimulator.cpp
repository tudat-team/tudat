//
// Created by ggarr on 18/10/2021.
//

#include "MultiArcSimulator.h"

#if(TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS)
extern template class MultiArcDynamicsSimulator< long double, double >;
extern template class MultiArcDynamicsSimulator< double, Time >;
extern template class MultiArcDynamicsSimulator< long double, Time >;
#endif


} // namespace propagators
