//
// Created by ggarr on 18/10/2021.
//

#include "HybridArcSimulator.h"


#if(TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS)
extern template class HybridArcDynamicsSimulator< long double, double >;
extern template class HybridArcDynamicsSimulator< double, Time >;
extern template class HybridArcDynamicsSimulator< long double, Time >;
#endif


} // namespace propagators
