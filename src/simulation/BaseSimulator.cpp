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

namespace tudat {

    namespace numerical_simulation {

        extern template
        class BaseSimulator<double, double>;

#if(TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS)
        extern template class BaseSimulator< long double, double >;
        extern template class BaseSimulator< double, Time >;
        extern template class BaseSimulator< long double, Time >;
#endif

    } // namespace numerical_simulation

} // namespace tudat

