/*", "Copyright (c) 2010-2019, Delft University of Technology
 *", "All rigths reserved
 *
 *", "This file is part of the Tudat. Redistribution and use in source and
 *", "binary forms, with or without modification, are permitted exclusively
 *", "under the terms of the Modified BSD license. You should have received
 *", "a copy of the license with this file. If not, please or visit:
 *", "http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_MPC_OBSERVATIONS_H
#define TUDAT_MPC_OBSERVATIONS_H

#include <vector>

#include <boost/bind.hpp>
#include <memory>
#include <functional>

#include <Eigen/Core>

#include "tudat/basics/basicTypedefs.h"
#include "tudat/basics/timeType.h"
#include "tudat/basics/tudatTypeTraits.h"
#include "tudat/basics/utilities.h"

#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/simulation/estimation_setup/observationOutput.h"

namespace tudat
{

namespace observation_models
{


static const std::map< std::string, std::string > mpcCatalogueNames = {
{ "a", "USNO-A1.0" },
{ "b", "USNO-SA1.0" },
{ "c", "USNO-A2.0" },
{ "d", "USNO-SA2.0" },
{ "e", "UCAC-1" },
{ "f", "Tycho-1" },
{ "g", "Tycho-2" },
{ "h", "GSC-1.0" },
{ "i", "GSC-1.1" },
{ "j", "GSC-1.2" },
{ "k", "GSC-2.2" },
{ "l", "ACT" },
{ "m", "GSC-ACT" },
{ "n", "SDSS-DR8" },
{ "o", "USNO-B1.0" },
{ "p", "PPM" },
{ "q", "UCAC-4" },
{ "r", "UCAC-2" },
{ "s", "USNO-B2.0" },
{ "t", "PPMXL" },
{ "u", "UCAC-3" },
{ "v", "NOMAD" },
{ "w", "CMC-14" },
{ "x", "Hipparcos 2" },
{ "y", "Hipparcos" },
{ "z", "GSC (version unspecified)" },
{ "A", "AC" },
{ "B", "SAO 1984" },
{ "C", "SAO" },
{ "D", "AGK 3" },
{ "E", "FK4" },
{ "F", "ACRS" },
{ "G", "Lick Gaspra Catalogue" },
{ "H", "Ida93 Catalogue" },
{ "I", "Perth 70" },
{ "J", "COSMOS/UKST Southern Sky Catalogue" },
{ "K", "Yale" },
{ "L", "2MASS" },
{ "M", "GSC-2.3" },
{ "N", "SDSS-DR7" },
{ "O", "SST-RC1" },
{ "P", "MPOSC3" },
{ "Q", "CMC-15" },
{ "R", "SST-RC4" },
{ "S", "URAT-1" },
{ "T", "URAT-2" },
{ "U", "Gaia-DR1" },
{ "V", "Gaia-DR2" },
{ "W", "Gaia-DR3" },
{ "X", "Gaia-EDR3" },
{ "Y", "UCAC-5" },
{ "Z", "ATLAS-2" },
{ "0", "IHW" },
{ "1", "PS1_DR1" },
{ "2", "PS1_DR2" },
{ "3", "Gaia_Int" },
{ "4", "GZ" },
{ "5", "USNO-UBAD" },
{ "6", "Gaia2016" }
};

} // namespace observation_models

} // namespace tudat

#endif // TUDAT_MPC_OBSERVATIONS_H
