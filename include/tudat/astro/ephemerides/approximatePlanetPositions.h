/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Standish, E.M. Keplerian Elements for Approximate Positions of the Major Planets,
 *          http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf, last accessed: 24 February, 2011.
 *
 */

#ifndef TUDAT_APPROXIMATE_PLANET_POSITIONS_H
#define TUDAT_APPROXIMATE_PLANET_POSITIONS_H

#include <memory>

#include "tudat/math/basic/mathematicalConstants.h"

#include "tudat/astro/basic_astro/convertMeanToEccentricAnomalies.h"
#include "tudat/astro/ephemerides/approximatePlanetPositionsBase.h"

#include "tudat/basics/basicTypedefs.h"

namespace tudat
{
namespace ephemerides
{

//! Ephemeris class using JPL "Approximate Positions of Major Planets".
/*!
 * Ephemeris class using JPL "Approximate Positions of Major Planets".
 */
class ApproximatePlanetPositions : public ApproximatePlanetPositionsBase
{
public:

    using Ephemeris::getCartesianState;

    //! Default constructor.
    /*!
     * Default constructor that initializes the class from the body for which the position is
     * approximated and the gravitational parameter of the Sun (default 1.32712440018e20). Other
     * members are initialized to default values (NAN).
     *
     * \param bodyWithEphemerisData The body for which the position is approximated.
     * \param sunGravitationalParameter The gravitational parameter of the Sun [m^3/s^2].
     * \sa BodiesWithEphemerisData, ApproximatePlanetPositionsBase.
     */
    ApproximatePlanetPositions( BodiesWithEphemerisData bodyWithEphemerisData,
                                const double sunGravitationalParameter = 1.32712440018e20 )
        : ApproximatePlanetPositionsBase( sunGravitationalParameter ),
          eccentricAnomalyAtGivenJulianDate_( TUDAT_NAN ),
          longitudeOfPerihelionAtGivenJulianDate_( TUDAT_NAN ),
          meanAnomalyAtGivenJulianDate_( TUDAT_NAN ),
          trueAnomalyAtGivenJulianData_( TUDAT_NAN )
    {
        setPlanet( bodyWithEphemerisData );
    }

    ApproximatePlanetPositions( const std::string& bodyName,
                                const double sunGravitationalParameter = 1.32712440018e20 )
        : ApproximatePlanetPositionsBase( sunGravitationalParameter ),
          eccentricAnomalyAtGivenJulianDate_( TUDAT_NAN ),
          longitudeOfPerihelionAtGivenJulianDate_( TUDAT_NAN ),
          meanAnomalyAtGivenJulianDate_( TUDAT_NAN ),
          trueAnomalyAtGivenJulianData_( TUDAT_NAN )
    {
        BodiesWithEphemerisData bodyWithEphemerisData = ApproximatePlanetPositionsBase::getBodiesWithEphemerisDataId(
                    bodyName );
        setPlanet( bodyWithEphemerisData );
    }

    //! Get cartesian state from ephemeris.
    /*!
     * Returns cartesian state from ephemeris.
     * \param secondsSinceEpoch Seconds since epoch.
     * \return State in Cartesian elements from ephemeris.
     */
    Eigen::Vector6d getCartesianState(
            const double secondsSinceEpoch );

    //! Get keplerian state from ephemeris.
    /*!
     * Returns keplerian state in from ephemeris.
     * \param secondsSinceEpoch Seconds since epoch.
     * \return State in Keplerian elements from ephemeris.
     */
    Eigen::Vector6d getKeplerianStateFromEphemeris(
            const double secondsSinceEpoch );

protected:

private:

    //! Eccentric anomaly at given Julian date.
    /*!
     * Eccentric anomaly of planet at given Julian date.
     */
    double eccentricAnomalyAtGivenJulianDate_;

    //! Longitude of perihelion at given Julian date.
    /*!
     * Longitude of perihelion of planet at given Julian date.
     */
    double longitudeOfPerihelionAtGivenJulianDate_;

    //! Mean anomaly at given Julian date.
    /*!
     * Mean anomaly of planet at given Julian date.
     */
    double meanAnomalyAtGivenJulianDate_;

    //! True anomaly at given Julian date.
    /*!
     * True anomaly of planet at given Julian date.
     */
    double trueAnomalyAtGivenJulianData_;
};

//! Typedef for shared-pointer to ApproximatePlanetPositions object.
typedef std::shared_ptr< ApproximatePlanetPositions > ApproximatePlanetPositionsPointer;

} // namespace ephemerides
} // namespace tudat

#endif // TUDAT_APPROXIMATE_PLANET_POSITIONS_H
