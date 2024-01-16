/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <algorithm>

#include "tudat/simulation/estimation_setup/processTrackingTxtFile.h"

namespace tudat
{
namespace observation_models
{

std::vector<ObservableType> findAvailableObservableTypes(std::vector<input_output::TrackingDataType> availableDataTypes)
{
  std::vector<ObservableType> availableObservableTypes;
  for (const auto& pair : dataTypeToObservableMap) {
    std::vector<input_output::TrackingDataType> requiredDataTypeSet = pair.first;
    if (utilities::containsAll(availableDataTypes, requiredDataTypeSet))
      availableObservableTypes.push_back(pair.second);
  }
  return availableObservableTypes;
}

void ProcessedTrackingTxtFileContents::updateObservations()
{
  const auto& observableTypes = getObservableTypes();
  linkEndsSet_ = utilities::vectorToSet(linkEndsVector_);
}

void ProcessedTrackingTxtFileContents::updateObservationTimes()
{

  observationTimes_.clear();
  const auto& dataMap = rawTrackingTxtFileContents_->getDoubleDataMap();
  const auto& numDataRows = rawTrackingTxtFileContents_->getNumRows();
  TimeRepresentation timeRepresentation = getTimeRepresentation();

  switch (timeRepresentation) {
    case tdb_seconds_j2000: {
      observationTimes_ = dataMap.at(input_output::TrackingDataType::tdb_time_j2000);
      break;
    }
    case calendar_day_time: {
      std::vector<double> observationJulianDaysSinceJ2000 = utilities::convertVectors(
          basic_astrodynamics::convertCalendarDateToJulianDaySinceJ2000<double>,
          dataMap.at(input_output::TrackingDataType::year),
          dataMap.at(input_output::TrackingDataType::month),
          dataMap.at(input_output::TrackingDataType::day),
          dataMap.at(input_output::TrackingDataType::hour),
          dataMap.at(input_output::TrackingDataType::minute),
          dataMap.at(input_output::TrackingDataType::second)
      );

      std::vector<double> observationTimesUtc;
      for (double julianDaySinceJ2000 : observationJulianDaysSinceJ2000) {
        observationTimesUtc.push_back(julianDaySinceJ2000 * physical_constants::JULIAN_DAY);
      }
      observationTimes_ = computeObservationTimesTdbFromJ2000(observationTimesUtc);

      break;
    }
    default: {
      throw std::runtime_error("Error while processing tracking txt file: Time representation not recognised or implemented.");
    }
  }
}

std::vector<double> ProcessedTrackingTxtFileContents::computeObservationTimesTdbFromJ2000(std::vector<double> observationTimesUtc)
{
  earth_orientation::TerrestrialTimeScaleConverter timeScaleConverter = earth_orientation::TerrestrialTimeScaleConverter();
  if (linkEndsVector_.size() != observationTimesUtc.size()) {
    throw std::runtime_error("Error while processing tracking data: vector of linkEnds and observationTimes not of equal size");
  }

  std::vector<Eigen::Vector3d> groundStationPositions;
  for (const auto& linkEnds : linkEndsVector_) {
    std::string currentGroundStation = linkEnds.at(receiver).getStationName(); // TODO: what if transmitter and receiver different?
    groundStationPositions.push_back(simulation_setup::getApproximateGroundStationPositionFromFile(currentGroundStation));
  }
  std::vector<double> observationTimesTdb = timeScaleConverter.getCurrentTimes(basic_astrodynamics::utc_scale,
                                                                               basic_astrodynamics::tdb_scale,
                                                                               observationTimesUtc,
                                                                               groundStationPositions);
  return observationTimesTdb;
}

void ProcessedTrackingTxtFileContents::updateLinkEnds()
{
  linkEndsVector_.clear();
  const auto& dataMap = rawTrackingTxtFileContents_->getDoubleDataMap();
  const auto& metaDataStrMap = rawTrackingTxtFileContents_->getMetaDataStrMap();
  const auto& numDataRows = rawTrackingTxtFileContents_->getNumRows();
  LinkEndsRepresentation linkEndsRepresentation = getLinkEndsRepresentation();

  switch (linkEndsRepresentation) {

    // TODO: make a cleaner implementation to allow adding different ways of providing the link ends easily
    case dsn_transmitting_receiving_station_nr: {
      const auto& dsnTransmitterIds = dataMap.at(input_output::TrackingDataType::dsn_transmitting_station_nr);
      const auto& dsnReceiverIds = dataMap.at(input_output::TrackingDataType::dsn_receiving_station_nr);

      for (size_t i = 0; i < numDataRows; ++i) {
        std::string transmitterName = simulation_setup::getGroundStationCodeFromFile(dsnTransmitterIds[i]);
        std::string receiverName = simulation_setup::getGroundStationCodeFromFile(dsnReceiverIds[i]);
        LinkEnds currentLinkEnds{
            {transmitter, LinkEndId("Earth", transmitterName)},
            {reflector, LinkEndId(spacecraftName_, "Antenna")},
            {receiver, LinkEndId("Earth", receiverName)},
        };
        linkEndsVector_.push_back(currentLinkEnds);
      }
      break;
    }

    case vlbi_station: {

      if (metaDataStrMap.count(input_output::TrackingDataType::vlbi_station_name)) {
        std::string vlbi_station_name = metaDataStrMap.at(input_output::TrackingDataType::vlbi_station_name);
        LinkEnds constantLinkEnds{
            {transmitter, LinkEndId(spacecraftName_, "Antenna")},
            {receiver, LinkEndId("Earth", vlbi_station_name)}, // FIXME!
        };
        for (size_t i = 0; i < numDataRows; ++i) {
          linkEndsVector_.push_back(constantLinkEnds);
        }
      } else if (dataMap.count(input_output::TrackingDataType::dsn_receiving_station_nr)) {
        for (size_t i = 0; i < numDataRows; ++i) {
          std::string vlbi_station_name = metaDataStrMap.at(input_output::TrackingDataType::vlbi_station_name);
          LinkEnds currentLinkEnds{
              {transmitter, LinkEndId(spacecraftName_, "Antenna")},
              {receiver, LinkEndId("Earth", vlbi_station_name)}, // FIXME!
          };
          linkEndsVector_.push_back(currentLinkEnds);
        }
      }
      break;
    }

    default: {
      throw std::runtime_error("Error while processing tracking txt file: LinkEnds representation not recognised or implemented.");
    }
  }

  // Creating a set with all the various linkEnds
  linkEndsSet_ = utilities::vectorToSet(linkEndsVector_);
}

ProcessedTrackingTxtFileContents::TimeRepresentation ProcessedTrackingTxtFileContents::getTimeRepresentation()
{
  auto const& availableDataTypes = rawTrackingTxtFileContents_->getDataColumnTypes();

  if (utilities::containsAll(availableDataTypes, std::vector<input_output::TrackingDataType>{input_output::TrackingDataType::tdb_time_j2000})) {
    return tdb_seconds_j2000;
  }

  if (utilities::containsAll(availableDataTypes,
                             std::vector<input_output::TrackingDataType>{
                                 input_output::TrackingDataType::year,
                                 input_output::TrackingDataType::month,
                                 input_output::TrackingDataType::day,
                                 input_output::TrackingDataType::hour,
                                 input_output::TrackingDataType::minute,
                                 input_output::TrackingDataType::second
                             })) {
    return calendar_day_time;
  }
  throw std::runtime_error("Error while processing tracking txt file: Time representation not recognised or implemented.");
}

ProcessedTrackingTxtFileContents::LinkEndsRepresentation ProcessedTrackingTxtFileContents::getLinkEndsRepresentation()
{
  auto const& availableDataTypes = rawTrackingTxtFileContents_->getAllAvailableDataTypes();

  if (utilities::containsAll(availableDataTypes,
                             std::vector<input_output::TrackingDataType>{
                                 input_output::TrackingDataType::dsn_transmitting_station_nr,
                                 input_output::TrackingDataType::dsn_receiving_station_nr
                             })) {
    return dsn_transmitting_receiving_station_nr;
  }

  if (utilities::containsAll(availableDataTypes, std::vector<input_output::TrackingDataType>{input_output::TrackingDataType::vlbi_station_name})) {
    return vlbi_station;
  }
  throw std::runtime_error("Error while processing tracking txt file: Link Ends representation not recognised or implemented.");
}

} // namespace observation_models
} // namespace tudat