// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/TrackParamsEstimationAlgorithm.hpp"

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Seeding/EstimateTrackParamsFromSeed.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <map>
#include <stdexcept>

ActsExamples::TrackParamsEstimationAlgorithm::TrackParamsEstimationAlgorithm(
    ActsExamples::TrackParamsEstimationAlgorithm::Config cfg,
    Acts::Logging::Level lvl)
    : ActsExamples::BareAlgorithm("TrackParamsEstimationAlgorithm", lvl),
      m_cfg(std::move(cfg)) {
  // Either seeds directly or proto tracks + space points
  if (m_cfg.inputSeeds.empty() and m_cfg.inputSpacePoints.empty()) {
    throw std::invalid_argument("Missing seeds or space point collection");
  }
  if (not m_cfg.inputSpacePoints.empty()) {
    for (const auto& i : m_cfg.inputSpacePoints) {
      if (i.empty()) {
        throw std::invalid_argument("Invalid space point input collection");
      }
    }
    if (m_cfg.inputProtoTracks.empty()) {
      throw std::invalid_argument("Missing proto tracks input collections");
    }
  }
  if (m_cfg.inputSourceLinks.empty()) {
    throw std::invalid_argument("Missing source links input collection");
  }
  if (m_cfg.outputTrackParameters.empty()) {
    throw std::invalid_argument("Missing track parameters output collection");
  }
  if (m_cfg.outputProtoTracks.empty()) {
    throw std::invalid_argument("Missing proto tracks output collections");
  }
  if (not m_cfg.trackingGeometry) {
    throw std::invalid_argument("Missing tracking geometry");
  }

  // Set up the track parameters covariance (the same for all tracks)
  m_covariance(Acts::eBoundLoc0, Acts::eBoundLoc0) =
      cfg.sigmaLoc0 * m_cfg.sigmaLoc0;
  m_covariance(Acts::eBoundLoc1, Acts::eBoundLoc1) =
      cfg.sigmaLoc1 * m_cfg.sigmaLoc1;
  m_covariance(Acts::eBoundPhi, Acts::eBoundPhi) =
      cfg.sigmaPhi * m_cfg.sigmaPhi;
  m_covariance(Acts::eBoundTheta, Acts::eBoundTheta) =
      cfg.sigmaTheta * m_cfg.sigmaTheta;
  m_covariance(Acts::eBoundQOverP, Acts::eBoundQOverP) =
      cfg.sigmaQOverP * m_cfg.sigmaQOverP;
  m_covariance(Acts::eBoundTime, Acts::eBoundTime) =
      m_cfg.sigmaT0 * m_cfg.sigmaT0;
}

ActsExamples::SimSeedContainer
ActsExamples::TrackParamsEstimationAlgorithm::createSeeds(
    const ActsExamples::ProtoTrackContainer& protoTracks,
    const ActsExamples::SimSpacePointContainer& spacePoints) const {
  SimSeedContainer seeds;
  seeds.reserve(protoTracks.size());
  for (std::size_t itrack = 0; itrack < protoTracks.size(); ++itrack) {
    // The list of hits and the initial start parameters
    const auto& protoTrack = protoTracks[itrack];
    if (protoTrack.size() < 3) {
      ACTS_WARNING("Proto track " << itrack << " size is less than 3.");
      continue;
    }
    // Space points on the proto track
    std::vector<const SimSpacePoint*> spacePointsOnTrack;
    spacePointsOnTrack.reserve(protoTrack.size());
    // Loop over the hit index on the proto track to find the space points
    for (const auto& hitIndex : protoTrack) {
      auto it =
          std::find_if(spacePoints.begin(), spacePoints.end(),
                       [&](const SimSpacePoint& spacePoint) {
                         return (spacePoint.measurementIndex() == hitIndex);
                       });
      if (it != spacePoints.end()) {
        spacePointsOnTrack.push_back(&(*it));
      }
    }
    // At least three space points are required
    if (spacePointsOnTrack.size() < 3) {
      continue;
    }
    // Sort the space points
    std::sort(spacePointsOnTrack.begin(), spacePointsOnTrack.end(),
              [](const SimSpacePoint* lhs, const SimSpacePoint* rhs) {
                return std::hypot(lhs->r(), lhs->z()) <
                       std::hypot(rhs->r(), rhs->z());
              });
    // Loop over the found space points to find seeds with simple selection
    // @todo add the check of deltaZ
    for (size_t ib = 0; ib < spacePointsOnTrack.size() - 2; ++ib) {
      for (size_t im = ib + 1; im < spacePointsOnTrack.size() - 1; ++im) {
        for (size_t it = im + 1; it < spacePointsOnTrack.size(); ++it) {
          double bmDeltaR = std::abs(spacePointsOnTrack[im]->r() -
                                     spacePointsOnTrack[ib]->r());
          double mtDeltaR = std::abs(spacePointsOnTrack[it]->r() -
                                     spacePointsOnTrack[im]->r());
          if (bmDeltaR >= m_cfg.deltaRMin and bmDeltaR <= m_cfg.deltaRMax and
              mtDeltaR >= m_cfg.deltaRMin and mtDeltaR <= m_cfg.deltaRMax) {
            seeds.emplace_back(*spacePointsOnTrack[ib], *spacePointsOnTrack[im],
                               *spacePointsOnTrack[it],
                               spacePointsOnTrack[im]->z());
          }
        }
      }
    }
  }
  return seeds;
}

ActsExamples::ProcessCode ActsExamples::TrackParamsEstimationAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // Read source links (necesary for retrieving the geometry identifer)
  const auto& sourceLinks =
      ctx.eventStore.get<IndexSourceLinkContainer>(m_cfg.inputSourceLinks);
  // Read seeds or create them from proto tracks and space points
  SimSeedContainer seeds;
  SimSpacePointContainer spacePoints;
  if (not m_cfg.inputSeeds.empty()) {
    seeds = ctx.eventStore.get<SimSeedContainer>(m_cfg.inputSeeds);
  } else {
    const auto& protoTracks =
        ctx.eventStore.get<ProtoTrackContainer>(m_cfg.inputProtoTracks);
    for (const auto& isp : m_cfg.inputSpacePoints) {
      const auto& sps = ctx.eventStore.get<SimSpacePointContainer>(isp);
      std::copy(sps.begin(), sps.end(), std::back_inserter(spacePoints));
    }
    seeds = createSeeds(protoTracks, spacePoints);
  }

  TrackParametersContainer trackParameters;
  ProtoTrackContainer tracks;
  trackParameters.reserve(seeds.size());
  tracks.reserve(seeds.size());

  // Loop over all found seeds to estimate track parameters
  for (size_t iseed = 0; iseed < seeds.size(); ++iseed) {
    const auto& seed = seeds[iseed];
    // Get the bottom space point and its reference surface
    const auto bottomSP = seed.sp().front();
    const auto hitIdx = bottomSP->measurementIndex();
    const auto sourceLink = sourceLinks.nth(hitIdx);
    auto geoId = sourceLink->geometryId();
    const Acts::Surface* surface = m_cfg.trackingGeometry->findSurface(geoId);
    if (surface == nullptr) {
      ACTS_WARNING("surface with geoID "
                   << geoId << " is not found in the tracking gemetry");
      continue;
    }

    // Get the magnetic field at the bottom space point
    Acts::Vector3 field = m_cfg.magneticField->getField(
        {bottomSP->x(), bottomSP->y(), bottomSP->z()});
    // Estimate the track parameters from seed
    auto optParams = Acts::estimateTrackParamsFromSeed(
        ctx.geoContext, seed.sp().begin(), seed.sp().end(), *surface, field,
        m_cfg.bFieldMin);
    if (not optParams.has_value()) {
      ACTS_WARNING("Estimation of track parameters for seed " << iseed
                                                              << " failed.");
      continue;
    } else {
      const auto& params = optParams.value();
      double charge = std::copysign(1, params[Acts::eBoundQOverP]);
      trackParameters.emplace_back(surface->getSharedPtr(), params, charge,
                                   m_covariance);
      // Create a proto track for this seed
      ProtoTrack track;
      track.reserve(3);
      for (const auto& sp : seed.sp()) {
        track.push_back(sp->measurementIndex());
      }
      tracks.emplace_back(track);
    }
  }

  ctx.eventStore.add(m_cfg.outputTrackParameters, std::move(trackParameters));
  ctx.eventStore.add(m_cfg.outputProtoTracks, std::move(tracks));
  return ActsExamples::ProcessCode::SUCCESS;
}
