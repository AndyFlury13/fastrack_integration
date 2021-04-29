// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Framework/BareAlgorithm.hpp"

namespace ActsExamples {

/// Select truth particles to be used as 'seeds' of reconstruction algorithms,
/// e.g. track fitting and track finding.
///
/// This pre-selection could help guarantee quality of the 'seeds', i.e. to
/// avoid empty proto track (no recorded hits for the particle). In addition, it
/// could help save unnecessary reconstruction time. For instance, when
/// investigating performance of CombinatorialKalmanFilter (CKF), we might be
/// interested in its performance for only truth particles with pT and number of
/// recorded hits (on sensitive detectors) safistying provided criteria (input
/// measurments of CKF are still recorded hits from all possible particles).
/// Then we could use particles only satistying provided criteria as the 'seeds'
/// of CKF instead of handling all the truth particles.
//
class TruthSeedSelector final : public BareAlgorithm {
 public:
  struct Config {
    /// The input truth particles that should be used to create proto tracks.
    std::string inputParticles;
    /// The input hit-particles map collection.
    std::string inputMeasurementParticlesMap;
    /// The output proto tracks collection.
    std::string outputParticles;
    /// Maximum distance from the origin in the transverse plane
    double rhoMax = std::numeric_limits<double>::max();
    /// Maximum absolute distance from the origin along z
    double absZMax = std::numeric_limits<double>::max();
    // Truth particle kinematic cuts
    double phiMin = std::numeric_limits<double>::lowest();
    double phiMax = std::numeric_limits<double>::max();
    double etaMin = std::numeric_limits<double>::lowest();
    double etaMax = std::numeric_limits<double>::max();
    double absEtaMin = std::numeric_limits<double>::lowest();
    double absEtaMax = std::numeric_limits<double>::max();
    double ptMin = 0.0;
    double ptMax = std::numeric_limits<double>::max();
    /// Keep neutral particles
    bool keepNeutral = false;
    /// Requirement on number of recorded hits
    //@TODO: implement detector-specific requirements
    size_t nHitsMin = 0;
    size_t nHitsMax = std::numeric_limits<size_t>::max();
  };

  TruthSeedSelector(const Config& cfg, Acts::Logging::Level lvl);

  ProcessCode execute(const AlgorithmContext& ctx) const override final;

 private:
  Config m_cfg;
};

}  // namespace ActsExamples
