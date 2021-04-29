// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/Framework/WriterT.hpp"

#include <memory>
#include <string>

namespace ActsExamples {

/// Write track finder performance measures.
///
/// Only considers the track finding itself, i.e. grouping of hits into tracks,
/// and computes relevant per-track and per-particles statistics.
class TrackFinderPerformanceWriter final : public WriterT<ProtoTrackContainer> {
 public:
  struct Config {
    /// Input reconstructed proto tracks collection.
    std::string inputProtoTracks;
    /// Input hit-particles map collection.
    std::string inputMeasurementParticlesMap;
    /// Input particles collection.
    std::string inputParticles;
    /// Output directory.
    std::string outputDir;
    /// Output filename.
    std::string outputFilename = "performance_track_finder.root";
  };

  TrackFinderPerformanceWriter(Config cfg, Acts::Logging::Level lvl);
  ~TrackFinderPerformanceWriter() final override;

  ProcessCode endRun() final override;

 private:
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const ProtoTrackContainer& tracks) final override;

  struct Impl;
  std::unique_ptr<Impl> m_impl;
};

}  // namespace ActsExamples
