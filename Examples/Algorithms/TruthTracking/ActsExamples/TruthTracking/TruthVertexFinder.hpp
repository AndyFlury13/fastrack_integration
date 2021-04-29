// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Framework/BareAlgorithm.hpp"

#include <string>

namespace ActsExamples {

/// Group particles into proto vertices using truth information.
class TruthVertexFinder final : public BareAlgorithm {
 public:
  struct Config {
    /// The input truth particles that should be used to create proto vertices.
    std::string inputParticles;
    /// The output proto vertices collection.
    std::string outputProtoVertices;
    /// Exclude secondary particles not originating from the primary vertex.
    bool excludeSecondaries = false;
    /// Build separate proto vertices for the secondary particles.
    bool separateSecondaries = false;
  };

  TruthVertexFinder(const Config& cfg, Acts::Logging::Level lvl);

  ProcessCode execute(const AlgorithmContext& ctx) const final;

 private:
  Config m_cfg;
};

}  // namespace ActsExamples
