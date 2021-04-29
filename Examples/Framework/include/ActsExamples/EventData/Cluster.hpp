// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsFatras/Digitization/Channelizer.hpp"

#include <vector>

namespace ActsExamples {

/// Simple struct holding cluster information.
struct Cluster {
  size_t sizeLoc0 = 0;
  size_t sizeLoc1 = 0;
  std::vector<ActsFatras::Channelizer::ChannelSegment> channels;
};

/// Clusters have a one-to-one relation with measurements
using ClusterContainer = std::vector<Cluster>;

}  // namespace ActsExamples
