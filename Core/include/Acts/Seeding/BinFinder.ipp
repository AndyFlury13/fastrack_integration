// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
template <typename external_spacepoint_t>
std::vector<size_t> Acts::BinFinder<external_spacepoint_t>::findBins(
    size_t phiBin, size_t zBin,
    const Acts::SpacePointGrid<external_spacepoint_t>* binnedSP) {
  return binnedSP->neighborHoodIndices({phiBin, zBin}).collect();
}
