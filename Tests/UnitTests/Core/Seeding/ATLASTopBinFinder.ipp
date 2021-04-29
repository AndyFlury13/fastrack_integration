// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// DEBUG: THIS REQUIRES THE BINS TO BE SET TO phi:41 z:11
template <typename SpacePoint>
std::set<size_t> Acts::ATLASTopBinFinder<SpacePoint>::findBins(
    size_t phiBin, size_t zBin,
    const Acts::SpacePointGrid<SpacePoint>* binnedSP) {
  std::set<size_t> neighbourBins =
      binnedSP->neighborHoodIndices({phiBin, zBin}, 1);
  if (zBin == 6) {
    return neighbourBins;
  }
  // <10 not necessary because grid doesn't return bins that don't exist (only
  // works for 11 zbins)
  if (zBin > 6) {
    neighbourBins.erase(binnedSP->getGlobalBinIndex({phiBin, zBin - 1}));
    neighbourBins.erase(binnedSP->getGlobalBinIndex({phiBin - 1, zBin - 1}));
    neighbourBins.erase(binnedSP->getGlobalBinIndex({phiBin + 1, zBin - 1}));
  } else {
    neighbourBins.erase(binnedSP->getGlobalBinIndex({phiBin, zBin + 1}));
    neighbourBins.erase(binnedSP->getGlobalBinIndex({phiBin - 1, zBin + 1}));
    neighbourBins.erase(binnedSP->getGlobalBinIndex({phiBin + 1, zBin + 1}));
  }
  return neighbourBins;
}
