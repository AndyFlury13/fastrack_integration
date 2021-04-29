// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/detail/Axis.hpp"

#include <memory>

template <typename SpacePoint>
std::unique_ptr<Acts::SpacePointGrid<SpacePoint>>
Acts::SpacePointGridCreator::createGrid(
    const Acts::SpacePointGridConfig& config) {
  // calculate circle intersections of helix and max detector radius
  float minHelixRadius = config.minPt / (300. * config.bFieldInZ);  // in mm
  float maxR2 = config.rMax * config.rMax;
  float xOuter = maxR2 / (2 * minHelixRadius);
  float yOuter = std::sqrt(maxR2 - xOuter * xOuter);
  float outerAngle = std::atan(xOuter / yOuter);
  // intersection of helix and max detector radius minus maximum R distance from
  // middle SP to top SP
  float innerAngle = 0;
  if (config.rMax > config.deltaRMax) {
    float innerCircleR2 =
        (config.rMax - config.deltaRMax) * (config.rMax - config.deltaRMax);
    float xInner = innerCircleR2 / (2 * minHelixRadius);
    float yInner = std::sqrt(innerCircleR2 - xInner * xInner);
    innerAngle = std::atan(xInner / yInner);
  }

  // FIXME: phibin size must include max impact parameters
  // divide 2pi by angle delta to get number of phi-bins
  // size is always 2pi even for regions of interest
  int phiBins = std::floor(2 * M_PI / (outerAngle - innerAngle));
  Acts::detail::Axis<detail::AxisType::Equidistant,
                     detail::AxisBoundaryType::Closed>
      phiAxis(-M_PI, M_PI, phiBins);

  // TODO: can probably be optimized using smaller z bins
  // and returning (multiple) neighbors only in one z-direction for forward
  // seeds
  // FIXME: zBinSize must include scattering
  float zBinSize = config.cotThetaMax * config.deltaRMax;
  int zBins = std::floor((config.zMax - config.zMin) / zBinSize);
  detail::Axis<detail::AxisType::Equidistant, detail::AxisBoundaryType::Bound>
      zAxis(config.zMin, config.zMax, zBins);
  return std::make_unique<Acts::SpacePointGrid<SpacePoint>>(
      std::make_tuple(phiAxis, zAxis));
}
