// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Material/AccumulatedVolumeMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/detail/AxisFwd.hpp"
#include "Acts/Utilities/detail/Grid.hpp"

#include <array>
#include <functional>
#include <utility>
#include <vector>

namespace Acts {

class MaterialSlab;

/// list of point used in the mapping of a volume
using RecordedMaterialVolumePoint =
    std::vector<std::pair<Acts::MaterialSlab, std::vector<Acts::Vector3>>>;

using EAxis = Acts::detail::EquidistantAxis;
using Grid2D =
    Acts::detail::Grid<Acts::AccumulatedVolumeMaterial, EAxis, EAxis>;
using Grid3D =
    Acts::detail::Grid<Acts::AccumulatedVolumeMaterial, EAxis, EAxis, EAxis>;
using MaterialGrid2D =
    Acts::detail::Grid<Acts::Material::ParametersVector, EAxis, EAxis>;
using MaterialGrid3D =
    Acts::detail::Grid<Acts::Material::ParametersVector, EAxis, EAxis, EAxis>;

/// @brief Helper method that creates the cache grid for the mapping. This
/// grid allows the collection of material at a the anchor points.
///
/// @param [in] gridAxis1 Axis data
/// @param [in] gridAxis2 Axis data
/// @note The data of the axes is given in the std::array as {minimum value,
/// maximum value, number of bins}
///
/// @return The grid
Grid2D createGrid(std::array<double, 3> gridAxis1,
                  std::array<double, 3> gridAxis2);

/// @brief Helper method that creates the cache grid for the mapping. This
/// grid allows the collection of material at a the anchor points.
///
/// @param [in] gridAxis1 Axis data
/// @param [in] gridAxis2 Axis data
/// @param [in] gridAxis3 Axis data
/// @note The data of the axes is given in the std::array as {minimum value,
/// maximum value, number of bins}
///
/// @return The grid
Grid3D createGrid(std::array<double, 3> gridAxis1,
                  std::array<double, 3> gridAxis2,
                  std::array<double, 3> gridAxis3);

/// @brief return a function that return the coordinate corresponding to type of
/// bin
///
/// @param [in] Type of bin
///
/// @return a coordinate transform function
std::function<double(Acts::Vector3)> globalToLocalFromBin(
    Acts::BinningValue& type);

/// @brief Create a 2DGrid using a BinUtility.
/// Also determine the coresponding global to local transform and grid mapping
/// function
///
/// @param [in] BinUtility of the volume to be mapped
/// @param [in] Global to local transform to be updated.
///
/// @return the 3D grid
Grid2D createGrid2D(
    const BinUtility& bins,
    std::function<Acts::Vector2(Acts::Vector3)>& transfoGlobalToLocal);

/// @brief Create a 3DGrid using a BinUtility.
/// Also determine the coresponding global to local transform and grid mapping
/// function
///
/// @param [in] BinUtility of the volume to be mapped
/// @param [in] Global to local transform to be updated.
///
/// @return the 3D grid
Grid3D createGrid3D(
    const BinUtility& bins,
    std::function<Acts::Vector3(Acts::Vector3)>& transfoGlobalToLocal);

/// @brief Concatenate a set of material at arbitrary space points on a set of
/// grid points and produces a grid containing the averaged material values.
///
/// @param [in] grid The material collecting grid
/// @param [in] mPoints The set of material at the space points
/// @param [in] transfoGlobalToLocal tranformation from local to local
/// coordinate
///
/// @return The average material grid decomposed into classification numbers
MaterialGrid2D mapMaterialPoints(
    Grid2D& grid, const Acts::RecordedMaterialVolumePoint& mPoints,
    std::function<Acts::Vector2(Acts::Vector3)>& transfoGlobalToLocal);

/// @brief Concatenate a set of material at arbitrary space points on a set of
/// grid points and produces a grid containing the averaged material values.
///
/// @param [in] grid The material collecting grid
/// @param [in] mPoints The set of material at the space points
/// @param [in] transfoGlobalToLocal tranformation from local to local
/// coordinate
///
/// @return The average material grid decomposed into classification numbers
MaterialGrid3D mapMaterialPoints(
    Grid3D& grid, const Acts::RecordedMaterialVolumePoint& mPoints,
    std::function<Acts::Vector3(Acts::Vector3)>& transfoGlobalToLocal);

}  // namespace Acts
