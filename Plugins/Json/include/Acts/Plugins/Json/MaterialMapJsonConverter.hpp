// This file is part of the Acts project.
//
// Copyright (C) 2017-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/IVolumeMaterial.hpp"
#include "Acts/Plugins/Json/GeometryHierarchyMapJsonConverter.hpp"
#include "Acts/Utilities/Logger.hpp"
#include <Acts/Geometry/TrackingVolume.hpp>
#include <Acts/Surfaces/Surface.hpp>

#include <map>

#include <nlohmann/json.hpp>

namespace Acts {

using SurfaceAndMaterial =
    std::pair<std::shared_ptr<const Acts::Surface>,
              std::shared_ptr<const Acts::ISurfaceMaterial>>;
using TrackingVolumeAndMaterial =
    std::pair<const Acts::TrackingVolume*,
              std::shared_ptr<const Acts::IVolumeMaterial>>;

/// @class MaterialMapJsonConverter
///
/// @brief read the material from Json
class MaterialMapJsonConverter {
 public:
  using SurfaceMaterialMap =
      std::map<GeometryIdentifier, std::shared_ptr<const ISurfaceMaterial>>;
  using VolumeMaterialMap =
      std::map<GeometryIdentifier, std::shared_ptr<const IVolumeMaterial>>;
  using DetectorMaterialMaps = std::pair<SurfaceMaterialMap, VolumeMaterialMap>;

  /// @class Config
  /// Configuration of the Converter
  class Config {
   public:
    /// The default logger
    std::shared_ptr<const Logger> logger;
    /// Default geometry context to extract surface tranforms
    GeometryContext context = GeometryContext();
    /// The name of the writer
    std::string name = "";

    /// Steering to handle sensitive data
    bool processSensitives = true;
    /// Steering to handle approach data
    bool processApproaches = true;
    /// Steering to handle representing data
    bool processRepresenting = true;
    /// Steering to handle boundary data
    bool processBoundaries = true;
    /// Steering to handle volume data
    bool processVolumes = true;
    /// Steering to handle volume data
    bool processDenseVolumes = false;
    /// Add proto material to all surfaces
    bool processNonMaterial = false;

    /// Constructor
    ///
    /// @param lname Name of the writer tool
    /// @param lvl The output logging level
    Config(const std::string& lname = "MaterialMapJsonConverter",
           Logging::Level lvl = Logging::INFO)
        : logger(getDefaultLogger(lname, lvl)), name(lname) {}
  };

  /// Constructor
  ///
  /// @param cfg configuration struct for the reader
  MaterialMapJsonConverter(const Config& cfg);

  /// Destructor
  ~MaterialMapJsonConverter() = default;

  /// Convert a json material map to a DetectorMaterialMaps
  ///
  /// @param materialmaps The json material
  DetectorMaterialMaps jsonToMaterialMaps(const nlohmann::json& materialmaps);

  /// Convert a DetectorMaterialMaps to json
  ///
  /// @param maps The material map collection
  nlohmann::ordered_json materialMapsToJson(const DetectorMaterialMaps& maps);

  /// Convert a tracking geometry to json.
  /// Can be used to initialise the material mapping process.
  ///
  /// @param tGeometry is the tracking geometry
  nlohmann::ordered_json trackingGeometryToJson(
      const TrackingGeometry& tGeometry);

  /// Go through a volume to find subvolume, layers and surfaces.
  /// Store volumes and surfaces in two vector used to initialised the geometry
  /// hierachy.
  ///
  /// @param volumeHierarchy is a vector of volume to be filled
  /// @param surfaceHierarchy is a vector of surfaces to be filled
  /// @param tVolume is a volume
  void convertToHierarchy(
      std::vector<std::pair<GeometryIdentifier,
                            Acts::TrackingVolumeAndMaterial>>& volumeHierarchy,
      std::vector<std::pair<GeometryIdentifier, Acts::SurfaceAndMaterial>>&
          surfaceHierarchy,
      const Acts::TrackingVolume* tVolume);

 private:
  /// The config class
  Config m_cfg;

  /// Name of the volume hierarchy
  std::string m_volumeName = "Material Volume Map";
  /// Geometry hierarchy writer for volume material.
  Acts::GeometryHierarchyMapJsonConverter<const IVolumeMaterial*>
      m_volumeMaterialConverter;
  /// Geometry hierarchy writer for tracking volume.
  Acts::GeometryHierarchyMapJsonConverter<Acts::TrackingVolumeAndMaterial>
      m_volumeConverter;

  /// Name of the surface hierarchy
  std::string m_surfaceName = "Material Surface Map";
  /// Geometry hierarchy writer for surface material.
  Acts::GeometryHierarchyMapJsonConverter<const ISurfaceMaterial*>
      m_surfaceMaterialConverter;
  /// Geometry hierarchy writer for surface.
  Acts::GeometryHierarchyMapJsonConverter<Acts::SurfaceAndMaterial>
      m_surfaceConverter;

  /// Private access to the logging instance
  const Logger& logger() const { return *m_cfg.logger; }
};

}  // namespace Acts
