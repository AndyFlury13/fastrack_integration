// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/DD4hep/ConvertDD4hepDetector.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/PassiveLayerBuilder.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/Geometry/TrackingVolumeArrayCreator.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Plugins/DD4hep/ActsExtension.hpp"
#include "Acts/Plugins/DD4hep/ConvertDD4hepMaterial.hpp"
#include "Acts/Plugins/DD4hep/DD4hepLayerBuilder.hpp"
#include "Acts/Plugins/DD4hep/DD4hepVolumeBuilder.hpp"
#include "Acts/Utilities/BinningData.hpp"

#include <list>
#include <stdexcept>

#include "TGeoManager.h"

namespace Acts {
std::unique_ptr<const TrackingGeometry> convertDD4hepDetector(
    dd4hep::DetElement worldDetElement, Logging::Level loggingLevel,
    BinningType bTypePhi, BinningType bTypeR, BinningType bTypeZ,
    double layerEnvelopeR, double layerEnvelopeZ, double defaultLayerThickness,
    const std::function<void(std::vector<dd4hep::DetElement>& detectors)>&
        sortSubDetectors,
    const Acts::GeometryContext& gctx,
    std::shared_ptr<const IMaterialDecorator> matDecorator) {
  // create local logger for conversion
  ACTS_LOCAL_LOGGER(Acts::getDefaultLogger("DD4hepConversion", loggingLevel));
  ACTS_INFO("Translating DD4hep geometry into Acts geometry");
  // get the sub detectors of the world detector e.g. beampipe, pixel detector,
  // strip detector
  std::vector<dd4hep::DetElement> subDetectors;
  // go through the detector hierarchies
  collectSubDetectors_dd4hep(worldDetElement, subDetectors);
  // sort to build detector from bottom to top
  sortSubDetectors(subDetectors);
  // the volume builders of the subdetectors
  std::list<std::shared_ptr<const ITrackingVolumeBuilder>> volumeBuilders;
  // the beam pipe volume builder needs special treatment and needs to be added
  // in the end (beampipe exceeds length of all other subdetectors)
  std::shared_ptr<const CylinderVolumeBuilder> beamPipeVolumeBuilder;
  // loop over the sub detectors
  for (auto& subDetector : subDetectors) {
    ACTS_INFO("Translating DD4hep sub detector: " << subDetector.name());
    // create volume builder
    auto volBuilder = volumeBuilder_dd4hep(
        subDetector, loggingLevel, bTypePhi, bTypeR, bTypeZ, layerEnvelopeR,
        layerEnvelopeZ, defaultLayerThickness);
    if (volBuilder) {
      // distinguish beam pipe
      if (volBuilder->getConfiguration().buildToRadiusZero) {
        // check if beam pipe is already present
        if (beamPipeVolumeBuilder) {
          throw std::logic_error(
              std::string("Beampipe has already been set! There can only "
                          "exist one beam pipe. Please check your "
                          "detector construction. Current volume name: ") +
              volBuilder->getConfiguration().volumeName +
              std::string(", name of volume, already set as beam pipe: ") +
              beamPipeVolumeBuilder->getConfiguration().volumeName);
        }
        // set the beam pipe
        beamPipeVolumeBuilder = volBuilder;
      } else {
        volumeBuilders.push_back(volBuilder);
      }
    }
  }
  // Finally add the beam pipe
  if (beamPipeVolumeBuilder) {
    volumeBuilders.push_back(beamPipeVolumeBuilder);
  }

  std::vector<std::function<std::shared_ptr<TrackingVolume>(
      const GeometryContext&, const TrackingVolumePtr&,
      const VolumeBoundsPtr&)>>
      volumeFactories;

  for (const auto& vb : volumeBuilders) {
    volumeFactories.push_back(
        [vb](const GeometryContext& vgctx,
             const std::shared_ptr<const TrackingVolume>& inner,
             const VolumeBoundsPtr&) {
          return vb->trackingVolume(vgctx, inner);
        });
  }

  // create cylinder volume helper
  auto volumeHelper = cylinderVolumeHelper_dd4hep();
  // hand over the collected volume builders
  Acts::TrackingGeometryBuilder::Config tgbConfig;
  tgbConfig.trackingVolumeHelper = volumeHelper;
  tgbConfig.materialDecorator = std::move(matDecorator);
  tgbConfig.trackingVolumeBuilders = std::move(volumeFactories);
  auto trackingGeometryBuilder =
      std::make_shared<const Acts::TrackingGeometryBuilder>(tgbConfig);
  return (trackingGeometryBuilder->trackingGeometry(gctx));
}

std::shared_ptr<const CylinderVolumeBuilder> volumeBuilder_dd4hep(
    dd4hep::DetElement subDetector, Logging::Level loggingLevel,
    BinningType bTypePhi, BinningType bTypeR, BinningType bTypeZ,
    double layerEnvelopeR, double layerEnvelopeZ,
    double defaultLayerThickness) {
  // create cylinder volume helper
  auto volumeHelper = cylinderVolumeHelper_dd4hep(loggingLevel);
  // create local logger for conversion
  ACTS_LOCAL_LOGGER(Acts::getDefaultLogger("D2A_Logger", loggingLevel));
  ACTS_VERBOSE("Processing detector element:  " << subDetector.name());

  Acts::ActsExtension* subDetExtension = nullptr;
  // at this stage not every DetElement needs to have an Extension attached
  try {
    subDetExtension = subDetector.extension<Acts::ActsExtension>();
  } catch (std::runtime_error& e) {
  }
  if (subDetector.type() == "compound") {
    ACTS_VERBOSE("Subdetector : '"
                 << subDetector.name()
                 << "' has no ActsExtension and has type compound ");
    ACTS_VERBOSE(
        "handling as a compound volume (a hierachy of a "
        "barrel-endcap structure) and resolving the "
        "subvolumes...");
    // Now create the Layerbuilders and Volumebuilder
    // the layers
    /// the dd4hep::DetElements of the layers of the negative volume
    std::vector<dd4hep::DetElement> negativeLayers;
    /// the dd4hep::DetElements of the layers of the central volume
    std::vector<dd4hep::DetElement> centralLayers;
    /// the dd4hep::DetElements of the layers of the positive volume
    std::vector<dd4hep::DetElement> positiveLayers;

    // the configuration object of the volume builder
    Acts::CylinderVolumeBuilder::Config cvbConfig;

    // go through sub volumes
    std::vector<dd4hep::DetElement> compounds;
    collectCompounds_dd4hep(subDetector, compounds);

    // get z position to distinguish positive & negative endcap
    double zPos = 0.;
    // flags to catch if sub volumes have been set already
    bool nEndCap = false;
    bool pEndCap = false;
    bool barrel = false;
    for (auto& volumeDetElement : compounds) {
      ACTS_VERBOSE("Volume : '"
                   << subDetector.name()
                   << "'is a compound volume -> resolve now the sub volumes");

      // get the dimensions of the volume
      TGeoShape* geoShape =
          volumeDetElement.placement().ptr()->GetVolume()->GetShape();
      // check if it has a shape (the other case should not happen)
      if (geoShape != nullptr) {
        zPos = volumeDetElement.placement()
                   .ptr()
                   ->GetMatrix()
                   ->GetTranslation()[2] *
               UnitConstants::cm;
      } else {
        throw std::logic_error(std::string("Volume of DetElement: ") +
                               volumeDetElement.name() +
                               std::string(" has no shape!"));
      }
      // check if it has a volume extension telling if it is a barrel or an
      // endcap
      ActsExtension* volumeExtension = nullptr;
      try {
        volumeExtension = volumeDetElement.extension<ActsExtension>();
      } catch (std::runtime_error& e) {
        throw std::logic_error(
            std::string("Current DetElement: ") + volumeDetElement.name() +
            std::string(" has no ActsExtension! At this stage it should be a "
                        "detector volume declared as Barrel or Endcap. Please"
                        "check your detector construction."));
      }

      if (volumeExtension->hasType("endcap", "detector")) {
        ACTS_VERBOSE(
            std::string("Subvolume : '") + volumeDetElement.name() +
            std::string("' is a disc volume -> handling as an endcap"));
        if (zPos < 0.) {
          if (nEndCap) {
            throw std::logic_error(
                "Negative Endcap was already given for this "
                "hierachy! Please create a new "
                "DD4hep_SubDetectorAssembly for the next "
                "hierarchy.");
          }
          nEndCap = true;
          ACTS_VERBOSE("      ->is negative endcap");
          collectLayers_dd4hep(volumeDetElement, negativeLayers);
          // Fill the volume material for barrel case
          if (volumeExtension->hasType("boundary_material")) {
            if (volumeExtension->hasValue("boundary_material_negative")) {
              cvbConfig.boundaryMaterial[2] = Acts::createProtoMaterial(
                  *volumeExtension, "boundary_material_negative",
                  {{"binPhi", Acts::closed}, {"binR", Acts::open}});
            }
            if (volumeExtension->hasValue("boundary_material_positive")) {
              cvbConfig.boundaryMaterial[3] = Acts::createProtoMaterial(
                  *volumeExtension, "boundary_material_positive",
                  {{"binPhi", Acts::closed}, {"binR", Acts::open}});
            }
          }
        } else {
          if (pEndCap) {
            throw std::logic_error(
                "Positive Endcap was already given for this "
                "hierachy! Please create a new "
                "DD4hep_SubDetectorAssembly for the next "
                "hierarchy.");
          }
          pEndCap = true;
          ACTS_VERBOSE("      ->is positive endcap");
          collectLayers_dd4hep(volumeDetElement, positiveLayers);
          // Fill the volume material for barrel case
          if (volumeExtension->hasType("boundary_material")) {
            if (volumeExtension->hasValue("boundary_material_negative")) {
              cvbConfig.boundaryMaterial[4] = Acts::createProtoMaterial(
                  *volumeExtension, "boundary_material_negative",
                  {{"binPhi", Acts::closed}, {"binR", Acts::open}});
            }
            if (volumeExtension->hasValue("boundary_material_positive")) {
              cvbConfig.boundaryMaterial[5] = Acts::createProtoMaterial(
                  *volumeExtension, "boundary_material_positive",
                  {{"binPhi", Acts::closed}, {"binR", Acts::open}});
            }
          }
        }
      } else if (volumeExtension->hasType("barrel", "detector")) {
        if (barrel) {
          throw std::logic_error(
              "Barrel was already given for this "
              "hierachy! Please create a new "
              "DD4hep_SubDetectorAssembly for the next "
              "hierarchy.");
        }
        barrel = true;
        ACTS_VERBOSE("Subvolume : "
                     << volumeDetElement.name()
                     << " is a cylinder volume -> handling as a barrel");
        collectLayers_dd4hep(volumeDetElement, centralLayers);
        // Fill the volume material for barrel case
        if (volumeExtension->hasType("boundary_material")) {
          if (volumeExtension->hasValue("boundary_material_negative")) {
            cvbConfig.boundaryMaterial[3] = Acts::createProtoMaterial(
                *volumeExtension, "boundary_material_negative",
                {{"binPhi", Acts::closed}, {"binR", Acts::open}});
          }
          if (volumeExtension->hasValue("boundary_material_positive")) {
            cvbConfig.boundaryMaterial[4] = Acts::createProtoMaterial(
                *volumeExtension, "boundary_material_positive",
                {{"binPhi", Acts::closed}, {"binR", Acts::open}});
          }
        }
      } else {
        throw std::logic_error(
            std::string("Current DetElement: ") + volumeDetElement.name() +
            std::string(
                " has wrong ActsExtension! At this stage it should be a "
                "detector volume declared as Barrel or Endcap. Please "
                "check your detector construction."));
      }

      // Fill the volume material for the inner / outer cover
      if (volumeExtension->hasType("boundary_material")) {
        if (volumeExtension->hasValue("boundary_material_inner")) {
          cvbConfig.boundaryMaterial[0] = Acts::createProtoMaterial(
              *volumeExtension, "boundary_material_inner",
              {{"binPhi", Acts::closed}, {"binZ", Acts::open}});
        }
        if (volumeExtension->hasValue("boundary_material_outer")) {
          cvbConfig.boundaryMaterial[1] = Acts::createProtoMaterial(
              *volumeExtension, "boundary_material_outer",
              {{"binPhi", Acts::closed}, {"binZ", Acts::open}});
        }
      }
    }

    if ((pEndCap && !nEndCap) || (!pEndCap && nEndCap)) {
      throw std::logic_error(
          "Only one Endcap is given for the current "
          "hierarchy! Endcaps should always occur in "
          "pairs. Please check your detector "
          "construction.");
    }

    // configure SurfaceArrayCreator
    auto surfaceArrayCreator =
        std::make_shared<const Acts::SurfaceArrayCreator>(
            Acts::getDefaultLogger("D2A_SAC", loggingLevel));
    // configure LayerCreator
    Acts::LayerCreator::Config lcConfig;
    lcConfig.surfaceArrayCreator = surfaceArrayCreator;
    auto layerCreator = std::make_shared<const Acts::LayerCreator>(
        lcConfig, Acts::getDefaultLogger("D2A_LAC", loggingLevel));
    // configure DD4hepLayerBuilder
    Acts::DD4hepLayerBuilder::Config lbConfig;
    lbConfig.configurationName = subDetector.name();
    lbConfig.layerCreator = layerCreator;
    lbConfig.negativeLayers = negativeLayers;
    lbConfig.centralLayers = centralLayers;
    lbConfig.positiveLayers = positiveLayers;
    lbConfig.bTypePhi = bTypePhi;
    lbConfig.bTypeR = bTypeR;
    lbConfig.bTypeZ = bTypeZ;
    lbConfig.defaultThickness = defaultLayerThickness;
    auto dd4hepLayerBuilder = std::make_shared<const Acts::DD4hepLayerBuilder>(
        lbConfig,
        Acts::getDefaultLogger(std::string("D2A_L:") + subDetector.name(),
                               loggingLevel));

    // Create the sub volume
    // Dimensions are created automatically by adding a tolerance to the
    // layer setup
    cvbConfig.layerEnvelopeR = std::make_pair(layerEnvelopeR, layerEnvelopeR);
    cvbConfig.layerEnvelopeZ = layerEnvelopeZ;
    cvbConfig.trackingVolumeHelper = volumeHelper;
    cvbConfig.volumeSignature = 0;
    cvbConfig.volumeName = subDetector.name();
    cvbConfig.layerBuilder = dd4hepLayerBuilder;
    auto cylinderVolumeBuilder =
        std::make_shared<const Acts::CylinderVolumeBuilder>(
            cvbConfig,
            Acts::getDefaultLogger(std::string("D2A_V:") + subDetector.name(),
                                   loggingLevel));
    return cylinderVolumeBuilder;

  } else if ((subDetExtension != nullptr) &&
             (subDetExtension->hasType("passive cylinder", "layer") ||
              subDetExtension->hasType("beampipe", "layer"))) {
    ACTS_VERBOSE("Subdetector : " << subDetector.name()
                                  << " - building a passive cylinder.");
    if (subDetExtension->hasType("beampipe", "layer")) {
      ACTS_VERBOSE("This is the beam pipe - will be built to r -> 0.");
    }

    // get the dimensions of the volume
    TGeoShape* geoShape =
        subDetector.placement().ptr()->GetVolume()->GetShape();
    TGeoTubeSeg* tube = dynamic_cast<TGeoTubeSeg*>(geoShape);
    if (tube == nullptr) {
      throw std::logic_error(
          "Cylinder has wrong shape - needs to be TGeoTubeSeg!");
    }
    // get the dimension of TGeo and convert lengths
    double rMin = tube->GetRmin() * UnitConstants::cm - layerEnvelopeR;
    double rMax = tube->GetRmax() * UnitConstants::cm + layerEnvelopeR;
    double halfZ = tube->GetDz() * UnitConstants::cm + layerEnvelopeZ;
    ACTS_VERBOSE(
        "Extracting cylindrical volume bounds ( rmin / rmax / "
        "halfZ )=  ( "
        << rMin << " / " << rMax << " / " << halfZ << " )");

    std::shared_ptr<Acts::ISurfaceMaterial> plMaterial = nullptr;
    if (subDetExtension->hasType("layer_material")) {
      // get the possible material of the surounding volume
      plMaterial = Acts::createProtoMaterial(
          *subDetExtension, "layer_material_representing",
          {{"binPhi", Acts::closed}, {"binZ", Acts::open}});
    }

    // configure the passive layer builder
    Acts::PassiveLayerBuilder::Config plbConfig;
    plbConfig.layerIdentification = subDetector.name();
    plbConfig.centralLayerRadii = std::vector<double>(1, 0.5 * (rMax + rMin));
    plbConfig.centralLayerHalflengthZ = std::vector<double>(1, halfZ);
    plbConfig.centralLayerThickness = std::vector<double>(1, fabs(rMax - rMin));
    plbConfig.centralLayerMaterial = {plMaterial};
    auto pcLayerBuilder = std::make_shared<const Acts::PassiveLayerBuilder>(
        plbConfig,
        Acts::getDefaultLogger(std::string("D2A_PL:") + subDetector.name(),
                               loggingLevel));

    // the configuration object of the volume builder
    Acts::CylinderVolumeBuilder::Config cvbConfig;
    cvbConfig.trackingVolumeHelper = volumeHelper;
    cvbConfig.volumeSignature = 0;
    cvbConfig.volumeName = subDetector.name();
    cvbConfig.layerBuilder = pcLayerBuilder;
    cvbConfig.layerEnvelopeR = {layerEnvelopeR, layerEnvelopeR};
    cvbConfig.layerEnvelopeZ = layerEnvelopeZ;
    cvbConfig.buildToRadiusZero = subDetExtension->hasType("beampipe", "layer");

    // beam pipe / passive cylinder volume builder
    auto pcVolumeBuilder = std::make_shared<const Acts::CylinderVolumeBuilder>(
        cvbConfig,
        Acts::getDefaultLogger(std::string("D2A_V:") + subDetector.name(),
                               loggingLevel));
    return pcVolumeBuilder;

  } else if ((subDetExtension != nullptr) &&
             subDetExtension->hasType("barrel", "detector")) {
    ACTS_VERBOSE("Subdetector: "
                 << subDetector.name()
                 << " is a (sensitive) Barrel volume - building barrel.");
    /// the dd4hep::DetElements of the layers of the central volume
    std::vector<dd4hep::DetElement> centralLayers, centralVolumes;
    collectLayers_dd4hep(subDetector, centralLayers);
    collectVolumes_dd4hep(subDetector, centralVolumes);

    // configure SurfaceArrayCreator
    auto surfaceArrayCreator =
        std::make_shared<const Acts::SurfaceArrayCreator>(
            Acts::getDefaultLogger("D2A_SAC", loggingLevel));
    // configure LayerCreator
    Acts::LayerCreator::Config lcConfig;
    lcConfig.surfaceArrayCreator = surfaceArrayCreator;
    auto layerCreator = std::make_shared<const Acts::LayerCreator>(
        lcConfig, Acts::getDefaultLogger("D2A_LAC", loggingLevel));
    // configure DD4hepLayerBuilder
    Acts::DD4hepLayerBuilder::Config lbConfig;
    lbConfig.configurationName = subDetector.name();
    lbConfig.layerCreator = layerCreator;
    lbConfig.centralLayers = centralLayers;
    lbConfig.bTypePhi = bTypePhi;
    lbConfig.bTypeZ = bTypeZ;
    lbConfig.defaultThickness = defaultLayerThickness;
    auto dd4hepLayerBuilder = std::make_shared<const Acts::DD4hepLayerBuilder>(
        lbConfig,
        Acts::getDefaultLogger(std::string("D2A_LB_") + subDetector.name(),
                               loggingLevel));

    // Configure DD4hepVolumeBuilder
    Acts::DD4hepVolumeBuilder::Config vbConfig;
    vbConfig.configurationName = subDetector.name();
    vbConfig.centralVolumes = centralVolumes;
    auto dd4hepVolumeBuilder =
        std::make_shared<const Acts::DD4hepVolumeBuilder>(
            vbConfig,
            Acts::getDefaultLogger(std::string("D2A_VB_") + subDetector.name(),
                                   loggingLevel));

    // the configuration object of the volume builder
    Acts::CylinderVolumeBuilder::Config cvbConfig;
    // get the dimensions of the volume
    TGeoShape* geoShape =
        subDetector.placement().ptr()->GetVolume()->GetShape();
    // this should not happen
    if (geoShape == nullptr) {
      throw std::logic_error(std::string("Volume of DetElement: ") +
                             subDetector.name() +
                             std::string(" has no a shape!"));
    }

    cvbConfig.layerEnvelopeR = std::make_pair(layerEnvelopeR, layerEnvelopeR);
    cvbConfig.layerEnvelopeZ = layerEnvelopeZ;
    cvbConfig.trackingVolumeHelper = volumeHelper;
    cvbConfig.volumeSignature = 0;
    cvbConfig.volumeName = subDetector.name();
    cvbConfig.layerBuilder = dd4hepLayerBuilder;
    cvbConfig.ctVolumeBuilder = dd4hepVolumeBuilder;
    auto cylinderVolumeBuilder =
        std::make_shared<const Acts::CylinderVolumeBuilder>(
            cvbConfig,
            Acts::getDefaultLogger(std::string("D2A_V:") + subDetector.name(),
                                   loggingLevel));
    return cylinderVolumeBuilder;

  } else {
    ACTS_INFO(
        "Subdetector with name : '"
        << subDetector.name()
        << "' has wrong ActsExtension for translation and is not of type "
           "'compound'. If you want to have this DetElement be translated "
           "into the tracking geometry you need add the right "
           "ActsExtension (at this stage the subvolume needs to be "
           "declared as beampipe or barrel) or if it is a compound "
           "DetElement (containing a barrel-endcap hierarchy), the type "
           "needs to be set to 'compound'.");
    return nullptr;
  }
}

std::shared_ptr<const Acts::CylinderVolumeHelper> cylinderVolumeHelper_dd4hep(
    Logging::Level loggingLevel) {
  // create cylindervolumehelper which can be used by all instances
  // hand over LayerArrayCreator
  Acts::LayerArrayCreator::Config lacConfig;
  auto layerArrayCreator = std::make_shared<const Acts::LayerArrayCreator>(
      lacConfig, Acts::getDefaultLogger(std::string("D2A_LAC"), loggingLevel));
  // tracking volume array creator
  Acts::TrackingVolumeArrayCreator::Config tvacConfig;
  auto trackingVolumeArrayCreator =
      std::make_shared<const Acts::TrackingVolumeArrayCreator>(
          tvacConfig,
          Acts::getDefaultLogger(std::string("D2A_TVAC"), loggingLevel));
  // configure the cylinder volume helper
  Acts::CylinderVolumeHelper::Config cvhConfig;
  cvhConfig.layerArrayCreator = layerArrayCreator;
  cvhConfig.trackingVolumeArrayCreator = trackingVolumeArrayCreator;
  auto cylinderVolumeHelper =
      std::make_shared<const Acts::CylinderVolumeHelper>(
          cvhConfig,
          Acts::getDefaultLogger(std::string("D2A_CVH"), loggingLevel));

  return cylinderVolumeHelper;
}

void collectCompounds_dd4hep(dd4hep::DetElement& detElement,
                             std::vector<dd4hep::DetElement>& compounds) {
  const dd4hep::DetElement::Children& children = detElement.children();
  for (auto& child : children) {
    dd4hep::DetElement childDetElement = child.second;
    Acts::ActsExtension* detExtension = nullptr;
    try {
      detExtension = childDetElement.extension<Acts::ActsExtension>();
    } catch (std::runtime_error& e) {
    }
    if ((detExtension != nullptr) &&
        (detExtension->hasType("barrel", "detector") ||
         detExtension->hasType("endcap", "detector"))) {
      compounds.push_back(childDetElement);
      continue;
    }
    collectCompounds_dd4hep(childDetElement, compounds);
  }
}

void collectSubDetectors_dd4hep(dd4hep::DetElement& detElement,
                                std::vector<dd4hep::DetElement>& subdetectors) {
  const dd4hep::DetElement::Children& children = detElement.children();
  for (auto& child : children) {
    dd4hep::DetElement childDetElement = child.second;
    Acts::ActsExtension* detExtension = nullptr;
    try {
      detExtension = childDetElement.extension<Acts::ActsExtension>();
    } catch (std::runtime_error& e) {
      if (childDetElement.type() == "compound") {
        subdetectors.push_back(childDetElement);
        continue;
      }
    }
    if ((detExtension != nullptr) &&
        (detExtension->hasType("barrel", "detector") ||
         detExtension->hasType("beampipe", "layer"))) {
      subdetectors.push_back(childDetElement);
      continue;
    }
    collectSubDetectors_dd4hep(childDetElement, subdetectors);
  }
}

void collectLayers_dd4hep(dd4hep::DetElement& detElement,
                          std::vector<dd4hep::DetElement>& layers) {
  const dd4hep::DetElement::Children& children = detElement.children();
  for (auto& child : children) {
    dd4hep::DetElement childDetElement = child.second;
    Acts::ActsExtension* detExtension = nullptr;
    try {
      detExtension = childDetElement.extension<Acts::ActsExtension>();
    } catch (std::runtime_error& e) {
    }
    if ((detExtension != nullptr) && detExtension->hasType("layer")) {
      layers.push_back(childDetElement);
      continue;
    }
    collectLayers_dd4hep(childDetElement, layers);
  }
}

void collectVolumes_dd4hep(dd4hep::DetElement& detElement,
                           std::vector<dd4hep::DetElement>& volumes) {
  const dd4hep::DetElement::Children& children = detElement.children();
  for (auto& child : children) {
    dd4hep::DetElement childDetElement = child.second;
    Acts::ActsExtension* detExtension = nullptr;
    try {
      detExtension = childDetElement.extension<Acts::ActsExtension>();
    } catch (std::runtime_error& e) {
    }
    if ((detExtension != nullptr) && detExtension->hasType("volume")) {
      volumes.push_back(childDetElement);
      continue;
    }
    collectVolumes_dd4hep(childDetElement, volumes);
  }
}
}  // End of namespace Acts
