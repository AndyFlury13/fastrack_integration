// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/DD4hepDetector/DD4hepGeometryService.hpp"

#include "Acts/Geometry/CylinderVolumeBuilder.hpp"
#include "Acts/Geometry/CylinderVolumeHelper.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolumeArrayCreator.hpp"
#include "Acts/Plugins/DD4hep/ConvertDD4hepDetector.hpp"
#include "Acts/Utilities/BinningType.hpp"

#include <stdexcept>

#include "DD4hep/Printout.h"

ActsExamples::DD4hep::DD4hepGeometryService::DD4hepGeometryService(
    const ActsExamples::DD4hep::DD4hepGeometryService::Config& cfg)
    : BareService("DD4hepGeometryService", cfg.logLevel), m_cfg(cfg) {
  if (m_cfg.xmlFileNames.empty()) {
    throw std::invalid_argument("Missing DD4hep XML filenames");
  }
}

ActsExamples::DD4hep::DD4hepGeometryService::~DD4hepGeometryService() {
  if (m_lcdd)
    m_lcdd->destroyInstance();
}

ActsExamples::ProcessCode
ActsExamples::DD4hep::DD4hepGeometryService::buildDD4hepGeometry() {
  switch (m_cfg.logLevel) {
    case Acts::Logging::Level::VERBOSE:
      dd4hep::setPrintLevel(dd4hep::PrintLevel::VERBOSE);
      break;
    case Acts::Logging::Level::DEBUG:
      dd4hep::setPrintLevel(dd4hep::PrintLevel::DEBUG);
      break;
    case Acts::Logging::Level::INFO:
      dd4hep::setPrintLevel(dd4hep::PrintLevel::INFO);
      break;
    case Acts::Logging::Level::WARNING:
      dd4hep::setPrintLevel(dd4hep::PrintLevel::WARNING);
      break;
    case Acts::Logging::Level::ERROR:
      dd4hep::setPrintLevel(dd4hep::PrintLevel::ERROR);
      break;
    case Acts::Logging::Level::FATAL:
      dd4hep::setPrintLevel(dd4hep::PrintLevel::FATAL);
      break;
    case Acts::Logging::Level::MAX:
      dd4hep::setPrintLevel(dd4hep::PrintLevel::ALWAYS);
      break;
  }
  m_lcdd = &(dd4hep::Detector::getInstance());
  for (auto& file : m_cfg.xmlFileNames) {
    m_lcdd->fromCompact(file.c_str());
  }
  m_lcdd->volumeManager();
  m_lcdd->apply("DD4hepVolumeManager", 0, 0);
  m_dd4hepGeometry = m_lcdd->world();

  return ActsExamples::ProcessCode::SUCCESS;
}

dd4hep::DetElement
ActsExamples::DD4hep::DD4hepGeometryService::dd4hepGeometry() {
  if (!m_dd4hepGeometry)
    buildDD4hepGeometry();
  return m_dd4hepGeometry;
}

dd4hep::Detector*
ActsExamples::DD4hep::DD4hepGeometryService::DD4hepGeometryService::lcdd() {
  if (!m_lcdd)
    buildDD4hepGeometry();
  return m_lcdd;
}

TGeoNode* ActsExamples::DD4hep::DD4hepGeometryService::tgeoGeometry() {
  if (!m_dd4hepGeometry)
    buildDD4hepGeometry();
  return m_dd4hepGeometry.placement().ptr();
}

ActsExamples::ProcessCode
ActsExamples::DD4hep::DD4hepGeometryService::buildTrackingGeometry(
    const Acts::GeometryContext& gctx) {
  // Set the tracking geometry
  m_trackingGeometry = Acts::convertDD4hepDetector(
      dd4hepGeometry(), m_cfg.logLevel, m_cfg.bTypePhi, m_cfg.bTypeR,
      m_cfg.bTypeZ, m_cfg.envelopeR, m_cfg.envelopeZ,
      m_cfg.defaultLayerThickness, m_cfg.sortDetectors, gctx,
      m_cfg.matDecorator);
  return ActsExamples::ProcessCode::SUCCESS;
}

std::unique_ptr<const Acts::TrackingGeometry>
ActsExamples::DD4hep::DD4hepGeometryService::trackingGeometry(
    const Acts::GeometryContext& gctx) {
  if (!m_trackingGeometry) {
    buildTrackingGeometry(gctx);
  }
  return std::move(m_trackingGeometry);
}
