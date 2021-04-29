// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/ContextualDetector/AlignedDetector.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/ContextualDetector/AlignedDetectorElement.hpp"
#include "ActsExamples/ContextualDetector/AlignmentDecorator.hpp"
#include "ActsExamples/Framework/IContextDecorator.hpp"
#include "ActsExamples/GenericDetector/BuildGenericDetector.hpp"
#include "ActsExamples/GenericDetector/GenericDetectorOptions.hpp"

#include <boost/program_options.hpp>

void AlignedDetector::addOptions(
    boost::program_options::options_description& opt) const {
  // Add the generic geometry options
  ActsExamples::Options::addGenericGeometryOptions(opt);
  // specify the rotation setp
  opt.add_options()(
      "align-seed",
      boost::program_options::value<size_t>()->default_value(1324354657),
      "Seed for the decorator random numbers.")(
      "align-iovsize",
      boost::program_options::value<size_t>()->default_value(100),
      "Size of a valid IOV.")(
      "align-flushsize",
      boost::program_options::value<size_t>()->default_value(200),
      "Span until garbage collection is active.")(
      "align-sigma-iplane",
      boost::program_options::value<double>()->default_value(100.),
      "Sigma of the in-plane misalignment in [um]")(
      "align-sigma-oplane",
      boost::program_options::value<double>()->default_value(50.),
      "Sigma of the out-of-plane misalignment in [um]")(
      "align-sigma-irot",
      boost::program_options::value<double>()->default_value(20.),
      "Sigma of the in-plane rotation misalignment in [mrad]")(
      "align-sigma-orot",
      boost::program_options::value<double>()->default_value(0.),
      "Sigma of the out-of-plane rotation misalignment in [mrad]")(
      "align-loglevel",
      boost::program_options::value<size_t>()->default_value(3),
      "Output log level of the alignment decorator.")(
      "align-firstnominal",
      boost::program_options::value<bool>()->default_value(false),
      "Keep the first iov batch nominal.");
}

auto AlignedDetector::finalize(
    const boost::program_options::variables_map& vm,
    std::shared_ptr<const Acts::IMaterialDecorator> mdecorator)
    -> std::pair<TrackingGeometryPtr, ContextDecorators> {
  // --------------------------------------------------------------------------------
  DetectorElement::ContextType nominalContext;

  auto buildLevel = vm["geo-generic-buildlevel"].template as<size_t>();
  // set geometry building logging level
  Acts::Logging::Level surfaceLogLevel =
      Acts::Logging::Level(vm["geo-surface-loglevel"].template as<size_t>());
  Acts::Logging::Level layerLogLevel =
      Acts::Logging::Level(vm["geo-layer-loglevel"].template as<size_t>());
  Acts::Logging::Level volumeLogLevel =
      Acts::Logging::Level(vm["geo-volume-loglevel"].template as<size_t>());

  bool buildProto =
      (vm["mat-input-type"].template as<std::string>() == "proto");

  /// return the generic detector - with aligned context decorator
  TrackingGeometryPtr aTrackingGeometry =
      ActsExamples::Generic::buildDetector<DetectorElement>(
          nominalContext, detectorStore, buildLevel, std::move(mdecorator),
          buildProto, surfaceLogLevel, layerLogLevel, volumeLogLevel);

  Acts::Logging::Level decoratorLogLevel =
      Acts::Logging::Level(vm["align-loglevel"].template as<size_t>());

  // Let's create a reandom number service
  ActsExamples::RandomNumbers::Config randomNumberConfig;
  randomNumberConfig.seed = vm["align-seed"].template as<size_t>();
  auto randomNumberSvc =
      std::make_shared<ActsExamples::RandomNumbers>(randomNumberConfig);

  // Alignment decorator service
  Decorator::Config agcsConfig;
  agcsConfig.detectorStore = detectorStore;
  agcsConfig.iovSize = vm["align-iovsize"].template as<size_t>();
  agcsConfig.flushSize = vm["align-flushsize"].template as<size_t>();

  // The misalingments
  double sigmaIp = vm["align-sigma-iplane"].template as<double>();
  double sigmaOp = vm["align-sigma-oplane"].template as<double>();
  double sigmaIr = vm["align-sigma-irot"].template as<double>();
  double sigmaOr = vm["align-sigma-orot"].template as<double>();
  agcsConfig.gSigmaX = sigmaIp * Acts::UnitConstants::um;
  agcsConfig.gSigmaY = sigmaIp * Acts::UnitConstants::um;
  agcsConfig.gSigmaZ = sigmaOp * Acts::UnitConstants::um;
  agcsConfig.aSigmaX = sigmaOr * 0.001;  // millirad
  agcsConfig.aSigmaY = sigmaOr * 0.001;  // millirad
  agcsConfig.aSigmaZ = sigmaIr * 0.001;  // millirad
  agcsConfig.randomNumberSvc = randomNumberSvc;
  agcsConfig.firstIovNominal = vm["align-firstnominal"].template as<bool>();

  // Now create the alignment decorator
  ContextDecorators aContextDecorators = {std::make_shared<Decorator>(
      agcsConfig,
      Acts::getDefaultLogger("AlignmentDecorator", decoratorLogLevel))};

  // return the pair of geometry and the alignment decorator(s)
  return std::make_pair<TrackingGeometryPtr, ContextDecorators>(
      std::move(aTrackingGeometry), std::move(aContextDecorators));
}
