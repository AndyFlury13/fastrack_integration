// This file is part of the Acts project.
//
// Copyright (C) 2019-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Detector/IBaseDetector.hpp"
#ifdef ACTS_PLUGIN_ONNX
#include "Acts/Plugins/Onnx/MLTrackClassifier.hpp"
#endif
#include "ActsExamples/Digitization/DigitizationOptions.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Geometry/CommonGeometry.hpp"
#include "ActsExamples/Io/Csv/CsvOptionsReader.hpp"
#include "ActsExamples/Io/Csv/CsvParticleReader.hpp"
#include "ActsExamples/Io/Csv/CsvSimHitReader.hpp"
#include "ActsExamples/Io/Performance/CKFPerformanceWriter.hpp"
#include "ActsExamples/Io/Performance/TrackFinderPerformanceWriter.hpp"
#include "ActsExamples/Io/Root/RootTrajectoryParametersWriter.hpp"
#include "ActsExamples/Io/Root/RootTrajectoryStatesWriter.hpp"
#include "ActsExamples/MagneticField/MagneticFieldOptions.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/TrackFinding/SeedingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/SpacePointMaker.hpp"
#include "ActsExamples/TrackFinding/SpacePointMakerOptions.hpp"
#include "ActsExamples/TrackFinding/TrackFindingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/TrackFindingOptions.hpp"
#include "ActsExamples/TrackFinding/TrackParamsEstimationAlgorithm.hpp"
#include "ActsExamples/TruthTracking/ParticleSmearing.hpp"
#include "ActsExamples/TruthTracking/TruthSeedSelector.hpp"
#include "ActsExamples/TruthTracking/TruthTrackFinder.hpp"
#include "ActsExamples/Utilities/Options.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include <Acts/Definitions/Units.hpp>

#include <memory>

#include <boost/filesystem.hpp>

#include "RecInput.hpp"


// Fastrack start

#include "fastrack/lib/ModelClass.h"


// Fastrack end 

using namespace Acts::UnitLiterals;
using namespace ActsExamples;
using namespace boost::filesystem;
using namespace std::placeholders;

void addRecCKFOptions(ActsExamples::Options::Description& desc) {
  using namespace ActsExamples;
  using boost::program_options::bool_switch;

  auto opt = desc.add_options();
  opt("ckf-truth-smeared-seeds", bool_switch(),
      "Use track parameters smeared from truth particles for steering CKF");
  opt("ckf-truth-estimated-seeds", bool_switch(),
      "Use track parameters estimated from truth tracks for steering CKF");
}

int runRecCKFTracks(int argc, char* argv[],
                    std::shared_ptr<ActsExamples::IBaseDetector> detector) {
  // setup and parse options
  auto desc = ActsExamples::Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addRandomNumbersOptions(desc);
  Options::addGeometryOptions(desc);
  Options::addMaterialOptions(desc);
  Options::addInputOptions(desc);
  Options::addOutputOptions(desc, OutputFormat::DirectoryOnly);
  detector->addOptions(desc);
  Options::addMagneticFieldOptions(desc);
  Options::addTrackFindingOptions(desc);
  addRecCKFOptions(desc);
  Options::addDigitizationOptions(desc);
  Options::addSpacePointMakerOptions(desc);

  auto vm = Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  Sequencer sequencer(Options::readSequencerConfig(vm));

  // Read some standard options
  auto logLevel = Options::readLogLevel(vm);
  auto inputDir = vm["input-dir"].as<std::string>();
  auto outputDir = ensureWritableDirectory(vm["output-dir"].as<std::string>());
  auto rnd = std::make_shared<ActsExamples::RandomNumbers>(
      Options::readRandomNumbersConfig(vm));
  bool truthSmearedSeeded = vm["ckf-truth-smeared-seeds"].template as<bool>();
  bool truthEstimatedSeeded =
      vm["ckf-truth-estimated-seeds"].template as<bool>();

  // Setup detector geometry
  auto geometry = Geometry::build(vm, *detector);
  auto trackingGeometry = geometry.first;
  // Add context decorators
  for (auto cdr : geometry.second) {
    sequencer.addContextDecorator(cdr);
  }
  // Setup the magnetic field
  auto magneticField = Options::readMagneticField(vm);

  // Read the sim hits
  auto simHitReaderCfg = setupSimHitReading(vm, sequencer);
  // Read the particles
  auto particleReader = setupParticleReading(vm, sequencer);

  // Run the sim hits smearing
  auto digiCfg = setupDigitization(vm, sequencer, rnd, trackingGeometry,
                                   simHitReaderCfg.outputSimHits);

  // Run the particle selection
  // The pre-selection will select truth particles satisfying provided criteria
  // from all particles read in by particle reader for further processing. It
  // has no impact on the truth hits read-in by the cluster reader.
  TruthSeedSelector::Config particleSelectorCfg;
  particleSelectorCfg.inputParticles = particleReader.outputParticles;
  particleSelectorCfg.inputMeasurementParticlesMap =
      digiCfg.outputMeasurementParticlesMap;
  particleSelectorCfg.outputParticles = "particles_selected";
  particleSelectorCfg.ptMin = 1_GeV;
  particleSelectorCfg.nHitsMin = 9;
  sequencer.addAlgorithm(
      std::make_shared<TruthSeedSelector>(particleSelectorCfg, logLevel));

  // The selected particles
  const auto& inputParticles = particleSelectorCfg.outputParticles;

  // Create starting parameters from either particle smearing or combined seed
  // finding and track parameters estimation
  std::string outputTrackParameters;
  if (truthSmearedSeeded) {
    // Run the particle smearing
    auto particleSmearingCfg =
        setupParticleSmearing(vm, sequencer, rnd, inputParticles);
    outputTrackParameters = particleSmearingCfg.outputTrackParameters;
  } else {
    // Create space points
    SpacePointMaker::Config spCfg = Options::readSpacePointMakerConfig(vm);
    spCfg.inputSourceLinks = digiCfg.outputSourceLinks;
    spCfg.inputMeasurements = digiCfg.outputMeasurements;
    spCfg.outputSpacePoints = "spacepoints";
    spCfg.trackingGeometry = trackingGeometry;
    sequencer.addAlgorithm(std::make_shared<SpacePointMaker>(spCfg, logLevel));

    // Create seeds (i.e. proto tracks) using either truth track finding or seed
    // finding algorithm
    std::string inputProtoTracks = "";
    std::string inputSeeds = "";
    if (truthEstimatedSeeded) {
      // Truth track finding algorithm
      TruthTrackFinder::Config trackFinderCfg;
      trackFinderCfg.inputParticles = inputParticles;
      trackFinderCfg.inputMeasurementParticlesMap =
          digiCfg.outputMeasurementParticlesMap;
      trackFinderCfg.outputProtoTracks = "prototracks";

      //// Fastrack start
      char geom_name[] = "geometry.bin";
      int geom_name_len = geom_name.length();
      char connections_name[] = "connetions.bin";
      int connections_name_len = connections_name.length();
      ModelClass* fastrack_model = new ModelClass("geometry.bin", geom_name_len, connections_name, connections_name_len);

      int nhits;
      int* hit_ids
      int* x;
      int* y;
      int* z;
      const float* px;
      const float* py;
      const float* pz;
      const unsigned long* particle;
      const float* w;
      fastrack_model::importHits(nhits, hit_ids, x, y, z, px, py, pz particle, w);

      int ncells;
      const int* input_hit_id;
      const int* ch0;
      const int* ch1;
      fastrack_modeL::importCells(ncells, input_hit_id, ch0, ch1);

      int* labels;
      fastrack_model::findTracks(labels);
      sequencer.addAlgorithm(std::make_shared<TruthTrackFinder>())
      //// Fastrack end


      sequencer.addAlgorithm(
          std::make_shared<TruthTrackFinder>(trackFinderCfg, logLevel));
      inputProtoTracks = trackFinderCfg.outputProtoTracks;
    } else {
      // Seeding algorithm
      SeedingAlgorithm::Config seedingCfg;
      seedingCfg.inputSpacePoints = {
          spCfg.outputSpacePoints,
      };
      seedingCfg.outputSeeds = "seeds";
      seedingCfg.outputProtoTracks = "prototracks";
      seedingCfg.rMax = 200.;
      seedingCfg.deltaRMax = 60.;
      seedingCfg.collisionRegionMin = -250;
      seedingCfg.collisionRegionMax = 250.;
      seedingCfg.zMin = -2000.;
      seedingCfg.zMax = 2000.;
      seedingCfg.maxSeedsPerSpM = 1;
      seedingCfg.cotThetaMax = 7.40627;  // 2.7 eta
      seedingCfg.sigmaScattering = 50;
      seedingCfg.radLengthPerSeed = 0.1;
      seedingCfg.minPt = 500.;
      seedingCfg.bFieldInZ = 0.00199724;
      seedingCfg.beamPosX = 0;
      seedingCfg.beamPosY = 0;
      seedingCfg.impactMax = 3.;
      sequencer.addAlgorithm(
          std::make_shared<SeedingAlgorithm>(seedingCfg, logLevel));
      inputProtoTracks = seedingCfg.outputProtoTracks;
      inputSeeds = seedingCfg.outputSeeds;
    }

    // write track finding/seeding performance
    TrackFinderPerformanceWriter::Config tfPerfCfg;
    tfPerfCfg.inputProtoTracks = inputProtoTracks;
    // using selected particles
    tfPerfCfg.inputParticles = inputParticles;
    tfPerfCfg.inputMeasurementParticlesMap =
        digiCfg.outputMeasurementParticlesMap;
    tfPerfCfg.outputDir = outputDir;
    tfPerfCfg.outputFilename = "performance_seeding_trees.root";
    sequencer.addWriter(
        std::make_shared<TrackFinderPerformanceWriter>(tfPerfCfg, logLevel));

    // Algorithm estimating track parameter from seed
    TrackParamsEstimationAlgorithm::Config paramsEstimationCfg;
    paramsEstimationCfg.inputSeeds = inputSeeds;
    paramsEstimationCfg.inputProtoTracks = inputProtoTracks;
    paramsEstimationCfg.inputSpacePoints = {
        spCfg.outputSpacePoints,
    };
    paramsEstimationCfg.inputSourceLinks = digiCfg.outputSourceLinks;
    paramsEstimationCfg.outputTrackParameters = "estimatedparameters";
    paramsEstimationCfg.outputProtoTracks = "prototracks_estimated";
    paramsEstimationCfg.trackingGeometry = trackingGeometry;
    paramsEstimationCfg.magneticField = magneticField;
    paramsEstimationCfg.bFieldMin = 0.1_T;
    paramsEstimationCfg.deltaRMax = 100._mm;
    paramsEstimationCfg.sigmaLoc0 = 25._um;
    paramsEstimationCfg.sigmaLoc1 = 100._um;
    paramsEstimationCfg.sigmaPhi = 0.005_degree;
    paramsEstimationCfg.sigmaTheta = 0.001_degree;
    paramsEstimationCfg.sigmaQOverP = 0.1 / 1._GeV;
    paramsEstimationCfg.sigmaT0 = 1400._s;
    sequencer.addAlgorithm(std::make_shared<TrackParamsEstimationAlgorithm>(
        paramsEstimationCfg, logLevel));

    outputTrackParameters = paramsEstimationCfg.outputTrackParameters;
  }

  // Setup the track finding algorithm with CKF
  // It takes all the source links created from truth hit smearing, seeds from
  // truth particle smearing and source link selection config
  auto trackFindingCfg = Options::readTrackFindingConfig(vm);
  trackFindingCfg.inputMeasurements = digiCfg.outputMeasurements;
  trackFindingCfg.inputSourceLinks = digiCfg.outputSourceLinks;
  trackFindingCfg.inputInitialTrackParameters = outputTrackParameters;
  trackFindingCfg.outputTrajectories = "trajectories";
  trackFindingCfg.findTracks = TrackFindingAlgorithm::makeTrackFinderFunction(
      trackingGeometry, magneticField);
  sequencer.addAlgorithm(
      std::make_shared<TrackFindingAlgorithm>(trackFindingCfg, logLevel));

  // write track states from CKF
  RootTrajectoryStatesWriter::Config trackStatesWriter;
  trackStatesWriter.inputTrajectories = trackFindingCfg.outputTrajectories;
  // @note The full particles collection is used here to avoid lots of warnings
  // since the unselected CKF track might have a majority particle not in the
  // filtered particle collection. Thsi could be avoided when a seperate track
  // selection algorithm is used.
  trackStatesWriter.inputParticles = particleReader.outputParticles;
  trackStatesWriter.inputSimHits = simHitReaderCfg.outputSimHits;
  trackStatesWriter.inputMeasurementParticlesMap =
      digiCfg.outputMeasurementParticlesMap;
  trackStatesWriter.inputMeasurementSimHitsMap =
      digiCfg.outputMeasurementSimHitsMap;
  trackStatesWriter.outputDir = outputDir;
  trackStatesWriter.outputFilename = "trackstates_ckf.root";
  trackStatesWriter.outputTreename = "trackstates_ckf";
  sequencer.addWriter(std::make_shared<RootTrajectoryStatesWriter>(
      trackStatesWriter, logLevel));

  // write track parameters from CKF
  RootTrajectoryParametersWriter::Config trackParamsWriter;
  trackParamsWriter.inputTrajectories = trackFindingCfg.outputTrajectories;
  // @note The full particles collection is used here to avoid lots of warnings
  // since the unselected CKF track might have a majority particle not in the
  // filtered particle collection. Thsi could be avoided when a seperate track
  // selection algorithm is used.
  trackParamsWriter.inputParticles = particleReader.outputParticles;
  trackParamsWriter.inputMeasurementParticlesMap =
      digiCfg.outputMeasurementParticlesMap;
  trackParamsWriter.outputDir = outputDir;
  trackParamsWriter.outputFilename = "trackparams_ckf.root";
  trackParamsWriter.outputTreename = "trackparams_ckf";
  sequencer.addWriter(std::make_shared<RootTrajectoryParametersWriter>(
      trackParamsWriter, logLevel));

  // Write CKF performance data
  CKFPerformanceWriter::Config perfWriterCfg;
  perfWriterCfg.inputParticles = inputParticles;
  perfWriterCfg.inputTrajectories = trackFindingCfg.outputTrajectories;
  perfWriterCfg.inputMeasurementParticlesMap =
      digiCfg.outputMeasurementParticlesMap;
  // The bottom seed on a pixel detector 'eats' one or two measurements?
  perfWriterCfg.nMeasurementsMin = particleSelectorCfg.nHitsMin - 2;
  perfWriterCfg.outputDir = outputDir;
#ifdef ACTS_PLUGIN_ONNX
  // Onnx plugin related options
  // Path to default demo ML model for track classification
  path currentFilePath(__FILE__);
  path parentPath = currentFilePath.parent_path();
  path demoModelPath =
      canonical(parentPath / "MLAmbiguityResolutionDemo.onnx").native();
  // Threshold probability for neural network to classify track as duplicate
  double decisionThreshProb = 0.5;
  // Initialize OnnxRuntime plugin
  Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "MLTrackClassifier");
  Acts::MLTrackClassifier neuralNetworkClassifier(env, demoModelPath.c_str());
  perfWriterCfg.duplicatedPredictor =
      std::bind(&Acts::MLTrackClassifier::isDuplicate, &neuralNetworkClassifier,
                std::placeholders::_1, decisionThreshProb);
#endif
  sequencer.addWriter(
      std::make_shared<CKFPerformanceWriter>(perfWriterCfg, logLevel));

  return sequencer.run();
}
