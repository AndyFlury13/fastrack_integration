// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/MagneticField/MagneticFieldOptions.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Options/Pythia8Options.hpp"
#include "ActsExamples/TruthTracking/ParticleSelector.hpp"
#include "ActsExamples/TruthTracking/ParticleSmearing.hpp"
#include "ActsExamples/Vertexing/IterativeVertexFinderAlgorithm.hpp"
#include "ActsExamples/Vertexing/VertexingOptions.hpp"

#include <memory>

using namespace Acts::UnitLiterals;
using namespace ActsExamples;

int main(int argc, char* argv[]) {
  // setup and parse options
  auto desc = Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addRandomNumbersOptions(desc);
  Options::addPythia8Options(desc);
  ParticleSelector::addOptions(desc);
  Options::addVertexingOptions(desc);
  Options::addMagneticFieldOptions(desc);
  Options::addOutputOptions(desc, OutputFormat::DirectoryOnly);
  auto vars = Options::parse(desc, argc, argv);
  if (vars.empty()) {
    return EXIT_FAILURE;
  }

  // basic setup
  auto logLevel = Options::readLogLevel(vars);
  auto rnd =
      std::make_shared<RandomNumbers>(Options::readRandomNumbersConfig(vars));
  Sequencer sequencer(Options::readSequencerConfig(vars));

  // Setup the magnetic field
  auto magneticField = Options::readMagneticField(vars);

  // setup event generator
  EventGenerator::Config evgen = Options::readPythia8Options(vars, logLevel);
  evgen.outputParticles = "particles_generated";
  evgen.randomNumbers = rnd;
  sequencer.addReader(std::make_shared<EventGenerator>(evgen, logLevel));

  // pre-select particles
  ParticleSelector::Config selectParticles = ParticleSelector::readConfig(vars);
  selectParticles.inputParticles = evgen.outputParticles;
  selectParticles.outputParticles = "particles_selected";
  // smearing only works with charge particles for now
  selectParticles.removeNeutral = true;
  selectParticles.absEtaMax = vars["vertexing-eta-max"].as<double>();
  selectParticles.rhoMax = vars["vertexing-rho-max"].as<double>() * 1_mm;
  selectParticles.ptMin = vars["vertexing-pt-min"].as<double>() * 1_MeV;
  sequencer.addAlgorithm(
      std::make_shared<ParticleSelector>(selectParticles, logLevel));

  // simulate track reconstruction by smearing truth track parameters
  ParticleSmearing::Config smearParticles;
  smearParticles.inputParticles = selectParticles.outputParticles;
  smearParticles.outputTrackParameters = "trackparameters";
  smearParticles.randomNumbers = rnd;
  sequencer.addAlgorithm(
      std::make_shared<ParticleSmearing>(smearParticles, logLevel));

  // find vertices
  IterativeVertexFinderAlgorithm::Config findVertices(magneticField);
  findVertices.inputTrackParameters = smearParticles.outputTrackParameters;
  findVertices.outputProtoVertices = "protovertices";
  sequencer.addAlgorithm(
      std::make_shared<IterativeVertexFinderAlgorithm>(findVertices, logLevel));

  return sequencer.run();
}
