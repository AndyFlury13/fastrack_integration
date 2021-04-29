// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootTrajectoryParametersWriter.hpp"

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/detail/TransformationBoundToFree.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "ActsExamples/EventData/AverageSimHits.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsExamples/Utilities/Range.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"

#include <ios>
#include <stdexcept>

#include <TFile.h>
#include <TTree.h>

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;
using Acts::VectorHelpers::theta;

ActsExamples::RootTrajectoryParametersWriter::RootTrajectoryParametersWriter(
    const ActsExamples::RootTrajectoryParametersWriter::Config& cfg,
    Acts::Logging::Level lvl)
    : WriterT(cfg.inputTrajectories, "RootTrajectoryParametersWriter", lvl),
      m_cfg(cfg),
      m_outputFile(cfg.rootFile) {
  // trajectories collection name is already checked by base ctor
  if (cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing particles input collection");
  }
  if (cfg.inputMeasurementParticlesMap.empty()) {
    throw std::invalid_argument("Missing hit-particles map input collection");
  }
  if (cfg.outputFilename.empty()) {
    throw std::invalid_argument("Missing output filename");
  }
  if (m_cfg.outputTreename.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  // Setup ROOT I/O
  if (m_outputFile == nullptr) {
    auto path = joinPaths(m_cfg.outputDir, m_cfg.outputFilename);
    m_outputFile = TFile::Open(path.c_str(), m_cfg.fileMode.c_str());
    if (m_outputFile == nullptr) {
      throw std::ios_base::failure("Could not open '" + path);
    }
  }
  m_outputFile->cd();
  m_outputTree =
      new TTree(m_cfg.outputTreename.c_str(), m_cfg.outputTreename.c_str());
  if (m_outputTree == nullptr)
    throw std::bad_alloc();
  else {
    // I/O parameters
    m_outputTree->Branch("event_nr", &m_eventNr);
    m_outputTree->Branch("multiTraj_nr", &m_multiTrajNr);
    m_outputTree->Branch("subTraj_nr", &m_subTrajNr);
    m_outputTree->Branch("t_barcode", &m_t_barcode, "t_barcode/l");
    m_outputTree->Branch("t_charge", &m_t_charge);
    m_outputTree->Branch("t_time", &m_t_time);
    m_outputTree->Branch("t_vx", &m_t_vx);
    m_outputTree->Branch("t_vy", &m_t_vy);
    m_outputTree->Branch("t_vz", &m_t_vz);
    m_outputTree->Branch("t_px", &m_t_px);
    m_outputTree->Branch("t_py", &m_t_py);
    m_outputTree->Branch("t_pz", &m_t_pz);
    m_outputTree->Branch("t_theta", &m_t_theta);
    m_outputTree->Branch("t_phi", &m_t_phi);
    m_outputTree->Branch("t_eta", &m_t_eta);
    m_outputTree->Branch("t_pT", &m_t_pT);

    m_outputTree->Branch("hasFittedParams", &m_hasFittedParams);
    m_outputTree->Branch("eLOC0_fit", &m_eLOC0_fit);
    m_outputTree->Branch("eLOC1_fit", &m_eLOC1_fit);
    m_outputTree->Branch("ePHI_fit", &m_ePHI_fit);
    m_outputTree->Branch("eTHETA_fit", &m_eTHETA_fit);
    m_outputTree->Branch("eQOP_fit", &m_eQOP_fit);
    m_outputTree->Branch("eT_fit", &m_eT_fit);
    m_outputTree->Branch("err_eLOC0_fit", &m_err_eLOC0_fit);
    m_outputTree->Branch("err_eLOC1_fit", &m_err_eLOC1_fit);
    m_outputTree->Branch("err_ePHI_fit", &m_err_ePHI_fit);
    m_outputTree->Branch("err_eTHETA_fit", &m_err_eTHETA_fit);
    m_outputTree->Branch("err_eQOP_fit", &m_err_eQOP_fit);
    m_outputTree->Branch("err_eT_fit", &m_err_eT_fit);
  }
}

ActsExamples::RootTrajectoryParametersWriter::
    ~RootTrajectoryParametersWriter() {
  if (m_outputFile) {
    m_outputFile->Close();
  }
}

ActsExamples::ProcessCode
ActsExamples::RootTrajectoryParametersWriter::endRun() {
  if (m_outputFile) {
    m_outputFile->cd();
    m_outputTree->Write();
    ACTS_INFO("Write parameters of trajectories to tree '"
              << m_cfg.outputTreename << "' in '"
              << joinPaths(m_cfg.outputDir, m_cfg.outputFilename) << "'");
  }
  return ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode ActsExamples::RootTrajectoryParametersWriter::writeT(
    const AlgorithmContext& ctx, const TrajectoriesContainer& trajectories) {
  using HitParticlesMap = IndexMultimap<ActsFatras::Barcode>;

  if (m_outputFile == nullptr)
    return ProcessCode::SUCCESS;

  // Read additional input collections
  const auto& particles =
      ctx.eventStore.get<SimParticleContainer>(m_cfg.inputParticles);
  const auto& hitParticlesMap =
      ctx.eventStore.get<HitParticlesMap>(m_cfg.inputMeasurementParticlesMap);

  // For each particle within a track, how many hits did it contribute
  std::vector<ParticleHitCount> particleHitCounts;

  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  // Get the event number
  m_eventNr = ctx.eventNumber;

  // Loop over the trajectories
  for (size_t itraj = 0; itraj < trajectories.size(); ++itraj) {
    const auto& traj = trajectories[itraj];

    if (traj.empty()) {
      ACTS_WARNING("Empty trajectories object " << itraj);
      continue;
    }

    // The trajectory index
    m_multiTrajNr = itraj;

    // The trajectory entry indices and the multiTrajectory
    const auto& trackTips = traj.tips();

    // Loop over the entry indices for the subtrajectories
    for (unsigned int isubtraj = 0; isubtraj < trackTips.size(); ++isubtraj) {
      // The subtrajectory index
      m_subTrajNr = isubtraj;
      // The entry index for this subtrajectory
      const auto& trackTip = trackTips[isubtraj];

      // Get the majority truth particle to this track
      identifyContributingParticles(hitParticlesMap, traj, trackTip,
                                    particleHitCounts);
      if (not particleHitCounts.empty()) {
        // Get the barcode of the majority truth particle
        m_t_barcode = particleHitCounts.front().particleId.value();
        // Find the truth particle via the barcode
        auto ip = particles.find(m_t_barcode);
        if (ip != particles.end()) {
          const auto& particle = *ip;
          ACTS_DEBUG("Find the truth particle with barcode = " << m_t_barcode);
          // Get the truth particle info at vertex
          const auto p = particle.absoluteMomentum();
          m_t_charge = particle.charge();
          m_t_time = particle.time();
          m_t_vx = particle.position().x();
          m_t_vy = particle.position().y();
          m_t_vz = particle.position().z();
          m_t_px = p * particle.unitDirection().x();
          m_t_py = p * particle.unitDirection().y();
          m_t_pz = p * particle.unitDirection().z();
          m_t_theta = theta(particle.unitDirection());
          m_t_phi = phi(particle.unitDirection());
          m_t_eta = eta(particle.unitDirection());
          m_t_pT = p * perp(particle.unitDirection());
        } else {
          ACTS_WARNING("Truth particle with barcode = " << m_t_barcode
                                                        << " not found!");
        }
      }

      // Get the fitted track parameter
      m_hasFittedParams = false;
      if (traj.hasTrackParameters(trackTip)) {
        m_hasFittedParams = true;
        const auto& boundParam = traj.trackParameters(trackTip);
        const auto& parameter = boundParam.parameters();
        const auto& covariance = *boundParam.covariance();
        m_eLOC0_fit = parameter[Acts::eBoundLoc0];
        m_eLOC1_fit = parameter[Acts::eBoundLoc1];
        m_ePHI_fit = parameter[Acts::eBoundPhi];
        m_eTHETA_fit = parameter[Acts::eBoundTheta];
        m_eQOP_fit = parameter[Acts::eBoundQOverP];
        m_eT_fit = parameter[Acts::eBoundTime];
        m_err_eLOC0_fit = sqrt(covariance(Acts::eBoundLoc0, Acts::eBoundLoc0));
        m_err_eLOC1_fit = sqrt(covariance(Acts::eBoundLoc1, Acts::eBoundLoc1));
        m_err_ePHI_fit = sqrt(covariance(Acts::eBoundPhi, Acts::eBoundPhi));
        m_err_eTHETA_fit =
            sqrt(covariance(Acts::eBoundTheta, Acts::eBoundTheta));
        m_err_eQOP_fit =
            sqrt(covariance(Acts::eBoundQOverP, Acts::eBoundQOverP));
        m_err_eT_fit = sqrt(covariance(Acts::eBoundTime, Acts::eBoundTime));
      }

      // fill the variables for one track to tree
      m_outputTree->Fill();
    }  // all subtrajectories
  }    // all trajectories

  return ProcessCode::SUCCESS;
}
