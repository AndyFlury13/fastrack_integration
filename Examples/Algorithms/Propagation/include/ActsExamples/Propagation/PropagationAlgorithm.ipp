// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <Acts/Utilities/Helpers.hpp>

#include <random>

template <typename propagator_t>
std::optional<Acts::BoundSymMatrix>
PropagationAlgorithm<propagator_t>::generateCovariance(
    ActsExamples::RandomEngine& rnd,
    std::normal_distribution<double>& gauss) const {
  if (m_cfg.covarianceTransport) {
    // We start from the correlation matrix
    Acts::BoundSymMatrix newCov(m_cfg.correlations);
    // Then we draw errors according to the error values
    Acts::BoundVector covs_smeared = m_cfg.covariances;
    for (size_t k = 0; k < size_t(covs_smeared.size()); ++k) {
      covs_smeared[k] *= gauss(rnd);
    }
    // and apply a double loop
    for (size_t i = 0; i < size_t(newCov.rows()); ++i) {
      for (size_t j = 0; j < size_t(newCov.cols()); ++j) {
        (newCov)(i, j) *= covs_smeared[i];
        (newCov)(i, j) *= covs_smeared[j];
      }
    }
    return newCov;
  }
  return std::nullopt;
}

template <typename propagator_t>
PropagationAlgorithm<propagator_t>::PropagationAlgorithm(
    const PropagationAlgorithm<propagator_t>::Config& cfg,
    Acts::Logging::Level loglevel)
    : BareAlgorithm("PropagationAlgorithm", loglevel), m_cfg(cfg) {}

/// Templated execute test method for
/// charged and netural particles
/// @param [in] context is the contextual data of this event
/// @param [in] startParameters the start parameters
/// @param [in] pathLength the maximal path length to go
template <typename propagator_t>
template <typename parameters_t>
PropagationOutput PropagationAlgorithm<propagator_t>::executeTest(
    const AlgorithmContext& context, const parameters_t& startParameters,
    double pathLength) const {
  ACTS_DEBUG("Test propagation/extrapolation starts");

  PropagationOutput pOutput;

  // This is the outside in mode
  if (m_cfg.mode == 0) {
    // The step length logger for testing & end of world aborter
    using MaterialInteractor = Acts::MaterialInteractor;
    using SteppingLogger = Acts::detail::SteppingLogger;
    using EndOfWorld = Acts::EndOfWorldReached;

    // Action list and abort list
    using ActionList = Acts::ActionList<SteppingLogger, MaterialInteractor>;
    using AbortList = Acts::AbortList<EndOfWorld>;
    using PropagatorOptions =
        Acts::DenseStepperPropagatorOptions<ActionList, AbortList>;

    PropagatorOptions options(context.geoContext, context.magFieldContext,
                              Acts::LoggerWrapper{logger()});
    options.pathLimit = pathLength;

    // Activate loop protection at some pt value
    options.loopProtection =
        (startParameters.transverseMomentum() < m_cfg.ptLoopers);

    // Switch the material interaction on/off & eventually into logging mode
    auto& mInteractor = options.actionList.get<MaterialInteractor>();
    mInteractor.multipleScattering = m_cfg.multipleScattering;
    mInteractor.energyLoss = m_cfg.energyLoss;
    mInteractor.recordInteractions = m_cfg.recordMaterialInteractions;

    // Switch the logger to sterile, e.g. for timing checks
    auto& sLogger = options.actionList.get<SteppingLogger>();
    sLogger.sterile = m_cfg.sterileLogger;
    // Set a maximum step size
    options.maxStepSize = m_cfg.maxStepSize;

    // Propagate using the propagator
    auto result = m_cfg.propagator.propagate(startParameters, options);
    if (result.ok()) {
      const auto& resultValue = result.value();
      auto steppingResults =
          resultValue.template get<SteppingLogger::result_type>();

      // Set the stepping result
      pOutput.first = std::move(steppingResults.steps);
      // Also set the material recording result - if configured
      if (m_cfg.recordMaterialInteractions) {
        auto materialResult =
            resultValue.template get<MaterialInteractor::result_type>();
        pOutput.second = std::move(materialResult);
      }
    }
  }
  return pOutput;
}

template <typename propagator_t>
ProcessCode PropagationAlgorithm<propagator_t>::execute(
    const AlgorithmContext& context) const {
  // Create a random number generator
  ActsExamples::RandomEngine rng =
      m_cfg.randomNumberSvc->spawnGenerator(context);

  // Standard gaussian distribution for covarianmces
  std::normal_distribution<double> gauss(0., 1.);

  // Setup random number distributions for some quantities
  std::uniform_real_distribution<double> phiDist(m_cfg.phiRange.first,
                                                 m_cfg.phiRange.second);
  std::uniform_real_distribution<double> etaDist(m_cfg.etaRange.first,
                                                 m_cfg.etaRange.second);
  std::uniform_real_distribution<double> ptDist(m_cfg.ptRange.first,
                                                m_cfg.ptRange.second);
  std::uniform_real_distribution<double> qDist(0., 1.);

  std::shared_ptr<const Acts::PerigeeSurface> surface =
      Acts::Surface::makeShared<Acts::PerigeeSurface>(
          Acts::Vector3(0., 0., 0.));

  // Output : the propagation steps
  std::vector<std::vector<Acts::detail::Step>> propagationSteps;
  propagationSteps.reserve(m_cfg.ntests);

  // Output (optional): the recorded material
  std::vector<RecordedMaterialTrack> recordedMaterial;
  if (m_cfg.recordMaterialInteractions) {
    recordedMaterial.reserve(m_cfg.ntests);
  }

  // loop over number of particles
  for (size_t it = 0; it < m_cfg.ntests; ++it) {
    /// get the d0 and z0
    double d0 = m_cfg.d0Sigma * gauss(rng);
    double z0 = m_cfg.z0Sigma * gauss(rng);
    double phi = phiDist(rng);
    double eta = etaDist(rng);
    double theta = 2 * atan(exp(-eta));
    double pt = ptDist(rng);
    double p = pt / sin(theta);
    double charge = qDist(rng) > 0.5 ? 1. : -1.;
    double qop = charge / p;
    double t = m_cfg.tSigma * gauss(rng);
    // parameters
    Acts::BoundVector pars;
    pars << d0, z0, phi, theta, qop, t;
    // some screen output

    Acts::Vector3 sPosition(0., 0., 0.);
    Acts::Vector3 sMomentum(0., 0., 0.);

    // The covariance generation
    auto cov = generateCovariance(rng, gauss);

    // execute the test for charged particles
    PropagationOutput pOutput;
    if (charge) {
      // charged extrapolation - with hit recording
      Acts::BoundTrackParameters startParameters(surface, std::move(pars),
                                                 std::move(cov));
      sPosition = startParameters.position(context.geoContext);
      sMomentum = startParameters.momentum();
      pOutput = executeTest(context, startParameters);
    } else {
      // execute the test for neutral particles
      Acts::NeutralBoundTrackParameters neutralParameters(
          surface, std::move(pars), std::move(cov));
      sPosition = neutralParameters.position(context.geoContext);
      sMomentum = neutralParameters.momentum();
      pOutput = executeTest(context, neutralParameters);
    }
    // Record the propagator steps
    propagationSteps.push_back(std::move(pOutput.first));
    if (m_cfg.recordMaterialInteractions &&
        pOutput.second.materialInteractions.size()) {
      // Create a recorded material track
      RecordedMaterialTrack rmTrack;
      // Start position
      rmTrack.first.first = std::move(sPosition);
      // Start momentum
      rmTrack.first.second = std::move(sMomentum);
      // The material
      rmTrack.second = std::move(pOutput.second);
      // push it it
      recordedMaterial.push_back(std::move(rmTrack));
    }
  }

  // Write the propagation step data to the event store
  context.eventStore.add(m_cfg.propagationStepCollection,
                         std::move(propagationSteps));

  // Write the recorded material to the event store
  if (m_cfg.recordMaterialInteractions) {
    context.eventStore.add(m_cfg.propagationMaterialCollection,
                           std::move(recordedMaterial));
  }

  return ProcessCode::SUCCESS;
}
