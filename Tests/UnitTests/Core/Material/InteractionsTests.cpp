// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/Material/Interactions.hpp"
#include "Acts/Tests/CommonHelpers/PredefinedMaterials.hpp"
#include "Acts/Utilities/PdgParticle.hpp"

namespace data = boost::unit_test::data;
using namespace Acts::UnitLiterals;

// fixed material
static const Acts::Material material = Acts::Test::makeSilicon();
// variable values for other parameters
// thickness
static const double valuesThickness[] = {200_um, 1_mm};
static auto thickness = data::make(valuesThickness);
// particle type, mass, and charge
static const int pdg[] = {Acts::eElectron, Acts::eMuon, Acts::ePionPlus,
                          Acts::eProton};
static const double mass[] = {511_keV, 105.7_MeV, 139.6_MeV, 938.3_MeV};
static const double charge[] = {-1_e, -1_e, 1_e, 1_e};
static const auto particle =
    data::make(pdg) ^ data::make(mass) ^ data::make(charge);
// momentum range
static const auto momentum_low = data::xrange(100_MeV, 10_GeV, 100_MeV);
static const auto momentum_med = data::xrange(10_GeV, 100_GeV, 10_GeV);
static const auto momentum_high = data::xrange(100_GeV, 10_TeV, 100_GeV);
static const auto momentum = momentum_low + momentum_med + momentum_high;

BOOST_AUTO_TEST_SUITE(interactions)

// consistency checks for the energy loss values
BOOST_DATA_TEST_CASE(energy_loss_consistency, thickness* particle* momentum, x,
                     i, m, q, p) {
  const auto slab = Acts::MaterialSlab(material, x);
  const auto qOverP = q / p;

  auto dEBethe = computeEnergyLossBethe(slab, i, m, qOverP, q);
  auto dELandau = computeEnergyLossLandau(slab, i, m, qOverP, q);
  auto dELandauSigma = computeEnergyLossLandauSigma(slab, i, m, qOverP, q);
  auto dELandauSigmaQOverP =
      computeEnergyLossLandauSigmaQOverP(slab, i, m, qOverP, q);
  auto dERad = computeEnergyLossRadiative(slab, i, m, qOverP, q);
  auto dEMean = computeEnergyLossMean(slab, i, m, qOverP, q);
  auto dEMode = computeEnergyLossMode(slab, i, m, qOverP, q);

  BOOST_CHECK_LT(0, dEBethe);
  BOOST_CHECK_LT(0, dELandau);
  BOOST_CHECK_LT(0, dELandauSigma);
  BOOST_CHECK_LT(0, dELandauSigmaQOverP);
  BOOST_CHECK_LE(dELandauSigma, dEBethe);
  // radiative terms only kick above some threshold -> can be zero
  BOOST_CHECK_LE(0, dERad);
  BOOST_CHECK_LT(0, dEMean);
  BOOST_CHECK_LT(0, dEMode);
  BOOST_CHECK_LE((dEBethe + dERad), dEMean);
  // TODO verify mode/mean relation for full energy loss
  // BOOST_CHECK_LE(dEMode, dEMean);
}

// consistency checks for multiple scattering
BOOST_DATA_TEST_CASE(multiple_scattering_consistency,
                     thickness* particle* momentum, x, i, m, q, p) {
  const auto slab = Acts::MaterialSlab(material, x);
  const auto slabDoubled = Acts::MaterialSlab(material, 2 * x);
  const auto qOverP = q / p;
  const auto qOver2P = q / (2 * p);

  auto t0 = computeMultipleScatteringTheta0(slab, i, m, qOverP, q);
  BOOST_CHECK_LT(0, t0);
  // use the anti-particle -> same scattering
  auto tanti = computeMultipleScatteringTheta0(slab, -i, m, -qOverP, -q);
  BOOST_CHECK_LT(0, tanti);
  BOOST_CHECK_EQUAL(t0, tanti);
  // double the material -> more scattering
  auto t2x = computeMultipleScatteringTheta0(slabDoubled, i, m, qOverP, q);
  BOOST_CHECK_LT(0, t2x);
  BOOST_CHECK_LT(t0, t2x);
  // double the momentum -> less scattering
  auto t2p = computeMultipleScatteringTheta0(slab, i, m, qOver2P, q);
  BOOST_CHECK_LT(0, t2p);
  BOOST_CHECK_LT(t2p, t0);
}

// no material -> no interactions
BOOST_DATA_TEST_CASE(vacuum, thickness* particle* momentum, x, i, m, q, p) {
  const auto vacuum = Acts::MaterialSlab(Acts::Material(), x);
  const auto qOverP = q / p;

  BOOST_CHECK_EQUAL(computeEnergyLossBethe(vacuum, i, m, qOverP, q), 0);
  BOOST_CHECK_EQUAL(computeEnergyLossLandau(vacuum, i, m, qOverP, q), 0);
  BOOST_CHECK_EQUAL(computeEnergyLossLandauSigma(vacuum, i, m, qOverP, q), 0);
  BOOST_CHECK_EQUAL(computeEnergyLossLandauSigmaQOverP(vacuum, i, m, qOverP, q),
                    0);
  BOOST_CHECK_EQUAL(computeEnergyLossRadiative(vacuum, i, m, qOverP, q), 0);
  BOOST_CHECK_EQUAL(computeEnergyLossMean(vacuum, i, m, qOverP, q), 0);
  BOOST_CHECK_EQUAL(computeEnergyLossMode(vacuum, i, m, qOverP, q), 0);
  BOOST_CHECK_EQUAL(computeMultipleScatteringTheta0(vacuum, i, m, qOverP, q),
                    0);
}

BOOST_AUTO_TEST_SUITE_END()
