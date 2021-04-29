// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Tests/CommonHelpers/GenerateParameters.hpp"
#include "Acts/Tests/CommonHelpers/TestSourceLink.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

#include <numeric>
#include <optional>
#include <random>
#include <tuple>

namespace {

using namespace Acts;
using namespace Acts::UnitLiterals;
using namespace Acts::Test;
namespace bd = boost::unit_test::data;

using ParametersVector = BoundTrackParameters::ParametersVector;
using CovarianceMatrix = BoundTrackParameters::CovarianceMatrix;
using Jacobian = BoundMatrix;

struct TestTrackState {
  std::shared_ptr<Surface> surface;
  TestSourceLink sourceLink;
  BoundTrackParameters predicted;
  BoundTrackParameters filtered;
  BoundTrackParameters smoothed;
  Jacobian jacobian;
  double chi2;
  double pathLength;

  // Generate a random TestTrackState.
  //
  // @param rng Random number generator
  // @param size_t nMeasurement either 1 or 2
  template <typename rng_t>
  TestTrackState(rng_t& rng, size_t nMeasurements)
      : surface(Surface::makeShared<PlaneSurface>(Vector3::Zero(),
                                                  Vector3::UnitZ())),
        // set bogus parameters first since they are not default-constructible
        predicted(surface, BoundVector::Zero()),
        filtered(surface, BoundVector::Zero()),
        smoothed(surface, BoundVector::Zero()),
        jacobian(Jacobian::Identity()),
        chi2(std::chi_squared_distribution<double>(nMeasurements)(rng)),
        pathLength(
            std::uniform_real_distribution<ActsScalar>(1_mm, 10_mm)(rng)) {
    // set a random geometry identifier to uniquely identify each surface
    auto geoId =
        std::uniform_int_distribution<GeometryIdentifier::Value>()(rng);
    surface->assignGeometryId(geoId);

    // create source link w/ inline 1d or 2d measurement data
    if (nMeasurements == 1u) {
      auto [par, cov] = generateParametersCovariance<ActsScalar, 1u>(rng);
      sourceLink = TestSourceLink(eBoundLoc0, par[0], cov(0, 0), geoId);
    } else if (nMeasurements == 2u) {
      auto [par, cov] = generateParametersCovariance<ActsScalar, 2u>(rng);
      sourceLink = TestSourceLink(eBoundLoc1, eBoundQOverP, par, cov, geoId);
    } else {
      throw std::runtime_error("invalid number of measurement dimensions");
    }

    // create track parameters
    auto [trkPar, trkCov] = generateBoundParametersCovariance(rng);
    // trkPar[eBoundPhi] = 45_degree;
    // trkPar[eBoundTheta] = 90_degree;
    // trkPar[eBoundQOverP] = 5.;
    // predicted
    predicted = BoundTrackParameters(surface, trkPar, trkCov);
    // filtered, modified q/p, reduced covariance
    // trkPar[eBoundQOverP] = 10.;
    filtered = BoundTrackParameters(surface, trkPar, 0.75 * trkCov);
    // smoothed, modified q/p, further reduced covariance
    // trkPar[eBoundQOverP] = 15.;
    smoothed = BoundTrackParameters(surface, trkPar, 0.5 * trkCov);

    // propagation jacobian is identity + corrections
    for (Eigen::Index c = 0; c < jacobian.cols(); ++c) {
      for (Eigen::Index r = 0; r < jacobian.rows(); ++r) {
        jacobian(c, r) +=
            std::uniform_real_distribution<ActsScalar>(-0.125, 0.125)(rng);
      }
    }
  }
};

// Fill a TrackStateProxy with values from a TestTrackState.
//
// @param[in] pc TestTrackState with the input values
// @param[in] mask Specifies which components are used/filled
// @param[out] ts TrackStateProxy which is filled
// @param [in] nMeasurements Dimension of the measurement
template <typename track_state_t>
void fillTrackState(const TestTrackState& pc, TrackStatePropMask mask,
                    track_state_t& ts) {
  // always set the reference surface
  ts.setReferenceSurface(pc.predicted.referenceSurface().getSharedPtr());

  if (ACTS_CHECK_BIT(mask, TrackStatePropMask::Predicted)) {
    ts.predicted() = pc.predicted.parameters();
    ts.predictedCovariance() = *(pc.predicted.covariance());
  }
  if (ACTS_CHECK_BIT(mask, TrackStatePropMask::Filtered)) {
    ts.filtered() = pc.filtered.parameters();
    ts.filteredCovariance() = *(pc.filtered.covariance());
  }
  if (ACTS_CHECK_BIT(mask, TrackStatePropMask::Smoothed)) {
    ts.smoothed() = pc.smoothed.parameters();
    ts.smoothedCovariance() = *(pc.smoothed.covariance());
  }
  if (ACTS_CHECK_BIT(mask, TrackStatePropMask::Jacobian)) {
    ts.jacobian() = pc.jacobian;
  }
  ts.chi2() = pc.chi2;
  ts.pathLength() = pc.pathLength;
  // source link defines the uncalibrated measurement
  if (ACTS_CHECK_BIT(mask, TrackStatePropMask::Uncalibrated)) {
    ts.uncalibrated() = pc.sourceLink;
  }
  // create calibrated measurements from source link
  if (ACTS_CHECK_BIT(mask, TrackStatePropMask::Calibrated)) {
    auto meas = TestSourceLinkCalibrator()(pc.sourceLink, nullptr);
    std::visit([&](const auto& m) { ts.setCalibrated(m); }, meas);
  }
}

const GeometryContext gctx;
// fixed seed for reproducible tests
std::default_random_engine rng(31415);

}  // namespace

BOOST_AUTO_TEST_SUITE(EventDataMultiTrajectory)

BOOST_AUTO_TEST_CASE(Build) {
  constexpr TrackStatePropMask kMask = TrackStatePropMask::Predicted;

  // construct trajectory w/ multiple components
  MultiTrajectory<TestSourceLink> t;
  auto i0 = t.addTrackState(kMask);
  // trajectory bifurcates here into multiple hypotheses
  auto i1a = t.addTrackState(kMask, i0);
  auto i1b = t.addTrackState(kMask, i0);
  auto i2a = t.addTrackState(kMask, i1a);
  auto i2b = t.addTrackState(kMask, i1b);

  // print each trajectory component
  std::vector<size_t> act;
  auto collect = [&](auto p) {
    act.push_back(p.index());
    BOOST_CHECK(!p.hasUncalibrated());
    BOOST_CHECK(!p.hasCalibrated());
    BOOST_CHECK(!p.hasFiltered());
    BOOST_CHECK(!p.hasSmoothed());
    BOOST_CHECK(!p.hasJacobian());
    BOOST_CHECK(!p.hasProjector());
  };

  std::vector<size_t> exp = {i2a, i1a, i0};
  t.visitBackwards(i2a, collect);
  BOOST_CHECK_EQUAL_COLLECTIONS(act.begin(), act.end(), exp.begin(), exp.end());

  act.clear();
  exp = {i2b, i1b, i0};
  t.visitBackwards(i2b, collect);
  BOOST_CHECK_EQUAL_COLLECTIONS(act.begin(), act.end(), exp.begin(), exp.end());

  act.clear();
  t.applyBackwards(i2b, collect);
  BOOST_CHECK_EQUAL_COLLECTIONS(act.begin(), act.end(), exp.begin(), exp.end());
}

BOOST_AUTO_TEST_CASE(ApplyWithAbort) {
  constexpr TrackStatePropMask kMask = TrackStatePropMask::Predicted;

  // construct trajectory with three components
  MultiTrajectory<TestSourceLink> t;
  auto i0 = t.addTrackState(kMask);
  auto i1 = t.addTrackState(kMask, i0);
  auto i2 = t.addTrackState(kMask, i1);

  size_t n = 0;
  t.applyBackwards(i2, [&](const auto&) {
    n++;
    return false;
  });
  BOOST_CHECK_EQUAL(n, 1u);

  n = 0;
  t.applyBackwards(i2, [&](const auto& ts) {
    n++;
    if (ts.index() == i1) {
      return false;
    }
    return true;
  });
  BOOST_CHECK_EQUAL(n, 2u);

  n = 0;
  t.applyBackwards(i2, [&](const auto&) {
    n++;
    return true;
  });
  BOOST_CHECK_EQUAL(n, 3u);
}

BOOST_AUTO_TEST_CASE(BitmaskOperators) {
  using PM = TrackStatePropMask;

  auto bs1 = PM::Uncalibrated;

  BOOST_CHECK(ACTS_CHECK_BIT(bs1, PM::Uncalibrated));
  BOOST_CHECK(!ACTS_CHECK_BIT(bs1, PM::Calibrated));

  auto bs2 = PM::Calibrated;

  BOOST_CHECK(!ACTS_CHECK_BIT(bs2, PM::Uncalibrated));
  BOOST_CHECK(ACTS_CHECK_BIT(bs2, PM::Calibrated));

  auto bs3 = PM::Calibrated | PM::Uncalibrated;

  BOOST_CHECK(ACTS_CHECK_BIT(bs3, PM::Uncalibrated));
  BOOST_CHECK(ACTS_CHECK_BIT(bs3, PM::Calibrated));

  BOOST_CHECK(ACTS_CHECK_BIT(PM::All, PM::Uncalibrated));
  BOOST_CHECK(ACTS_CHECK_BIT(PM::All, PM::Calibrated));

  auto bs4 = PM::Predicted | PM::Jacobian | PM::Uncalibrated;
  BOOST_CHECK(ACTS_CHECK_BIT(bs4, PM::Predicted));
  BOOST_CHECK(ACTS_CHECK_BIT(bs4, PM::Uncalibrated));
  BOOST_CHECK(ACTS_CHECK_BIT(bs4, PM::Jacobian));
  BOOST_CHECK(!ACTS_CHECK_BIT(bs4, PM::Calibrated));
  BOOST_CHECK(!ACTS_CHECK_BIT(bs4, PM::Filtered));
  BOOST_CHECK(!ACTS_CHECK_BIT(bs4, PM::Smoothed));

  auto cnv = [](auto a) -> std::bitset<8> {
    return static_cast<std::underlying_type<PM>::type>(a);
  };

  BOOST_CHECK(cnv(PM::All).all());    // all ones
  BOOST_CHECK(cnv(PM::None).none());  // all zeros

  // test orthogonality
  std::array<PM, 6> values{PM::Predicted, PM::Filtered,     PM::Smoothed,
                           PM::Jacobian,  PM::Uncalibrated, PM::Calibrated};
  for (size_t i = 0; i < values.size(); i++) {
    for (size_t j = 0; j < values.size(); j++) {
      PM a = values[i];
      PM b = values[j];

      if (i == j) {
        BOOST_CHECK(cnv(a & b).count() == 1);
      } else {
        BOOST_CHECK(cnv(a & b).none());
      }
    }
  }

  BOOST_CHECK(cnv(PM::Predicted ^ PM::Filtered).count() == 2);
  BOOST_CHECK(cnv(PM::Predicted ^ PM::Predicted).none());
  BOOST_CHECK(~(PM::Predicted | PM::Calibrated) ==
              (PM::All ^ PM::Predicted ^ PM::Calibrated));

  PM base = PM::None;
  BOOST_CHECK(cnv(base) == 0);

  base &= PM::Filtered;
  BOOST_CHECK(cnv(base) == 0);

  base |= PM::Filtered;
  BOOST_CHECK(base == PM::Filtered);

  base |= PM::Calibrated;
  BOOST_CHECK(base == (PM::Filtered | PM::Calibrated));

  base ^= PM::All;
  BOOST_CHECK(base == ~(PM::Filtered | PM::Calibrated));
}

BOOST_AUTO_TEST_CASE(AddTrackStateWithBitMask) {
  using PM = TrackStatePropMask;

  MultiTrajectory<TestSourceLink> t;

  auto ts = t.getTrackState(t.addTrackState(PM::All));
  BOOST_CHECK(ts.hasPredicted());
  BOOST_CHECK(ts.hasFiltered());
  BOOST_CHECK(ts.hasSmoothed());
  BOOST_CHECK(ts.hasUncalibrated());
  BOOST_CHECK(ts.hasCalibrated());
  BOOST_CHECK(ts.hasProjector());
  BOOST_CHECK(ts.hasJacobian());

  ts = t.getTrackState(t.addTrackState(PM::None));
  BOOST_CHECK(!ts.hasPredicted());
  BOOST_CHECK(!ts.hasFiltered());
  BOOST_CHECK(!ts.hasSmoothed());
  BOOST_CHECK(!ts.hasUncalibrated());
  BOOST_CHECK(!ts.hasCalibrated());
  BOOST_CHECK(!ts.hasProjector());
  BOOST_CHECK(!ts.hasJacobian());

  ts = t.getTrackState(t.addTrackState(PM::Predicted));
  BOOST_CHECK(ts.hasPredicted());
  BOOST_CHECK(!ts.hasFiltered());
  BOOST_CHECK(!ts.hasSmoothed());
  BOOST_CHECK(!ts.hasUncalibrated());
  BOOST_CHECK(!ts.hasCalibrated());
  BOOST_CHECK(!ts.hasProjector());
  BOOST_CHECK(!ts.hasJacobian());

  ts = t.getTrackState(t.addTrackState(PM::Filtered));
  BOOST_CHECK(!ts.hasPredicted());
  BOOST_CHECK(ts.hasFiltered());
  BOOST_CHECK(!ts.hasSmoothed());
  BOOST_CHECK(!ts.hasUncalibrated());
  BOOST_CHECK(!ts.hasCalibrated());
  BOOST_CHECK(!ts.hasProjector());
  BOOST_CHECK(!ts.hasJacobian());

  ts = t.getTrackState(t.addTrackState(PM::Smoothed));
  BOOST_CHECK(!ts.hasPredicted());
  BOOST_CHECK(!ts.hasFiltered());
  BOOST_CHECK(ts.hasSmoothed());
  BOOST_CHECK(!ts.hasUncalibrated());
  BOOST_CHECK(!ts.hasCalibrated());
  BOOST_CHECK(!ts.hasProjector());
  BOOST_CHECK(!ts.hasJacobian());

  ts = t.getTrackState(t.addTrackState(PM::Uncalibrated));
  BOOST_CHECK(!ts.hasPredicted());
  BOOST_CHECK(!ts.hasFiltered());
  BOOST_CHECK(!ts.hasSmoothed());
  BOOST_CHECK(ts.hasUncalibrated());
  BOOST_CHECK(!ts.hasCalibrated());
  BOOST_CHECK(!ts.hasProjector());
  BOOST_CHECK(!ts.hasJacobian());

  ts = t.getTrackState(t.addTrackState(PM::Calibrated));
  BOOST_CHECK(!ts.hasPredicted());
  BOOST_CHECK(!ts.hasFiltered());
  BOOST_CHECK(!ts.hasSmoothed());
  BOOST_CHECK(!ts.hasUncalibrated());
  BOOST_CHECK(ts.hasCalibrated());
  BOOST_CHECK(ts.hasProjector());
  BOOST_CHECK(!ts.hasJacobian());

  ts = t.getTrackState(t.addTrackState(PM::Jacobian));
  BOOST_CHECK(!ts.hasPredicted());
  BOOST_CHECK(!ts.hasFiltered());
  BOOST_CHECK(!ts.hasSmoothed());
  BOOST_CHECK(!ts.hasUncalibrated());
  BOOST_CHECK(!ts.hasCalibrated());
  BOOST_CHECK(!ts.hasProjector());
  BOOST_CHECK(ts.hasJacobian());
}

// assert expected "cross-talk" between trackstate proxies
BOOST_AUTO_TEST_CASE(TrackStateProxyCrossTalk) {
  TestTrackState pc(rng, 2u);

  // multi trajectory w/ a single, fully set, track state
  MultiTrajectory<TestSourceLink> traj;
  size_t index = traj.addTrackState();
  {
    auto ts = traj.getTrackState(index);
    fillTrackState(pc, TrackStatePropMask::All, ts);
  }
  // get two TrackStateProxies that reference the same data
  auto tsa = traj.getTrackState(index);
  auto tsb = traj.getTrackState(index);
  // then modify one and check that the other was modified as well
  {
    auto [par, cov] = generateBoundParametersCovariance(rng);
    tsb.predicted() = par;
    tsb.predictedCovariance() = cov;
    BOOST_CHECK_EQUAL(tsa.predicted(), par);
    BOOST_CHECK_EQUAL(tsa.predictedCovariance(), cov);
    BOOST_CHECK_EQUAL(tsb.predicted(), par);
    BOOST_CHECK_EQUAL(tsb.predictedCovariance(), cov);
  }
  {
    auto [par, cov] = generateBoundParametersCovariance(rng);
    tsb.filtered() = par;
    tsb.filteredCovariance() = cov;
    BOOST_CHECK_EQUAL(tsa.filtered(), par);
    BOOST_CHECK_EQUAL(tsa.filteredCovariance(), cov);
    BOOST_CHECK_EQUAL(tsb.filtered(), par);
    BOOST_CHECK_EQUAL(tsb.filteredCovariance(), cov);
  }
  {
    auto [par, cov] = generateBoundParametersCovariance(rng);
    tsb.smoothed() = par;
    tsb.smoothedCovariance() = cov;
    BOOST_CHECK_EQUAL(tsa.smoothed(), par);
    BOOST_CHECK_EQUAL(tsa.smoothedCovariance(), cov);
    BOOST_CHECK_EQUAL(tsb.smoothed(), par);
    BOOST_CHECK_EQUAL(tsb.smoothedCovariance(), cov);
  }
  {
    // create a new (invalid) source link
    TestSourceLink invalid;
    BOOST_CHECK_NE(tsa.uncalibrated(), invalid);
    BOOST_CHECK_NE(tsb.uncalibrated(), invalid);
    tsb.uncalibrated() = invalid;
    BOOST_CHECK_EQUAL(tsa.uncalibrated(), invalid);
    BOOST_CHECK_EQUAL(tsb.uncalibrated(), invalid);
  }
  {
    // reset measurements w/ full parameters
    auto [measPar, measCov] = generateBoundParametersCovariance(rng);
    tsb.calibrated() = measPar;
    tsb.calibratedCovariance() = measCov;
    BOOST_CHECK_EQUAL(tsa.calibrated(), measPar);
    BOOST_CHECK_EQUAL(tsa.calibratedCovariance(), measCov);
    BOOST_CHECK_EQUAL(tsb.calibrated(), measPar);
    BOOST_CHECK_EQUAL(tsb.calibratedCovariance(), measCov);
  }
  {
    // reset only the effective measurements
    auto [measPar, measCov] = generateBoundParametersCovariance(rng);
    size_t nMeasurements = tsb.effectiveCalibrated().rows();
    auto effPar = measPar.head(nMeasurements);
    auto effCov = measCov.topLeftCorner(nMeasurements, nMeasurements);
    tsb.effectiveCalibrated() = effPar;
    tsb.effectiveCalibratedCovariance() = effCov;
    BOOST_CHECK_EQUAL(tsa.effectiveCalibrated(), effPar);
    BOOST_CHECK_EQUAL(tsa.effectiveCalibratedCovariance(), effCov);
    BOOST_CHECK_EQUAL(tsb.effectiveCalibrated(), effPar);
    BOOST_CHECK_EQUAL(tsb.effectiveCalibratedCovariance(), effCov);
  }
  {
    Jacobian jac = Jacobian::Identity();
    BOOST_CHECK_NE(tsa.jacobian(), jac);
    BOOST_CHECK_NE(tsb.jacobian(), jac);
    tsb.jacobian() = jac;
    BOOST_CHECK_EQUAL(tsa.jacobian(), jac);
    BOOST_CHECK_EQUAL(tsb.jacobian(), jac);
  }
  {
    tsb.chi2() = 98.0;
    BOOST_CHECK_EQUAL(tsa.chi2(), 98.0);
    BOOST_CHECK_EQUAL(tsb.chi2(), 98.0);
  }
  {
    tsb.pathLength() = 66.0;
    BOOST_CHECK_EQUAL(tsa.pathLength(), 66.0);
    BOOST_CHECK_EQUAL(tsb.pathLength(), 66.0);
  }
}

BOOST_AUTO_TEST_CASE(TrackStateReassignment) {
  TestTrackState pc(rng, 1u);

  MultiTrajectory<TestSourceLink> t;
  size_t index = t.addTrackState();
  auto ts = t.getTrackState(index);
  fillTrackState(pc, TrackStatePropMask::All, ts);

  // assert contents of original measurement (just to be safe)
  BOOST_CHECK_EQUAL(ts.calibratedSize(), 1u);
  BOOST_CHECK_EQUAL(ts.effectiveCalibrated(),
                    (pc.sourceLink.parameters.head<1>()));
  BOOST_CHECK_EQUAL(ts.effectiveCalibratedCovariance(),
                    (pc.sourceLink.covariance.topLeftCorner<1, 1>()));

  // use temporary measurement to reset calibrated data
  TestTrackState ttsb(rng, 2u);
  auto meas = TestSourceLinkCalibrator()(ttsb.sourceLink, nullptr);
  auto m2 = std::get<Measurement<TestSourceLink, BoundIndices, 2u>>(meas);
  ts.setCalibrated(m2);

  BOOST_CHECK_EQUAL(ts.calibratedSize(), m2.size());
  BOOST_CHECK_EQUAL(ts.effectiveCalibrated(), m2.parameters());
  BOOST_CHECK_EQUAL(ts.effectiveCalibratedCovariance(), m2.covariance());
  BOOST_CHECK_EQUAL(ts.effectiveProjector(), m2.projector());

  // check that the overallocated parts are zeroed
  ParametersVector mParFull = ParametersVector::Zero();
  CovarianceMatrix mCovFull = CovarianceMatrix::Zero();
  ActsMatrix<MultiTrajectory<TestSourceLink>::MeasurementSizeMax, eBoundSize>
      projFull;
  mParFull.head<2>() = m2.parameters();
  mCovFull.topLeftCorner<2, 2>() = m2.covariance();
  projFull.setZero();
  projFull.topLeftCorner<2, eBoundSize>() = m2.projector();
  BOOST_CHECK_EQUAL(ts.calibrated(), mParFull);
  BOOST_CHECK_EQUAL(ts.calibratedCovariance(), mCovFull);
  BOOST_CHECK_EQUAL(ts.projector(), projFull);
}

BOOST_DATA_TEST_CASE(TrackStateProxyStorage, bd::make({1u, 2u}),
                     nMeasurements) {
  TestTrackState pc(rng, nMeasurements);

  // create trajectory with a single fully-filled random track state
  MultiTrajectory<TestSourceLink> t;
  size_t index = t.addTrackState();
  auto ts = t.getTrackState(index);
  fillTrackState(pc, TrackStatePropMask::All, ts);

  // check that the surface is correctly set
  BOOST_CHECK_EQUAL(&ts.referenceSurface(), pc.surface.get());
  BOOST_CHECK_EQUAL(ts.referenceSurface().geometryId(),
                    pc.sourceLink.geometryId());

  // check that the track parameters are set
  BOOST_CHECK(ts.hasPredicted());
  BOOST_CHECK_EQUAL(ts.predicted(), pc.predicted.parameters());
  BOOST_CHECK_EQUAL(ts.predictedCovariance(), *pc.predicted.covariance());
  BOOST_CHECK(ts.hasFiltered());
  BOOST_CHECK_EQUAL(ts.filtered(), pc.filtered.parameters());
  BOOST_CHECK_EQUAL(ts.filteredCovariance(), *pc.filtered.covariance());
  BOOST_CHECK(ts.hasSmoothed());
  BOOST_CHECK_EQUAL(ts.smoothed(), pc.smoothed.parameters());
  BOOST_CHECK_EQUAL(ts.smoothedCovariance(), *pc.smoothed.covariance());

  // check that the jacobian is set
  BOOST_CHECK(ts.hasJacobian());
  BOOST_CHECK_EQUAL(ts.jacobian(), pc.jacobian);
  BOOST_CHECK_EQUAL(ts.pathLength(), pc.pathLength);
  // check that chi2 is set
  BOOST_CHECK_EQUAL(ts.chi2(), pc.chi2);

  // check that the uncalibrated source link is set
  BOOST_CHECK(ts.hasUncalibrated());
  BOOST_CHECK_EQUAL(ts.uncalibrated(), pc.sourceLink);

  // check that the calibrated measurement is set
  BOOST_CHECK(ts.hasCalibrated());
  BOOST_CHECK_EQUAL(ts.calibratedSourceLink(), pc.sourceLink);
  BOOST_CHECK_EQUAL(ts.effectiveCalibrated(),
                    pc.sourceLink.parameters.head(nMeasurements));
  BOOST_CHECK_EQUAL(
      ts.effectiveCalibratedCovariance(),
      pc.sourceLink.covariance.topLeftCorner(nMeasurements, nMeasurements));
  {
    ParametersVector mParFull = ParametersVector::Zero();
    CovarianceMatrix mCovFull = CovarianceMatrix::Zero();
    mParFull.head(nMeasurements) = pc.sourceLink.parameters.head(nMeasurements);
    mCovFull.topLeftCorner(nMeasurements, nMeasurements) =
        pc.sourceLink.covariance.topLeftCorner(nMeasurements, nMeasurements);
    BOOST_CHECK_EQUAL(ts.calibrated(), mParFull);
    BOOST_CHECK_EQUAL(ts.calibratedCovariance(), mCovFull);
  }

  BOOST_CHECK(ts.hasProjector());
  ActsMatrix<MultiTrajectory<TestSourceLink>::MeasurementSizeMax, eBoundSize>
      fullProj;
  fullProj.setZero();
  {
    // create a temporary measurement to extract the projector matrix
    auto meas = TestSourceLinkCalibrator()(pc.sourceLink, nullptr);
    std::visit(
        [&](const auto& m) {
          fullProj.topLeftCorner(nMeasurements, eBoundSize) = m.projector();
        },
        meas);
  }
  BOOST_CHECK_EQUAL(ts.effectiveProjector(),
                    fullProj.topLeftCorner(nMeasurements, eBoundSize));
  BOOST_CHECK_EQUAL(ts.projector(), fullProj);
}

BOOST_AUTO_TEST_CASE(TrackStateProxyAllocations) {
  TestTrackState pc(rng, 2u);

  // this should allocate for all components in the trackstate, plus filtered
  MultiTrajectory<TestSourceLink> t;
  size_t i = t.addTrackState(
      TrackStatePropMask::Predicted | TrackStatePropMask::Filtered |
      TrackStatePropMask::Uncalibrated | TrackStatePropMask::Jacobian);
  auto tso = t.getTrackState(i);
  fillTrackState(pc, TrackStatePropMask::Predicted, tso);
  fillTrackState(pc, TrackStatePropMask::Filtered, tso);
  fillTrackState(pc, TrackStatePropMask::Uncalibrated, tso);
  fillTrackState(pc, TrackStatePropMask::Jacobian, tso);

  BOOST_CHECK(tso.hasPredicted());
  BOOST_CHECK(tso.hasFiltered());
  BOOST_CHECK(!tso.hasSmoothed());
  BOOST_CHECK(tso.hasUncalibrated());
  BOOST_CHECK(!tso.hasCalibrated());
  BOOST_CHECK(tso.hasJacobian());
}

BOOST_AUTO_TEST_CASE(TrackStateProxyGetMask) {
  using PM = TrackStatePropMask;

  std::array<PM, 6> values{PM::Predicted, PM::Filtered,     PM::Smoothed,
                           PM::Jacobian,  PM::Uncalibrated, PM::Calibrated};
  PM all = std::accumulate(values.begin(), values.end(), PM::None,
                           [](auto a, auto b) { return a | b; });

  MultiTrajectory<TestSourceLink> mj;
  {
    auto ts = mj.getTrackState(mj.addTrackState(PM::All));
    BOOST_CHECK(ts.getMask() == all);
  }
  {
    auto ts = mj.getTrackState(mj.addTrackState(PM::Filtered | PM::Calibrated));
    BOOST_CHECK(ts.getMask() == (PM::Filtered | PM::Calibrated));
  }
  {
    auto ts = mj.getTrackState(
        mj.addTrackState(PM::Filtered | PM::Smoothed | PM::Predicted));
    BOOST_CHECK(ts.getMask() == (PM::Filtered | PM::Smoothed | PM::Predicted));
  }
  {
    for (PM mask : values) {
      auto ts = mj.getTrackState(mj.addTrackState(mask));
      BOOST_CHECK(ts.getMask() == mask);
    }
  }
}

BOOST_AUTO_TEST_CASE(TrackStateProxyCopy) {
  using PM = TrackStatePropMask;

  std::array<PM, 6> values{PM::Predicted, PM::Filtered,     PM::Smoothed,
                           PM::Jacobian,  PM::Uncalibrated, PM::Calibrated};

  MultiTrajectory<TestSourceLink> mj;
  auto mkts = [&](PM mask) { return mj.getTrackState(mj.addTrackState(mask)); };

  // orthogonal ones
  for (PM a : values) {
    for (PM b : values) {
      auto tsa = mkts(a);
      auto tsb = mkts(b);
      // doesn't work
      if (a != b) {
        BOOST_CHECK_THROW(tsa.copyFrom(tsb), std::runtime_error);
        BOOST_CHECK_THROW(tsb.copyFrom(tsa), std::runtime_error);
      } else {
        tsa.copyFrom(tsb);
        tsb.copyFrom(tsa);
      }
    }
  }

  auto ts1 = mkts(PM::Filtered | PM::Predicted);  // this has both
  ts1.filtered().setRandom();
  ts1.filteredCovariance().setRandom();
  ts1.predicted().setRandom();
  ts1.predictedCovariance().setRandom();

  // ((src XOR dst) & src) == 0
  auto ts2 = mkts(PM::Predicted);
  ts2.predicted().setRandom();
  ts2.predictedCovariance().setRandom();

  // they are different before
  BOOST_CHECK(ts1.predicted() != ts2.predicted());
  BOOST_CHECK(ts1.predictedCovariance() != ts2.predictedCovariance());

  // ts1 -> ts2 fails
  BOOST_CHECK_THROW(ts2.copyFrom(ts1), std::runtime_error);
  BOOST_CHECK(ts1.predicted() != ts2.predicted());
  BOOST_CHECK(ts1.predictedCovariance() != ts2.predictedCovariance());

  // ts2 -> ts1 is ok
  ts1.copyFrom(ts2);
  BOOST_CHECK(ts1.predicted() == ts2.predicted());
  BOOST_CHECK(ts1.predictedCovariance() == ts2.predictedCovariance());

  size_t i0 = mj.addTrackState();
  size_t i1 = mj.addTrackState();
  ts1 = mj.getTrackState(i0);
  ts2 = mj.getTrackState(i1);
  TestTrackState rts1(rng, 1u);
  TestTrackState rts2(rng, 2u);
  fillTrackState(rts1, TrackStatePropMask::All, ts1);
  fillTrackState(rts2, TrackStatePropMask::All, ts2);

  auto ots1 = mkts(PM::All);
  auto ots2 = mkts(PM::All);
  // make full copy for later. We prove full copy works right below
  ots1.copyFrom(ts1);
  ots2.copyFrom(ts2);

  BOOST_CHECK_NE(ts1.predicted(), ts2.predicted());
  BOOST_CHECK_NE(ts1.predictedCovariance(), ts2.predictedCovariance());
  BOOST_CHECK_NE(ts1.filtered(), ts2.filtered());
  BOOST_CHECK_NE(ts1.filteredCovariance(), ts2.filteredCovariance());
  BOOST_CHECK_NE(ts1.smoothed(), ts2.smoothed());
  BOOST_CHECK_NE(ts1.smoothedCovariance(), ts2.smoothedCovariance());

  BOOST_CHECK_NE(ts1.uncalibrated(), ts2.uncalibrated());

  BOOST_CHECK_NE(ts1.calibratedSourceLink(), ts2.calibratedSourceLink());
  BOOST_CHECK_NE(ts1.calibrated(), ts2.calibrated());
  BOOST_CHECK_NE(ts1.calibratedCovariance(), ts2.calibratedCovariance());
  BOOST_CHECK_NE(ts1.calibratedSize(), ts2.calibratedSize());
  BOOST_CHECK_NE(ts1.projector(), ts2.projector());

  BOOST_CHECK_NE(ts1.jacobian(), ts2.jacobian());
  BOOST_CHECK_NE(ts1.chi2(), ts2.chi2());
  BOOST_CHECK_NE(ts1.pathLength(), ts2.pathLength());
  BOOST_CHECK_NE(&ts1.referenceSurface(), &ts2.referenceSurface());

  ts1.copyFrom(ts2);

  BOOST_CHECK_EQUAL(ts1.predicted(), ts2.predicted());
  BOOST_CHECK_EQUAL(ts1.predictedCovariance(), ts2.predictedCovariance());
  BOOST_CHECK_EQUAL(ts1.filtered(), ts2.filtered());
  BOOST_CHECK_EQUAL(ts1.filteredCovariance(), ts2.filteredCovariance());
  BOOST_CHECK_EQUAL(ts1.smoothed(), ts2.smoothed());
  BOOST_CHECK_EQUAL(ts1.smoothedCovariance(), ts2.smoothedCovariance());

  BOOST_CHECK_EQUAL(ts1.uncalibrated(), ts2.uncalibrated());

  BOOST_CHECK_EQUAL(ts1.calibratedSourceLink(), ts2.calibratedSourceLink());
  BOOST_CHECK_EQUAL(ts1.calibrated(), ts2.calibrated());
  BOOST_CHECK_EQUAL(ts1.calibratedCovariance(), ts2.calibratedCovariance());
  BOOST_CHECK_EQUAL(ts1.calibratedSize(), ts2.calibratedSize());
  BOOST_CHECK_EQUAL(ts1.projector(), ts2.projector());

  BOOST_CHECK_EQUAL(ts1.jacobian(), ts2.jacobian());
  BOOST_CHECK_EQUAL(ts1.chi2(), ts2.chi2());
  BOOST_CHECK_EQUAL(ts1.pathLength(), ts2.pathLength());
  BOOST_CHECK_EQUAL(&ts1.referenceSurface(), &ts2.referenceSurface());

  // full copy proven to work. now let's do partial copy
  ts2 = mkts(PM::Predicted | PM::Jacobian | PM::Calibrated);
  ts2.copyFrom(ots2, PM::Predicted | PM::Jacobian | PM::Calibrated);
  // copy into empty ts, only copy some
  ts1.copyFrom(ots1);  // reset to original
  // is different again
  BOOST_CHECK_NE(ts1.predicted(), ts2.predicted());
  BOOST_CHECK_NE(ts1.predictedCovariance(), ts2.predictedCovariance());

  BOOST_CHECK_NE(ts1.calibratedSourceLink(), ts2.calibratedSourceLink());
  BOOST_CHECK_NE(ts1.calibrated(), ts2.calibrated());
  BOOST_CHECK_NE(ts1.calibratedCovariance(), ts2.calibratedCovariance());
  BOOST_CHECK_NE(ts1.calibratedSize(), ts2.calibratedSize());
  BOOST_CHECK_NE(ts1.projector(), ts2.projector());

  BOOST_CHECK_NE(ts1.jacobian(), ts2.jacobian());
  BOOST_CHECK_NE(ts1.chi2(), ts2.chi2());
  BOOST_CHECK_NE(ts1.pathLength(), ts2.pathLength());
  BOOST_CHECK_NE(&ts1.referenceSurface(), &ts2.referenceSurface());

  ts1.copyFrom(ts2);

  // some components are same now
  BOOST_CHECK_EQUAL(ts1.predicted(), ts2.predicted());
  BOOST_CHECK_EQUAL(ts1.predictedCovariance(), ts2.predictedCovariance());

  BOOST_CHECK_EQUAL(ts1.calibratedSourceLink(), ts2.calibratedSourceLink());
  BOOST_CHECK_EQUAL(ts1.calibrated(), ts2.calibrated());
  BOOST_CHECK_EQUAL(ts1.calibratedCovariance(), ts2.calibratedCovariance());
  BOOST_CHECK_EQUAL(ts1.calibratedSize(), ts2.calibratedSize());
  BOOST_CHECK_EQUAL(ts1.projector(), ts2.projector());

  BOOST_CHECK_EQUAL(ts1.jacobian(), ts2.jacobian());
  BOOST_CHECK_EQUAL(ts1.chi2(), ts2.chi2());              // always copied
  BOOST_CHECK_EQUAL(ts1.pathLength(), ts2.pathLength());  // always copied
  BOOST_CHECK_EQUAL(&ts1.referenceSurface(),
                    &ts2.referenceSurface());  // always copied
}

BOOST_AUTO_TEST_SUITE_END()
