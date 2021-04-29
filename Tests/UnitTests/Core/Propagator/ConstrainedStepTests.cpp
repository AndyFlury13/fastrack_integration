// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Propagator/ConstrainedStep.hpp"

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;

namespace Acts {

namespace Test {

// This tests the implementation of the AbortList
// and the standard aborters
BOOST_AUTO_TEST_CASE(ConstrainedStepTest) {
  // forward stepping test
  ConstrainedStep stepSize_p = 0.25;

  // All of the types should be 0.25 now
  BOOST_CHECK_EQUAL(stepSize_p.values[ConstrainedStep::accuracy],
                    std::numeric_limits<double>::max());
  BOOST_CHECK_EQUAL(stepSize_p.values[ConstrainedStep::actor],
                    std::numeric_limits<double>::max());
  BOOST_CHECK_EQUAL(stepSize_p.values[ConstrainedStep::aborter],
                    std::numeric_limits<double>::max());
  BOOST_CHECK_EQUAL(stepSize_p.values[ConstrainedStep::user], 0.25);

  // Check the cast operation to double
  BOOST_CHECK_EQUAL(stepSize_p, 0.25);

  // now we update the accuracy
  stepSize_p.update(0.1, ConstrainedStep::accuracy);
  BOOST_CHECK_EQUAL(stepSize_p, 0.1);

  // now we update the actor to smaller
  stepSize_p.update(0.05, ConstrainedStep::actor);
  BOOST_CHECK_EQUAL(stepSize_p, 0.05);
  // we increase the actor, but do not release the step size
  stepSize_p.update(0.15, ConstrainedStep::actor, false);
  BOOST_CHECK_EQUAL(stepSize_p, 0.05);
  // we increase the actor, but now DO release the step size
  // it falls back to the accuracy
  stepSize_p.update(0.15, ConstrainedStep::actor, true);
  BOOST_CHECK_EQUAL(stepSize_p, 0.1);

  // now set two and update them
  stepSize_p.update(0.05, ConstrainedStep::user);
  stepSize_p.update(0.03, ConstrainedStep::accuracy);
  BOOST_CHECK_EQUAL(stepSize_p, 0.03);

  // now we release the accuracy - to the highest available value
  stepSize_p.release(ConstrainedStep::accuracy);
  BOOST_CHECK_EQUAL(stepSize_p.values[ConstrainedStep::accuracy],
                    std::numeric_limits<double>::max());
  BOOST_CHECK_EQUAL(stepSize_p, 0.05);

  // backward stepping test
  ConstrainedStep stepSize_n = -0.25;

  // All of the types should be 0.25 now
  BOOST_CHECK_EQUAL(stepSize_n.values[ConstrainedStep::accuracy],
                    -std::numeric_limits<double>::max());
  BOOST_CHECK_EQUAL(stepSize_n.values[ConstrainedStep::actor],
                    -std::numeric_limits<double>::max());
  BOOST_CHECK_EQUAL(stepSize_n.values[ConstrainedStep::aborter],
                    -std::numeric_limits<double>::max());
  BOOST_CHECK_EQUAL(stepSize_n.values[ConstrainedStep::user], -0.25);

  // Check the cast operation to double
  BOOST_CHECK_EQUAL(stepSize_n, -0.25);

  // now we update the accuracy
  stepSize_n.update(-0.1, ConstrainedStep::accuracy);
  BOOST_CHECK_EQUAL(stepSize_n, -0.1);

  // now we update the actor to smaller
  stepSize_n.update(-0.05, ConstrainedStep::actor);
  BOOST_CHECK_EQUAL(stepSize_n, -0.05);
  // we increase the actor and accuracy is smaller again, without reset
  stepSize_n.update(-0.15, ConstrainedStep::actor, false);
  BOOST_CHECK_EQUAL(stepSize_n, -0.05);
  // we increase the actor and accuracy is smaller again, with reset
  stepSize_n.update(-0.15, ConstrainedStep::actor, true);
  BOOST_CHECK_EQUAL(stepSize_n, -0.1);

  // now set two and update them
  stepSize_n.update(-0.05, ConstrainedStep::user);
  stepSize_n.update(-0.03, ConstrainedStep::accuracy);
  BOOST_CHECK_EQUAL(stepSize_n, -0.03);

  // now we release the accuracy - to the highest available value
  stepSize_n.release(ConstrainedStep::accuracy);
  BOOST_CHECK_EQUAL(stepSize_n.values[ConstrainedStep::accuracy],
                    -std::numeric_limits<double>::max());
  BOOST_CHECK_EQUAL(stepSize_n, -0.05);
}

}  // namespace Test
}  // namespace Acts
