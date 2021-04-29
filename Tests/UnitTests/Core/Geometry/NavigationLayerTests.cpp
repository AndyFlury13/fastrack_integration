// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/NavigationLayer.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

#include "LayerStub.hpp"

using boost::test_tools::output_test_stream;
namespace utf = boost::unit_test;

namespace Acts {

namespace Test {

// Create a test context
GeometryContext tgContext = GeometryContext();

namespace Layers {
BOOST_AUTO_TEST_SUITE(Layers)

/// Unit test for creating compliant/non-compliant NavigationLayer object
BOOST_AUTO_TEST_CASE(NavigationLayerConstruction) {
  // default constructor, copy and assignment are all deleted
  std::shared_ptr<const Surface> pSurface = Surface::makeShared<SurfaceStub>();
  auto pNavigationLayer = NavigationLayer::create(std::move(pSurface));
  BOOST_CHECK_EQUAL(pNavigationLayer->layerType(), LayerType::navigation);
  // next level: with thickness
  const double thickness = 0.1;
  auto pSurface2 = Surface::makeShared<SurfaceStub>();
  auto pThickNavigationLayer =
      NavigationLayer::create(std::move(pSurface2), thickness);
  BOOST_CHECK_EQUAL(pThickNavigationLayer->layerType(), LayerType::navigation);
}

/// Unit test for testing NavigationLayer properties
BOOST_AUTO_TEST_CASE(NavigationLayerProperties, *utf::expected_failures(1)) {
  const double thickness = 0.1;
  std::shared_ptr<const Surface> pSurface = Surface::makeShared<SurfaceStub>();
  auto rawSurfacePtr = pSurface.get();
  auto pNavigationLayer =
      NavigationLayer::create(std::move(pSurface), thickness);
  BinningValue b{BinningValue::binZ};
  Vector3 origin{0., 0., 0.};
  // binningPosition(), needs a better test
  BOOST_CHECK_EQUAL(pNavigationLayer->binningPosition(tgContext, b), origin);
  // surfaceRepresentation() [looks dangerous]
  BOOST_CHECK_EQUAL(rawSurfacePtr,
                    &(pNavigationLayer->surfaceRepresentation()));
  // isOnLayer()
  BOOST_CHECK(pNavigationLayer->isOnLayer(tgContext, origin, true));
  // isOnLayer()
  Vector3 crazyPosition{1000., 10000., std::nan("")};
  BOOST_CHECK(!pNavigationLayer->isOnLayer(tgContext, crazyPosition, true));
  // resolve()
  BOOST_CHECK(!pNavigationLayer->resolve(true, true, true));
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace Layers
}  // namespace Test

}  // namespace Acts
