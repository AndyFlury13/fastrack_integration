// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <bitset>

using namespace Acts::VectorHelpers;

namespace Acts {
namespace Test {

BOOST_AUTO_TEST_SUITE(Utilities)

BOOST_AUTO_TEST_CASE(bitset_to_matrix_to_bitset) {
  Eigen::Matrix<int, 4, 3> mat;
  mat << 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0;

  std::bitset<4 * 3> act = matrixToBitset(mat);
  std::bitset<4 * 3> exp{"101001011010"};

  BOOST_CHECK_EQUAL(exp, act);

  Eigen::Matrix<int, 4, 3> cnv;
  cnv = bitsetToMatrix<decltype(cnv)>(act);

  BOOST_CHECK_EQUAL(mat, cnv);
}

struct MyStruct {
  double phi() const { return 42; }
};

BOOST_AUTO_TEST_CASE(phi_helper_test) {
  Vector3 v(0, 1, 0);
  CHECK_CLOSE_ABS(phi(v), M_PI / 2., 1e-9);

  MyStruct s;
  BOOST_CHECK_EQUAL(phi(s), 42u);
}

BOOST_AUTO_TEST_CASE(perp_helper_test) {
  Vector3 v(1, 2, 3);
  CHECK_CLOSE_ABS(perp(v), std::sqrt(1 + 2 * 2), 1e-9);
}

BOOST_AUTO_TEST_CASE(theta_eta_test_helper) {
  Vector3 v(1, 2, 3);
  CHECK_CLOSE_ABS(theta(v), 0.640522, 1e-5);
  CHECK_CLOSE_ABS(eta(v), 1.10359, 1e-5);
}

BOOST_AUTO_TEST_CASE(cross_test_helper) {
  {
    Vector3 v(1, 2, 3);
    ActsMatrix<3, 3> mat;
    mat << 1, 2, 3, 4, 5, 6, 7, 8, 9;

    ActsMatrix<3, 3> act = cross(mat, v);
    ActsMatrix<3, 3> exp;
    exp << -2, -1, 0, 4, 2, 0, -2, -1, 0;

    CHECK_CLOSE_ABS(act, exp, 1e-9);
  }
}

BOOST_AUTO_TEST_CASE(toString_test_helper) {
  ActsMatrix<3, 3> mat;
  mat << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  std::string out;
  out = toString(mat);
  BOOST_CHECK(out.size() > 0);

  Translation3 trl{Vector3{1, 2, 3}};
  out = toString(trl);
  BOOST_CHECK(out.size() > 0);

  Transform3 trf;
  trf = trl;
  out = toString(trf);
  BOOST_CHECK(out.size() > 0);
}

BOOST_AUTO_TEST_CASE(shared_vector_helper_test) {
  {
    std::vector<std::shared_ptr<int>> vec;
    vec = {std::make_shared<int>(5), std::make_shared<int>(9),
           std::make_shared<int>(26), std::make_shared<int>(18473)};

    std::vector<int*> unpacked = unpack_shared_vector(vec);

    std::vector<int*> exp = {
        vec[0].get(),
        vec[1].get(),
        vec[2].get(),
        vec[3].get(),
    };

    BOOST_CHECK_EQUAL_COLLECTIONS(unpacked.begin(), unpacked.end(), exp.begin(),
                                  exp.end());
  }

  // same for const
  {
    std::vector<std::shared_ptr<const int>> vec;
    vec = {std::make_shared<const int>(5), std::make_shared<const int>(9),
           std::make_shared<const int>(26), std::make_shared<const int>(18473)};

    std::vector<const int*> unpacked = unpack_shared_vector(vec);

    std::vector<const int*> exp = {
        vec[0].get(),
        vec[1].get(),
        vec[2].get(),
        vec[3].get(),
    };

    BOOST_CHECK_EQUAL_COLLECTIONS(unpacked.begin(), unpacked.end(), exp.begin(),
                                  exp.end());
  }
}

BOOST_AUTO_TEST_CASE(VectorHelpersPosition) {
  Vector4 pos4 = Vector4::Constant(-1);
  pos4[ePos0] = 1;
  pos4[ePos1] = 2;
  pos4[ePos2] = 3;
  BOOST_CHECK_EQUAL(position(pos4), Vector3(1, 2, 3));

  FreeVector params = FreeVector::Constant(-1);
  params[eFreePos0] = 1;
  params[eFreePos1] = 2;
  params[eFreePos2] = 3;
  BOOST_CHECK_EQUAL(position(params), Vector3(1, 2, 3));
}

template <size_t I>
struct functor {
  static constexpr size_t invoke() { return I * I * I; }
};

BOOST_AUTO_TEST_CASE(test_matrix_dimension_switch) {
  constexpr size_t imax = 20;
  for (size_t i = 0; i < imax; i++) {
    size_t val = template_switch<functor, 0, imax>(i);
    BOOST_CHECK_EQUAL(val, i * i * i);
  }
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace Test
}  // namespace Acts
