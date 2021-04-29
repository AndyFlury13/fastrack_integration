// This file is part of the Acts project.
//
// Copyright (C) 2017-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/Helpers.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/MagneticField/MagneticField.hpp"
#include "ActsExamples/MagneticField/MagneticFieldOptions.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"

#include <string>

#include <boost/program_options.hpp>

#include "BFieldWritingBase.hpp"

namespace po = boost::program_options;

template <class T>
struct always_false : std::false_type {};

/// The main executable
///
/// Creates an InterpolatedBFieldMap from a txt or csv file and writes out the
/// grid points and values of the map into root format. The Field can then be
/// displayed using the root script printBField.cpp

int main(int argc, char* argv[]) {
  using boost::program_options::value;

  // setup and parse options
  auto desc = ActsExamples::Options::makeDefaultOptions();
  ActsExamples::Options::addMagneticFieldOptions(desc);
  desc.add_options()("bf-file-out",
                     value<std::string>()->default_value("BFieldOut.root"),
                     "Set this name for an output root file.")(
      "bf-map-out", value<std::string>()->default_value("bField"),
      "Set this name for the tree in the out file.")(
      "bf-out-rz", value<bool>()->default_value(false),
      "Please set this flag to true, if you want to print out the field map in "
      "cylinder coordinates (r,z). The default are cartesian coordinates "
      "(x,y,z). ")(
      "bf-rRange", value<ActsExamples::Options::Reals<2>>(),
      "[optional] range which the bfield map should be written out in either r "
      "(cylinder "
      "coordinates) or x/y (cartesian coordinates)  in [mm]. In case no value "
      "is handed over the whole map will be written out. Please "
      "hand over by simply seperating the values by space")(
      "bf-zRange", value<ActsExamples::Options::Reals<2>>(),
      "[optional] range which the bfield map should be written out in z in "
      "[mm].In case no value is handed over for 'bf-rRange' and 'bf-zRange the "
      "whole map will be written out. "
      "Please hand over by simply seperating the values by space")(
      "bf-rBins", value<size_t>()->default_value(200),
      "[optional] The number of bins in r. This parameter only needs to be "
      "specified if 'bf-rRange' and 'bf-zRange' are given.")(
      "bf-ZBins", value<size_t>()->default_value(300),
      "[optional] The number of bins in z. This parameter only needs to be "
      "specified if 'bf-rRange' and 'bf-zRange' are given.")(
      "bf-PhiBins", value<size_t>()->default_value(100),
      "[optional] The number of bins in phi. This parameter only needs to be "
      "specified if 'bf-rRange' and 'bf-zRange' are given and 'bf-out-rz' is "
      "turned on.");
  auto vm = ActsExamples::Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  auto bFieldVar = ActsExamples::Options::readMagneticField(vm);

  if (auto bField2D = std::dynamic_pointer_cast<
          const ActsExamples::detail::InterpolatedMagneticField2>(bFieldVar);
      bField2D) {
    ActsExamples::BField::writeField(vm, *bField2D);
    return EXIT_SUCCESS;
  } else if (auto bField3D = std::dynamic_pointer_cast<
                 const ActsExamples::detail::InterpolatedMagneticField3>(
                 bFieldVar);
             bField3D) {
    ActsExamples::BField::writeField(vm, *bField3D);
    return EXIT_SUCCESS;
  } else {
    std::cout << "Bfield map could not be read. Exiting." << std::endl;
    return EXIT_FAILURE;
  }
}
