// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// System include(s)
#include <string>

// SYCL include
#include <CL/sycl.hpp>

namespace Acts::Sycl {
/// @brief Custom device selector that refuses to select NVIDIA OpenCL backends.
///
/// It is also possible to the tell the selector explicitly which device we want
/// to use by providing a substring of the preferred device's name.
struct DeviceSelector : public cl::sycl::device_selector {
  DeviceSelector(const std::string& deviceName = "");

  int operator()(const cl::sycl::device& d) const;

 private:
  /// Fallback device selector
  cl::sycl::default_selector m_defaultSelector;

  /// Substring of the preferred device's name
  std::string m_deviceName;
};
}  // namespace Acts::Sycl
