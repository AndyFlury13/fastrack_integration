// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Surfaces/Surface.hpp"

#include <sstream>

namespace Acts {

/// Simple struct to select surfaces
struct SurfaceSelector {
  bool selectSensitive = true;
  bool selectMaterial = false;
  bool selectPassive = false;

  /// SurfaceSelector with options
  ///
  /// @param sSensitive is the directive to select sensitive surfaces
  /// @param sMaterial is the directive to select material surfaces
  /// @param sPassive is the directive to select passive surfaces
  SurfaceSelector(bool sSensitive = true, bool sMaterial = false,
                  bool sPassive = false)
      : selectSensitive(sSensitive),
        selectMaterial(sMaterial),
        selectPassive(sPassive) {}

  /// Call operator to check if a surface should be selected
  ///
  /// @param surface is the test surface
  bool operator()(const Acts::Surface& surface) const {
    if (selectSensitive && surface.associatedDetectorElement() != nullptr) {
      return true;
    }
    if (selectMaterial && surface.surfaceMaterial() != nullptr) {
      return true;
    }
    if (selectPassive) {
      return true;
    }
    return false;
  }
};

/// The information to be writtern out per hit surface
struct SurfaceHit {
  const Surface* surface = nullptr;
  Vector3 position;
  Vector3 direction;
};

/// A Surface Collector struct
/// templated with a Selector type
///
/// Whenever a surface is passed in the propagation
/// that satisfies the selector, it is recorded
/// for further usage in the flow.
template <typename Selector = SurfaceSelector>
struct SurfaceCollector {
  /// The selector used for this surface
  Selector selector;

  /// Simple result struct to be returned
  /// It has all the SurfaceHit objects that
  /// are collected (and thus have been selected)
  struct this_result {
    std::vector<SurfaceHit> collected;
  };

  using result_type = this_result;

  /// Collector action for the ActionList of the Propagator
  /// It checks if the propagator state has a current surface,
  /// in which case the action is performed:
  /// - it records the surface given the configuration
  ///
  /// @tparam propagator_state_t is the type of Propagator state
  /// @tparam stepper_t Type of the stepper used for the propagation
  ///
  /// @param [in,out] state is the mutable stepper state object
  /// @param [in] stepper The stepper in use
  /// @param [in,out] result is the mutable result object
  template <typename propagator_state_t, typename stepper_t>
  void operator()(propagator_state_t& state, const stepper_t& stepper,
                  result_type& result) const {
    const auto& logger = state.options.logger;
    // The current surface has been assigned by the navigator
    if (state.navigation.currentSurface &&
        selector(*state.navigation.currentSurface)) {
      // Create for recording
      SurfaceHit surface_hit;
      surface_hit.surface = state.navigation.currentSurface;
      surface_hit.position = stepper.position(state.stepping);
      surface_hit.direction = stepper.direction(state.stepping);
      // Save if in the result
      result.collected.push_back(surface_hit);
      // Screen output
      ACTS_VERBOSE("Collect surface  "
                   << state.navigation.currentSurface->geometryId());
    }
  }

  /// Pure observer interface
  /// - this does not apply to the surface collector
  template <typename propagator_state_t, typename stepper_t>
  void operator()(propagator_state_t& /*state*/,
                  const stepper_t& /*unused*/) const {}
};

}  // namespace Acts
