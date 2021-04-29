// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/MaterialSlab.hpp"

#include <iosfwd>

namespace Acts {

/// @class HomogeneousSurfaceMaterial
///
/// It extends the ISurfaceMaterial virutal base class to describe
/// a simple homogeneous material on a surface
class HomogeneousSurfaceMaterial : public ISurfaceMaterial {
 public:
  /// Default Constructor - defaulted
  HomogeneousSurfaceMaterial() = default;

  /// Explicit constructor
  ///
  /// @param full are the full material properties
  /// @param splitFactor is the split for pre/post update
  HomogeneousSurfaceMaterial(const MaterialSlab& full, double splitFactor = 1.);

  /// Copy Constructor
  ///
  /// @param hsm is the source material
  HomogeneousSurfaceMaterial(const HomogeneousSurfaceMaterial& hsm) = default;

  /// Copy Move Constructor
  ///
  /// @param hsm is the source material
  HomogeneousSurfaceMaterial(HomogeneousSurfaceMaterial&& hsm) = default;

  /// Destructor
  ~HomogeneousSurfaceMaterial() override = default;

  /// Assignment operator
  ///
  /// @param hsm is the source material
  HomogeneousSurfaceMaterial& operator=(const HomogeneousSurfaceMaterial& hsm) =
      default;

  /// Assignment Move operator
  ///
  /// @param hsm is the source material
  HomogeneousSurfaceMaterial& operator=(HomogeneousSurfaceMaterial&& hsm) =
      default;

  /// Scale operator
  /// - it is effectively a thickness scaling
  ///
  /// @param scale is the scale factor
  HomogeneousSurfaceMaterial& operator*=(double scale) final;

  /// Equality operator
  ///
  /// @param hsm is the source material
  bool operator==(const HomogeneousSurfaceMaterial& hsm) const;

  /// @copydoc SurfaceMaterial::materialSlab(const Vector2&)
  ///
  /// @note the input parameter is ignored
  const MaterialSlab& materialSlab(const Vector2& lp) const final;

  /// @copydoc SurfaceMaterial::materialSlab(const Vector3&)
  ///
  /// @note the input parameter is ignored
  const MaterialSlab& materialSlab(const Vector3& gp) const final;

  /// @copydoc SurfaceMaterial::materialSlab(size_t, size_t)
  ///
  /// @param ib0 The bin at local 0 for retrieving the material
  /// @param ib1 The bin at local 1 for retrieving the material
  ///
  /// @note the input parameter is ignored
  const MaterialSlab& materialSlab(size_t ib0, size_t ib1) const final;

  /// The inherited methods - for MaterialSlab access
  using ISurfaceMaterial::materialSlab;

  /// The interited methods - for scale access
  using ISurfaceMaterial::factor;

  /// Output Method for std::ostream
  ///
  /// @param sl The outoput stream
  std::ostream& toStream(std::ostream& sl) const final;

 private:
  /// The five different MaterialSlab
  MaterialSlab m_fullMaterial = MaterialSlab();
};

inline const MaterialSlab& HomogeneousSurfaceMaterial::materialSlab(
    const Vector2& /*lp*/) const {
  return (m_fullMaterial);
}

inline const MaterialSlab& HomogeneousSurfaceMaterial::materialSlab(
    const Vector3& /*gp*/) const {
  return (m_fullMaterial);
}

inline const MaterialSlab& HomogeneousSurfaceMaterial::materialSlab(
    size_t /*ib0*/, size_t /*ib1*/) const {
  return (m_fullMaterial);
}

inline bool HomogeneousSurfaceMaterial::operator==(
    const HomogeneousSurfaceMaterial& hsm) const {
  return (m_fullMaterial == hsm.m_fullMaterial);
}

}  // namespace Acts
