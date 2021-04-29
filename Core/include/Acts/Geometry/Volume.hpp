// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryObject.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/BoundingBox.hpp"

#include <iosfwd>
#include <memory>

namespace Acts {

class VolumeBounds;
using VolumeBoundsPtr = std::shared_ptr<const VolumeBounds>;

/// @class Volume
///
/// It inhertis of GeometryObject for TDD identification
///
/// Base class for all volumes inside the tracking realm, it defines
/// the interface for inherited Volume classes
/// regarding the geometrical information.

class Volume : public virtual GeometryObject {
 public:
  using BoundingBox = AxisAlignedBoundingBox<Volume, ActsScalar, 3>;

  /// Explicit constructor with shared arguments
  ///
  /// @param transform is the transform to position the volume in 3D space
  /// @param volbounds is the volume boundary definitions
  /// @note This will automatically build an oriented bounding box with an
  /// envelope value of (0.05, 0.05, 0.05)mm
  Volume(const Transform3& transform, VolumeBoundsPtr volbounds);

  /// Copy Constructor - with optional shift
  ///
  /// @param vol is the source volume for the copy
  /// @param shift is the optional shift applied as : shift * vol.transform()
  /// @note This will automatically build an oriented bounding box with an
  /// envelope value of (0.05, 0.05, 0.05)mm
  Volume(const Volume& vol, const Transform3& shift = Transform3::Identity());

  Volume() = delete;
  virtual ~Volume() = default;

  /// Assignment operator
  ///
  /// @param vol is the source volume to be copied
  Volume& operator=(const Volume& vol);

  /// Return methods for geometry transform
  const Transform3& transform() const;

  /// Returns the inverted transform of this volume.
  const Transform3& itransform() const;

  /// returns the center of the volume
  const Vector3& center() const;

  /// returns the volumeBounds()
  const VolumeBounds& volumeBounds() const;

  /// Construct bounding box for this shape
  /// @param envelope Optional envelope to add / subtract from min/max
  /// @return Constructed bounding box pointing to this volume
  BoundingBox boundingBox(const Vector3& envelope = {0, 0, 0}) const;

  /// Construct oriented bounding box for this shape
  /// @return Constructed oriented bounding box pointing to this volume
  const BoundingBox& orientedBoundingBox() const;

  /// Inside() method for checks
  ///
  /// @param gpos is the position to be checked
  /// @param tol is the tolerance parameter
  ///
  /// @return boolean indicator if the position is inside
  bool inside(const Vector3& gpos, double tol = 0.) const;

  /// The binning position method
  /// - as default the center is given, but may be overloaded
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param bValue is the binning value schema
  ///
  /// @return vector 3D that can be used for the binning
  Vector3 binningPosition(const GeometryContext& gctx,
                          BinningValue bValue) const override;

 protected:
  Transform3 m_transform;
  Transform3 m_itransform;
  Vector3 m_center;
  VolumeBoundsPtr m_volumeBounds;
  BoundingBox m_orientedBoundingBox;
};

inline const Transform3& Volume::transform() const {
  return m_transform;
}

inline const Transform3& Volume::itransform() const {
  return m_itransform;
}

inline const Vector3& Volume::center() const {
  return m_center;
}

inline const VolumeBounds& Volume::volumeBounds() const {
  return (*(m_volumeBounds.get()));
}

/**Overload of << operator for std::ostream for debug output*/
std::ostream& operator<<(std::ostream& sl, const Volume& vol);

}  // namespace Acts
