// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/detail/CovarianceEngine.hpp"

#include "Acts/EventData/detail/TransformationBoundToFree.hpp"
#include "Acts/EventData/detail/TransformationFreeToBound.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

namespace Acts {
namespace {
/// Some type defs
using Jacobian = BoundMatrix;

using BoundState = std::tuple<BoundTrackParameters, Jacobian, double>;
using CurvilinearState =
    std::tuple<CurvilinearTrackParameters, Jacobian, double>;

/// @brief Evaluate the projection Jacobian from free to curvilinear parameters
///
/// @param [in] direction Normalised direction vector
///
/// @return Projection Jacobian
FreeToBoundMatrix freeToCurvilinearJacobian(const Vector3& direction) {
  // Optimized trigonometry on the propagation direction
  const double x = direction(0);  // == cos(phi) * sin(theta)
  const double y = direction(1);  // == sin(phi) * sin(theta)
  const double z = direction(2);  // == cos(theta)
  // can be turned into cosine/sine
  const double cosTheta = z;
  const double sinTheta = sqrt(x * x + y * y);
  const double invSinTheta = 1. / sinTheta;
  const double cosPhi = x * invSinTheta;
  const double sinPhi = y * invSinTheta;
  // prepare the jacobian to curvilinear
  FreeToBoundMatrix jacToCurv = FreeToBoundMatrix::Zero();
  if (std::abs(cosTheta) < s_curvilinearProjTolerance) {
    // We normally operate in curvilinear coordinates defined as follows
    jacToCurv(0, 0) = -sinPhi;
    jacToCurv(0, 1) = cosPhi;
    jacToCurv(1, 0) = -cosPhi * cosTheta;
    jacToCurv(1, 1) = -sinPhi * cosTheta;
    jacToCurv(1, 2) = sinTheta;
  } else {
    // Under grazing incidence to z, the above coordinate system definition
    // becomes numerically unstable, and we need to switch to another one
    const double c = sqrt(y * y + z * z);
    const double invC = 1. / c;
    jacToCurv(0, 1) = -z * invC;
    jacToCurv(0, 2) = y * invC;
    jacToCurv(1, 0) = c;
    jacToCurv(1, 1) = -x * y * invC;
    jacToCurv(1, 2) = -x * z * invC;
  }
  // Time parameter
  jacToCurv(5, 3) = 1.;
  // Directional and momentum parameters for curvilinear
  jacToCurv(2, 4) = -sinPhi * invSinTheta;
  jacToCurv(2, 5) = cosPhi * invSinTheta;
  jacToCurv(3, 4) = cosPhi * cosTheta;
  jacToCurv(3, 5) = sinPhi * cosTheta;
  jacToCurv(3, 6) = -sinTheta;
  jacToCurv(4, 7) = 1.;

  return jacToCurv;
}

/// @brief This function calculates the full jacobian from local parameters at
/// the start surface to bound parameters at the final surface
///
/// @note Modifications of the jacobian related to the
/// projection onto a surface is considered. Since a variation of the start
/// parameters within a given uncertainty would lead to a variation of the end
/// parameters, these need to be propagated onto the target surface. This an
/// approximated approach to treat the (assumed) small change.
///
/// @param [in] geoContext The geometry Context
/// @param [in] freeParameters Free, nominal parametrisation
/// @param [in] boundToFreeJacobian The projection jacobian from start local
/// to start free parameters
/// @param [in] freeTransportJacobian The transport jacobian from start free to
/// final free parameters
/// @param [in] freeToPathDerivatives Path length derivatives of the final free
/// parameters
/// @param [in, out] fullTransportJacobian The full jacobian from start local to
/// bound parameters at the final surface
/// @param [in] surface The final surface onto which the projection should be
/// performed
void boundToBoundJacobian(const GeometryContext& geoContext,
                          const FreeVector& freeParameters,
                          const BoundToFreeMatrix& boundToFreeJacobian,
                          const FreeMatrix& freeTransportJacobian,
                          const FreeVector& freeToPathDerivatives,
                          BoundMatrix& fullTransportJacobian,
                          const Surface& surface) {
  // Calculate the derivative of path length at the final surface or the
  // point-of-closest approach w.r.t. free parameters
  const FreeToPathMatrix freeToPath =
      surface.freeToPathDerivative(geoContext, freeParameters);
  // Calculate the jacobian from free to bound at the final surface
  FreeToBoundMatrix freeToBoundJacobian =
      surface.freeToBoundJacobian(geoContext, freeParameters);
  // Calculate the full jacobian from the local/bound parameters at the start
  // surface to local/bound parameters at the final surface
  // @note jac(locA->locB) = jac(gloB->locB)*(1+
  // pathCorrectionFactor(gloB))*jacTransport(gloA->gloB) *jac(locA->gloA)
  fullTransportJacobian =
      freeToBoundJacobian *
      (FreeMatrix::Identity() + freeToPathDerivatives * freeToPath) *
      freeTransportJacobian * boundToFreeJacobian;
}

/// @brief This function calculates the full jacobian from local parameters at
/// the start surface to final curvilinear parameters
///
/// @note Modifications of the jacobian related to the
/// projection onto a curvilinear surface is considered. Since a variation of
/// the start parameters within a given uncertainty would lead to a variation of
/// the end parameters, these need to be propagated onto the target surface.
/// This is an approximated approach to treat the (assumed) small change.
///
/// @param [in] direction Normalised direction vector
/// @param [in] boundToFreeJacobian The projection jacobian from local start
/// to global final parameters
/// @param [in] freeTransportJacobian The transport jacobian from start free to
/// final free parameters
/// @param [in] freeToPathDerivatives Path length derivatives of the final free
/// parameters
/// @param [in, out] jacFull The full jacobian from start local to curvilinear
/// parameters
///
/// @note The parameter @p surface is only required if projected to bound
/// parameters. In the case of curvilinear parameters the geometry and the
/// position is known and the calculation can be simplified
void boundToCurvilinearJacobian(const Vector3& direction,
                                const BoundToFreeMatrix& boundToFreeJacobian,
                                const FreeMatrix& freeTransportJacobian,
                                const FreeVector& freeToPathDerivatives,
                                BoundMatrix& fullTransportJacobian) {
  // Calculate the derivative of path length at the the curvilinear surface
  // w.r.t. free parameters
  FreeToPathMatrix freeToPath = FreeToPathMatrix::Zero();
  freeToPath.segment<3>(eFreePos0) = -1.0 * direction;
  // Calculate the jacobian from global to local at the curvilinear surface
  FreeToBoundMatrix freeToBoundJacobian = freeToCurvilinearJacobian(direction);
  // Calculate the full jocobian from the local parameters at the start surface
  // to curvilinear parameters
  // @note jac(locA->locB) = jac(gloB->locB)*(1+
  // pathCorrectionFactor(gloB))*jacTransport(gloA->gloB) *jac(locA->gloA)
  fullTransportJacobian =
      freeToBoundJacobian *
      (FreeMatrix::Identity() + freeToPathDerivatives * freeToPath) *
      freeTransportJacobian * boundToFreeJacobian;
}

/// @brief This function reinitialises the state members required for the
/// covariance transport
///
/// @param [in] geoContext The geometry context
/// @param [in, out] freeTransportJacobian The transport jacobian from start
/// free to final free parameters
/// @param [in, out] freeToPathDerivatives Path length derivatives of the free,
/// nominal parameters
/// @param [in, out] boundToFreeJacobian Projection jacobian of the last bound
/// parametrisation to free parameters
/// @param [in] freeParameters Free, nominal parametrisation
/// @param [in] surface The reference surface of the local parametrisation
Result<void> reinitializeJacobians(const GeometryContext& geoContext,
                                   FreeMatrix& freeTransportJacobian,
                                   FreeVector& freeToPathDerivatives,
                                   BoundToFreeMatrix& boundToFreeJacobian,
                                   const FreeVector& freeParameters,
                                   const Surface& surface) {
  using VectorHelpers::phi;
  using VectorHelpers::theta;

  // Reset the jacobians
  freeTransportJacobian = FreeMatrix::Identity();
  freeToPathDerivatives = FreeVector::Zero();

  // Get the local position
  const Vector3 position = freeParameters.segment<3>(eFreePos0);
  const Vector3 direction = freeParameters.segment<3>(eFreeDir0);
  auto lpResult = surface.globalToLocal(geoContext, position, direction);
  if (not lpResult.ok()) {
    ACTS_LOCAL_LOGGER(
        Acts::getDefaultLogger("CovarianceEngine", Logging::INFO));
    ACTS_FATAL(
        "Inconsistency in global to local transformation during propagation.")
  }
  // Transform from free to bound parameters
  Result<BoundVector> boundParameters = detail::transformFreeToBoundParameters(
      freeParameters, surface, geoContext);
  if (!boundParameters.ok()) {
    return boundParameters.error();
  }
  // Reset the jacobian from local to global
  boundToFreeJacobian =
      surface.boundToFreeJacobian(geoContext, *boundParameters);
  return Result<void>::success();
}

/// @brief This function reinitialises the state members required for the
/// covariance transport
///
/// @param [in, out] freeTransportJacobian The transport jacobian from start
/// free to final free parameters
/// @param [in, out] derivatives Path length derivatives of the free, nominal
/// parameters
/// @param [in, out] boundToFreeJacobian Projection jacobian of the last bound
/// parametrisation to free parameters
/// @param [in] direction Normalised direction vector
void reinitializeJacobians(FreeMatrix& freeTransportJacobian,
                           FreeVector& freeToPathDerivatives,
                           BoundToFreeMatrix& boundToFreeJacobian,
                           const Vector3& direction) {
  // Reset the jacobians
  freeTransportJacobian = FreeMatrix::Identity();
  freeToPathDerivatives = FreeVector::Zero();
  boundToFreeJacobian = BoundToFreeMatrix::Zero();

  // Optimized trigonometry on the propagation direction
  const double x = direction(0);  // == cos(phi) * sin(theta)
  const double y = direction(1);  // == sin(phi) * sin(theta)
  const double z = direction(2);  // == cos(theta)
  // can be turned into cosine/sine
  const double cosTheta = z;
  const double sinTheta = sqrt(x * x + y * y);
  const double invSinTheta = 1. / sinTheta;
  const double cosPhi = x * invSinTheta;
  const double sinPhi = y * invSinTheta;

  boundToFreeJacobian(eFreePos0, eBoundLoc0) = -sinPhi;
  boundToFreeJacobian(eFreePos0, eBoundLoc1) = -cosPhi * cosTheta;
  boundToFreeJacobian(eFreePos1, eBoundLoc0) = cosPhi;
  boundToFreeJacobian(eFreePos1, eBoundLoc1) = -sinPhi * cosTheta;
  boundToFreeJacobian(eFreePos2, eBoundLoc1) = sinTheta;
  boundToFreeJacobian(eFreeTime, eBoundTime) = 1;
  boundToFreeJacobian(eFreeDir0, eBoundPhi) = -sinTheta * sinPhi;
  boundToFreeJacobian(eFreeDir0, eBoundTheta) = cosTheta * cosPhi;
  boundToFreeJacobian(eFreeDir1, eBoundPhi) = sinTheta * cosPhi;
  boundToFreeJacobian(eFreeDir1, eBoundTheta) = cosTheta * sinPhi;
  boundToFreeJacobian(eFreeDir2, eBoundTheta) = -sinTheta;
  boundToFreeJacobian(eFreeQOverP, eBoundQOverP) = 1;
}
}  // namespace

namespace detail {

Result<BoundState> boundState(const GeometryContext& geoContext,
                              BoundSymMatrix& covarianceMatrix,
                              BoundMatrix& jacobian,
                              FreeMatrix& transportJacobian,
                              FreeVector& derivatives,
                              BoundToFreeMatrix& boundToFreeJacobian,
                              const FreeVector& parameters, bool covTransport,
                              double accumulatedPath, const Surface& surface) {
  // Covariance transport
  std::optional<BoundSymMatrix> cov = std::nullopt;
  if (covTransport) {
    // Initialize the jacobian from start local to final local
    jacobian = BoundMatrix::Identity();
    // Calculate the jacobian and transport the covarianceMatrix to final local.
    // Then reinitialize the transportJacobian, derivatives and the
    // boundToFreeJacobian
    transportCovarianceToBound(geoContext, covarianceMatrix, jacobian,
                               transportJacobian, derivatives,
                               boundToFreeJacobian, parameters, surface);
  }
  if (covarianceMatrix != BoundSymMatrix::Zero()) {
    cov = covarianceMatrix;
  }

  // Create the bound parameters
  Result<BoundVector> bv =
      detail::transformFreeToBoundParameters(parameters, surface, geoContext);
  if (!bv.ok()) {
    return bv.error();
  }
  // Create the bound state
  return std::make_tuple(
      BoundTrackParameters(surface.getSharedPtr(), *bv, std::move(cov)),
      jacobian, accumulatedPath);
}

CurvilinearState curvilinearState(BoundSymMatrix& covarianceMatrix,
                                  BoundMatrix& jacobian,
                                  FreeMatrix& transportJacobian,
                                  FreeVector& derivatives,
                                  BoundToFreeMatrix& boundToFreeJacobian,
                                  const FreeVector& parameters,
                                  bool covTransport, double accumulatedPath) {
  const Vector3& direction = parameters.segment<3>(eFreeDir0);

  // Covariance transport
  std::optional<BoundSymMatrix> cov = std::nullopt;
  if (covTransport) {
    // Initialize the jacobian from start local to final local
    jacobian = BoundMatrix::Identity();
    // Calculate the jacobian and transport the covarianceMatrix to final local.
    // Then reinitialize the transportJacobian, derivatives and the
    // boundToFreeJacobian
    transportCovarianceToCurvilinear(covarianceMatrix, jacobian,
                                     transportJacobian, derivatives,
                                     boundToFreeJacobian, direction);
  }
  if (covarianceMatrix != BoundSymMatrix::Zero()) {
    cov = covarianceMatrix;
  }

  // Create the curvilinear parameters
  Vector4 pos4 = Vector4::Zero();
  pos4[ePos0] = parameters[eFreePos0];
  pos4[ePos1] = parameters[eFreePos1];
  pos4[ePos2] = parameters[eFreePos2];
  pos4[eTime] = parameters[eFreeTime];
  CurvilinearTrackParameters curvilinearParams(
      pos4, direction, parameters[eFreeQOverP], std::move(cov));
  // Create the curvilinear state
  return std::make_tuple(std::move(curvilinearParams), jacobian,
                         accumulatedPath);
}

void transportCovarianceToBound(
    const GeometryContext& geoContext, BoundSymMatrix& boundCovariance,
    BoundMatrix& fullTransportJacobian, FreeMatrix& freeTransportJacobian,
    FreeVector& freeToPathDerivatives, BoundToFreeMatrix& boundToFreeJacobian,
    const FreeVector& freeParameters, const Surface& surface) {
  // Calculate the full jacobian from local parameters at the start surface to
  // current bound parameters
  boundToBoundJacobian(geoContext, freeParameters, boundToFreeJacobian,
                       freeTransportJacobian, freeToPathDerivatives,
                       fullTransportJacobian, surface);

  // Apply the actual covariance transport to get covariance of the current
  // bound parameters
  boundCovariance = fullTransportJacobian * boundCovariance *
                    fullTransportJacobian.transpose();

  // Reinitialize jacobian components:
  // ->The transportJacobian is reinitialized to Identity
  // ->The derivatives is reinitialized to Zero
  // ->The boundToFreeJacobian is initialized to that at the current surface
  reinitializeJacobians(geoContext, freeTransportJacobian,
                        freeToPathDerivatives, boundToFreeJacobian,
                        freeParameters, surface);
}

void transportCovarianceToCurvilinear(BoundSymMatrix& boundCovariance,
                                      BoundMatrix& fullTransportJacobian,
                                      FreeMatrix& freeTransportJacobian,
                                      FreeVector& freeToPathDerivatives,
                                      BoundToFreeMatrix& boundToFreeJacobian,
                                      const Vector3& direction) {
  // Calculate the full jacobian from local parameters at the start surface to
  // current curvilinear parameters
  boundToCurvilinearJacobian(direction, boundToFreeJacobian,
                             freeTransportJacobian, freeToPathDerivatives,
                             fullTransportJacobian);

  // Apply the actual covariance transport to get covariance of the current
  // curvilinear parameters
  boundCovariance = fullTransportJacobian * boundCovariance *
                    fullTransportJacobian.transpose();

  // Reinitialize jacobian components:
  // ->The free transportJacobian is reinitialized to Identity
  // ->The path derivatives is reinitialized to Zero
  // ->The boundToFreeJacobian is reinitialized to that at the current
  // curvilinear surface
  reinitializeJacobians(freeTransportJacobian, freeToPathDerivatives,
                        boundToFreeJacobian, direction);
}

}  // namespace detail
}  // namespace Acts
