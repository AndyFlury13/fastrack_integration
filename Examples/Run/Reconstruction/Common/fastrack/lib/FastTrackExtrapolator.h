#ifndef __FAST_TRACK_EXTRAPOLATOR_H__
#define __FAST_TRACK_EXTRAPOLATOR_H__

#include "AlgorithmParameters.h"
#include "Geometry.h"
#include "FieldModel.h"
#include "DataStructures.h"
#include "RecoEvent.h"

typedef class FastTrackExtrapolator {
 public:

 FastTrackExtrapolator(const ALGORITHM_PARAMETERS& p, const GEOMETRY& g, const HITS* h) : m_params(p), m_geo(g), m_hits(h) {
    m_fModel.initialize(p);
  }

  ~FastTrackExtrapolator() {};

  int fastExtrapolate(TRACK_STATE*, int, float, TRAJECTORY&) const;
  void extrapolateOutsideFast(const TRACK_STATE*, TRAJECTORY&) const;
  void extrapolateInsideFast(const TRACK_STATE*,  TRAJECTORY&) const;
  int fastExtrapolator(TRACK_STATE*, int, bool checkBounds = false) const;
  

 protected:

  inline void crossProduct(float const *, float const *, float*) const;
  inline void fastCrossProduct(float const *, float const *, float*) const;
  inline void fastCrossProduct3(float const *, float const *, float const *, float const *,
				float*, float*, float*) const;

  int RungeKutta3(float*, float*, const DETECTOR_ELEMENT*, float minStep = 0.0) const;

  FIELD_MODEL m_fModel;
  const ALGORITHM_PARAMETERS& m_params;
  const GEOMETRY& m_geo;

  const HITS* m_hits;

} FTE;


#endif
