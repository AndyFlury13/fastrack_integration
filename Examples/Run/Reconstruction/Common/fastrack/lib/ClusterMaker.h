#ifndef __CLUSTER_MAKER__
#define __CLUSTER_MAKER__

#include<vector>
#include<map>

#include "Geometry.h"
#include "DataStructures.h"
#include "Event.h"
#include "AlgorithmParameters.h"

typedef class ClusterMaker {
 public:
  ClusterMaker(const HITS*, const GEOMETRY&);
  ~ClusterMaker();

  void findClusters(const ALGORITHM_PARAMETERS&, HIT_LAYER*);

  const HITS* m_hits;
  const GEOMETRY& m_geo;

} CLUSTER_MAKER;


#endif
