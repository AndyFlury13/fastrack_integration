#ifndef __FAST_TRACK_FITTER_H__
#define __FAST_TRACK_FITTER_H__

#include<iostream>

#include "AlgorithmParameters.h"
#include "Geometry.h"
#include "Event.h"
#include "RecoEvent.h"
#include "DataStructures.h"
#include "FastTrackExtrapolator.h"

#include<vector>

#define MAX_CLONE 1000

typedef struct FitHit {
public:

FitHit(int h): m_h(h), m_status(0), m_dchi2(0.0) {};
  
  ~FitHit() {};

  int m_h;//hit index
  int m_status;//outlier, weight, etc
  float m_dchi2;
  float m_detr;
  
} FIT_HIT;


typedef class FastTrackFitter {
 public:
 FastTrackFitter(const ALGORITHM_PARAMETERS& p, const GEOMETRY& g, const HITS* h, const FTE& fte, int id = -1) : m_params(p), m_geo(g), m_hits(h), m_fte(fte), m_fitterId(id), m_stage(-1) {};
  ~FastTrackFitter() {};

  bool fitTrackCandidate(NEW_TRACK*, int) const;
  void fitTrack(NEW_TRACK*, int);

  void forwardPropagation(NEW_TRACK*, TRAJECTORY&);
  void backwardPropagation(NEW_TRACK*, TRAJECTORY&);
  bool getTrajectory(NEW_TRACK*, TRAJECTORY&);
  
  int addHits(NEW_TRACK*, std::vector<int>&);
  void addHitsBackward(NEW_TRACK*, std::vector<int>&);
 
  

 protected:

  int getTriplet(NEW_TRACK*, int*) const;
  int getLastTriplet(NEW_TRACK*, int*) const;

  bool initialize(NEW_TRACK*) const;
  
  void update(TRACK_STATE*, FIT_HIT*, float) const;

  void compareTrackStates(TRACK_STATE*, TRACK_STATE*) const;

  float calculateLikelihood(float, int, int*, HIT_INFO*) const;
  
  const ALGORITHM_PARAMETERS& m_params;
  
  const GEOMETRY& m_geo;
 
  const HITS* m_hits;

  const FTE& m_fte;//the extrapolator

  int m_fitterId;

  TRACK_CLONE m_tClones[MAX_CLONE];

  mutable int m_stage;
  
} FTF;

#endif
