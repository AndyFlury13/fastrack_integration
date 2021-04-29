#ifndef __RECO_EVENT_H__
#define __RECO_EVENT_H__

#include<iostream>
#include<cmath>
#include<cstring>
#include<vector>
#include<list>
#include "Event.h"
#include "Geometry.h"
#include "DataStructures.h"

#define MAX_HIT_ON_TRACK 19

typedef struct HitInfo {
public:

  HitInfo() {};

HitInfo(float dc, int stat) : m_dchi2(dc), m_status(stat) {};
  
HitInfo(const HitInfo& hi) : m_dchi2(hi.m_dchi2), m_status(hi.m_status) {};

  ~HitInfo() {};

  float m_dchi2;
  int m_status;

} HIT_INFO;


typedef struct TrackState {
public:
TrackState() : m_pDE(0), m_chi2(0.0), m_ndof(-5), m_initialized(false) {};
  
  TrackState(float*, const DETECTOR_ELEMENT*);
TrackState(const TrackState& ts) : m_pDE(ts.m_pDE), m_chi2(ts.m_chi2), m_ndof(ts.m_ndof), m_initialized(true) {

  memcpy(&m_Rk[0],&ts.m_Rk[0],sizeof(m_Rk));
  memcpy(&m_Re[0],&ts.m_Re[0],sizeof(m_Re));

  memcpy(&m_Ck[0],&ts.m_Ck[0],sizeof(m_Ck));
  memcpy(&m_Ce[0],&ts.m_Ce[0],sizeof(m_Ce));
  
}
  
  ~TrackState() {};

  void initialize(float*, const DETECTOR_ELEMENT*);
  void initialize(const TrackState& ts) {
    m_pDE = ts.m_pDE;
    m_chi2 = ts.m_chi2;
    m_ndof = ts.m_ndof;
    m_initialized = true;

    memcpy(&m_Rk[0],&ts.m_Rk[0],sizeof(m_Rk));
    memcpy(&m_Re[0],&ts.m_Re[0],sizeof(m_Re));

    memcpy(&m_Ck[0],&ts.m_Ck[0],sizeof(m_Ck));
    memcpy(&m_Ce[0],&ts.m_Ce[0],sizeof(m_Ce));
  }
  
  void clone(const TrackState&);

  void resetCovariance();

  void report() {
    std::cout<<"xLoc="<<m_Rk[0]<<" yLoc="<<m_Rk[1]<<" phi="
	     <<m_Rk[2]<<" theta="<<m_Rk[3]<<" mom="<<1/m_Rk[4]<<" chi2="
	     <<m_chi2<<" ndof="<<m_ndof<<std::endl;
    std::cout<<"------- Covariance: ----------"<<std::endl;

    int idx=0;
    
    for(int i=0;i<5;i++) {
      for(int j=0;j<=i;j++,idx++) {
	std::cout<<m_Ce[idx]<<" ";
      }
      std::cout<<std::endl;
    }
    
  }

  bool valid();

  float m_Rk[5];
  const DETECTOR_ELEMENT* m_pDE;

  float m_Re[5];
  float m_chi2;
  int m_ndof;
  float m_Ck[15], m_Ce[15];//for the fast implementation of the extrapolator
  float m_initialized;
  
} TRACK_STATE;


typedef struct NewTrack {
struct CompareLength {
    bool operator()(const struct NewTrack* pT1, const struct NewTrack* pT2) {
      return pT1->m_nlayers > pT2->m_nlayers;
    }
  };

  struct CompareLikelihood {
    bool operator()(const struct NewTrack* pT1, const struct NewTrack* pT2) {
      return pT1->m_ll > pT2->m_ll;
    }
  };
  
  struct CompareHits {
    bool operator()(const struct NewTrack* pT1, const struct NewTrack* pT2) {
      return pT1->m_nHits > pT2->m_nHits;
    }
  };
  
  struct CompareQuality {
    bool operator()(const struct NewTrack* pT1, const struct NewTrack* pT2) {
      return pT1->m_Q > pT2->m_Q;
    }
  };

  
  struct SortHitsByDistance {
  SortHitsByDistance(const HITS* p) : m_p(p) {};
    bool operator()(int h1, int h2) {
      return m_p->m_distance[h1] < m_p->m_distance[h2];
    }
    const HITS* m_p;
  };
  
  
NewTrack() : m_nHits(0), m_ll(-1e8), m_status(0) {};

  void initialize(int, const int*, int, float);
  
  ~NewTrack() {

  }

  void print(const HITS*);
  void update(const HITS*);
  
  void sortHits(const HITS*);
  void maskHits(int, const HITS*) const;
  void add(std::vector<int>&, const HITS*);
  void getLayers2(std::vector<int>&, const HITS*) const;
  void getLayers(std::vector<int>&, const HITS*) const;
  void getModules(std::vector<int>&, const HITS*) const;
  
  void removeDuplicates();

  bool valid() {
    return (m_endTrackState.valid() && m_beginTrackState.valid());
  }
  
  int m_trackId;

  int m_nHits;

  int m_hits[MAX_HIT_ON_TRACK];//hit indices

  unsigned int m_nlayers;

  int m_cloneFlag;

  float m_Q;

  float m_ll;

  float m_chi2;
  int m_ndof;

  int m_status;
  
  TRACK_STATE m_endTrackState;
  TRACK_STATE m_beginTrackState;

} NEW_TRACK;

/*
typedef struct Track {
public:

  struct CompareLength {
    bool operator()(const struct Track* pT1, const struct Track* pT2) {
      return pT1->m_nlayers > pT2->m_nlayers;
    }
  };

  struct CompareLikelihood {
    bool operator()(const struct Track* pT1, const struct Track* pT2) {
      return pT1->m_ll > pT2->m_ll;
    }
  };
  
  struct CompareHits {
    bool operator()(const struct Track* pT1, const struct Track* pT2) {
      return pT1->m_nHits > pT2->m_nHits;
    }
  };
  
  struct CompareQuality {
    bool operator()(const struct Track* pT1, const struct Track* pT2) {
      return pT1->m_Q > pT2->m_Q;
    }
  };

  
  struct SortHitsByDistance {
  SortHitsByDistance(const HITS* p) : m_p(p) {};
    bool operator()(int h1, int h2) {
      return m_p->m_distance[h1] < m_p->m_distance[h2];
    }
    const HITS* m_p;
  };
  
  
  Track(int, const int*, int, float, int);
  ~Track();
  
  void print(const HITS*);
  void update(const HITS*);
  
  void sortHits(const HITS*);
  void maskHits(int, const HITS*);
  void add(std::vector<int>&, const HITS*);
  float calculateLikelihood(float, const HITS*);
  void removeDuplicates();
  int m_trackId;

  int m_nHits;

  int m_hits[MAX_HIT_ON_TRACK];//hit indices

  std::vector<HIT_INFO> m_info; //clone marker, variance, isOutlier, dchi2, etc.

  unsigned  int m_nlayers;

  int m_cloneFlag;

  float m_Q;

  float m_ll;

  float m_chi2;
  int m_ndof;

  int m_status;

  int m_stage;
  
  TRACK_STATE* m_endTrackState;
  TRACK_STATE* m_beginTrackState;

} TRACK;
*/

typedef struct TrackClone {

  struct CompareQuality {
    bool operator()(const struct TrackClone* pT1, const struct TrackClone* pT2) {
      return pT1->m_Q > pT2->m_Q;
    }
  };
  
TrackClone() : m_Q(0.0), m_nHits(0) {};

TrackClone(TRACK_STATE* pTS) : m_Q(0.0), m_nHits(0) {
  m_ts.clone(*pTS);
}
  
TrackClone(const TrackClone& tc) : m_Q(tc.m_Q), m_nHits(tc.m_nHits) {
  memcpy(m_hits, tc.m_hits, m_nHits*sizeof(int));
  m_ts.clone(tc.m_ts);
}
  ~TrackClone() {};

  void initialize(TRACK_STATE* pTS) {
    m_Q = 0.0;
    m_nHits = 0;
    m_ts.clone(*pTS);
  }

  void clone(const TrackClone& tc) {
    m_Q = tc.m_Q;
    m_nHits = tc.m_nHits;
    memcpy(&m_hits[0], &tc.m_hits[0], m_nHits*sizeof(int));
    m_ts.clone(tc.m_ts);
  }

  void addHit(int h, float dQ) {
    if(m_nHits < MAX_HIT_ON_TRACK) {
      m_hits[m_nHits++] = h;
      m_Q += dQ;
    }
  }
  
  void addHole(float dQ) {
    m_Q += dQ;
  }

  float m_Q;
  int m_nHits;
  int m_hits[MAX_HIT_ON_TRACK];
  
  TRACK_STATE m_ts;
} TRACK_CLONE;


typedef struct TrackGroup {
public:

  struct CompareSize {
    bool operator()(const struct TrackGroup* g1, const struct TrackGroup* g2) {
      return g1->m_size < g2->m_size;
    }
  };

TrackGroup(int gid, const std::vector<NEW_TRACK*>& v) : m_gid(gid), m_vTracks(v), m_size(v.size()) {};
  ~TrackGroup() {};


  int m_gid;

  std::vector<NEW_TRACK*> m_vTracks;
  
  unsigned int m_size;

} TRACK_GROUP;


#endif

