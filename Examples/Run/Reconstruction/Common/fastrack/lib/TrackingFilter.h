#ifndef __TRACKING_FILTER__
#define __TRACKING_FILTER__

#include<cstring>
#include<vector>
#include<algorithm>

#include "AlgorithmParameters.h"
#include "Event.h"
#include "DataStructures.h"
#include "RecoEvent.h"

#define MAX_TS 2500
/*
typedef struct S_Clone {
public:

  struct CompareJ {
    bool operator()(const struct S_Clone* ts1, const struct S_Clone* ts2) {
      return ts1->m_J > ts2->m_J;
    }
  };
  

  void initialize(const SEGMENT* pS, int t) {
    m_vtype = t; 
    m_J = 0.0;

    m_vn.clear();
    
    m_vn.push_back(pS->m_n2);
    
    m_J += 1.0;//number of nodes
  
    //initialize using the segment parameters
  
    m_X[0] = pS->m_n2->m_phi;

    double dphi = pS->m_n1->m_phi - pS->m_n2->m_phi;

    if(dphi < -M_PI) dphi += 2*M_PI;
    if(dphi >  M_PI) dphi -= 2*M_PI;

    double dr = pS->m_n1->m_r - pS->m_n2->m_r;
    double dz = pS->m_n1->m_z - pS->m_n2->m_z;
    
    m_X[1] = dphi/dr;
    
    m_X[2] = 0.0;//omega in O-U model
    
    memset(&m_Cx[0][0], 0, sizeof(m_Cx));
    
    m_Cx[0][0] = 0.001;
    m_Cx[1][1] = 0.001;
    m_Cx[2][2] = 0.001;
    
    if(m_vtype != 0) {
      m_refCoord = pS->m_n2->m_z;
      m_Y[0] = pS->m_n2->m_r;
      m_Y[1] = dr/dz;
    }
    else {
      m_refCoord = pS->m_n2->m_r;
      m_Y[0] = pS->m_n2->m_z;
      m_Y[1] = dz/dr;
    }
    
    m_Cy[0][0] = 1.0;
    m_Cy[1][1] = 0.001;
    m_Cy[0][1] = m_Cy[1][0] = 0.0; 
  }

  void clone(const TS_Clone& ts) {

    m_vtype = ts.m_vtype;
    m_J = ts.m_J;
    m_refCoord = ts.m_refCoord; 
    m_chi2x = ts.m_chi2x;
    m_chi2y = ts.m_chi2y;
    
    memcpy(&m_Y[0], &ts.m_Y[0], sizeof(m_Y));
    memcpy(&m_Cy[0][0], &ts.m_Cy[0][0], sizeof(m_Cy));  
    memcpy(&m_X[0], &ts.m_X[0], sizeof(m_X));
    memcpy(&m_Cx[0][0], &ts.m_Cx[0][0], sizeof(m_Cx));

    m_vs.clear();
    m_vs.reserve(ts.m_vs.size());
    std::copy(ts.m_vs.begin(), ts.m_vs.end(), std::back_inserter(m_vs));
    m_vn.clear();
    m_vn.reserve(ts.m_vn.size());
    std::copy(ts.m_vn.begin(), ts.m_vn.end(), std::back_inserter(m_vn));
  }

  inline void addSegment(const SEGMENT* p) {
    m_vs.push_back(p);
  }

  inline void addNode(const NODE* p) {
    m_vn.push_back(p);
  }
  
  int m_vtype;
  double m_J;

  double m_refCoord;

  double m_chi2x, m_chi2y;
  
  double  m_X[3];
  double m_Cx[3][3];//r-phi plane
  
  double m_Y[2];
  double m_Cy[2][2];//r-z plane
  
  std::vector<const SEGMENT*> m_vs;
  std::vector<const NODE*> m_vn;
  
} S_CLONE;
*/

typedef struct TS_Clone {
public:

  struct CompareJ {
    bool operator()(const struct TS_Clone* ts1, const struct TS_Clone* ts2) {
      return ts1->m_J > ts2->m_J;
    }
  };
  

  void initialize(const SEGMENT* pS, int t) {
    m_vtype = t; 
    m_J = 0.0;

    m_vn.clear();
    
    m_vn.push_back(pS->m_n2);
    
    m_J += 1.0;//number of nodes
  
    //initialize using the segment parameters
  
    m_X[0] = pS->m_n2->m_phi;

    double dphi = pS->m_n1->m_phi - pS->m_n2->m_phi;

    if(dphi < -M_PI) dphi += 2*M_PI;
    if(dphi >  M_PI) dphi -= 2*M_PI;

    double dr = pS->m_n1->m_r - pS->m_n2->m_r;
    double dz = pS->m_n1->m_z - pS->m_n2->m_z;
    
    m_X[1] = dphi/dr;
    
    m_X[2] = 0.0;//omega in O-U model
    
    memset(&m_Cx[0][0], 0, sizeof(m_Cx));
    
    m_Cx[0][0] = 0.001;
    m_Cx[1][1] = 0.001;
    m_Cx[2][2] = 0.001;
    
    if(m_vtype != 0) {
      m_refCoord = pS->m_n2->m_z;
      m_Y[0] = pS->m_n2->m_r;
      m_Y[1] = dr/dz;
    }
    else {
      m_refCoord = pS->m_n2->m_r;
      m_Y[0] = pS->m_n2->m_z;
      m_Y[1] = dz/dr;
    }
    
    m_Cy[0][0] = 1.0;
    m_Cy[1][1] = 0.001;
    m_Cy[0][1] = m_Cy[1][0] = 0.0; 
  }

  void clone(const TS_Clone& ts) {

    m_vtype = ts.m_vtype;
    m_J = ts.m_J;
    m_refCoord = ts.m_refCoord; 
    m_chi2x = ts.m_chi2x;
    m_chi2y = ts.m_chi2y;
    
    memcpy(&m_Y[0], &ts.m_Y[0], sizeof(m_Y));
    memcpy(&m_Cy[0][0], &ts.m_Cy[0][0], sizeof(m_Cy));  
    memcpy(&m_X[0], &ts.m_X[0], sizeof(m_X));
    memcpy(&m_Cx[0][0], &ts.m_Cx[0][0], sizeof(m_Cx));

    m_vs.clear();
    m_vs.reserve(ts.m_vs.size());
    std::copy(ts.m_vs.begin(), ts.m_vs.end(), std::back_inserter(m_vs));
    m_vn.clear();
    m_vn.reserve(ts.m_vn.size());
    std::copy(ts.m_vn.begin(), ts.m_vn.end(), std::back_inserter(m_vn));
  }

  inline void addSegment(const SEGMENT* p) {
    m_vs.push_back(p);
  }

  inline void addNode(const NODE* p) {
    m_vn.push_back(p);
  }
  
  int m_vtype;
  double m_J;

  double m_refCoord;

  double m_chi2x, m_chi2y;
  
  double  m_X[3];
  double m_Cx[3][3];//r-phi plane
  
  double m_Y[2];
  double m_Cy[2][2];//r-z plane
  
  std::vector<const SEGMENT*> m_vs;
  std::vector<const NODE*> m_vn;
  
} TS_CLONE;


typedef class TrackingFilter {
 public:
  TrackingFilter(const ALGORITHM_PARAMETERS&, const HITS*, NEW_SEG_BANK* sb[N_SEG_BANKS]);
  ~TrackingFilter(){};

  float followTrack(const SEGMENT*, std::vector<const SEGMENT*>&, std::vector<const NODE*>&);

 protected:
  
  void propagate(const SEGMENT*, TS_CLONE&);

  void update(const NODE*, TS_CLONE&);

  int volumeType(const NODE*);
  
  std::vector<TS_CLONE*> m_tsVec;

  const ALGORITHM_PARAMETERS& m_params;
  
  const HITS* m_hits;
  NEW_SEG_BANK* m_segBank[N_SEG_BANKS];
  
  bool m_debug;

  TS_CLONE m_tClones[MAX_TS];

  int m_globalCloneCounter;

} TRACKING_FILTER;

/*
typedef class NewTrackingFilter {
 public:
  NewTrackingFilter(const ALGORITHM_PARAMETERS&, const HITS*, NEW_SEG_BANK* sb[N_SEG_BANKS]);
  ~NewTrackingFilter(){};

  float followTrack(const SEGMENT*, std::vector<const SEGMENT*>&, std::vector<const NODE*>&);

 protected:
  
  void propagate(const SEGMENT*, S_CLONE&);

  void update(const SEGMENT*, S_CLONE&);

  int volumeType(const NODE*);
  
  std::vector<S_CLONE*> m_sVec;

  const ALGORITHM_PARAMETERS& m_params;
  
  const HITS* m_hits;
  NEW_SEG_BANK* m_segBank[N_SEG_BANKS];
  
  bool m_debug;

  S_CLONE m_sClones[MAX_TS];

  int m_globalCloneCounter;

} NEW_TRACKING_FILTER;
*/

#endif
