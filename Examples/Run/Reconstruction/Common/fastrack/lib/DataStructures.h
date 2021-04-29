#ifndef __DATA_STRUCTURES_H__
#define __DATA_STRUCTURES_H__

#include<iostream>
#include<cmath>
#include<cstring>
#include<cstdint>

#include<vector>
#include<map>
#include<set>
#include<algorithm>

#include "Event.h"
#include "Geometry.h"

#define N_THREADS 3

#define N_SEGMENTS_PER_BANK 700000 //per bank

#define N_SEG_CONNS  6

#define N_SEG_BANKS 4
#define SEG_MASK_SHIFT  30
#define SEG_BANK_MASK   0xC0000000
#define SEG_INDEX_MASK  0x3FFFFFFF


#define MAX_OUT_SEGS 170

#define MAX_SEG_PER_NODE 160

#define N_TRACKS 100000

typedef struct PhiNode {
public:
PhiNode(float f, struct Node* p) : m_phi(f), m_pN(p) {};
PhiNode(const PhiNode& pn) : m_phi(pn.m_phi), m_pN(pn.m_pN) {};
  float m_phi;//extended phi-range beyond [-pi, +pi]
  struct Node* m_pN;
} PHI_NODE;

typedef struct Node {
public:

  struct ComparePhi {
    bool operator()(const struct Node* n1, const struct Node* n2) {
      return n1->m_phi < n2->m_phi;
    }
  };

  struct isEmpty {
    bool operator()(const struct Node* n) {
      return n->m_h == 0;
    }
  };

Node() : m_h(0), m_type(0), m_inOverlap(false) {
  m_minCutOnPar0 = -100.0;
  m_maxCutOnPar0 = 100.0;
}
  
Node(int hitIndex, const HITS* hits, int t) : m_h(hitIndex), m_type(t), m_inOverlap(false) {

  m_hits.push_back(hitIndex);
  m_x = hits->m_x[hitIndex];
  m_y = hits->m_y[hitIndex];
  m_z = hits->m_z[hitIndex];
  m_phi = hits->m_phi[hitIndex];
  m_r = hits->m_r[hitIndex];
  m_nCellsT = hits->m_clusterHeight[m_h];
  m_nCells = hits->m_clusterLength[m_h];
  m_minCutOnPar0 = -100.0;
  m_maxCutOnPar0 = 100.0;
}

  void initialize(int hitIndex, const HITS* hits, int t) {
    m_h = hitIndex;
    m_type = t;
    m_inOverlap = false;
    m_hits.push_back(hitIndex);
    m_x = hits->m_x[hitIndex];
    m_y = hits->m_y[hitIndex];
    m_z = hits->m_z[hitIndex];
    m_phi = hits->m_phi[hitIndex];
    m_r = hits->m_r[hitIndex];
    m_nCellsT = hits->m_clusterHeight[m_h];
    m_nCells = hits->m_clusterLength[m_h];
  }
  
  inline void addIn(int i) {
    if(m_in.size()<MAX_SEG_PER_NODE) {
      m_in.push_back(i);
    }
  }

  inline void addOut(int i) {
    if(m_out.size()<MAX_SEG_PER_NODE) {
      m_out.push_back(i);
    }
  }

  void calculateCellShape(const HITS* hits) {

    int nCellsT = hits->m_clusterHeight[m_h];
    if(m_hits.size()>1) {
      for(std::vector<int>::iterator hIt=m_hits.begin();hIt!=m_hits.end();++hIt) {
	int nC = hits->m_clusterHeight[*hIt];
	if(nC>nCellsT) nCellsT = nC;
      }
    }
    
    m_nCellsT = nCellsT;
    
    int nCells = hits->m_clusterLength[m_h];
    
    if(m_hits.size()>1) {
      for(std::vector<int>::iterator hIt=m_hits.begin();hIt!=m_hits.end();++hIt) {
	int nC = hits->m_clusterLength[*hIt];
	if(nC>nCells) nCells = nC;
      }
    }
    m_nCells = nCells;
  }
  
  void updateCoordinates(const HITS* hits) {
    
    m_x = 0.0;m_y = 0.0;m_z = 0.0;
    
    for(std::vector<int>::iterator it = m_hits.begin();it!=m_hits.end();++it) {
      m_x += hits->m_x[*it];
      m_y += hits->m_y[*it];
      m_z += hits->m_z[*it];
    }
    
    m_x /= m_hits.size();m_y /= m_hits.size();m_z /= m_hits.size();
    m_phi = atan2(m_y, m_x);
    m_r = sqrt(m_x*m_x + m_y*m_y);
  }

  int checkModset(const HITS* hits) {

    std::set<unsigned int> modSet;
    for(std::vector<int>::iterator it = m_hits.begin();it!=m_hits.end();++it) {
      modSet.insert(hits->m_mod_id[*it]);
    }
    
    return m_hits.size()-modSet.size();

  }

  int m_h;
  int m_type;

  float m_phi, m_r, m_x, m_y, m_z;

  std::vector<unsigned int> m_in;//indices of the segments in the segment bank
  std::vector<unsigned int> m_out;
  
  std::vector<int> m_hits;//possibly clusterized hits
  
  bool m_inOverlap;

  int m_nCells, m_nCellsT;
  float m_minCutOnPar0, m_maxCutOnPar0;

  ~Node() {};
  
} NODE;

typedef struct HitCluster {
public:
  
HitCluster(unsigned int i1, unsigned int i2) : m_idx1(i1), m_idx2(i2) {};
  
HitCluster(const HitCluster& c) : m_idx1(c.m_idx1), m_idx2(c.m_idx2) {};
  
  unsigned int m_idx1;
  unsigned int m_idx2;
  
} HIT_CLUSTER;

typedef struct HitBin {
public:

  HitBin() {};
  ~HitBin() {
    clear();
  }

  

  void sort() {
    std::sort(m_nodes.begin(), m_nodes.end(), NODE::ComparePhi());
  }

  void clear();
  void deleteUsedNodes(const HITS*);
  
  void selectByPhi(float, float, std::vector<NODE*>&);

  void searchPhiRange(float, float, int&, int&);

  inline void updatePhiRangeStart(float phi_min, int& j0) {
    while(m_l2phiIndices[j0].m_phi < phi_min && j0<(int)m_l2phiIndices.size()) j0++;
  }

  inline void updatePhiRangeEnd(float phi_max, int& j1) {
    while(m_l2phiIndices[j1].m_phi < phi_max && j1<(int)m_l2phiIndices.size()) j1++;
  }

  void removeEmptyNodes();

  void updateCoordinates(const HITS*);
  void calculateCellShape(const HITS*);

  int checkModset(const HITS*);

  void fillIndices(bool);

  std::vector<NODE*> m_nodes;

  std::vector<PHI_NODE> m_phiIndices;

  std::vector<NODE*> m_l1nodes;

  std::vector<PHI_NODE> m_l2phiIndices;
  
} HIT_BIN;


typedef struct HitLayer {
public:

  struct CompareAbsRefs {
    bool operator()(const struct HitLayer* pl1, const struct HitLayer* pl2) {
      return fabs(pl1->m_info.m_refCoord) > fabs(pl2->m_info.m_refCoord);
    }
  };

  HitLayer(const LAYER_SUMMARY& ls, float);
  ~HitLayer();

  void sort();
  void fillIndices();
  void clear();
  void deleteUsedNodes(const HITS*);
  void importNewHit(int, const HITS*);
  bool getBinRange(float, float, int&, int&);
  unsigned int numberOfHits();
  
  int getBinIndex(float);
  float getBinRadius(int);
  void getBinTaus(float, float, int, float&, float&);
  
  void prepareL1Nodes();
  
  const LAYER_SUMMARY& m_info;

  int m_nBins;
  float m_minEta, m_maxEta, m_deltaEta;
  HIT_BIN* m_bins;
  float* m_binCoords;
  float* m_binRads;
  float* m_binBoundMin;
  float* m_binBoundMax;

} HIT_LAYER;

typedef struct BinTable {
public:
  std::vector<int> m_binStart;
  std::vector<int> m_binEnd;
  std::vector<std::vector<float> > m_phiWidth;
} BIN_TABLE;


typedef class HitBank {
 public:

  HitBank(const GEOMETRY&, float);
  ~HitBank();
  
  void clear();

  void importHits(const HITS*, bool useMasked = false);

  void sort();

  void fillIndices();

  void getLayersFromVolume(int, std::vector<HIT_LAYER*>&);
  void getLayers(std::vector<HIT_LAYER*>&);
  
  HIT_LAYER* getLayer(unsigned int);
  
  unsigned int numberOfHits();

  void collectHits(const std::vector<int>&, const TRAJECTORY&, std::vector<int>&, const HITS*);

  void collectHitsInward(const TRAJECTORY&, std::vector<int>&);
  
  void findClosestHits(const std::vector<int>&, const TRAJECTORY&, std::vector<int>&, const HITS*);
  
  void getConnectingNodes(std::vector<NODE*>&);  

  void deleteUsedNodes(const HITS*);
 
 protected:

  void addNewLayer(unsigned int, unsigned int, float);
  
  void importNewHit(int);
    
  const GEOMETRY& m_geo;

  std::map<unsigned int, HIT_LAYER*> m_layMap;

  const HITS* m_hits;//SoA

  HIT_LAYER* m_layArray[128];
 
} HIT_BANK;

typedef struct Segment {
public:

  inline void initialize(NODE* n1, NODE* n2) {
    m_n1 = n1; 
    m_n2 = n2;
    m_level = 1;
    m_next = 1;
    m_nNei = 0;
  }


  NODE* m_n1;
  NODE* m_n2;
  
  int8_t m_level, m_next;

  uint8_t m_nNei;
  float m_p[4];
  
  unsigned int m_vNei[N_SEG_CONNS];//global indices of the connected segments

} SEGMENT;


typedef class NewSegmentBank {
 public:

 NewSegmentBank() : m_nSegments(0) {

  }

  ~NewSegmentBank() {

  }

  unsigned int numberOfSegments() {
    return m_nSegments;
  }
  
  void clear() {
    m_nSegments = 0;
  }
  
  SEGMENT m_S[N_SEGMENTS_PER_BANK];

  unsigned int m_nSegments;

} NEW_SEG_BANK;


#endif
