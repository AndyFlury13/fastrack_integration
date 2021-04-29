#ifndef __LAYER_LINKER_H__
#define __LAYER_LINKER_H__

#include<fstream>
#include<vector>

typedef struct LayerLink {
public:

  struct CompareFlow {
    bool operator()(const struct LayerLink& l1, const struct LayerLink& l2) {
      return l1.m_flow > l2.m_flow;
    }
  };

  struct CompareProb {
    bool operator()(const struct LayerLink& l1, const struct LayerLink& l2) {
      return l1.m_prob > l2.m_prob;
    }
  };
  
LayerLink(int idx, unsigned int d, unsigned int s, float p, float f) : m_index(idx), m_src(s), m_dst(d), m_prob(p), m_flow(f) {};
  ~LayerLink() {};

  int m_index;
  unsigned int m_src, m_dst;
  float m_prob, m_flow;
  
} LAYER_LINK;

typedef struct LinkStat {
public:
  struct CompareStats {
    bool operator()(const struct LinkStat& l1, const struct LinkStat& l2) {
      return l1.m_fSeg > l2.m_fSeg;
    }
  };
    
LinkStat(float f, int i) : m_fSeg(f), m_index(i) {};

  float m_fSeg;
  
  int m_index;
  
} LINK_STAT;


typedef class LayerLinker {
 public:
  LayerLinker(std::ifstream&);
  ~LayerLinker() {};

  //void getLongStrips(int, std::vector<int>&);

  std::vector<LAYER_LINK> m_links;

} LAYER_LINKER;

#endif
