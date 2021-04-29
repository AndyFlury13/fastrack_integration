#ifndef __EVENT_H__
#define __EVENT_H__

typedef struct Hits {
public:
Hits() : m_phi(0), m_r(0), m_mask(0), m_m1(0), m_m2(0), m_eta(0), m_distance(0),
    m_clusterLength(0), m_clusterHeight(0) {};
  ~Hits() {
    delete[] m_phi;
    delete[] m_r;
    delete[] m_mask;
    delete[] m_m1;
    delete[] m_m2;
    delete[] m_eta;
    delete[] m_distance;
    delete[] m_clusterLength;
    delete[] m_clusterHeight;
    
  }
  int m_nHits;
  int m_nCells;
  
  const int* m_hit_index;
  const int* m_vol_id;
  const int* m_lay_id;
  const int* m_mod_id;
  const float* m_x;
  const float* m_y;
  const float* m_z;
  float* m_phi;
  float* m_r;
  int*   m_mask;
  float* m_m1;
  float* m_m2;
  float* m_eta;
  float* m_distance;

  char*  m_clusterLength;
  char*  m_clusterHeight;

  const int* m_true_hit_id;
  const float* m_true_x;
  const float* m_true_y;
  const float* m_true_z;
  const float* m_true_px;
  const float* m_true_py;
  const float* m_true_pz;
  const unsigned long* m_particle;
  const float* m_weight;

  int m_eventId;
  
} HITS;

typedef struct Particle {
public:
  Particle() {};
  ~Particle() {};

  int m_particleType;
  float m_vx, m_vy, m_vz;
  float m_px, m_py, m_pz;
  int m_q;
  int m_nhits;

} PARTICLE;
#endif
