#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__

#include<iostream>
#include<fstream>
#include<cstring>
#include<map>
#include<vector>

#define MAX_POINT_ON_TRACK 150

typedef struct Point3d {
public:
Point3d(float x, float y, float z, unsigned int k = 0) : m_x(x), m_y(y), m_z(z), m_key(k) {};
Point3d(const Point3d& p) : m_x(p.m_x), m_y(p.m_y), m_z(p.m_z), m_key(p.m_key) {};
  Point3d() {};
  
  ~Point3d() {};
  
  float m_x, m_y, m_z;
  unsigned int m_key;
  
} POINT_3D;

typedef struct Trajectory {

public:

Trajectory() : m_nPoints(0), m_nPointsOnTrack(0) {};

  inline void addPoint(float x, float y, float z) {
    if(m_nPoints >= MAX_POINT_ON_TRACK) return;
    m_p[m_nPoints].m_x = x;
    m_p[m_nPoints].m_y = y;
    m_p[m_nPoints].m_z = z;
    m_nPoints++;
  }

  inline void addPoint(const POINT_3D& p) {
    if(m_nPoints >= MAX_POINT_ON_TRACK) return;
    m_p[m_nPoints].m_x = p.m_x;
    m_p[m_nPoints].m_y = p.m_y;
    m_p[m_nPoints].m_z = p.m_z;
    m_nPoints++;
  }

  inline void addPoint(float x, float y, float z, unsigned int k) {
    if(m_nPoints >= MAX_POINT_ON_TRACK) return;
    m_p[m_nPoints].m_x = x;
    m_p[m_nPoints].m_y = y;
    m_p[m_nPoints].m_z = z;
    m_p[m_nPoints].m_key = k;
    m_nPoints++;
  }
  
  void clear() {
    m_nPoints = 0;
  }

  int m_nPoints, m_nPointsOnTrack;
  POINT_3D m_p[MAX_POINT_ON_TRACK];

private:

Trajectory(const Trajectory&t) : m_nPoints(t.m_nPoints) {};
  
} TRAJECTORY;


typedef struct DetectorElement {
public:
DetectorElement(unsigned int v, unsigned int l, unsigned int m) : m_volume_id(v), m_layer_id(l), m_module_id(m) {};
  ~DetectorElement(){};

  unsigned int m_volume_id , m_layer_id, m_module_id;
  float m_cx, m_cy, m_cz;
  float m_rot_xu, m_rot_xv, m_rot_xw;
  float m_rot_yu, m_rot_yv, m_rot_yw;
  float m_rot_zu, m_rot_zv, m_rot_zw;
  float m_module_t, m_module_minhu, m_module_maxhu, m_module_hv, m_pitch_u, m_pitch_v;

  float m_LG[3][3];
  float m_GL[3][3];

  float m_P[4];
  
  bool m_isRect;

  float m_R;

  void prepare();
  void toLocal(const float*, float*) const;
  void toGlobal(const float*, float*) const;

  inline void rotateToLocal(const float* X, float* Y) const {
    Y[0] = m_GL[0][0]*X[0] + m_GL[0][1]*X[1] + m_GL[0][2]*X[2]; 
    Y[1] = m_GL[1][0]*X[0] + m_GL[1][1]*X[1] + m_GL[1][2]*X[2];
    Y[2] = m_GL[2][0]*X[0] + m_GL[2][1]*X[1] + m_GL[2][2]*X[2];
  }

  inline float localVz(const float* X) const {
    return (m_GL[2][0]*X[0] + m_GL[2][1]*X[1] + m_GL[2][2]*X[2]);
  }

  void rotateToGlobal(const float*, float*) const;
  
  inline void getPlaneParameters(float* P) const {
    memcpy(P, &m_P[0], sizeof(m_P));
  }
  
  void getTransformLG(float m[3][3]) const;
  void getTransformGL(float g[3][3]) const;
  bool inOverlap(const float*) const;
  bool isInSide(const float*, float tol = 0.0) const;
  void getResidualTransform(float, float, float*) const;
  void print() const;

  void transformTrackStateToLocal(float*, const float*, const float*, float*) const;
  
} DETECTOR_ELEMENT;

typedef struct LayerSummary {
public:
  unsigned int m_vol_id, m_lay_id;
  int m_type;
  float m_refCoord;
  float m_minCoord;
  float m_maxCoord;

  void print() {
    std::cout<<"V"<<m_vol_id<<" L"<<m_lay_id<<" refCoord="<<m_refCoord<<" min "<<m_minCoord<<" max="<<m_maxCoord<<std::endl;
  }

} LAYER_SUMMARY;

typedef struct Layer {
public:

  struct Compare {
    bool operator()(const struct Layer* l1, const struct Layer* l2) {
      return l1->m_summary.m_refCoord < l2->m_summary.m_refCoord;
    }

  };
  
  Layer(int t, unsigned int v, unsigned int l) {
    m_summary.m_type = t;
    m_summary.m_vol_id = v;
    m_summary.m_lay_id = l;
  }
  
  ~Layer() {};

  void addNewDetectorElement(const DETECTOR_ELEMENT* p) {
    m_deVec.push_back(p);
  }

  void finalize();
  
  std::vector<const DETECTOR_ELEMENT*> m_deVec;
  LAYER_SUMMARY m_summary;
  
} LAYER;

typedef struct Volume {
public:
Volume(int t);
  ~Volume();

  void addNewDetectorElement(const DETECTOR_ELEMENT*);
  void finalize();
  void getLayerIds(std::vector<unsigned int>&) const;
  void getLayerSummary(std::vector<LAYER_SUMMARY>&) const;
  const LAYER* getLayer(unsigned int) const;
  int getLayerIntersection(float, float, float, float, const POINT_3D&, const POINT_3D&, TRAJECTORY&) const;

  void initializeBounds(float minR, float maxR, float minZ, float maxZ) {
    m_minR = minR;
    m_maxR = maxR;
    m_minZ = minZ;
    m_maxZ = maxZ;
  }
  
  int m_type;
  
  std::map<unsigned int, LAYER*> m_layMap;

  std::vector<const LAYER*> m_layVec;

  float m_minR, m_maxR, m_minZ, m_maxZ;

} VOLUME;

typedef class Geometry {
 public:
  Geometry(std::ifstream&);
  ~Geometry();
  
  unsigned long calculateKey(unsigned int, unsigned int, unsigned int) const;
  const DETECTOR_ELEMENT* getDetectorElement(unsigned int, unsigned int, unsigned int) const;
  
  void getVolumeIds(std::vector<unsigned int>&) const;

  const VOLUME* getVolume(unsigned int) const;

  bool getLayerIntersection(const POINT_3D&, const POINT_3D&, TRAJECTORY&) const;
 
 protected:
  std::map<unsigned long, const DETECTOR_ELEMENT*> m_deMap;
  bool addNewDetectorElement(const DETECTOR_ELEMENT*);
  void finalize();//calculate summaries, etc.
  
  int getVolumeId(float, float) const;


  std::map<unsigned int, VOLUME*> m_volMap;
  
  DETECTOR_ELEMENT* m_deArray[84][4000];
  
} GEOMETRY;


#endif
