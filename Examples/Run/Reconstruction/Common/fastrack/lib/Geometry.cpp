#include<iostream>
#include<cstring>
#include<cmath>

#include "Geometry.h"

#include<algorithm>

void DetectorElement::prepare() {

  m_LG[0][0] = m_rot_xu;m_LG[0][1] = m_rot_xv;m_LG[0][2] = m_rot_xw;
  m_LG[1][0] = m_rot_yu;m_LG[1][1] = m_rot_yv;m_LG[1][2] = m_rot_yw;
  m_LG[2][0] = m_rot_zu;m_LG[2][1] = m_rot_zv;m_LG[2][2] = m_rot_zw;

  m_GL[0][0] =  m_LG[1][1]*m_LG[2][2] - m_LG[1][2]*m_LG[2][1];
  m_GL[1][0] = -m_LG[1][0]*m_LG[2][2] + m_LG[1][2]*m_LG[2][0];
  m_GL[2][0] =  m_LG[1][0]*m_LG[2][1] - m_LG[1][1]*m_LG[2][0];
  m_GL[0][1] = -m_LG[0][1]*m_LG[2][2] + m_LG[0][2]*m_LG[2][1];
  m_GL[1][1] =  m_LG[0][0]*m_LG[2][2] - m_LG[0][2]*m_LG[2][0];
  m_GL[2][1] = -m_LG[0][0]*m_LG[2][1] + m_LG[0][1]*m_LG[2][0];
  m_GL[0][2] =  m_LG[0][1]*m_LG[1][2] - m_LG[0][2]*m_LG[1][1];
  m_GL[1][2] = -m_LG[0][0]*m_LG[1][2] + m_LG[0][2]*m_LG[1][0];
  m_GL[2][2] =  m_LG[0][0]*m_LG[1][1] - m_LG[0][1]*m_LG[1][0];

  m_isRect = fabs(m_module_maxhu - m_module_minhu) < 0.001;

  m_R = 1e8;

  if(!m_isRect) {
    float eps = m_module_maxhu/m_module_minhu;
    m_R = m_module_hv*(1+eps)/(eps-1);
  }

  //normal

  float n[3] = {m_LG[0][2], m_LG[1][2], m_LG[2][2]};
  
  for(int i=0;i<3;i++) m_P[i] = n[i];

  m_P[3] = -m_cx*n[0] - m_cy*n[1] - m_cz*n[2];
  
  //test
  /*
  float T[3][3];
  memset(&T[0][0], 0, sizeof(T));
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++) {
      for(int k=0;k<3;k++) T[i][j] += m_LG[i][k]*m_GL[k][j];
    }
  for(int i=0;i<3;i++) std::cout<<T[i][0]<<" "<<T[i][1]<<" "<<T[i][2]<<" "<<std::endl;
  */
  
}

void DetectorElement::getResidualTransform(float mu, float mv, float* sc) const {

  if(m_isRect) {
    sc[0] = 0.0;sc[1] = 1.0;
    return;
  }

  float mR = mv + m_R;
  float S = 1/sqrt(mu*mu+mR*mR);

  sc[0] = mu*S;sc[1] = mR*S;

}

void DetectorElement::toLocal(const float* X, float* Y) const {
  float r[3];

  r[0] = X[0] - m_cx;r[1] = X[1] - m_cy;r[2] = X[2] - m_cz;

  Y[0] = m_GL[0][0]*r[0] + m_GL[0][1]*r[1] + m_GL[0][2]*r[2]; 
  Y[1] = m_GL[1][0]*r[0] + m_GL[1][1]*r[1] + m_GL[1][2]*r[2];
  Y[2] = m_GL[2][0]*r[0] + m_GL[2][1]*r[1] + m_GL[2][2]*r[2];
}

void DetectorElement::toGlobal(const float* X, float* Y) const {

  float r[3] = {m_cx, m_cy, m_cz};

  r[0] += m_LG[0][0]*X[0] + m_LG[0][1]*X[1];
  r[1] += m_LG[1][0]*X[0] + m_LG[1][1]*X[1];
  r[2] += m_LG[2][0]*X[0] + m_LG[2][1]*X[1];

  Y[0] = r[0];Y[1] = r[1];Y[2] = r[2];

}

void DetectorElement::rotateToGlobal(const float* X, float* Y) const {

  for(int i=0;i<3;i++) {
    Y[i] = 0.0;
    for(int k=0;k<3;k++) Y[i] += m_LG[i][k]*X[k];
  }
}


void DetectorElement::getTransformLG(float m[3][3]) const {
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++) m[i][j] = m_LG[i][j];
}

void DetectorElement::getTransformGL(float g[3][3]) const {
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++) g[i][j] = m_GL[i][j];
}


void DetectorElement::print() const {
  std::cout<<" ("<<m_volume_id<<","<<m_layer_id<<","<<m_module_id<<") at "<<m_cx<<" "<<m_cy<<" "<<m_cz<<" ";
}


void DetectorElement::transformTrackStateToLocal(float* Re, const float* P, const float* G, float* J) const {

  double r[3] = {P[0]-m_cx,P[1]-m_cy,P[2]-m_cz};

  Re[0] = r[0]*m_GL[0][0]+r[1]*m_GL[0][1]+r[2]*m_GL[0][2];
  Re[1] = r[0]*m_GL[1][0]+r[1]*m_GL[1][1]+r[2]*m_GL[1][2];
  Re[2] = atan2f(P[4],P[3]);
  Re[3] = acosf(P[5]);
  Re[4] = P[6];

  float alpha[3];//direction vector in the local coordinate system

  rotateToLocal(P+3,alpha);

  alpha[2] = 1/alpha[2];
  
  float dV0[3], dV1[3], dV2[3], dV3[3], dV4[3];

  //transforming the jacobian columns
  
  rotateToLocal(G,   dV0);
  rotateToLocal(G+7, dV1);
  rotateToLocal(G+14,dV2);
  rotateToLocal(G+21,dV3);
  rotateToLocal(G+28,dV4);

  float L0 = alpha[0]*alpha[2];
 
  J[0] = dV0[0] - dV0[2]*L0;
  J[1] = dV1[0] - dV1[2]*L0;
  J[2] = dV2[0] - dV2[2]*L0;
  J[3] = dV3[0] - dV3[2]*L0;
  J[4] = dV4[0] - dV4[2]*L0;

  float L1 = alpha[1]*alpha[2];
 
  J[5] = dV0[1] - dV0[2]*L1;
  J[6] = dV1[1] - dV1[2]*L1;
  J[7] = dV2[1] - dV2[2]*L1;
  J[8] = dV3[1] - dV3[2]*L1;
  J[9] = dV4[1] - dV4[2]*L1;

  float Vt = 1/(1 - P[5]*P[5]);

  float a =  P[3]*Vt;
  float b = -P[4]*Vt;
  
  alpha[2] *= P[6];
  
  float T = (a*G[36]+b*G[35])*alpha[2];//36, 35
  
  J[12] = a*G[18] + b*G[17] - dV2[2]*T;//18, 17
  J[13] = a*G[25] + b*G[24] - dV3[2]*T;//25, 24
  J[14] = a*G[32] + b*G[31] - dV4[2]*T;//32, 31
  
  float sqV = -sqrt(Vt);
  T = sqV*G[37]*alpha[2];//37
  
  J[17] = sqV*G[19] - dV2[2]*T;//19
  J[18] = sqV*G[26] - dV3[2]*T;//26
  J[19] = sqV*G[33] - dV4[2]*T;//33

  J[20] = 1.0;

}

bool DetectorElement::inOverlap(const float* pos) const {

  float cut_v = 0.9;
  float cut_u = 0.9;

    
  if(m_volume_id == 8) {
    cut_v = 0.7;
    cut_u = 0.4;
  }

  if(m_volume_id == 13) {
    cut_v = 0.7;
    cut_u = 0.35;
  }

  if(m_volume_id == 7 || m_volume_id == 9) {
    cut_v = 0.3;
    cut_u = 0.2;
  }

  bool inSide = fabs(pos[1]) < cut_v*m_module_hv;

  float h = 0.0;

  if(m_volume_id == 7 || m_volume_id == 9) {

    float r_top = 0.0;
    float r_bot = 0.2;

    float minhu = (1-2*r_bot)*m_module_minhu;
    float maxhu = (1-2*r_top)*m_module_maxhu;
    h = minhu + (maxhu - minhu)*(pos[1] - m_module_hv)/(2*m_module_hv);

  } else {
    h = m_module_minhu + (m_module_maxhu - m_module_minhu)*(pos[1] - m_module_hv)/(2*m_module_hv);
  }
  inSide  = inSide && (fabs(pos[0]) < cut_u*h);

  return !inSide;
}

bool DetectorElement::isInSide(const float* pos, float tol) const {

  bool inSide = fabs(pos[1]) < (1-tol)*m_module_hv;

  if(!inSide) return false;

  float h = m_module_minhu + (m_module_maxhu - m_module_minhu)*(pos[1] - m_module_hv)/(2*m_module_hv);

  inSide  = (fabs(pos[0]) < (1-tol)*h);

  return inSide;
}

void Layer::finalize() {

  const float border = 1.0;
  
  int type = m_summary.m_type;
  
  if(type == -1 || type == 1){//endcap disks

    float zSum = 0.0;
    float minR = 100000.0;
    float maxR = 0.0;
    
    int nModules = 0;
    
    for(std::vector<const DETECTOR_ELEMENT*>::iterator it=m_deVec.begin();it!=m_deVec.end();++it) {
      zSum += (*it)->m_cz;

      float rc = sqrt((*it)->m_cx*(*it)->m_cx + (*it)->m_cy*(*it)->m_cy);
      float rmin = rc - (*it)->m_module_hv - border;
      float rmax = rc + (*it)->m_module_hv + border;
      if(minR > rmin) minR = rmin;
      if(maxR < rmax) maxR = rmax;
      nModules++;
    }
    m_summary.m_refCoord = zSum/nModules;
    m_summary.m_minCoord = minR;
    m_summary.m_maxCoord = maxR;
  }
  if(type == 0) {//barrel cylinders
    
    float rSum = 0.0;
    float minZ = 1000000.0;
    float maxZ = -1000000.0;

    int nModules = 0;
    
    for(std::vector<const DETECTOR_ELEMENT*>::iterator it=m_deVec.begin();it!=m_deVec.end();++it) {
      rSum += sqrt((*it)->m_cx*(*it)->m_cx + (*it)->m_cy*(*it)->m_cy);
      float zmin = (*it)->m_cz - (*it)->m_module_hv - border;
      float zmax = (*it)->m_cz + (*it)->m_module_hv + border;
      if(minZ > zmin) minZ = zmin;
      if(maxZ < zmax) maxZ = zmax;
      nModules++;
    }
    m_summary.m_refCoord = rSum/nModules;
    m_summary.m_minCoord = minZ;
    m_summary.m_maxCoord = maxZ;
  }

}

Volume::Volume(int t) : m_type(t) {

}

Volume::~Volume() {
  for(std::map<unsigned int, LAYER*>::iterator it=m_layMap.begin();it!=m_layMap.end();++it) {
    delete (*it).second;
  }
  m_layMap.clear();
}

void Volume::addNewDetectorElement(const DETECTOR_ELEMENT* p) {

  std::map<unsigned int, LAYER*>::iterator it=m_layMap.find(p->m_layer_id);

  LAYER* pL = 0;
  if(it == m_layMap.end()) {
    pL = new LAYER(m_type, p->m_volume_id, p->m_layer_id);
    m_layMap.insert(std::pair<unsigned int, LAYER*>(p->m_layer_id, pL));
  } else pL = (*it).second;
      
  pL->addNewDetectorElement(p);

}

void Volume::finalize() {
  m_layVec.clear();
  for(std::map<unsigned int, LAYER*>::iterator it=m_layMap.begin();it!=m_layMap.end();++it) {
    (*it).second->finalize();
    m_layVec.push_back((*it).second);
  }
  std::sort(m_layVec.begin(), m_layVec.end(), LAYER::Compare());
}

void Volume::getLayerIds(std::vector<unsigned int>& v) const {
  v.clear();
  for(std::map<unsigned int, LAYER*>::const_iterator it=m_layMap.begin();it!=m_layMap.end();++it) {
    v.push_back((*it).first);
  }
}

void Volume::getLayerSummary(std::vector<LAYER_SUMMARY>& v) const {
  v.clear();
  for(std::map<unsigned int, LAYER*>::const_iterator it=m_layMap.begin();it!=m_layMap.end();++it) {
    v.push_back((*it).second->m_summary);
  }
}

const LAYER* Volume::getLayer(unsigned int id) const {
  std::map<unsigned int, LAYER*>::const_iterator it=m_layMap.find(id);
  if(it==m_layMap.end()) return 0;
  return (*it).second;
}

int Volume::getLayerIntersection(float r1, float z1, float r2, float z2,
				 const POINT_3D& p1, const POINT_3D& p2, TRAJECTORY& vxp) const {
  
  const int INSIDE = 0; // 0000
  const int LEFT = 1;   // 0001
  const int RIGHT = 2;  // 0010
  const int BOTTOM = 4; // 0100
  const int TOP = 8;    // 1000

  float r[2] = {r1, r2};
  float z[2] = {z1, z2};
  
  int codes[2] = {INSIDE, INSIDE};

  for(int i=0;i<2;i++) {
    if(z[i] < m_minZ) codes[i] |= LEFT;
    else
      if(z[i] > m_maxZ) codes[i] |= RIGHT;

    if(r[i] < m_minR) codes[i] |= BOTTOM;
    else
      if(r[i] > m_maxR) codes[i] |= TOP;
  }

  if (codes[0] & codes[1]) return 0;
  int rc = 1;
  if(!(codes[0] | codes[1])) rc = 2;
  
  float L2 = 1.0/((z2-z1)*(z2-z1) + (r2-r1)*(r2-r1));
  float A = L2*(z2-z1);
  float B = L2*(r2-r1);

  float tz = (z2-z1)/(r2-r1);
  float tr = (r2-r1)/(z2-z1);

  int vol_id = (*m_layVec.begin())->m_summary.m_vol_id;
  unsigned int layerKey = 1000*vol_id;
  
  for(std::vector<const LAYER*>::const_iterator it=m_layVec.begin();it!=m_layVec.end();++it) {

    const LAYER* pL = (*it);

    int type = pL->m_summary.m_type;
    int lay_id = pL->m_summary.m_lay_id;
    
    if(type == 0){//barrel

      float R = pL->m_summary.m_refCoord;

      if(r1 < R && r2 < R) return rc;
      if(r1 > R && r2 > R) continue;

      float zp = z1 + (R-r1)*tz;
      float Zmin = pL->m_summary.m_minCoord;
      float Zmax = pL->m_summary.m_maxCoord;
      if(zp>=Zmin && zp<=Zmax) {

	float delta_x = zp - z1;
	float delta_y = R - r1;
	float lambda = delta_x*A + delta_y*B;
	
	float x,y,z;
	x = p1.m_x + lambda*(p2.m_x - p1.m_x);
	y = p1.m_y + lambda*(p2.m_y - p1.m_y);
	z = p1.m_z + lambda*(p2.m_z - p1.m_z);
	vxp.addPoint(x, y, z, layerKey + lay_id);
      }
    }
    else {//endcap

      float Z = pL->m_summary.m_refCoord;

      if(z1 < Z && z2 < Z) return rc;
      if(z1 > Z && z2 > Z) continue;
      
      float rp = r1 + (Z-z1)*tr;
      float Rmin = pL->m_summary.m_minCoord;
      float Rmax = pL->m_summary.m_maxCoord;

      if(rp>=Rmin && rp<=Rmax) {

	float delta_x = Z - z1;
	float delta_y = rp - r1;
	float lambda = delta_x*A + delta_y*B;

	float x,y,z;
	x = p1.m_x + lambda*(p2.m_x - p1.m_x);
	y = p1.m_y + lambda*(p2.m_y - p1.m_y);
	z = p1.m_z + lambda*(p2.m_z - p1.m_z);
	vxp.addPoint(x, y, z, layerKey + lay_id);
      }

    }
  }

  return rc;
}


Geometry::Geometry(std::ifstream& inFile) {
  
  while(!inFile.eof()) {
    unsigned int vol_id, lay_id, mod_id;
    inFile.read((char*)&vol_id, sizeof(vol_id));
    if(inFile.eof()) break;
    inFile.read((char*)&lay_id, sizeof(lay_id));
    inFile.read((char*)&mod_id, sizeof(mod_id));

    DETECTOR_ELEMENT* newDE = new DETECTOR_ELEMENT(vol_id, lay_id, mod_id);

    int address = (vol_id-7)*7 + (lay_id >> 1);
    m_deArray[address][mod_id] = newDE;

    inFile.read((char*)&newDE->m_cx, sizeof(newDE->m_cx));
    inFile.read((char*)&newDE->m_cy, sizeof(newDE->m_cy));
    inFile.read((char*)&newDE->m_cz, sizeof(newDE->m_cz));

    inFile.read((char*)&newDE->m_rot_xu, sizeof(newDE->m_rot_xu));
    inFile.read((char*)&newDE->m_rot_xv, sizeof(newDE->m_rot_xv));
    inFile.read((char*)&newDE->m_rot_xw, sizeof(newDE->m_rot_xw));
    
    inFile.read((char*)&newDE->m_rot_yu, sizeof(newDE->m_rot_yu));
    inFile.read((char*)&newDE->m_rot_yv, sizeof(newDE->m_rot_yv));
    inFile.read((char*)&newDE->m_rot_yw, sizeof(newDE->m_rot_yw));
    
    inFile.read((char*)&newDE->m_rot_zu, sizeof(newDE->m_rot_zu));
    inFile.read((char*)&newDE->m_rot_zv, sizeof(newDE->m_rot_zv));
    inFile.read((char*)&newDE->m_rot_zw, sizeof(newDE->m_rot_zw));


    inFile.read((char*)&newDE->m_module_t, sizeof(newDE->m_module_t));
    inFile.read((char*)&newDE->m_module_minhu, sizeof(newDE->m_module_minhu));
    inFile.read((char*)&newDE->m_module_maxhu, sizeof(newDE->m_module_maxhu));

    inFile.read((char*)&newDE->m_module_hv, sizeof(newDE->m_module_hv));
    inFile.read((char*)&newDE->m_pitch_u, sizeof(newDE->m_pitch_u));
    inFile.read((char*)&newDE->m_pitch_v, sizeof(newDE->m_pitch_v));

    newDE->prepare();
    
    if(!addNewDetectorElement(newDE)) {
      std::cout<<"ERROR in reading Geometry"<<std::endl;
      break;
    }

  }
  
  std::cout<<"Geometry has "<<m_deMap.size()<<" detector elements"<<std::endl;
  
  inFile.close();

  finalize();

  for(std::map<unsigned int, VOLUME*>::iterator vIt = m_volMap.begin();vIt!=m_volMap.end();++vIt) {
    VOLUME* pV = (*vIt).second;

    switch ((*vIt).first) {

    case 7 : pV->initializeBounds(29, 177, -1550, -550); break;
    case 8 : pV->initializeBounds(29, 172, -520,   520); break;
    case 9 : pV->initializeBounds(29, 177,  550,   1550); break;
    case 12: pV->initializeBounds(239, 710, -2960, -1200); break;
    case 13: pV->initializeBounds(260, 660, -1100, 1100); break;
    case 14: pV->initializeBounds(239, 710, 1200, 2960); break;
    case 16: pV->initializeBounds(720, 1051, -2960, -1200); break;
    case 17: pV->initializeBounds(800, 1050, -1100, 1100); break;
    case 18: pV->initializeBounds(720, 1051, 1200, 2960); break;
    default: std::cout<<"setting bounds error: no such volume "<<(*vIt).first<<std::endl;
      
    }
  }
  
  /*
  for(std::map<unsigned int, VOLUME*>::iterator vIt = m_volMap.begin();vIt!=m_volMap.end();++vIt) {
    VOLUME* pV = (*vIt).second;

    std::vector<LAYER_SUMMARY> vs;
    pV->getLayerSummary(vs);

    for(std::vector<LAYER_SUMMARY>::iterator it = vs.begin();it!=vs.end();++it) {
      (*it).print();
    }
  }
  */
 
}

Geometry::~Geometry() {
  for(std::map<unsigned long, const DETECTOR_ELEMENT*>::iterator it = m_deMap.begin();it != m_deMap.end();++it) {
    delete (*it).second;
  }
  m_deMap.clear();
  for(std::map<unsigned int, VOLUME*>::iterator it = m_volMap.begin();it!=m_volMap.end();++it) {
    delete (*it).second;
  }
  m_volMap.clear();
}

void Geometry::finalize() {
  for(std::map<unsigned int, VOLUME*>::iterator it = m_volMap.begin();it!=m_volMap.end();++it) {
    (*it).second->finalize();
  }
}

unsigned long Geometry::calculateKey(unsigned int v, unsigned int l, unsigned int m) const {

  unsigned long key = m;
  key += 10000l*l;
  key += 1000000l*v;
  return key;
}

bool Geometry::addNewDetectorElement(const DETECTOR_ELEMENT* p) {

  const int vtypes[19] = {2,2,2,2,2,2,2,-1,0,1,2,2,-1,0,1,2,-1,0,1};
  
  unsigned long key = calculateKey(p->m_volume_id, p->m_layer_id, p->m_module_id);
  std::map<unsigned long, const DETECTOR_ELEMENT*>::iterator it = m_deMap.find(key);
  if(it!=m_deMap.end()) return false;
  m_deMap.insert(std::pair<unsigned long, const DETECTOR_ELEMENT*>(key, p));

  std::map<unsigned int, VOLUME*>::iterator vIt = m_volMap.find(p->m_volume_id);
  VOLUME* pV = 0;
  if(vIt==m_volMap.end()) {
    pV = new VOLUME(vtypes[p->m_volume_id]);
    m_volMap.insert(std::pair<unsigned int, VOLUME*>(p->m_volume_id, pV));
  } else pV = (*vIt).second;

  pV->addNewDetectorElement(p);


  return true;
}

const DETECTOR_ELEMENT* Geometry::getDetectorElement(unsigned int v, unsigned int l, unsigned int m) const {
  int address = (v-7)*7 + (l >> 1);
  return m_deArray[address][m];
}

void Geometry::getVolumeIds(std::vector<unsigned int>& v) const {
  for(std::map<unsigned int, VOLUME*>::const_iterator it = m_volMap.begin();it!=m_volMap.end();++it) {
    v.push_back((*it).first);
  }
}

const VOLUME* Geometry::getVolume(unsigned int id) const {
  std::map<unsigned int, VOLUME*>::const_iterator it = m_volMap.find(id);
  if(it == m_volMap.end()) return 0;
  return (*it).second;
}

bool Geometry::getLayerIntersection(const POINT_3D& p1, const POINT_3D& p2, TRAJECTORY& v) const {

  const float pixBoundMin = 20.0;

  const float lgsBoundMax = 1020.0;
  
  const float zBound = 3000.0;

  //layer crossing in r-z plane

  int r1 = sqrt(p1.m_x*p1.m_x + p1.m_y*p1.m_y);
  int r2 = sqrt(p2.m_x*p2.m_x + p2.m_y*p2.m_y);

  int z1 = p1.m_z;
  int z2 = p2.m_z;

  if(r1 > lgsBoundMax && r2 > lgsBoundMax) return false;
  if(r1 < pixBoundMin && r2 < pixBoundMin) return false;

  if(fabs(z1) > zBound && fabs(z2) > zBound) return false;


  for(std::map<unsigned int, VOLUME*>::const_iterator it = m_volMap.begin();it!=m_volMap.end();++it) {

    if((*it).second->getLayerIntersection(r1, z1, r2, z2, p1, p2, v) == 2) break;

  }

  return (v.m_nPoints != 0);
}
