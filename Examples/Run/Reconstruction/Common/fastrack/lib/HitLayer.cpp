#include "DataStructures.h"

HitLayer::HitLayer(const LAYER_SUMMARY& ls, float deta) : m_info(ls) {

  float min_eta, max_eta, tau;

  if(ls.m_type == 0) {//barrel
    tau = ls.m_maxCoord/ls.m_refCoord;
    max_eta = -log(sqrt(1+tau*tau)-tau);
    tau = ls.m_minCoord/ls.m_refCoord;
    min_eta = -log(sqrt(1+tau*tau)-tau);
  }
  else {
    tau = ls.m_refCoord/ls.m_minCoord;
    max_eta = -log(sqrt(1+tau*tau)-tau);
    tau = ls.m_refCoord/ls.m_maxCoord;
    min_eta = -log(sqrt(1+tau*tau)-tau);
    if(min_eta>max_eta) std::swap(min_eta, max_eta);
  }
  m_minEta = min_eta;
  m_maxEta = max_eta;
  m_nBins = 1+(int)((m_maxEta-m_minEta)/deta);
  m_deltaEta = (m_maxEta-m_minEta)/m_nBins;

  m_bins = new HIT_BIN[m_nBins];
  m_binCoords = new float[m_nBins];
  m_binRads   = new float[m_nBins];
  m_binBoundMin  = new float[m_nBins];
  m_binBoundMax  = new float[m_nBins];
  
  for(int i=0;i<m_nBins;i++) {

    float eta = m_minEta+(i+0.5)*m_deltaEta;
    float thetah = 2*atan(exp(-eta));

    float eta1 = eta - 0.5*m_deltaEta;
    float eta2 = eta + 0.5*m_deltaEta;

    float theta1 = 2*atan(exp(-eta1));
    float theta2 = 2*atan(exp(-eta2));
    
    if(ls.m_type == 0) {//barrel
      m_binCoords[i] = ls.m_refCoord*cos(thetah)/sin(thetah);//z-coordinate

      m_binBoundMin[i] = ls.m_refCoord*cos(theta1)/sin(theta1);
      m_binBoundMax[i] = ls.m_refCoord*cos(theta2)/sin(theta2);
      
      m_binRads[i] = ls.m_refCoord;
    }
    else {
      m_binCoords[i] = ls.m_refCoord*sin(thetah)/cos(thetah);//r-coordinate

      m_binBoundMin[i] = ls.m_refCoord*sin(theta1)/cos(theta1);
      m_binBoundMax[i] = ls.m_refCoord*sin(theta2)/cos(theta2);
      
      m_binRads[i] = m_binCoords[i];
    }
  }

}

HitLayer::~HitLayer() {
  delete[] m_bins;
  delete[] m_binCoords;
  delete[] m_binRads;
  delete[] m_binBoundMin;
  delete[] m_binBoundMax;
}

void HitLayer::sort() {
  for(int i=0;i<m_nBins;i++) m_bins[i].sort();
}

void HitLayer::clear() {
  for(int i=0;i<m_nBins;i++) m_bins[i].clear();
}

void HitLayer::fillIndices() {

  bool isPixelEndcap = (m_info.m_type != 0) && (m_info.m_vol_id == 7 || m_info.m_vol_id == 9);
  
  for(int i=0;i<m_nBins;i++) m_bins[i].fillIndices(isPixelEndcap);
}

void HitLayer::deleteUsedNodes(const HITS* hits) {
  for(int i=0;i<m_nBins;i++) m_bins[i].deleteUsedNodes(hits);
}

void HitLayer::prepareL1Nodes() {

  const float Table1[6]      = {0.0, 1.0265, 1.509946, 2.12928, 2.12928, 1.509946};
  const float Table2[6]      = {40.0, 40.0, 30.1, 9.05956, 9.05956, 9.05956};
  const float absTauCut[13]  = {0.0, 0.453, 0.73363, 1.14457, 1.69838, 2.01427, 2.79041, 3.26816, 3.62686, 4.23419, 4.98758, 5.24827, 5.57847};
  const float absTauCut2[12] = {0.0, 0.0,   0.001,   0.211547, 0.221779, 0.226903, 0.237169, 0.247458, 0.283673, 0.299297, 0.399962, 0.410752};

  int vol1 = m_info.m_vol_id;

  bool isPixelBarrel  = vol1 == 8;
  bool isShsBarrel    = vol1 == 13;
  bool isShsEndcap    = vol1 == 12 || vol1 == 14;  
  bool isPixelEndcap  = vol1 == 7  || vol1 == 9;
  
  for(int i=0;i<m_nBins;i++) {

    std::vector<NODE*>& nodes = m_bins[i].m_nodes;
    std::vector<NODE*>& l1v = m_bins[i].m_l1nodes;
    l1v.clear();

    for(std::vector<NODE*>::iterator n1It = nodes.begin();n1It != nodes.end();++n1It) {

      NODE& n1 = *(*n1It);
      
      if(n1.m_nCellsT > 3) continue;
      
      if(isPixelEndcap && n1.m_nCells > 7) continue;
      if(isShsBarrel && n1.m_nCells > 3) continue;
      if(isShsEndcap && n1.m_nCells > 2) continue;
      
      l1v.push_back(*n1It);
      
      if((*n1It)->m_minCutOnPar0 > -100) continue;//lazy initialization
      
      int nCells  = n1.m_nCells;
      
      (*n1It)->m_minCutOnPar0 = 0.0;
      (*n1It)->m_maxCutOnPar0 = 40.0;

      if(isPixelEndcap) {
	if(nCells<6) {
	  (*n1It)->m_minCutOnPar0 = Table1[nCells];
	  (*n1It)->m_maxCutOnPar0 = Table2[nCells];
	}
	else {
	  (*n1It)->m_minCutOnPar0 = 0.521095;
	  (*n1It)->m_maxCutOnPar0 = 5.466229;
	}
      }
      
      if(isShsBarrel) {
	if(nCells == 1) {
	  (*n1It)->m_maxCutOnPar0 = 4.732;
	}
	if(nCells == 2) {
	  (*n1It)->m_maxCutOnPar0 = 4.93696;
	}
	if(nCells == 3) {
	  (*n1It)->m_maxCutOnPar0 = 4.93696;
	}
      }
      else {
	if(nCells == 1) {
	  (*n1It)->m_minCutOnPar0 = 1.0275;
	}
	
	if(nCells == 2) {
	  (*n1It)->m_minCutOnPar0 = 1.509946;
	  (*n1It)->m_maxCutOnPar0 = 30.1;
	}
      }
	  
      if(isPixelBarrel) {
	
	(*n1It)->m_minCutOnPar0 = nCells<11 ? absTauCut2[nCells] : absTauCut2[11];

	if(nCells<13) {
	  (*n1It)->m_maxCutOnPar0 = absTauCut[nCells];
	}
	else {
	  if(nCells<15) {
	    (*n1It)->m_maxCutOnPar0 = 5.466229;
	  }
	  else {
	    if(nCells < 25) {
	      (*n1It)->m_maxCutOnPar0 = 10.01787;
	    }
	  }
	}
      }
    }
  }

}

/*
int HitLayer::getBinIndex(float r, float z) {

  float tau, eta;

  tau = z/r;
  eta = -log(sqrt(1+tau*tau)-tau);

  if(eta<m_minEta) return 0;
  if(eta>m_maxEta) return m_nBins-1;

  int bin = (int)((eta-m_minEta)/m_deltaEta);

  //std::cout<<"bin="<<bin<<" r="<<r<<" z="<<z<<" nBins="<<m_nBins<<std::endl;
  
  return bin;
}
*/

int HitLayer::getBinIndex(float eta) {

  
  if(eta<m_minEta) return 0;
  if(eta>m_maxEta) return m_nBins-1;

  int bin = (int)((eta-m_minEta)/m_deltaEta);

  //std::cout<<"bin="<<bin<<" r="<<r<<" z="<<z<<" nBins="<<m_nBins<<std::endl;
  
  return bin;
}

float HitLayer::getBinRadius(int idx) {
  return m_binRads[idx];
}

void HitLayer::getBinTaus(float r1, float z1, int binIndex, float& tmin, float& tmax) {

  if(m_info.m_type == 0) {//barrel

    float zb1 = m_binBoundMin[binIndex];
    float zb2 = m_binBoundMax[binIndex];
    float rb = m_binRads[binIndex]; 
    float delta = 1/(rb-r1);
    tmin = fabs((zb1-z1)*delta); 
    tmax = fabs((zb2-z1)*delta); 

  }
  else {

    float rb1 = m_binBoundMin[binIndex];
    float rb2 = m_binBoundMax[binIndex];
    float zb = m_info.m_refCoord;
    tmax = fabs((zb-z1)/(rb1-r1)); 
    tmin = fabs((zb-z1)/(rb2-r1)); 
    
  }

  if(tmax < tmin) std::swap(tmax, tmin);

}

//void HitLayer::importNewHit(int hitIndex, const HITS* hits) {
//  NODE* p = new NODE(hitIndex, hits, m_info.m_type);
//  float eta = hits->m_eta[hitIndex];
//  int b = getBinIndex(eta);
//  m_bins[b].importNewNode(p);  
//}


bool HitLayer::getBinRange(float t1, float t2, int& i1, int& i2) {

  float emin = -log(sqrt(1+t1*t1)-t1);
  float emax = -log(sqrt(1+t2*t2)-t2);
  
  i1 = -1; i2 = -1;

  if(std::isnan(emax) || std::isnan(emin)) return false;

  if(emax < m_minEta || emin > m_maxEta) return false;

  if(emax > m_maxEta) emax = m_maxEta;
  if(emin < m_minEta) emin = m_minEta;
  
  i1 = int((emin-m_minEta)/m_deltaEta);
  i2 = int((emax-m_minEta)/m_deltaEta);

  if(i1<0) i1 = 0;

  if(i2>=m_nBins) i2 = m_nBins-1;
  if(i1>=m_nBins) i1 = m_nBins-1;
  
  return true;
}

unsigned int HitLayer::numberOfHits() {

  unsigned int n=0;
    
  for(int i=0;i<m_nBins;i++) n += m_bins[i].m_nodes.size();

  return n;
}
