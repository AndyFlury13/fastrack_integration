#include<iostream>
#include<cstring>
#include<set> 

#include "RecoEvent.h"

#include<algorithm>

TrackState::TrackState(float* P, const DETECTOR_ELEMENT* p) : m_pDE(p) {

  memcpy(&m_Rk[0], P, sizeof(m_Rk));
  memcpy(&m_Re[0], P, sizeof(m_Re));

  memset(&m_Ck[0], 0, sizeof(m_Ck));
  
  m_Ck[0] = 4.8f;
  m_Ck[2] = 4.8f;
  m_Ck[5] = 0.001f;
  m_Ck[9] = 0.001f;
  m_Ck[14] = 1.1f;

  memcpy(&m_Ce[0], &m_Ck[0], sizeof(m_Ce));


  m_chi2 = 0.0;
  m_ndof = -5;

}

void TrackState::initialize(float* P, const DETECTOR_ELEMENT* p) {

  m_pDE = p;
  m_initialized = true;
  memcpy(&m_Rk[0], P, sizeof(m_Rk));
  memcpy(&m_Re[0], P, sizeof(m_Re));

  memset(&m_Ck[0], 0, sizeof(m_Ck));
  
  m_Ck[0] = 4.8f;
  m_Ck[2] = 4.8f;
  m_Ck[5] = 0.001f;
  m_Ck[9] = 0.001f;
  m_Ck[14] = 1.1f;

  memcpy(&m_Ce[0], &m_Ck[0], sizeof(m_Ce));


  m_chi2 = 0.0;
  m_ndof = -5;
}

void TrackState::clone(const TrackState& ts) {

  m_pDE = ts.m_pDE;

  m_initialized = true;
  
  memcpy(&m_Rk[0], &ts.m_Rk[0], sizeof(m_Rk));
  memcpy(&m_Re[0], &ts.m_Re[0], sizeof(m_Re));

  memcpy(&m_Ck[0], &ts.m_Ck[0], sizeof(m_Ck));
  memcpy(&m_Ce[0], &ts.m_Ce[0], sizeof(m_Ce));

  m_chi2 = ts.m_chi2;
  m_ndof = ts.m_ndof;

}



void TrackState::resetCovariance() {

  const float factor = 2.0;

  m_Ck[0] *= factor;
  m_Ck[2] *= factor;
  m_Ck[5] *= factor;
  m_Ck[9] *= factor;
  m_Ck[14] *= factor;

  memcpy(&m_Ce[0], &m_Ck[0], sizeof(m_Ce));

  m_chi2 = 0.0;
  m_ndof = -5;
}

bool TrackState::valid() {
  if(!m_initialized) return false;
  if(std::isnan(m_chi2)) return false;
  if(m_chi2 < 0.0) return false;
  
  return true;
}

void NewTrack::initialize(int id, const int* v, int n, float q) {

  m_trackId = id;
  m_nHits = n;
  m_Q = q;
  m_ll = -1e8;
  m_status = 1;
  if(m_nHits >= MAX_HIT_ON_TRACK) m_nHits = MAX_HIT_ON_TRACK-1;
  memcpy(&m_hits[0], v, m_nHits*sizeof(int));

  m_nlayers = 0;
  m_beginTrackState.m_initialized = false;
  m_endTrackState.m_initialized = false;
}

void NewTrack::print(const HITS* hits) {
  std::cout<<"Track "<<m_trackId<<" with "<<m_nlayers<<" layers, LL="<<m_ll<<std::endl;
  
  for(int idx = 0;idx<m_nHits;idx++) {
    int h = m_hits[idx];
    std::cout<<"Hit "<<h<<" at "<<hits->m_vol_id[h]<<":"<<hits->m_lay_id[h]<<":"<<hits->m_mod_id[h]
	     <<" x="<<hits->m_x[h]<<" y="<<hits->m_y[h]<<" z="
	     <<hits->m_z[h]<<" phi="<<hits->m_phi[h]<< std::endl;
  }
  
  
}

void NewTrack::update(const HITS* hits) {

  unsigned int vl[256];
  memset(&vl[0],0,sizeof(vl));

  int nL = 0;
  
  for(int idx=0;idx<m_nHits;idx++) {

    int h = m_hits[idx];

    unsigned int layerKey = 7*(hits->m_vol_id[h]-7) + (hits->m_lay_id[h]>>1);
    if(layerKey>=255) {
      std::cout<<"Layer key = "<<layerKey<<std::endl;
    }
    if(vl[layerKey] > 0) continue;
    nL++;
    vl[layerKey] = 1;
  }
  
  m_nlayers = nL;
}

void NewTrack::sortHits(const HITS* hits) {
  std::sort(m_hits, m_hits+m_nHits, SortHitsByDistance(hits));
}

void NewTrack::maskHits(int m, const HITS* hits) const {

  for(int idx=0;idx<m_nHits;idx++) {
    int h = m_hits[idx];
    hits->m_mask[h] = m;
 }
}

void NewTrack::add(std::vector<int>& v, const HITS* hits) {

  if(m_nHits>=MAX_HIT_ON_TRACK) return;

  //merging

  std::vector<unsigned int> modVec;

  modVec.resize(m_nHits);
  
  int nHits0 = m_nHits;
  
  for(int hitIdx=0;hitIdx<nHits0;hitIdx++) {
    int h = m_hits[hitIdx];
    modVec[hitIdx] = 1000000*hits->m_vol_id[h] + 10000*hits->m_lay_id[h] + hits->m_mod_id[h];
  }

  for(std::vector<int>::iterator hIt = v.begin();hIt!=v.end();++hIt) {

    int h = *hIt;
    
    unsigned int key = 1000000*hits->m_vol_id[h] + 10000*hits->m_lay_id[h] + hits->m_mod_id[h];

    if(std::find(modVec.begin(), modVec.end(), key) != modVec.end()) continue;
	
    m_hits[m_nHits++] = h;
    if(m_nHits>=MAX_HIT_ON_TRACK) break;
  }
}

void NewTrack::getLayers2(std::vector<int>& vLayers, const HITS* hits) const {

  vLayers.clear();
  
  unsigned int vl[256];

  memset(&vl[0],0,sizeof(vl));
  
  for(int idx=0;idx<m_nHits;idx++) {

    int h = m_hits[idx];

    int vol_id = hits->m_vol_id[h];
    int lay_id = hits->m_lay_id[h];

    unsigned int th = vol_id > 10 ? 1 : 2;
    
    unsigned int layerKey = 7*(vol_id-7) + (lay_id>>1);
    
    vl[layerKey]++;
    if(vl[layerKey]>th) vLayers.push_back(1000*vol_id + lay_id);
  }
}

void NewTrack::getLayers(std::vector<int>& vLayers, const HITS* hits) const {

  vLayers.clear();
  
  unsigned int vl[256];

  memset(&vl[0],0,sizeof(vl));
  
  for(int idx=0;idx<m_nHits;idx++) {

    int h = m_hits[idx];

    int vol_id = hits->m_vol_id[h];
    int lay_id = hits->m_lay_id[h];

    unsigned int th = vol_id > 10 ? 1 : 2;
    
    unsigned int layerKey = 7*(vol_id-7) + (lay_id>>1);
    
    vl[layerKey]++;
    if(vl[layerKey]>th) vLayers.push_back(1000*vol_id + lay_id);
  }
}

void NewTrack::getModules(std::vector<int>& vM, const HITS* hits) const {
  vM.clear();
  for(int idx=0;idx<m_nHits;idx++) {

    int h = m_hits[idx];

    int vol_id = hits->m_vol_id[h];
    int lay_id = hits->m_lay_id[h];
    int mod_id = hits->m_lay_id[h];

    vM.push_back(1000000*vol_id + 10000*lay_id + mod_id);
  }
}

void NewTrack::removeDuplicates() {

  std::vector<int> hitVec;
  
  hitVec.reserve(m_nHits);
  
  for(int idx = 0;idx<m_nHits;idx++) {
    int h = m_hits[idx];

    if(std::find(hitVec.begin(), hitVec.end(), h) != hitVec.end()) continue;
  
    hitVec.push_back(h);
  }

  if((int)hitVec.size() == m_nHits) return;

  m_nHits = 0;

  for(std::vector<int>::iterator it = hitVec.begin();it!=hitVec.end();++it) {
    m_hits[m_nHits++] = (*it);
  }
}
