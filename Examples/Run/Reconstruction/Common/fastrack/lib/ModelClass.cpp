#include<iostream>
#include<fstream>
#include<cstring>
#include<cmath>
#include<omp.h>

#include "ModelClass.h"

#include<algorithm>
#include<iterator>


ModelClass::ModelClass(const char* c1, int n1, const char* c2, int n2) : m_geo(0), m_link(0), m_tf(0) {
  std::cout<<"initializing TrackML Model ver 1.4"<<std::endl;

  std::string geoName;
  geoName.assign(c1, n1);

  std::ifstream geoFile(geoName.c_str(),std::ios::binary);
  m_geo = new GEOMETRY(geoFile);
  
  std::string linkName;
  linkName.assign(c2, n2);

  std::ifstream linkFile(linkName.c_str(),std::ios::binary);
  m_link = new LAYER_LINKER(linkFile);

  m_params = new ALGORITHM_PARAMETERS;
  
  m_tf = new TRACK_FINDER(*m_params, *m_geo, *m_link);

  m_eventCounter = 0;
}

ModelClass::~ModelClass() {
  delete m_params;
  delete m_tf;
  delete m_geo;
  delete m_link;
  std::cout<<"Model instance deleted"<<std::endl;
}

void ModelClass::importHits(int nHits, const int* id, const float* x, const float* y, const float* z,
			    const int* v, const int* l, const int* m) {
  
  //std::cout<<"received "<<nHits<<" nHits"<<std::endl;

  m_H.m_eventId = m_eventCounter++;
  
  m_H.m_nHits = nHits;  
  m_H.m_x = x;
  m_H.m_y = y;
  m_H.m_z = z;
  m_H.m_hit_index = id;
  m_H.m_vol_id = v;
  m_H.m_lay_id = l;
  m_H.m_mod_id = m;

  
  delete[] m_H.m_r;
  delete[] m_H.m_phi;
  delete[] m_H.m_mask;
  delete[] m_H.m_m1;
  delete[] m_H.m_m2;
  delete[] m_H.m_eta;
  delete[] m_H.m_distance;
  
  //delete[] m_H.m_key;
  
  m_H.m_r = new float[nHits];
  m_H.m_phi = new float[nHits];
  m_H.m_mask = new int[nHits];
  
  m_H.m_m1 = new float[nHits];
  m_H.m_m2 = new float[nHits];
  m_H.m_eta = new float[nHits];
  m_H.m_distance = new float[nHits];

  memset(&m_H.m_mask[0], 0, nHits*sizeof(int));

  int nBlocks = N_THREADS;

  omp_set_num_threads(nBlocks);

  int tid;
  
#pragma omp parallel private(tid) shared(nHits, nBlocks)
  {

    tid = omp_get_thread_num();
  

    for(int i=tid;i<nHits;i+=nBlocks) {

      m_H.m_phi[i] = atan2(y[i], x[i]);
      
      float r2 = x[i]*x[i] + y[i]*y[i];
      
      m_H.m_r[i] = sqrt(r2);
      
      float tau;
      
      tau = z[i]/m_H.m_r[i];
      m_H.m_eta[i] = -log(sqrt(1+tau*tau)-tau);
      
      m_H.m_distance[i] = sqrt(r2 + z[i]*z[i]);
      
      const DETECTOR_ELEMENT* pDE = m_geo->getDetectorElement(v[i],l[i],m[i]);
      
      float m[3];//pre-calculate measurements in the local coordinate system
      float X[3] = {x[i], y[i], z[i]};
      
      pDE->toLocal(X, m);
      
      m_H.m_m1[i] = m[0]; m_H.m_m2[i] = m[1];
      
    }
  }

}

void ModelClass::importCells(int nCells, const int* input_hit_id, const int* ch0, const int* ch1) {//, const float* value) {
  
  //std::cout<<"received "<<nCells<<" nCells"<<std::endl;
  
  m_H.m_nCells = nCells;

  delete[] m_H.m_clusterLength;
  delete[] m_H.m_clusterHeight;
  
  int* cellBegin = new int[m_H.m_nHits];
  int* cellEnd   = new int[m_H.m_nHits];

  m_H.m_clusterLength = new char[m_H.m_nHits];
  m_H.m_clusterHeight = new char[m_H.m_nHits];
    
  memset(cellBegin, 0, sizeof(int)*m_H.m_nHits);
  memset(cellEnd, 0, sizeof(int)*m_H.m_nHits);

  int prev_hit_id = 0;

  for(int idx=0;idx<nCells;idx++) {

    int hit_id = input_hit_id[idx];
    
    if(hit_id != prev_hit_id) {
      cellBegin[hit_id-1] = idx;
    }
    if(prev_hit_id != 0) {
      cellEnd[prev_hit_id-1] = idx;
    }
    prev_hit_id = hit_id;
  }

  cellEnd[prev_hit_id-1] = nCells;

  int nBlocks = N_THREADS;

  omp_set_num_threads(nBlocks);

  int nHits = m_H.m_nHits;

  int tid;
  
#pragma omp parallel private(tid) shared(nHits, nBlocks)
  {

    tid = omp_get_thread_num();

  for(int i=tid;i<nHits;i+=nBlocks) {

    int nC = cellEnd[i]-cellBegin[i];
    
    if(nC==1) {
      m_H.m_clusterLength[i] = 1;
      m_H.m_clusterHeight[i] = 1;
      continue;
    }
    int k=cellBegin[i];
    int min_ch0 = ch0[k];
    int max_ch0 = ch0[k];
    int min_ch1 = ch1[k];
    int max_ch1 = ch1[k];
    for(;k<cellEnd[i];k++) {
      if(min_ch0 > ch0[k]) min_ch0 = ch0[k];
      if(max_ch0 < ch0[k]) max_ch0 = ch0[k];
      if(min_ch1 > ch1[k]) min_ch1 = ch1[k];
      if(max_ch1 < ch1[k]) max_ch1 = ch1[k];
    }
    m_H.m_clusterLength[i] = max_ch1-min_ch1+1;
    m_H.m_clusterHeight[i] = max_ch0-min_ch0+1;    
  }

  }

  delete[] cellBegin;
  delete[] cellEnd;

}

void ModelClass::findTracks(int* tid) {

  tid[0] = 0;
  tid[m_H.m_nHits-1] = 0;
  
  std::vector<NEW_TRACK*> vTracks, vTracks1, vTracks2, vTracks3;
  
  m_tf->prepare(&m_H);
  m_tf->run(1, &m_H, vTracks1, 7);
  m_tf->run(2, &m_H, vTracks2, 5);
  m_tf->run(4, &m_H, vTracks3, 2);
  m_tf->clean();
  
  
  vTracks.insert(vTracks.end(), std::make_move_iterator(vTracks1.begin()),
		 std::make_move_iterator(vTracks1.end()));
  vTracks.insert(vTracks.end(), std::make_move_iterator(vTracks2.begin()),
		 std::make_move_iterator(vTracks2.end()));
  vTracks.insert(vTracks.end(), std::make_move_iterator(vTracks3.begin()),
		 std::make_move_iterator(vTracks3.end()));

  //m_tf->groupTracks(vTracks);
  
  //std::cout<<"found "<<vTracks.size()<<" tracks"<<std::endl;

  generateHitAssignment(vTracks, tid);

  vTracks.clear();
  
  //for(std::vector<TRACK*>::iterator it=vTracks.begin();it!=vTracks.end();++it) delete (*it);
}

void ModelClass::generateHitAssignment(std::vector<NEW_TRACK*>& vTracks, int* H2T) {

  //create unique hit-to-track association
 
  int dummyTrack = vTracks.size() + 1;

  std::fill_n(H2T, m_H.m_nHits, dummyTrack);
  
  int trackId = 0;
  
  for(std::vector<NEW_TRACK*>::iterator it = vTracks.begin();it!=vTracks.end();++it) {
    
    trackId++;
    
    for(int k=0;k<(*it)->m_nHits;k++) {

    int h = (*it)->m_hits[k];
    
      int tid = H2T[h];
      if(tid < dummyTrack) continue;
      
      H2T[h] = trackId;
      
    }
  }
}

void ModelClass::importParticles(int nPart, const unsigned long* pid, const int* type, 
				 const float* vx, const float* vy, const float* vz, 
				 const float* px, const float* py, const float* pz, 
				 const int* q, const int* nhits) {


  m_P.clear();
  
  for(int k=0;k<nPart;k++) {

    PARTICLE p;
    
    p.m_particleType = type[k];

    p.m_vx = vx[k];
    p.m_vy = vy[k];
    p.m_vz = vz[k];
    p.m_px = px[k];
    p.m_py = py[k];
    p.m_pz = pz[k];
    p.m_q = q[k];
    p.m_nhits = nhits[k];

    m_P.insert(std::pair<unsigned long, PARTICLE>(pid[k], p));

  }
}


void ModelClass::importTruth(int nHits, const int* hit_id, const float* x, const float* y, const float* z, const float* px, const float* py, const float* pz, const unsigned long* particle, const float* w) {

  if(nHits != m_H.m_nHits) {
    std::cout<<"nHits mismatch in importTruth"<<std::endl;
  }

  m_H.m_true_hit_id = hit_id;
  m_H.m_true_x = x;
  m_H.m_true_y = y;
  m_H.m_true_z = z;
  m_H.m_true_px = px;
  m_H.m_true_py = py;
  m_H.m_true_pz = pz;
  m_H.m_particle = particle;
  m_H.m_weight = w;
}

void ModelClass::train() {
  
  m_tf->train(&m_H, m_P);

}
