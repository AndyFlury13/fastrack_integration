#include<iostream>

#include "ClusterMaker.h"
#include "Event.h"
#include<map>

ClusterMaker::ClusterMaker(const HITS* h, const GEOMETRY& g) : m_hits(h), m_geo(g) {

}

ClusterMaker::~ClusterMaker() {

}

void ClusterMaker::findClusters(const ALGORITHM_PARAMETERS& params, HIT_LAYER* pL) {

  float dnCellsCut = params.cm_f[0];
  float dPhiR_max  = params.cm_f[1];
  float maxZ0_disk = params.cm_f[2];
  float maxZ0_barr = params.cm_f[3];

  bool isPixel = pL->m_info.m_vol_id < 10;
  
  float maxZ0 = pL->m_info.m_type == 0 ? maxZ0_barr : maxZ0_disk;

  if(pL->m_info.m_vol_id == 8) {
    if(pL->m_info.m_lay_id == 2) {
      dPhiR_max = params.cm_f[4]*dPhiR_max;//p4
      maxZ0 = params.cm_f[5];//p5
    }
    if(pL->m_info.m_lay_id == 4) {
      dPhiR_max = params.cm_f[6]*dPhiR_max;//p6
      maxZ0 = params.cm_f[7];//p7
    }
    if(pL->m_info.m_lay_id == 6) {
      dPhiR_max = params.cm_f[8]*dPhiR_max;//p8
      maxZ0 = params.cm_f[9];//p9
    }
    if(pL->m_info.m_lay_id == 8) {
      dPhiR_max = params.cm_f[10]*dPhiR_max;//p10
      maxZ0 = params.cm_f[11];//p11
    }
  }

  if(pL->m_info.m_vol_id == 7 || pL->m_info.m_vol_id == 9) {
    maxZ0 = params.cm_f[12];//p12
    dPhiR_max = params.cm_f[13];//p13
  }

  if(pL->m_info.m_vol_id > 10) {//short strips

    if(pL->m_info.m_type == 0) {
      dPhiR_max = params.cm_f[14];//p14
    }
    else {
      dPhiR_max = params.cm_f[15];//p15
    }
  }

  if(pL->m_info.m_vol_id == 13) {
    maxZ0 = params.cm_f[16];//p16
  }


  float surfaceRad = 0.0;
  if(pL->m_info.m_type == 0) surfaceRad = pL->m_info.m_refCoord;

  int curr_state[5000];
  int next_state[5000];
  
  for(int binIdx=0;binIdx<pL->m_nBins;binIdx++) {

    HIT_BIN& B = pL->m_bins[binIdx];

    if(B.m_nodes.size() < 2) continue;

    std::vector<HIT_CLUSTER> vClusters;

    if(pL->m_info.m_type != 0) surfaceRad = pL->m_binCoords[binIdx]; 

    for(unsigned int n1=0;n1<B.m_nodes.size()-1;n1++) {

      if(!B.m_nodes[n1]->m_inOverlap) continue;

      float z1 = B.m_nodes[n1]->m_z;

      int h1 = B.m_nodes[n1]->m_h;

      float phi1 = B.m_nodes[n1]->m_phi;
      
      float r1 = B.m_nodes[n1]->m_r;

      int mod1 = m_hits->m_mod_id[h1];
      
      unsigned int n2 = n1 + 1;

      std::vector<int> modVec;
      
      float nCells1 = m_hits->m_clusterLength[h1];

      while(n2<B.m_nodes.size()) {

	if(!B.m_nodes[n2]->m_inOverlap) break;

	int h2 = B.m_nodes[n2]->m_h;
	int mod2 = m_hits->m_mod_id[h2];
	
	if(mod1 == mod2) break;

	float z2 = B.m_nodes[n2]->m_z;

	float phi2 = B.m_nodes[n2]->m_phi;
	
	float r2 = B.m_nodes[n2]->m_r;

	float dphi = phi2 - phi1;

	if(isPixel) {
	  float nCells2 = m_hits->m_clusterLength[h2];

	  float ncDiff = fabs(nCells2-nCells1)/std::max(nCells2,nCells1);
	  
	  if(ncDiff > dnCellsCut) {
	    n2++;
	    continue;
	  }
	}
	
	if(surfaceRad*fabs(dphi) <= dPhiR_max) {

	  float z0 = (r2*z1-r1*z2)/(r2-r1);

	  if(fabs(z0)<maxZ0 && (std::find(modVec.begin(), modVec.end(), mod2) == modVec.end()) ) {//new cluster detected
	    vClusters.push_back(HIT_CLUSTER(n1, n2));
	    modVec.push_back(mod2);
	  }
	  
	}
	n2++;
      }
    }

    if(vClusters.empty()) continue;

    //merging using CCA

    unsigned int nNodes = B.m_nodes.size();

    for(unsigned int nodeIdx=0;nodeIdx< nNodes;nodeIdx++) {
      curr_state[nodeIdx] = next_state[nodeIdx] = nodeIdx;
    }

    int nChanges = vClusters.size();
    
    while(nChanges !=0) {

      for(std::vector<HIT_CLUSTER>::iterator cIt=vClusters.begin();cIt!=vClusters.end();++cIt) {

	unsigned int idx1 = (*cIt).m_idx1;
	unsigned int idx2 = (*cIt).m_idx2;

	if(curr_state[idx1] < curr_state[idx2]) next_state[idx2] = curr_state[idx1];
	else next_state[idx1] = curr_state[idx2];
      }

      nChanges = 0;
      for(unsigned int nodeIdx=0;nodeIdx < nNodes;nodeIdx++) {
	if(curr_state[nodeIdx] != next_state[nodeIdx]) nChanges++;
	curr_state[nodeIdx] = next_state[nodeIdx];
      }
    }

    for(unsigned int idx1=0;idx1 < nNodes-1;idx1++) {

      int tag1 = curr_state[idx1];
      if(tag1 == -1) continue;

      NODE* node1 = B.m_nodes.at(idx1);
      
      for(unsigned int idx2=idx1+1;idx2 < nNodes;idx2++) {
	int tag2 = curr_state[idx2];
	if(tag2 == -1) continue;

	if(tag1 == tag2) {//merging
	  NODE* node2 = B.m_nodes.at(idx2);
	  //1. transfer hit(s) from node2 to node1

	  std::copy(node2->m_hits.begin(), node2->m_hits.end(), std::back_inserter(node1->m_hits));

	  //2. mark node2 as empty

	  node2->m_h = 0;
	  node2->m_hits.clear();
	  curr_state[idx2] = -1;
	}
      }
    }
    
  }
  
  //finally, remove nodes with zero m_h, update merged coordinates
    
  for(int binIdx=0;binIdx<pL->m_nBins;binIdx++) {

    HIT_BIN& B = pL->m_bins[binIdx];
    
    if(B.m_nodes.empty()) continue;
    
    B.removeEmptyNodes();
    
    B.updateCoordinates(m_hits);
    B.calculateCellShape(m_hits);

  }

}
