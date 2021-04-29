#include<iostream>
#include<fstream>
#include<cstring>

#include "DataStructures.h"

#include<set>

#include<omp.h>

HitBank::HitBank(const GEOMETRY& g, float deta) : m_geo(g), m_hits(0)  {

  std::cout<<"HitBank initialization ..."<<std::endl;

  memset(&m_layArray[0], 0, sizeof(m_layArray));
 
  std::vector<unsigned int> volIds;
  
  m_geo.getVolumeIds(volIds);

  for(std::vector<unsigned int>::iterator vIt=volIds.begin();vIt!=volIds.end();++vIt) {

    unsigned int vId = (*vIt);
    
    const VOLUME* pV = m_geo.getVolume(vId);

    std::vector<unsigned int> vL;
    pV->getLayerIds(vL);
  
    //std::cout<<"Volume "<<vId<<" has "<<vL.size()<<" layers"<<std::endl;

    for(unsigned int lIdx=0;lIdx<vL.size();lIdx++) {
      addNewLayer(vId, vL.at(lIdx), deta);
    }
  }
 
  std::cout<<"done"<<std::endl;
}

HitBank::~HitBank() {
  for(std::map<unsigned int, HIT_LAYER*>::iterator it = m_layMap.begin();it!=m_layMap.end();++it) {
    delete (*it).second;
  }
  m_layMap.clear();
}

void HitBank::addNewLayer(unsigned int vId, unsigned int lId, float deta) {

  unsigned int layerKey = 1000*vId + lId;

  const VOLUME* pV = m_geo.getVolume(vId);
  const LAYER* pL = pV->getLayer(lId);

  HIT_LAYER* pHL = new HIT_LAYER(pL->m_summary, deta);
  
  m_layMap.insert(std::pair<unsigned int, HIT_LAYER*>(layerKey, pHL));

  m_layArray[(vId-7)*7+(lId>>1)] = pHL;

}

void HitBank::getLayersFromVolume(int vId, std::vector<HIT_LAYER*>& vL) {
  for(std::map<unsigned int, HIT_LAYER*>::iterator it = m_layMap.begin();it!=m_layMap.end();++it) {
    int layerKey = (*it).first;
    int vol_id = layerKey / 1000;
    if(vol_id == vId) vL.push_back((*it).second);
  }
}

void HitBank::getLayers(std::vector<HIT_LAYER*>& vL) {
  for(std::map<unsigned int, HIT_LAYER*>::iterator it = m_layMap.begin();it!=m_layMap.end();++it) {
    vL.push_back((*it).second);
  }
}

HIT_LAYER* HitBank::getLayer(unsigned int layerKey) {
  return (*m_layMap.find(layerKey)).second;
}


void HitBank::clear() {

   for(std::map<unsigned int, HIT_LAYER*>::iterator it = m_layMap.begin();it!=m_layMap.end();++it) {
     (*it).second->clear();
   }
}

void HitBank::deleteUsedNodes(const HITS* hits) {
   
  for(std::map<unsigned int, HIT_LAYER*>::iterator it = m_layMap.begin();it!=m_layMap.end();++it) {
     (*it).second->deleteUsedNodes(hits);
  }
  
}

void HitBank::importHits(const HITS* pE, bool useMasked) {
  
  m_hits = pE;
  
  for(int i=0;i<pE->m_nHits;i++) {
    
    if(!useMasked) {
      if(pE->m_mask[i] != 0) continue;
    }
    
    importNewHit(i);

  }

  //check hits in overlapping areas

  for(std::map<unsigned int, HIT_LAYER*>::iterator it = m_layMap.begin();it!=m_layMap.end();++it) {

    HIT_LAYER* pL = (*it).second;

    for(int i=0;i<pL->m_nBins;i++) {

      HIT_BIN& B = pL->m_bins[i];

      for(unsigned int n1=0;n1<B.m_nodes.size();n1++) {

	float x1 = B.m_nodes[n1]->m_x;
	float y1 = B.m_nodes[n1]->m_y;
	float z1 = B.m_nodes[n1]->m_z;
	int h1 = B.m_nodes[n1]->m_h;
	
	float gP1[3] = {x1, y1, z1};
	
	const DETECTOR_ELEMENT* pE1 = m_geo.getDetectorElement(m_hits->m_vol_id[h1],
							       m_hits->m_lay_id[h1],
							       m_hits->m_mod_id[h1]);
	
	float lP1[3];
	
	pE1->toLocal(gP1, lP1);
	
	B.m_nodes[n1]->m_inOverlap = pE1->inOverlap(lP1);
      }
    }
  }


}

void HitBank::importNewHit(int i) {
  
  int vol_id = m_hits->m_vol_id[i];
  int lay_id = m_hits->m_lay_id[i];
  
  HIT_LAYER* pL = m_layArray[(vol_id-7)*7 + (lay_id>>1)];
  
  if(pL == 0) return;
  
  NODE* p = new NODE(i, m_hits, pL->m_info.m_type);
  
  float eta = m_hits->m_eta[i];

  int b = pL->getBinIndex(eta);

  pL->m_bins[b].m_nodes.push_back(p); 
 
}

void HitBank::sort() {
  for(std::map<unsigned int, HIT_LAYER*>::iterator it = m_layMap.begin();it!=m_layMap.end();++it) {
    (*it).second->sort();
  }
}

void HitBank::fillIndices() {
  for(std::map<unsigned int, HIT_LAYER*>::iterator it = m_layMap.begin();it!=m_layMap.end();++it) {
    (*it).second->fillIndices();
  }
}

unsigned int HitBank::numberOfHits() {

  unsigned int n=0;

  for(std::map<unsigned int, HIT_LAYER*>::iterator it = m_layMap.begin();it!=m_layMap.end();++it) {
    n += (*it).second->numberOfHits();
  }

  return n;
}

void HitBank::getConnectingNodes(std::vector<NODE*>& vn) {
  
  vn.reserve(numberOfHits());
  
  for(std::map<unsigned int, HIT_LAYER*>::iterator it = m_layMap.begin();it!=m_layMap.end();++it) {

    HIT_LAYER* pL = (*it).second;

    for(int ib=0;ib<pL->m_nBins;ib++) {
     
      HIT_BIN& b = pL->m_bins[ib];

      for(std::vector<NODE*>::iterator nIt = b.m_nodes.begin();nIt!=b.m_nodes.end();++nIt) {
	if((*nIt)->m_in.empty()) continue;
	if((*nIt)->m_out.empty()) continue;
	vn.push_back(*nIt);
      }
    }
  }
}


void HitBank::findClosestHits(const std::vector<int>& layVec, const TRAJECTORY& vP, std::vector<int>& vHits, const HITS* hits) {

  float delta_z_m = 9.0;//was 7.0
  float delta_z_p = 9.0;

  float delta_R_m = 6.0;//was 9.0
  float delta_R_p = 6.0;

  float phiWidth = 0.010;//was 0.010

  float chi2Cut = 200.0;//was 200.0

  const double sigmaU_pix = 0.088;//was 0.085
  const double sigmaU_pix2 = 1/(sigmaU_pix*sigmaU_pix);
  
  const double sigmaV_pix = 0.16;//was 0.16
  const double sigmaV_pix2 = 1/(sigmaV_pix*sigmaV_pix);
  
  const double sigmaU_shs = 0.25;//was 0.25
  const double sigmaU_shs2 = 1/(sigmaU_shs*sigmaU_shs);

  const double sigmaV_shs = 0.80;//was 0.79
  const double sigmaV_shs2 = 1/(sigmaV_shs*sigmaV_shs);

  const double sigmaU_lgs = 0.33;//was 0.33
  const double sigmaU_lgs2 = 1/(sigmaU_lgs*sigmaU_lgs);
  
  const double sigmaV_lgs = 2.8;//was 2.5
  const double sigmaV_lgs2 = 1/(sigmaV_lgs*sigmaV_lgs);
  
  const float tol = 0.0;
  
  vHits.clear();

  if(vP.m_nPoints<2) return;

  int it1 = 0;
  int it2 = 1;

  TRAJECTORY vOut;
  
  for(;it2<vP.m_nPoints;it1++,it2++) {
    m_geo.getLayerIntersection(vP.m_p[it1], vP.m_p[it2], vOut);
  }

  for(int impactIdx=0;impactIdx<vOut.m_nPoints;impactIdx++) {

    unsigned int key = vOut.m_p[impactIdx].m_key;

    if(std::find(layVec.begin(), layVec.end(), (int)key) != layVec.end()) continue;

    POINT_3D& imp = vOut.m_p[impactIdx];
    
    HIT_LAYER* pL = getLayer(key);
    unsigned int vol_id = key / 1000;
    unsigned int lay_id = key % 1000;

    float sigmaV2 = 0.0;
    float sigmaU2 = 0.0;
    
    if(vol_id < 10) {//pixel
      sigmaU2 = sigmaU_pix2;
      sigmaV2 = sigmaV_pix2;

    } 
    else if(vol_id < 16) {
      sigmaU2 = sigmaU_shs2;
      sigmaV2 = sigmaV_shs2;
    } else {
      sigmaU2 = sigmaU_lgs2;
      sigmaV2 = sigmaV_lgs2;
    }

    float Phi = atan2(imp.m_y, imp.m_x); 
    float Rz = 0.0;
    if(pL->m_info.m_type == 0) {//barrel
      Rz = imp.m_z;
    } else {
      Rz = sqrt(imp.m_x*imp.m_x + imp.m_y*imp.m_y);
    }

    float refCoord = pL->m_info.m_refCoord;

    float cm, cp;
    
    if(pL->m_info.m_type == 0) {//barrel

      cm = (Rz - delta_z_m)/refCoord;
      cp = (Rz + delta_z_p)/refCoord;
    }
    else {
      cm = refCoord/(Rz - delta_R_m);
      cp = refCoord/(Rz + delta_R_p);
    }

    if(cm > cp) std::swap(cm, cp);
    
    int bin_start, bin_end;

    bool valid = pL->getBinRange(cm, cp, bin_start, bin_end);
 
    if(!valid) continue;
    
    if(bin_end<bin_start) {
      std::cout<<"ERROR: collectHits bin range : "<<bin_start<<" "<<bin_end<<" z_range :"<<cm<<" "<<cp<<std::endl;
      continue;
    }
    //loop over bins in range

    std::vector<int> vTmp;
 
    for(int bin = bin_start; bin<=bin_end;bin++) {

      HIT_BIN& B = pL->m_bins[bin];

      //find compatible hits
	
      std::vector<NODE*> vN;

      B.selectByPhi(Phi, phiWidth, vN);

      if(vN.empty()) continue;
	
      for(std::vector<NODE*>::iterator nIt = vN.begin();nIt!=vN.end();++nIt) {

	for(std::vector<int>::iterator hIt=(*nIt)->m_hits.begin();hIt!=(*nIt)->m_hits.end();++hIt) {

	  int h = (*hIt);
      
	  int nCellsT = hits->m_clusterHeight[h];

	  if(nCellsT > 4) continue;

	  vTmp.push_back(*hIt);
	}
      }
    }

    std::vector<int> modVec;
    std::vector<int> hitArray[100];

    for(std::vector<int>::iterator hIt=vTmp.begin();hIt!=vTmp.end();++hIt) {

      int h = (*hIt);
	
      int mod_id = hits->m_mod_id[h];

      int hitCollIdx = -1;
	
      for(int k=0;k<(int)modVec.size();k++) {
	if(modVec[k] == mod_id) {
	  hitCollIdx = k;break;
	}
      }

      if(hitCollIdx == -1) {
	if(modVec.size()<99) {
	  modVec.push_back(mod_id);
	}
	hitCollIdx = (int)modVec.size() - 1;
	hitArray[hitCollIdx].reserve(10);
      }
      hitArray[hitCollIdx].push_back(h);
    }

    for(int modIdx=0;modIdx<(int)modVec.size();modIdx++)  {

      int mod_id = modVec[modIdx];
	
      const DETECTOR_ELEMENT* pDE = m_geo.getDetectorElement(vol_id, lay_id, mod_id);
      
      if(pDE==0) {
	std::cout<<"NOT FOUND "<<vol_id<<" "<<lay_id<<" "<<mod_id<<std::endl;
	continue;
      }
      
      float lP[3];
      float gP[3] = {imp.m_x, imp.m_y, imp.m_z};
      
      pDE->toLocal(gP, lP);
      
      if(!pDE->isInSide(lP, tol)) continue;
      
      std::vector<int>& vH = hitArray[modIdx];
      
      int bestIdx=-1;
      float minDChi2=1e9;
      
      for(unsigned int idx=0;idx<vH.size();idx++) {
	
	int h = vH[idx];
	float r[2] = {hits->m_m1[h]-lP[0], hits->m_m2[h]-lP[1]};
	
	if(!pDE->m_isRect) {
	  
	  float sc[2], d[2];
	  pDE->getResidualTransform(hits->m_m1[h], hits->m_m2[h], sc);
	  
	  d[0] = r[0]*sc[1] - r[1]*sc[0];
	  d[1] = r[0]*sc[0] + r[1]*sc[1];
	  r[0] = d[0];r[1] = d[1];
	  
	}
	
	float dchi2 = r[0]*r[0]*sigmaU2 + r[1]*r[1]*sigmaV2;
	
	if(dchi2 > chi2Cut) continue;
	
	if(dchi2 < minDChi2) {
	  minDChi2 = dchi2;
	  bestIdx = h;
	}
      }
      
      if(bestIdx!=-1 && minDChi2<chi2Cut) {
	vHits.push_back(bestIdx);
      }      
    }
  }
}

void HitBank::collectHits(const std::vector<int>& layVec, const TRAJECTORY& vP, std::vector<int>& vH, const HITS* hits) {

  float delta_Rz_m = 19.0;//was 18.0
  float delta_Rz_p = 19.0;
  float valGate_pix = 2.0;//was 2.0
  float valGate_shs = 18.0;//was 18.0
  float valGate_lgs = 22.0;//was 22.0
  float phiGate_pix = 0.032;//was 0.041
  float phiGate_shs = 0.064;//was 0.064
  float phiGate_lgs = 0.055;//was 0.054

  vH.clear();

  if(vP.m_nPoints<2) return;

  int it1 = 0;
  int it2 = 1;

  TRAJECTORY vOut;
  
  for(;it2<vP.m_nPoints;it1++,it2++) {
    m_geo.getLayerIntersection(vP.m_p[it1], vP.m_p[it2], vOut);
    if(vOut.m_nPoints>0) {
      unsigned int layerKey = vOut.m_p[vOut.m_nPoints-1].m_key;
      if(layerKey == 16002 || layerKey == 17004 || layerKey == 18012 || layerKey == 14012 || layerKey == 12002) break;
    }
    
  }

  for(int impactIdx=0;impactIdx<vOut.m_nPoints;impactIdx++) {

    unsigned int key = vOut.m_p[impactIdx].m_key;

    if(std::find(layVec.begin(), layVec.end(), (int)key) != layVec.end()) continue;
    
    POINT_3D& imp = vOut.m_p[impactIdx];

    HIT_LAYER* pL = getLayer(key);
    unsigned int vol_id = key / 1000;

    float gateWidth = 0.0;
    float phiGateWidth = 0.0;
    
    if(vol_id < 10) {//pixel
      gateWidth = valGate_pix;
      phiGateWidth = phiGate_pix;
    } 
    else if(vol_id < 16) {
      gateWidth = valGate_shs;
      phiGateWidth = phiGate_shs;
    } else {
      gateWidth = valGate_lgs;
      phiGateWidth = phiGate_lgs;
    }

    float phiWidth = phiGateWidth;

    
    float Rz = 0.0;
    if(pL->m_info.m_type == 0) {//barrel
      Rz = imp.m_z;
    } else {
      Rz = sqrt(imp.m_x*imp.m_x + imp.m_y*imp.m_y);
    }

    float refCoord = pL->m_info.m_refCoord;

    float cm, cp;
    
    if(pL->m_info.m_type == 0) {//barrel

      cm = (Rz - delta_Rz_m)/refCoord;
      cp = (Rz + delta_Rz_p)/refCoord;
    }
    else {

      cm = refCoord/(Rz - delta_Rz_m);
      cp = refCoord/(Rz + delta_Rz_p);

    }

    if(cm > cp) std::swap(cm, cp);
    
    int bin_start, bin_end;

    bool valid = pL->getBinRange(cm, cp, bin_start, bin_end);
 
    if(!valid) continue;
    
    if(bin_end<bin_start) {
      std::cout<<"ERROR: collectHits bin range : "<<bin_start<<" "<<bin_end<<" tau_range :"<<cm<<" "<<cp<<std::endl;
      continue;
    }
    //loop over bins in range

    std::vector<int> vTmp;

    float Phi = atan2(imp.m_y, imp.m_x);
 
    for(int bin = bin_start; bin<=bin_end;bin++) {

      HIT_BIN& B = pL->m_bins[bin];

	//find compatible hits
	
	std::vector<NODE*> vN;
	B.selectByPhi(Phi, phiWidth, vN);

	if(vN.empty()) continue;
	
	for(std::vector<NODE*>::iterator nIt = vN.begin();nIt!=vN.end();++nIt) {
	  std::copy((*nIt)->m_hits.begin(), (*nIt)->m_hits.end(), std::back_inserter(vTmp));
	}
    }

    //copy with validation
    if(pL->m_info.m_type == 0) {//barrel
	  
      for(std::vector<int>::iterator hIt=vTmp.begin();hIt!=vTmp.end();++hIt) {
	    
	if(fabs(Rz-hits->m_z[*hIt]) > gateWidth) continue;
	
	vH.push_back(*hIt);
	
      }
    }
    else {

      for(std::vector<int>::iterator hIt=vTmp.begin();hIt!=vTmp.end();++hIt) {
	
	float r = hits->m_r[*hIt];
	if(fabs(Rz-r) > gateWidth) continue;
	
	vH.push_back(*hIt);
	
      }
    }
  }
}

void HitBank::collectHitsInward(const TRAJECTORY& vP, std::vector<int>& vH) {

  const float delta_Rz_m_2 = 2.15;//was 2.2
  const float delta_Rz_p_2 = 2.15;

  const float phiWidth_2 = 0.068;//was 0.068

  const float delta_Rz_m_4 = 3.9;
  const float delta_Rz_p_4 = 3.9;

  const float phiWidth_4 = 0.030;//was 0.030

  float valGate_pix = 2.7;//was 2.7
  float valGate_shs = 2.5;//was 2.5

  float phiGate_pix = 0.034;//was 0.034
  float phiGate_shs = 0.020;//was 0.020

  vH.clear();

  if(vP.m_nPoints<2) return;

  int it1 = vP.m_nPointsOnTrack;
  int it2 = it1 + 1;

  TRAJECTORY vOut;
  
  for(;it2<vP.m_nPoints;it1++,it2++) m_geo.getLayerIntersection(vP.m_p[it1], vP.m_p[it2], vOut);

  for(int impactIdx=0;impactIdx<vOut.m_nPoints;impactIdx++) {

    unsigned int key = vOut.m_p[impactIdx].m_key;

    POINT_3D& imp = vOut.m_p[impactIdx];

    unsigned int vol_id = key / 1000;
    unsigned int lay_id = key % 1000;

    float delta_Rz_m = 0.0;
    float delta_Rz_p = 0.0;

    float phiWidth = 0.0;

    if(vol_id == 8 && lay_id == 2) {
      delta_Rz_m = delta_Rz_m_2;
      delta_Rz_p = delta_Rz_p_2;
      phiWidth = phiWidth_2;
    }
    else if(vol_id == 8 && lay_id == 4) {
      delta_Rz_m = delta_Rz_m_4;
      delta_Rz_p = delta_Rz_p_4;
      phiWidth = phiWidth_4;
    } else {

      if(vol_id < 10) {//pixel
	delta_Rz_m = valGate_pix;
	delta_Rz_p = valGate_pix;
	phiWidth = phiGate_pix;
      }
      else if(vol_id < 16) {
	delta_Rz_m = valGate_shs;
	delta_Rz_p = valGate_shs;
	phiWidth = phiGate_shs;
      }

    }
    
    //std::cout<<"volume = "<<vol_id<<" layer = "<<lay_id<<std::endl;
    
    HIT_LAYER* pL = getLayer(key);

    //std::cout<<"nHits="<<pL->numberOfHits()<<" layer type "<<pL->m_info.m_type<<std::endl;

    float Phi = atan2(imp.m_y, imp.m_x);
 
    float Rz = 0.0;
    if(pL->m_info.m_type == 0) {//barrel
      Rz = imp.m_z;
    } else {
      Rz = sqrt(imp.m_x*imp.m_x + imp.m_y*imp.m_y);
    }

    if(std::isnan(Rz)) continue;

    float refCoord = pL->m_info.m_refCoord;

    float cm, cp;
    
    if(pL->m_info.m_type == 0) {//barrel

      cm = (Rz - delta_Rz_m)/refCoord;
      cp = (Rz + delta_Rz_p)/refCoord;
    }
    else {
      cm = refCoord/(Rz - delta_Rz_m);
      cp = refCoord/(Rz + delta_Rz_p);
    }

    if(cm > cp) std::swap(cm, cp);

    int bin_start, bin_end;

    bool valid = pL->getBinRange(cm, cp, bin_start, bin_end);
 
    if(!valid) continue;

    if(bin_end<bin_start) {
      std::cout<<"ERROR: collectHitsInward bin range : "<<bin_start<<" "<<bin_end<<" z_range :"<<cm<<" "<<cp<<std::endl;
      continue;
    }
    //loop over bins in range

    for(int bin = bin_start; bin<=bin_end;bin++) {

      HIT_BIN& B = pL->m_bins[bin];

      //find compatible hits
	
      std::vector<NODE*> vN;

      B.selectByPhi(Phi, phiWidth, vN);

      if(vN.empty()) continue;
	
      for(std::vector<NODE*>::iterator nIt = vN.begin();nIt!=vN.end();++nIt) {
	
	std::copy((*nIt)->m_hits.begin(), (*nIt)->m_hits.end(), std::back_inserter(vH));
      }
    }
  }
  
}

