#include<iostream>
#include<cmath>
#include<algorithm>

#include "DataStructures.h"

void HitBin::clear() {
  for(std::vector<NODE*>::iterator it = m_nodes.begin();it!=m_nodes.end();++it) {
    delete (*it);
  }
  m_nodes.clear();
}

void HitBin::deleteUsedNodes(const HITS* hits) {
  
  std::vector<NODE*> vTmp;

  for(std::vector<NODE*>::iterator it=m_nodes.begin();it!= m_nodes.end();++it) {

    if((*it)->m_hits.empty()) {
      if(hits->m_mask[(*it)->m_h] != 0) {
	delete (*it);
	continue;
      }
    }
    else {

      unsigned int nMasked = 0;
      
      for(std::vector<int>::const_iterator hit = (*it)->m_hits.begin();hit!=(*it)->m_hits.end();++hit) {
	if(hits->m_mask[*hit] != 0) nMasked++;
      }

      if(nMasked == (*it)->m_hits.size()) {
	delete (*it);
	continue;
      }
      
      if(nMasked != 0) {
	
	std::vector<int> vClear;
      
	for(std::vector<int>::const_iterator hit = (*it)->m_hits.begin();hit!=(*it)->m_hits.end();++hit) {
	  if(hits->m_mask[*hit] == 0) vClear.push_back(*hit); 
	}
	if(vClear.empty()) {
	  delete (*it);
	  continue;
	}
      
	(*it)->m_hits = std::move(vClear);
      
	(*it)->m_h = (*(*it)->m_hits.begin());
      
	(*it)->updateCoordinates(hits);
	(*it)->calculateCellShape(hits);
      }
    }

    unsigned int s1 = (*it)->m_in.size();
    (*it)->m_in.clear();
    (*it)->m_in.reserve(s1);
    unsigned int s2 = (*it)->m_out.size();
    (*it)->m_out.clear();
    (*it)->m_out.reserve(s2);
    
    vTmp.push_back(*it);
  }

  m_nodes = std::move(vTmp);
}

void HitBin::selectByPhi(float phi, float width, std::vector<NODE*>& vn) {

  vn.clear();

  float phi_min = phi-width;
  float phi_max = phi+width;

  for(std::vector<PHI_NODE>::iterator it=m_phiIndices.begin();it!=m_phiIndices.end();++it) {
    
    float phi = (*it).m_phi;
    if(phi<phi_min) continue;
    if(phi<=phi_max) {
      vn.push_back((*it).m_pN);
    }
    if(phi>phi_max) break;

  }

}

void HitBin::searchPhiRange(float phi_min, float phi_max, int& j0, int& j1) {
  while(j0<(int)m_phiIndices.size() && m_phiIndices[j0].m_phi < phi_min) j0++;
  while(j1<(int)m_phiIndices.size() && m_phiIndices[j1].m_phi < phi_max) j1++;
}




void HitBin::removeEmptyNodes() {

  std::vector<NODE*> vTmp;

  for(std::vector<NODE*>::iterator it=m_nodes.begin();it!= m_nodes.end();++it) {
    if((*it)->m_h == 0) {
      delete (*it);
      continue;
    }
    vTmp.push_back(*it);
  }

  m_nodes = std::move(vTmp);
}

void HitBin::updateCoordinates(const HITS* hits) {

  for(unsigned int nIdx=0;nIdx<m_nodes.size();nIdx++) {
    NODE* pN = m_nodes.at(nIdx);
    if(pN->m_hits.size() == 1) continue;
    pN->updateCoordinates(hits);
  }
 
}

void HitBin::calculateCellShape(const HITS* hits) {

  for(unsigned int nIdx=0;nIdx<m_nodes.size();nIdx++) {
    NODE* pN = m_nodes.at(nIdx);
    if(pN->m_hits.size() == 1) continue;
    pN->calculateCellShape(hits);
  }
 
}

int HitBin::checkModset(const HITS* hits) {

  int nErrors = 0;

  for(unsigned int nIdx=0;nIdx<m_nodes.size();nIdx++) {
    NODE* pN = m_nodes.at(nIdx);
    if(pN->m_hits.size() == 1) continue;
    nErrors += pN->checkModset(hits);
  }
  return nErrors;
}

void HitBin::fillIndices(bool isPixeEndcap) {

  const float phiMargin = 0.23;
  const int T2Cut = 3;
  
  float two_pi = 2*M_PI;

  m_phiIndices.clear();

  for(unsigned int nIdx=0;nIdx<m_nodes.size();nIdx++) {
    NODE* pN = m_nodes[nIdx];
    float phi = pN->m_phi;
    if(phi > M_PI-phiMargin) {
      m_phiIndices.push_back(PHI_NODE(phi-two_pi,pN));
    }
  }
  for(unsigned int nIdx=0;nIdx<m_nodes.size();nIdx++) {
    NODE* pN = m_nodes[nIdx];
    m_phiIndices.push_back(PHI_NODE(pN->m_phi,pN)); 
  }

  for(unsigned int nIdx=0;nIdx<m_nodes.size();nIdx++) {
    NODE* pN = m_nodes[nIdx];
    float phi = pN->m_phi;
    if(phi < -M_PI+phiMargin) {
      m_phiIndices.push_back(PHI_NODE(phi+two_pi,pN)); 
    } else break;
  }

  m_l2phiIndices.clear();
  if(isPixeEndcap) {
    for(unsigned int nIdx=0;nIdx<m_nodes.size();nIdx++) {
      NODE* pN = m_nodes[nIdx];
      int nCellsT2 = pN->m_nCellsT;
      if(nCellsT2 > T2Cut) continue;
      if(pN->m_nCells > 7) continue;
      float phi = pN->m_phi;
      if(phi > M_PI-phiMargin) {
	m_l2phiIndices.push_back(PHI_NODE(phi-two_pi,pN)); 
      }
    }
    for(unsigned int nIdx=0;nIdx<m_nodes.size();nIdx++) {
      NODE* pN = m_nodes[nIdx];
      int nCellsT2 = pN->m_nCellsT;
      if(nCellsT2 > T2Cut) continue;
      if(pN->m_nCells > 7) continue;
      m_l2phiIndices.push_back(PHI_NODE(pN->m_phi,pN)); 
    }
    
    for(unsigned int nIdx=0;nIdx<m_nodes.size();nIdx++) {
      NODE* pN = m_nodes[nIdx];
      int nCellsT2 = pN->m_nCellsT;
      if(nCellsT2 > T2Cut) continue;
      if(pN->m_nCells > 7) continue;
      float phi = pN->m_phi;
      if(phi < -M_PI+phiMargin) {
	m_l2phiIndices.push_back(PHI_NODE(phi+two_pi,pN)); 
      } else break;
    }
  }
  else {
    for(unsigned int nIdx=0;nIdx<m_nodes.size();nIdx++) {
      NODE* pN = m_nodes[nIdx];
      int nCellsT2 = pN->m_nCellsT;
      if(nCellsT2 > T2Cut) continue;
      float phi = pN->m_phi;
      if(phi > M_PI-phiMargin) {
	m_l2phiIndices.push_back(PHI_NODE(phi-two_pi,pN)); 
      }
    }
    for(unsigned int nIdx=0;nIdx<m_nodes.size();nIdx++) {
      NODE* pN = m_nodes[nIdx];
      int nCellsT2 = pN->m_nCellsT;
      if(nCellsT2 > T2Cut) continue;
      m_l2phiIndices.push_back(PHI_NODE(pN->m_phi,pN)); 
    }
    
    for(unsigned int nIdx=0;nIdx<m_nodes.size();nIdx++) {
      NODE* pN = m_nodes[nIdx];
      int nCellsT2 = pN->m_nCellsT;
      if(nCellsT2 > T2Cut) continue;
      float phi = pN->m_phi;
      if(phi < -M_PI+phiMargin) {
	m_l2phiIndices.push_back(PHI_NODE(phi+two_pi,pN)); 
      } else break;
    }
  }
  
}
