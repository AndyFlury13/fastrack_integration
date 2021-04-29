#include<iostream>
#include<cstring>
#include<cmath>
#include <math.h>

#include "FastTrackFitter.h"

#include<algorithm>

float FastTrackFitter::calculateLikelihood(float w, int nHits, int* hArray, HIT_INFO* info) const {

  float chi2_sum = 0;

  unsigned int vl[256];
  unsigned int layerKeys[128];
  
  memset(&vl[0],0,sizeof(vl));

  int nL = 0;
  
  for(int hitIdx = 0;hitIdx<nHits;hitIdx++) {

    if(info[hitIdx].m_status != 0) continue;

    int h = hArray[hitIdx];

    unsigned int layerKey = 7*(m_hits->m_vol_id[h]-7) + (m_hits->m_lay_id[h]>>1);

    if(vl[layerKey] == 0) nL++;
    vl[layerKey]++;
    layerKeys[hitIdx] = layerKey;
  }
  
  for(int hitIdx=0;hitIdx<nHits;hitIdx++) {
    
    if(info[hitIdx].m_status != 0) continue;

    int nc = vl[layerKeys[hitIdx]];

    chi2_sum += info[hitIdx].m_dchi2/nc;
  }

  float L = w*nL - chi2_sum;
  
  return L;

}

int FastTrackFitter::getTriplet(NEW_TRACK* pT, int* v) const {

  int it = 0;

  int h = pT->m_hits[it];
  
  v[0] = h;

  ++it;

  int hitCounter = 1;
  
  while(it < pT->m_nHits) {

    int vol_id = m_hits->m_vol_id[h];
    int lay_id = m_hits->m_lay_id[h];

    int h2 = pT->m_hits[it];

    int vol_id_2 = m_hits->m_vol_id[h2];
    int lay_id_2 = m_hits->m_lay_id[h2];
    
    if(vol_id == vol_id_2 && lay_id == lay_id_2) {
      ++it;continue;
    }
    h = h2;
    v[hitCounter++] = h;
    if(hitCounter == 3) break;
    ++it;
  }
  return hitCounter;
}

int FastTrackFitter::getLastTriplet(NEW_TRACK* pT, int* v) const {
  
  int it = pT->m_nHits-1;

  int h = pT->m_hits[it];
  v[0] = h;

  --it;

  int hitCounter = 0;

  while(it >=0) {
    int vol_id = m_hits->m_vol_id[h];
    int lay_id = m_hits->m_lay_id[h];
    
    int h2 = pT->m_hits[it];

    int vol_id_2 = m_hits->m_vol_id[h2];
    int lay_id_2 = m_hits->m_lay_id[h2];
    
    if(vol_id == vol_id_2 && lay_id == lay_id_2) {
      --it;continue;
    }
    h = h2;
    v[hitCounter++] = h;
    if(hitCounter == 3) break;
    --it;
  }
  return hitCounter;
}

bool FastTrackFitter::initialize(NEW_TRACK* pT) const {

  const double Cs = 0.000299997;
  
  //select 3 hits

  int vh3[3];

  int n3 = getTriplet(pT, vh3);

  if(n3<3) return false;

  int h1 = vh3[0];
  int h2 = vh3[1];
  int h3 = vh3[2];

  float x0 = m_hits->m_x[h1];
  float y0 = m_hits->m_y[h1];
  float z0 = m_hits->m_z[h1];

  float dx1 = m_hits->m_x[h2]-x0;
  float dy1 = m_hits->m_y[h2]-y0;
  
  float dx2 = m_hits->m_x[h3]-x0;
  float dy2 = m_hits->m_y[h3]-y0;

  float dz2 = m_hits->m_z[h3]-z0; 

  //conformal mapping

  double rn = dx2*dx2+dy2*dy2;
  double r2 = 1/rn;
  double u1 = 1/sqrt(dx1*dx1+dy1*dy1);

  double K = u1*r2;

  double a  = dx1*K;
  double b  = dy1*K;

  double u2 = a*dx2 + b*dy2;
  double v2 = a*dy2 - b*dx2;

  double A  = v2/(u2-u1);
  double B  = v2-A*u2;

  double Curv  = B/sqrt(1+A*A);//curvature

  float ctg  = dz2*sqrt(r2);//cot(theta)
  
  const DETECTOR_ELEMENT* pDE = m_geo.getDetectorElement(m_hits->m_vol_id[h1],
							 m_hits->m_lay_id[h1],
							 m_hits->m_mod_id[h1]);

  float Corr = m_params.tf_f[0]; //correction factor due to the non-uniform magnetic field

  float P[5];

  P[0] = m_hits->m_m1[h1];
  P[1] = m_hits->m_m2[h1];

  P[2] = atan2(b+a*A, a-b*A);//phi0
  P[3] = atan(1/ctg);//theta

  if(P[3]<0.0) P[3] += M_PI;

  float fC = Cs*Corr;

  float pt_inv = -Curv/fC;// 1/pT in GeV

  P[4] = pt_inv/sqrt(1+ctg*ctg);//full inverse momentum

  pT->m_endTrackState.initialize(P, pDE);
  
  return true;

}

bool FastTrackFitter::fitTrackCandidate(NEW_TRACK* pT, int stage) const {

  float maxDChi2 = 27.0;
  float maxFrac  = 0.10;
  float llWeight = 5.0;

  if(stage<1 || stage>5) return false;
  
  m_stage = stage;

  int offset = 3*stage-3;
  
  maxDChi2 = m_params.tf_f[1+offset];
  maxFrac  = m_params.tf_f[2+offset];
  llWeight = m_params.tf_f[3+offset];
  
  //1. create track, sort hits

  pT->sortHits(m_hits);

  if(!initialize(pT)) return false;
  
  TRACK_STATE* pTS = &pT->m_endTrackState;

  bool result = true;
  pT->m_status = 0;
  int nOutliers=0;

  int nHits0 = pT->m_nHits;
  float OutlThr = maxFrac*nHits0;

  int hit_status[100];
  
  for(int hitIdx=0;hitIdx<nHits0;hitIdx++) {
    
    int h = pT->m_hits[hitIdx];

    int code = m_fte.fastExtrapolator(pTS, h, false);

    if(code == -20) {//not in bounds
      hit_status[hitIdx] = -1;
      nOutliers++;
    
      if(nOutliers >= OutlThr) {
	result = false;
	break;
      }
      memcpy(&pTS->m_Rk[0], &pTS->m_Re[0], sizeof(pTS->m_Re));
      memcpy(&pTS->m_Ck[0], &pTS->m_Ce[0], sizeof(pTS->m_Ce));
      continue;
    }

    if(code!=0) {
      result = false;
      break;
    }

    FIT_HIT fh(h);
    
    update(pTS, &fh, maxDChi2);

    if(!pTS->valid()) {
      result = false;
      break;
    }
 
    hit_status[hitIdx] = fh.m_status;
    
    if(fh.m_status!=0) nOutliers++;
    
    if(nOutliers >= OutlThr) {
      result = false;
      break;
    }
    
  }
  
  if(result) {
    
    if(nOutliers>0 && nOutliers < OutlThr) {
      pT->m_nHits = 0;
      for(int hitIdx=0;hitIdx<nHits0;hitIdx++) {	
	if(hit_status[hitIdx] == 0) {
	  pT->m_hits[pT->m_nHits++] = pT->m_hits[hitIdx];
	}
      }
      pT->update(m_hits);
    }    
  }
  
  if(result) {

    float L = llWeight*pT->m_nlayers - pTS->m_chi2;
    
    pT->m_ll = L;

    pT->m_chi2 = pTS->m_chi2;
    pT->m_ndof = pTS->m_ndof;
    pT->m_status = 1;
  }
  else {
    pT->m_endTrackState.m_initialized = false;

  }
  
  if(!result) return false;

  return true;
  
}


void FastTrackFitter::fitTrack(NEW_TRACK* pT, int stage) {

  float maxDChi2 = 100.0;
  float llWeight = 60.0;
  float fracOutl = 0.33;

  if( !((stage>=1 && stage <= 4) || (stage>=11 && stage <= 14))) {
    pT->m_status = 0;
    return;
  }
    
  int offset = stage < 11 ? 16+3*(stage-1) : 28+3*(stage-11);

  maxDChi2 = m_params.tf_f[offset];
  llWeight = m_params.tf_f[offset+1];
  fracOutl = m_params.tf_f[offset+2];
  
  //1. create track, sort hits
  
  pT->sortHits(m_hits);

  pT->m_status = 0;
  
  if(pT->m_beginTrackState.valid()) {
    pT->m_endTrackState.initialize(pT->m_beginTrackState);
  }
  else {
    if(!initialize(pT)) return;
  }
 
  TRACK_STATE* pTS = &pT->m_endTrackState;

  bool result = true;

  pT->m_status = 0;
  
  int nOutliers=0;

  int nHits0 = pT->m_nHits;

  HIT_INFO info[MAX_HIT_ON_TRACK];

  for(int hitIdx=0;hitIdx<nHits0;hitIdx++) {
    
    int h = pT->m_hits[hitIdx];

    int code = m_fte.fastExtrapolator(pTS, h, false);

    if(code == -20) {//not in bounds
      nOutliers++;
      info[hitIdx].m_dchi2 = 0.0;
      info[hitIdx].m_status = -1;
      memcpy(&pTS->m_Rk[0], &pTS->m_Re[0], sizeof(pTS->m_Re));
      memcpy(&pTS->m_Ck[0], &pTS->m_Ce[0], sizeof(pTS->m_Ce));
      continue;
    }

    if(code!=0) {
      result = false;
      //std::cout<<"extrapolation failure, code = "<<code<<std::endl;
      break;
    }

    FIT_HIT fh(h);

    update(pTS, &fh, maxDChi2);

    if(!pTS->valid()) {
      result = false;
      break;
    }

    if(fh.m_status!=0 || std::isnan(fh.m_dchi2) ) nOutliers++;
    
    info[hitIdx].m_dchi2 = fh.m_dchi2;
    info[hitIdx].m_status = fh.m_status;
    
  }
  
  if(result) {

    //pTS->print();
    
    //update track here ....
    
    if(nOutliers>0 && nOutliers<fracOutl*nHits0) {

      std::vector<int> newHits;
      std::vector<HIT_INFO> newInfo;
      
      for(int hitIdx=0;hitIdx<nHits0;hitIdx++) {
	
	if(info[hitIdx].m_status != 0) continue;

	newInfo.push_back(info[hitIdx]);
	newHits.push_back(pT->m_hits[hitIdx]);
      }

      pT->m_nHits=0;
      for(unsigned int k=0;k<newHits.size();k++) {
	info[pT->m_nHits] = newInfo[k];
      	pT->m_hits[pT->m_nHits++] = newHits[k];
      }

      pT->update(m_hits);
    }

    pT->m_ll = calculateLikelihood(llWeight, pT->m_nHits, pT->m_hits, info);

    pT->m_status = 1;
  }
  
}

bool FastTrackFitter::getTrajectory(NEW_TRACK* pT, TRAJECTORY& vPoints) {

  float maxDChi2 = m_params.tf_f[40];
  float minStep  = m_params.tf_f[41];
  float maxFrac  = m_params.tf_f[42];
  
  bool result = true;

  vPoints.m_nPoints = 0;

  if(!pT->m_endTrackState.valid()) {
    return false;
  }

  int nHits0 = pT->m_nHits;

  if(nHits0 == 0) return false;
  
  pT->sortHits(m_hits);

  pT->m_beginTrackState.initialize(pT->m_endTrackState); 
  
  TRACK_STATE* pInit = &pT->m_beginTrackState;
    
  pInit->resetCovariance();

  TRACK_STATE* pTS = pInit;

  int nOutliers=0;
  
  float maxOutl = maxFrac*nHits0;

  int kIdx = nHits0-1;

  const DETECTOR_ELEMENT* pDE0 = pT->m_endTrackState.m_pDE;
  
  for(;kIdx>0;kIdx--) {
    int h = pT->m_hits[kIdx];
    const DETECTOR_ELEMENT* pNE = m_geo.getDetectorElement(m_hits->m_vol_id[h], m_hits->m_lay_id[h], m_hits->m_mod_id[h]);
    if(pNE == pDE0) break;
  }

  if(kIdx == 0) kIdx = nHits0-1;

  FIT_HIT fh(pT->m_hits[kIdx]);

  update(pTS, &fh, maxDChi2);
  
  if(fh.m_status!=0) nOutliers++;

  kIdx--;
    
  for(;kIdx>=0;kIdx--) {
    
    int h = pT->m_hits[kIdx];

    int code = m_fte.fastExtrapolate(pTS, h, minStep, vPoints);
        
    if(code!=0) {
      result = false;
      break;
    }

    FIT_HIT fh(h);
    
    update(pTS, &fh, maxDChi2);

    if(!pTS->valid()) {
      result = false;
      break;
    }
    
    if(fh.m_status!=0) nOutliers++;

    if(nOutliers >= maxOutl) {
      result = false;
      break;
    }

  }

  if(result) {
    
    int firstHit = pT->m_hits[0];

    bool searchInside = true;
    
    if(m_hits->m_vol_id[firstHit] == 8 && m_hits->m_lay_id[firstHit] == 2) searchInside = false;
     
    if(searchInside) m_fte.extrapolateInsideFast(pTS, vPoints);
      
  }

  return result;
}

void FastTrackFitter::backwardPropagation(NEW_TRACK* pT, TRAJECTORY& vPoints) {

  float maxDChi2 = m_params.tf_f[43];
  float maxFrac  = m_params.tf_f[44];
  
  //1. create track state
  
  if(!pT->m_endTrackState.valid()) {
    return;
  }

  int nHits0 = pT->m_nHits;

  if(nHits0 == 0) return;

  pT->m_beginTrackState.initialize(pT->m_endTrackState);
  
  TRACK_STATE* pInit = &pT->m_beginTrackState;
    
  pInit->resetCovariance();

  TRACK_STATE* pTS = pInit;

  bool result = true;
  
  int nOutliers=0;

  float maxOutl = maxFrac*nHits0;
  
  int kIdx = nHits0-1;
  const DETECTOR_ELEMENT* pDE0 = pT->m_endTrackState.m_pDE;
  for(;kIdx>0;kIdx--) {
    int h = pT->m_hits[kIdx];
    const DETECTOR_ELEMENT* pNE = m_geo.getDetectorElement(m_hits->m_vol_id[h], m_hits->m_lay_id[h], m_hits->m_mod_id[h]);
    if(pNE == pDE0) break;
  }

  if(kIdx == 0) kIdx = nHits0-1;
 
  FIT_HIT fh(pT->m_hits[kIdx]);
  update(pTS, &fh, maxDChi2);
  
  if(fh.m_status!=0) nOutliers++;

  kIdx--;
    
  for(;kIdx>=0;kIdx--) {
    
    int h = pT->m_hits[kIdx];

    int code = m_fte.fastExtrapolate(pTS, h, 0.0, vPoints);
    
    if(code!=0) {
      result = false;
      //std::cout<<"extrapolation failure, code = "<<code<<std::endl;
      break;
    }

    FIT_HIT fh(h);
    
    update(pTS, &fh, maxDChi2);
    if(!pTS->valid()) {
      result = false;
      break;
    }
    
    if(fh.m_status!=0) nOutliers++;

    if(nOutliers >= maxOutl) {
      result = false;
      break;
    }
  }
  
  vPoints.m_nPointsOnTrack = vPoints.m_nPoints;

  if(result) {

    int firstHit = pT->m_hits[0];

    bool searchInside = true;
    
    if(m_hits->m_vol_id[firstHit] == 8 && m_hits->m_lay_id[firstHit] == 2) searchInside = false;
     
    if(searchInside) m_fte.extrapolateInsideFast(pTS, vPoints);

  }

}


void FastTrackFitter::addHitsBackward(NEW_TRACK* pT, std::vector<int>& vHits) {

  const int nMaxNoAdd = 3;
  
  float dQ_hole      = m_params.tf_f[45];
  float dQ_hit       = m_params.tf_f[46];
  float minQ         = m_params.tf_f[47];
  float maxDChi2_upd = m_params.tf_f[48];
  int nBest          = m_params.tf_i[0];
  
  int globalBranchCounter = 0;

  if(pT->m_nHits >= MAX_HIT_ON_TRACK) return;

  if(!pT->m_beginTrackState.valid()) {
    std::cout<<"No valid begin track state found, skipping the extension refit !"<<std::endl;
    for(std::vector<int>::iterator it=vHits.begin(); it!=vHits.end();it++) {
      pT->m_hits[pT->m_nHits++] = (*it);
      if(pT->m_nHits >= MAX_HIT_ON_TRACK) break;
    }
    return;
  }

  //do branching exploration here ...

  //arrange hits in layers

  std::vector<int> hitArray[100];

  std::vector<unsigned int> keyArray;

  for(std::vector<int>::iterator hIt=vHits.begin();hIt!=vHits.end();++hIt) {

    int h = (*hIt);

    if(m_hits->m_mask[h] == 2) continue;

    int vol_id = m_hits->m_vol_id[h];
    int lay_id = m_hits->m_lay_id[h];

    unsigned int layerKey = (vol_id-7) + (lay_id>>1);

    if(std::find(keyArray.begin(), keyArray.end(), layerKey) == keyArray.end()) {
      keyArray.push_back(layerKey);
    }
    hitArray[layerKey].push_back(h);
  }


  //std::cout<<"Processing new track extension ....."<<std::endl;
  
  std::vector<TRACK_CLONE*> vBranches;

  TRACK_CLONE* pClone = &m_tClones[globalBranchCounter++];

  pClone->initialize(&pT->m_beginTrackState);
  vBranches.push_back(pClone);
  
  std::vector<TRACK_CLONE*> newBranches;

  int lIdx=0;

  bool errorFlag = false;

  int nNoAdd = 0;
  
  for(;lIdx<(int)keyArray.size();lIdx++) {

    std::vector<int>& vHL = hitArray[keyArray[lIdx]];
      
    unsigned int nHits = vHL.size();

    int nNewTotal=0;
    
    for(std::vector<TRACK_CLONE*>::iterator bIt=vBranches.begin();bIt!=vBranches.end();++bIt) {

      int bCounter = 0;

      //1. add "skip-all", "no-detect-in-this layer" branch

      TRACK_CLONE* nodBr = &m_tClones[globalBranchCounter++];

      if(globalBranchCounter >= MAX_CLONE) {
	errorFlag = true;
	break;
      }
      
      nodBr->clone(**bIt);
	
      nNewTotal++;
      nodBr->addHole(dQ_hole);
      
      //std::cout<<"added hole, Q="<<nodBr->m_Q<<std::endl;

      newBranches.push_back(nodBr);bCounter++;

      nNoAdd++;
      
      //2. investigate hits

      //arrange hits into 1- or 2-hit clusters -> each successfully fitted (!) cluster results in a new branch
      
      
      for(unsigned int hIdx1=0;hIdx1<nHits;hIdx1++) {

	int h1 = vHL.at(hIdx1);

	TRACK_CLONE* newBr1 = &m_tClones[globalBranchCounter++];

	if(globalBranchCounter >= MAX_CLONE) {
	  errorFlag = true;
	  break;
	}

	newBr1->clone(**bIt);

	int code = m_fte.fastExtrapolator(&newBr1->m_ts, h1);

	if(code!=0) {
	  newBr1 = 0;
	  continue;
	}

	FIT_HIT fh(h1);

	update(&newBr1->m_ts, &fh, maxDChi2_upd);
	
	if(fh.m_status!=0 || std::isnan(fh.m_dchi2)) {
	  newBr1 = 0;
	  continue;
	} 
	  
	//std::cout<<"new hit added"<<std::endl;

	float dQ = dQ_hit - fh.m_dchi2;// - log(2*3.1415926*fh.m_detr);
	  
	newBr1->addHit(h1, dQ);
	newBranches.push_back(newBr1);bCounter++;
	nNewTotal++;
	nNoAdd = 0;
	
	int nNewBr = 0;

	float best_so_far = newBr1->m_Q;

	//std::cout<<"propagating a newly added branch .."<<std::endl;

	//try to add another hit
	for(unsigned int hIdx2=hIdx1+1;hIdx2<nHits;hIdx2++) {

	  int h2 = vHL.at(hIdx2);

	  if(m_hits->m_mod_id[h2] == m_hits->m_mod_id[h1]) continue;

	  if(newBr1 == 0) {
	    std::cout<<"ZERO BRANCH !!!!"<<std::endl;
	  }

	  TRACK_CLONE* newBr12 = &m_tClones[globalBranchCounter++];

	  if(globalBranchCounter >= MAX_CLONE) {
	    errorFlag = true;
	    break;
	  }

	  newBr12->clone(*newBr1);

	  int code = m_fte.fastExtrapolator(&newBr12->m_ts, h2, true);

	  if(code!=0) {
	    newBr12 = 0;
	    continue;
	  }

	  FIT_HIT fh2(h2);
	
	  update(&newBr12->m_ts, &fh2, maxDChi2_upd);
	
	  if(fh2.m_status!=0 || std::isnan(fh2.m_dchi2)) {
	    newBr12 = 0;
	    continue;
	  } 
	  
	  float dQ = dQ_hit - fh2.m_dchi2;// - log(2*3.1415926*fh2.m_detr);
	    
	  float newQ = dQ + newBr12->m_Q;

	  if(best_so_far > newQ) {
	    newBr12 = 0;
	    continue;
	  } 
	  else {
	    best_so_far = newQ;
	  }

	  newBr12->addHit(h2, dQ);
	  newBranches.push_back(newBr12);bCounter++;
	  nNewTotal++;
	  nNewBr++;

	  //try to add the third hit

	  for(unsigned int hIdx3=hIdx2+1;hIdx3<nHits;hIdx3++) {

	    int h3 = vHL.at(hIdx3);

	    if(m_hits->m_mod_id[h3] == m_hits->m_mod_id[h2]) continue;
	    if(m_hits->m_mod_id[h3] == m_hits->m_mod_id[h1]) continue;
	    if(newBr12 == 0) {
	      std::cout<<"ZERO BRANCH !!!!"<<std::endl;
	    }

	    TRACK_CLONE* newBr123 = &m_tClones[globalBranchCounter++];
	    if(globalBranchCounter>=MAX_CLONE) {
	      errorFlag = true;
	      break;
	    }
	    newBr123->clone(*newBr12);

	    int code = m_fte.fastExtrapolator(&newBr123->m_ts, h3, false);

	    if(code!=0) {
	      newBr123 = 0;
	      continue;
	    }

	    FIT_HIT fh3(h3);
	
	    update(&newBr123->m_ts, &fh3, maxDChi2_upd);
	
	    if(fh3.m_status!=0 || std::isnan(fh3.m_dchi2)) {
	      newBr123 = 0;
	      continue;
	    } 
	  
	    float dQ = dQ_hit - fh3.m_dchi2;// - log(2*3.1415926*fh2->m_detr);

	    newQ = dQ + newBr123->m_Q;

	    if(best_so_far > newQ) {
	      newBr123 = 0;
	      continue;
	    } 
	    else {
	      best_so_far = newQ;
	    }
	    

	    newBr123->addHit(h3, dQ);
	    newBranches.push_back(newBr123);
	    
	    bCounter++;
	    nNewTotal++;
	    nNewBr++;
	  }
	}
	if(errorFlag) break;
      }
      if(nNoAdd > nMaxNoAdd) {
	errorFlag = true;
      }
      if(errorFlag) break;
    }

    if(errorFlag) break;

    std::sort(newBranches.begin(), newBranches.end(), TRACK_CLONE::CompareQuality());

    //delete old branches

    vBranches.clear();

    int brIdx=0;

    for(std::vector<TRACK_CLONE*>::iterator bIt=newBranches.begin();bIt!=newBranches.end();++bIt, brIdx++) {
      if(brIdx <= nBest) vBranches.push_back(*bIt);
    }
    newBranches.clear();
    if(vBranches.empty()) break;
  }

  if(!errorFlag && !vBranches.empty()) {

    TRACK_CLONE* bestBr = (*vBranches.begin());
    if(bestBr->m_Q > minQ) {
      for(int k=0;k<bestBr->m_nHits;k++) {
	pT->m_hits[pT->m_nHits++] = bestBr->m_hits[k];
	if(pT->m_nHits >= MAX_HIT_ON_TRACK) break;
      }
      pT->m_beginTrackState.initialize(bestBr->m_ts);

    }
  }

}


void FastTrackFitter::forwardPropagation(NEW_TRACK* pT, TRAJECTORY& vPoints) {
  
  TRACK_STATE* pTS = &pT->m_endTrackState;
  
  int lastHit = pT->m_hits[pT->m_nHits-1];

  float theta = pTS->m_Rk[3];

  bool searchOutside = true;

  int vol = m_hits->m_vol_id[lastHit];
  int lay = m_hits->m_lay_id[lastHit];
  
  if(vol == 12 && lay == 2) return;
  if(vol == 14 && lay == 12) return;
  if(vol == 9 && lay == 14) {
    if(theta < 0.08) return;
  }
  if(vol == 7 && lay == 2) {
    if(theta > M_PI-0.08) return;
  }

  if(searchOutside) m_fte.extrapolateOutsideFast(pTS, vPoints);
  
}

int FastTrackFitter::addHits(NEW_TRACK* pT, std::vector<int>& vHits) {
 
  float dQ_hole = m_params.tf_f[49];
  float dQ_hit  = m_params.tf_f[50];
  float maxDChi2_upd = m_params.tf_f[51];
  float minQ = m_params.tf_f[52];

  int nBest = m_params.tf_i[1];

  int globalBranchCounter = 0;

  if(pT->m_nHits >= MAX_HIT_ON_TRACK) return 0;

  if(!pT->m_endTrackState.valid()) {
    std::cout<<"No valid end track state found, skipping the extension refit !"<<std::endl;
    for(std::vector<int>::iterator it=vHits.begin(); it!=vHits.end();it++) {
      pT->m_hits[pT->m_nHits++] = (*it);
      if(pT->m_nHits >= MAX_HIT_ON_TRACK) break;
    }
    return 0;
  }

  int nAdded = 0;

  //do branching exploration here ...

  //arrange hits in layers along the trajectory

  std::vector<int> hitArray[100];

  std::vector<unsigned int> keyArray;

  int nM = 0;

  for(std::vector<int>::iterator hIt=vHits.begin();hIt!=vHits.end();++hIt) {

    int h = (*hIt);

    if(m_hits->m_mask[h] == 2) continue;

    int vol_id = m_hits->m_vol_id[h];
    int lay_id = m_hits->m_lay_id[h];

    unsigned int layerKey = (vol_id-7) + (lay_id>>1);

    if(std::find(keyArray.begin(), keyArray.end(), layerKey) == keyArray.end()) {
      keyArray.push_back(layerKey);
    }
    hitArray[layerKey].push_back(h);
    nM++;
  }

  if(nM==0) return 0;

  std::vector<TRACK_CLONE*> vBranches;

  TRACK_CLONE* pClone = &m_tClones[globalBranchCounter++];
  pClone->initialize(&pT->m_endTrackState);
  vBranches.push_back(pClone);

  std::vector<TRACK_CLONE*> newBranches;

  int lIdx=0;

  bool errorFlag = false;
  
  for(;lIdx<(int)keyArray.size();lIdx++) {

    std::vector<int>& vHL = hitArray[keyArray[lIdx]];

    unsigned int nHits = vHL.size();

    int nNewTotal=0;

    for(std::vector<TRACK_CLONE*>::iterator bIt=vBranches.begin();bIt!=vBranches.end();++bIt) {

      int bCounter = 0;

      //1. add "skip-all", "no-detect-in-this layer" branch

      TRACK_CLONE* nodBr = &m_tClones[globalBranchCounter++];
      
      if(globalBranchCounter>=MAX_CLONE) {
	errorFlag = true;
	break;
      }

      nodBr->clone(**bIt);
      
      nNewTotal++;
      nodBr->addHole(dQ_hole);
      
      //std::cout<<"added hole, Q="<<nodBr->m_Q<<std::endl;

      newBranches.push_back(nodBr);

      bCounter++;
      
      float best_so_far = nodBr->m_Q;

      //2. investigate hits

      //arrange hits into 1- , 2- or 3-hit clusters -> each successfully fitted (!) cluster results in a new branch

      for(unsigned int hIdx1=0;hIdx1<nHits;hIdx1++) {

	int h1 = vHL.at(hIdx1);
	
	TRACK_CLONE* newBr1 = &m_tClones[globalBranchCounter++];
	if(globalBranchCounter>=MAX_CLONE) {
	  errorFlag = true;
	  break;
	}
	newBr1->clone(**bIt);

	int code = m_fte.fastExtrapolator(&newBr1->m_ts, h1, true);
	
	if(code!=0) {
	  newBr1 = 0;
	  continue;//just forget about this branch
	}
	
	FIT_HIT fh(h1);
	
	update(&newBr1->m_ts, &fh, maxDChi2_upd);
	
	if(fh.m_status!=0 || std::isnan(fh.m_dchi2)) {
	  newBr1 = 0;
	  continue;
	} 
	  
	//std::cout<<"new hit added"<<std::endl;

	float dQ = dQ_hit - fh.m_dchi2;// - log(2*3.1415926*fh->m_detr);
	
	newBr1->addHit(h1, dQ);

	best_so_far = newBr1->m_Q;
	
	newBranches.push_back(newBr1);
	bCounter++;
	nNewTotal++;

	int nNewBr = 0;

	//try to add another hit
	for(unsigned int hIdx2=hIdx1+1;hIdx2<nHits;hIdx2++) {

	  int h2 = vHL.at(hIdx2);

	  if(m_hits->m_mod_id[h2] == m_hits->m_mod_id[h1]) continue;

	  if(newBr1 == 0) {
	    std::cout<<"ZERO BRANCH !!!!"<<std::endl;
	  }

	  TRACK_CLONE* newBr12 = &m_tClones[globalBranchCounter++];
	  if(globalBranchCounter>=MAX_CLONE) {
	    errorFlag = true;
	    break;
	  }
	  newBr12->clone(*newBr1);

	  int code = m_fte.fastExtrapolator(&newBr12->m_ts, h2, true);

	  if(code!=0) {
	    newBr12 = 0;
	    continue;
	  }

	  FIT_HIT fh2(h2);
	
	  update(&newBr12->m_ts, &fh2, maxDChi2_upd);
	
	  if(fh2.m_status!=0 || std::isnan(fh2.m_dchi2)) {
	    newBr12 = 0;
	    continue;
	  } 
	  
	  float dQ = dQ_hit - fh2.m_dchi2;// - log(2*3.1415926*fh2->m_detr);
	    
	  float newQ = dQ + newBr12->m_Q;

	  if(best_so_far > newQ) {
	    newBr12 = 0;
	    continue;
	  } 
	  else {
	    best_so_far = newQ;
	  }


	  newBr12->addHit(h2, dQ);
	  newBranches.push_back(newBr12);

	  bCounter++;
	  nNewTotal++;
	  nNewBr++;

	  //try to add the third hit
	  for(unsigned int hIdx3=hIdx2+1;hIdx3<nHits;hIdx3++) {

	    int h3 = vHL.at(hIdx3);

	    if(m_hits->m_mod_id[h3] == m_hits->m_mod_id[h2]) continue;
	    if(m_hits->m_mod_id[h3] == m_hits->m_mod_id[h1]) continue;
	    if(newBr12 == 0) {
	      std::cout<<"ZERO BRANCH !!!!"<<std::endl;
	    }

	    TRACK_CLONE* newBr123 = &m_tClones[globalBranchCounter++];
	    if(globalBranchCounter>=MAX_CLONE) {
	      errorFlag = true;
	      break;
	    }
	    newBr123->clone(*newBr12);

	    int code = m_fte.fastExtrapolator(&newBr123->m_ts, h3, true);

	    if(code!=0) {
	      newBr123 = 0;
	      continue;
	    }

	    FIT_HIT fh3(h3);
	
	    update(&newBr123->m_ts, &fh3, maxDChi2_upd);
	
	    if(fh3.m_status!=0 || std::isnan(fh3.m_dchi2)) {
	      newBr123 = 0;
	      continue;
	    } 
	  
	    float dQ = dQ_hit - fh3.m_dchi2;// - log(2*3.1415926*fh2->m_detr);

	    newQ = dQ + newBr123->m_Q;

	    if(best_so_far > newQ) {
	      newBr123 = 0;
	      continue;
	    } 
	    else {
	      best_so_far = newQ;
	    }
	    

	    newBr123->addHit(h3, dQ);
	    newBranches.push_back(newBr123);
	    
	    bCounter++;
	    nNewTotal++;
	    nNewBr++;
	  }
	  if(errorFlag) break;
	}
	if(errorFlag) break;
      }

      if(errorFlag) break;
    } //loop over branches

    if(errorFlag) break;

    std::sort(newBranches.begin(), newBranches.end(), TRACK_CLONE::CompareQuality());

    //drop old branches

    vBranches.clear();

    int brIdx=0;

    for(std::vector<TRACK_CLONE*>::iterator bIt=newBranches.begin();bIt!=newBranches.end();++bIt, brIdx++) {
      if(brIdx <= nBest && (*bIt)->m_Q > 2*minQ) vBranches.push_back(*bIt);
    }
    newBranches.clear();
    if(vBranches.empty()) break;
  }

  if(!errorFlag && !vBranches.empty()) {
  
    TRACK_CLONE* bestBr = (*vBranches.begin());
    if(bestBr->m_Q > minQ) {

      for(int k=0;k<bestBr->m_nHits;k++) {
	pT->m_hits[pT->m_nHits++] = bestBr->m_hits[k];
	nAdded++;
	if(pT->m_nHits >= MAX_HIT_ON_TRACK) break;
      }
      pT->m_endTrackState.initialize(bestBr->m_ts);
    }
  }
  return nAdded;
}

void FastTrackFitter::update(TRACK_STATE* pTS, FIT_HIT* fh, float maxDChi2) const {
  
  int h = fh->m_h;
  int vol_id = m_hits->m_vol_id[h];
  if(vol_id == 8) {//pixel

    float cost = fabs(pTS->m_Re[3]-0.5*M_PI);
    int nCells = m_hits->m_clusterLength[h];

    bool drop = false;
    if(nCells < 7) drop = cost > m_params.tf_table1[nCells];
    else drop = cost > m_params.tf_table1[8];
    
    if(drop) {//reject this hit

      fh->m_status = -1;
      
      memcpy(&pTS->m_Rk[0], &pTS->m_Re[0], sizeof(pTS->m_Re));
      memcpy(&pTS->m_Ck[0], &pTS->m_Ce[0], sizeof(pTS->m_Ce));
      
      return;
    }
  }

  float m[3] = {m_hits->m_m1[h], m_hits->m_m2[h], 0.0};

  double R[2][2];

  R[0][1] = R[1][0] = 0.0;
  R[0][0] = m_params.tf_d[4];
  R[1][1] = m_params.tf_d[5];
  
  if(vol_id < 10) {//pixel
    R[0][0] = m_params.tf_d[0];
    R[1][1] = m_params.tf_d[1];
  } 
  else if(vol_id < 16) {
    R[0][0] = m_params.tf_d[2];
    R[1][1] = m_params.tf_d[3];
  }
  
  float Ce[16];

  memcpy(&Ce[0], &pTS->m_Ce[0], sizeof(pTS->m_Ce));

  
  double CHT[5][2];
  double D[2][2];
  float resid[2];

  if(pTS->m_pDE->m_isRect) {
    
    CHT[0][0] = Ce[0];
    CHT[0][1] = Ce[1];

    CHT[1][0] = Ce[1];
    CHT[1][1] = Ce[2];

    CHT[2][0] = Ce[3];
    CHT[2][1] = Ce[4];

    CHT[3][0] = Ce[6];
    CHT[3][1] = Ce[7];

    CHT[4][0] = Ce[10];
    CHT[4][1] = Ce[11];

    D[0][0] = Ce[0] + R[0][0];
    D[0][1] = D[1][0] = Ce[1] + R[0][1];
    D[1][1] = Ce[2] + R[1][1];

    resid[0] = m[0] - pTS->m_Re[0];
    resid[1] = m[1] - pTS->m_Re[1];
  }
  else {

    float sc[2];
    pTS->m_pDE->getResidualTransform(m[0], m[1], sc);

    float p1 = sc[0]*Ce[1];
    float p2 = sc[1]*Ce[1];

    CHT[0][0] = sc[1]*Ce[0] - p1;
    CHT[0][1] = sc[0]*Ce[0] + p2;

    CHT[1][0] = p2 - sc[0]*Ce[2];
    CHT[1][1] = p1 + sc[1]*Ce[2];

    CHT[2][0] = sc[1]*Ce[3] - sc[0]*Ce[4];
    CHT[2][1] = sc[0]*Ce[3] + sc[1]*Ce[4];
    
    CHT[3][0] = sc[1]*Ce[6] - sc[0]*Ce[7];
    CHT[3][1] = sc[0]*Ce[6] + sc[1]*Ce[7];
    
    CHT[4][0] = sc[1]*Ce[10] - sc[0]*Ce[11];
    CHT[4][1] = sc[0]*Ce[10] + sc[1]*Ce[11];


    float sc12 = sc[1]*sc[1];
    float sc02 = sc[0]*sc[0];
    float scx  = sc[0]*sc[1];

    float scG = 2*scx*Ce[1];

    double RG00 = R[0][0]+Ce[0];
    double RG11 = R[1][1]+Ce[2];
    
    D[0][0] = sc12*RG00 + sc02*RG11 - scG;
    D[0][1] = scx*(RG00-RG11);
    D[1][1] = sc02*RG00 + sc12*RG11 + scG;
    
    D[0][1] += (sc12-sc02)*Ce[1];
    D[1][0] = D[0][1];
    
    float d[2] = {m[0] - pTS->m_Re[0], m[1] - pTS->m_Re[1]};
    resid[0] = sc[1]*d[0] - sc[0]*d[1];
    resid[1] = sc[0]*d[0] + sc[1]*d[1];
    
  }

  double detr = 1.0/(D[0][0]*D[1][1] - D[0][1]*D[1][0]);

  //std::cout<<"detr="<<detr<<std::endl;

  double D_1[2][2];
  
  D_1[0][0] = detr*D[1][1];
  D_1[0][1] =-detr*D[1][0];
  D_1[1][0] = D_1[0][1];
  D_1[1][1] = detr*D[0][0];

  float dchi2 = resid[0]*resid[0]*D_1[0][0] + resid[1]*(2*resid[0]*D_1[0][1] + resid[1]*D_1[1][1]);

  fh->m_dchi2 = dchi2;  
  fh->m_detr  = detr;
  
  if(dchi2>maxDChi2) {//outlier

    fh->m_status = -1;

    memcpy(&pTS->m_Rk[0], &pTS->m_Re[0], sizeof(pTS->m_Re));
    memcpy(&pTS->m_Ck[0], &Ce[0], sizeof(pTS->m_Ce));
    return;
  }
  
  double K[5][2];

  for(int i=0;i<5;i++)
    for(int j=0;j<2;j++) {
      K[i][j] = CHT[i][0]*D_1[0][j] + CHT[i][1]*D_1[1][j];
    }

  pTS->m_chi2 += dchi2;
  pTS->m_ndof += 2;
  fh->m_status = 0;

  for(int i=0;i<5;i++) pTS->m_Rk[i] = pTS->m_Re[i] + K[i][0]*resid[0] + K[i][1]*resid[1];

  int idx = 0;
  for(int i=0;i<5;i++) {
    for(int j=0;j<=i;j++, idx++) {
      pTS->m_Ck[idx] = Ce[idx] - K[i][0]*CHT[j][0] - K[i][1]*CHT[j][1];
    }
  }

  float& a = pTS->m_Rk[2];

  if(a>M_PI) 
    a=a-2*M_PI; 
  else {
    if(a<-M_PI) 
      a=a+2*M_PI;
  }

  float& b = pTS->m_Rk[3];

  if(b<0)  b+=M_PI;
  if(b>M_PI) b-=M_PI;

}

void FastTrackFitter::compareTrackStates(TRACK_STATE* p1, TRACK_STATE* p2) const {

  std::cout<<"estimated track params:"<<std::endl;
  for(int i=0;i<5;i++) {
    std::cout<<p1->m_Rk[i]<<"               "<<p2->m_Rk[i]<<std::endl;
  }
  std::cout<<"extrapolated track params:"<<std::endl;
  for(int i=0;i<5;i++) {
    std::cout<<p1->m_Re[i]<<"               "<<p2->m_Re[i]<<std::endl;
  }
  std::cout<<"estimated covariance :"<<std::endl;
  int idx1=0, idx2=0;

  for(int i=0;i<5;i++) {
    for(int j=0;j<=i;j++, idx1++) {
      std::cout<<p1->m_Ck[idx1]<<" ";
    }
    std::cout<<"     vs    ";
    for(int j=0;j<=i;j++, idx2++) {
      std::cout<<p2->m_Ck[idx2]<<" ";
    }
    std::cout<<std::endl;
  }

  std::cout<<"extrapolated covariance :"<<std::endl;
  idx1=0; 
  idx2=0;

  for(int i=0;i<5;i++) {
    for(int j=0;j<=i;j++, idx1++) {
      std::cout<<p1->m_Ce[idx1]<<" ";
    }
    std::cout<<"     vs    ";
    for(int j=0;j<=i;j++, idx2++) {
      std::cout<<p2->m_Ce[idx2]<<" ";
    }
    std::cout<<std::endl;
  }
}

