#include<iostream>
#include<fstream>
#include<iomanip>
#include<cstring>
#include<cmath>
#include<sstream>

#include <omp.h>

#include<algorithm>

#include "TrackFinder.h"

#include "ClusterMaker.h"
#include "TrackingFilter.h"
#include "FastTrackExtrapolator.h"
#include "FastTrackFitter.h"

#include<set>
#include<bitset>

#include <chrono>

TrackFinder::TrackFinder(const ALGORITHM_PARAMETERS& p, const GEOMETRY& g, const LAYER_LINKER& ll) : m_params(p), m_geo(g), m_linker(ll), m_H2T(0), m_trackBank(0), m_eventCounter(0) {

  const float deltaEta = m_params.sf_f[0];


  m_hitBank = new HIT_BANK(m_geo, deltaEta);

  for(int i=0;i<N_SEG_BANKS;i++) {
    m_segBank[i] = new NEW_SEG_BANK;
  }

  m_globalTrackCounter = 0;

  m_trackBank = new NEW_TRACK[N_TRACKS];

  for(int i=0;i<5;i++) m_binTable[i] = 0;

  initializeCuts();
  m_binTable[1] = createBinTable(1);
  m_binTable[2] = createBinTable(2);
  m_binTable[3] = createBinTable(3);
  m_binTable[4] = createBinTable(4);

  memset(&m_time[0], 0, sizeof(m_time));
  memset(&m_segStat[0], 0, sizeof(m_segStat));
}

TrackFinder::~TrackFinder() {
  delete m_hitBank;
  for(int i=0;i<N_SEG_BANKS;i++) {
    delete m_segBank[i];
  }

  delete[] m_trackBank;

  for(int i=0;i<5;i++) delete[] m_binTable[i];
  delete[] m_H2T;

  /*
  double sum = 0.0;

  std::string timerNames[15] = {"fillIndices", "prepareL1Nodes", "createSegments", "connectSegments",
				"evolveNetwork", "collectTracks", "prefitTracks", "removeClones",
				"fitTracks", "extendTracks", "reassignHits",
				"maskHits", "removeUsed", "prepare", "clean"};

  for(int i=0;i<15;i++) {
    std::cout<<"Timer"<<i<<" "<<timerNames[i]<<" "<<m_time[i]/m_eventCounter<<" ms "<<std::endl;
    sum += m_time[i];
  }

  std::cout<<"Processing time per event  "<<sum/m_eventCounter<<" ms "<<std::endl;
  */
}

void TrackFinder::prepare(const HITS* pE) {

  auto start = std::chrono::high_resolution_clock::now();

  m_pH = pE;
  m_hitBank->importHits(pE);
  m_hitBank->sort();
  findClusters();
  m_globalTrackCounter = 0;
  m_H2T = new int[m_pH->m_nHits+1];

  auto finish = std::chrono::high_resolution_clock::now();

  m_time[13] += std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count();
}

void TrackFinder::clean() {

  auto start = std::chrono::high_resolution_clock::now();

  m_hitBank->clear();
  delete[] m_H2T;
  m_H2T = 0;

  auto finish = std::chrono::high_resolution_clock::now();

  m_time[14] += std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count();

  m_eventCounter++;
}


void TrackFinder::run(int stage, const HITS* pE, std::vector<NEW_TRACK*>& vTracks, int thr) {

  m_stage = stage;
  m_pH = pE;

  auto start = std::chrono::high_resolution_clock::now();

  m_hitBank->fillIndices();

  auto finish = std::chrono::high_resolution_clock::now();
  m_time[0] += std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count();
  start = finish;

  prepareL1Nodes();//Layer 1 nodes for track segment creation

  finish = std::chrono::high_resolution_clock::now();
  m_time[1] += std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count();
  start = finish;

  createSegments();

  finish = std::chrono::high_resolution_clock::now();

  auto deltaT = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count();

  m_time[2] += deltaT;

  m_time[19+m_stage] += deltaT;

  start = finish;

  connectSegments();

  finish = std::chrono::high_resolution_clock::now();
  m_time[3] += std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count();
  start = finish;

  runNetworkEvolution();

  finish = std::chrono::high_resolution_clock::now();
  m_time[4] += std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count();
  start = finish;

  vTracks.clear();

  collectTracks(vTracks);

  finish = std::chrono::high_resolution_clock::now();
  m_time[5] += std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count();
  start = finish;

  prefitTracks(vTracks);

  std::sort(vTracks.begin(), vTracks.end(), NEW_TRACK::CompareLikelihood());

  finish = std::chrono::high_resolution_clock::now();
  m_time[6] += std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count();
  start = finish;

  removeClones(vTracks);

  finish = std::chrono::high_resolution_clock::now();
  m_time[7] += std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count();
  start = finish;

  fitTracks(vTracks, m_stage);

  std::sort(vTracks.begin(), vTracks.end(), NEW_TRACK::CompareLikelihood());

  finish = std::chrono::high_resolution_clock::now();
  m_time[8] += std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count();
  start = finish;

  runTrackExtensions(vTracks);

  finish = std::chrono::high_resolution_clock::now();
  m_time[9] += std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count();
  start = finish;

  reAssignHits(vTracks);


  finish = std::chrono::high_resolution_clock::now();
  m_time[10] += std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count();
  start = finish;

  removeClones(vTracks);

  finish = std::chrono::high_resolution_clock::now();
  m_time[7] += std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count();
  start = finish;

  if(thr!=-2) maskHits(vTracks, thr);

  for(int i=0;i<N_SEG_BANKS;i++) {
    m_segBank[i]->clear();
  }

  finish = std::chrono::high_resolution_clock::now();
  m_time[11] += std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count();
  start = finish;


  if(m_stage < 4) m_hitBank->deleteUsedNodes(m_pH);

  finish = std::chrono::high_resolution_clock::now();

  m_time[12] += std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count();

}

void TrackFinder::fitTracks(std::vector<NEW_TRACK*>& vTracks, int stage) {

  int nBlocks = N_THREADS;

  int nTracks = vTracks.size();

  omp_set_num_threads(nBlocks);

  int tid;

#pragma omp parallel private(tid) shared(vTracks, nTracks, nBlocks)
  {

    tid = omp_get_thread_num();

    FTE trackExtrapolator(m_params, m_geo, m_pH);
    FTF trackFit(m_params, m_geo, m_pH, trackExtrapolator);

    for(int trackIdx=tid;trackIdx<nTracks;trackIdx += nBlocks) {

      NEW_TRACK* pT = vTracks.at(trackIdx);

      if(stage == 1) {
      	if(pT->m_nHits<=4) continue;
      }

      trackFit.fitTrack(pT, stage);

    }
  }
}

void TrackFinder::runTrackExtensions(std::vector<NEW_TRACK*>& vTracks) {

  int maxLength_inside = 200;
  int maxLength_outside = 200;

  int nBlocks = N_THREADS;

  int nTracks = vTracks.size();

  omp_set_num_threads(nBlocks);

  int tid;

  //tracks are split into groups, groups processed sequentially while in-group tracks are processed in parallel

  int groupSize = 10*nBlocks;
  int firstTrack = 0;

  while(firstTrack < nTracks) {
    int lastTrack  = firstTrack + groupSize;
    if(lastTrack > nTracks) lastTrack = nTracks;

#pragma omp parallel private(tid) shared(vTracks, firstTrack, lastTrack, nBlocks)
    {

      tid = omp_get_thread_num();

      FTE trackExtrapolator(m_params, m_geo, m_pH);
      FTF trackFitter(m_params, m_geo, m_pH, trackExtrapolator);

      for(int trackIdx=firstTrack+tid;trackIdx<lastTrack;trackIdx += nBlocks) {

	NEW_TRACK* pT = vTracks.at(trackIdx);

	TRAJECTORY vPointsIn;

	vPointsIn.clear();

	bool proceed = true;

	//1. extend inside

	int h = pT->m_hits[0];

	if(m_pH->m_vol_id[h] == 8 && m_pH->m_lay_id[h] == 2) proceed = false;

	if(proceed) {
	  trackFitter.backwardPropagation(pT, vPointsIn);
	  if(vPointsIn.m_nPoints > maxLength_inside) proceed = false;
	  if(proceed) {
	    std::vector<int> vHits;

	    m_hitBank->collectHitsInward(vPointsIn, vHits);

	    if(!vHits.empty()) {
	      trackFitter.addHitsBackward(pT, vHits);
	      pT->update(m_pH);
	      pT->removeDuplicates();
	    }
	  }
	}

	//2. extend outside

	TRAJECTORY vPointsOut;
	vPointsOut.clear();

	trackFitter.forwardPropagation(pT, vPointsOut);

	proceed = true;

	if(vPointsOut.m_nPoints < 2 || vPointsOut.m_nPoints > maxLength_outside) proceed = false;

	if(proceed) {
	  std::vector<int> vHits;
	  std::vector<int> vLayers;
	  pT->getLayers(vLayers, m_pH);
	  m_hitBank->collectHits(vLayers, vPointsOut, vHits, m_pH);//mask them later
	  if(!vHits.empty()) {
	    int nAdded = trackFitter.addHits(pT, vHits);
	    if(nAdded != 0) {
	      pT->update(m_pH);
	      pT->removeDuplicates();
	    }
	  }
	}
      }
    }//end of parallel section
    //mask all hits on the tracks in the group
    for(int trackId=firstTrack;trackId<lastTrack;trackId++) {
      NEW_TRACK* pT = vTracks.at(trackId);
      pT->maskHits(2, m_pH);
    }
    firstTrack += groupSize;
  }

  //restore the masks

  memset(m_pH->m_mask, 0, m_pH->m_nHits*sizeof(int));
}

void TrackFinder::reAssignHits(std::vector<NEW_TRACK*>& vTracks) {

  int nBlocks = N_THREADS;

  int nTracks = vTracks.size();

  omp_set_num_threads(nBlocks);

  int tid;

#pragma omp parallel private(tid) shared(vTracks, nTracks, nBlocks)
  {

    tid = omp_get_thread_num();

    FTE trackExtrapolator(m_params, m_geo, m_pH);
    FTF trackFitter(m_params, m_geo, m_pH, trackExtrapolator);

    int nhits_cut = 2;

    if(m_stage == 1) nhits_cut = 5;
    if(m_stage == 2) nhits_cut = 4;

    TRAJECTORY vPoints;

    for(int trackIdx=tid;trackIdx<nTracks;trackIdx+=nBlocks) {

      NEW_TRACK* pT = vTracks.at(trackIdx);

      if(pT->m_nHits > nhits_cut) {

	vPoints.m_nPoints = 0;

	bool result = trackFitter.getTrajectory(pT, vPoints);

	if(result) {

	  std::vector<int> vHits;

	  std::vector<int> vLayers;

	  pT->getLayers(vLayers, m_pH);

	  m_hitBank->findClosestHits(vLayers, vPoints, vHits, m_pH);

	  if(!vHits.empty()) {
	    pT->add(vHits, m_pH);
	  }
	}
      }

    }

  }

}

void TrackFinder::removeClones(std::vector<NEW_TRACK*>& vTracks) {

  float maxPhiDist_barr = m_params.sf_f[1];
  float maxRzDist_barr  = m_params.sf_f[2];
  float maxPhiDist_all  = m_params.sf_f[3]*maxPhiDist_barr;
  float maxRzDist_all   = m_params.sf_f[4]*maxRzDist_barr;

  assignCloneFlags(vTracks);

  std::vector<NEW_TRACK*> vTmp;

  vTmp.clear();
  vTmp.reserve(vTracks.size());

  for(std::vector<NEW_TRACK*>::iterator it = vTracks.begin();it!=vTracks.end();++it) {
    if((*it)->m_cloneFlag == 0) vTmp.push_back(*it);
  }

  for(unsigned int t = 0; t<vTracks.size();t++) {

    NEW_TRACK* pClone = vTracks[t];

    if(pClone->m_cloneFlag <= 0) continue;

    int tid = pClone->m_cloneFlag;

    NEW_TRACK* pDst = vTracks[tid-1];

    if(pDst->m_nHits >= MAX_HIT_ON_TRACK) continue;

    if(pDst->m_cloneFlag != 0) continue;

    //merging (*it) -> pDst

    std::vector<int> hitVec;

    std::vector<unsigned int> modVec;

    modVec.resize(pDst->m_nHits);


    for(int i=0;i<pDst->m_nHits;i++) {
      int h = pDst->m_hits[i];
      modVec[i] = 1000000*m_pH->m_vol_id[h] + 10000*m_pH->m_lay_id[h] + m_pH->m_mod_id[h];
    }

    //source track

    for(int i=0;i<pClone->m_nHits;i++) {

      int h = pClone->m_hits[i];

      unsigned int key = 1000000*m_pH->m_vol_id[h] + 10000*m_pH->m_lay_id[h] + m_pH->m_mod_id[h];

      if(std::find(modVec.begin(), modVec.end(), key) != modVec.end()) continue;

      hitVec.push_back(h);
    }

    if(hitVec.empty()) continue;

    //validate distance

    std::vector<int> vNew;

    for(std::vector<int>::iterator kIt=hitVec.begin(); kIt!=hitVec.end();++kIt) {

      int h = (*kIt);

      int vol_id = m_pH->m_vol_id[h];
      int lay_id = m_pH->m_lay_id[h];

      float isBarrel = vol_id == 8 || vol_id == 13;

      float z = m_pH->m_z[h];
      float phi = m_pH->m_phi[h];
      float r = m_pH->m_r[h];

      int nHL = 0;

      float dPhiMin = 10000.0;
      float dRzMin = 10000.0;

      for(int dstIdx=0;dstIdx<pDst->m_nHits;dstIdx++) {

	int dstH = pDst->m_hits[dstIdx];

	int vol_id_2 = m_pH->m_vol_id[dstH];

	if(vol_id != vol_id_2) continue;

	int lay_id_2 = m_pH->m_lay_id[dstH];

	if(lay_id != lay_id_2) continue;

	nHL++;

	float dPhi = m_pH->m_phi[dstH] - phi;

	if(dPhi<-M_PI) dPhi += 2*M_PI;
	if(dPhi> M_PI) dPhi -= 2*M_PI;

	if(fabs(dPhi)<dPhiMin) dPhiMin = fabs(dPhi);

	if(isBarrel) {
	  float dz = m_pH->m_z[dstH] - z;
	  if(fabs(dz) < dRzMin) dRzMin = fabs(dz);
	}
	else {
	  float dr = m_pH->m_r[dstH] - r;
	  if(fabs(dr) < dRzMin) dRzMin = fabs(dr);
	}
      }

      if(nHL==0) vNew.push_back(h);

      else {
	if(isBarrel) {
	  if(dPhiMin < maxPhiDist_barr && dRzMin < maxRzDist_barr)  vNew.push_back(h);
	}
	else {
	  if(dPhiMin < maxPhiDist_all && dRzMin < maxRzDist_all)  vNew.push_back(h);
	}
      }

    }


    for(unsigned int k=0;k<vNew.size();k++) {
      pDst->m_hits[pDst->m_nHits++] = vNew[k];
      if(pDst->m_nHits >= MAX_HIT_ON_TRACK) break;
    }

    pClone->m_cloneFlag = -1;//merged, to be deleted

  }

  //store tracks

  vTracks = std::move(vTmp);
}


void TrackFinder::prefitTracks(std::vector<NEW_TRACK*>& vTracks) {

  float chi2cut = 5.5;//was 5.6

  if(m_stage == 1) chi2cut = 5.4;
  if(m_stage == 2) chi2cut = 4.4;
  if(m_stage == 3) chi2cut = 4.9;
  if(m_stage == 4) chi2cut = 5.7;

  int nBlocks = N_THREADS;

  std::vector<NEW_TRACK*> vBigVec[N_THREADS];

  for(int i=0;i<N_THREADS;i++) vBigVec[i].reserve(vTracks.size());

  int vecId=0;

  for(std::vector<NEW_TRACK*>::iterator it=vTracks.begin();it!=vTracks.end();++it) {
    vBigVec[vecId].push_back(*it);
    vecId++;
    if(vecId == N_THREADS) vecId = 0;
  }

  omp_set_num_threads(nBlocks);

  int tid;

#pragma omp parallel private(tid) shared(vBigVec)
  {

    tid = omp_get_thread_num();
    std::vector<NEW_TRACK*>& vTr = vBigVec[tid];
    FTE trackExtrapolator(m_params, m_geo, m_pH);
    FTF trackFit(m_params, m_geo, m_pH, trackExtrapolator,tid);

    for(unsigned int trackIdx = 0;trackIdx<vTr.size();trackIdx++) {

      NEW_TRACK* pT = vTr.at(trackIdx);

      pT->m_status = 1;
      bool result = trackFit.fitTrackCandidate(pT, m_stage);

      if(!result || pT->m_nHits<3 ) {
	pT->m_status = 0;
	continue;
      }
      if(pT->m_chi2<0 || pT->m_ndof<=0) {
	pT->m_status = 0;
	continue;
      }

      float chi2n = pT->m_chi2/pT->m_ndof;

      if(chi2n>chi2cut) {
	pT->m_status = 0;
      }
    }
  } // the end of parallel section

  std::vector<NEW_TRACK*> vTmp = std::move(vTracks);

  vTracks.clear();

  for(std::vector<NEW_TRACK*>::iterator it=vTmp.begin();it!=vTmp.end();++it) {
    if((*it)->m_status == 1) vTracks.push_back(*it);
  }

}



void TrackFinder::collectTracks(std::vector<NEW_TRACK*>& vTracks) {

  int minLevel = 4;

  if(m_stage == 1) minLevel = 4;//min requirement 4 segments on a track
  if(m_stage == 2) minLevel = 3;
  if(m_stage == 3) minLevel = 3;
  if(m_stage == 4) minLevel = 2;

  std::vector<const SEGMENT*> vBigSeg[N_THREADS];

  int vecId = 0;

  int nBlocks = N_THREADS;

  for(int bankId=0;bankId<N_SEG_BANKS;bankId++) {
    for(unsigned int segmentIndex=0;segmentIndex<m_segBank[bankId]->m_nSegments;segmentIndex++) {
      const SEGMENT* pS = &m_segBank[bankId]->m_S[segmentIndex];
      int level = pS->m_level;
      if(level<minLevel) continue;
      vBigSeg[vecId].push_back(pS);
      vecId++;
      if(vecId == nBlocks) vecId = 0;
    }
  }

  omp_set_num_threads(nBlocks);

  int tid;

  std::vector<int> vA[N_THREADS];

#pragma omp parallel private(tid) shared(vBigSeg, vA)
  {

    tid = omp_get_thread_num();

    std::vector<int>& newTracks = vA[tid];
    newTracks.reserve(N_TRACKS);

    std::vector<const SEGMENT*>& vSeg = vBigSeg[tid];

    float Qcut = -100.0;

    if(m_stage == 1) Qcut = -1.7;
    if(m_stage == 2) Qcut = -2.0;
    if(m_stage == 3) Qcut = -2.0;
    if(m_stage == 4) Qcut = -2.0;

    TRACKING_FILTER tFilter(m_params, m_pH, m_segBank);

    for(std::vector<const SEGMENT*>::iterator sIt = vSeg.begin();sIt!=vSeg.end();++sIt) {

      const SEGMENT* pS = (*sIt);

      std::vector<const SEGMENT*> chain;
      std::vector<const NODE*> track;
      float Q;

      Q = tFilter.followTrack(pS, chain, track);

      if(chain.size() < 2) continue;
      if(Q < Qcut) continue;

      int vH[100];
      int nHits = 0;

      for(std::vector<const NODE*>::iterator nIt=track.begin();nIt!=track.end();++nIt) {
	for(unsigned int s=0;s<(*nIt)->m_hits.size();s++) {
	  vH[nHits++] = (*nIt)->m_hits[s];
	  if(nHits>=100) break;
	}
      }

      int newTrackIdx;

#pragma omp critical
      {
	newTrackIdx = m_globalTrackCounter;
	m_globalTrackCounter++;
      }
      if(newTrackIdx < N_TRACKS) {

	newTracks.push_back(newTrackIdx);

	NEW_TRACK* pT = &m_trackBank[newTrackIdx];

	pT->initialize(newTrackIdx, vH, nHits, Q);
	pT->update(m_pH);
      }
    }
  }

  vTracks.clear();

  if(m_globalTrackCounter > N_TRACKS) {
    std::cout<<"too many tracks "<<m_globalTrackCounter<<std::endl;
    m_globalTrackCounter = N_TRACKS;
  }

  for(int b=0;b<nBlocks;b++) {
    for(std::vector<int>::iterator it=vA[b].begin();it!=vA[b].end();it++) {
      vTracks.push_back(&m_trackBank[*it]);
    }
  }

}

void TrackFinder::runNetworkEvolution() {

  const int maxIter = 30;

  int iter = 0;

  std::vector<SEGMENT*> uSegs_old[N_THREADS];
  std::vector<SEGMENT*> uSegs_new[N_THREADS];

  int vecId = 0;

  for(int bankId=0;bankId<N_SEG_BANKS; bankId++) {
    for(unsigned int segmentIndex=0;segmentIndex<m_segBank[bankId]->m_nSegments;segmentIndex++) {
      SEGMENT* pS = &m_segBank[bankId]->m_S[segmentIndex];
      if(pS->m_nNei == 0) continue;
      uSegs_old[vecId].push_back(pS);
      vecId++;
      if(vecId == N_THREADS) vecId = 0;
    }
  }

  omp_set_num_threads(N_THREADS);

  int nBlocks = N_THREADS;

  int tid;

  for(;iter<maxIter;iter++) {

    //generate proposals

#pragma omp parallel private(tid) shared(nBlocks, uSegs_old, uSegs_new)
    {

      tid = omp_get_thread_num();

      uSegs_new[tid].clear();

      for(std::vector<SEGMENT*>::iterator sIt = uSegs_old[tid].begin();sIt != uSegs_old[tid].end();++sIt) {
	SEGMENT* pS = (*sIt);
	int next_level = pS->m_level;

	for(int nIdx=0;nIdx<pS->m_nNei;nIdx++) {

	  unsigned int segIdx = pS->m_vNei[nIdx];
	  unsigned int nextBankId = (segIdx & SEG_BANK_MASK) >> SEG_MASK_SHIFT;
	  unsigned int nextSegmentIdx = (segIdx & SEG_INDEX_MASK);

	  SEGMENT* pN = &(m_segBank[nextBankId]->m_S[nextSegmentIdx]);

	  if(pS->m_level == pN->m_level) {
	    next_level = pS->m_level + 1;
	    uSegs_new[tid].push_back(pS);
	    break;
	  }
	}
	pS->m_next = next_level;
      }
    }

    //update

    int nChanges[N_THREADS];

#pragma omp parallel private(tid) shared(nBlocks, nChanges)
    {

      tid = omp_get_thread_num();

      nChanges[tid] = 0;

      for(std::vector<SEGMENT*>::iterator it=uSegs_new[tid].begin();it!=uSegs_new[tid].end();++it) {

	SEGMENT* pS = (*it);
	if(pS->m_next != pS->m_level) {
	  nChanges[tid]++;
	  pS->m_level = pS->m_next;
	}
      }
    }

    int nChTotal = 0;

    for(int idx=0;idx<nBlocks;idx++) {
      nChTotal += nChanges[idx];
    }
    //std::cout<<"iter = "<<iter<<" nChanges = "<<nChTotal<<std::endl;
    if(nChTotal == 0) break;

    for(int i=0;i<nBlocks;i++) {
      uSegs_old[i] = std::move(uSegs_new[i]);
      uSegs_new[i].clear();
    }
  }
}



void TrackFinder::findClusters() {

  std::vector<HIT_LAYER*> vL, vTmp;

  m_hitBank->getLayers(vTmp);


  for(std::vector<HIT_LAYER*>::iterator lIt = vTmp.begin();lIt!=vTmp.end();++lIt) {

    HIT_LAYER* pL = (*lIt);

    if(pL->m_info.m_vol_id >= 16) continue;
    vL.push_back(pL);
  }

  int nBlocks  = N_THREADS;

  int nLayers = vL.size();

  omp_set_num_threads(nBlocks);

  int tid;

#pragma omp parallel private(tid) shared(vL, nLayers, nBlocks)
  {

    tid = omp_get_thread_num();

    CLUSTER_MAKER cm(m_pH, m_geo);

    for(int layIdx=tid;layIdx<nLayers;layIdx+=nBlocks) {

      HIT_LAYER* pL = vL.at(layIdx);

      cm.findClusters(m_params, pL);
    }
  }

}

void TrackFinder::prepareL1Nodes() {

  std::vector<HIT_LAYER*> vL, vTmp;

  m_hitBank->getLayers(vTmp);


  for(std::vector<HIT_LAYER*>::iterator lIt = vTmp.begin();lIt!=vTmp.end();++lIt) {

    HIT_LAYER* pL = (*lIt);

    if(pL->m_info.m_vol_id >= 16) continue;



    vL.push_back(pL);
  }

  int nBlocks  = N_THREADS;

  int nLayers = vL.size();

  omp_set_num_threads(nBlocks);

  int tid;

#pragma omp parallel private(tid) shared(vL, nLayers, nBlocks)
  {

    tid = omp_get_thread_num();

    for(int layIdx=tid;layIdx<nLayers;layIdx+=nBlocks) {

      HIT_LAYER* pL = vL.at(layIdx);
      pL->prepareL1Nodes();
    }
  }

}


void TrackFinder::createSegments() {


  //layer pairs allocation for load balancing with 4 segment banks

  std::vector<int> lT0 = {33, 37, 48, 54, 19, 7,  46, 27, 44, 13, 15, 39, 53, 57, 18};
  std::vector<int> lT1 = {0,  2,  32, 30, 58, 49, 22, 23, 47, 12, 28, 40, 52, 43, 17};
  std::vector<int> lT2 = {3,  1,  31, 36, 35, 59, 21, 24, 9,  11, 26, 50, 51, 38, 56};
  std::vector<int> lT3 = {4,  5,  6,  29, 34, 55, 20, 8,  10, 25, 45, 14, 16, 41, 42};
  std::vector<int>* loadTable[4] = {&lT0, &lT1, &lT2, &lT3};


  //selecting layer pairs

  std::vector<const LAYER_LINK*> vGoodLinks;

  const std::vector<LAYER_LINK>& LL = m_linker.m_links;

  int pairIdx = 0;

  for(std::vector<LAYER_LINK>::const_iterator it=LL.begin();it!=LL.end();++it, pairIdx++) {

    if((*it).m_src >= 16000 || (*it).m_dst >= 16000) continue;//skipping the long strips
    vGoodLinks.push_back(&(*it));
  }

  int nBlocks = N_SEG_BANKS;

  omp_set_num_threads(nBlocks);

  int nLinkPairs = vGoodLinks.size();

  nLinkPairs -= 2;//drop the last 2 the least important links

  m_maxLink = nLinkPairs;

  int tid;

#pragma omp parallel private(tid) shared(vGoodLinks, nLinkPairs)
  {

    tid = omp_get_thread_num();

    for(std::vector<int>::iterator it=loadTable[tid]->begin();it!=loadTable[tid]->end();++it) {
      int pairIdx = (*it);
      if(pairIdx == -1) continue;
      const LAYER_LINK* lll = vGoodLinks.at(pairIdx);
      createSegments(*lll, tid);
    }

    /*
    for(int pairIdx=tid;pairIdx<nLinkPairs;pairIdx+=nBlocks) {
      const LAYER_LINK* lll = vGoodLinks.at(pairIdx);
      createSegments(*lll, tid);
    }
    */
  }

  createLinks();

}

void TrackFinder::createLinks() {

  for(int bankId=0;bankId<N_SEG_BANKS;bankId++) {

    unsigned int bankMask = bankId << SEG_MASK_SHIFT;

    for(unsigned int segmentIndex=0;segmentIndex<m_segBank[bankId]->m_nSegments;segmentIndex++) {

      unsigned int segId = bankMask | segmentIndex;

      m_segBank[bankId]->m_S[segmentIndex].m_n1->addIn(segId);
      m_segBank[bankId]->m_S[segmentIndex].m_n2->addOut(segId);
    }
  }

}



typedef struct SegIndex {

  struct CompareTau {
    bool operator()(const struct SegIndex& s1, const struct SegIndex& s2) {
      return s1.tau < s2.tau;
    }
  };
  struct CompareTauReverse {
    bool operator()(const struct SegIndex& s1, const struct SegIndex& s2) {
      return s1.tau > s2.tau;
    }
  };

  int index;
  float tau;
  SEGMENT* pS;
} SEG_INDEX;

void TrackFinder::connectSegments() {

  float etaPreCut = 0.06;

  float minCurv = 0.3/110.0;

  int nBlocks  = N_THREADS;

  omp_set_num_threads(nBlocks);

  int tid;

  const PATTERN_CUTS* cuts = &m_cuts[m_stage];

  const float min_deta[4] = {1.0, expf(cuts->dEta_min_pix), expf(cuts->dEta_min_mix), expf(cuts->dEta_min_shs)};
  const float max_deta[4] = {1.0, expf(cuts->dEta_max_pix), expf(cuts->dEta_max_mix), expf(cuts->dEta_max_shs)};

  const float min_dphi[4] = {0.0, cuts->dPhi_min_pix, cuts->dPhi_min_mix, cuts->dPhi_min_shs};
  const float max_dphi[4] = {0.0, cuts->dPhi_max_pix, cuts->dPhi_max_mix, cuts->dPhi_max_shs};

  const float max_curv    = cuts->dCurv_max;

  if(m_stage == 1) {
    etaPreCut    = 0.06;
    minCurv = 0.3/120.0;
  }

  const float min_EtaPreCut = exp(-etaPreCut);
  const float max_EtaPreCut = exp( etaPreCut);


  //1. get all relevant nodes

  std::vector<NODE*> vNodes;

  m_hitBank->getConnectingNodes(vNodes);

  //std::cout<<"found "<<vNodes.size()<<" connections"<<std::endl;

  int nNodes = vNodes.size();

#pragma omp parallel private(tid) shared(cuts, nNodes, nBlocks)
  {

    tid = omp_get_thread_num();

    float minCurv2 = minCurv*minCurv;

    SEG_INDEX vSO[MAX_OUT_SEGS];
    SEG_INDEX vSI[MAX_OUT_SEGS];

    for(int nodeIdx=tid;nodeIdx<nNodes;nodeIdx+=nBlocks) {

      NODE* pN = vNodes.at(nodeIdx);

      //2. sort the segments to accelerate the inner loop

      int idxO=0;

      for(unsigned int k=0;k<pN->m_out.size();k++) {

	if(k>=MAX_OUT_SEGS) break;

	unsigned int nn = pN->m_out[k];
	unsigned int nextBankId = (nn & SEG_BANK_MASK) >> SEG_MASK_SHIFT;
	unsigned int nextSegmentIdx = (nn & SEG_INDEX_MASK);

	SEGMENT* pNS = &(m_segBank[nextBankId]->m_S[nextSegmentIdx]);

	vSO[idxO].index = k;
	vSO[idxO].tau = pNS->m_p[0];
	vSO[idxO].pS = pNS;
	idxO++;
	if(idxO==MAX_OUT_SEGS) break;
      }

      if(idxO==0) continue;

      std::sort(vSO, vSO+idxO, SEG_INDEX::CompareTau()); //sort segments by tau

      int idxI=0;

      for(unsigned int k=0;k<pN->m_in.size();k++) {

	if(k>=MAX_OUT_SEGS) break;

	unsigned int nn = pN->m_in[k];
	unsigned int nextBankId = (nn & SEG_BANK_MASK) >> SEG_MASK_SHIFT;
	unsigned int nextSegmentIdx = (nn & SEG_INDEX_MASK);

	SEGMENT* pNS = &(m_segBank[nextBankId]->m_S[nextSegmentIdx]);

	vSI[idxI].index = k;
	vSI[idxI].tau = pNS->m_p[0];
	vSI[idxI].pS = pNS;
	idxI++;
	if(idxI==MAX_OUT_SEGS) break;
      }

      if(idxI==0) continue;

      std::sort(vSI, vSI+idxI, SEG_INDEX::CompareTau()); //sort segments by tau

      //3. loop over the incoming segments

      int segStart = 0;
      int segEnd = 0;

      int v2 = m_pH->m_vol_id[pN->m_h];
      int nCells1 = pN->m_nCells;

      for(int inIdx = 0;inIdx<idxI;inIdx++) {

	SEGMENT* pS = vSI[inIdx].pS;
	float tau1 = vSI[inIdx].tau;
	float etaPreCutM = tau1*min_EtaPreCut;
	float etaPreCutP = tau1*max_EtaPreCut;

	while(vSO[segStart].tau < etaPreCutM && segStart<idxO) segStart++;
	segEnd = segStart;
	while(vSO[segEnd].tau < etaPreCutP && segEnd<idxO) segEnd++;

	if(segEnd - segStart < 1) continue;

	pS->m_nNei = 0;

	int v1 = m_pH->m_vol_id[pS->m_n2->m_h];

	bool isPix1 = (v1 == 8);

	bool pred1 = v2 < 10 && v1 < 10;
	bool pred2 = v2 > 10 && v1 > 10;

	float curv1 = pS->m_p[1];

	float tau_cut_min[4];
	float tau_cut_max[4];

	for(int i=1;i<4;i++) {
	  tau_cut_min[i] = tau1*min_deta[i];
	  tau_cut_max[i] = tau1*max_deta[i];
	}

	NODE& n2 = *pS->m_n1;
	NODE& n3 = *pS->m_n2;

	float dx = n3.m_x - n2.m_x;
	float dy = n3.m_y - n2.m_y;

	float kappa2Cut = minCurv2*(dx*dx+dy*dy);

	float curvCutM   = -max_curv + curv1;
	float curvCutP   =  max_curv + curv1;

	for(int outIdx=segStart;outIdx<segEnd;outIdx++) {

	  SEGMENT* pNS = vSO[outIdx].pS;

	  float curv2 = pNS->m_p[1];

	  if(curv2 < curvCutM || curv2 > curvCutP) continue;

	  bool isPix2 = m_pH->m_vol_id[pNS->m_n1->m_h] < 10;
	  NODE* n_2 = pNS->m_n1;

	  int nCells2 = n_2->m_nCells;

	  if(isPix1 && isPix2) {
	    if(abs(nCells2-nCells1) > 1) continue;
	  }

	  int code = 2;
	  if(pred1 && isPix2) code = 1;
	  else if(pred2 && !isPix2) code = 3;

	  float tau2 = vSO[outIdx].tau;

	  if(tau2 < tau_cut_min[code]) continue;
	  if(tau2 > tau_cut_max[code]) break;

	  float dphi_min = min_dphi[code];
	  float dphi_max = max_dphi[code];

	  float dPhi =  pNS->m_p[3] - pS->m_p[2];

	  if(dPhi<-M_PI) dPhi += 2*M_PI;
	  else if(dPhi>M_PI) dPhi -= 2*M_PI;

	  if(dPhi < dphi_min || dPhi > dphi_max) continue;

	  //curvature test

	  NODE& n1 = *pNS->m_n1;

	  float dy21 = n2.m_y-n1.m_y;
	  float dx21 = n2.m_x-n1.m_x;
	  float dy31 = n3.m_y-n1.m_y;
	  float dx31 = n3.m_x-n1.m_x;

	  float D = dy31*dx21-dx31*dy21;
	  float D2 = D*D;
	  float N = (dx21*dx21+dy21*dy21)*(dx31*dx31+dy31*dy31);

	  float kappa2 = D2/N;

	  if(kappa2>kappa2Cut) continue;

	  pS->m_vNei[pS->m_nNei++] = pN->m_out[vSO[outIdx].index];

	  if(pS->m_nNei >= N_SEG_CONNS) break;

	}

      }

    }
  }

}

void TrackFinder::createSegments(const LAYER_LINK& ll, int bankId) {

  float minCurv = 0.3/134.0;
  float phiWidth = 0.187;

  const PATTERN_CUTS& cuts = m_cuts[m_stage];

  if(m_stage == 2) {
    phiWidth = 0.299;
    minCurv = 0.3/68.0;
  }

  if(m_stage == 3) {
    phiWidth = 0.24;
    minCurv = 0.3/110.0;
  }

  if(m_stage == 4) {
    phiWidth = 0.297;
    minCurv = 0.3/24.0;
  }

  float minCurv2 = minCurv*minCurv;

  const float absTauCut_CL2[13]  = {0.0, 0.566629, 0.861533, 1.30254, 1.9043, 2.24961, 3.10129, 3.62686, 4.02186, 4.69117, 5.52207, 5.80969, 6.17407};
  const float absTauCut2[12] = {0.0,0.0,0.0,0.211547, 0.221779, 0.226903, 0.237169, 0.247458, 0.283673, 0.299297, 0.399962, 0.410752};

  const float PhiCut[5] = {0.0, 0.49, 0.51, 0.53, 0.4};

  int tableIndex = ll.m_index;

  //std::cout<<"Link index="<<tableIndex<<std::endl;

  std::vector<int>& bin2Start = m_binTable[m_stage][tableIndex].m_binStart;
  std::vector<int>& bin2End   = m_binTable[m_stage][tableIndex].m_binEnd;

  std::vector<std::vector<float> >& bin2Width = m_binTable[m_stage][tableIndex].m_phiWidth;

  HIT_LAYER* pL1 = m_hitBank->getLayer(ll.m_dst);
  HIT_LAYER* pL2 = m_hitBank->getLayer(ll.m_src);

  int vol1 = pL1->m_info.m_vol_id;
  int vol2 = pL2->m_info.m_vol_id;

  //assuming that n1 is always closer to the IP than n2

  bool isPixelBarrel  = vol1 == 8;
  bool isPixelBarrel2 = vol2 == 8;

  bool twoPix = vol1 == 8 && vol2 == 8;
  bool isPixelEndcap = vol1 == 7 || vol1 == 9;
  bool isPixelEndcap2 = vol2 == 7 || vol2 == 9;

  bool barrelSeeds = false;
  barrelSeeds = pL1->m_info.m_type == 0 && pL2->m_info.m_type == 0;

  float cutOnZ0 = barrelSeeds ? cuts.maxZ0_barr : cuts.maxZ0_all;

  //loop over bins/hits in L1

  for(int binIdx=0;binIdx<pL1->m_nBins;binIdx++) {

    int bin_start = bin2Start[binIdx];
    int bin_end   = bin2End[binIdx];

    if(bin_start == -1) continue;

    std::vector<float>& phiWidthVec = bin2Width[binIdx];

    const std::vector<NODE*>& vn1Good = pL1->m_bins[binIdx].m_l1nodes;

    if(vn1Good.empty()) continue;

    float binPhiWidth = phiWidth;

    if(twoPix) {
      binPhiWidth = m_stage == 1 ? 0.11 : 0.19;
    }

    std::vector<int> nSegs;
    nSegs.resize(vn1Good.size(), 0);

    int nb2 = bin_end-bin_start+1;
    if(m_stage == 1 && nb2 > 6) {
      bin_start+=1;
      bin_end-=1;
    }

    for(int bin2 = bin_start; bin2<=bin_end;bin2++) {

      HIT_BIN& B2 = pL2->m_bins[bin2];

      if(B2.m_nodes.empty()) continue;
      if(B2.m_l2phiIndices.empty()) continue;

      if(!twoPix && m_stage != 1) {
	binPhiWidth = 1.75*phiWidthVec.at(bin2-bin_start);
	if(binPhiWidth > phiWidth)
	  binPhiWidth = phiWidth;
      }

      //loop over hits in Bin2

      int j_start = 0;
      int j_end = j_start;

      for(unsigned int n1Idx=0;n1Idx<vn1Good.size();n1Idx++) {

	NODE* n1 = vn1Good[n1Idx];

	if(nSegs[n1Idx]>70) continue;

	float minTau = n1->m_minCutOnPar0;
	float maxTau = n1->m_maxCutOnPar0;

	float r1 = n1->m_r;
	float z1 = n1->m_z;

	if(minTau > -100.0) {

	  //check if bin2 compatible with the minimax "measurement"

	  float binTauMin, binTauMax;

	  pL2->getBinTaus(r1, z1, bin2, binTauMin, binTauMax);

	  if(binTauMax < minTau || binTauMin > maxTau) continue;
	}

	int nCellsT = n1->m_nCellsT;

	float Phi1 = n1->m_phi;

	//find compatible hits

	float minPhi = Phi1 - binPhiWidth;
	float maxPhi = Phi1 + binPhiWidth;

	B2.updatePhiRangeStart(minPhi, j_start);
	j_end = j_start;
	B2.updatePhiRangeEnd(maxPhi, j_end);

	if(j_end-j_start<1) continue;

	for(int phi_node_idx=j_start;phi_node_idx<j_end;phi_node_idx++) {//sliding window

	  NODE* n2 = B2.m_l2phiIndices[phi_node_idx].m_pN;

	  int nCellsT2 = n2->m_nCellsT;

	  if(isPixelEndcap && isPixelEndcap2) {
	    if(fabs(nCellsT2-nCellsT) > 3) continue;
	  }

	  float r2 = n2->m_r;

	  float dr = r2-r1;

	  if(dr<=2.5) continue;

	  float z2 = n2->m_z;
	  float dz = z2-z1;

	  float tau = dz/dr;

	  float absTau = fabs(tau);

	  if(absTau>29.0) continue;//detector acceptance

	  if(absTau < minTau) continue;
	  if(absTau > maxTau) continue;

	  float z0 = z1 - r1*tau;
	  if(fabs(z0) > cutOnZ0) continue;

	  if(isPixelEndcap2) {

	    int nCells2 = n2->m_nCells;
	    if(nCells2 == 1 && absTau < 2.376) continue;
	    if(nCells2 == 2) {
	      if(absTau > 30.1) continue;
	      if(absTau < 1.904) continue;
	    }
	  }

	  if(isPixelBarrel2) {
	    int nCells2 = n2->m_nCells;

	    if(nCells2 < 25 && absTau > 9.05956) continue;//barrel acceptance
	    if(nCells2<13 && absTau > absTauCut_CL2[nCells2]) continue;

	    if(nCells2<11) {
	      if(absTau < absTauCut2[nCells2]-0.015) continue;
	    }
	    else {
	      if(absTau < absTauCut2[11]+0.01) continue;
	    }

	  }
	  float dx = n2->m_x - n1->m_x;
	  float dy = n2->m_y - n1->m_y;
	  float L2 = 1/(dx*dx+dy*dy);

	  float D = (n2->m_y*n1->m_x - n1->m_y*n2->m_x)/(r1*r2);

	  float kappa = D*D*L2;

	  if(kappa>minCurv2) continue;

	  float curv = D*sqrt(L2);//curvature

	  //float dPhi2 = asinf(curv*r2);

	  float df = curv*r2;
	  float df2 = df*df;
	  float dPhi2 = df*(1 + df2*(0.1667 + 0.075*df2));//asinf

	  if(isPixelBarrel2) {
	    if(fabs(dPhi2) > PhiCut[nCellsT2]) continue;
	  }

	  //float dPhi1 = asinf(curv*r1);

	  float ef = curv*r1;
	  float ef2 = ef*ef;
	  float dPhi1 = ef*(1 + ef2*(0.1667 + 0.075*ef2));//asinf

	  if(isPixelBarrel) {
	    if(fabs(dPhi1) > PhiCut[nCellsT]) continue;
	  }

	  float Phi2 = n2->m_phi;

	  int segIndex  = m_segBank[bankId]->m_nSegments;
	  float* params = m_segBank[bankId]->m_S[segIndex].m_p;//exp(-eta), curvature, phi1, phi2
	  params[0] = sqrt(1+tau*tau)-tau;
	  params[1] = curv;
	  params[2] = Phi1 + dPhi1;
	  params[3] = Phi2 + dPhi2;


	  m_segBank[bankId]->m_S[segIndex].initialize(n1, n2);

	  if(m_segBank[bankId]->m_nSegments<N_SEGMENTS_PER_BANK-1) {
	    m_segBank[bankId]->m_nSegments++;
	    nSegs[n1Idx]++;
	  }
	}
      }
    }
  }

}


BIN_TABLE* TrackFinder::createBinTable(int cut_idx) {

  float maxD0 = 30.0;

  if(cut_idx == 1) maxD0 = 10.0;
  if(cut_idx == 2) maxD0 = 130.0;
  if(cut_idx == 3) maxD0 = 10.0;
  if(cut_idx == 4) maxD0 = 140.0;

  float zvMax = m_cuts[cut_idx].zvMax;
  float gammaMS = m_cuts[cut_idx].gammaMS;

  const std::vector<LAYER_LINK>& LL = m_linker.m_links;

  BIN_TABLE* binTable = new BIN_TABLE[LL.size()];

  //std::ofstream binPhi("phi_table.txt");

  for(std::vector<LAYER_LINK>::const_iterator llIt=LL.begin();llIt!=LL.end();++llIt) {

    int tableIndex = (*llIt).m_index;

    //std::cout<<"creating binTable for index "<<tableIndex<<" stage "<<cut_idx<<std::endl;

    HIT_LAYER* pL1 = m_hitBank->getLayer((*llIt).m_dst);
    HIT_LAYER* pL2 = m_hitBank->getLayer((*llIt).m_src);

    binTable[tableIndex].m_binStart.resize(pL1->m_nBins,-1);
    binTable[tableIndex].m_binEnd.resize(pL1->m_nBins,-1);
    binTable[tableIndex].m_phiWidth.resize(pL1->m_nBins);

    for(int binIdx=0;binIdx<pL1->m_nBins;binIdx++) {

      float r_b1, z_b1;//bin1 coordinates

      if(pL1->m_info.m_type == 0) {//barrel
	z_b1 = pL1->m_binCoords[binIdx];
	r_b1 = pL1->m_info.m_refCoord;
      }
      else {//disk
	r_b1 = pL1->m_binCoords[binIdx];
	z_b1 = pL1->m_info.m_refCoord;
      }

      float cp, cm, path;

      if(pL2->m_info.m_type == 0) {//barrel

	float r_2 = pL2->m_info.m_refCoord;
	path = r_2 - r_b1;
	float rr = r_2/r_b1;

	cp = z_b1*rr - zvMax*(1-rr) + gammaMS*path;
	cm = z_b1*rr + zvMax*(1-rr) - gammaMS*path;

      } else {//disk

	float z_2 = pL2->m_info.m_refCoord;
	path = z_2-z_b1;

	cp = r_b1*(z_2-zvMax)/(z_b1-zvMax) + gammaMS*path;
	cm = r_b1*(z_2+zvMax)/(z_b1+zvMax) - gammaMS*path;
      }


      float refCoord2 = pL2->m_info.m_refCoord;

      if(pL2->m_info.m_type == 0) {//barrel

	cm = cm/refCoord2;
	cp = cp/refCoord2;
      }
      else {
	cm = refCoord2/cm;
	cp = refCoord2/cp;
      }

      if(cp<cm) {
	std::swap(cm, cp);
      }

      int bin_start, bin_end;

      bool valid = pL2->getBinRange(cm, cp, bin_start, bin_end);

      if(valid && bin_end>=bin_start) {

	binTable[tableIndex].m_binStart[binIdx] = bin_start;
	binTable[tableIndex].m_binEnd[binIdx] = bin_end;

	for(int idx=bin_start;idx<=bin_end;idx++) {
	  float r_b2 = pL2->getBinRadius(idx);
	  float cosPhi = (r_b1*r_b1 - r_b2*maxD0)/(r_b1*(r_b2-maxD0));
	  float Phi=100.0;
	  if(cosPhi<1.0 && cosPhi>0.0) {
	    Phi = acos(cosPhi);
	  }
	  //binPhi<<pL1->m_info.m_vol_id<<":"<<pL1->m_info.m_lay_id<<" -> "<<pL2->m_info.m_vol_id<<":"<<pL2->m_info.m_lay_id<<" ";
	  //binPhi<<binIdx<<" "<<idx<<" "<<Phi<<std::endl;
	  binTable[tableIndex].m_phiWidth[binIdx].push_back(Phi);
	}

      } else {

	binTable[tableIndex].m_binStart[binIdx] = -1;
	binTable[tableIndex].m_binEnd[binIdx] = -1;

      }

    }
  }
  //binPhi.close();

  return binTable;
}

void TrackFinder::assignCloneFlags(std::vector<NEW_TRACK*>& vTracks) {

  float shareThreshold = m_params.sf_f[5];
  float addThreshold   = m_params.sf_f[6];

  //assign new labels

  unsigned int trackId = 1;

  for(std::vector<NEW_TRACK*>::iterator it=vTracks.begin();it!=vTracks.end();++it,trackId++) {
    (*it)->m_trackId = trackId;
    (*it)->m_cloneFlag = 0;//accept by default
  }

  int maxHit = m_pH->m_nHits;
  memset(m_H2T, 0, (maxHit+1)*sizeof(int));


  for(std::vector<NEW_TRACK*>::iterator it = vTracks.begin();it!=vTracks.end();++it) {

    for(int k=0;k<(*it)->m_nHits;k++) {

      int h = (*it)->m_hits[k];

      unsigned int hit_id = h + 1;

      int tid = m_H2T[hit_id];

      if(tid == 0 || tid > (*it)->m_trackId) {
	m_H2T[hit_id] = (*it)->m_trackId;
      }
    }
  }

  int nTracks = vTracks.size();

  int nBlocks  = N_THREADS;

  omp_set_num_threads(nBlocks);

  int tid;

#pragma omp parallel private(tid) shared(nTracks, vTracks, nBlocks)
  {

    tid = omp_get_thread_num();

    for(int trackIdx = tid; trackIdx<nTracks; trackIdx += nBlocks) {

      NEW_TRACK* pT = vTracks.at(trackIdx);

      int nTotal = 0;
      int nOther = 0;
      int trackId = pT->m_trackId;

      std::map<int, int> aM;//hit asssociation map

      for(int k=0;k<pT->m_nHits;k++) {
	nTotal++;
	int h = pT->m_hits[k];
	unsigned int hit_id = h + 1;
	int tid = m_H2T[hit_id];
	if(tid != trackId) {
	  nOther++;
	  std::map<int, int>::iterator aIt = aM.find(tid);
	  if(aIt == aM.end()) {
	    aM.insert(std::pair<int, int>(tid, 1));
	  } else (*aIt).second++;
	}
      }

      int nContrib = (int)aM.size();

      if(nContrib >= 2) {
	pT->m_cloneFlag = -1;//reject
	continue;
      }

      if(nOther > shareThreshold*nTotal) {


	int nLargest = 0;
	int largest_tid = -1;

	for(std::map<int, int>::iterator aIt = aM.begin();aIt!=aM.end();++aIt) {
	  if((*aIt).second > nLargest) {
	    nLargest = (*aIt).second;
	    largest_tid = (*aIt).first;
	  }
	}
	pT->m_cloneFlag = -1;//reject

	int n1 = vTracks.at(largest_tid-1)->m_nHits;
	int n2 = nTotal;//this track
	int ns = nLargest;
	int nAdd = n2-ns;

	if(nContrib == 1) {
	  if(ns >= 2 && nAdd < addThreshold*n1) pT->m_cloneFlag = largest_tid;//merge with this track
	}
      }
    }
  }

  //  delete[] H2T;
}

void TrackFinder::groupTracks(std::vector<NEW_TRACK*>& vTracks) {

  //const float llWeight = 19.0;//was 20.0

  int nHits = m_pH->m_nHits;

  std::vector<TRACK_GROUP*> allGroups;

  createGroups(nHits, vTracks, allGroups);

  std::sort(allGroups.begin(), allGroups.end(), TRACK_GROUP::CompareSize());

  for(std::vector<TRACK_GROUP*>::iterator it=allGroups.begin();it!=allGroups.end();++it) {

    if((*it)->m_size == 1) continue;

    std::sort((*it)->m_vTracks.begin(), (*it)->m_vTracks.end(), NEW_TRACK::CompareLikelihood());

  }

  vTracks.clear();

  //for(std::vector<TRACK_GROUP*>::reverse_iterator it=allGroups.rbegin();it!=allGroups.rend();++it) {
  //  std::copy((*it)->m_vTracks.begin(), (*it)->m_vTracks.end(), std::back_inserter(vTracks));
  //}

  for(std::vector<TRACK_GROUP*>::iterator it=allGroups.begin();it!=allGroups.end();++it) {
    std::copy((*it)->m_vTracks.begin(), (*it)->m_vTracks.end(), std::back_inserter(vTracks));
  }


  for(std::vector<TRACK_GROUP*>::iterator it=allGroups.begin();it!=allGroups.end();++it) {
    delete *it;
  }

}
void TrackFinder::createGroups(int nHits, std::vector<NEW_TRACK*>& vTracks, std::vector<TRACK_GROUP*>& vG) {

  int* Usage = new int[nHits+1];//hit id starts with 1
  memset(Usage, 0, (nHits+1)*sizeof(int));

  //1. assign unique track ids

  int trackId = 1;

  for(std::vector<NEW_TRACK*>::iterator it = vTracks.begin();it!=vTracks.end();++it, trackId++) {
    (*it)->m_trackId = trackId;
  }

  //2. iterate

  for(int iter=0;iter<100;iter++) {

    int nFlops = 0;

    for(std::vector<NEW_TRACK*>::iterator it = vTracks.begin();it!=vTracks.end();++it) {

      int newTrackId = (*it)->m_trackId;

      for(int k=0;k<(*it)->m_nHits;k++) {

	int h = (*it)->m_hits[k];

	unsigned int hit_id = m_pH->m_hit_index[h];
	int tid = Usage[hit_id];
	if(tid == 0) {
	  Usage[hit_id] = (*it)->m_trackId;
	  continue;
	}
	else {
	  if(tid>=(*it)->m_trackId) Usage[hit_id] = (*it)->m_trackId;
	  else {
	    if(newTrackId > tid) newTrackId = tid;
	  }
	}
      }

      if(newTrackId != (*it)->m_trackId) nFlops++;
      (*it)->m_trackId = newTrackId;
    }
    if(nFlops==0) break;
  }

  //counting groups

  std::set<int> groupSet;

  for(std::vector<NEW_TRACK*>::iterator it = vTracks.begin();it!=vTracks.end();++it) {
    int tid = (*it)->m_trackId;
    if(groupSet.find(tid) == groupSet.end()) groupSet.insert(tid);
  }

  //  std::cout<<"Found "<<groupSet.size()<<" disjoint track groups, nTracks="<<vTracks.size()<<std::endl;

  std::map<int, std::vector<NEW_TRACK*> > gMap;

  for(std::vector<NEW_TRACK*>::iterator it = vTracks.begin();it!=vTracks.end();++it) {
    int tid = (*it)->m_trackId;
    std::map<int, std::vector<NEW_TRACK*> >::iterator gIt = gMap.find(tid);
    if(gIt==gMap.end()) {
      std::vector<NEW_TRACK*> v(1, *it);
      gMap.insert(std::pair<int, std::vector<NEW_TRACK*> >(tid, v));
    } else (*gIt).second.push_back(*it);
  }

  for(std::map<int, std::vector<NEW_TRACK*> >::iterator gIt = gMap.begin();gIt!=gMap.end();++gIt) {

    TRACK_GROUP* pG = new TRACK_GROUP((*gIt).first, (*gIt).second);
    vG.push_back(pG);

  }

  delete[] Usage;

}

void TrackFinder::maskHits(std::vector<NEW_TRACK*>& vTracks, int threshold) {

  //create unique hit-to-track association

  memset(m_H2T, 0, (m_pH->m_nHits+1)*sizeof(int));

  int trackId = 1;

  std::vector<NEW_TRACK*> vTmp;

  for(std::vector<NEW_TRACK*>::iterator it = vTracks.begin();it!=vTracks.end();++it, trackId++) {

    //if(!(*it)->valid()) continue;

    std::vector<int> vHits;//uniquely associated to this track

    for(int k=0;k<(*it)->m_nHits;k++) {

      int h = (*it)->m_hits[k];

      unsigned int hit_id = m_pH->m_hit_index[h];

      int tid = m_H2T[hit_id];
      if(tid == 0) {
	m_H2T[hit_id] = trackId;
	vHits.push_back(h);
      }
    }

    if((int)vHits.size() > threshold) {
      //mask
      vTmp.push_back(*it);

      if(m_stage < 4) {

	for(std::vector<int>::iterator hIt = vHits.begin();hIt!=vHits.end();++hIt) {
	  m_pH->m_mask[*hIt] = 1;
	}
      }
    }
  }

  vTracks = std::move(vTmp);
}

void TrackFinder::initializeCuts() {

  m_cuts[1].zvMax = 150.0;//was 150.0
  m_cuts[1].gammaMS = 0.16;//0.16

  m_cuts[2].zvMax = 92.0;//was 92
  m_cuts[2].gammaMS = 0.18;//was 0.18

  m_cuts[3].zvMax = 100.0;//was 92
  m_cuts[3].gammaMS = 0.22;//was 0.18

  m_cuts[4].zvMax = 140.0;//was 140
  m_cuts[4].gammaMS = 0.155;//was 0.155

  m_cuts[1].maxZ0_barr = 185.0;//was 185.0
  m_cuts[1].maxZ0_all  = 230.0;//was 230

  m_cuts[2].maxZ0_barr = 190.0;//was 180.0
  m_cuts[2].maxZ0_all  = 200.0;//was 200.0

  m_cuts[3].maxZ0_barr = 180.0;//was 180.0
  m_cuts[3].maxZ0_all  = 300.0;//was 100.0

  m_cuts[4].maxZ0_barr = 180.0;//was 180.0
  m_cuts[4].maxZ0_all  = 900.0;//was 900.0

  m_cuts[1].dEta_min_pix = -0.031;//was -0.031
  m_cuts[1].dEta_max_pix = 0.030;//was 0.030

  m_cuts[1].dEta_min_shs = -0.058;//was -0.058
  m_cuts[1].dEta_max_shs = 0.060;//was 0.058

  m_cuts[1].dEta_min_mix = -0.049;//was -0.049
  m_cuts[1].dEta_max_mix = 0.051;//was 0.051

  m_cuts[1].dPhi_min_pix = -0.028;//was -0.029
  m_cuts[1].dPhi_max_pix = 0.0295;//was 0.029

  m_cuts[1].dPhi_min_shs = -0.044;//was -0.044-
  m_cuts[1].dPhi_max_shs = 0.043;//was 0.043

  m_cuts[1].dPhi_min_mix = -0.026;//was -0.026
  m_cuts[1].dPhi_max_mix = 0.026;//was 0.026

  m_cuts[1].dCurv_max = 0.00075;//was 0.0008

  m_cuts[2].dEta_min_pix = -0.039;//was -0.034
  m_cuts[2].dEta_max_pix = 0.039;//was 0.036

  m_cuts[2].dEta_min_shs = -0.063;//was -0.061
  m_cuts[2].dEta_max_shs = 0.061;//was 0.061

  m_cuts[2].dEta_min_mix = -0.060;//was -0.051
  m_cuts[2].dEta_max_mix = 0.051;//was 0.052

  m_cuts[2].dPhi_min_pix = -0.092;//was -0.092
  m_cuts[2].dPhi_max_pix = 0.092;//was 0.092

  m_cuts[2].dPhi_min_shs = -0.075;//was -0.075
  m_cuts[2].dPhi_max_shs = 0.075;//was 0.075

  m_cuts[2].dPhi_min_mix = -0.088;//was -0.088
  m_cuts[2].dPhi_max_mix = 0.088;//was 0.088

  m_cuts[2].dCurv_max = 0.0020;//was 0.0035

  m_cuts[3].dEta_min_pix = -0.038;//was -0.038
  m_cuts[3].dEta_max_pix = 0.038;//was 0.038

  m_cuts[3].dEta_min_shs = -0.063;//was -0.063
  m_cuts[3].dEta_max_shs = 0.061;//was 0.061

  m_cuts[3].dEta_min_mix = -0.053;//was -0.051
  m_cuts[3].dEta_max_mix = 0.052;//was 0.052

  m_cuts[3].dPhi_min_pix = -0.095;//was -0.092
  m_cuts[3].dPhi_max_pix = 0.095;//was 0.092

  m_cuts[3].dPhi_min_shs = -0.077;//was -0.077
  m_cuts[3].dPhi_max_shs = 0.077;//was 0.077

  m_cuts[3].dPhi_min_mix = -0.089;//was -0.089
  m_cuts[3].dPhi_max_mix = 0.089;//was 0.089

  m_cuts[3].dCurv_max = 0.0006;//was 0.003


  m_cuts[4].dEta_min_pix = -0.034;//was -0.034
  m_cuts[4].dEta_max_pix = 0.033;//was 0.033

  m_cuts[4].dEta_min_shs = -0.056;//was -0.057
  m_cuts[4].dEta_max_shs = 0.057;//was 0.059

  m_cuts[4].dEta_min_mix = -0.049;//was -0.049
  m_cuts[4].dEta_max_mix = 0.052;//was 0.051

  m_cuts[4].dPhi_min_pix = -0.192;//was -0.191
  m_cuts[4].dPhi_max_pix = 0.192;//was 0.191

  m_cuts[4].dPhi_min_shs = -0.084;//was -0.084
  m_cuts[4].dPhi_max_shs = 0.084;//was 0.084

  m_cuts[4].dPhi_min_mix = -0.081;//was -0.081
  m_cuts[4].dPhi_max_mix = 0.081;//was 0.081

  m_cuts[4].dCurv_max = 0.004;//was 0.004

}
void TrackFinder::train(const HITS* pE, const std::map<unsigned long, PARTICLE>&) {

  std::cout<<"In train ..."<<std::endl;

  m_pH = pE;

  m_hitBank->importHits(pE);

  m_hitBank->sort();

  std::vector<HIT_LAYER*> vL, vTmp;

  m_hitBank->getLayers(vTmp);

  for(std::vector<HIT_LAYER*>::iterator lIt = vTmp.begin();lIt!=vTmp.end();++lIt) {

    HIT_LAYER* pL = (*lIt);

    if(pL->m_info.m_vol_id >= 16) continue;
    vL.push_back(pL);
  }

  int nLayers = vL.size();

  int layStart = 0;
  int layEnd = nLayers;

  CLUSTER_MAKER cm(m_pH, m_geo);

  for(int layIdx=layStart;layIdx<layEnd;layIdx++) {

    HIT_LAYER* pL = vL.at(layIdx);

    cm.findClusters(m_params, pL);

    //evaluate clusters

    int nTotal = 0;
    int nBad = 0;

    std::map<unsigned long, std::set<int> > partMap;

    int nodeIdx=0;

    for(int binIdx=0;binIdx<pL->m_nBins;binIdx++) {

      HIT_BIN& B = pL->m_bins[binIdx];

      for(unsigned int n1=0;n1<B.m_nodes.size();n1++, nodeIdx++) {

	std::vector<int>& vHits = B.m_nodes[n1]->m_hits;

	if(vHits.size() > 1) nTotal++;

	std::set<unsigned long> partSet;

	for(unsigned int k=0;k<vHits.size();k++) {

	  int h1 = vHits.at(k);

	  //float phi1 = B.m_nodes[n1]->m_phi;
	  //float r1 = B.m_nodes[n1]->m_r;

	  float w1 = m_pH->m_weight[h1];

	  if(w1==0.0) continue;

	  unsigned long part1 = m_pH->m_particle[h1];

	  partSet.insert(part1);

	  if(part1 == 0) continue;

	  if(partMap.find(part1) == partMap.end()) {
	    std::set<int> s;
	    s.insert(nodeIdx);
	    partMap.insert(std::pair<unsigned long, std::set<int> >(part1, s));
	  } else (*partMap.find(part1)).second.insert(nodeIdx);

	}

	if(partSet.size() > 1) nBad++;

      }
    }

    int nSplit = 0;

    for(std::map<unsigned long, std::set<int> >::iterator it = partMap.begin();it!=partMap.end();++it) {
      if((*it).second.size() > 1) nSplit++;
    }

    std::cout<<"L"<<layIdx<<" nBad="<<nBad<<" nTotal="<<nTotal<<" nSplit="<<nSplit<<" "<<partMap.size()<<std::endl;
  }

  //int nLeft = m_hitBank->numberOfHits();

  m_hitBank->clear();

}
