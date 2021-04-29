#ifndef __TRACK_FINDER_H__
#define __TRACK_FINDER_H__

#include<vector>
#include<set>
#include<map>

#include "AlgorithmParameters.h"
#include "Event.h"
#include "Geometry.h"
#include "LayerLinker.h"
#include "RecoEvent.h"
#include "DataStructures.h"
#include "PatternCuts.h"


typedef class TrackFinder {
 public:
  TrackFinder(const ALGORITHM_PARAMETERS&, const GEOMETRY&, const LAYER_LINKER&);
  ~TrackFinder(); 

  void prepare(const HITS*);
  void clean();

  void run(int, const HITS*, std::vector<NEW_TRACK*>&, int);

  void train(const HITS*, const std::map<unsigned long, PARTICLE>&);

  void groupTracks(std::vector<NEW_TRACK*>&);

 protected:

  BIN_TABLE* createBinTable(int);
  
  void findClusters();
  
  void createSegments();
  void connectSegments();

  void runNetworkEvolution();

  void collectTracks(std::vector<NEW_TRACK*>&);
  
  void prefitTracks(std::vector<NEW_TRACK*>&);
  
  void fitTracks(std::vector<NEW_TRACK*>&, int);
  void removeClones(std::vector<NEW_TRACK*>&);

  void reAssignHits(std::vector<NEW_TRACK*>&);

  void runTrackExtensions(std::vector<NEW_TRACK*>&);

  void maskHits(std::vector<NEW_TRACK*>&, int);
  
  void createSegments(const LAYER_LINK&, int);
 
  void assignCloneFlags(std::vector<NEW_TRACK*>&);

  void createGroups(int, std::vector<NEW_TRACK*>&, std::vector<TRACK_GROUP*>&);

  void prepareL1Nodes();
  void createLinks();
  
  void initializeCuts();

  const ALGORITHM_PARAMETERS& m_params;
  
  const GEOMETRY& m_geo;
  const LAYER_LINKER& m_linker;
  
  HIT_BANK* m_hitBank;

  NEW_SEG_BANK* m_segBank[N_SEG_BANKS];

  const HITS* m_pH;

  int m_stage;

  PATTERN_CUTS m_cuts[5];

  BIN_TABLE* m_binTable[5];

  int* m_H2T;

  NEW_TRACK* m_trackBank;

  int m_globalTrackCounter;

  int m_eventCounter;

  double m_time[25];

  float m_segStat[100];
  unsigned int m_maxLink;
  
} TRACK_FINDER;

#endif
