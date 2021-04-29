#ifndef __MODEL_CLASS_H__
#define __MODEL_CLASS_H__

#include "ModelInterface.h"
#include "AlgorithmParameters.h"
#include "Event.h"
#include "Geometry.h"
#include "LayerLinker.h"
#include "TrackFinder.h"

#include<map>

class ModelClass : public ModelInterface {

 public:
  ModelClass(const char*, int, const char*, int);
  virtual ~ModelClass();

  virtual void importHits(int, const int*, const float*, const float*, const float*, const int*, const int*, const int*);
  //virtual void importCells(int, const int*, const int*, const int*, const float*);

  virtual void importCells(int, const int*, const int*, const int*);
  
  virtual void importParticles(int, const unsigned long*, const int*, 
			       const float*, const float*, const float*,  
			       const float*, const float*, const float*, 
			       const int*, const int*);

  virtual void importTruth(int, const int*, const float*, const float*, const float*, const float*, const float*, const float*, const unsigned long*, const float*);



  virtual void findTracks(int*);
  virtual void train();
  
 protected:

  void generateHitAssignment(std::vector<NEW_TRACK*>& vTracks, int*);
  
  GEOMETRY* m_geo;
  LAYER_LINKER* m_link;
  ALGORITHM_PARAMETERS* m_params;
  
  TRACK_FINDER* m_tf;
  
  HITS m_H;

  int m_eventCounter;
  
  std::map<unsigned long, PARTICLE> m_P;
};

#endif
