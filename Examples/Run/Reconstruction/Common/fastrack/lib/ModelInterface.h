#ifndef __MODEL_INTERFACE_H__
#define __MODEL_INTERFACE_H__

class ModelInterface {
 public:
  ModelInterface() {};
  virtual ~ModelInterface() {};

  virtual void importHits(int, const int*, const float*, const float*, const float*, const int*, const int*, const int*) = 0;
  //virtual void importCells(int, const int*, const int*, const int*, const float*) = 0;

  virtual void importCells(int, const int*, const int*, const int*) = 0;
 
  virtual void importParticles(int, const unsigned long*, const int*, 
			       const float*, const float*, const float*,  
			       const float*, const float*, const float*, 
			       const int*, const int*) = 0;

  virtual void importTruth(int, const int*, const float*, const float*, const float*, const float*, const float*, const float*, const unsigned long*, const float*) = 0;

  virtual void findTracks(int*) = 0; 

  virtual void train() = 0;
 
};

#endif

