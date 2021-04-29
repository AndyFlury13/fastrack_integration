#include "ModelInterface.h"
#include "ModelClass.h"

extern "C" {

  ModelInterface* createModel(const char* s1, int n1, const char* s2, int n2) {
    ModelInterface* pModel = new ModelClass(s1, n1, s2, n2);
    return pModel;
  }

  void deleteModel(ModelInterface* t)  {
    delete t;
  }
  
  void importHits(ModelInterface* t, int n, int* id, float* x, float* y, float* z, int* v, int* l, int* m) {
    t->importHits(n, id, x, y, z, v, l, m);
  }
  /*
  void importCells(ModelInterface* t, int n, int* id, int* ch1, int* ch2, float* v) {
    t->importCells(n, id, ch1, ch2, v);
  }
  */
  void importCells(ModelInterface* t, int n, int* id, int* ch1, int* ch2) {
    t->importCells(n, id, ch1, ch2);
  }
  
  void importParticles(ModelInterface* t, int n, unsigned long* pid, int* type, 
				 float* vx, float* vy, float* vz, 
				 float* px, float* py, float* pz, 
		       int* q, int* nhits) {
    t->importParticles(n, pid, type, vx, vy, vz, px, py, pz, q, nhits);
  }

  void importTruth(ModelInterface* t, int n, int* hit_id, float* x, float* y, float* z, float* px, float* py, float* pz, unsigned long* particle, float* w) {

    t->importTruth(n, hit_id, x, y, z, px, py, pz, particle, w);
  }

  void findTracks(ModelInterface* t, int* trackId) {
    t->findTracks(trackId);
  }

  void train(ModelInterface* t) {
    t->train();
  }

}
