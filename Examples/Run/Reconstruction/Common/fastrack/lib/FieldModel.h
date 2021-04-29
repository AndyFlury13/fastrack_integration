#ifndef __FIELD_MODEL__
#define __FIELD_MODEL__

#include "AlgorithmParameters.h"

typedef class FieldModel {
 public:
  FieldModel();
  ~FieldModel();

  void initialize(const ALGORITHM_PARAMETERS&);

  void getField(const float*, float*) const;
  void getFieldFast(const float*, float*) const;
  
 protected:

  float m_C1, m_dC, m_C1_fast, m_dC_fast;
  
} FIELD_MODEL;

#endif
