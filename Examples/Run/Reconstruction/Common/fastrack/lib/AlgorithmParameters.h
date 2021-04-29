#ifndef __ALGORITHM_PARAMETERS_H__
#define __ALGORITHM_PARAMETERS_H__

typedef struct AlgoParams {
public:

  AlgoParams();

  //1. ClusterMaker

  float  cm_f[17];  

  //2. Magnetic field model

  float  mf_f[3];

  //3. Track extrapolator

  int    te_i[5];
  float  te_f[25];

  //4. Track fitter

  int    tf_i[5];
  float  tf_f[60];
  double tf_d[20];
  float  tf_table1[9];


  //5. Track following filter
  
  float ou_d[20];

  //6. Track finder settings

  float sf_f[20];

} ALGORITHM_PARAMETERS;



#endif
