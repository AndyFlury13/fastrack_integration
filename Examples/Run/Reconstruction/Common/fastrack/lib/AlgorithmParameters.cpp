#include "AlgorithmParameters.h"

AlgoParams::AlgoParams() {
  
  //1. ClusterMaker

  cm_f[ 0] = 0.5;
  cm_f[ 1] = 0.49;
  cm_f[ 2] = 200.0;
  cm_f[ 3] = 90.0;
  cm_f[ 4] = 0.23;
  cm_f[ 5] = 120.0;
  cm_f[ 6] = 0.68;
  cm_f[ 7] = 125.0;
  cm_f[ 8] = 1.04;
  cm_f[ 9] = 100.0;
  cm_f[10] = 1.26;
  cm_f[11] = 135.0;
  cm_f[12] = 190.0;
  cm_f[13] = 0.3;
  cm_f[14] = 1.6;
  cm_f[15] = 1.6;
  cm_f[16] = 150.0;
  
  //2. Magnetic field model
  
  mf_f[0] = 0.301;
  mf_f[1] = 2200.0;
  mf_f[2] = 21.95;
  
  //3. Track extrapolator 
  
  //3a. detector material model parameters:
  
  te_f[0] = 0.0092;
  te_f[1] = 0.0305;
  te_f[2] = 0.12;
  te_f[3] = 0.21;
  
  //3b. steps and detector boundaries
  
  te_i[ 0] = 60;
  te_i[ 1] = 80;
  
  te_f[ 4] = 150.0;//was 150.0
  te_f[ 5] = 0.8;
  te_f[ 6] = 0.085;
  te_f[ 7] = 1200.0;
  te_f[ 8] = 3200.0;
  te_f[ 9] = -22.0;//was -22.0
  te_f[10] = 20.0;
  
  //3c. Runge-Kutta parameters
  
  te_f[11] = 1.0;//was 1.0
  te_f[12] = 20.0;//was 10.0
  te_f[13] = 2000.0;
  te_f[14] = 15.0;
  te_f[15] = 300.0;
  
  //4. Track fitter
  
  tf_f[0] = 1.1179;
  
  //4b. Prefitter maxDChi2, maxFrac of outliers, likelihood layer weight
  
  tf_f[ 1] = 22.1;
  tf_f[ 2] = 0.39;
  tf_f[ 3] = 15.0;
  
  tf_f[ 4] = 25.0;
  tf_f[ 5] = 0.38;
  tf_f[ 6] = 17.0;
  
  tf_f[ 7] = 24.0;
  tf_f[ 8] = 0.39;
  tf_f[ 9] = 18.0;
  
  tf_f[10] = 26.1;
  tf_f[11] = 0.39;
  tf_f[12] = 18.0;
  
  tf_f[13] = 26.5;
  tf_f[14] = 0.42;
  tf_f[15] = 22.0;

  //4c. Fitter maxDChi2, likelihood layer weight, maxFrac of outliers

  tf_f[16] = 25.0;
  tf_f[17] = 11.0;
  tf_f[18] = 0.2;
  
  tf_f[19] = 22.0;
  tf_f[20] = 12.0;
  tf_f[21] = 0.1;
  
  tf_f[22] = 20.0;
  tf_f[23] = 12.0;
  tf_f[24] = 0.15;
  
  tf_f[25] = 26.0;
  tf_f[26] = 12.0;
  tf_f[27] = 0.1;
  
  tf_f[28] = 480.0;
  tf_f[29] = 35.0;
  tf_f[30] = 0.28;
  
  tf_f[31] = 20.5;
  tf_f[32] = 10.0;
  tf_f[33] = 0.05;
  
  tf_f[34] = 25.0;
  tf_f[35] = 10.0;
  tf_f[36] = 0.05;
  
  tf_f[37] = 21.0;
  tf_f[38] = 12.0;
  tf_f[39] = 0.05;

  //4d. Hit pick up stage
  
  tf_f[40] = 48.0;//was 44.0
  tf_f[41] = 70.0;//was 40.0
  tf_f[42] = 0.34;
  
  //4e. Backward propagation and hit search

  tf_f[43] = 50.0;
  tf_f[44] = 0.35;
  
  tf_f[45] = -17.0;//hole on track penalty
  tf_f[46] = 50.0;//hit on track premium
  tf_f[47] = 44.0;//minimum acceptable quality of a branch
  tf_f[48] = 23.5;// maxDChi2 28.5
  
  //4f. Forward hit search
  
  tf_f[49] = -65.0;//hole on track penalty
  tf_f[50] = 112.0;//hit on track premium
  tf_f[51] = 280.0;// maxDChi2
  tf_f[52] = -50.0;//minimum acceptable quality of a branch
  
  tf_i[0] = 1;//max num of propagated branches : backward 1
  tf_i[1] = 3;//max num of propagated branches : forward 3

  //4g. Measurement errors and cluster_length(theta) calibration

  const double sigmaU_pix = 0.082;
  const double sigmaV_pix = 0.155;
  const double sigmaU_shs = 0.118;
  const double sigmaV_shs = 0.69;
  const double sigmaU_lgs = 0.085;
  const double sigmaV_lgs = 4.2;
  
  tf_d[0] = sigmaU_pix*sigmaU_pix;
  tf_d[1] = sigmaV_pix*sigmaV_pix;
  tf_d[2] = sigmaU_shs*sigmaU_shs;
  tf_d[3] = sigmaV_shs*sigmaV_shs;
  tf_d[4] = sigmaU_lgs*sigmaU_lgs;
  tf_d[5] = sigmaV_lgs*sigmaV_lgs;

  float cellCut[9] = {0.0, 0.42, 0.69, 0.905, 0.98, 1.12, 1.22, 1.34, 1.475};

  for(int i=0;i<9;i++) tf_table1[i] = cellCut[i];

  //5. Track following filter    

  ou_d[ 0] = 2.9;
  ou_d[ 1] = 0.0195;
  ou_d[ 2] = 0.01;
  ou_d[ 3] = 0.50;
  ou_d[ 4] = 0.002;
  ou_d[ 5] = 0.0006;
  ou_d[ 6] = 0.00008;
  ou_d[ 7] = 0.008;
  ou_d[ 8] = 4.2e-8;
  ou_d[ 9] = 0.9e-4;
  ou_d[10] = 0.5;
  ou_d[11] = 0.5;
  ou_d[12] = 11.0;
  ou_d[13] = 0.6;

  //6. Track finder settings

  sf_f[0] = 0.2;//eta bin width

  //6b. clone removal
  
  sf_f[1] = 0.016;
  sf_f[2] = 28.0;
  sf_f[3] = 0.92;
  sf_f[4] = 0.2;
  
  sf_f[5] = 0.49;
  sf_f[6] = 0.85;
  
}
