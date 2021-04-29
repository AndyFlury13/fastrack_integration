#include "FastTrackExtrapolator.h"
#include<x86intrin.h>

void FastTrackExtrapolator::crossProduct(float const * B, float const * V, float* A) const {
  A[0] = -B[1]*V[2] + B[2]*V[1];
  A[1] =  B[0]*V[2] - B[2]*V[0];
  A[2] = -B[0]*V[1] + B[1]*V[0];
}

void FastTrackExtrapolator::fastCrossProduct(float const * B, float const * V, float* A) const {
  __m128 a = _mm_load_ps(B);
  __m128 b = _mm_load_ps(V);

  __m128 a_yzx = _mm_shuffle_ps(a, a, _MM_SHUFFLE(3, 0, 2, 1));

  __m128 b_yzx = _mm_shuffle_ps(b, b, _MM_SHUFFLE(3, 0, 2, 1));

  __m128 c = _mm_sub_ps(_mm_mul_ps(a, b_yzx), _mm_mul_ps(a_yzx, b));

  _mm_store_ps(A, _mm_shuffle_ps(c, c, _MM_SHUFFLE(3, 0, 2, 1)));

}

void FastTrackExtrapolator::fastCrossProduct3(float const * B,
				     float const * V1, float const * V2, float const * V3,
				     float* A1, float* A2, float* A3) const {
  __m128 a = _mm_load_ps(B);
  __m128 a_yzx = _mm_shuffle_ps(a, a, _MM_SHUFFLE(3, 0, 2, 1));

  __m128 b1 = _mm_load_ps(V1);
  __m128 b2 = _mm_load_ps(V2);
  __m128 b3 = _mm_load_ps(V3);

  __m128 b_yzx = _mm_shuffle_ps(b1, b1, _MM_SHUFFLE(3, 0, 2, 1));
  __m128 c = _mm_sub_ps(_mm_mul_ps(a, b_yzx), _mm_mul_ps(a_yzx, b1));
  _mm_store_ps(A1, _mm_shuffle_ps(c, c, _MM_SHUFFLE(3, 0, 2, 1)));

  b_yzx = _mm_shuffle_ps(b2, b2, _MM_SHUFFLE(3, 0, 2, 1));
  c = _mm_sub_ps(_mm_mul_ps(a, b_yzx), _mm_mul_ps(a_yzx, b2));
  _mm_store_ps(A2, _mm_shuffle_ps(c, c, _MM_SHUFFLE(3, 0, 2, 1)));

  b_yzx = _mm_shuffle_ps(b3, b3, _MM_SHUFFLE(3, 0, 2, 1));
  c = _mm_sub_ps(_mm_mul_ps(a, b_yzx), _mm_mul_ps(a_yzx, b3));
  _mm_store_ps(A3, _mm_shuffle_ps(c, c, _MM_SHUFFLE(3, 0, 2, 1)));

}

void  FastTrackExtrapolator::extrapolateOutsideFast(const TRACK_STATE* pTS, TRAJECTORY& vP) const {

  float Step    = m_params.te_f[4];
  int nMaxStep  = m_params.te_i[0];
  float stepAdj = m_params.te_f[5];
  float epsilon = m_params.te_f[6];
  float Rbound  = m_params.te_f[7];
  float Zbound  = m_params.te_f[8];

  float Re[5];
  memcpy(&Re[0], &pTS->m_Rk[0], sizeof(Re));

  float gV[3], gP[3], lP[3];

  float sinf, cosf, sint, cost;

  __sincosf(Re[2], &sinf, &cosf);
  __sincosf(Re[3], &sint, &cost);

  gV[0]=sint*cosf;gV[1]=sint*sinf;gV[2]=cost;

  lP[0] = pTS->m_Re[0];
  lP[1] = pTS->m_Re[1];
  lP[2] = 0.0;

  pTS->m_pDE->toGlobal(lP, gP);

  float B[3];

  m_fModel.getFieldFast(gP, B);

  vP.addPoint(gP[0], gP[1], gP[2]);

  float H = Step;

  double Cq = pTS->m_Re[4];

  float Y[6] = {gP[0], gP[1], gP[2], gV[0], gV[1], gV[2]};

  for(int iStep = 0;iStep<nMaxStep;iStep++) {

    float B2[3], B3[3];

    float H2 = 0.5*H;
    float H6 = H2/3;
    float CqH  = Cq*H;
    float CqH2 = Cq*H2;
    float CqH6 = Cq*H6;


    float YB1[3] = {Y[4]*B[2] - Y[5]*B[1], Y[5]*B[0] - Y[3]*B[2], Y[3]*B[1] - Y[4]*B[0]};

    //second point

    float Y2[6];

    for(int i=0;i<3;i++) Y2[i] = Y[i] + H2*Y[i+3];
    for(int i=3;i<6;i++) Y2[i] = Y[i] + CqH2*YB1[i-3];

    m_fModel.getFieldFast(Y2, B2);

    float YB2[3] = {Y2[4]*B2[2] - Y2[5]*B2[1], Y2[5]*B2[0] - Y2[3]*B2[2], Y2[3]*B2[1] - Y2[4]*B2[0]};

    //third point

    float Y3[6];

    for(int i=0;i<3;i++) Y3[i] = Y[i] + H2*Y2[i+3];
    for(int i=3;i<6;i++) Y3[i] = Y[i] + CqH2*YB2[i-3];

    m_fModel.getFieldFast(Y3, B3);

    float YB3[3] = {Y3[4]*B3[2] - Y3[5]*B3[1], Y3[5]*B3[0] - Y3[3]*B3[2], Y3[3]*B3[1] - Y3[4]*B3[0]};

    //last point

    float Y4[6];

    float B4[3];

    for(int i=0;i<3;i++) Y4[i] = Y[i] + H*Y3[i+3];
    for(int i=3;i<6;i++) Y4[i] = Y[i] + CqH*YB3[i-3];

    m_fModel.getFieldFast(Y4, B4);

    float YB4[3] = {Y4[4]*B4[2] - Y4[5]*B4[1], Y4[5]*B4[0] - Y4[3]*B4[2], Y4[3]*B4[1] - Y4[4]*B4[0]};


    float dd=0.0;
    for(int i=0;i<3;i++) {
      dd += fabs((YB1[i]+YB4[i])-(YB2[i]+YB3[i]));
    }
    dd = fabs(CqH*H)*dd;

    if(dd>epsilon && H > 25.0) {
      H *= stepAdj;

      continue;
    }

    H = Step;


    float Y1[6];

    for(int i=3;i<6;i++) Y1[i-3]  = Y[i-3] + H6*(Y[i] + 2*Y2[i] + 2*Y3[i] + Y4[i]);
    for(int i=0;i<3;i++) Y1[i+3]  = Y[i+3] + CqH6*(YB1[i] + 2*YB2[i] + 2*YB3[i] + YB4[i]);

    if(fabs(Y1[5])>1) break;

    if(fabs(Y1[2]) > Zbound) break;

    float R = sqrtf(Y1[0]*Y1[0]+Y1[1]*Y1[1]);
    if(R > Rbound) break;

    m_fModel.getFieldFast(Y1, B);

    vP.addPoint(Y1[0], Y1[1], Y1[2]);

    memcpy(&Y[0], &Y1[0], sizeof(Y));
  }

}

void FastTrackExtrapolator::extrapolateInsideFast(const TRACK_STATE* pTS, TRAJECTORY& vP) const {

  const float Step    = m_params.te_f[9];
  const int nMaxStep  = m_params.te_i[1];
  const float Rmin    = m_params.te_f[10];

  float Re[5];

  memcpy(&Re[0], &pTS->m_Rk[0], sizeof(Re));

  float gV[3], gP[3], lP[3];

  float sinf, cosf;
  float sint, cost;

  __sincosf(Re[2], &sinf, &cosf);
  __sincosf(Re[3], &sint, &cost);

  gV[0]=sint*cosf;gV[1]=sint*sinf;gV[2]=cost;

  lP[0] = pTS->m_Re[0];
  lP[1] = pTS->m_Re[1];
  lP[2] = 0.0;

  pTS->m_pDE->toGlobal(lP, gP);

  float R0 = sqrt(gP[0]*gP[0]+gP[1]*gP[1]);

  float B[3];

  m_fModel.getFieldFast(gP, B);

  vP.addPoint(gP[0], gP[1], gP[2]);

  float mom = pTS->m_Re[4];

  float Y[6] = {gP[0], gP[1], gP[2], gV[0], gV[1], gV[2]};

  float H = Step;
  float H23 = 2*H/3;
  float H4 = 0.25*H;
  float H23mom = mom*H23;
  float H4mom =  mom*H4;

  for(int iStep = 0;iStep<nMaxStep;iStep++) {

    float YB[3];

    crossProduct(B, Y+3, YB);

    //mid point

    float Y2[6];

    for(int i=0;i<3;i++) Y2[i] = Y[i] + H23*Y[i+3];
    for(int i=3;i<6;i++) Y2[i] = Y[i] + H23mom*YB[i-3];

    float B2[3];
    m_fModel.getFieldFast(Y2, B2);

    float YB2[3];

    crossProduct(B2, Y2+3, YB2);

    float Y1[6];

    for(int i=3;i<6;i++) Y1[i-3]  = Y[i-3] + H4*(Y[i] + 3*Y2[i]);
    for(int i=0;i<3;i++) Y1[i+3]  = Y[i+3] + H4mom*(YB[i] + 3*YB2[i]);

    if(fabs(Y1[5])>1) break;

    float R = sqrt(Y1[0]*Y1[0]+Y1[1]*Y1[1]);

    if(R < Rmin) break;
    if(R > R0) break;
    if(R < 50.0) {
      float dot_prod = Y1[3]*Y1[0]+Y1[4]*Y1[1];
      if(dot_prod<0.0) break;
    }

    m_fModel.getFieldFast(Y1, B);

    vP.addPoint(Y2[0], Y2[1], Y2[2]);
    vP.addPoint(Y1[0], Y1[1], Y1[2]);

    memcpy(&Y[0], &Y1[0], sizeof(Y));
  }

}

int FastTrackExtrapolator::fastExtrapolate(TRACK_STATE* pTS, int h, float minStep, TRAJECTORY& vP) const {

  const DETECTOR_ELEMENT* pNE = m_geo.getDetectorElement(m_hits->m_vol_id[h],
							 m_hits->m_lay_id[h],
							 m_hits->m_mod_id[h]);


  memcpy(&pTS->m_Re[0], &pTS->m_Rk[0], sizeof(pTS->m_Re));

  if(pTS->m_pDE == pNE) {//the same element
    memcpy(&pTS->m_Ce[0], &pTS->m_Ck[0], sizeof(pTS->m_Ce));
    return 0;
  }

  float Ce[15], Re[5];

  memcpy(&Ce[0], &pTS->m_Ck[0], sizeof(Ce));
  memcpy(&Re[0], &pTS->m_Rk[0], sizeof(Re));



  float gV[3];

  float sinf, cosf;
  float sint, cost;

  __sincosf(Re[2], &sinf, &cosf);
  __sincosf(Re[3], &sint, &cost);

  gV[0]=cosf*sint;gV[1]=sinf*sint;gV[2]=cost;

  //1. apply material effects

  float lV = pTS->m_pDE->localVz(gV);

  float rl = m_params.te_f[2]*pTS->m_pDE->m_module_t/fabs(lV);
  float pt2 = Re[4]*Re[4]*rl;

  //energy loss corrections

  Re[4]    += m_params.te_f[1]*Re[4]*rl*(1.0-0.5*rl);


  float ms2=0.00018496*pt2;
  float sinThetaLoc=sqrt(1-lV*lV);
  float invSin=1/sinThetaLoc;
  float mis = ms2*invSin;
  float eLoss = m_params.te_f[0]*pt2*(0.415-0.744*rl);

  Ce[5] += mis*invSin;
  Ce[9] += ms2;
  Ce[8] += mis;
  //energy loss corrections
  Ce[14] += eLoss;

  float P[8];//parameters + path
  float Jm[40];//jacobian

  memset(&P[0],0,sizeof(P));
  memset(&Jm[0],0,sizeof(Jm));

  //2.

  float lP[3], gP[3];

  lP[0] = pTS->m_Re[0];
  lP[1] = pTS->m_Re[1];
  lP[2] = 0.0;

  pTS->m_pDE->toGlobal(lP, gP);

  //parameters

  P[0] = gP[0];
  P[1] = gP[1];
  P[2] = gP[2];
  P[3] = gV[0];
  P[4] = gV[1];
  P[5] = gV[2];
  P[6] = Re[4];
  P[7] = 0.0;

  //Jacobian

  Jm[17] = -gV[1];//17
  Jm[18] = gV[0];//18
  Jm[24] = cosf*cost;//24
  Jm[25] = sinf*cost;//25
  Jm[26] =   -sint;//26
  for(int i=0;i<3;i++) {
    Jm[i]   = pTS->m_pDE->m_LG[i][0];
    Jm[7+i] = pTS->m_pDE->m_LG[i][1];
  }

  Jm[34] = 1.0;//34

  //3.

  int code = RungeKutta3(P, Jm, pNE, minStep);

  if(code!=0) return code;

  //4.

  float J[21];

  memset(&J[0],0,sizeof(J));

  pNE->transformTrackStateToLocal(pTS->m_Re, P, Jm, J);

  vP.addPoint(P[0], P[1], P[2]);

  //5. J*C*J^T

  float V1[5] = {Ce[ 0], Ce[ 1], Ce[ 3], Ce[ 6], Ce[10]};
  float V2[5] = {Ce[ 1], Ce[ 2], Ce[ 4], Ce[ 7], Ce[11]};
  float V3[5] = {Ce[ 3], Ce[ 4], Ce[ 5], Ce[ 8], Ce[12]};
  float V4[5] = {Ce[ 6], Ce[ 7], Ce[ 8], Ce[ 9], Ce[13]};
  float V5[5] = {Ce[10], Ce[11], Ce[12], Ce[13], Ce[14]};

  float A1[5], A2[5], A3[5], A4[5], A5[5];

  memset(&A1[0],0, sizeof(A1));
  memset(&A2[0],0, sizeof(A2));
  memcpy(&A5[0],&Ce[10],sizeof(A5));

  for(int i=0;i<5;i++) {

    float mult = J[i];

    A1[0] += mult*V1[i];
    A1[1] += mult*V2[i];
    A1[2] += mult*V3[i];
    A1[3] += mult*V4[i];
    A1[4] += mult*V5[i];

  }

  pTS->m_Ce[ 0] =  A1[0]*J[ 0] + A1[1]*J[ 1] + A1[2]*J[ 2] + A1[3]*J[ 3] + A1[4]*J[ 4];

  for(int i=0;i<5;i++) {

    float mult = J[i+5];

    A2[0] += mult*V1[i];
    A2[1] += mult*V2[i];
    A2[2] += mult*V3[i];
    A2[3] += mult*V4[i];
    A2[4] += mult*V5[i];
  }

  pTS->m_Ce[ 1] =  A2[0]*J[ 0] + A2[1]*J[ 1] + A2[2]*J[ 2] + A2[3]*J[ 3] + A2[4]*J[ 4];
  pTS->m_Ce[ 2] =  A2[0]*J[ 5] + A2[1]*J[ 6] + A2[2]*J[ 7] + A2[3]*J[ 8] + A2[4]*J[ 9];

  A3[0] = J[12]*V1[2] + J[13]*V1[3] + J[14]*V1[4];
  A3[1] = J[12]*V2[2] + J[13]*V2[3] + J[14]*V2[4];
  A3[2] = J[12]*V3[2] + J[13]*V3[3] + J[14]*V3[4];
  A3[3] = J[12]*V4[2] + J[13]*V4[3] + J[14]*V4[4];
  A3[4] = J[12]*V5[2] + J[13]*V5[3] + J[14]*V5[4];

  pTS->m_Ce[ 3] =  A3[0]*J[ 0] + A3[1]*J[ 1] + A3[2]*J[ 2] + A3[3]*J[ 3] + A3[4]*J[ 4];
  pTS->m_Ce[ 4] =  A3[0]*J[ 5] + A3[1]*J[ 6] + A3[2]*J[ 7] + A3[3]*J[ 8] + A3[4]*J[ 9];
  pTS->m_Ce[ 5] =                              A3[2]*J[12] + A3[3]*J[13] + A3[4]*J[14];

  A4[0] = J[17]*V1[2] + J[18]*V1[3] + J[19]*V1[4];
  A4[1] = J[17]*V2[2] + J[18]*V2[3] + J[19]*V2[4];
  A4[2] = J[17]*V3[2] + J[18]*V3[3] + J[19]*V3[4];
  A4[3] = J[17]*V4[2] + J[18]*V4[3] + J[19]*V4[4];
  A4[4] = J[17]*V5[2] + J[18]*V5[3] + J[19]*V5[4];

  pTS->m_Ce[ 6] =  A4[0]*J[ 0] + A4[1]*J[ 1] + A4[2]*J[ 2] + A4[3]*J[ 3] + A4[4]*J[ 4];
  pTS->m_Ce[ 7] =  A4[0]*J[ 5] + A4[1]*J[ 6] + A4[2]*J[ 7] + A4[3]*J[ 8] + A4[4]*J[ 9];
  pTS->m_Ce[ 8] =                              A4[2]*J[12] + A4[3]*J[13] + A4[4]*J[14];
  pTS->m_Ce[ 9] =                              A4[2]*J[17] + A4[3]*J[18] + A4[4]*J[19];



  pTS->m_Ce[10] =  A5[0]*J[ 0] + A5[1]*J[ 1] + A5[2]*J[ 2] + A5[3]*J[ 3] + A5[4]*J[ 4];
  pTS->m_Ce[11] =  A5[0]*J[ 5] + A5[1]*J[ 6] + A5[2]*J[ 7] + A5[3]*J[ 8] + A5[4]*J[ 9];
  pTS->m_Ce[12] =                              A5[2]*J[12] + A5[3]*J[13] + A5[4]*J[14];
  pTS->m_Ce[13] =                              A5[2]*J[17] + A5[3]*J[18] + A5[4]*J[19];

  pTS->m_Ce[14] = Ce[14];

  pTS->m_pDE = pNE;

  return 0;
}

int FastTrackExtrapolator::fastExtrapolator(TRACK_STATE* pTS, int h, bool checkBounds) const {

  const DETECTOR_ELEMENT* pNE = m_geo.getDetectorElement(m_hits->m_vol_id[h],
							 m_hits->m_lay_id[h],
							 m_hits->m_mod_id[h]);


  memcpy(&pTS->m_Re[0], &pTS->m_Rk[0], sizeof(pTS->m_Re));

  if(pTS->m_pDE == pNE) {//the same element
    memcpy(&pTS->m_Ce[0], &pTS->m_Ck[0], sizeof(pTS->m_Ce));
    return 0;
  }

  float Ce[15], Re[5];

  memcpy(&Ce[0], &pTS->m_Ck[0], sizeof(Ce));
  memcpy(&Re[0], &pTS->m_Rk[0], sizeof(Re));

  //1. apply material effects

  float gV[3];

  float sinf, cosf;
  float sint, cost;

  __sincosf(Re[2], &sinf, &cosf);
  __sincosf(Re[3], &sint, &cost);

  gV[0]=cosf*sint;gV[1]=sinf*sint;gV[2]=cost;

  float lV = pTS->m_pDE->localVz(gV);

  float rl = m_params.te_f[2]*pTS->m_pDE->m_module_t/fabs(lV);
  float pt2 = Re[4]*Re[4]*rl;

  //energy loss corrections

  Re[4]    += m_params.te_f[1]*Re[4]*rl*(1.0-0.5*rl);


  float ms2=0.00018496*pt2;
  float sinThetaLoc=sqrt(1-lV*lV);
  float invSin=1/sinThetaLoc;
  float mis = ms2*invSin;
  float eLoss = m_params.te_f[0]*pt2*(0.415-0.744*rl);

  Ce[5] += mis*invSin;
  Ce[9] += ms2;
  Ce[8] += mis;
  //energy loss corrections
  Ce[14] += eLoss;

  float P[8];//parameters + path
  float Jm[40];//jacobian

  memset(&P[0],0,sizeof(P));
  memset(&Jm[0],0,sizeof(Jm));

  //2.

  float lP[3], gP[3];

  lP[0] = pTS->m_Re[0];
  lP[1] = pTS->m_Re[1];
  lP[2] = 0.0;

  pTS->m_pDE->toGlobal(lP, gP);

  P[0] = gP[0];
  P[1] = gP[1];
  P[2] = gP[2];
  P[3] = gV[0];
  P[4] = gV[1];
  P[5] = gV[2];
  P[6] = Re[4];
  P[7] = 0.0;

  Jm[17] = -gV[1];//17
  Jm[18] = gV[0];//18
  Jm[24] = cosf*cost;//24
  Jm[25] = sinf*cost;//25
  Jm[26] =   -sint;//26
  for(int i=0;i<3;i++) {
    Jm[i]   = pTS->m_pDE->m_LG[i][0];
    Jm[7+i] = pTS->m_pDE->m_LG[i][1];
  }

  Jm[34] = 1.0;//34

  //3.

  int code = RungeKutta3(P, Jm, pNE);

  if(code!=0) return code;

  //4.

  float J[21];

  memset(&J[0],0,sizeof(J));

  pNE->transformTrackStateToLocal(pTS->m_Re, P, Jm, J);

  if(checkBounds) {
    float tol = -1.0;
    float fP[3] = {pTS->m_Re[0], pTS->m_Re[1], 0.0};
    if(!pNE->isInSide(fP, tol)) return -20;
  }

  //5. J*C*J^T

  float V1[5] = {Ce[ 0], Ce[ 1], Ce[ 3], Ce[ 6], Ce[10]};
  float V2[5] = {Ce[ 1], Ce[ 2], Ce[ 4], Ce[ 7], Ce[11]};
  float V3[5] = {Ce[ 3], Ce[ 4], Ce[ 5], Ce[ 8], Ce[12]};
  float V4[5] = {Ce[ 6], Ce[ 7], Ce[ 8], Ce[ 9], Ce[13]};
  float V5[5] = {Ce[10], Ce[11], Ce[12], Ce[13], Ce[14]};

  float A1[5], A2[5], A3[5], A4[5], A5[5];

  memset(&A1[0],0, sizeof(A1));
  memset(&A2[0],0, sizeof(A2));
  memcpy(&A5[0],&Ce[10],sizeof(A5));

  for(int i=0;i<5;i++) {

    float mult = J[i];

    A1[0] += mult*V1[i];
    A1[1] += mult*V2[i];
    A1[2] += mult*V3[i];
    A1[3] += mult*V4[i];
    A1[4] += mult*V5[i];

  }

  pTS->m_Ce[ 0] =  A1[0]*J[ 0] + A1[1]*J[ 1] + A1[2]*J[ 2] + A1[3]*J[ 3] + A1[4]*J[ 4];

  for(int i=0;i<5;i++) {

    float mult = J[i+5];

    A2[0] += mult*V1[i];
    A2[1] += mult*V2[i];
    A2[2] += mult*V3[i];
    A2[3] += mult*V4[i];
    A2[4] += mult*V5[i];
  }

  pTS->m_Ce[ 1] =  A2[0]*J[ 0] + A2[1]*J[ 1] + A2[2]*J[ 2] + A2[3]*J[ 3] + A2[4]*J[ 4];
  pTS->m_Ce[ 2] =  A2[0]*J[ 5] + A2[1]*J[ 6] + A2[2]*J[ 7] + A2[3]*J[ 8] + A2[4]*J[ 9];

  A3[0] = J[12]*V1[2] + J[13]*V1[3] + J[14]*V1[4];
  A3[1] = J[12]*V2[2] + J[13]*V2[3] + J[14]*V2[4];
  A3[2] = J[12]*V3[2] + J[13]*V3[3] + J[14]*V3[4];
  A3[3] = J[12]*V4[2] + J[13]*V4[3] + J[14]*V4[4];
  A3[4] = J[12]*V5[2] + J[13]*V5[3] + J[14]*V5[4];

  A4[0] = J[17]*V1[2] + J[18]*V1[3] + J[19]*V1[4];
  A4[1] = J[17]*V2[2] + J[18]*V2[3] + J[19]*V2[4];
  A4[2] = J[17]*V3[2] + J[18]*V3[3] + J[19]*V3[4];
  A4[3] = J[17]*V4[2] + J[18]*V4[3] + J[19]*V4[4];
  A4[4] = J[17]*V5[2] + J[18]*V5[3] + J[19]*V5[4];

  pTS->m_Ce[ 3] =  A3[0]*J[ 0] + A3[1]*J[ 1] + A3[2]*J[ 2] + A3[3]*J[ 3] + A3[4]*J[ 4];
  pTS->m_Ce[ 4] =  A3[0]*J[ 5] + A3[1]*J[ 6] + A3[2]*J[ 7] + A3[3]*J[ 8] + A3[4]*J[ 9];
  pTS->m_Ce[ 5] =                              A3[2]*J[12] + A3[3]*J[13] + A3[4]*J[14];

  pTS->m_Ce[ 6] =  A4[0]*J[ 0] + A4[1]*J[ 1] + A4[2]*J[ 2] + A4[3]*J[ 3] + A4[4]*J[ 4];
  pTS->m_Ce[ 7] =  A4[0]*J[ 5] + A4[1]*J[ 6] + A4[2]*J[ 7] + A4[3]*J[ 8] + A4[4]*J[ 9];
  pTS->m_Ce[ 8] =                              A4[2]*J[12] + A4[3]*J[13] + A4[4]*J[14];
  pTS->m_Ce[ 9] =                              A4[2]*J[17] + A4[3]*J[18] + A4[4]*J[19];

  pTS->m_Ce[10] =  A5[0]*J[ 0] + A5[1]*J[ 1] + A5[2]*J[ 2] + A5[3]*J[ 3] + A5[4]*J[ 4];
  pTS->m_Ce[11] =  A5[0]*J[ 5] + A5[1]*J[ 6] + A5[2]*J[ 7] + A5[3]*J[ 8] + A5[4]*J[ 9];
  pTS->m_Ce[12] =                              A5[2]*J[12] + A5[3]*J[13] + A5[4]*J[14];
  pTS->m_Ce[13] =                              A5[2]*J[17] + A5[3]*J[18] + A5[4]*J[19];

  pTS->m_Ce[14] = Ce[14];

  pTS->m_pDE = pNE;

  return 0;
}


int FastTrackExtrapolator::RungeKutta3(float* P, float* J, const DETECTOR_ELEMENT* pNE, float exStep) const {

  float min_step         = m_params.te_f[11];
  float const_field_step = m_params.te_f[12];
  float maxPath          = m_params.te_f[13];
  float minQp            = m_params.te_f[14];
  float minRad           = m_params.te_f[15];

  int maxIter = 10;

  if(fabs(P[6]) > minQp) return -1;

  float planeParams[4];//target surface parameters

  pNE->getPlaneParameters(planeParams);//normal and (n,r0)

  float a = 0.0;
  float Sum = planeParams[3];

  for(int i=0;i<3;i++) {
    a   += planeParams[i]*P[i+3];
    Sum += planeParams[i]*P[i];
  }

  if(a==0.0) return -2;

  float Step = -Sum/a;

  float absStep   = fabs(Step);

  if(absStep <= min_step) {
    for(int i=0;i<3;i++) P[i] += Step*P[i+3];
    P[7] += Step;
    return 0;
  }

  if(fabs(P[6]*Step) > minRad) {
    Step = (Step > 0.0 ? minRad : -minRad)/fabs(P[6]);
  }

  int nFlips = 0;

  float mom = P[6];

  float Y[6];

  memcpy(&Y[0], P, sizeof(Y));

  float B[3];
  m_fModel.getFieldFast(Y, B);

  if(exStep != 0) {
    if(absStep > fabs(exStep))
      Step = exStep;
  }

  for(int iter=0;iter<maxIter;iter++) {

    bool useConstField = fabs(Step) < const_field_step;

    if(!useConstField) {
      m_fModel.getFieldFast(Y, B);
    }

    float B2[3], B3[3];

    float H = Step;
    float H3 = H/3;
    float H23 = 2*H3;
    float H4 = 0.25*H;
    float H34 = 3*H4;

    float H3mom =  mom*H3;
    float H23mom = mom*H23;
    float H4mom =  mom*H4;
    float H34mom = mom*H34;

    float YB[3];

    crossProduct(B, Y+3, YB);

    //second point

    float Y2[6];

    for(int i=0;i<3;i++) Y2[i] = Y[i] + H3*Y[i+3];
    for(int i=3;i<6;i++) Y2[i] = Y[i] + H3mom*YB[i-3];

    float YB2[3];

    if(useConstField) {
      crossProduct(B, Y2+3, YB2);
    }
    else {
      m_fModel.getFieldFast(Y2, B2);
      crossProduct(B2, Y2+3, YB2);
    }

    //last point

    float Y3[6];

    for(int i=0;i<3;i++) Y3[i] = Y[i] + H23*Y2[i+3];
    for(int i=3;i<6;i++) Y3[i] = Y[i] + H23mom*YB2[i-3];

    float YB3[3];

    if(useConstField) {
      crossProduct(B, Y3+3, YB3);
    }
    else {
      m_fModel.getFieldFast(Y3, B3);
      crossProduct(B3, Y3+3, YB3);
    }

    float Y1[6];

    for(int i=3;i<6;i++) Y1[i-3]  = Y[i-3] + H4*(Y[i] + 3*Y3[i]);
    for(int i=0;i<3;i++) Y1[i+3]  = Y[i+3] + H4mom*(YB[i] + 3*YB3[i]);

    if(fabs(Y1[5])>1) return -10;

    float J1C[9], L2C[9], J2C[9], L3C[9], J3C[9];

    float CqB3H34[3];
    float CqB2H23[3];
    float CqBH3[3];

    if(!useConstField) {
      for(int i=0;i<3;i++) CqBH3[i]   = H3mom*B[i];
      for(int i=0;i<3;i++) CqB2H23[i] = H23mom*B2[i];
      for(int i=0;i<3;i++) CqB3H34[i] = H34mom*B3[i];
    }
    else {
      for(int i=0;i<3;i++) CqBH3[i]   = H3mom*B[i];
      for(int i=0;i<3;i++) CqB2H23[i] = H23mom*B[i];
      for(int i=0;i<3;i++) CqB3H34[i] = H34mom*B[i];
    }

    crossProduct(CqBH3, J+17, J1C);
    crossProduct(CqBH3, J+24, J1C+3);
    crossProduct(CqBH3, J+31, J1C+6);

    J1C[6] += H3*YB[0];
    J1C[7] += H3*YB[1];
    J1C[8] += H3*YB[2];

    L2C[0] = J[17] + J1C[0];
    L2C[1] = J[18] + J1C[1];
    L2C[2] = J[19] + J1C[2];

    L2C[3] = J[24] + J1C[3];
    L2C[4] = J[25] + J1C[4];
    L2C[5] = J[26] + J1C[5];

    L2C[6] = J[31] + J1C[6];
    L2C[7] = J[32] + J1C[7];
    L2C[8] = J[33] + J1C[8];

    crossProduct(CqB2H23, L2C,   J2C);
    crossProduct(CqB2H23, L2C+3, J2C+3);
    crossProduct(CqB2H23, L2C+6, J2C+6);

    J2C[6] += H23*YB2[0];
    J2C[7] += H23*YB2[1];
    J2C[8] += H23*YB2[2];

    L3C[0] = J[17] + J2C[0];
    L3C[1] = J[18] + J2C[1];
    L3C[2] = J[19] + J2C[2];

    L3C[3] = J[24] + J2C[3];
    L3C[4] = J[25] + J2C[4];
    L3C[5] = J[26] + J2C[5];

    L3C[6] = J[31] + J2C[6];
    L3C[7] = J[32] + J2C[7];
    L3C[8] = J[33] + J2C[8];

    crossProduct(CqB3H34, L3C,   J3C);
    crossProduct(CqB3H34, L3C+3, J3C+3);
    crossProduct(CqB3H34, L3C+6, J3C+6);

    J3C[6] += H34*YB3[0];
    J3C[7] += H34*YB3[1];
    J3C[8] += H34*YB3[2];

    for(int i=0;i<9;i++) J1C[i] = 0.75*J1C[i] + J3C[i];
    for(int i=0;i<9;i++) J2C[i] *= H34;

    J[14] += H*J[17];
    J[15] += H*J[18];
    J[16] += H*J[19];

    J[21] += H*J[24];
    J[22] += H*J[25];
    J[23] += H*J[26];

    J[28] += H*J[31];
    J[29] += H*J[32];
    J[30] += H*J[33];

    J[14] += J2C[0];
    J[15] += J2C[1];
    J[16] += J2C[2];
    J[21] += J2C[3];
    J[22] += J2C[4];
    J[23] += J2C[5];
    J[28] += J2C[6];
    J[29] += J2C[7];
    J[30] += J2C[8];

    J[17] += J1C[0];
    J[18] += J1C[1];
    J[19] += J1C[2];
    J[24] += J1C[3];
    J[25] += J1C[4];
    J[26] += J1C[5];
    J[31] += J1C[6];
    J[32] += J1C[7];
    J[33] += J1C[8];

    P[7] += Step;

    if(fabs(P[7]) > maxPath) return -3;

    float norm = 1/sqrtf(Y1[3]*Y1[3]+Y1[4]*Y1[4]+Y1[5]*Y1[5]);
    Y1[3]*=norm; Y1[4]*=norm; Y1[5]*=norm;

    float a = 0.0;
    float Sum = planeParams[3];

    for(int i=0;i<3;i++) {
      a   += planeParams[i]*Y1[i+3];
      Sum += planeParams[i]*Y1[i];
    }

    if(a==0.0) return -4;

    float newStep = -Sum/a;

    float absNewStep = fabs(newStep);

    if(absNewStep <= min_step) {//straight line
      if(!useConstField) {
	crossProduct(B3, Y1+3, J+35);
      }
      else {
	crossProduct(B, Y1+3, J+35);
      }
      for(int i=0;i<3;i++) {
	P[i+3] = Y1[i+3];
	P[i] = Y1[i] + newStep*Y1[i+3];
      }
      P[7] += newStep;

      return 0;
    }

    float absStep = fabs(Step);

    if(Step*newStep < 0.0) {//overshot
      if(++nFlips > 2) return -5;//oscillations
      Step =  absNewStep < absStep ? newStep : -Step;
    }
    else
      if(absNewStep < absStep) Step = newStep;

    for(int i=0;i<6;i++) Y[i] = Y1[i];
  }

  return -1;
}
