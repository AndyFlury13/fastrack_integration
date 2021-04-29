#include<iostream>
#include<cmath>

#include "FieldModel.h"

FieldModel::FieldModel() {

  const double Cs = 0.00002999975;
  
  //solenoid parameters

  float ratio = 0.301;// R/L 0.3
  
  float Ls = 2*2218.0;//length 2x2218
  float Rs = ratio*Ls;//radius

  float Ls2 = Ls*Ls;
  float Rs2 = Rs*Rs;
  
  float C1 = Ls/sqrt(Ls2/4 + Rs2);
  float C2 = 0.5*C1/(Ls2/4 + Rs2);
  float C3 = 0.25*C2*Ls2/(Ls2/4 + Rs2);

  float dC = C2 - C3;
  
  float B0 = 21.95;//Bz in the center 22.0
  float alpha = B0/C1;//calibration const

  m_dC = dC*alpha;
  m_C1 = C1*alpha;

  m_dC_fast = m_dC*Cs;
  m_C1_fast = m_C1*Cs;
  
}

FieldModel::~FieldModel() {

}

void FieldModel::initialize(const ALGORITHM_PARAMETERS& p) {

  const double Cs = 0.00002999975;
  
  //solenoid parameters

  //for more detail see
  //arXiv:1003.3720v2 [physics.atom-ph]

  float ratio = p.mf_f[0];// solenoid aspect ratio: Radius/Half-Length
  
  float Ls = 2*p.mf_f[1];// solenoid length

  float B0 = p.mf_f[2];//Bz in the center

  float Rs = ratio*Ls;//radius

  float Ls2 = Ls*Ls;
  float Rs2 = Rs*Rs;
  
  float C1 = Ls/sqrt(Ls2/4 + Rs2);
  float C2 = 0.5*C1/(Ls2/4 + Rs2);
  float C3 = 0.25*C2*Ls2/(Ls2/4 + Rs2);

  float dC = C2 - C3;
  
  float alpha = B0/C1;//calibration const

  m_dC = dC*alpha;
  m_C1 = C1*alpha;

  m_dC_fast = m_dC*Cs;
  m_C1_fast = m_C1*Cs;

}


void FieldModel::getField(const float* P, float* B) const {

  float rt2 = P[0]*P[0] + P[1]*P[1];
  float rt_1 = 1.0/sqrt(rt2);

  float z   = P[2];
  float r2 = rt2 + z*z;
  float r   = sqrt(r2);
  
  float cost = z/r;    

  float cost2 = cost*cost;

  float sint = sqrt(1-cost2);

  float coeff = 1.5*r2*m_dC;
  float coeff2 = 2*rt_1*coeff*cost*sint;

  B[0] = coeff2*P[0];B[1] = coeff2*P[1];B[2] = m_C1 + coeff*(1-3*cost2);
}

void FieldModel::getFieldFast(const float* P, float* B) const {

  //returns field premultiplied by the speed of light to speed up the extrapolator code

  float rt2 = P[0]*P[0] + P[1]*P[1];
  float rt_1 = 1.0/sqrt(rt2);

  float z   = P[2];
  float r2 = rt2 + z*z;
  float r   = sqrt(r2);
  
  float cost = z/r;    

  float cost2 = cost*cost;

  float sint = sqrt(1-cost2);

  float coeff = 1.5*r2*m_dC_fast;
  float coeff2 = 2*rt_1*coeff*cost*sint;

  B[0] = coeff2*P[0];B[1] = coeff2*P[1];B[2] = m_C1_fast + coeff*(1-3*cost2);
}
