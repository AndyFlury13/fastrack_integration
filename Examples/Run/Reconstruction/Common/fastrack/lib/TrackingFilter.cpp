#include<iostream>
#include<cmath>
#include<cstring>

#include "TrackingFilter.h"

#include<list>

TrackingFilter::TrackingFilter(const ALGORITHM_PARAMETERS& p, const HITS* h, NEW_SEG_BANK* sb[N_SEG_BANKS]) :
  m_params(p), m_hits(h), m_debug(false) {
  for(int i=0;i<N_SEG_BANKS;i++) m_segBank[i] = sb[i];
}


float TrackingFilter::followTrack(const SEGMENT* pS, std::vector<const SEGMENT*>& chain, std::vector<const NODE*>& track) {

  

  m_debug = false;
  
  if(pS->m_level == -1) return -999.0;//already collected

  m_globalCloneCounter = 0;

  //create track state

  int vtype = volumeType(pS->m_n2);

  TS_CLONE* pInitState = &m_tClones[m_globalCloneCounter++];
  
  pInitState->initialize(pS, vtype);
  
  m_tsVec.clear();

  //recursive branching and propagation

  //std::cout<<"at propagate ..."<<std::endl;
  
  propagate(pS, *pInitState);

  //std::cout<<"finished propagation, sorting now ..."<<std::endl;

  std::sort(m_tsVec.begin(), m_tsVec.end(), TS_CLONE::CompareJ());

  TS_CLONE* best = (*m_tsVec.begin());

  //if(m_debug) {
  //  std::cout<<"best path J="<<best->m_J<<" chi2="<<best->m_chi2y<<" tau="<<best->m_Y[1]<<" r="<<best->m_Y[0]<<std::endl;
  //}

  chain = std::move(best->m_vs);
  track = std::move(best->m_vn);

  m_globalCloneCounter = 0;
  
  return best->m_J;
}

void TrackingFilter::propagate(const SEGMENT* pS, TS_CLONE& ts) {

  double dropThreshold = m_params.ou_d[0];

  if(m_globalCloneCounter >= MAX_TS) return;
  
  TS_CLONE* p_new_ts = &m_tClones[m_globalCloneCounter++];
  
  TS_CLONE& new_ts = *p_new_ts;
  new_ts.clone(ts);

  new_ts.addSegment(pS);
  
  //1. update new_ts with pS->m_n1 == n2 of the possible neighbours

  update(pS->m_n1, new_ts);
  
  int level = pS->m_level;
  
  if(!m_tsVec.empty()) {
    float best_so_far = (*m_tsVec.begin())->m_J;
    float bestJ = level - 1 + ts.m_J;
    if(best_so_far > dropThreshold*bestJ) return;
  }

  std::list<const SEGMENT*> lCont;

  for(int nIdx=0;nIdx<pS->m_nNei;nIdx++) {
    unsigned int segIdx = pS->m_vNei[nIdx];
    unsigned int nextBankId = (segIdx & SEG_BANK_MASK) >> SEG_MASK_SHIFT;
    unsigned int nextSegmentIdx = (segIdx & SEG_INDEX_MASK);
    
    const SEGMENT* pN = &(m_segBank[nextBankId]->m_S[nextSegmentIdx]);

    if(pN->m_level == -1) continue;

    if(pN->m_level == level - 1) {

      float phi = ts.m_X[0];

      float phiN = pN->m_p[2];
      float dPhi = phi - phiN;
      if(dPhi < -M_PI) dPhi += 2*M_PI;
      if(dPhi >  M_PI) dPhi -= 2*M_PI;

      if(fabs(dPhi)>0.5) continue;
      
      lCont.push_back(pN);
    }
  }

  if(lCont.empty()) {//the end of chain

    //store in the vector
    if(m_globalCloneCounter < MAX_TS) {

      if(m_tsVec.empty()) {
	TS_CLONE* p = &m_tClones[m_globalCloneCounter++];
	p->clone(new_ts);
	m_tsVec.push_back(p);
      }
      else {
	float best_so_far = (*m_tsVec.begin())->m_J;
	if(new_ts.m_J > best_so_far) {
	  TS_CLONE* p = &m_tClones[m_globalCloneCounter++];
	  p->clone(new_ts);
	  m_tsVec.push_back(p);
	}
      }
    }
  } else {
    int nBranches = 0;
    for(std::list<const SEGMENT*>::iterator sIt = lCont.begin();sIt!=lCont.end();++sIt, nBranches++) {
      //      if(nBranches > 10) break;
      propagate((*sIt), new_ts);//recursive call
    }
  }
 
}

void TrackingFilter::update(const NODE* tn, TS_CLONE& ts) {

  const double sigmaMS      = m_params.ou_d[1];
  const double sigmaR_pix   = m_params.ou_d[2];
  const double sigmaR_shs   = m_params.ou_d[3];
  const double sigmaPhi_pix = m_params.ou_d[4];
  const double sigmaPhi_shs = m_params.ou_d[5];
  const double sigma_t      = m_params.ou_d[6];
  const double sigma_ms     = m_params.ou_d[7];
  const double sigma_w      = m_params.ou_d[8];
  const double alpha        = m_params.ou_d[9];
  const double chi2y_coeff  = m_params.ou_d[10];
  const double chi2x_coeff  = m_params.ou_d[11];
  const double maxDChi2_y   = m_params.ou_d[12];
  const double maxDChi2_x   = m_params.ou_d[13];

  const double sw2 = sigma_w*sigma_w;
  const double st2 = sigma_t*sigma_t;
  const double sigma_ms2 = sigma_ms*sigma_ms;

  //r-z update

  int type = ts.m_vtype;
  int ntype = volumeType(tn);

  //extrapolation
  
  //add noise

  double t2 = 1.0 + ts.m_Y[1]*ts.m_Y[1];
  double s1 = sigmaMS*t2;
  double s2 = s1*s1;

  s2 *= sqrt(t2);

  ts.m_Cy[1][1] += s2;

  ts.m_vtype = ntype;

  double Y[2];
  double Cy[2][2];
  
  double ds = 0.0;

  if(type !=0) {
    if(ntype != 0) {//D -> D

      double dz = tn->m_z - ts.m_refCoord;
      ds = dz*ts.m_Y[1];
      double dC = dz*ts.m_Cy[1][1];

      Y[0] = ts.m_Y[0] + ds;
      Y[1] = ts.m_Y[1];
      Cy[0][0] = ts.m_Cy[0][0] + dz*(2*ts.m_Cy[0][1] + dC);
      Cy[0][1] = Cy[1][0] = ts.m_Cy[0][1] + dC;
      Cy[1][1] = ts.m_Cy[1][1];

      ts.m_refCoord = tn->m_z;
    }
    else {          //D -> C
      
      //tan -> cot

      double r2 = tn->m_r;
      double ctg = 1.0/ts.m_Y[1];
      ds = r2 - ts.m_Y[0];
      double dCoord = ds*ctg;

      Y[0] = ts.m_refCoord + dCoord;
      Y[1] = ctg;

      double J[2][2];

      J[0][0] = -ctg; J[0][1] = ctg*dCoord;
      J[1][0] = 0.0;J[1][1] = -ctg*ctg;
      
      for(int i0=0;i0<2;i0++) {
	double JCy[2];
	JCy[0] = J[i0][0]*ts.m_Cy[0][0] + J[i0][1]*ts.m_Cy[1][0];
	JCy[1] = J[i0][0]*ts.m_Cy[0][1] + J[i0][1]*ts.m_Cy[1][1];
	
	for(int j2=0;j2<=i0;j2++) {
	  Cy[i0][j2] = 0.0;
	  for(int k=0;k<2;k++) Cy[i0][j2] += JCy[k]*J[j2][k];
	  Cy[j2][i0] = Cy[i0][j2];
	}
      }

      ts.m_refCoord = r2;
    }
  } 
  else {
    if(ntype == 0) {//C -> C
      
      //Y2 is cot, Y1 is z
      
      double r1 = ts.m_refCoord;
      double r2 = tn->m_r;
      double dr = r2-r1;

      ds = dr;
      double dC = dr*ts.m_Cy[1][1];

      Y[0] = ts.m_Y[0] + ts.m_Y[1]*dr;
      Y[1] = ts.m_Y[1];
      
      Cy[0][0] = ts.m_Cy[0][0] + dr*(2*ts.m_Cy[0][1] + dC);
      Cy[0][1] = Cy[1][0] = ts.m_Cy[0][1] + dC;
      Cy[1][1] = ts.m_Cy[1][1];
      
      ts.m_refCoord = r2;
    }
    else {//C->D
      
      //Y2 is tan, Y1 is r
      
      double z2 = tn->m_z;
      double r1 = ts.m_refCoord;
      double tg = 1.0/ts.m_Y[1];
      double dCoord = (z2 - ts.m_Y[0])*tg;

      Y[1] = tg;
      Y[0] = r1 + dCoord;
      
      ds = Y[0] - r1;
      
      double J[2][2];
      
      J[0][0] = -tg; J[0][1] = tg*dCoord;
      J[1][0] = 0.0;J[1][1] = -tg*tg;
      
      for(int i0=0;i0<2;i0++) {
	double JCy[2];
	JCy[0] = J[i0][0]*ts.m_Cy[0][0] + J[i0][1]*ts.m_Cy[1][0];
	JCy[1] = J[i0][0]*ts.m_Cy[0][1] + J[i0][1]*ts.m_Cy[1][1];
	for(int j2=0;j2<=i0;j2++) {
	  Cy[i0][j2] = 0.0;
	  for(int k=0;k<2;k++) Cy[i0][j2] += JCy[k]*J[j2][k];
	  Cy[j2][i0] = Cy[i0][j2];
	}
      }

      ts.m_refCoord = z2;
    }
  }

  double e1 = exp(-fabs(ds)*alpha);
  double f1 = (1.0-e1)/alpha;

  double F[3][3];

  memset(&F[0][0],0,sizeof(F));

  F[0][0] = 1.0; F[0][1] = ds;  F[0][2] = (fabs(ds)-f1)/alpha; 
                 F[1][1] = 1.0; F[1][2] = f1;
                                F[2][2] = e1;

  double X[3];
  double Cx[3][3];

  X[0] = ts.m_X[0] + F[0][1]*ts.m_X[1] + F[0][2]*ts.m_X[2];
  X[1] = ts.m_X[1] + F[1][2]*ts.m_X[2];
  X[2] = F[2][2]*ts.m_X[2];

  //double CFt[3][3];
  
  double CFt[3][3];

  CFt[0][0] = ts.m_Cx[0][0] + F[0][1]*ts.m_Cx[1][0] + F[0][2]*ts.m_Cx[2][0];
  CFt[0][1] = ts.m_Cx[0][1] + F[0][1]*ts.m_Cx[1][1] + F[0][2]*ts.m_Cx[2][1];
  CFt[0][2] = ts.m_Cx[0][2] + F[0][1]*ts.m_Cx[1][2] + F[0][2]*ts.m_Cx[2][2];
  
  CFt[1][0] = ts.m_Cx[1][0] + F[1][2]*ts.m_Cx[2][0];
  CFt[1][1] = ts.m_Cx[1][1] + F[1][2]*ts.m_Cx[2][1];
  CFt[1][2] = ts.m_Cx[1][2] + F[1][2]*ts.m_Cx[2][2];

  CFt[2][0] = F[2][2]*ts.m_Cx[2][0];
  CFt[2][1] = F[2][2]*ts.m_Cx[2][1];
  CFt[2][2] = F[2][2]*ts.m_Cx[2][2];

  
  Cx[0][0] = CFt[0][0] + CFt[0][1]*F[0][1] + CFt[0][2]*F[0][2];

  Cx[1][0] = CFt[1][0] + CFt[1][1]*F[0][1] + CFt[1][2]*F[0][2];
  Cx[1][1] = CFt[1][1] + CFt[1][2]*F[1][2];

  Cx[2][0] = CFt[2][0] + CFt[2][1]*F[0][1] + CFt[2][2]*F[0][2];
  Cx[2][1] = CFt[2][1] + CFt[2][2]*F[1][2];
  Cx[2][2] = CFt[2][2]*F[2][2];

  Cx[0][1] = Cx[1][0];Cx[0][2] = Cx[2][0];Cx[1][2] = Cx[2][1];

  double Q[3][3];
   
  double ds2 = ds*ds;
  double dsw2 = ds2*sw2;

  Q[0][0] = ds2*(st2 + 0.25*dsw2) + sigma_ms2;
  Q[0][1] = Q[1][0] = ds*(st2 + 0.5*dsw2);
  Q[0][2] = Q[2][0] = 0.5*dsw2;
  Q[1][1] = st2 + dsw2;
  Q[1][2] = Q[2][1] = ds*sw2;
  Q[2][2] = sw2;
  
  for(int i=0;i<3;i++) {
    for(int j=0;j<3;j++) {
      Cx[i][j] += Q[i][j];
    }
  }

  int vol_id = m_hits->m_vol_id[tn->m_h];
  
  double sigmaR   = (vol_id >= 10) ? sigmaR_shs   : sigmaR_pix;

  double sigmaPhi = (vol_id >= 10) ? sigmaPhi_shs : sigmaPhi_pix;
  

  double meas_y = (ntype == 0) ? tn->m_z : tn->m_r;
  
  double meas_x = tn->m_phi; 

  double detr_y = 1.0/(sigmaR*sigmaR + Cy[0][0]);
  
  double detr_x = 1.0/(sigmaPhi*sigmaPhi + Cx[0][0]);
  
  //std::cout<<"measurement "<<tn->m_r<<" expected "<<Y[0]<<std::endl;
  
  double resid_y = meas_y - Y[0];

  double dchi2y = resid_y*resid_y*detr_y;

  double resid_x = meas_x - X[0];

  if(resid_x < -M_PI) resid_x += 2*M_PI;
  if(resid_x >  M_PI) resid_x -= 2*M_PI;

  double dchi2x = resid_x*resid_x*detr_x;
  
  if(dchi2y < maxDChi2_y && dchi2x < maxDChi2_x) {// && !tn->isMasked(m_hits)) {
  
    ts.m_J += 1.0 - chi2y_coeff*dchi2y - chi2x_coeff*dchi2x;

    ts.addNode(tn);//add a node

    ts.m_chi2y += dchi2y;
    ts.m_chi2x += dchi2x;

    //std::cout<<"dchi2="<<dchi2<<std::endl;

    double Ky[2], Kx[3];

    Ky[0] = detr_y*Cy[0][0];
    Ky[1] = detr_y*Cy[0][1];
    Kx[0] = detr_x*Cx[0][0];
    Kx[1] = detr_x*Cx[0][1];
    Kx[2] = detr_x*Cx[0][2];


    ts.m_Y[0] = Y[0] + Ky[0]*resid_y;
    ts.m_Y[1] = Y[1] + Ky[1]*resid_y;

    ts.m_X[0] = X[0] + Kx[0]*resid_x;
    ts.m_X[1] = X[1] + Kx[1]*resid_x;
    ts.m_X[2] = X[2] + Kx[2]*resid_x;

    ts.m_Cy[0][0] = Cy[0][0] - Ky[0]*Cy[0][0];
    
    ts.m_Cy[0][1] = Cy[0][1] - Ky[0]*Cy[0][1];
    ts.m_Cy[1][0] = ts.m_Cy[0][1];
    
    ts.m_Cy[1][1] = Cy[1][1] - Ky[1]*Cy[0][1];
    
    ts.m_Cx[0][0] = Cx[0][0] - Kx[0]*Cx[0][0];
    ts.m_Cx[0][1] = Cx[0][1] - Kx[0]*Cx[0][1];
    ts.m_Cx[1][0] = ts.m_Cx[0][1];
 
    ts.m_Cx[0][2] = Cx[0][2] - Kx[0]*Cx[0][2];
    ts.m_Cx[2][0] = ts.m_Cx[0][2];

    ts.m_Cx[1][1] = Cx[1][1] - Kx[1]*Cx[0][1];
    ts.m_Cx[1][2] = Cx[1][2] - Kx[1]*Cx[0][2];
    ts.m_Cx[2][1] = ts.m_Cx[1][2];

    ts.m_Cx[2][2] = Cx[2][2] - Kx[2]*Cx[0][2];
  }
  else {//just propagate 

    ts.m_Y[0] = Y[0];
    ts.m_Y[1] = Y[1];

    memcpy(&ts.m_X[0], &X[0], sizeof(ts.m_X));
    memcpy(&ts.m_Cy[0][0], &Cy[0][0], sizeof(ts.m_Cy));
    memcpy(&ts.m_Cx[0][0], &Cx[0][0], sizeof(ts.m_Cx));

    //if(tn->isMasked()) ts.m_J -= 100.0; 
  }

  if(ts.m_X[0] < -M_PI) ts.m_X[0] += 2*M_PI;
  if(ts.m_X[0] >  M_PI) ts.m_X[0] -= 2*M_PI;
  
  //std::cout<<"Updated: "<<ts.m_Y[0]<<" "<<ts.m_Y[1]<<std::endl;
  //std::cout<<"Covariance: "<<ts.m_Cy[0][0]<<" "<<ts.m_Cy[0][1]<<" "<<ts.m_Cy[1][1]<<std::endl;
}

int TrackingFilter::volumeType(const NODE* n) {

  int type = 1;//disks

  int vol_id = m_hits->m_vol_id[n->m_h];

  if(vol_id == 8 || vol_id == 13 || vol_id == 17) type = 0;//barrels

  return type;
}
