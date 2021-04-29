#ifndef __PATTERN_CUTS_H__
#define __PATTERN_CUTS_H__

typedef struct patternCuts {
public:

  float zvMax, gammaMS;
  
  float maxZ0_barr, maxZ0_all;

  float dEta_min_pix, dEta_max_pix, dEta_min_shs, dEta_max_shs, dEta_min_mix, dEta_max_mix;
  float dPhi_min_pix, dPhi_max_pix, dPhi_min_shs, dPhi_max_shs, dPhi_min_mix, dPhi_max_mix;
  float dCurv_max;

} PATTERN_CUTS;


#endif
