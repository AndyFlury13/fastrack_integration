#include<iostream>
#include<fstream>
#include<algorithm>

#include "LayerLinker.h"

LayerLinker::LayerLinker(std::ifstream& linkFile) {
  
  unsigned int nSources;

  linkFile.read((char*)&nSources, sizeof(nSources));

  int linkIndex=0;
  
  for(unsigned int idx=0;idx<nSources;idx++) {

    unsigned int src;
    float totalFlow;
    unsigned int nLinks;
    
    linkFile.read((char*)&src, sizeof(src));
    linkFile.read((char*)&totalFlow, sizeof(totalFlow));
    linkFile.read((char*)&nLinks, sizeof(nLinks));
    
    for(unsigned int k=0;k<nLinks;k++) {

      unsigned int dst;
      float prob, flow;
      
      linkFile.read((char*)&dst, sizeof(dst));
      linkFile.read((char*)&prob, sizeof(prob));
      linkFile.read((char*)&flow, sizeof(flow));

      m_links.push_back(LAYER_LINK(linkIndex++,src, dst, prob, flow));
      //std::cout<<"Source "<<src<<" Dest "<<dst<<std::endl;
    }
  }

  linkFile.close();
  
  std::sort(m_links.begin(), m_links.end(), LAYER_LINK::CompareFlow());
  /*
  std::ofstream lf("layer_links.csv");
  
  lf<<"to,from,score_flow"<<std::endl;
  
  for(std::vector<LAYER_LINK>::iterator it=m_links.begin();it!=m_links.end();++it) {

    if((*it).m_src >= 16000 || (*it).m_dst >= 16000) continue;//skipping the long strips
    
    lf<<(*it).m_src<<","<<(*it).m_dst<<","<<(*it).m_flow<<std::endl;
  }

  lf.close();
  */
}
