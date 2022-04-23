#ifndef SFLS_h_
#define SFLS_h_


// std
#include <list>

class CSFLS
{
public:
  typedef CSFLS Self;

  // The node is the index of the pixel, not coordinates, to avoid computing index from coordinates
  typedef long NodeType;
  typedef std::list< NodeType > CSFLSLayer; 

  // ctor
  CSFLS()  {}
  
  CSFLSLayer m_lz;
  CSFLSLayer m_ln1;
  CSFLSLayer m_ln2;
  CSFLSLayer m_lp1;
  CSFLSLayer m_lp2;
};

#endif
