#ifndef HeapNodeWithPre_hxx
#define HeapNodeWithPre_hxx

// local
#include "HeapNodeWithPre.h"

namespace gth818n
{
  void HeapNodeWithPre::Print()
  {
    FibHeapNode::Print();
    //    mexPrintf( "%f (%d)" , N , IndexV );
  }

  void HeapNodeWithPre::operator =(double NewKeyVal)
  {
    HeapNodeWithPre Tmp;
    Tmp.N = N = NewKeyVal;
    FHN_Assign(Tmp);
  }

  void HeapNodeWithPre::operator =(FibHeapNode& RHS) {
    FHN_Assign(RHS);
    N = ((HeapNodeWithPre&) RHS).N;
  }

  int HeapNodeWithPre::operator ==(FibHeapNode& RHS) {
    if (FHN_Cmp(RHS)) return 0;
    return N == ((HeapNodeWithPre&) RHS).N ? 1 : 0;
  }

  int HeapNodeWithPre::operator <(FibHeapNode& RHS) {
    int X;
    if ((X=FHN_Cmp(RHS)) != 0)
      return X < 0 ? 1 : 0;
    return N < ((HeapNodeWithPre&) RHS).N ? 1 : 0;
  }

} ///< namespace FGC


#endif
