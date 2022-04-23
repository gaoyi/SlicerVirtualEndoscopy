#ifndef LJHeapNode_h
#define LJHeapNode_h

#include "fibheap.h"

namespace gth818n
{
  /***************************************************************************
   * class LJHeapNode
   ***************************************************************************/
  class LJHeapNode : public FibHeapNode
  {
    float   N; // The key value in the heap. It's the "distance value" when computing shortest path.
    long IndexV; // The image linear index of the node with smallest key

  public:
    LJHeapNode() : FibHeapNode() { N = 0; }

    virtual void operator =(FibHeapNode& RHS)
    {
      FHN_Assign(RHS);
      N = ((LJHeapNode&) RHS).N;
    }

    virtual int  operator ==(FibHeapNode& RHS)
    {
      if (FHN_Cmp(RHS)) return 0;
      return N == ((LJHeapNode&) RHS).N ? 1 : 0;
    }

    virtual int  operator <(FibHeapNode& RHS)
    {
      int X;
      if ((X=FHN_Cmp(RHS)) != 0)
        return X < 0 ? 1 : 0;
      return N < ((LJHeapNode&) RHS).N ? 1 : 0;
    }

    virtual void operator =(double NewKeyVal )
    {
      LJHeapNode Tmp;
      Tmp.N = N = NewKeyVal;
      FHN_Assign(Tmp);
    }

    virtual void Print() { FibHeapNode::Print(); }

    double GetKeyValue() { return N; }
    void SetKeyValue(double n) { N = n; }
    long int GetIndexValue() { return IndexV; }
    void SetIndexValue( long int v) { IndexV = v; }
  };

} ///< namespace gth818n

#endif
