#ifndef HeapNodeWithPre_h
#define HeapNodeWithPre_h

#include "fibheap.h"

namespace gth818n
{
  /***************************************************************************
   * class HeapNodeWithPre
   ***************************************************************************/
  class HeapNodeWithPre : public FibHeapNode
  {
    float   N;
    long IndexV;
    long IndexOfPredecessor;

  public:
    HeapNodeWithPre() : FibHeapNode() { N = 0; }
    virtual void operator =(double NewKeyVal );
    virtual void operator =(FibHeapNode& RHS);

    virtual int  operator ==(FibHeapNode& RHS);

    virtual int  operator <(FibHeapNode& RHS);

    virtual void Print();

    double GetKeyValue() { return N; }
    void SetKeyValue(double n) { N = n; }

    long GetIndexValue() { return IndexV; }
    void SetIndexValue( long v) { IndexV = v; }

    long GetPredecessorIndexValue() { return IndexOfPredecessor; }
    void SetPredecessorIndexValue( long v) { IndexOfPredecessor = v; }
  };

} ///< namespace gth818n

#include "HeapNodeWithPre.hxx"

#endif
