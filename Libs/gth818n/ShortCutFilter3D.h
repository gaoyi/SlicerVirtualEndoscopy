#ifndef ShortCutFilter3D_h
#define ShortCutFilter3D_h

#include <math.h>
#include <queue>
#include <set>
#include <vector>
#include <stdlib.h>
#include <fstream>
#include <iterator>

#include "LJHeapNode.h"

namespace gth818n
{
  //template<typename SrcPixelType, typename LabPixelType>
  template<typename SourceImageType, typename LabelImageType>
  class ShortCutFilter3D
  {
  public:
    ShortCutFilter3D();
    ~ShortCutFilter3D();

    void SetSourceImage(const SourceImageType* img);
    void SetSeedlImage(const LabelImageType* seedImage);
    void SetBackgroundLabel(typename LabelImageType::PixelType backgroundLabel) {m_backgroundLabel = backgroundLabel;}

    void SetWorkMode(bool bSegInitialized = false);

    void SetBackgroundDistance(float dist);


    void DoSegmentation();

    typename LabelImageType::Pointer GetLabeImage();

  private:
    float DIST_INF;
    float DIST_EPSION;
    unsigned char NNGBH;

    typedef float FPixelType;

    void InitializationAHP();
    void DijkstraBasedClassificationAHP();

    // // dbg
    // void checkLabelIntegrity();
    // // dbg, end

    typedef typename SourceImageType::PixelType SourcePixelType;
    typedef typename LabelImageType::PixelType LabelPixelType;

    const SourceImageType* m_img;
    const SourcePixelType* m_imgBufferPointer;

    const LabelImageType* m_seedImage;
    const LabelPixelType* m_seedImageBufferPointer;

    LabelPixelType m_backgroundLabel; ///< This label will be set to 0 in the output image

    std::vector<LabelPixelType> m_imLabPre;
    std::vector<FPixelType> m_imDistPre;
    std::vector<FPixelType> m_imDist;

    /// This is the default distance between two nodes regardless of
    /// the image intensity difference between two nodes. If this
    /// value is very high, then the process will basically ignore the
    /// image intensity and the result will be like the Voronoi
    /// partition of the space given the input seeds.
    float m_backgroundDistance;


    typename LabelImageType::Pointer m_labelImage;
    LabelPixelType* m_imLab;
    LabelPixelType* m_labelImageBufferPointer;

    //std::vector<long> m_imSize;
    long m_DIMX, m_DIMY, m_DIMZ, m_DIMXY, m_DIMXYZ;
    std::vector<int> m_indOff;
    std::vector<double> m_spacingToNeighbors;
    std::vector<unsigned char>  m_NBSIZE;

    FibHeap *m_heap;
    LJHeapNode *m_hpNodes;
    bool m_bSegInitialized;
  };

} ///< namespace gth818n

#include "ShortCutFilter3D.hxx"

#endif
