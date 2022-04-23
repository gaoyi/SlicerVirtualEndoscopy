#ifndef SeededOptimalPathFilter3D_h
#define SeededOptimalPathFilter3D_h

#include <math.h>
#include <queue>
#include <set>
#include <vector>
#include <stdlib.h>
#include <fstream>
#include <iterator>

#include "itkImage.h"


#include "HeapNodeWithPre.h"

namespace gth818n
{
  template<typename SourceImageType, typename LabelImageType>
  class SeededOptimalPathFilter3D
  {
  public:
    SeededOptimalPathFilter3D();
    ~SeededOptimalPathFilter3D();

    void SetSourceImage(const SourceImageType* img);
    void SetSeedlImage(const LabelImageType* seedImage);

    typename LabelImageType::Pointer GetBackTraceLabelImage() {return m_backTraceLabelImage;}

    typedef itk::Image<float, 3> FloatImageType;
    FloatImageType::Pointer GetDistanceImage() {return m_distanceImage;}

    void update();

  private:
    float DIST_INF;
    float DIST_EPSION;
    unsigned char NNGBH;

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

    typename LabelImageType::Pointer m_backTraceLabelImage;
    typename LabelImageType::PixelType* m_backTraceLabelImageBufferPointer;

    typename LabelImageType::PixelType m_startSeedLabelValue;
    typename LabelImageType::PixelType m_endSeedLabelValue;

    std::vector<long> m_startSeedIndexList;
    std::vector<long> m_endSeedIndexList;

    FloatImageType::Pointer m_distanceImage;
    FloatImageType::PixelType* m_distanceImageBufferPointer;

    long m_DIMX, m_DIMY, m_DIMZ, m_DIMXY, m_DIMXYZ;
    std::vector<int> m_indOff;
    std::vector<double> m_spacingToNeighbors;
    std::vector<unsigned char>  m_NBSIZE;

    FibHeap *m_heap;
    HeapNodeWithPre *m_hpNodes;

    void SetupBeforeRun();

    void DoSegmentation();

    void BackTrace();

    void scanSeedLabelImage();
  };

} ///< namespace gth818n

#include "SeededOptimalPathFilter3D.hxx"

#endif
