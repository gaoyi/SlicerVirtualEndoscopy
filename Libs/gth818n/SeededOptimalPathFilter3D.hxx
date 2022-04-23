#ifndef SeededOptimalPathFilter3D_hxx
#define SeededOptimalPathFilter3D_hxx

#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkNumericTraits.h"
#include "itkRegionOfInterestImageFilter.h"

// local
#include "SeededOptimalPathFilter3D.h"


namespace gth818n
{
  template<typename SourceImageType, typename LabelImageType>
  SeededOptimalPathFilter3D<SourceImageType, LabelImageType>::SeededOptimalPathFilter3D()
  {
    m_heap = NULL;
    m_hpNodes = NULL;

    DIST_INF = std::numeric_limits<float>::max();
    DIST_EPSION = 1e-3;
    NNGBH = 26;

    m_startSeedLabelValue = 1;
    m_endSeedLabelValue = 2;

    return;
  }

  template<typename SourceImageType, typename LabelImageType>
  SeededOptimalPathFilter3D<SourceImageType, LabelImageType>::~SeededOptimalPathFilter3D()
  {
    if(m_heap != NULL)
      {
        delete m_heap;
        m_heap = NULL;
      }

    if(m_hpNodes != NULL)
      {
        delete []m_hpNodes;
        m_hpNodes = NULL;
      }

    //std::cout << "SeededOptimalPathFilter3D destroyed\n";

    return;
  }

  template<typename SourceImageType, typename LabelImageType>
  void
  SeededOptimalPathFilter3D<SourceImageType, LabelImageType>::SetSourceImage(const SourceImageType* img)
  {
    m_img = img;

    return;
  }

  template<typename SourceImageType, typename LabelImageType>
  void
  SeededOptimalPathFilter3D<SourceImageType, LabelImageType>::SetSeedlImage(const LabelImageType* seedImage)
  {
    m_seedImage = seedImage;

    return;
  }

  template<typename SourceImageType, typename LabelImageType>
  void
  SeededOptimalPathFilter3D<SourceImageType, LabelImageType>::InitializationAHP()
  {
    if((m_heap = new FibHeap) == NULL || (m_hpNodes = new HeapNodeWithPre[m_DIMXYZ+1]) == NULL)
      {
        std::cerr << "Memory allocation failed-- ABORTING.\n";
        abort();
      }
    m_heap->ClearHeapOwnership();


    long i,j,k, index;

    m_distanceImage = FloatImageType::New();
    m_distanceImage->SetRegions( m_img->GetLargestPossibleRegion() );
    m_distanceImage->Allocate();
    m_distanceImage->CopyInformation(m_img);
    m_distanceImage->FillBuffer(DIST_INF);

    m_distanceImageBufferPointer = m_distanceImage->GetBufferPointer();

    // Compute index offset
    m_indOff.clear();
    m_spacingToNeighbors.clear();

    long ix,iy,iz;

    typename SourceImageType::SpacingType spacing = m_img->GetSpacing();

    for(ix = -1; ix <= 1; ix++)
      {
        for(iy = -1; iy <= 1; iy++)
          {
            for(iz = -1; iz <= 1; iz++)
              {
                if(!(ix == 0 && iy == 0 && iz == 0))
                  {
                    m_indOff.push_back(ix + iy*m_DIMX + iz*m_DIMXY);
                    m_spacingToNeighbors.push_back(1.0/sqrt((ix*spacing[0])*(ix*spacing[0]) + (iy*spacing[1])*(iy*spacing[1]) + (iz*spacing[2])*(iz*spacing[2]))); ///< This is in fact inverse of spacing, so that in later computing the gradient, we use multiplication, not division
                  }
              }
          }
      }

    // Determine neighborhood size at each *non-boundary* vertice
    m_NBSIZE = std::vector<unsigned char>(m_DIMXYZ, 0);
    for(i = 1; i < m_DIMX - 1; i++)
      {
        for(j = 1; j < m_DIMY - 1; j++)
          {
            for(k = 1; k < m_DIMZ - 1; k++)
              {
                index = i + j*m_DIMX + k*m_DIMXY;
                m_NBSIZE[index] = NNGBH;
              }
          }
      }

    for(index = 0; index < m_DIMXYZ; index++)
      {
        if(m_seedImageBufferPointer[index] != m_startSeedLabelValue)
          {
            m_hpNodes[index] = (float)DIST_INF;
            //m_distanceImageBufferPointer[index] = DIST_INF;
          }
        else
          {
            m_hpNodes[index] = (float)DIST_EPSION;
            m_hpNodes[index].SetPredecessorIndexValue(-1);
            m_distanceImageBufferPointer[index] = DIST_EPSION;
          }

        m_heap->Insert(&m_hpNodes[index]);
        m_hpNodes[index].SetIndexValue(index);
      }


    return;
  }

  template<typename SourceImageType, typename LabelImageType>
  void
  SeededOptimalPathFilter3D<SourceImageType, LabelImageType>::DijkstraBasedClassificationAHP()
  {
    HeapNodeWithPre *hnMin;
    HeapNodeWithPre hnTmp;
    long i, index;

    LabelPixelType labSrc;
    SourcePixelType pixCenter;

    // Insert 0 then extract it, which will balance heap
    m_heap->Insert(&hnTmp);
    m_heap->ExtractMin();

    long numberOfEndSeedReached = 0;

    // Normal Dijkstra to be used in initializing the segmenter for the current image
    while(!m_heap->IsEmpty())
      {
        hnMin = (HeapNodeWithPre *) m_heap->ExtractMin();
        index = hnMin->GetIndexValue();
        float tSrc = hnMin->GetKeyValue();
        m_distanceImageBufferPointer[index] = tSrc;

        if (m_seedImageBufferPointer[index] == m_endSeedLabelValue)
          {
            ++numberOfEndSeedReached;

            std::cout<<numberOfEndSeedReached<<" end seed reached. total of "<<m_endSeedIndexList.size()<<"\n"<<std::flush;

            if (numberOfEndSeedReached == static_cast<long>(m_endSeedIndexList.size()))
              {
                std::cout<<"all end seeds reached, exit\n"<<std::flush;
                break;
              }
          }

        // Update neighbors
        pixCenter = m_imgBufferPointer[index];
        for(i = 0; i < m_NBSIZE[index]; i++)
          {
            long indexNgbh = index + m_indOff[i];
            float tOri = m_distanceImageBufferPointer[indexNgbh];

            float t = (pixCenter + m_imgBufferPointer[indexNgbh])/2.0 + tSrc;

            if(tOri > t)
              {
                m_distanceImageBufferPointer[indexNgbh] = t;

                hnTmp = m_hpNodes[indexNgbh];
                hnTmp.SetKeyValue(t);
                m_hpNodes[indexNgbh].SetPredecessorIndexValue(index);
                m_heap->DecreaseKey(&m_hpNodes[indexNgbh], hnTmp);
              }
          }
      }
  }

  template<typename SourceImageType, typename LabelImageType>
  void
  SeededOptimalPathFilter3D<SourceImageType, LabelImageType>::SetupBeforeRun()
  {
    if (!m_seedImage || !m_img)
      {
        std::cerr<<"!m_seedImage || !m_img\n";
        abort();
      }

    typedef typename SourceImageType::RegionType RegionType;
    RegionType region = m_img->GetLargestPossibleRegion();

    m_imgBufferPointer = m_img->GetBufferPointer();
    m_seedImageBufferPointer = m_seedImage->GetBufferPointer();

    m_DIMX = region.GetSize()[0];
    m_DIMY = region.GetSize()[1];
    m_DIMZ = region.GetSize()[2];
    m_DIMXY = m_DIMX*m_DIMY;
    m_DIMXYZ = m_DIMXY*m_DIMZ;


    scanSeedLabelImage();

    return;
  }


  template<typename SourceImageType, typename LabelImageType>
  void
  SeededOptimalPathFilter3D<SourceImageType, LabelImageType>::DoSegmentation()
  {
    InitializationAHP();

    DijkstraBasedClassificationAHP();
  }

  template<typename SourceImageType, typename LabelImageType>
  void
  SeededOptimalPathFilter3D<SourceImageType, LabelImageType>::update()
  {
    SetupBeforeRun();

    DoSegmentation();

    BackTrace();

    return;
  }

  template<typename SourceImageType, typename LabelImageType>
  void
  SeededOptimalPathFilter3D<SourceImageType, LabelImageType>::scanSeedLabelImage()
  {
    for(long index = 0; index < m_DIMXYZ; index++)
      {
        LabelPixelType l = m_seedImageBufferPointer[index];

        if(l == m_endSeedLabelValue)
          {
            m_endSeedIndexList.push_back(index);
          }
        else if (l == m_startSeedLabelValue)
          {
            m_startSeedIndexList.push_back(index);
          }
      }

    return;
  }

  template<typename SourceImageType, typename LabelImageType>
  void
  SeededOptimalPathFilter3D<SourceImageType, LabelImageType>::BackTrace()
  {
    m_backTraceLabelImage = LabelImageType::New();
    m_backTraceLabelImage->SetRegions( m_img->GetLargestPossibleRegion() );
    m_backTraceLabelImage->Allocate();
    m_backTraceLabelImage->CopyInformation(m_img);
    m_backTraceLabelImage->FillBuffer(0);

    m_backTraceLabelImageBufferPointer = m_backTraceLabelImage->GetBufferPointer();


    HeapNodeWithPre hnTmp;
    //std::cout<<"m_endSeedIndexList.size() = "<<m_endSeedIndexList.size()<<std::endl<<std::flush;

    for (std::size_t it = 0; it < m_endSeedIndexList.size(); ++it)
      {
        //std::cout<<"back tracing from the "<<it<<"-th end seed.\n"<<std::flush;

        long thisEndSeedIndex = m_endSeedIndexList[it];
        m_backTraceLabelImageBufferPointer[thisEndSeedIndex] = m_endSeedLabelValue;

        //std::cout<<"it is "<<thisEndSeedIndex<<"\n"<<std::flush;

        hnTmp = m_hpNodes[thisEndSeedIndex];

        long previousSeedIndex = hnTmp.GetPredecessorIndexValue();

        //std::cout<<"previous is "<<previousSeedIndex<<"\n"<<std::flush;

        while (m_seedImageBufferPointer[previousSeedIndex] != m_startSeedLabelValue)
          {
            //std::cout<<previousSeedIndex<<", "<<std::flush;

            m_backTraceLabelImageBufferPointer[previousSeedIndex] = m_endSeedLabelValue;

            hnTmp = m_hpNodes[previousSeedIndex];
            previousSeedIndex = hnTmp.GetPredecessorIndexValue();
          }
        //std::cout<<"\n"<<std::flush;
      }

    return;
  }


} ///< namespace FGC


#endif
