#ifndef ShortCutFilter3D_hxx
#define ShortCutFilter3D_hxx

#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkNumericTraits.h"
#include "itkRegionOfInterestImageFilter.h"

// local
#include "ShortCutFilter3D.h"

#include "ShortCutFilter3D.h"

namespace gth818n
{
  template<typename SourceImageType, typename LabelImageType>
  ShortCutFilter3D<SourceImageType, LabelImageType>::ShortCutFilter3D()
  {
    m_heap = NULL;
    m_hpNodes = NULL;
    m_bSegInitialized = false;

    DIST_INF = std::numeric_limits<float>::max();
    DIST_EPSION = 1e-3;
    NNGBH = 26;

    m_backgroundDistance = 1e-2;

    m_backgroundLabel = 0;

    return;
  }

  template<typename SourceImageType, typename LabelImageType>
  ShortCutFilter3D<SourceImageType, LabelImageType>::~ShortCutFilter3D()
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

    std::cout << "ShortCutFilter3D destroyed\n";

    return;
  }

  template<typename SourceImageType, typename LabelImageType>
  void
  ShortCutFilter3D<SourceImageType, LabelImageType>::SetBackgroundDistance(float dist)
  {
    if (dist <= 0.0f)
      {
        std::cerr<<"Error: dist should > 0, but got "<<dist<<std::endl;
      }

    m_backgroundDistance = dist;

    return;
  }


  template<typename SourceImageType, typename LabelImageType>
  void
  ShortCutFilter3D<SourceImageType, LabelImageType>::SetSourceImage(const SourceImageType* img)
  {
    m_img = img;

    return;
  }

  template<typename SourceImageType, typename LabelImageType>
  void
  ShortCutFilter3D<SourceImageType, LabelImageType>::SetSeedlImage(const LabelImageType* seedImage)
  {
    m_seedImage = seedImage;

    return;
  }

  template<typename SourceImageType, typename LabelImageType>
  void
  ShortCutFilter3D<SourceImageType, LabelImageType>::SetWorkMode(bool bSegUnInitialized )
  {
    m_bSegInitialized = bSegUnInitialized;

    return;
  }

  template<typename SourceImageType, typename LabelImageType>
  void
  ShortCutFilter3D<SourceImageType, LabelImageType>::InitializationAHP()
  {
    if((m_heap = new FibHeap) == NULL || (m_hpNodes = new LJHeapNode[m_DIMXYZ+1]) == NULL)
      {
        std::cerr << "Memory allocation failed-- ABORTING.\n";
        abort();
      }
    m_heap->ClearHeapOwnership();

    long i,j,k, index;
    if(!m_bSegInitialized)
      {
        m_labelImage = LabelImageType::New();
        m_labelImage->SetRegions( m_img->GetLargestPossibleRegion() );
        m_labelImage->Allocate();
        m_labelImage->CopyInformation(m_img);
        m_labelImageBufferPointer = m_labelImage->GetBufferPointer();

        m_imDist.resize(m_DIMXYZ);
        m_imLabPre.resize(m_DIMXYZ);
        m_imDistPre.resize(m_DIMXYZ);

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

        // for (int it = 0; it < m_spacingToNeighbors.size(); ++it)
        //   {
        //     std::cout<<m_spacingToNeighbors[it]<<std::endl;
        //   }
        //std::cout<<m_spacingToNeighbors.size()<<std::endl;

        // Determine neighborhood size at each vertice
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
            m_labelImageBufferPointer[index] = m_seedImageBufferPointer[index];
            if(m_labelImageBufferPointer[index] == 0)
              {
                m_hpNodes[index] = (float)DIST_INF;
                m_imDist[index] = DIST_INF;
              }
            else
              {
                m_hpNodes[index] = (float)DIST_EPSION;
                m_imDist[index] = DIST_EPSION;
              }

            m_heap->Insert(&m_hpNodes[index]);
            m_hpNodes[index].SetIndexValue(index);
          }
      }
    else
      {
        for(index = 0; index < m_DIMXYZ; index++)
          {
            if(m_seedImageBufferPointer[index] != 0 && m_seedImageBufferPointer[index] != m_imLabPre[index])
              {
                //             if(m_seedImageBufferPointer[index] != 0 && (m_imDistPre[index] != 0 ||  (m_imDistPre[index] == 0 && m_seedImageBufferPointer[index] != m_imLabPre[index]))) {
                m_hpNodes[index] = (float)DIST_EPSION;
                m_imDist[index] = DIST_EPSION;
                m_labelImageBufferPointer[index] = m_seedImageBufferPointer[index];
              }
            else
              {
                m_hpNodes[index] = (float)DIST_INF;
                m_imDist[index] = DIST_INF;
                m_labelImageBufferPointer[index] = 0;
              }
            m_heap->Insert(&m_hpNodes[index]);
            m_hpNodes[index].SetIndexValue(index);
          }
      }
  }

  template<typename SourceImageType, typename LabelImageType>
  void
  ShortCutFilter3D<SourceImageType, LabelImageType>::DijkstraBasedClassificationAHP()
  {
    LJHeapNode *hnMin, hnTmp;
    float t, tOri, tSrc;
    long i, index, indexNgbh;

    LabelPixelType labSrc;
    SourcePixelType pixCenter;

    // Insert 0 then extract it, which will balance heap
    m_heap->Insert(&hnTmp); m_heap->ExtractMin();

    long k = 0;

    // Adpative Dijkstra
    if(m_bSegInitialized)
      {
        while(!m_heap->IsEmpty())
          {
            hnMin = (LJHeapNode *) m_heap->ExtractMin();
            index = hnMin->GetIndexValue();
            tSrc = hnMin->GetKeyValue();

            // stop propagation when the new distance is larger than the previous one
            if(tSrc == DIST_INF) break;
            if(tSrc > m_imDistPre[index])
              {
                m_imDist[index] = m_imDistPre[index];
                m_labelImageBufferPointer[index] = m_imLabPre[index];
                continue;
              }

            labSrc = m_labelImageBufferPointer[index];
            m_imDist[index] = tSrc;

            // Update neighbors
            pixCenter = m_imgBufferPointer[index];
            for(i = 0; i < m_NBSIZE[index]; i++)
              {
                indexNgbh = index + m_indOff[i];
                tOri = m_imDist[indexNgbh];
                //t = 1.0;
                //t = std::abs(pixCenter - m_imgBufferPointer[indexNgbh])*m_spacingToNeighbors[i] + tSrc;
                t = (std::fabs(pixCenter - m_imgBufferPointer[indexNgbh]) + m_backgroundDistance)*m_spacingToNeighbors[i] + tSrc;
                //t = std::abs(pixCenter - m_imgBufferPointer[indexNgbh] + 1.0)*m_spacingToNeighbors[i] + tSrc;
                //t = std::abs(pixCenter - m_imgBufferPointer[indexNgbh] + 1000.0)*m_spacingToNeighbors[i] + tSrc;
                //t = std::abs(pixCenter - m_imgBufferPointer[indexNgbh]) + tSrc;
                if(tOri > t)
                  {
                    m_imDist[indexNgbh] = t;
                    m_labelImageBufferPointer[indexNgbh] = labSrc;

                    hnTmp = m_hpNodes[indexNgbh];
                    hnTmp.SetKeyValue(t);
                    m_heap->DecreaseKey(&m_hpNodes[indexNgbh], hnTmp);
                  }
              }
            k++;

            // // dbg
            // if (k % 10 == 0)
            //   {
            //     //checkLabelIntegrity();
            //   }
            // // dbg, end
          }

        // Update previous labels and distance information
        for(long i = 0; i < m_DIMXYZ; i++)
          {
            if(m_imDist[i] < DIST_INF)
              {
                m_imLabPre[i] = m_labelImageBufferPointer[i];
                m_imDistPre[i] = m_imDist[i];
              }
          }
        std::memcpy(m_labelImageBufferPointer, &m_imLabPre[0], m_DIMXYZ*sizeof(LabelPixelType));
        std::copy(m_imDistPre.begin(), m_imDistPre.end(), m_imDist.begin());

        //        m_labelImageBufferPointer = m_imLabPre;
        //        m_imDist = m_imDistPre;
      }
    // Normal Dijkstra to be used in initializing the segmenter for the current image
    else
      {
        while(!m_heap->IsEmpty())
          {
            hnMin = (LJHeapNode *) m_heap->ExtractMin();
            index = hnMin->GetIndexValue();
            tSrc = hnMin->GetKeyValue();
            labSrc = m_labelImageBufferPointer[index];
            m_imDist[index] = tSrc;

            // Update neighbors
            pixCenter = m_imgBufferPointer[index];
            for(i = 0; i < m_NBSIZE[index]; i++)
              {

                indexNgbh = index + m_indOff[i];
                tOri = m_imDist[indexNgbh];

                //t = m_spacingToNeighbors[i] + tSrc; ///< result is like voronoi diagram from the seeds.
                //t = std::abs(pixCenter - m_imgBufferPointer[indexNgbh])*m_spacingToNeighbors[i] + tSrc;
                t = (std::fabs(pixCenter - m_imgBufferPointer[indexNgbh]) + m_backgroundDistance)*m_spacingToNeighbors[i] + tSrc;
                //t = std::abs(pixCenter - m_imgBufferPointer[indexNgbh] + 1000.0)*m_spacingToNeighbors[i] + tSrc;
                //t = std::abs(pixCenter - m_imgBufferPointer[indexNgbh]) + tSrc;

                //t = std::abs(pixCenter - m_imgBufferPointer[indexNgbh]) + tSrc;
                if(tOri > t)
                  {
                    m_imDist[indexNgbh] = t;
                    m_labelImageBufferPointer[indexNgbh] = labSrc;

                    hnTmp = m_hpNodes[indexNgbh];
                    hnTmp.SetKeyValue(t);
                    m_heap->DecreaseKey(&m_hpNodes[indexNgbh], hnTmp);
                  }
              }
            k++;
          }
      }

    std::memcpy(&m_imLabPre[0], m_labelImageBufferPointer, m_DIMXYZ*sizeof(LabelPixelType));
    std::copy(m_imDist.begin(), m_imDist.end(), m_imDistPre.begin());

    // Release memory
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
  }

  template<typename SourceImageType, typename LabelImageType>
  void
  ShortCutFilter3D<SourceImageType, LabelImageType>::DoSegmentation()
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


    InitializationAHP();
    DijkstraBasedClassificationAHP();
  }


  // template<typename SourceImageType, typename LabelImageType>
  // void ShortCutFilter3D<SourceImageType, LabelImageType>
  // ::checkLabelIntegrity()
  // {
  //   for (std::size_t it = 0; it < m_labelImageBufferPointer.size(); ++it)
  //     {
  //       if (m_labelImageBufferPointer[it] > 2000)
  //         {
  //           std::cerr<<"it = "<<m_labelImageBufferPointer[it]<<std::endl;
  //           abort();
  //         }
  //     }

  //   return;
  // }


  template<typename SourceImageType, typename LabelImageType>
  typename LabelImageType::Pointer
  ShortCutFilter3D<SourceImageType, LabelImageType>::GetLabeImage()
  {
    if(m_backgroundLabel)
      {
        LabelPixelType* ptr = m_labelImage->GetBufferPointer();
        long numPixel = m_labelImage->GetLargestPossibleRegion().GetNumberOfPixels();

        for (long it = 0; it < numPixel; ++it)
          {
            if (ptr[it] == m_backgroundLabel)
              {
                ptr[it] = 0;
              }
          }
      }

    return m_labelImage;
  }


} ///< namespace FGC


#endif
