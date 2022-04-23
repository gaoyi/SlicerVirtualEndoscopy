#ifndef SFLSRobustStatSegmentor3DLabelMap_txx_
#define SFLSRobustStatSegmentor3DLabelMap_txx_

#include "SFLSRobustStatSegmentor3DLabelMap.h"

#include <algorithm>
#include <ctime>

#include <limits>

// //debug//
#include <fstream>
// //DEBUG//

#include <omp.h>


#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionConstIterator.h"
#include "itkMinimumMaximumImageCalculator.h"


/* ============================================================   */
template< typename TPixel >
void
CSFLSRobustStatSegmentor3DLabelMap< TPixel >
::basicInit()
{
  SuperClassType::basicInit();

  m_statNeighborX = 1;
  m_statNeighborY = 1;
  m_statNeighborZ = 1;

  m_kernelWidthFactor = 10.0;

  m_inputImageIntensityMin = 0;
  m_inputImageIntensityMax = 0;

  return;
}


/* ============================================================  */
template< typename TPixel >
void
CSFLSRobustStatSegmentor3DLabelMap< TPixel >
::setInputLabelImage(TLabelImagePointer l)
{
  m_inputLabelImage = l;
  m_inputLabelImage_buffer_ptr = m_inputLabelImage->GetBufferPointer();

  TSize size = m_inputLabelImage->GetLargestPossibleRegion().GetSize();

  TIndex start = m_inputLabelImage->GetLargestPossibleRegion().GetIndex();
  TIndex origin = {{0, 0, 0}};
  if (start != origin)
    {
      std::cout<<"Warrning: Force mask start to be (0, 0, 0)\n";

      TRegion region = m_inputLabelImage->GetLargestPossibleRegion();
      region.SetIndex(origin);

      m_inputLabelImage->SetRegions(region);
    }


  if (this->m_nx + this->m_ny + this->m_nz == 0)
    {
      this->m_nx = size[0];
      this->m_ny = size[1];
      this->m_nz = size[2];
    }
  else if ( this->m_nx != (long)size[0] || this->m_ny != (long)size[1] || this->m_nz != (long)size[2] )
    {
      std::cerr<<"Error: image sizes do not match with label image size.\n";
      raise(SIGABRT);
    }

  return;
}


/* ============================================================  */
template< typename TPixel >
void
CSFLSRobustStatSegmentor3DLabelMap< TPixel >
::computeForce()
{
  double fmax = std::numeric_limits<double>::min();
  double kappaMax = std::numeric_limits<double>::min();

  long n = this->m_lz.size();
  double* kappaOnZeroLS = new double[ n ];
  double* cvForce = new double[ n ];

  std::vector<typename CSFLSLayer::iterator> m_lzIterVct( n );
  {
    long iiizzz = 0;
    for (typename CSFLSLayer::iterator itz = this->m_lz.begin(); itz != this->m_lz.end(); ++itz)
      m_lzIterVct[iiizzz++] = itz;
  }

#pragma omp parallel
  {
    double fmaxOfThisThread = std::numeric_limits<double>::min();
    double kappaMaxOfThisThread = std::numeric_limits<double>::min();

#pragma omp for
    for (long i = 0; i < n; ++i)
      {
        typename CSFLSLayer::iterator itz = m_lzIterVct[i];

        long idx = *itz;

        kappaOnZeroLS[i] = this->computeKappa(idx);

        std::vector<double> f(m_numberOfFeature);

        computeFeatureAt(idx, f);

        double a = -kernelEvaluationUsingPDF(f);

        fmaxOfThisThread = fmaxOfThisThread>fabs(a)?fmaxOfThisThread:fabs(a);
        kappaMaxOfThisThread = kappaMaxOfThisThread>fabs(kappaOnZeroLS[i])?kappaMaxOfThisThread:fabs(kappaOnZeroLS[i]);

        cvForce[i] = a;
      }

#pragma omp critical
    {
      fmax = fmax>fmaxOfThisThread?fmax:fmaxOfThisThread;
      kappaMax = kappaMax>kappaMaxOfThisThread?kappaMax:kappaMaxOfThisThread;
    }
  }


  this->m_force.resize(n);
  for (long i = 0; i < n; ++i)
    {
      //this->m_force.push_back(cvForce[i]/(fmax + 1e-10) +  (this->m_curvatureWeight)*kappaOnZeroLS[i]);
      this->m_force[i] = (1 - (this->m_curvatureWeight))*cvForce[i]/(fmax + 1e-10) + (this->m_curvatureWeight)*kappaOnZeroLS[i]/(kappaMax + 1e-10);
    }


  delete[] kappaOnZeroLS;
  delete[] cvForce;
}

/* ============================================================  */
template< typename TPixel >
void
CSFLSRobustStatSegmentor3DLabelMap< TPixel >
::inputLableImageToSeeds()
{
  for (long idx = 0; idx < this->m_nAll; ++idx)
    {
      if (m_inputLabelImage_buffer_ptr[idx])
        {
          m_seeds.push_back(idx);
        }
    }


  /*--------------------------------------------------
    Sub-sample seeds in case there are too many seeds */
  std::size_t n = m_seeds.size();
  std::cout<<"There are "<<n<<" seeds\n"<<std::flush;
  // if (n > )
  // std::random_shuffle( m_seeds.begin(), m_seeds.end() );



  return;
}

/* ============================================================  */
template< typename TPixel >
void
CSFLSRobustStatSegmentor3DLabelMap< TPixel >
::getThingsReady()
{
  /*
   1. Generate mp_mask from seeds
   2. Compute feature at each point
   3. Extract feature at/around the seeds
  */
  inputLableImageToSeeds();

  seedToMask();

  //dialteSeeds();

  initFeatureComputedImage();
  initFeatureImage();


  //computeFeature();
  getFeatureAroundSeeds();
  estimateFeatureStdDevs();

  estimatePDFs();

  return;
}


/* ============================================================ */
template< typename TPixel >
void
CSFLSRobustStatSegmentor3DLabelMap< TPixel >
::initFeatureImage()
{
  if (!(this->mp_img))
    {
      std::cerr<<"Error: set input image first.\n";
      raise(SIGABRT);
    }


  for (long ifeature = 0; ifeature < m_numberOfFeature; ++ifeature)
    {
      TFloatImagePointer fimg = TFloatImage::New();
      fimg->SetRegions(this->mp_img->GetLargestPossibleRegion() );
      fimg->Allocate();
      fimg->CopyInformation(this->mp_img);

      m_featureImageList.push_back(fimg);
      m_featureImage_buffer_ptr_list.push_back(fimg->GetBufferPointer());
    }

  return;
}


/* ============================================================ */
template< typename TPixel >
void
CSFLSRobustStatSegmentor3DLabelMap< TPixel >
::initFeatureComputedImage()
{
  if (!(this->mp_img))
    {
      std::cerr<<"Error: set input image first.\n";
      raise(SIGABRT);
    }

  m_featureComputed = TLabelImage::New();
  m_featureComputed->SetRegions(this->mp_img->GetLargestPossibleRegion());
  m_featureComputed->Allocate();
  m_featureComputed->CopyInformation(this->mp_img);
  m_featureComputed->FillBuffer(0);
  m_featureComputed_buffer_ptr = m_featureComputed->GetBufferPointer();

  return;
}


/* ============================================================ */
template< typename TPixel >
void
CSFLSRobustStatSegmentor3DLabelMap< TPixel >
::computeFeatureAt(long idx, std::vector<double>& f)
{
  f.resize(m_numberOfFeature);

  if (m_featureComputed_buffer_ptr[idx])
    {
      // the feature at this pixel is computed, just retrive

      for (long i = 0; i < m_numberOfFeature; ++i)
        {
          f[i] = *(m_featureImage_buffer_ptr_list[i] + idx);
        }
    }
  else
    {
      // compute the feature
      std::vector< double > neighborIntensities;

     long ix, iy, iz;
     this->ind2sub(idx, ix, iy, iz);

     for (long iiz = iz - m_statNeighborZ; iiz <= iz + m_statNeighborZ; ++iiz)
       {
         for (long iiy = iy - m_statNeighborY; iiy <= iy + m_statNeighborY; ++iiy)
           {
             for (long iix = ix - m_statNeighborX; iix <= ix + m_statNeighborX; ++iix)
               {
                 if (0 <= iix && iix < this->m_nx && 0 <= iiy && iiy < this->m_ny && 0 <= iiz && iiz < this->m_nz)
                   {
                     long idxx = this->sub2ind(iix, iiy, iiz);
                     neighborIntensities.push_back(this->mp_img_buffer_ptr[idxx]);
                   }
               }
           }
       }

      getRobustStatistics(neighborIntensities, f);

      for (long ifeature = 0; ifeature < m_numberOfFeature; ++ifeature)
        {
          *(m_featureImage_buffer_ptr_list[ifeature] + idx) = f[ifeature];
            //m_featureImageList[ifeature]->SetPixel(idx, f[ifeature]);
        }

      //m_featureComputed->SetPixel(idx, 1); // mark as computed
      m_featureComputed_buffer_ptr[idx] = 1;
    }

  return;
}


/* ============================================================  */
template< typename TPixel >
void
CSFLSRobustStatSegmentor3DLabelMap< TPixel >
::doSegmenation()
{
  double startingTime = clock();

  getThingsReady();

  /*============================================================
   * From the initial mask, generate: 1. SFLS, 2. mp_label and
   * 3. mp_phi.
   */
  //this->initializeSFLS();

  m_probabilityImage = TFloatImage::New();
  m_probabilityImage->SetRegions(this->mp_img->GetLargestPossibleRegion());
  m_probabilityImage->Allocate();
  m_probabilityImage->CopyInformation(this->mp_img);
  m_probabilityImage->FillBuffer(0);

  typename TFloatImage::PixelType* m_probabilityImagePtr = m_probabilityImage->GetBufferPointer();
  long n = m_probabilityImage->GetLargestPossibleRegion().GetNumberOfPixels();

  for (long idx = 0; idx < n; ++idx)
    {
      std::vector<double> featureHere(m_numberOfFeature);
      computeFeatureAt(idx, featureHere);

      m_probabilityImagePtr[idx] = kernelEvaluationUsingPDF(featureHere);
    }

  this->m_done = true;

  return;
}



/* ============================================================ */
template< typename TPixel >
void
CSFLSRobustStatSegmentor3DLabelMap< TPixel >
::getRobustStatistics(std::vector<double>& samples, std::vector<double>& robustStat)
{
  /* note, sample is sorted, so the order is changed */
  robustStat.resize(m_numberOfFeature);

  std::sort(samples.begin(), samples.end() );

  double n = samples.size();

  double q1 = n/4.0;
  double q1_floor;
  double l1 = modf(q1, &q1_floor);

  double q2 = n/2.0;
  double q2_floor;
  double l2 = modf(q2, &q2_floor);

  double q3 = 3.0*n/4.0;
  double q3_floor;
  double l3 = modf(q3, &q3_floor);

  double median = (1 - l2)*samples[static_cast<long>(q2_floor)] + l2*samples[static_cast<long>(q2_floor) + 1];

  double iqr = ( (1 - l3)*samples[static_cast<long>(q3_floor)] + l3*samples[static_cast<long>(q3_floor) + 1] ) \
    - ( (1 - l1)*samples[static_cast<long>(q1_floor)] + l1*samples[static_cast<long>(q1_floor) + 1] );

  robustStat[0] = median;
  robustStat[1] = iqr;

  /* next compute MAD */
  long nn = samples.size();
  std::vector<double> samplesDeMedian(nn);
  for (long i = 0; i < nn; ++i)
    {
      samplesDeMedian[i] = fabs(samples[i] - median);
    }

  std::sort(samplesDeMedian.begin(), samplesDeMedian.end() );

  double mad = (1 - l2)*samplesDeMedian[static_cast<long>(q2_floor)] + l2*samplesDeMedian[static_cast<long>(q2_floor) + 1];
  robustStat[2] = mad;


  return;
}


/* ============================================================ */
template< typename TPixel >
void
CSFLSRobustStatSegmentor3DLabelMap< TPixel >
::seedToMask()
{
  if (!(this->mp_img))
    {
      std::cerr<<"Error: set input image first.\n";
      raise(SIGABRT);
    }

  if (this->mp_mask)
    {
      /* Sometimes, the mask is not corresponding to the seed, like
         using the mean shape as the mask and just some other seed as
         the samples for feature. In such cases, do not touch mask
         and just return. */

      return;
    }


  long n = m_seeds.size();
  if (n == 0)
    {
      std::cerr << "Error: No seeds specified." << std::endl;
      raise(SIGABRT);
    }


  this->mp_mask = TMaskImage::New();
  this->mp_mask->SetRegions(this->mp_img->GetLargestPossibleRegion());
  this->mp_mask->Allocate();
  this->mp_mask->CopyInformation(this->mp_img);
  this->mp_mask->FillBuffer(0);
  this->mp_mask_buffer_ptr = this->mp_mask->GetBufferPointer();

  for (long i = 0; i < n; ++i)
    {
      // if (3 != m_seeds[i].size())
      //   {
      //     std::cerr<<"Error: 3 != m_seeds[i].size()\n";
      //     raise(SIGABRT);
      //   }
      long idx = m_seeds[i];

      long ix = 0;
      long iy = 0;
      long iz = 0;

      this->ind2sub(idx, ix, iy, iz);

      for (long iiz = iz - 1; iiz <= iz + 1; ++iiz)
        {
          for (long iiy = iy - 1; iiy <= iy + 1; ++iiy)
            {
              for (long iix = ix - 1; iix <= ix + 1; ++iix)
                {
                  if (0 <= iix && iix < this->m_nx && 0 <= iiy && iiy < this->m_ny && 0 <= iiz && iiz < this->m_nz)
                    {
                      //TIndex idx = {{iix, iiy, iiz}};
                      long idx = this->sub2ind(iix, iiy, iiz);

                      this->mp_mask_buffer_ptr[idx] = 1;
                    }
                }
            }
        }
    }

  return;
}


// /* ============================================================ */
// template< typename TPixel >
// void
// CSFLSRobustStatSegmentor3DLabelMap< TPixel >
// ::dialteSeeds()
// {
//   /* For each seed, add its 26 neighbors into the seed list. */

//   if (!(this->mp_img))
//     {
//       std::cerr<<"Error: set input image first.\n";
//       raise(SIGABRT);
//     }


//   long n = m_seeds.size();
//   std::vector<std::vector<long> > newSeeds;

//   if (n == 0)
//     {
//       std::cerr << "Error: No seeds specified." << std::endl;
//       raise(SIGABRT);
//     }


//   for (long i = 0; i < n; ++i)
//     {
//       long ix = m_seeds[i][0];
//       long iy = m_seeds[i][1];
//       long iz = m_seeds[i][2];

//       for (long iiz = iz - 1; iiz <= iz + 1; ++iiz)
//         {
//           for (long iiy = iy - 1; iiy <= iy + 1; ++iiy)
//             {
//               for (long iix = ix - 1; iix <= ix + 1; ++iix)
//                 {
//                   if (0 <= iix && iix < this->m_nx    \
//                       && 0 <= iiy && iiy < this->m_ny    \
//                       && 0 <= iiz && iiz < this->m_nz)
//                     {
//                       /* Some locations may be added multiple times,
//                          if the original seeds are close. But I think
//                          this is fine */

//                       std::vector<long> s(3);
//                       s[0] = iix;
//                       s[1] = iiy;
//                       s[2] = iiz;

//                       newSeeds.push_back(s);
//                     }
//                 }
//             }
//         }
//     }

//   m_seeds.assign(newSeeds.begin(), newSeeds.end() );

//   return;
// }


/* ============================================================  */
template< typename TPixel >
void
CSFLSRobustStatSegmentor3DLabelMap< TPixel >
::getFeatureAroundSeeds()
{
  if (!m_featureImageList[m_numberOfFeature-1])
    {
      // last feature image is not constructed
      std::cerr<<"Error: construct feature images first.\n";
      raise(SIGABRT);
    }

  long n = m_seeds.size();
  if (n == 0)
    {
      std::cerr << "Error: No seeds specified." << std::endl;
      raise(SIGABRT);
    }

  short ax = 0;
  short ay = 0;
  short az = 0;

  for (long i = 0; i < n; ++i)
    {
      long idx = m_seeds[i];

      long ix, iy, iz;
      this->ind2sub(idx, ix, iy, iz);

      for (long iiz = iz - az; iiz <= iz + az; ++iiz)
        {
          for (long iiy = iy - ay; iiy <= iy + ay; ++iiy)
            {
              for (long iix = ix - ax; iix <= ix + ax; ++iix)
                {
                  if (0 <= iix && iix < this->m_nx    \
                      && 0 <= iiy && iiy < this->m_ny    \
                      && 0 <= iiz && iiz < this->m_nz)
                    {
                      //TIndex idx = {{iix, iiy, iiz}};

                      long idxx = this->sub2ind(iix, iiy, iiz);

                      std::vector<double> featureHere(m_numberOfFeature);
                      computeFeatureAt(idxx, featureHere);

                      m_featureAtTheSeeds.push_back(featureHere);
                    }
                }
            }
        }
    }


  return;
}

/* ============================================================ */
template< typename TPixel >
void
CSFLSRobustStatSegmentor3DLabelMap< TPixel >
::estimateFeatureStdDevs()
{
  m_kernelStddev.assign(m_numberOfFeature, 0.0);

  long n = m_seeds.size(); // == m_featureAtTheSeeds.size()

  for (long i = 0; i < m_numberOfFeature; ++i)
    {
      double m = 0;
      for (long ii = 0; ii < n; ++ii)
        {
          m += m_featureAtTheSeeds[ii][i];
        }
      m /= n;

      for (long ii = 0; ii < n; ++ii)
        {
          m_kernelStddev[i] += (m_featureAtTheSeeds[ii][i] - m)*(m_featureAtTheSeeds[ii][i] - m);
        }

      m_kernelStddev[i] /= (n-1);
      m_kernelStddev[i] = sqrt(m_kernelStddev[i]);
    }

  return;
}

template< typename TPixel >
double
CSFLSRobustStatSegmentor3DLabelMap< TPixel >
::kernelEvaluationUsingPDF(const std::vector<double>& newFeature)
{
  double p = 1;

  for (long i = 0; i < m_numberOfFeature; ++i)
    {
      long idx = static_cast<long>(newFeature[i] - m_inputImageIntensityMin);

      double probOfThisFeature = m_PDFlearnedFromSeeds[i][idx];

      p *= probOfThisFeature;
    }

  return p;
}

/* ============================================================  */
template< typename TPixel >
double
CSFLSRobustStatSegmentor3DLabelMap< TPixel >
::kernelEvaluation(const std::vector<double>& newFeature)
{
  long n = m_seeds.size(); // == m_featureAtTheSeeds.size()

  double p = 1;
  //double p = 0;

  for (long i = 0; i < m_numberOfFeature; ++i)
    {
      double pp = 0.0;

      double stdDev = m_kernelStddev[i]/m_kernelWidthFactor; // /10 as in Eric's appendix

      double var2 = -1.0/(2*stdDev*stdDev);
      double c = 1.0/sqrt(2*(vnl_math::pi))/stdDev;

      for (long ii = 0; ii < n; ++ii)
        {
          pp += exp(var2*(newFeature[i] - m_featureAtTheSeeds[ii][i])*(newFeature[i] - m_featureAtTheSeeds[ii][i]));
        }

      pp *= c;
      pp /= n;

      p *= pp;
      //p = p>pp?p:pp;
    }

  return p;
}


/* ============================================================  */
template< typename TPixel >
void
CSFLSRobustStatSegmentor3DLabelMap< TPixel >
::setKernelWidthFactor(double f)
{
  if (f < 0.3)
    {
      m_kernelWidthFactor = 0.3;
    }

  if (f > 30.0)
    {
      m_kernelWidthFactor = 30.0;
    }

  m_kernelWidthFactor = f;


//   std::ofstream fil("/tmp/d.txt", std::ios_base::app);
//   fil<<"m_kernelWidthFactor = "<<m_kernelWidthFactor<<std::endl;
//   fil.close();


  return;
}


/* ============================================================  */
template< typename TPixel >
void
CSFLSRobustStatSegmentor3DLabelMap< TPixel >
::setIntensityHomogeneity(double h)
{
//   std::ofstream fil("/tmp/d.txt", std::ios_base::app);
//   fil<<"intensity homogeneity = "<<h<<std::endl;
//   fil.close();


  double f = h*(30.0 - 0.3) + 0.3;

  setKernelWidthFactor(f);

  return;
}

/* ============================================================  */
template< typename TPixel >
void
CSFLSRobustStatSegmentor3DLabelMap< TPixel >
::estimatePDFs()
{
  m_PDFlearnedFromSeeds.clear();

  computeMinMax(); // so we have the range of all pdfs

  long n = m_seeds.size();

  for (long ifeature = 0; ifeature < m_numberOfFeature; ++ifeature)
    {
      std::vector<double> thisPDF(m_inputImageIntensityMax - m_inputImageIntensityMin + 1);
      // assumption: TPixel are of integer types.

      double stdDev = m_kernelStddev[ifeature]/m_kernelWidthFactor; // /10 as in Eric's appendix
      double var2 = -1.0/(2*stdDev*stdDev);
      double c = 1.0/sqrt(2*(vnl_math::pi))/stdDev;

#pragma omp parallel for
      for (TPixel a = m_inputImageIntensityMin; a <= m_inputImageIntensityMax; ++a)
        {
          long ia = static_cast<long>(a - m_inputImageIntensityMin);

          double pp = 0.0;
          for (long ii = 0; ii < n; ++ii)
            {
              pp += exp(var2*(a - m_featureAtTheSeeds[ii][ifeature])*(a - m_featureAtTheSeeds[ii][ifeature]));
            }

          pp *= c;
          pp /= n;

          thisPDF[ia] = pp;
        }


      m_PDFlearnedFromSeeds.push_back(thisPDF);
    }

  return;
}


/* ============================================================  */
template< typename TPixel >
void
CSFLSRobustStatSegmentor3DLabelMap< TPixel >
::computeMinMax()
{
  if (!(this->mp_img))
    {
      std::cerr<<"Error: set input image first.\n";
      raise(SIGABRT);
    }

  typedef itk::Image<TPixel, 3> itkImage_t;

  typedef itk::MinimumMaximumImageCalculator <itkImage_t> ImageCalculatorFilterType;

  typename ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New();
  imageCalculatorFilter->SetImage(this->mp_img);
  imageCalculatorFilter->Compute();

  m_inputImageIntensityMax = imageCalculatorFilter->GetMaximum();
  m_inputImageIntensityMin = imageCalculatorFilter->GetMinimum();

  // typedef itk::ImageRegionConstIterator<itkImage_t> itkImageRegionConstIterator_t;

  // itkImageRegionConstIterator_t it((this->mp_img), (this->mp_img)->GetLargestPossibleRegion() );
  // it.GoToBegin();

  // m_inputImageIntensityMin = std::numeric_limits<unsigned>::max(); // yes, it's twisted so easity to compute.
  // m_inputImageIntensityMax = std::numeric_limits<unsigned>::min();

  // for (; !it.IsAtEnd(); ++it)
  //   {
  //     TPixel v = it.Get();

  //     m_inputImageIntensityMin = m_inputImageIntensityMin<v?m_inputImageIntensityMin:v;
  //     m_inputImageIntensityMax = m_inputImageIntensityMax>v?m_inputImageIntensityMax:v;
  //   }

  return;
}


#endif
