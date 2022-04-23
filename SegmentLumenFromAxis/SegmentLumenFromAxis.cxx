#include <algorithm>
#include <set>

#include "ShortCutFilter3D.h"

#include "itkIdentityTransform.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkResampleImageFilter.h"
//#include "itkTransformType.h"


#include "itkOtsuThresholdImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkBinaryFillholeImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"

#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"

#include "itkImageRegionIterator.h"

#include "itkLaplacianSegmentationLevelSetImageFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"

#include "itkImageFileWriter.h"

#include "SFLSRobustStatSegmentor3DLabelMap.h"
#include "itkDiscreteGaussianImageFilter.h"

#include "itkPluginUtilities.h"

#include "SegmentLumenFromAxisCLP.h"

// Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
//
namespace
{
  //--------------------------------------------------------------------------------
  // basic typedef and const
  const int ImageDimension = 3;

  typedef float FloatPixelType;
  typedef itk::Image<FloatPixelType, ImageDimension> FloatImageType;
  typedef FloatImageType ImageType;

  typedef short LabelPixelType;
  typedef itk::Image<LabelPixelType, ImageDimension> LabelImageType;

  typedef itk::Image<short, ImageDimension> ShortImageType;

  //--------------------------------------------------------------------------------
  // functions

  // based on the axis lable image and the ct image, compute the label image to be used by shortcut to segment the vessel
  LabelImageType::Pointer computeShortcutLabelImageFromAxisLabelImage(LabelImageType::Pointer axisLabelImage);
  LabelImageType::Pointer computeShortcutLabelImageFromAxisLabelImage1(LabelImageType::Pointer axisLabelImage, ImageType::Pointer ctImage);

  // based on the axis lable image and the ct image, compute the vessel using shortcut
  LabelImageType::Pointer segmentVesselFromAxisLabelImage(ImageType::Pointer ctImage, LabelImageType::Pointer axisLabelImage);

  LabelImageType::Pointer adjustLabel(LabelImageType::Pointer label, ImageType::Pointer ctImage);
  //ImageType::Pointer removeCalcifiedRegion(ImageType::Pointer ctImage);

  // based on the axis lable image and the ct image, compute the vessel using shortcut
  LabelImageType::Pointer segmentVesselFromAxisLabelImageOutputProb(ImageType::Pointer ctImage, LabelImageType::Pointer axisLabelImage);


  // upsample images to 0.1mm isotropic
  LabelImageType::Pointer superResolutionSegmentation(ImageType::Pointer ctImage, LabelImageType::Pointer lumenSegmentationLabelImageInNativeResolution);


  // template< typename TImageType >
  // typename TImageType::Pointer
  // upsampleImage(typename TImageType::Pointer img, double targetIsotrpicResolution, typename TImageType::PixelType fillInValue);

  // template< typename TImageType, typename TLabelImageType >
  // typename TLabelImageType::Pointer
  // laplaceLevelsetSegmentationAtHighResolution(typename TImageType::Pointer img, \
  //                                             typename TLabelImageType::Pointer initialBinaryLabelImage, \
  //                                             int numiter);

}

namespace gth818n
{
  template<typename image_t>
  typename image_t::RegionType
  computeNonZeroRegion(typename image_t::Pointer img)
  {
    /**
     * Given the img, compute the region where outside this region,
     * the image is all zero.
     *
     * The minx, y, z are initialized as sizeX, y, z; then, whenever
     * encounter an non-zero voxel, the minx, y, z are updated
     * accordingly. Similar for maxX, y, z except that they are
     * intialized to 0, 0, 0
     */
    typedef typename image_t::RegionType imageRegion_t;
    typedef typename image_t::IndexType imageIndex_t;
    typedef typename image_t::SizeType imageSize_t;


    imageRegion_t entireRegion = img->GetLargestPossibleRegion();

    long minX = entireRegion.GetSize()[0];
    long minY = entireRegion.GetSize()[1];
    long minZ = entireRegion.GetSize()[2];

    long maxX = 0;
    long maxY = 0;
    long maxZ = 0;

    //    std::cout<<"hahaha = "<<minX<<'\t'<<minY<<'\t'<<minZ<<'\t'<<maxX<<'\t'<<maxX<<'\t'<<maxX<<'\n';

    typedef itk::ImageRegionConstIteratorWithIndex< image_t > itkImageRegionConstIteratorWithIndex_t;

    itkImageRegionConstIteratorWithIndex_t it(img, entireRegion);

    char foundNonZero = 0;

    {
      imageIndex_t idx;
      for (it.GoToBegin(); !it.IsAtEnd(); ++it)
        {
          if (it.Get() != 0)
            {
              foundNonZero = 1;

              idx = it.GetIndex();

              minX = minX<idx[0]?minX:idx[0];
              minY = minY<idx[1]?minY:idx[1];
              minZ = minZ<idx[2]?minZ:idx[2];

              maxX = maxX>idx[0]?maxX:idx[0];
              maxY = maxY>idx[1]?maxY:idx[1];
              maxZ = maxZ>idx[2]?maxZ:idx[2];
            }
        }
    }

    imageRegion_t nonZeroRegion;

    if (1 == foundNonZero)
      {
        imageIndex_t startIdx;
        startIdx[0] = minX;
        startIdx[1] = minY;
        startIdx[2] = minZ;

        imageSize_t size;
        size[0] = maxX - minX;
        size[1] = maxY - minY;
        size[2] = maxZ - minZ;

        nonZeroRegion.SetSize( size );
        nonZeroRegion.SetIndex( startIdx );
      }
    else
      {
        imageIndex_t startIdx;
        startIdx[0] = 0;
        startIdx[1] = 0;
        startIdx[2] = 0;

        imageSize_t size;
        size[0] = entireRegion.GetSize()[0];
        size[1] = entireRegion.GetSize()[1];
        size[2] = entireRegion.GetSize()[2];

        nonZeroRegion.SetSize( size );
        nonZeroRegion.SetIndex( startIdx );
      }


    return nonZeroRegion;
  }


  template<typename image_t>
  typename image_t::Pointer
  extractROI(const image_t* img, typename image_t::RegionType region)
  {
    typedef itk::RegionOfInterestImageFilter<image_t, image_t> itkRegionOfInterestImageFilter_t;

    typename itkRegionOfInterestImageFilter_t::Pointer ROIfilter = itkRegionOfInterestImageFilter_t::New();
    ROIfilter->SetInput( img );
    ROIfilter->SetRegionOfInterest( region );
    ROIfilter->Update();

    return ROIfilter->GetOutput();
  }


  template<typename input_image_t, typename output_image_t>
  typename output_image_t::Pointer
  binarilizeImage(typename input_image_t::Pointer input,                \
                  typename input_image_t::PixelType inclusiveLowerT,    \
                  typename input_image_t::PixelType inclusiveUpperT,    \
                  typename output_image_t::PixelType insideValue,       \
                  typename output_image_t::PixelType outsideValue)
  {
    /**
     * O(x) :=    I(x) \in [lowerT, upperT] ? insideValue : outsideValue
     */

    //tst
    //   std::cout<<lowerT<<std::endl;
    //   std::cout<<upperT<<std::endl;
    //tst//

    typedef itk::BinaryThresholdImageFilter<input_image_t, output_image_t> binaryThresholdImageFilter_t;

    typename binaryThresholdImageFilter_t::Pointer thlder = binaryThresholdImageFilter_t::New();
    thlder->SetInput(input);
    thlder->SetInsideValue(insideValue);
    thlder->SetOutsideValue(outsideValue);
    thlder->SetUpperThreshold(inclusiveUpperT);
    thlder->SetLowerThreshold(inclusiveLowerT);
    thlder->Update();

    return thlder->GetOutput();
  }



  template<typename input_image_t, typename output_image_t>
  typename output_image_t::Pointer
  binarilizeImage(typename input_image_t::Pointer input,                \
                  typename input_image_t::PixelType inclusiveThld,      \
                  typename output_image_t::PixelType insideValue)
  {
    typename input_image_t::PixelType lowerT = inclusiveThld;
    typename input_image_t::PixelType upperT = static_cast<typename input_image_t::PixelType>(itk::NumericTraits< typename input_image_t::PixelType >::max() );
    typename output_image_t::PixelType outsideValue = 0;

    return binarilizeImage<input_image_t, output_image_t>(input, lowerT, upperT, insideValue, outsideValue);
  }

  template<typename image_t>
  typename image_t::RegionType
  enlargeNonZeroRegion(typename image_t::Pointer img, typename image_t::RegionType nonZeroRegion, double enlargeRatio)
  {
    typedef typename image_t::RegionType imageRegion_t;
    typedef typename image_t::IndexType imageIndex_t;
    typedef typename image_t::SizeType imageSize_t;

    imageRegion_t entireRegion = img->GetLargestPossibleRegion();
    imageSize_t entireSize = entireRegion.GetSize();

    imageIndex_t start = nonZeroRegion.GetIndex();
    imageSize_t sz = nonZeroRegion.GetSize();

    for (unsigned int id = 0; id < image_t::ImageDimension; ++id)
      {
        start[id] = std::max(0l, static_cast<long>(start[id] - static_cast<float>(sz[id])*enlargeRatio));

        sz[id] = std::min(static_cast<double>(entireSize[id] - start[id]), (1.0 + 2.0*enlargeRatio)*static_cast<float>(sz[id])*enlargeRatio);
      }

    /**********************************************************************************
    {
    //tst
      std::cout<<"\t\t start =    "<<start<<std::endl<<std::flush;
      std::cout<<"\t\t entireSize =    "<<entireSize<<std::endl<<std::flush;
      std::cout<<"\t\t entireSize[1] - start[1], 7*sz[1]/5   "<<entireSize[1] - start[1]<<'\t'<<7*sz[1]/5<<'\t'<<sz[1]<<std::endl<<std::flush;
      //tst//
    }
    **********************************************************************************/

    imageRegion_t largerRegion;
    largerRegion.SetSize( sz );
    largerRegion.SetIndex( start );

    return largerRegion;
  }


  template< typename TImageType >
  typename TImageType::Pointer
  upsampleImage(typename TImageType::Pointer img, double targetIsotrpicResolution, typename TImageType::PixelType fillInValue)
  {
    // Typedefs for the different (numerous!) elements of the "resampling"

    // Identity transform.
    // We don't want any transform on our image except rescaling which is not
    // specified by a transform but by the input/output spacing as we will see
    // later.
    // So no transform will be specified.
    typedef itk::IdentityTransform<double, TImageType::ImageDimension> TransformType;

    // If ITK resampler determines there is something to interpolate which is
    // usually the case when upscaling (!) then we must specify the interpolation
    // algorithm. In our case, we want bicubic interpolation. One way to implement
    // it is with a third order b-spline. So the type is specified here and the
    // order will be specified with a method call later on.
    typedef itk::BSplineInterpolateImageFunction<TImageType, double, double> InterpolatorType;

    // The resampler type itself.
    typedef itk::ResampleImageFilter<TImageType, TImageType> ResampleFilterType;


    // Prepare the resampler.
    // Instantiate the transform and specify it should be the id transform.
    typename TransformType::Pointer idTransform = TransformType::New();
    idTransform->SetIdentity();

    // Instantiate the b-spline interpolator and set it as the third order
    // for bicubic.
    typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
    interpolator->SetSplineOrder(3);

    // Instantiate the resampler. Wire in the transform and the interpolator.
    typename ResampleFilterType::Pointer resizeFilter = ResampleFilterType::New();
    resizeFilter->SetTransform(idTransform);
    resizeFilter->SetInterpolator(interpolator);
    resizeFilter->SetOutputDirection( img->GetDirection() );
    resizeFilter->SetDefaultPixelValue( fillInValue );
    resizeFilter->SetOutputOrigin( img->GetOrigin() );

    const typename TImageType::RegionType& inputRegion = img->GetLargestPossibleRegion();
    const typename TImageType::SizeType& vnInputSize = inputRegion.GetSize();

    const typename TImageType::SpacingType& vfInputSpacing = img->GetSpacing();

    typename TImageType::SpacingType outputSpacing;
    typename TImageType::SizeType outputSize;
    for (unsigned int it = 0; it < TImageType::ImageDimension; ++it)
      {
        double factor = vfInputSpacing[it]/targetIsotrpicResolution;
        outputSize[it] = static_cast<unsigned int>(round(vnInputSize[it]*factor));
        outputSpacing[it] = targetIsotrpicResolution;
      }

    resizeFilter->SetOutputSpacing(outputSpacing);
    resizeFilter->SetSize(outputSize);
    resizeFilter->SetInput(img);
#ifdef USE_DisableITKMultithreading
    resizeFilter->SetNumberOfThreads(1);
#endif
    resizeFilter->Update();

    return resizeFilter->GetOutput();
  }


  template< typename TImageType, typename TLabelImageType >
  typename TLabelImageType::Pointer
  laplaceLevelsetSegmentationAtHighResolution(typename TImageType::Pointer img, \
                                              typename TLabelImageType::Pointer initialBinaryLabelImage, \
                                              int numiter)
  {
    typedef itk::LaplacianSegmentationLevelSetImageFilter< TLabelImageType, TImageType > LaplacianSegmentationLevelSetImageFilterType;
    typename LaplacianSegmentationLevelSetImageFilterType::Pointer laplacianSegmentation = LaplacianSegmentationLevelSetImageFilterType::New();
    laplacianSegmentation->SetCurvatureScaling( 100.0 ); //10000.0
    laplacianSegmentation->SetPropagationScaling( -1.0 );
    laplacianSegmentation->SetMaximumRMSError( 1e-6 );
    laplacianSegmentation->SetNumberOfIterations( numiter );

    laplacianSegmentation->SetInput( initialBinaryLabelImage ); // input initla mask image
    laplacianSegmentation->SetIsoSurfaceValue( 0.5 );

    laplacianSegmentation->SetFeatureImage( img );

    laplacianSegmentation->Update();


    std::cout << "Max. no. iterations: " << laplacianSegmentation->GetNumberOfIterations() << std::endl;
    std::cout << "Max. RMS error: " << laplacianSegmentation->GetMaximumRMSError() << std::endl;
    std::cout << "No. elpased iterations: " << laplacianSegmentation->GetElapsedIterations() << std::endl;
    std::cout << "RMS change: " << laplacianSegmentation->GetRMSChange() << std::endl;
    std::cout << std::endl <<std::flush;


    typename TLabelImageType::Pointer seg                              \
      = binarilizeImage<typename LaplacianSegmentationLevelSetImageFilterType::OutputImageType, TLabelImageType>(laplacianSegmentation->GetOutput(), 0.0, 1);

    return seg;

    // itk::ImageFileWriter< InternalImageType >::Pointer speedWriter = itk::ImageFileWriter< InternalImageType >::New();
    // speedWriter->SetInput( laplacianSegmentation->GetSpeedImage() );
    // speedWriter->SetFileName( "speedImage.nrrd" );
    // speedWriter->Update();
  }


}

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  //--------------------------------------------------------------------------------
  // read CT image
  typedef itk::ImageFileReader<ImageType>  ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputVolume.c_str() );
  reader->Update();
  ImageType::Pointer image = reader->GetOutput();


  typedef itk::ImageFileReader<LabelImageType>  LabelReaderType;
  LabelReaderType::Pointer readerL = LabelReaderType::New();
  readerL->SetFileName( inputAxisLabelVolume.c_str() );
  readerL->Update();
  LabelImageType::Pointer caAxiaLabelImage = readerL->GetOutput();

  std::cout<<"Segment vessle lumen from axis label image.\n"<<std::flush;
  //ImageType::Pointer imageWithOutCalcification = removeCalcifiedRegion(image);
  LabelImageType::Pointer vesselLabelImage = segmentVesselFromAxisLabelImage(image, caAxiaLabelImage);
  //LabelImageType::Pointer vesselLabelImage = segmentVesselFromAxisLabelImageOutputProb(image, caAxiaLabelImage);

  // remove pixel whose intenstiy is below this value
  LabelImageType::PixelType* vesselLabelImagePtr = vesselLabelImage->GetBufferPointer();
  const ImageType::PixelType* imagePtr = image->GetBufferPointer();

  long np = image->GetLargestPossibleRegion().GetNumberOfPixels();

  for (long it = 0; it < np; ++it)
    {
      if (imagePtr[it] < lowerThreshold)
        {
          vesselLabelImagePtr[it] = 0;
        }

      if (imagePtr[it] > calcificationThreshold)
        {
          vesselLabelImagePtr[it] = 0;
        }
    }

  // //--------------------------------------------------------------------------------
  // // write output
  // typedef itk::ImageFileWriter<LabelImageType> WriterType;
  // //typedef itk::ImageFileWriter<ImageType> WriterType;
  // typename WriterType::Pointer writer = WriterType::New();
  // writer->SetFileName( outputLumenMaskVolume.c_str() );
  // writer->SetInput( vesselLabelImage );
  // writer->SetUseCompression(1);
  // writer->Update();

  typedef itk::BinaryFillholeImageFilter< LabelImageType > BinaryFillholeImageFilterType;
  BinaryFillholeImageFilterType::Pointer filler = BinaryFillholeImageFilterType::New();
  filler->SetInput( vesselLabelImage );
  filler->SetForegroundValue( 1 );
  filler->Update();

  // double gaussianSigma = 1.0;
  // typedef itk::DiscreteGaussianImageFilter<LabelImageType, FloatImageType> GaussianType;
  // GaussianType::Pointer gaussianFilter = GaussianType::New();
  // gaussianFilter->SetInput( filler->GetOutput() );
  // gaussianFilter->SetVariance( gaussianSigma * gaussianSigma );
  // gaussianFilter->ReleaseDataFlagOn();
  // gaussianFilter->Update();

  // typedef itk::BinaryThresholdImageFilter<FloatImageType, LabelImageType>       OutputThresholdType;
  // OutputThresholdType::Pointer outputThresholder = OutputThresholdType::New();
  // outputThresholder->SetInput( gaussianFilter->GetOutput() );
  // outputThresholder->SetLowerThreshold( 0.2 ); // 0.5 gives too thin vessel
  // outputThresholder->SetInsideValue( 1 );
  // outputThresholder->SetOutsideValue( 0 );
  // outputThresholder->ReleaseDataFlagOn();
  // outputThresholder->Update();


  LabelImageType::Pointer finalMask = filler->GetOutput();

  if (superResolution)
    {
      std::cout<<"I'm in ********************************************************************************\n"<<std::flush;
      finalMask = superResolutionSegmentation(image, finalMask);
      std::cout<<"I'm done ********************************************************************************\n"<<std::flush;
    }

  //--------------------------------------------------------------------------------
  // write output
  typedef itk::ImageFileWriter<LabelImageType> WriterType;
  //typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputLumenMaskVolume.c_str() );
  writer->SetInput( finalMask );
  writer->SetUseCompression(1);
  writer->Update();


  return EXIT_SUCCESS;
}

namespace
{
  LabelImageType::Pointer computeShortcutLabelImageFromAxisLabelImage(LabelImageType::Pointer axisLabelImage)
  {
    //--------------------------------------------------------------------------------
    // allocate sc label image. default 2 to be rejection region
    LabelImageType::Pointer scLabelImage = LabelImageType::New();
    scLabelImage->SetRegions(axisLabelImage->GetLargestPossibleRegion());
    scLabelImage->Allocate();
    scLabelImage->CopyInformation(axisLabelImage);
    scLabelImage->FillBuffer(2);

    //--------------------------------------------------------------------------------
    // dilate axis to be over-segment: so everything else must be ouside
    typedef itk::BinaryBallStructuringElement<LabelImageType::PixelType, LabelImageType::ImageDimension> StructuringElementType;
    StructuringElementType structuringElement;
    structuringElement.SetRadius(6);
    structuringElement.CreateStructuringElement();

    typedef itk::BinaryDilateImageFilter <LabelImageType, LabelImageType, StructuringElementType> BinaryDilateImageFilterType;
    BinaryDilateImageFilterType::Pointer dilateFilter = BinaryDilateImageFilterType::New();
    dilateFilter->SetInput(axisLabelImage);
    dilateFilter->SetForegroundValue(1);
    dilateFilter->SetKernel(structuringElement);
    dilateFilter->Update();

    //--------------------------------------------------------------------------------
    // put pieces in
    long np = static_cast<long>(axisLabelImage->GetLargestPossibleRegion().GetNumberOfPixels());

    LabelImageType::Pointer dilatedAxisLabelImage = dilateFilter->GetOutput();
    const LabelImageType::PixelType* dilatedAxisLabelImageBufferPointer = dilatedAxisLabelImage->GetBufferPointer();

    const LabelImageType::PixelType* axisLabelImageBufferPointer = axisLabelImage->GetBufferPointer();

    LabelImageType::PixelType* scLabelImageBufferPointer = scLabelImage->GetBufferPointer();

    for (long it = 0; it < np; ++it)
      {
        if (dilatedAxisLabelImageBufferPointer[it] > 0)
          {
            scLabelImageBufferPointer[it] = 0;
          }
        if (axisLabelImageBufferPointer[it] > 0)
          {
            scLabelImageBufferPointer[it] = 1;
          }
      }

    return scLabelImage;
  }


  LabelImageType::Pointer computeShortcutLabelImageFromAxisLabelImage1(LabelImageType::Pointer axisLabelImage, ImageType::Pointer ctImage)
  {
    //--------------------------------------------------------------------------------
    // allocate sc label image. default 2 to be rejection region
    LabelImageType::Pointer scLabelImage = LabelImageType::New();
    scLabelImage->SetRegions(axisLabelImage->GetLargestPossibleRegion());
    scLabelImage->Allocate();
    scLabelImage->CopyInformation(axisLabelImage);
    scLabelImage->FillBuffer(2);

    //--------------------------------------------------------------------------------
    // dilate axis to be over-segment: so everything else must be ouside
    typedef itk::BinaryBallStructuringElement<LabelImageType::PixelType, LabelImageType::ImageDimension> StructuringElementType;
    StructuringElementType structuringElement;
    structuringElement.SetRadius(15);
    structuringElement.CreateStructuringElement();

    typedef itk::BinaryDilateImageFilter <LabelImageType, LabelImageType, StructuringElementType> BinaryDilateImageFilterType;
    BinaryDilateImageFilterType::Pointer dilateFilter = BinaryDilateImageFilterType::New();
    dilateFilter->SetInput(axisLabelImage);
    dilateFilter->SetForegroundValue(1);
    dilateFilter->SetKernel(structuringElement);
    dilateFilter->Update();

    //--------------------------------------------------------------------------------
    // put pieces in
    long np = static_cast<long>(axisLabelImage->GetLargestPossibleRegion().GetNumberOfPixels());

    LabelImageType::Pointer dilatedAxisLabelImage = dilateFilter->GetOutput();
    const LabelImageType::PixelType* dilatedAxisLabelImageBufferPointer = dilatedAxisLabelImage->GetBufferPointer();

    const ImageType::PixelType* ctImageBufferPointer = ctImage->GetBufferPointer();

    const LabelImageType::PixelType* axisLabelImageBufferPointer = axisLabelImage->GetBufferPointer();

    LabelImageType::PixelType* scLabelImageBufferPointer = scLabelImage->GetBufferPointer();

    for (long it = 0; it < np; ++it)
      {
        if (dilatedAxisLabelImageBufferPointer[it] > 0)
          {
            scLabelImageBufferPointer[it] = 0;
          }
        if (axisLabelImageBufferPointer[it] > 0)
          {
            scLabelImageBufferPointer[it] = 1;
          }

        /// removing the calcification
        // if (ctImageBufferPointer[it] > 600)
        //   {
        //     scLabelImageBufferPointer[it] = 2;
        //   }
      }

    return scLabelImage;
  }


  // based on the axis lable image and the ct image, compute the vessel using shortcut
  LabelImageType::Pointer segmentVesselFromAxisLabelImage(ImageType::Pointer ctImage, LabelImageType::Pointer axisLabelImage)
  {
    //--------------------------------------------------------------------------------
    // as fn name indicates
    //LabelImageType::Pointer shortcutLabelImage = computeShortcutLabelImageFromAxisLabelImage(axisLabelImage);
    LabelImageType::Pointer shortcutLabelImage = computeShortcutLabelImageFromAxisLabelImage1(axisLabelImage, ctImage);

    //--------------------------------------------------------------------------------
    // as fn name indicates
    gth818n::ShortCutFilter3D<ImageType, LabelImageType> sc;
    sc.SetSourceImage(ctImage);
    sc.SetSeedlImage(shortcutLabelImage);
    sc.SetBackgroundLabel(2);

    // Do Segmentation
    sc.DoSegmentation();

    return sc.GetLabeImage();
  }


  LabelImageType::Pointer segmentVesselFromAxisLabelImageOutputProb(ImageType::Pointer ctImage, LabelImageType::Pointer axisLabelImage)
  {
    typedef short PixelType;
    typedef CSFLSRobustStatSegmentor3DLabelMap< PixelType > SFLSRobustStatSegmentor3DLabelMap_c;

    typedef itk::CastImageFilter< ImageType, ShortImageType > itkCastFilter_t;
    typename itkCastFilter_t::Pointer caster = itkCastFilter_t::New();
    caster->SetInput( ctImage );
    caster->Update();

    SFLSRobustStatSegmentor3DLabelMap_c seg;
    seg.setImage(caster->GetOutput());
    seg.setInputLabelImage(axisLabelImage);

    seg.setNumIter(100);

    seg.setIntensityHomogeneity(0.6);
    seg.setCurvatureWeight(0.1);
    seg.doSegmenation();

    ImageType::Pointer probImage = seg.m_probabilityImage;

    typedef itk::OtsuThresholdImageFilter<ImageType, LabelImageType> FilterType;
    typename FilterType::Pointer otsuFilter = FilterType::New();
    otsuFilter->SetInput(probImage);
    otsuFilter->SetInsideValue(1);
    otsuFilter->Update(); // To compute threshold
    LabelImageType::Pointer label = otsuFilter->GetOutput();
    LabelImageType::PixelType* labelPtr = label->GetBufferPointer();
    long n = label->GetLargestPossibleRegion().GetNumberOfPixels();
    for (long i = 0; i < n; ++i)
      {
        labelPtr[i] = 1 - labelPtr[i];
      }


    typedef itk::ConnectedComponentImageFilter <LabelImageType, LabelImageType > ConnectedComponentImageFilterType;
    typename ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New ();
    connected->SetInput(label);
    //connected->SetForegroundValue(label);
    connected->Update();

    label = connected->GetOutput();
    labelPtr = label->GetBufferPointer();

    std::set<LabelImageType::PixelType> labelAlongAxis;
    const LabelImageType::PixelType* axisLabelImagePtr = axisLabelImage->GetBufferPointer();


    for (long i = 0; i < n; ++i)
      {
        if (0 != axisLabelImagePtr[i])
          {
            if (0 != labelPtr[i])
              {
                labelAlongAxis.insert(labelPtr[i]);
                std::cout<<int(labelPtr[i])<<std::endl;
              }
          }
      }


    for (long i = 0; i < n; ++i)
      {
        if (0 == labelAlongAxis.count(labelPtr[i]))
          {
            labelPtr[i] = 0;
          }
        else if (labelPtr[i] != 0)
          {
            labelPtr[i] = 1;
          }

        if (0 != axisLabelImagePtr[i])
          {
            labelPtr[i] = 1;
          }
      }

    return label;
  }


  ImageType::Pointer removeCalcifiedRegion(ImageType::Pointer ctImage)
  {
    ImageType::Pointer ctImageWithOutCalcification = ImageType::New();
    ctImageWithOutCalcification->SetRegions(ctImage->GetLargestPossibleRegion());
    ctImageWithOutCalcification->Allocate();
    ctImageWithOutCalcification->CopyInformation(ctImage);
    ctImageWithOutCalcification->FillBuffer(0);

    long n = ctImage->GetLargestPossibleRegion().GetNumberOfPixels();

    const ImageType::PixelType* ctImageBufferPtr = ctImage->GetBufferPointer();
    ImageType::PixelType* ctImageWithOutCalcificationBufferPtr = ctImageWithOutCalcification->GetBufferPointer();

    ImageType::PixelType calcificationThreshold = 600;
    for (long i = 0; i < n; ++i)
      {
        ctImageWithOutCalcificationBufferPtr[i] = ctImageBufferPtr[i] > calcificationThreshold?calcificationThreshold:ctImageBufferPtr[i];
      }

    return ctImageWithOutCalcification;
  }

  LabelImageType::Pointer adjustLabel(LabelImageType::Pointer label, ImageType::Pointer ctImage)
  {
    LabelImageType::Pointer labelImageWithOutCalcification = LabelImageType::New();
    labelImageWithOutCalcification->SetRegions(label->GetLargestPossibleRegion());
    labelImageWithOutCalcification->Allocate();
    labelImageWithOutCalcification->CopyInformation(label);
    labelImageWithOutCalcification->FillBuffer(0);

    long n = label->GetLargestPossibleRegion().GetNumberOfPixels();

    const ImageType::PixelType* ctImageBufferPtr = ctImage->GetBufferPointer();
    const LabelImageType::PixelType* labelBufferPtr = label->GetBufferPointer();
    LabelImageType::PixelType* labelImageWithOutCalcificationBufferPtr = labelImageWithOutCalcification->GetBufferPointer();

    ImageType::PixelType calcificationThreshold = 600;
    for (long i = 0; i < n; ++i)
      {
        labelImageWithOutCalcificationBufferPtr[i] = ctImageBufferPtr[i] > calcificationThreshold? 0:labelBufferPtr[i];
      }

    return labelImageWithOutCalcification;
  }

  LabelImageType::Pointer superResolutionSegmentation(ImageType::Pointer ctImage, LabelImageType::Pointer lumenSegmentationLabelImageInNativeResolution)
  {
    //--------------------------------------------------------------------------------
    // Check inpute label image is 0/1
    // bool isImage01 = isImage01Binary<LabelImageType>(lumenSegmentationLabelImageInNativeResolution);
    // if (!isImage01)
    //   {
    //     std::cerr<<"Error: input label image must be 0/1.\n";

    //     writeImage<LabelImageType>(lumenSegmentationLabelImageInNativeResolution, "lumenSegmentationLabelImageInNativeResolution.nrrd");

    //     abort();
    //   }


    //--------------------------------------------------------------------------------
    // Crop to only around the target (label 1) region
    typedef typename ImageType::RegionType RegionType;
    RegionType nonzeroRegion = gth818n::computeNonZeroRegion<LabelImageType>(lumenSegmentationLabelImageInNativeResolution);
    //RegionType largerNonzeroRegion = gth818n::enlargeNonZeroRegion<LabelImageType>(lumenSegmentationLabelImageInNativeResolution, nonzeroRegion, 0.1);
    RegionType largerNonzeroRegion = nonzeroRegion;
    std::cout<<"nonzeroRegion = "<<nonzeroRegion<<std::endl<<std::flush;
    std::cout<<"largerNonzeroRegion = "<<largerNonzeroRegion<<std::endl<<std::flush;

    typename ImageType::Pointer subImageInRegion = gth818n::extractROI<ImageType>(ctImage, largerNonzeroRegion);
    typename LabelImageType::Pointer subBinaryLabelImage = gth818n::extractROI<LabelImageType>(lumenSegmentationLabelImageInNativeResolution, largerNonzeroRegion);

    std::cout<<"image size before upsample = "<<subImageInRegion->GetLargestPossibleRegion().GetSize()<<std::endl<<std::flush;

    //--------------------------------------------------------------------------------
    // Smooth binary label image to continuous image by Gaussian
    // with sigma 1.0. So in the next step of up-sampling the binary
    // image, we don't have tha alias. Then, the smoothed upsampled
    // image is thresholded with threshold 0.46. This value is
    // obtained in the test at
    // "gth818n/method/itkVtk/itkExampmles/testBinaryGaussianThenThreshold.cxx". I'm
    // loading a binary image of the hippocampus generated by this
    // program, then Guassian smooth with 1.0 sigma and then
    // threshold. 0.46 seems giving the highest Dice with the
    // original binary image. The optimal value is in [0.4, 0.5]
    typedef float RealType;
    typedef itk::Image<RealType, ImageType::ImageDimension> RealImageType;
    double gaussianSigma = 1.0;

    typedef itk::SmoothingRecursiveGaussianImageFilter<LabelImageType, RealImageType>  GaussianFilterType;
    typename GaussianFilterType::Pointer gaussianFilter = GaussianFilterType::New();
    gaussianFilter->SetInput( subBinaryLabelImage );
    gaussianFilter->SetSigma( gaussianSigma );
    gaussianFilter->Update();
    typename RealImageType::Pointer smoothSubBinaryLabelImage = gaussianFilter->GetOutput();

    std::cout<<"label image size before upsample = "<<smoothSubBinaryLabelImage->GetLargestPossibleRegion().GetSize()<<std::endl<<std::flush;
    std::cout<<"********************************************************************************\n"<<std::flush;

    //--------------------------------------------------------------------------------
    // Super dense interpolate
    double targetIsotrpicResolution = 0.2;
    //double factor = 10;

    typename ImageType::Pointer highResolutioinSubImageInRegion = gth818n::upsampleImage<ImageType>(subImageInRegion, targetIsotrpicResolution, 0);
    typename RealImageType::Pointer highResolutioinSubProbabilityImageInRegion = gth818n::upsampleImage<RealImageType>(smoothSubBinaryLabelImage, targetIsotrpicResolution, 0);
    std::cout<<"image size after upsample = "<<highResolutioinSubImageInRegion->GetLargestPossibleRegion().GetSize()<<std::endl<<std::flush;
    std::cout<<"label image size after upsample = "<<highResolutioinSubProbabilityImageInRegion->GetLargestPossibleRegion().GetSize()<<std::endl<<std::flush;

    //--------------------------------------------------------------------------------
    // Level set seg in dense space.
    // see above for the reason for 0.46
    //typename LabelImageType::Pointer initMaskHighResolution = binarilizeImage<RealImageType, LabelImageType>(highResolutioinSubProbabilityImageInRegion, 0.46); // over seg
    //typename LabelImageType::Pointer initMaskHighResolution = binarilizeImage<RealImageType, LabelImageType>(highResolutioinSubProbabilityImageInRegion, 0.7); // not bad
    typename LabelImageType::Pointer initMaskHighResolution = gth818n::binarilizeImage<RealImageType, LabelImageType>(highResolutioinSubProbabilityImageInRegion, 0.3, 1); // 0.1 over seg, 0.5 under seg
    unsigned long numberOfLevelSetIterations = 10;
    typename LabelImageType::Pointer finalMask = gth818n::laplaceLevelsetSegmentationAtHighResolution< ImageType, LabelImageType >(highResolutioinSubImageInRegion, initMaskHighResolution, numberOfLevelSetIterations);

    return finalMask;
  }




}
