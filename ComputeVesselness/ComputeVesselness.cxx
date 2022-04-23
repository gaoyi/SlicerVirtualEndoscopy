#include <algorithm>

#include "itkImageFileWriter.h"

#include "itkMinimumMaximumImageCalculator.h"

#include "itkHessian3DToVesselnessMeasureImageFilter.h"
#include "itkHessianRecursiveGaussianImageFilter.h"

#include "itkPluginUtilities.h"

#include "ComputeVesselnessCLP.h"

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

  //--------------------------------------------------------------------------------
  // functions
  ImageType::Pointer flipImageIntensity(ImageType::Pointer image);

  // Compute the vesselness. Assuming the vessle is BRIGHT wrt surrounding
  ImageType::Pointer compueVesselnessImage1(ImageType::Pointer ctImage, float sigma, float alpha1, float alpha2, double calcificationThreshold);
  ImageType::Pointer compueVesselnessMetricImage1(ImageType::Pointer ctImage, float sigma, float alpha1, float alpha2, double calcificationThreshold);
}


// end of anonymous namespace

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

  if (!vesselIsBrighter)
    {
      std::cout<<"---------flipImageIntensity(image)\n--------------"<<std::flush;

      image = flipImageIntensity(image);
    }

  //--------------------------------------------------------------------------------
  // compute vesselness image from CT
  std::cout<<"compute vesselness metric image ..."<<std::flush;
  //ImageType::Pointer vesselnessImage = compueVesselnessImage1(image, sigma, alpha1, alpha2, calcificationThreshold);
  ImageType::Pointer vesselnessMetricImage = compueVesselnessMetricImage1(image, sigma, alpha1, alpha2, calcificationThreshold);
  std::cout<<"done.\n"<<std::flush;

  //--------------------------------------------------------------------------------
  // write output
  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputVesselnessVolume.c_str() );
  writer->SetInput(vesselnessMetricImage);
  writer->SetUseCompression(0);
  writer->Update();

  return EXIT_SUCCESS;
}

namespace
{
  ImageType::Pointer compueVesselnessImage1(ImageType::Pointer ctImage, float sigma, float alpha1, float alpha2, double calcificationThreshold)
  {
    typedef itk::HessianRecursiveGaussianImageFilter< ImageType > HessianFilterType;
    typedef itk::Hessian3DToVesselnessMeasureImageFilter< FloatPixelType > VesselnessMeasureFilterType;

    HessianFilterType::Pointer hessianFilter = HessianFilterType::New();
    hessianFilter->SetInput( ctImage );
    hessianFilter->SetSigma( sigma );

    VesselnessMeasureFilterType::Pointer vesselnessFilter = VesselnessMeasureFilterType::New();
    vesselnessFilter->SetInput( hessianFilter->GetOutput() );
    vesselnessFilter->SetAlpha1( alpha1 );
    vesselnessFilter->SetAlpha2( alpha2 );
    vesselnessFilter->Update();

    //return vesselnessFilter->GetOutput();
    ImageType::Pointer vesselnessImage = vesselnessFilter->GetOutput();

    ImageType::PixelType* vesselnessImageBufferPointer = vesselnessImage->GetBufferPointer();

    long np = vesselnessImage->GetLargestPossibleRegion().GetNumberOfPixels();

    const ImageType::PixelType* ctImageBufferPtr = ctImage->GetBufferPointer();

    ImageType::PixelType v;
    ImageType::PixelType m;
    for (long it = 0; it < np; ++it)
      {
        if (ctImageBufferPtr[it] > calcificationThreshold)
          {
            vesselnessImageBufferPointer[it] = 0;
          }
      }

    return vesselnessImage;
  }

  ImageType::Pointer compueVesselnessMetricImage1(ImageType::Pointer ctImage, float sigma, float alpha1, float alpha2, double calcificationThreshold)
  {
    typedef itk::HessianRecursiveGaussianImageFilter< ImageType > HessianFilterType;
    typedef itk::Hessian3DToVesselnessMeasureImageFilter< FloatPixelType > VesselnessMeasureFilterType;

    HessianFilterType::Pointer hessianFilter = HessianFilterType::New();
    hessianFilter->SetInput( ctImage );
    hessianFilter->SetSigma( sigma );

    VesselnessMeasureFilterType::Pointer vesselnessFilter = VesselnessMeasureFilterType::New();
    vesselnessFilter->SetInput( hessianFilter->GetOutput() );
    vesselnessFilter->SetAlpha1( alpha1 );
    vesselnessFilter->SetAlpha2( alpha2 );
    vesselnessFilter->Update();

    //return vesselnessFilter->GetOutput();
    ImageType::Pointer vesselnessImage = vesselnessFilter->GetOutput();

    //return vesselnessImage;

    //--------------------------------------------------------------------------------
    // from vesselness image, compute metric image to be used in short cut.
    ImageType::Pointer metricImage = ImageType::New();
    metricImage->SetRegions(vesselnessImage->GetLargestPossibleRegion() );
    metricImage->Allocate();
    metricImage->CopyInformation(vesselnessImage);
    metricImage->FillBuffer(0);

    ImageType::PixelType* vesselnessImageBufferPointer = vesselnessImage->GetBufferPointer();
    long np = vesselnessImage->GetLargestPossibleRegion().GetNumberOfPixels();
    ImageType::PixelType maxVesselness = *std::max_element(vesselnessImageBufferPointer, vesselnessImageBufferPointer+np);

    const ImageType::PixelType* ctImageBufferPtr = ctImage->GetBufferPointer();
    ImageType::PixelType* metricImageBufferPointer = metricImage->GetBufferPointer();

    ImageType::PixelType v;
    ImageType::PixelType m;
    for (long it = 0; it < np; ++it)
      {
        v = vesselnessImageBufferPointer[it];

        //m = maxVesselness - v;

        m = exp(-(v - 13.0)*(v - 13.0)/(2*3.14*80.0));
        //m = exp(-(v - 13.0)*(v - 13.0));
        //m = atan2(30.0 - v, 10.0)/3.1415926 + 0.5;
        //m = m<0?0:m;
        m = 1.0 - m; // flip s.t. vessel has small (close to 0) value
        m *= 100; // scale to 0-100
        m += 0.01; // add small background mass

        //m = (v<20 && v > 10)?1:0;

        metricImageBufferPointer[it] = m;

        // if (ctImageBufferPtr[it] > calcificationThreshold)
        //   {
        //     metricImageBufferPointer[it] = 80;
        //   }
      }

    return metricImage;
  }


  ImageType::Pointer flipImageIntensity(ImageType::Pointer image)
  {
    typedef itk::MinimumMaximumImageCalculator<ImageType> ImageCalculatorFilterType;

    ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New ();
    imageCalculatorFilter->SetImage(image);
    imageCalculatorFilter->Compute();

    ImageType::PixelType maxValue = imageCalculatorFilter->GetMaximum();

    ImageType::Pointer flipImage = ImageType::New();
    flipImage->SetRegions(image->GetLargestPossibleRegion() );
    flipImage->Allocate();
    flipImage->CopyInformation(image);
    flipImage->FillBuffer(0);

    long np = image->GetLargestPossibleRegion().GetNumberOfPixels();

    const ImageType::PixelType* imagePtr = image->GetBufferPointer();
    ImageType::PixelType *flipImagePtr = flipImage->GetBufferPointer();

    for (long it = 0; it < np; ++it)
      {
        flipImagePtr[it] = maxValue - imagePtr[it];
      }

    return flipImage;
  }
}
