#include <algorithm>

#include <omp.h>

#include "SeededOptimalPathFilter3D.h"
#include "ShortCutFilter3D.h"

#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"

#include "itkImageRegionIterator.h"

#include "itkImageFileWriter.h"

#include "itkPluginUtilities.h"

#include "ComputeAxisFromVesselnessCLP.h"

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


  //--------------------------------------------------------------------------------
  // functions

  //--------------------------------------------------------------------------------
  // From fiducial list along the same vessel. Compute a list of label images.
  //
  // This is needed coz only two ponits is not enough to find a
  // complicated vessel and we need n points. Then, for every point i,
  // we trace back (using optimal path) to i-1. Then, we connected all
  // the sections. So the optimal path is run for n-1 times.
  //
  // But here in this function, we only constructed the label image
  // from the fiducial point list so that each i-1, i point pair forms
  // a single label image with i-1 point being 1 (starting point) and
  // i point being 2 (end point)
  std::vector<LabelImageType::Pointer> getLabelImageListFromFiducialList(ImageType::Pointer ctImage, const std::vector<std::vector<float> >& fiducialList);

  // Given an image and a label map where source is labeled as 1 and
  // end is labeled as 2, compute the optimal path
  LabelImageType::Pointer computeOptimalPathLabelImage(ImageType::Pointer metricImage, LabelImageType::Pointer inputSeedLabelImage);

  // Merge each optimal path between i-1 and i points
  LabelImageType::Pointer mergeLabelImages(std::vector<LabelImageType::Pointer> labelImageList);

  // based on the axis lable image and the ct image, compute the label image to be used by shortcut to segment the vessel
  LabelImageType::Pointer computeShortcutLabelImageFromAxisLabelImage(LabelImageType::Pointer axisLabelImage);

  // based on the axis lable image and the ct image, compute the vessel using shortcut
  LabelImageType::Pointer segmentVesselFromAxisLabelImage(ImageType::Pointer ctImage, LabelImageType::Pointer axisLabelImage);


  // based on the center line label image and a set of points
  // (DIFFERNT FROM above where all points are along the same vessle,
  // here, each point belongs to a distinct branch), find the optimal
  // path back to the center axis.
  //
  // Yes, this way, we can't have complecated branches because there
  // is only one seed point on the branch and we want to trace back to
  // the main axis. This confines that the branch is relatively not
  // too curved.
  LabelImageType::Pointer computeBranchesFromAxisAndBranchSeedPoints(ImageType::Pointer metricImage, LabelImageType::Pointer axisLabelImage, const std::vector<std::vector<float> >& fiducialList);
}

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  //--------------------------------------------------------------------------------
  // read CT image
  typedef itk::ImageFileReader<ImageType>  ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputVesselnessVolume.c_str() );
  reader->Update();

  ImageType::Pointer vesselnessMetricImage = reader->GetOutput();

  //--------------------------------------------------------------------------------
  // Get label image from fiducial list
  std::vector<LabelImageType::Pointer> seedImageList = getLabelImageListFromFiducialList(vesselnessMetricImage, fiducialsAlongCA);

  //--------------------------------------------------------------------------------
  // for each separated seed label image, trace in that one, in parallel
  std::vector<LabelImageType::Pointer> optimalPathLabelImageList(seedImageList.size());
#pragma omp parallel for
  for (std::size_t it = 0; it < seedImageList.size(); ++it)
    {
      std::cout<<"tracing in the "<<it<<" of "<<seedImageList.size()<<" images\n"<<std::flush;
      LabelImageType::Pointer thisSeedImage = seedImageList[it];

      // {
      //   char name[1000];
      //   sprintf(name, "seed_%d.nrrd", it);
      //   typedef itk::ImageFileWriter<LabelImageType> WriterType;
      //   typename WriterType::Pointer writer = WriterType::New();
      //   writer->SetFileName( name );
      //   writer->SetInput( thisSeedImage );
      //   writer->SetUseCompression(1);
      //   writer->Update();
      // }

      optimalPathLabelImageList[it] = computeOptimalPathLabelImage(vesselnessMetricImage, thisSeedImage);
    }


  //--------------------------------------------------------------------------------
  // merge resulting optimal path label images
  std::cout<<"merge sections of axis ..."<<std::flush;
  LabelImageType::Pointer caAxiaLabelImage = mergeLabelImages(optimalPathLabelImageList);
  std::cout<<"done.\n"<<std::flush;

  //--------------------------------------------------------------------------------
  // write output
  typedef itk::ImageFileWriter<LabelImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputAxisMaskVolume.c_str() );
  writer->SetInput( caAxiaLabelImage );
  writer->SetUseCompression(1);
  writer->Update();


  return EXIT_SUCCESS;
}

namespace
{
  std::vector<LabelImageType::Pointer> getLabelImageListFromFiducialList(ImageType::Pointer ctImage, const std::vector<std::vector<float> >& fiducialList)
  {
    std::vector<LabelImageType::Pointer> seedImageList(fiducialList.size() - 1);

    for( std::size_t i = 0; i < fiducialList.size() - 1; i++ )
      {
        typename LabelImageType::PointType lpsPointStart;
        lpsPointStart[0] = fiducialList[i][0];
        lpsPointStart[1] = fiducialList[i][1];
        lpsPointStart[2] = fiducialList[i][2];
        typename LabelImageType::IndexType indexStart;
        ctImage->TransformPhysicalPointToIndex(lpsPointStart, indexStart);


        typename LabelImageType::PointType lpsPointEnd;
        typename LabelImageType::IndexType indexEnd;
        lpsPointEnd[0] = fiducialList[i+1][0];
        lpsPointEnd[1] = fiducialList[i+1][1];
        lpsPointEnd[2] = fiducialList[i+1][2];
        ctImage->TransformPhysicalPointToIndex(lpsPointEnd, indexEnd);

        std::cout<<"lpsPointStart = "<<lpsPointStart<<std::endl;
        std::cout<<"indexStart = "<<indexStart<<std::endl;
        std::cout<<"lpsPointEnd = "<<lpsPointEnd<<std::endl;
        std::cout<<"indexEnd = "<<indexEnd<<std::endl<<std::flush;

        LabelImageType::Pointer thisSeedLabelImage = LabelImageType::New();
        thisSeedLabelImage->SetRegions(ctImage->GetLargestPossibleRegion() );
        thisSeedLabelImage->Allocate();
        thisSeedLabelImage->CopyInformation(ctImage);
        thisSeedLabelImage->FillBuffer(0);

        thisSeedLabelImage->SetPixel(indexStart, 1);
        thisSeedLabelImage->SetPixel(indexEnd, 2);

        seedImageList[i] = thisSeedLabelImage;
      }

    return seedImageList;
  }


  LabelImageType::Pointer computeOptimalPathLabelImage(ImageType::Pointer metricImage, LabelImageType::Pointer inputSeedLabelImage)
  {
    //----------------------------------------------------------------------
    // Use short cut to find LCA axis
    gth818n::SeededOptimalPathFilter3D<ImageType, LabelImageType> sc;
    sc.SetSourceImage(metricImage);
    sc.SetSeedlImage(inputSeedLabelImage);

    // Do Segmentation
    sc.update();

    LabelImageType::Pointer optimalPathLabelImage = sc.GetBackTraceLabelImage();

    // typedef gth818n::SeededOptimalPathFilter3D<ImageType, LabelImageType>::FloatImageType FloatImageType;
    // FloatImageType::Pointer distanceImage = sc.GetDistanceImage();
    // gth818n::writeImage<FloatImageType>(distanceImage, "distanceMapLCA.nrrd", true);

    return optimalPathLabelImage;
  }

  LabelImageType::Pointer mergeLabelImages(std::vector<LabelImageType::Pointer> labelImageList)
  {
    //--------------------------------------------------------------------------------
    // check if all images in the list share the same region
    std::size_t n = labelImageList.size();
    std::vector<const LabelImageType::PixelType*> bufferPointerList(n);
    for (std::size_t it = 0; it < n - 1; ++it)
      {
        //std::cout<<labelImageList[it]->GetLargestPossibleRegion().GetSize()<<std::endl;

        if (labelImageList[it]->GetLargestPossibleRegion() != labelImageList[it+1]->GetLargestPossibleRegion())
          {
            std::cerr<<"labelImageList[it]->GetLargestPossibleRegion() != labelImageList[it+1]->GetLargestPossibleRegion()\n";
            abort();
          }
      }

    for (std::size_t it = 0; it < n; ++it)
      {
        bufferPointerList[it] = labelImageList[it]->GetBufferPointer();
      }

    long np = static_cast<long>(labelImageList[0]->GetLargestPossibleRegion().GetNumberOfPixels());

    //----------------------------------------------------------------------
    // final output label image
    LabelImageType::Pointer mergedLabelImage = LabelImageType::New();
    mergedLabelImage->SetRegions(labelImageList[0]->GetLargestPossibleRegion() );
    mergedLabelImage->Allocate();
    mergedLabelImage->CopyInformation(labelImageList[0]);
    mergedLabelImage->FillBuffer(0);

    LabelImageType::PixelType* mergedLabelImageBufferPointer = mergedLabelImage->GetBufferPointer();

    std::cout<<"starting merging...........................\n"<<std::flush;

#pragma omp parallel for
    for (long iPixel = 0; iPixel < np; ++iPixel)
      {
        for (std::size_t it = 0; it < n; ++it)
          {
            if (bufferPointerList[it][iPixel] > 0)
              {
                mergedLabelImageBufferPointer[iPixel] = 1;
                //continue;
                break;
              }
          }
      }

    return mergedLabelImage;
  }
}
