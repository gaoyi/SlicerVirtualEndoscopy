#ifndef SFLSSegmentor3D_h_
#define SFLSSegmentor3D_h_

#include "SFLS.h"

#include <list>
#include <vector>

//itk
#include "itkImage.h"

template< typename TPixel >
class CSFLSSegmentor3D : public CSFLS
{
public:
  typedef CSFLSSegmentor3D< TPixel > Self;

  typedef CSFLS SuperClassType;

  typedef SuperClassType::NodeType NodeType;
  typedef SuperClassType::CSFLSLayer CSFLSLayer;

  typedef itk::Image<TPixel, 3> TImage;
  typedef itk::Image<float, 3> TFloatImage;
//  typedef itk::Image<double, 3> TDoubleImage;
  typedef itk::Image<char, 3> TCharImage;
  typedef itk::Image<short, 3> TShortImage;
  typedef itk::Image<unsigned char, 3> TUCharImage;

  typedef float RealType;
  typedef itk::Image<RealType, 3> RealImageType;


  typedef TImage ImageType;
  typedef TFloatImage LSImageType;
  typedef TCharImage LabelImageType;
  typedef TUCharImage MaskImageType;
  //typedef TUCharImage NeighborLabelImageType;
  typedef TUCharImage NeighborLabelImageType;

  typedef typename TImage::IndexType TIndex;
  typedef typename TImage::SizeType TSize;
  typedef typename TImage::RegionType TRegion;

  CSFLSSegmentor3D();
  virtual ~CSFLSSegmentor3D() {}

  /* ============================================================
   * functions         */
  void basicInit();

  void setNumIter(unsigned long n);

  void setImage(typename ImageType::Pointer img);
  void setMask(typename MaskImageType::Pointer mask);

  virtual void computeForce() = 0;

  void normalizeForce();

  bool getPhiOfTheNbhdWhoIsClosestToZeroLevelInLayerCloserToZeroLevel(long idx, double& thePhi);

  void oneStepLevelSetEvolution();

  void initializeSFLS() { initializeSFLSFromMask(); }
  void initializeSFLSFromMask(); // m_insideVoxelCount is first computed here

  void initializeLabel();
  void initializePhi();

  virtual void doSegmenation() = 0;

  // geometry
  double computeKappa(long idx);


  void setMaxVolume(double v); // v is in mL
  void setMaxRunningTime(double t); // t in min


  // about evolution history
  void keepZeroLayerHistory(bool b) {m_keepZeroLayerHistory = b;}
  CSFLSLayer getZeroLayerAtIteration(unsigned long i);
  void writeZeroLayerAtIterationToFile(unsigned long i, const char* name);
  void writeZeroLayerToFile(const char* namePrefix);
    
  void setCurvatureWeight(double a);

  LSImageType::Pointer getLevelSetFunction();

  MaskImageType::Pointer getFinalMask();

  /* ============================================================
   * data     */
  typename ImageType::Pointer mp_img;
  typename ImageType::PixelType* mp_img_buffer_ptr;

  typename LabelImageType::Pointer mp_label;
  typename LabelImageType::PixelType* mp_label_buffer_ptr;

  typename MaskImageType::Pointer mp_mask; // 0, non-0 mask for object
  typename MaskImageType::PixelType* mp_mask_buffer_ptr;

  typename LSImageType::Pointer mp_phi;
  typename LSImageType::PixelType* mp_phi_buffer_ptr;

  typename NeighborLabelImageType::Pointer mp_neighborLabel;
  typename NeighborLabelImageType::PixelType* mp_neighborLabel_buffer_ptr;

  std::vector< double > m_force;

  double m_timeStep;

  unsigned long m_numIter;

protected:
  double m_curvatureWeight;

  bool m_done;

  long m_nx;
  long m_ny;
  long m_nz;
  long m_nAll;

  long m_incrementX;
  long m_incrementY;
  long m_incrementZ;

  double m_dx; // in mm
  double m_dy; // in mm
  double m_dz; // in mm

  long m_insideVoxelCount;
  double m_insideVolume;

  double m_maxVolume; // max physical volume, in mm^3
  double m_maxRunningTime; // in sec


  /*----------------------------------------------------------------------
    These two record the pts which change status 

    Because they are created and visited sequentially, and when not
    needed, are clear-ed as a whole. No random insertion or removal is
    needed. So use vector is faster than list.  */ 
  CSFLSLayer m_lIn2out;
  CSFLSLayer m_lOut2in;


  void updateInsideVoxelCount();

  void initneighborLabelImage(long nx, long ny, long nz);

  inline bool doubleEqual(double a, double b, double eps = 1e-10)
  {
    return (a-b < eps && b-a < eps);
  }


  bool m_keepZeroLayerHistory;
  std::vector< CSFLSLayer > m_zeroLayerHistory;

  static const NeighborLabelImageType::PixelType m_neighborMinusX = 1; // has minus-x neighbor
  static const NeighborLabelImageType::PixelType m_neighborPlusX = 2; // has plus-x neighbor
  static const NeighborLabelImageType::PixelType m_neighborMinusY = 4; // has minus-y neighbor
  static const NeighborLabelImageType::PixelType m_neighborPlusY = 8; // has minus-y neighbor
  static const NeighborLabelImageType::PixelType m_neighborMinusZ = 16;
  static const NeighborLabelImageType::PixelType m_neighborPlusZ = 32;

  // This function is just for help and debug, at the end of the
  // testing, this should be rare in the real code. Everything should
  // directly using linear index.
  void ind2sub(long ind, long& ix, long& iy, long& iz);
  long sub2ind(long ix, long iy, long iz);

};

#include "SFLSSegmentor3D.hxx"

#endif
