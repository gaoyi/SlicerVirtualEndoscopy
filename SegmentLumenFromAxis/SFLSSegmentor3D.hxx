#ifndef SFLSSegmentor3D_hpp_
#define SFLSSegmentor3D_hpp_

#include "SFLSSegmentor3D.h"

#include <algorithm>
#include <cmath>

#include <csignal>

#include <cassert>

#include <fstream>

template< typename TPixel >
CSFLSSegmentor3D< TPixel >
::CSFLSSegmentor3D() : CSFLS()
{
  basicInit();
}


/* ============================================================
   basicInit    */
template< typename TPixel >
void
CSFLSSegmentor3D< TPixel >
::basicInit()
{
  m_numIter = 100;
  m_timeStep = 1.0;

  m_nx = 0;
  m_ny = 0;
  m_nz = 0;
  m_nAll = 0;

  m_incrementX = 0;
  m_incrementY = 0;
  m_incrementZ = 0;

  m_dx = 1.0;
  m_dy = 1.0;
  m_dz = 1.0;

  mp_img_buffer_ptr = 0;
  mp_label_buffer_ptr = 0;
  mp_mask_buffer_ptr = 0;
  mp_phi_buffer_ptr = 0;

  m_curvatureWeight = 0.0;

  m_insideVoxelCount = 0;
  m_insideVolume = 0;

  m_maxVolume = 1e10; // in mm^3
  m_maxRunningTime = 3600; // in sec

  m_keepZeroLayerHistory = false;

  m_done = false;
}

/* ============================================================
   setNumIter    */
template< typename TPixel >
void
CSFLSSegmentor3D< TPixel >
::setNumIter(unsigned long n)
{
  m_numIter = n;
}

/* ============================================================
   setImage    */
template< typename TPixel >
void
CSFLSSegmentor3D< TPixel >
::setImage(typename ImageType::Pointer img)
{
  mp_img = img;
  mp_img_buffer_ptr = mp_img->GetBufferPointer();

  TIndex start = mp_img->GetLargestPossibleRegion().GetIndex();
  TIndex origin = {{0, 0, 0}};
  if (start != origin)
    {
      std::cout<<"Warrning: Force image start to be (0, 0, 0)\n";

      TRegion region = mp_img->GetLargestPossibleRegion();
      region.SetIndex(origin);

      mp_img->SetRegions(region);
    }

  TSize size = img->GetLargestPossibleRegion().GetSize();

  if (m_nx + m_ny + m_nz == 0)
    {
      m_nx = size[0];
      m_ny = size[1];
      m_nz = size[2];

      typename ImageType::SpacingType spc = img->GetSpacing();
      m_dx = spc[0];
      m_dy = spc[1];
      m_dz = spc[2];
    }
  else if ( m_nx != (long)size[0] || m_ny != (long)size[1] || m_nz != (long)size[2] )
    {
      std::cerr<<"image sizes do not match. abort\n";
      raise(SIGABRT);
    }

  m_nAll = m_nx*m_ny*m_nz;

  m_incrementX = 1;
  m_incrementY = m_nx;
  m_incrementZ = m_nx*m_ny;

  return;
}


/* ============================================================
   setMaxVolume    */
template< typename TPixel >
void
CSFLSSegmentor3D< TPixel >
::setMaxVolume(double v)
{
  if (v <= 0)
    {
      std::cerr<<"Error: max volume >= 0\n";
      raise(SIGABRT);
    }

  m_maxVolume = v*1000; // v is in mL, m_maxVolume is in mm^3

  return;
}

/* ============================================================ */
template< typename TPixel >
void
CSFLSSegmentor3D< TPixel >
::setMaxRunningTime(double t)
{
  if (t <= 0)
    {
      std::cerr<<"Error: t <= 0\n";
      raise(SIGABRT);
    }

  m_maxRunningTime = t*60; // t is in min, m_maxRunningTime is in second

  return;
}


/* ============================================================
   setCurvatureWeight    */
template< typename TPixel >
void
CSFLSSegmentor3D< TPixel >
::setCurvatureWeight(double a)
{
  if (a < 0)
    {
      std::cerr<<"Error: curvature weight < 0\n";
      raise(SIGABRT);
    }

  m_curvatureWeight = a;

  return;
}



/* ============================================================
   setMask    */
template< typename TPixel >
void
CSFLSSegmentor3D< TPixel >
::setMask(typename MaskImageType::Pointer mask)
{
  mp_mask = mask;
  mp_mask_buffer_ptr = mp_mask->GetBufferPointer();

  TSize size = mask->GetLargestPossibleRegion().GetSize();


  TIndex start = mp_mask->GetLargestPossibleRegion().GetIndex();
  TIndex origin = {{0, 0, 0}};
  if (start != origin)
    {
      std::cout<<"Warrning: Force mask start to be (0, 0, 0)\n";

      TRegion region = mp_mask->GetLargestPossibleRegion();
      region.SetIndex(origin);

      mp_mask->SetRegions(region);
    }

  if (m_nx + m_ny + m_nz == 0)
    {
      m_nx = size[0];
      m_ny = size[1];
      m_nz = size[2];
    }
  else if ( m_nx != (long)size[0] || m_ny != (long)size[1] || m_nz != (long)size[2] )
    {
      std::cerr<<"image sizes do not match. abort\n";
      raise(SIGABRT);
    }

  return;
}



template< typename TPixel >
bool
CSFLSSegmentor3D< TPixel >
::getPhiOfTheNbhdWhoIsClosestToZeroLevelInLayerCloserToZeroLevel(long idx, double& thePhi)
{
  /*--------------------------------------------------
   *
   * Look in all the neighbors, to find the phi value of the nbhd:
   * this nbhd should satisfy: 1. its layer is strictly closer to
   * the zero layer hence its value is thought to be updated. 2. If
   * there are several nbhd's belonging to the same layer, choose
   * the one whose phi value has the smallest abs value.  If (ix,
   * iy) is outside, go through all nbhd who is in the layer of
   * label = mylevel-1 pick the SMALLEST phi. If (ix, iy) is inside,
   * go through all nbhd who is in the layer of label = mylevel+1
   * pick the LARGEST phi.
   */
  LabelImageType::PixelType mylevel = mp_label_buffer_ptr[idx];
  NeighborLabelImageType::PixelType neighborLabel = mp_neighborLabel_buffer_ptr[idx];

  bool foundNbhd = false;

  LSImageType::PixelType itsPhi = 0.0;

  if (mylevel > 0)
    {
      // find the SMALLEST phi
      thePhi = 10000;

      if ( (neighborLabel & m_neighborPlusX) && (mp_label_buffer_ptr[idx + m_incrementX] == mylevel - 1) )
        {
          itsPhi = mp_phi_buffer_ptr[idx + m_incrementX];
          thePhi = thePhi<itsPhi?thePhi:itsPhi;

          foundNbhd = true;
        }

      if ( (neighborLabel & m_neighborMinusX) && (mp_label_buffer_ptr[idx - m_incrementX] == mylevel - 1) )
        {
          itsPhi = mp_phi_buffer_ptr[idx - m_incrementX];
          thePhi = thePhi<itsPhi?thePhi:itsPhi;

          foundNbhd = true;
        }

      if ( (neighborLabel & m_neighborPlusY) && (mp_label_buffer_ptr[idx + m_incrementY] == mylevel - 1) )
        {
          itsPhi = mp_phi_buffer_ptr[idx + m_incrementY];
          thePhi = thePhi<itsPhi?thePhi:itsPhi;

          foundNbhd = true;
        }

      if ( (neighborLabel & m_neighborMinusY) && (mp_label_buffer_ptr[idx - m_incrementY] == mylevel - 1) )
        {
          itsPhi = mp_phi_buffer_ptr[idx - m_incrementY];
          thePhi = thePhi<itsPhi?thePhi:itsPhi;

          foundNbhd = true;
        }

      if ( (neighborLabel & m_neighborPlusZ) && (mp_label_buffer_ptr[idx + m_incrementZ] == mylevel - 1) )
        {
          itsPhi = mp_phi_buffer_ptr[idx + m_incrementZ];
          thePhi = thePhi<itsPhi?thePhi:itsPhi;

          foundNbhd = true;
        }

      if ( (neighborLabel & m_neighborMinusZ) && (mp_label_buffer_ptr[idx - m_incrementZ] == mylevel - 1) )
        {
          itsPhi = mp_phi_buffer_ptr[idx - m_incrementZ];
          thePhi = thePhi<itsPhi?thePhi:itsPhi;

          foundNbhd = true;
        }
    }
  else
    {
      // find the LARGEST phi
      thePhi = -10000;

      if ( (neighborLabel & m_neighborPlusX) && (mp_label_buffer_ptr[idx + m_incrementX] == mylevel + 1) )
        {
          itsPhi = mp_phi_buffer_ptr[idx + m_incrementX];
          thePhi = thePhi>itsPhi?thePhi:itsPhi;

          foundNbhd = true;
        }

      if ( (neighborLabel & m_neighborMinusX) && (mp_label_buffer_ptr[idx - m_incrementX] == mylevel + 1) )
        {
          itsPhi = mp_phi_buffer_ptr[idx - m_incrementX];
          thePhi = thePhi>itsPhi?thePhi:itsPhi;

          foundNbhd = true;
        }

      if ( (neighborLabel & m_neighborPlusY) && (mp_label_buffer_ptr[idx + m_incrementY] == mylevel + 1) )
        {
          itsPhi = mp_phi_buffer_ptr[idx + m_incrementY];
          thePhi = thePhi>itsPhi?thePhi:itsPhi;

          foundNbhd = true;
        }

      if ( (neighborLabel & m_neighborMinusY) && (mp_label_buffer_ptr[idx - m_incrementY] == mylevel + 1) )
        {
          itsPhi = mp_phi_buffer_ptr[idx - m_incrementY];
          thePhi = thePhi>itsPhi?thePhi:itsPhi;

          foundNbhd = true;
        }

      if ( (neighborLabel & m_neighborPlusZ) && (mp_label_buffer_ptr[idx + m_incrementZ] == mylevel + 1) )
        {
          itsPhi = mp_phi_buffer_ptr[idx + m_incrementZ];
          thePhi = thePhi>itsPhi?thePhi:itsPhi;

          foundNbhd = true;
        }

      if ( (neighborLabel & m_neighborMinusZ) && (mp_label_buffer_ptr[idx - m_incrementZ] == mylevel + 1) )
        {
          itsPhi = mp_phi_buffer_ptr[idx - m_incrementZ];
          thePhi = thePhi>itsPhi?thePhi:itsPhi;

          foundNbhd = true;
        }
    }

  return foundNbhd;
}


/* ============================================================
   normalizeForce
   Normalize m_force s.t. max(abs(m_force)) < 0.5 */
template< typename TPixel >
void
CSFLSSegmentor3D< TPixel >
::normalizeForce()
{
  unsigned long nLz = m_lz.size();

  if (m_force.size() != nLz )
    {
      std::cerr<<"m_force.size() = "<<m_force.size()<<std::endl;
      std::cerr<<"nLz = "<<nLz<<std::endl;

      std::cerr<<"m_force.size() != nLz, abort.\n";
      raise(SIGABRT);
    }

  double fMax = fabs( m_force.front() );

  //for (std::list<double>::const_iterator itf = m_force.begin(); itf != m_force.end(); ++itf)
  {
    long nf = m_force.size();
    //for (std::list<double>::const_iterator itf = m_force.begin(); itf != m_force.end(); ++itf)
    for (long itf = 0; itf < nf; ++itf)
      {
        double v = fabs(m_force[itf]);
        fMax = fMax>v?fMax:v;
      }
  }
  fMax /= 0.49;

  {
    long nf = m_force.size();

    //for (std::list<double>::iterator itf = m_force.begin(); itf != m_force.end(); ++itf)
    for (long itf = 0; itf < nf; ++itf)
      {
        //(*itf) /= (fMax + 1e-10);
        m_force[itf] /= (fMax + 1e-10);
      }
  }
}


/* ============================================================
   updateInsideVoxelCount    */
template< typename TPixel >
void
CSFLSSegmentor3D< TPixel >
::updateInsideVoxelCount()
{
  m_insideVoxelCount -= m_lIn2out.size();
  m_insideVoxelCount += m_lOut2in.size();

  m_insideVolume = m_insideVoxelCount*m_dx*m_dy*m_dz;

    //dbg
    std::cout<<"m_insideVolume = "<<m_insideVolume<<std::endl;
    //dbg, end

  return;
}



/* ============================================================
   oneStepLevelSetEvolution    */
template< typename TPixel >
void
CSFLSSegmentor3D< TPixel >
::oneStepLevelSetEvolution()
{
  assert( m_nAll && m_incrementX && m_incrementY && m_incrementZ);

  // create 'changing status' lists
  CSFLSLayer Sz;
  CSFLSLayer Sn1;
  CSFLSLayer Sp1;
  CSFLSLayer Sn2;
  CSFLSLayer Sp2;

  m_lIn2out.clear();
  m_lOut2in.clear();

  /*--------------------------------------------------
    1. add F to phi(Lz), create Sn1 & Sp1
    scan Lz values [-2.5 -1.5)[-1.5 -.5)[-.5 .5](.5 1.5](1.5 2.5]
    ========                */
  {
    long itf = 0;
    for (CSFLSLayer::iterator itz = m_lz.begin(); itz != m_lz.end(); ++itf) // advance of itz is inside the loop, do NOT ++itz here.
      {
        long idx = (*itz);

        double phi_old = mp_phi_buffer_ptr[idx];
        double phi_new = phi_old + m_force[itf];

        /*----------------------------------------------------------------------
          Update the lists of pt who change the state, for faster
          energy fnal computation. */
        if ( phi_old <= 0 && phi_new > 0 )
          {
            m_lIn2out.push_back(idx);
          }

        if( phi_old > 0  && phi_new <= 0)
          {
            m_lOut2in.push_back(idx);
          }

        //           // DEBUG
        //           if (phi_new > 3.1 || phi_new < -3.1)
        //             {
        //               std::cout<<"phi_old = "<<phi_old<<std::endl;
        //               std::cout<<"its lbl = "<<(int)mp_label->get(ix, iy)<<std::endl;

        //               std::cerr<<"phi_new > 3.1 || phi_new < -3.1\n";
        //               raise(SIGABRT);
        //             }

        mp_phi_buffer_ptr[idx] = phi_new;

        if(phi_new > 0.5)
          {
            Sp1.push_back(idx);
            itz = m_lz.erase(itz);
          }
        else if (phi_new < -0.5)
          {
            Sn1.push_back(idx);
            itz = m_lz.erase(itz);
          }
        else
          {
            ++itz;
          }
        /*--------------------------------------------------
          NOTE, mp_label are (should) NOT update here. They should
          be updated with Sz, Sn/p's
          --------------------------------------------------*/
      }
  }


  //     // debug
  //     labelsCoherentCheck1();


  /*--------------------------------------------------
    2. update Ln1,Lp1,Lp2,Lp2, ****in that order****

    2.1 scan Ln1 values [-2.5 -1.5)[-1.5 -.5)[-.5 .5](.5 1.5](1.5 2.5]
    ==========                     */
  for (CSFLSLayer::iterator itn1 = m_ln1.begin(); itn1 != m_ln1.end(); )
    {
      long idx = (*itn1);

      double thePhi;
      bool found = getPhiOfTheNbhdWhoIsClosestToZeroLevelInLayerCloserToZeroLevel(idx, thePhi);

      if (found)
        {
          double phi_new = thePhi-1;
          mp_phi_buffer_ptr[idx] = phi_new;

          if (phi_new >= -0.5)
            {
              Sz.push_back(idx);
              itn1 = m_ln1.erase(itn1);
            }
          else if (phi_new < -1.5)
            {
              Sn2.push_back(idx);
              itn1 = m_ln1.erase(itn1);
            }
          else
            {
              ++itn1;
            }
        }
      else
        {
          /*--------------------------------------------------
            No nbhd in inner (closer to zero contour) layer, so
            should go to Sn2. And the phi shold be further -1
          */
          Sn2.push_back(idx);
          itn1 = m_ln1.erase(itn1);

          mp_phi_buffer_ptr[idx] -= 1.0;
        }
    }



  //     // debug
  //     labelsCoherentCheck1();


  /*--------------------------------------------------
    2.2 scan Lp1 values [-2.5 -1.5)[-1.5 -.5)[-.5 .5](.5 1.5](1.5 2.5]
    ========          */
  for (CSFLSLayer::iterator itp1 = m_lp1.begin(); itp1 != m_lp1.end();)
    {
      long idx = (*itp1);

      double thePhi;
      bool found = getPhiOfTheNbhdWhoIsClosestToZeroLevelInLayerCloserToZeroLevel(idx, thePhi);

      if (found)
        {
          double phi_new = thePhi+1;
          mp_phi_buffer_ptr[idx] = phi_new;

          if (phi_new <= 0.5)
            {
              Sz.push_back(idx);
              itp1 = m_lp1.erase(itp1);
            }
          else if (phi_new > 1.5)
            {
              Sp2.push_back(idx);
              itp1 = m_lp1.erase(itp1);
            }
          else
            {
              ++itp1;
            }
        }
      else
        {
          /*--------------------------------------------------
            No nbhd in inner (closer to zero contour) layer, so
            should go to Sp2. And the phi shold be further +1
          */

          Sp2.push_back(idx);
          itp1 = m_lp1.erase(itp1);

          //mp_phi->SetPixel(idx, mp_phi->GetPixel(idx) + 1);
          mp_phi_buffer_ptr[idx] += 1.0;
        }
    }


  //     // debug
  //     labelsCoherentCheck1();



  /*--------------------------------------------------
    2.3 scan Ln2 values [-2.5 -1.5)[-1.5 -.5)[-.5 .5](.5 1.5](1.5 2.5]
    ==========                                      */
  for (CSFLSLayer::iterator itn2 = m_ln2.begin(); itn2 != m_ln2.end(); )
    {
      long idx = (*itn2);

      double thePhi;
      bool found = getPhiOfTheNbhdWhoIsClosestToZeroLevelInLayerCloserToZeroLevel(idx, thePhi);

      if (found)
        {
          double phi_new = thePhi-1;
          mp_phi_buffer_ptr[idx] = phi_new;

          if (phi_new >= -1.5)
            {
              Sn1.push_back(idx);
              itn2 = m_ln2.erase(itn2);
            }
          else if (phi_new < -2.5)
            {
              itn2 = m_ln2.erase(itn2);

              mp_phi_buffer_ptr[idx] = -3.0;
              mp_label_buffer_ptr[idx] = -3;
            }
          else
            {
              ++itn2;
            }
        }
      else
        {
          itn2 = m_ln2.erase(itn2);

          mp_phi_buffer_ptr[idx] = -3.0;
          mp_label_buffer_ptr[idx] = -3;
        }
    }


  //     // debug
  //     labelsCoherentCheck1();


  /*--------------------------------------------------
    2.4 scan Lp2 values [-2.5 -1.5)[-1.5 -.5)[-.5 .5](.5 1.5](1.5 2.5]
    ========= */
  for (CSFLSLayer::iterator itp2 = m_lp2.begin(); itp2 != m_lp2.end(); )
    {
      long idx = (*itp2);

      double thePhi;
      bool found = getPhiOfTheNbhdWhoIsClosestToZeroLevelInLayerCloserToZeroLevel(idx, thePhi);

      if (found)
        {
          double phi_new = thePhi+1.0;
          mp_phi_buffer_ptr[idx] = phi_new;

          if (phi_new <= 1.5)
            {
              Sp1.push_back(idx);
              itp2 = m_lp2.erase(itp2);
            }
          else if (phi_new > 2.5)
            {
              itp2 = m_lp2.erase(itp2);

              mp_phi_buffer_ptr[idx] = 3.0;
              mp_label_buffer_ptr[idx] = 3;
            }
          else
            {
              ++itp2;
            }
        }
      else
        {
          itp2 = m_lp2.erase(itp2);

          mp_phi_buffer_ptr[idx] = 3.0;
          mp_label_buffer_ptr[idx] = 3;
        }
    }


  //     // debug
  //     labelsCoherentCheck1();


  /*--------------------------------------------------
    3. Deal with S-lists Sz,Sn1,Sp1,Sn2,Sp2
    3.1 Scan Sz */
  for (CSFLSLayer::iterator itSz = Sz.begin(); itSz != Sz.end(); ++itSz)
    {
      long idx = (*itSz);

      m_lz.push_back(idx);
      mp_label_buffer_ptr[idx] = 0;
    }
  //     // debug
  //     labelsCoherentCheck1();


  /*--------------------------------------------------
    3.2 Scan Sn1     */
  for (CSFLSLayer::iterator itSn1 = Sn1.begin(); itSn1 != Sn1.end(); ++itSn1)
    {
      long idx = (*itSn1);

      m_ln1.push_back(idx);

      mp_label_buffer_ptr[idx] = -1;

      LSImageType::PixelType phiHere = mp_phi_buffer_ptr[idx];
      NeighborLabelImageType::PixelType neighborLabel = mp_neighborLabel_buffer_ptr[idx];

      if ( (neighborLabel & m_neighborPlusX) && doubleEqual(mp_phi_buffer_ptr[idx + m_incrementX], -3.0) )
        {
          Sn2.push_back(idx + m_incrementX);

          mp_phi_buffer_ptr[idx + m_incrementX] = phiHere - 1.0;
        }

      if ( (neighborLabel & m_neighborMinusX) && doubleEqual(mp_phi_buffer_ptr[idx - m_incrementX], -3.0) )
        {
          Sn2.push_back(idx - m_incrementX);
          mp_phi_buffer_ptr[idx - m_incrementX] = phiHere - 1.0;
        }

      if ( (neighborLabel & m_neighborPlusY) && doubleEqual(mp_phi_buffer_ptr[idx + m_incrementY], -3.0) )
        {
          Sn2.push_back(idx + m_incrementY);
          mp_phi_buffer_ptr[idx + m_incrementY] = phiHere - 1.0;
        }

      if ( (neighborLabel & m_neighborMinusY) && doubleEqual(mp_phi_buffer_ptr[idx - m_incrementY], -3.0) )
        {
          Sn2.push_back(idx - m_incrementY );
          mp_phi_buffer_ptr[idx - m_incrementY] = phiHere - 1.0;
        }

      if ( (neighborLabel & m_neighborPlusZ) && doubleEqual(mp_phi_buffer_ptr[idx + m_incrementZ], -3.0) )
        {
          Sn2.push_back(idx + m_incrementZ);
          mp_phi_buffer_ptr[idx + m_incrementZ] = phiHere - 1.0;
        }

      if ( (neighborLabel & m_neighborMinusZ) && doubleEqual(mp_phi_buffer_ptr[idx - m_incrementZ], -3.0) )
        {
          Sn2.push_back(idx - m_incrementZ );
          mp_phi_buffer_ptr[idx - m_incrementZ] = phiHere - 1.0;
        }
    }

  //     // debug
  //     labelsCoherentCheck1();


  /*--------------------------------------------------
    3.3 Scan Sp1     */
  for (CSFLSLayer::iterator itSp1 = Sp1.begin(); itSp1 != Sp1.end(); ++itSp1)
    {
      long idx = (*itSp1);

      m_lp1.push_back(idx);
      mp_label_buffer_ptr[idx] = 1;

      LSImageType::PixelType phiHere = mp_phi_buffer_ptr[idx];
      NeighborLabelImageType::PixelType neighborLabel = mp_neighborLabel_buffer_ptr[idx];

      if ( (neighborLabel & m_neighborPlusX) && doubleEqual(mp_phi_buffer_ptr[idx + m_incrementX], 3.0) )
        {
          Sp2.push_back(idx + m_incrementX);
          mp_phi_buffer_ptr[idx + m_incrementX] = phiHere + 1.0;
        }

      if ( (neighborLabel & m_neighborMinusX) && doubleEqual(mp_phi_buffer_ptr[idx - m_incrementX], 3.0) )
        {
          Sp2.push_back(idx - m_incrementX);
          mp_phi_buffer_ptr[idx - m_incrementX] = phiHere + 1.0;
        }

      if ( (neighborLabel & m_neighborPlusY) && doubleEqual(mp_phi_buffer_ptr[idx + m_incrementY], 3.0) )
        {
          Sp2.push_back(idx + m_incrementY);
          mp_phi_buffer_ptr[idx + m_incrementY] = phiHere + 1.0;
        }

      if ( (neighborLabel & m_neighborMinusY) && doubleEqual(mp_phi_buffer_ptr[idx - m_incrementY], 3.0) )
        {
          Sp2.push_back(idx - m_incrementY );
          mp_phi_buffer_ptr[idx - m_incrementY] = phiHere + 1.0;
        }

      if ( (neighborLabel & m_neighborPlusZ) && doubleEqual(mp_phi_buffer_ptr[idx + m_incrementZ], 3.0) )
        {
          Sp2.push_back(idx + m_incrementZ);
          mp_phi_buffer_ptr[idx + m_incrementZ] = phiHere + 1.0;
        }

      if ( (neighborLabel & m_neighborMinusZ) && doubleEqual(mp_phi_buffer_ptr[idx - m_incrementZ], 3.0) )
        {
          Sp2.push_back(idx - m_incrementZ );
          mp_phi_buffer_ptr[idx - m_incrementZ] = phiHere + 1.0;
        }
    }

  //     // debug
  //     labelsCoherentCheck1();


  /*--------------------------------------------------
    3.4 Scan Sn2     */
  {
    //debug
    int aaa = 0;
    for (CSFLSLayer::iterator itSn2 = Sn2.begin(); itSn2 != Sn2.end(); ++itSn2, ++aaa)
      {
        long idx = *itSn2;

        m_ln2.push_back(idx);

        mp_label_buffer_ptr[idx] = -2;
        //           // debug
        //           labelsCoherentCheck1();
      }
  }



  /*--------------------------------------------------
    3.5 Scan Sp2     */
  for (CSFLSLayer::iterator itSp2 = Sp2.begin(); itSp2 != Sp2.end(); ++itSp2)
    {
      long idx = (*itSp2);

      m_lp2.push_back(idx);

      mp_label_buffer_ptr[idx] = 2;
    }

  //     // debug
  //     labelsCoherentCheck1();
}

/*================================================================================
  initializeLabel*/
template< typename TPixel >
void
CSFLSSegmentor3D< TPixel >
::initializeLabel()
{
  if (m_nx + m_ny + m_nz == 0)
    {
      std::cerr<<"set mp_img first.\n";
      raise(SIGABRT);
    }

  //find interface and mark as 0, create Lz
  char defaultLabel = 0;

  mp_label = LabelImageType::New();
  TRegion region = mp_img->GetLargestPossibleRegion();

  mp_label->SetRegions( region );
  mp_label->Allocate();
  mp_label->CopyInformation(mp_img);
  mp_label->FillBuffer(defaultLabel);

  mp_label_buffer_ptr = mp_label->GetBufferPointer();

  return;
}


/*================================================================================
  initializePhi*/
template< typename TPixel >
void
CSFLSSegmentor3D< TPixel >
::initializePhi()
{
  if (m_nx + m_ny + m_nz == 0)
    {
      std::cerr<<"set mp_img first.\n";
      raise(SIGABRT);
    }

  double arbitraryInitPhi = 1000;

  mp_phi = LSImageType::New();
  TRegion region = mp_img->GetLargestPossibleRegion();

  mp_phi->SetRegions( region );
  mp_phi->Allocate();
  mp_phi->CopyInformation(mp_img);
  mp_phi->FillBuffer(arbitraryInitPhi);

  mp_phi_buffer_ptr = mp_phi->GetBufferPointer();

  // After phi is set-up, init neighbor label. This should be put here to
  // enfore neighbor label image is init-ed. Should NOT be put, say, after
  // calling of initializePhi(). Coz that way I may forget.
  initneighborLabelImage(m_nx, m_ny, m_nz);

  return;
}


/* ============================================================
   initializeSFLSFromMask    */
template< typename TPixel >
void
CSFLSSegmentor3D< TPixel >
::initializeSFLSFromMask()
{
  if (!mp_mask)
    {
      std::cerr<<"set mp_mask first.\n";
      raise(SIGABRT);
    }

  initializePhi();
  initializeLabel();

  for (long idx = 0; idx < m_nAll; ++idx)
    {
      typename NeighborLabelImageType::PixelType neighborLabel = mp_neighborLabel_buffer_ptr[idx];

      //mark the inside and outside of label and phi
      if( mp_mask_buffer_ptr[idx] == 0 )
        {
          mp_label_buffer_ptr[idx] = 3;
          mp_phi_buffer_ptr[idx] = 3.0;
        }
      else
        {
          mp_label_buffer_ptr[idx] = -3;
          mp_phi_buffer_ptr[idx] = -3.0;

          ++m_insideVoxelCount;

          if ( ( (neighborLabel & m_neighborMinusX) && (mp_mask_buffer_ptr[idx - m_incrementX] == 0) ) \
               || ( (neighborLabel & m_neighborPlusX) && (mp_mask_buffer_ptr[idx + m_incrementX] == 0) ) \
               || ( (neighborLabel & m_neighborMinusY) && (mp_mask_buffer_ptr[idx - m_incrementY] == 0) ) \
               || ( (neighborLabel & m_neighborPlusY) && (mp_mask_buffer_ptr[idx + m_incrementY] == 0) ) \
               || ( (neighborLabel & m_neighborMinusZ) && (mp_mask_buffer_ptr[idx - m_incrementZ] == 0) ) \
               || ( (neighborLabel & m_neighborPlusZ) && (mp_mask_buffer_ptr[idx + m_incrementZ] == 0) ) )
            {
              // long ix, iy, iz;
              // ind2sub(idx, ix, iy, iz);
              m_lz.push_back(idx);

              mp_label_buffer_ptr[idx] = 0;
              mp_phi_buffer_ptr[idx] = 0.0;
            }
        }
    }

  m_insideVolume = m_insideVoxelCount*m_dx*m_dy*m_dz;

//   //debug//
//   gth818n::writeImage3<char>(mp_label, "lable.mha");
//   gth818n::writeImage3<unsigned char>(mp_mask, "mask.mha");
//   gth818n::writeImage3<double>(mp_phi, "phi.mha");
//   exit(0);
//   //DEBUG//


  //scan Lz to create Ln1 and Lp1
  for (CSFLSLayer::const_iterator it = m_lz.begin(); it != m_lz.end(); ++it)
    {
      // long ix = (*it)[0];
      // long iy = (*it)[1];
      // long iz = (*it)[2];

      long idx = (*it);//sub2ind(ix, iy, iz);

      NeighborLabelImageType::PixelType neighborLabel = mp_neighborLabel_buffer_ptr[idx];

      if((neighborLabel & m_neighborPlusX))
        {
          //TIndex idx = {{ix+1, iy, iz}};
          long idx1 = idx + m_incrementX;

          if (mp_label_buffer_ptr[idx1] == 3)
            {
              mp_label_buffer_ptr[idx1] = 1;
              mp_phi_buffer_ptr[idx1] = 1.0;

              m_lp1.push_back( idx1 );
            }
          else if ( mp_label_buffer_ptr[idx1] == -3 )
            {
              mp_label_buffer_ptr[idx1] = -1;
              mp_phi_buffer_ptr[idx1] = -1.0;

              m_ln1.push_back( idx1 );
            }
        }

      if((neighborLabel & m_neighborMinusX))
        {
          //TIndex idx = {{ix-1, iy, iz}};
          long idx1 = idx - m_incrementX;

          if (mp_label_buffer_ptr[idx1] == 3)
            {
              mp_label_buffer_ptr[idx1] = 1;
              mp_phi_buffer_ptr[idx1] = 1.0;

              m_lp1.push_back( idx1 );
            }
          else if ( mp_label_buffer_ptr[idx1] == -3 )
            {
              mp_label_buffer_ptr[idx1] = -1;
              mp_phi_buffer_ptr[idx1] = -1.0;

              m_ln1.push_back( idx1 );
            }
        }

      if((neighborLabel & m_neighborPlusY))
        {
          //TIndex idx = {{ix, iy+1, iz}};
          long idx1 = idx + m_incrementY;

          if (mp_label_buffer_ptr[idx1] == 3)
            {
              mp_label_buffer_ptr[idx1] = 1;
              mp_phi_buffer_ptr[idx1] = 1.0;

              m_lp1.push_back(idx1 );
            }
          else if ( mp_label_buffer_ptr[idx1] == -3 )
            {
              mp_label_buffer_ptr[idx1] = -1;
              mp_phi_buffer_ptr[idx1] = -1.0;

              m_ln1.push_back( idx1 );
            }
        }

      if((neighborLabel & m_neighborMinusY))
        {
          //TIndex idx = {{ix, iy-1, iz}};
          long idx1 = idx - m_incrementY;

          if (mp_label_buffer_ptr[idx1] == 3)
            {
              mp_label_buffer_ptr[idx1] = 1;
              mp_phi_buffer_ptr[idx1] = 1.0;

              m_lp1.push_back( idx1 );
            }
          else if ( mp_label_buffer_ptr[idx1] == -3 )
            {
              mp_label_buffer_ptr[idx1] = -1;
              mp_phi_buffer_ptr[idx1] = -1.0;

              m_ln1.push_back( idx1 );
            }
        }

      if((neighborLabel & m_neighborPlusZ))
        {
          //TIndex idx = {{ix, iy, iz+1}};
          long idx1 = idx + m_incrementZ;

          if (mp_label_buffer_ptr[idx1] == 3)
            {
              mp_label_buffer_ptr[idx1] = 1;
              mp_phi_buffer_ptr[idx1] = 1.0;

              m_lp1.push_back(idx1 );
            }
          else if ( mp_label_buffer_ptr[idx1] == -3 )
            {
              mp_label_buffer_ptr[idx1] = -1;
              mp_phi_buffer_ptr[idx1] = -1.0;

              m_ln1.push_back( idx1 );
            }
        }

      if((neighborLabel & m_neighborMinusZ))
        {
          //TIndex idx = {{ix, iy, iz-1}};
          long idx1 = idx - m_incrementZ;

          if (mp_label_buffer_ptr[idx1] == 3)
            {
              mp_label_buffer_ptr[idx1] = 1;
              mp_phi_buffer_ptr[idx1] = 1.0;

              m_lp1.push_back( idx1 );
            }
          else if ( mp_label_buffer_ptr[idx1] == -3 )
            {
              mp_label_buffer_ptr[idx1] = -1;
              mp_phi_buffer_ptr[idx1] = -1.0;

              m_ln1.push_back( idx1 );
            }
        }
    }

  //scan Ln1 to create Ln2
  for (CSFLSLayer::const_iterator it = m_ln1.begin(); it != m_ln1.end(); ++it)
    {
      // long ix = (*it)[0];
      // long iy = (*it)[1];
      // long iz = (*it)[2];

      long idx = *it;

      NeighborLabelImageType::PixelType neighborLabel = mp_neighborLabel_buffer_ptr[idx];

      //TIndex idx1 = {{ix+1, iy, iz}};
      long idx1 = idx + m_incrementX;
      if((neighborLabel & m_neighborPlusX) && mp_label_buffer_ptr[idx1] == -3 )
        {
          mp_label_buffer_ptr[idx1] = -2;
          mp_phi_buffer_ptr[idx1] = -2.0;

          m_ln2.push_back( idx1 );
        }

      long idx2 = idx - m_incrementX;
      if((neighborLabel & m_neighborMinusX) && mp_label_buffer_ptr[idx2] == -3 )
        {
          mp_label_buffer_ptr[idx2] = -2;
          mp_phi_buffer_ptr[idx2] = -2.0;

          m_ln2.push_back( idx2 );
        }

      //TIndex idx3 = {{ix, iy+1, iz}};
      long idx3 = idx + m_incrementY;
      if((neighborLabel & m_neighborPlusY) && mp_label_buffer_ptr[idx3] == -3 )
        {
          mp_label_buffer_ptr[idx3] = -2;
          mp_phi_buffer_ptr[idx3] = -2.0;

          m_ln2.push_back( idx3 );
        }

      //TIndex idx4 = {{ix, iy-1, iz}};
      long idx4 = idx - m_incrementY;
      if((neighborLabel & m_neighborMinusY) && mp_label_buffer_ptr[idx4] == -3 )
        {
          mp_label_buffer_ptr[idx4] = -2;
          mp_phi_buffer_ptr[idx4] = -2.0;

          m_ln2.push_back( idx4 );
        }

      //TIndex idx5 = {{ix, iy, iz+1}};
      long idx5 = idx + m_incrementZ;
      if((neighborLabel & m_neighborPlusZ) && mp_label_buffer_ptr[idx5] == -3 )
        {
          mp_label_buffer_ptr[idx5] = -2;
          mp_phi_buffer_ptr[idx5] = -2.0;

          m_ln2.push_back( idx5 );
        }

      //TIndex idx6 = {{ix, iy, iz-1}};
      long idx6 = idx - m_incrementZ;
      if((neighborLabel & m_neighborMinusZ) && mp_label_buffer_ptr[idx6] == -3 )
        {
          mp_label_buffer_ptr[idx6] = -2;
          mp_phi_buffer_ptr[idx6] = -2.0;

          m_ln2.push_back( idx6 );
        }
    }

  //scan Lp1 to create Lp2
  for (CSFLSLayer::const_iterator it = m_lp1.begin(); it != m_lp1.end(); ++it)
    {
      // long ix = (*it)[0];
      // long iy = (*it)[1];
      // long iz = (*it)[2];

      long idx = *it;//sub2ind(ix, iy, iz);
      NeighborLabelImageType::PixelType neighborLabel = mp_neighborLabel_buffer_ptr[idx];

      //TIndex idx1 = {{ix+1, iy, iz}};
      long idx1 = idx + m_incrementX;
      if((neighborLabel & m_neighborPlusX) && mp_label_buffer_ptr[idx1] == 3 )
        {
          mp_label_buffer_ptr[idx1] = 2;
          mp_phi_buffer_ptr[idx1] = 2.0;

          m_lp2.push_back( idx1 );
        }

      //TIndex idx2 = {{ix-1, iy, iz}};
      long idx2 = idx - m_incrementX;
      if((neighborLabel & m_neighborMinusX) && mp_label_buffer_ptr[idx2] == 3 )
        {
          mp_label_buffer_ptr[idx2] = 2;
          mp_phi_buffer_ptr[idx2] = 2.0;

          m_lp2.push_back( idx2 );
        }

      //TIndex idx3 = {{ix, iy+1, iz}};
      long idx3 = idx + m_incrementY;
      if((neighborLabel & m_neighborPlusY) && mp_label_buffer_ptr[idx3] == 3 )
        {
          mp_label_buffer_ptr[idx3] = 2;
          mp_phi_buffer_ptr[idx3] = 2.0;

          m_lp2.push_back( idx3 );
        }

      //TIndex idx4 = {{ix, iy-1, iz}};
      long idx4 = idx - m_incrementY;
      if((neighborLabel & m_neighborMinusY) && mp_label_buffer_ptr[idx4] == 3 )
        {
          mp_label_buffer_ptr[idx4] = 2;
          mp_phi_buffer_ptr[idx4] = 2.0;

          m_lp2.push_back( idx4 );
        }

      //TIndex idx5 = {{ix, iy, iz+1}};
      long idx5 = idx + m_incrementZ;
      if((neighborLabel & m_neighborPlusZ) && mp_label_buffer_ptr[idx5] == 3 )
        {
          mp_label_buffer_ptr[idx5] = 2;
          mp_phi_buffer_ptr[idx5] = 2.0;

          m_lp2.push_back( idx5 );
        }

      //TIndex idx6 = {{ix, iy, iz-1}};
      long idx6 = idx - m_incrementZ;
      if((neighborLabel & m_neighborMinusZ) && mp_label_buffer_ptr[idx6] == 3 )
        {
          mp_label_buffer_ptr[idx6] = 2;
          mp_phi_buffer_ptr[idx6] = 2.0;

          m_lp2.push_back( idx6 );
        }
    }
}


// /* ============================================================
//    doSegmenation    */
// template< typename TPixel >
// void
// CSFLSSegmentor3D< TPixel >
// ::doSegmenation()
// {
// //   double arbitraryInitPhi = 1000;
// //   mp_phi.reset(new cArray3D< double >(m_nx, m_ny, m_nz, arbitraryInitPhi) );


//   // gth818n::saveAsImage3< double >(mp_phi, "init0.nrrd");


//   /*============================================================
//    * From the initial mask, generate: 1. SFLS, 2. mp_label and
//    * 3. mp_phi.
//    */
//   initializeSFLS();

//   //gth818n::saveAsImage3< double >(mp_phi, "initPhi.nrrd");

//   for (unsigned int it = 0; it < m_numIter; ++it)
//     {
//       std::cout<<"iteration "<<it<<"\n"<<std::flush;

//       /*--------------------------------------------------
//         Compute the force on the zero level set, NOT on the whole domain.
//         This is NOT implemented in this base class.

//         This function will compute the m_force. m_force has the same
//         size as the m_ln, indicating the change at each pixel on the
//         zero level set.
//       */
//       computeForce();


//       normalizeForce();

//       //         // debug
//       //         for (std::list< double >::const_iterator itf = this->m_force.begin(); itf != this->m_force.end(); ++itf)
//       //           {
//       //             std::cout<<(*itf)<<", ";
//       //           }
//       //         std::cout<<std::endl<<it<<std::endl<<std::endl;



//       //         //debug//
//       //         labelsCoherentCheck1();

//       oneStepLevelSetEvolution();


//       //         //debug//
//       //         std::cout<<"-----------------------"<<it<<"---------------------------"<<std::endl;
//       //         std::cout<<"lz \t ln1 \t ln2 \t lp1 \t lp2 \n";
//       //         std::cout<<m_lz.size()<<"\t"<<m_ln1.size()<<"\t"<<m_ln2.size()<<"\t"<<m_lp1.size()<<"\t"<<m_lp2.size()<<std::endl;
//       //         std::cout<<"--------------------------------------------------"<<std::endl;


//       //         // debug
//       //         labelsCoherentCheck1();


//       //        gth818n::saveAsImage3< double >(mp_phi, "temp.nrrd");

// updateInsideVoxelCount();
//     }
// }


/* getLevelSetFunction */
template< typename TPixel >
itk::Image<float, 3>::Pointer
CSFLSSegmentor3D< TPixel >
::getLevelSetFunction()
{
//   if (!m_done)
//     {
//       std::cerr<<"Error: not done.\n";
//       raise(SIGABRT);
//     }

  return mp_phi;
}


/*============================================================
  computeKappa

  Compute kappa at a point in the zero level set  */
template< typename TPixel >
double
CSFLSSegmentor3D< TPixel >
::computeKappa(long idx)
{
  double dx = 0;
  double dy = 0;
  double dz = 0;

  double dxx = 0;
  double dyy = 0;
  double dzz = 0;

  double dx2 = 0;
  double dy2 = 0;
  double dz2 = 0;

  double dxy = 0;
  double dxz = 0;
  double dyz = 0;

  char xok = 0;
  char yok = 0;
  char zok = 0;

  //long idx = sub2ind(ix, iy, iz);
  NeighborLabelImageType::PixelType neighborLabel = mp_neighborLabel_buffer_ptr[idx];

  if( (neighborLabel & m_neighborMinusX) && (neighborLabel & m_neighborPlusX) )
    {
      xok = 1;
    }

  if( (neighborLabel & m_neighborMinusY) && (neighborLabel & m_neighborPlusY) )
    {
      yok = 1;
    }

  if( (neighborLabel & m_neighborMinusZ) && (neighborLabel & m_neighborPlusZ) )
    {
      zok = 1;
    }

  float tmp = mp_phi_buffer_ptr[idx];

  /**
     test
  if (isnan(tmp))
    {
      gth818n::writeImage<float, 3>(this->mp_phi, "mp_phi.nrrd");
      raise(SIGABRT);
    }
  */

  assert(!isnan(tmp));

  if (xok)
    {
      float tmp1 = mp_phi_buffer_ptr[idx - m_incrementX];
      float tmp2 = mp_phi_buffer_ptr[idx + m_incrementX];

      dx  = (tmp2 - tmp1 )/(2.0*m_dx);
      dxx = (tmp2 - 2.0*tmp + tmp1)/(m_dx*m_dx);
      dx2 = dx*dx;
    }

  if (yok)
    {
      float tmp3 = mp_phi_buffer_ptr[idx - m_incrementY];
      float tmp4 = mp_phi_buffer_ptr[idx + m_incrementY];

      dy  = ((tmp4 - tmp3 ))/(2.0*m_dy);
      dyy = (tmp4 - 2*tmp + tmp3)/(m_dy*m_dy);
      dy2 = dy*dy;
    }

  if (zok)
    {
      float tmp5 = mp_phi_buffer_ptr[idx - m_incrementZ];
      float tmp6 = mp_phi_buffer_ptr[idx + m_incrementZ];

      dz  = ((tmp6 - tmp5 ))/(2.0*m_dz);
      dzz = (tmp6 - 2.0*tmp + tmp5)/(m_dz*m_dz);
      dz2 = dz*dz;
    }


  if(xok && yok)
    {
      dxy = 0.25*(mp_phi_buffer_ptr[idx + m_incrementX + m_incrementY]  \
                  + mp_phi_buffer_ptr[idx - m_incrementX - m_incrementY] \
                  - mp_phi_buffer_ptr[idx + m_incrementX - m_incrementY] \
                  - mp_phi_buffer_ptr[idx - m_incrementX + m_incrementY]) /(m_dx*m_dy);
    }

  if(xok && zok)
    {
      dxz = 0.25*(mp_phi_buffer_ptr[idx + m_incrementX + m_incrementZ]  \
                  + mp_phi_buffer_ptr[idx - m_incrementX - m_incrementZ] \
                  - mp_phi_buffer_ptr[idx + m_incrementX - m_incrementZ] \
                  - mp_phi_buffer_ptr[idx - m_incrementX + m_incrementZ]) /(m_dx*m_dz);
    }

  if(yok && zok)
    {
      dyz = 0.25*(mp_phi_buffer_ptr[idx + m_incrementY + m_incrementZ]  \
                  + mp_phi_buffer_ptr[idx - m_incrementY - m_incrementZ] \
                  - mp_phi_buffer_ptr[idx + m_incrementY - m_incrementZ] \
                  - mp_phi_buffer_ptr[idx - m_incrementY + m_incrementZ]) /(m_dy*m_dz);
    }

  double k = (dxx*(dy2 + dz2) + dyy*(dx2 + dz2) + dzz*(dx2 + dy2) - 2*dx*dy*dxy - 2*dx*dz*dxz - 2*dy*dz*dyz) \
    /(dx2 + dy2 + dz2 + vnl_math::eps);

  assert(!isnan(k));

  return k;
}


template <typename TPixel>
typename CSFLSSegmentor3D<TPixel>::CSFLSLayer
CSFLSSegmentor3D<TPixel>
::getZeroLayerAtIteration(unsigned long i)
{
  if (!m_keepZeroLayerHistory)
    {
      std::cerr<<"Error: no history stored.";
      std::cerr<<"By default, they are not. Set keepZeroLayerHistory(true) *before* doSegmentation to keep history.\n";
      raise(SIGABRT);
    }
  else if(i >= m_numIter)
    {
      std::cerr<<"Error: history requested, "<<i<<", exceeds number of records, "<<m_numIter<<std::endl;
      raise(SIGABRT);
    }
  else
    {
      // return the i-th zero-layer
      return m_zeroLayerHistory[i];
    }
}


template< typename TPixel >
void
CSFLSSegmentor3D< TPixel >
::writeZeroLayerAtIterationToFile(unsigned long i, const char* name)
{
  if (!m_keepZeroLayerHistory)
    {
      std::cerr<<"Error: no history stored.";
      std::cerr<<"By default, they are not. Set keepZeroLayerHistory(true) *before* doSegmentation to keep history.\n";
      raise(SIGABRT);
    }
  else if(i >= m_numIter)
    {
      std::cerr<<"Error: history requested, "<<i<<", exceeds number of records, "<<m_numIter<<std::endl;
      raise(SIGABRT);
    }
  else
    {
      // write it out
      const CSFLSLayer& requestedHistory = m_zeroLayerHistory[i];

      std::ofstream f(name);

      typename ImageType::PointType physicalPoint;
      typename ImageType::IndexType index;

      for (CSFLSLayer::const_iterator itz = requestedHistory.begin(); itz != requestedHistory.end(); ++itz)
        {
//           /* output ijk */
//           f<<*itz;
//           /* _output ijk */


          /* output physical points */
          long idx = (*itz);

          ind2sub(idx, index[0], index[1], index[2]);

          // index[0] = (*itz)[0];
          // index[1] = (*itz)[1];
          // index[2] = (*itz)[2];

          mp_img->TransformIndexToPhysicalPoint(index, physicalPoint);
          // the returned physical coord is in physicalPoint, but in
          // RAS, but Slicer displays LPS, so add - at the first two coords:
          f<<-physicalPoint[0]<<" "<<-physicalPoint[1]<<" "<<physicalPoint[2]<<"\n";
          /* _output physical points */

        }
      f.close();
    }

  return;
}

template< typename TPixel >
void
CSFLSSegmentor3D< TPixel >
::writeZeroLayerToFile(const char* namePrefix)
{
  for (unsigned long i = 0; i < m_numIter; ++i)
    {
      char thisName[1000];
      sprintf(thisName, "%s_%ld.layer", namePrefix, i);
      writeZeroLayerAtIterationToFile(i, thisName);
    }

  return;
}

template< typename TPixel >
void
CSFLSSegmentor3D< TPixel >
::ind2sub(long ind, long& ix, long& iy, long& iz)
{
  iz = ind / m_incrementZ;
  iy = (ind % m_incrementZ)/m_incrementY;
  ix = ind % m_incrementY;

  return;
}

template< typename TPixel >
long
CSFLSSegmentor3D< TPixel >
::sub2ind(long ix, long iy, long iz)
{
  return iz*m_incrementZ + iy*m_incrementY + ix;
}

template< typename TPixel >
void
CSFLSSegmentor3D< TPixel >
::initneighborLabelImage(long nx, long ny, long nz)
{
  if (nx == 0 || ny == 0 || nz == 0 || (!mp_phi))
    {
      std::cerr<<"Error: nx == 0 || ny == 0 || nz == 0 || (!mp_phi)\n";
      raise(SIGABRT);
    }

  mp_neighborLabel = NeighborLabelImageType::New();
  TRegion region = mp_phi->GetLargestPossibleRegion();

  mp_neighborLabel->SetRegions( region );
  mp_neighborLabel->Allocate();
  mp_neighborLabel->CopyInformation(mp_phi);
  mp_neighborLabel->FillBuffer(-1); // all-1 coz pixel type is unsigned. this is not a good way coz it's too tricky. should use numetric-traits max-value of this type or some thing alike

  mp_neighborLabel_buffer_ptr = mp_neighborLabel->GetBufferPointer();

  long idx = 0;
  for (long iz = 0; iz < nz; ++iz)
    {
      for (long iy = 0; iy < ny; ++iy)
        {
          for (long ix = 0; ix < nx; ++ix)
            {
              if (ix == 0)
                {
                  mp_neighborLabel_buffer_ptr[idx] &= (!m_neighborMinusX); // XNeg
                }

              if (ix == m_nx-1)
                {
                  mp_neighborLabel_buffer_ptr[idx] &= (!m_neighborPlusX); // XPos
                }

              if (iy == 0)
                {
                  mp_neighborLabel_buffer_ptr[idx] &= (!m_neighborMinusY); // YNeg
                }

              if (iy == m_ny-1)
                {
                  mp_neighborLabel_buffer_ptr[idx] &= (!m_neighborPlusY); // YPos
                }

              if (iz == 0)
                {
                  mp_neighborLabel_buffer_ptr[idx] &= (!m_neighborMinusZ); // ZNeg
                }

              if (iz == m_nz-1)
                {
                  mp_neighborLabel_buffer_ptr[idx] &= (!m_neighborPlusZ); // ZPos
                }

              ++idx;
            }
        }
    }

  return;
}


/* getLevelSetFunction */
template< typename TPixel >
typename CSFLSSegmentor3D< TPixel >::MaskImageType::Pointer
CSFLSSegmentor3D< TPixel >
::getFinalMask()
{
  if (!m_done)
    {
      std::cerr<<"Error: not done.\n";
      abort();
    }

  typename MaskImageType::Pointer finalMask = MaskImageType::New();
  finalMask->SetRegions( mp_phi->GetLargestPossibleRegion() );
  finalMask->Allocate();
  finalMask->FillBuffer(0);
  finalMask->CopyInformation(mp_phi);

  const typename LSImageType::PixelType* phiBufferPointer = mp_phi->GetBufferPointer();
  typename MaskImageType::PixelType* maskBufferPointer = finalMask->GetBufferPointer();

  for (std::size_t it = 0; it < mp_phi->GetLargestPossibleRegion().GetNumberOfPixels(); ++it)
    {
      maskBufferPointer[it] = phiBufferPointer[it]<0?1:0;
    }

  return finalMask;
}


#endif
