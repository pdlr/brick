/**
***************************************************************************
* @file brick/computerVision/cameraIntrinsicsPinhole.hh
*
* Include file defining a inline and template functions declared in
* CameraIntrinsicsPinhole.hh
*
* Copyright (C) 2007,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_CAMERAINTRINSICSPINHOLE_IMPL_HH
#define BRICK_COMPUTERVISION_CAMERAINTRINSICSPINHOLE_IMPL_HH

// This file is included by cameraIntrinsicsPinhole.hh, and should not be
// directly included by user code, so no need to include
// cameraIntrinsicsPinhole.hh here.
// 
// #include <brick/computerVision/cameraIntrinsicsPinhole.hh>

#include <iomanip>
#include <brick/common/expect.hh>

using namespace brick::numeric;
using namespace brick::geometry;

namespace brick {

  namespace computerVision {

    // The default constructor initializes the CameraIntrinsicsPinhole
    // instance to a consistent (but not terribly useful) state.
    template <class FloatType>
    CameraIntrinsicsPinhole<FloatType>::
    CameraIntrinsicsPinhole()
      : CameraIntrinsics<FloatType>(),
        m_centerU(50),
        m_centerV(50),
        m_focalLength(1.0),
        m_kX(1.0),
        m_kY(1.0),
        m_numPixelsX(100),
        m_numPixelsY(100)
    {
      // Empty.
    }
      

    // This constructor allows the caller to explicitly set the
    // camera intrinsic parameters.
    template <class FloatType>
    CameraIntrinsicsPinhole<FloatType>::
    CameraIntrinsicsPinhole(unsigned int numPixelsX,
                            unsigned int numPixelsY,
                            FloatType focalLength,
                            FloatType pixelSizeX,
                            FloatType pixelSizeY,
                            FloatType centerU,
                            FloatType centerV)
      : CameraIntrinsics<FloatType>(),
        m_centerU(centerU),
        m_centerV(centerV),
        m_focalLength(focalLength),
        m_kX(focalLength / pixelSizeX),
        m_kY(focalLength / pixelSizeY),
        m_numPixelsX(numPixelsX),
        m_numPixelsY(numPixelsY)
    {
      // Empty.
    }
    
    
    // This member function returns a coordinate transform that
    // "matches" *this.
    template <class FloatType>
    brick::numeric::Array2D<FloatType>
    CameraIntrinsicsPinhole<FloatType>::
    getProjectionMatrix() const
    {
      brick::numeric::Array2D<FloatType> result(3, 4);
      result[0] = m_kX; result[1] = 0.0; result[2] = m_centerU; result[3] = 0.0;
      result[4] = 0.0; result[5] = m_kY; result[6] = m_centerV; result[7] = 0.0;
      result[8] = 0.0; result[9] = 0.0; result[10] = 1.0; result[11] = 0.0;
      return result;
    }


    // This member function takes a point in 3D camera coordinates
    // and projects it into pixel coordinates.
    template <class FloatType>
    Vector2D<FloatType>
    CameraIntrinsicsPinhole<FloatType>::
    project(const brick::numeric::Vector3D<FloatType>& point) const
    {
      return Vector2D<FloatType>(m_kX * point.x() / point.z() + m_centerU,
                              m_kY * point.y() / point.z() + m_centerV);
    }
    

    // This member function sets the calibration from an input
    // stream.
    template <class FloatType>
    std::istream&
    CameraIntrinsicsPinhole<FloatType>::
    readFromStream(std::istream& stream)
    {
      // If stream is in a bad state, we can't read from it.
      if (!stream){
        return stream;
      }

      // We'll silently skip whitespace.
      common::Expect::FormatFlag flags = common::Expect::SkipWhitespace;

      // Read input data into temporary variables.
      FloatType centerU, centerV, kX, kY;
      unsigned int numpixelsX, numpixelsY;
      stream >> common::Expect("CameraIntrinsicsPinhole", flags);
      stream >> common::Expect("{", flags);
      stream >> numpixelsX;
      stream >> common::Expect(",", flags);
      stream >> numpixelsY;
      stream >> common::Expect(",", flags);
      stream >> kX;
      stream >> common::Expect(",", flags);
      stream >> kY;
      stream >> common::Expect(",", flags);
      stream >> centerU;
      stream >> common::Expect(",", flags);
      stream >> centerV;
      stream >> common::Expect("}", flags);

      // If all went well, update self.
      if(stream) {
        m_centerU = centerU;
        m_centerV = centerV;
        m_kX = kX;
        m_kY = kY;
        m_numPixelsX = numpixelsX;
        m_numPixelsY = numpixelsY;
      }
      return stream;
    }


    // This member function takes a point in 2D pixel coordinates
    // and returns a ray in 3D camera coordinates passing through
    // all of the 3D points that project to the specified 2D
    // position.
    template <class FloatType>
    geometry::Ray3D<FloatType>
    CameraIntrinsicsPinhole<FloatType>::
    reverseProject(const Vector2D<FloatType>& pixelPosition,
                   bool normalize) const
    {
      // For pinhole camera model, assume 3D point [x_cam, y_cam,
      // z_cam]^T, projection matrix
      //
      //   [[k_x, 0, x_c, 0],
      //    [0, k_y, y_c, 0],
      //    [0,   0,   1, 0]]],
      //
      // and pixel coordinate [u, v]^T.  We can hold z_cam = 1.0 and
      // write:
      //
      //   x_pix = k_x * x_cam + x_c, and
      //   y_pix = k_y * y_cam + x_c.
      // 
      // This is trivial to invert:
      //
      //   x_cam = (x_pix - x_c) / k_x, and
      //   y_cam = (y_pix - y_c) / k_y.
      //   z_cam = 1.0
      return Ray3D<FloatType>(
        Vector3D<FloatType>(0.0, 0.0, 0.0),
        Vector3D<FloatType>((pixelPosition.x() - m_centerU) / m_kX,
                         (pixelPosition.y() - m_centerV) / m_kY,
                         1.0), normalize);
    }


    // This member function writes the calibration to an
    // outputstream in a format which is compatible with member
    // function readFromStream().
    template <class FloatType>
    std::ostream&
    CameraIntrinsicsPinhole<FloatType>::
    writeToStream(std::ostream& stream) const
    {
      stream << "CameraIntrinsicsPinhole {"
             << std::fixed << std::setw(14)
             << m_numPixelsX << ", "
             << m_numPixelsY << ", "
             << m_kX << ", "
             << m_kY << ", "
             << m_centerU << ", "
             << m_centerV << "}";
      return stream;
    }

    
  } // namespace computerVision
  
} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_CAMERAINTRINSICSPINHOLE_IMPL_HH */
