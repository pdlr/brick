/**
***************************************************************************
* @file brick/computerVision/cameraIntrinsicsPinhole.cc
*
* Source file defining a CameraIntrinsics subclass for pinhole
* cameras.
*
* Copyright (C) 2007,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#include <iomanip>
#include <brick/common/expect.hh>
#include <brick/computerVision/cameraIntrinsicsPinhole.hh>

using namespace brick::numeric;
using namespace brick::geometry;

namespace brick {

  namespace computerVision {

    // The default constructor initializes the CameraIntrinsicsPinhole
    // instance to a consistent (but not terribly useful) state.
    CameraIntrinsicsPinhole::
    CameraIntrinsicsPinhole()
      : CameraIntrinsics(),
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
    CameraIntrinsicsPinhole::
    CameraIntrinsicsPinhole(size_t numPixelsX,
                            size_t numPixelsY,
                            double focalLength,
                            double pixelSizeX,
                            double pixelSizeY,
                            double centerU,
                            double centerV)
      : CameraIntrinsics(),
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
    brick::numeric::Array2D<double>
    CameraIntrinsicsPinhole::
    getProjectionMatrix() const
    {
      brick::numeric::Array2D<double> result(3, 4);
      result[0] = m_kX; result[1] = 0.0; result[2] = m_centerU; result[3] = 0.0;
      result[4] = 0.0; result[5] = m_kY; result[6] = m_centerV; result[7] = 0.0;
      result[8] = 0.0; result[9] = 0.0; result[10] = 1.0; result[11] = 0.0;
      return result;
    }


    // This member function takes a point in 3D camera coordinates
    // and projects it into pixel coordinates.
    Vector2D<double>
    CameraIntrinsicsPinhole::
    project(const brick::numeric::Vector3D<double>& point) const
    {
      return Vector2D<double>(m_kX * point.x() / point.z() + m_centerU,
                              m_kY * point.y() / point.z() + m_centerV);
    }
    

    // This member function sets the calibration from an input
    // stream.
    std::istream&
    CameraIntrinsicsPinhole::
    readFromStream(std::istream& stream)
    {
      // If stream is in a bad state, we can't read from it.
      if (!stream){
        return stream;
      }

      // We'll silently skip whitespace.
      common::Expect::FormatFlag flags = common::Expect::SkipWhitespace;

      // Read input data into temporary variables.
      double centerU, centerV, kX, kY;
      size_t numpixelsX, numpixelsY;
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
    geometry::Ray3D<double>
    CameraIntrinsicsPinhole::
    reverseProject(const Vector2D<double>& pixelPosition,
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
      return Ray3D<double>(
        Vector3D<double>(0.0, 0.0, 0.0),
        Vector3D<double>((pixelPosition.x() - m_centerU) / m_kX,
                         (pixelPosition.y() - m_centerV) / m_kY,
                         1.0), normalize);
    }


    // This member function writes the calibration to an
    // outputstream in a format which is compatible with member
    // function readFromStream().
    std::ostream&
    CameraIntrinsicsPinhole::
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
