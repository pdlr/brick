/**
***************************************************************************
* @file brick/computerVision/cameraIntrinsics.hh
*
* Header file declaring a parent class from which to derive classes
* that represent camera intrinsic parameters.
*
* Copyright (C) 2007,2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_CAMERAINTRINSICS_HH
#define BRICK_COMPUTERVISION_CAMERAINTRINSICS_HH

#include <brick/geometry/ray3D.hh>
#include <brick/numeric/index2D.hh>
#include <brick/numeric/vector2D.hh>
#include <brick/numeric/vector3D.hh>

namespace brick {

  namespace computerVision {

    /**
     ** This abstract base class defines an interface for classes that
     ** describe camera projection and distortion parameters.  There
     ** are two coordinate systems referenced in the documentation for
     ** this class:
     **
     **   The "camera coordinate system" has its origin at the camera
     **   focus, its Z axis pointing in the direction of camera gaze,
     **   its X axis pointing to the right (from the camera's
     **   perspective), and its Y axis pointing down.
     **
     **   The "pixel coordinate system" describes positions in the
     **   image.  The upper left corner of the image has pixel
     **   coordinates (0, 0).  The center of the upper left pixel of
     **   the image has pixel coordinates (0.5, 0.5).  The upper left
     **   corner of the next-to-leftmost pixel in the top row is at
     **   pixel coordinates (1, 0).
     **/
    template <class FloatType = double>
    class CameraIntrinsics {
    public:

      /** 
       * The default constructor currently does nothing.
       */
      CameraIntrinsics() {};


      /**
       * The destructor currently does nothing.
       */
      virtual
      ~CameraIntrinsics() {}


      /** 
       * This function should be overridden by derived classes so that
       * it takes points in 3D camera coordinates and projects them
       * into image pixel coordinates.
       * 
       * @param point This argument is the 3D point to be projected.
       * 
       * @return The return value is the resulting pixel coordinate.
       */
      virtual brick::numeric::Vector2D<FloatType>
      project(const brick::numeric::Vector3D<FloatType>& point) const = 0;


      /** 
       * This function returns a ray in 3D camera coordinates starting
       * at the camera focus, and passing through the center of the
       * specified pixel.
       * 
       * @param pixelPosition This argument is the pixel through which
       * the ray should pass.
       * 
       * @param normalize This argument indicates whether the ray
       * should be normalized to unit length before being returned.
       * 
       * @return The return value is the resulting ray.
       */
      virtual inline geometry::Ray3D<FloatType>
      reverseProject(const brick::numeric::Index2D& pixelPosition,
                     bool normalize = true) const;


      /** 
       * This function should be overridden by derived classes to
       * return a ray in 3D camera coordinates starting at the camera
       * focus and passing through the specified pixel position.
       * 
       * @param pixelPosition This argument is the point in pixel
       * coordinates through which the returned ray should pass.
       * 
       * @param normalize This argument indicates whether the ray
       * should be normalized to unit length before being returned.
       * 
       * @return The return value is the resulting ray.
       */
      virtual geometry::Ray3D<FloatType>
      reverseProject(const brick::numeric::Vector2D<FloatType>& pixelPosition,
                     bool normalize = true) const = 0;


    protected:

    };

  } // namespace computerVision
  
} // namespace brick


/* ============ Definitions of inline & template functions ============ */


#include <cmath>

namespace brick {

  namespace computerVision {

    // This function returns a ray in 3D camera coordinates starting
    // at the camera focus, and passing through the center of the
    // specified pixel.
    template <class FloatType>
    inline brick::geometry::Ray3D<FloatType>
    CameraIntrinsics<FloatType>::
    reverseProject(const brick::numeric::Index2D& pixelPosition,
                   bool normalize) const
    {
      return this->reverseProject(
        brick::numeric::Vector2D<FloatType>(pixelPosition.getColumn() + 0.5,
                                            pixelPosition.getRow() + 0.5),
        normalize);
    }

  } // namespace computerVision
  
} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_CAMERAINTRINSICS_HH */
