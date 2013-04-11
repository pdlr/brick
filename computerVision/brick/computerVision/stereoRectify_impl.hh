/**
***************************************************************************
* @file brick/computerVision/stereoRectify_impl.hh
*
* Header file defining inline and template functions declared in
* stereoRectify.hh.
*
* Copyright (C) 2009,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_STEREORECTIFY_IMPL_HH
#define BRICK_COMPUTERVISION_STEREORECTIFY_IMPL_HH

// This file is included by stereoRectify.hh, and should not be
// directly included by user code, so no need to include
// stereoRectify.hh here.
// 
// #include <brick/computerVision/stereoRectify.hh>

#include <brick/linearAlgebra/linearAlgebra.hh>
#include <brick/numeric/utilities.hh>


namespace brick {

  namespace computerVision {

    namespace privateCode {

      template <class FloatType>
      brick::numeric::Array2D<FloatType>
      extractUpperLeftBlock(
        CameraIntrinsicsPinhole<FloatType> const& intrinsics,
        brick::numeric::Transform3D<FloatType> const& cameraTworld)
      {
        brick::numeric::Array2D<FloatType> result(3, 3);
#if 0
        result(0, 0) = (intrinsics.getKx() * cameraTworld(0, 0)
                        + intrinsics.getCenterU() * cameraTworld(2, 0));
        result(0, 1) = (intrinsics.getKx() * cameraTworld(0, 1)
                        + intrinsics.getCenterU() * cameraTworld(2, 1));
        result(0, 2) = (intrinsics.getKx() * cameraTworld(0, 2)
                        + intrinsics.getCenterU() * cameraTworld(2, 2));
        result(1, 0) = (intrinsics.getKy() * cameraTworld(1, 0)
                        + intrinsics.getCenterV() * cameraTworld(2, 0));
        result(1, 1) = (intrinsics.getKy() * cameraTworld(1, 1)
                        + intrinsics.getCenterV() * cameraTworld(2, 1));
        result(1, 2) = (intrinsics.getKy() * cameraTworld(1, 2)
                        + intrinsics.getCenterV() * cameraTworld(2, 2));
        result(2, 0) = cameraTworld(2, 0);
        result(2, 1) = cameraTworld(2, 1);
        result(2, 2) = cameraTworld(2, 2);
#else
        result(0, 0) = (intrinsics.getKx() * cameraTworld(0, 0)
                        + intrinsics.getCenterU() * cameraTworld(2, 0));
        result(0, 1) = (intrinsics.getKx() * cameraTworld(0, 1)
                        + intrinsics.getCenterU() * cameraTworld(2, 1));
        result(0, 2) = (intrinsics.getKx() * cameraTworld(0, 2)
                        + intrinsics.getCenterU() * cameraTworld(2, 2));

        FloatType ky = intrinsics.getKy();
        FloatType c10 =  cameraTworld(1, 0);
        FloatType dummy = ky * c10;
        dummy += intrinsics.getCenterV() * cameraTworld(2, 0);
        result(1, 0) = dummy;
        result(1, 1) = (intrinsics.getKy() * cameraTworld(1, 1)
                        + intrinsics.getCenterV() * cameraTworld(2, 1));
        result(1, 2) = (intrinsics.getKy() * cameraTworld(1, 2)
                        + intrinsics.getCenterV() * cameraTworld(2, 2));
        result(2, 0) = cameraTworld(2, 0);
        result(2, 1) = cameraTworld(2, 1);
        result(2, 2) = cameraTworld(2, 2);
#endif
        return result;
      }

  
    }; // namespace privateCode
    
    template <class FloatType>
    void
    stereoRectify(CameraIntrinsicsPinhole<FloatType> const& intrinsics0,
                  CameraIntrinsicsPinhole<FloatType> const& intrinsics1,
                  numeric::Transform3D<FloatType> const& camera0Tworld,
                  numeric::Transform3D<FloatType> const& camera1Tworld,
                  CameraIntrinsicsPinhole<FloatType>& rectifiedIntrinsics0,
                  CameraIntrinsicsPinhole<FloatType>& rectifiedIntrinsics1,
                  numeric::Transform3D<FloatType>& rcamera0Tworld,
                  numeric::Transform3D<FloatType>& rcamera1Tworld,
                  numeric::Transform2D<FloatType>& image0Trimage0,
                  numeric::Transform2D<FloatType>& image1Trimage1)
    {
      // Hack(xxx): find a principled way to set this threshold.
      FloatType localEpsilon = 1.0E-12;
      
      // Recover optical centers (focal points) of the two cameras in
      // world coords.  These won't change with rectification.
      brick::numeric::Transform3D<FloatType> worldTcamera0 =
        camera0Tworld.invert();
      brick::numeric::Transform3D<FloatType> worldTcamera1 =
        camera1Tworld.invert();
      brick::numeric::Vector3D<FloatType> focalPoint0(
        worldTcamera0(0, 3), worldTcamera0(1, 3), worldTcamera0(2, 3));
      brick::numeric::Vector3D<FloatType> focalPoint1(
        worldTcamera1(0, 3), worldTcamera1(1, 3), worldTcamera1(2, 3));

      // Recover new X axis (expressed in world coordinates) for both
      // cameras.  This is parallel to the stereo baseline.
      brick::numeric::Vector3D<FloatType> xAxis = focalPoint1 - focalPoint0;
      FloatType xMagnitude = brick::numeric::magnitude<FloatType>(xAxis);
      if(xMagnitude < localEpsilon) {
        BRICK_THROW(brick::common::ValueException, "rectify()",
                    "Optical centers of camera1 and camera0 are too close "
                    "together.");
      }
      xAxis /= xMagnitude;

      // New Y axis (expressed in world coordinates), again for both
      // cameras, is orthogonal to X, and orthogonal to Z axis of
      // unrectified camera0.  Note that this fails if camera1 lies on
      // the optical axis of camera0.
      brick::numeric::Vector3D<FloatType> oldZ(
        camera0Tworld(2, 0), camera0Tworld(2, 1), camera0Tworld(2, 2));
      brick::numeric::Vector3D<FloatType> yAxis =
        brick::numeric::cross<FloatType>(oldZ, xAxis);
      FloatType yMagnitude = brick::numeric::magnitude<FloatType>(yAxis);
      if(yMagnitude < localEpsilon) {
        BRICK_THROW(brick::common::ValueException, "rectify()",
                    "Camera1 lies too close to the optical axis of camera0.");
      }
      yAxis /= yMagnitude;

      // New Z axis (expressed in world coordinates) follows directly
      // from X and Y.
      brick::numeric::Vector3D<FloatType> zAxis =
        brick::numeric::cross<FloatType>(xAxis, yAxis);
      FloatType zMagnitude = brick::numeric::magnitude<FloatType>(zAxis);
      if(zMagnitude < localEpsilon) {
        BRICK_THROW(brick::common::LogicException, "rectify()",
                    "Error in computation of z axis.");
      }
      zAxis /= zMagnitude;
      
      // The rectified cameras share the same orientation.  The
      // locations of the optical centers (i.e., the translation
      // component of the cameraTworld transforms) remain unchanged.
      brick::numeric::Transform3D<FloatType> worldTrcamera0(
        xAxis.x(), yAxis.x(), zAxis.x(), worldTcamera0(0, 3),
        xAxis.y(), yAxis.y(), zAxis.y(), worldTcamera0(1, 3),
        xAxis.z(), yAxis.z(), zAxis.z(), worldTcamera0(2, 3),
        0.0, 0.0, 0.0, 1.0);
      brick::numeric::Transform3D<FloatType> worldTrcamera1(
        xAxis.x(), yAxis.x(), zAxis.x(), worldTcamera1(0, 3),
        xAxis.y(), yAxis.y(), zAxis.y(), worldTcamera1(1, 3),
        xAxis.z(), yAxis.z(), zAxis.z(), worldTcamera1(2, 3),
        0.0, 0.0, 0.0, 1.0);

      rcamera0Tworld = worldTrcamera0.invert();
      rcamera1Tworld = worldTrcamera1.invert();
      
      // The rectified cameras will share the same intrinsic
      // parameters.  Note that these parameters are not uniquely
      // determined.  Choosing them poorly just means that the input
      // images will project into the rectified images suboptimally.
      //
      // Note(xxx): eventually we'll want to do something smart here,
      // but for now we just take the easy path of averaging values
      // from the two input cameras, and taking image size from
      // camera0.
      unsigned int numPixelsX = intrinsics0.getNumPixelsX();
      unsigned int numPixelsY = intrinsics0.getNumPixelsY();
      FloatType focalLength =
        (intrinsics0.getFocalLength() + intrinsics1.getFocalLength()) / 2.0;
      FloatType pixelSizeX = 
        (intrinsics0.getPixelSizeX() + intrinsics1.getPixelSizeX()) / 2.0;
      FloatType pixelSizeY = 
        (intrinsics0.getPixelSizeY() + intrinsics1.getPixelSizeY()) / 2.0;
      FloatType centerU = 
        (intrinsics0.getCenterU() + intrinsics1.getCenterU()) / 2.0;
      FloatType centerV = 
        (intrinsics0.getCenterV() + intrinsics1.getCenterV()) / 2.0;
      rectifiedIntrinsics0 = CameraIntrinsicsPinhole<FloatType>(
        numPixelsX, numPixelsY, focalLength, pixelSizeX, pixelSizeY,
        centerU, centerV);
      rectifiedIntrinsics1 = rectifiedIntrinsics0;

      // Finally, 2D transformations to take coordinates in the
      // rectified images and return matching coordinates in the
      // unrectified image.
      brick::numeric::Array2D<FloatType> oldQ0 =
        privateCode::extractUpperLeftBlock(intrinsics0, camera0Tworld);
      brick::numeric::Array2D<FloatType> oldQ1 =
        privateCode::extractUpperLeftBlock(intrinsics1, camera1Tworld);
      brick::numeric::Array2D<FloatType> newQ0 =
        privateCode::extractUpperLeftBlock(
          rectifiedIntrinsics0, rcamera0Tworld);
      brick::numeric::Array2D<FloatType> newQ1 =
        privateCode::extractUpperLeftBlock(
          rectifiedIntrinsics1, rcamera1Tworld);

      brick::numeric::Array2D<FloatType> rectArray0 =
        brick::numeric::matrixMultiply<FloatType>(
          oldQ0, brick::linearAlgebra::inverse(newQ0));
      brick::numeric::Array2D<FloatType> rectArray1 =
        brick::numeric::matrixMultiply<FloatType>(
          oldQ1, brick::linearAlgebra::inverse(newQ1));
      image0Trimage0.setTransform(
        rectArray0(0, 0), rectArray0(0, 1), rectArray0(0, 2),
        rectArray0(1, 0), rectArray0(1, 1), rectArray0(1, 2),
        rectArray0(2, 0), rectArray0(2, 1), rectArray0(2, 2));
      image1Trimage1.setTransform(
        rectArray1(0, 0), rectArray1(0, 1), rectArray1(0, 2),
        rectArray1(1, 0), rectArray1(1, 1), rectArray1(1, 2),
        rectArray1(2, 0), rectArray1(2, 1), rectArray1(2, 2));
    }


  } // namespace computerVision

} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_STEREORECTIFY_IMPL_HH */
