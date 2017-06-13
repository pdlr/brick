/**
***************************************************************************
* @file brick/computerVision/stereoRectify_impl.hh
*
* Header file defining inline and template functions declared in
* stereoRectify.hh.
*
* Copyright (C) 2009-2017 David LaRose, dlr@cs.cmu.edu
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
#include <brick/numeric/mathFunctions.hh>
#include <brick/numeric/numericTraits.hh>
#include <brick/numeric/utilities.hh>


namespace brick {

  namespace computerVision {

    namespace privateCode {

      template <class FloatType>
      bool
      checkEqual(FloatType const& x0, FloatType const& x1)
      {
        return (brick::numeric::absoluteValue(x1 - x0)
                <= brick::numeric::NumericTraits<FloatType>::epsilon());
      }

      
      template <class FloatType>
      brick::numeric::Array2D<FloatType>
      extractUpperLeftBlock(
        CameraIntrinsicsPinhole<FloatType> const& intrinsics,
        brick::numeric::Transform3D<FloatType> const& cameraTworld)
      {
        brick::numeric::Array2D<FloatType> result(3, 3);
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

        return result;
      }

    } // namespace privateCode


    // Given rectified intrinsics and a stereo baseline for a camera pair,
    // calculate the reprojection matrix Q that transforms image
    // coordinates and disparities into 3D coordinates.
    template <class FloatType>
    brick::numeric::Transform3D<FloatType>
    getReprojectionMatrix(
      CameraIntrinsicsPinhole<FloatType> const& intrinsics0,
      CameraIntrinsicsPinhole<FloatType> const& intrinsics1,
      FloatType baseline)
    {
      // First check arguments.
      if(!(privateCode::checkEqual(intrinsics0.getFocalLengthX(),
                                   intrinsics1.getFocalLengthX())
           && privateCode::checkEqual(intrinsics0.getFocalLengthY(),
                                      intrinsics1.getFocalLengthY())
           && privateCode::checkEqual(intrinsics0.getCenterV(),
                                      intrinsics1.getCenterV()))) {
        BRICK_THROW(brick::common::ValueException,
                    "calculateReprojectionMatrix()",
                    "Intrinsics instances may differ only in the value "
                    "of the projection center U coordinate.");
      }

      // In the comments below, we assume intrinsics0 refers to
      // the left camera of a stereo pair, and intrinsics1 refers
      // to the right camera.
      //
      // Projection equations for the left camera are:
      //
      // @verbatim
      //       |u|   |f_x,   0, c_u, 0|   |x|
      //   a * |v| = |0,   f_y, c_v, 0| * |y|
      //       |1|   |0,     0,   1, 0|   |z|
      //                                  |1|
      // @endverbatim
      //
      // where a is an arbitrary scale factor.
      //
      // Equations for the right camera are very similar.
      // 
      // @verbatim
      //       |u - d|   |f_x,   0, c_u', 0|   |x - b|
      //   a * |v    | = |0,   f_y,  c_v, 0| * |y    |
      //       |1    |   |0,     0,    1, 0|   |z    |
      //                                       |1    |
      // @endverbatim
      //
      // where d is the stereo disparity (the u coordinate in the
      // left camera vs the u coordinate in right camera), b is
      // the stereo baseline, and c_u' is the only right camera
      // projection parameter permitted to differ from the left
      // camera.
      //
      // The right camera projection can be rearranged to move b
      // into the projection matrix.
      // 
      // @verbatim
      //       |u - d|   |f_x,   0, c_u', -b*f_x|   |x|
      //   a * |v    | = |0,   f_y,  c_v,      0| * |y|
      //       |1    |   |0,     0,    1,      0|   |z|
      //                                            |1|
      // @endverbatim
      // 
      // Only the top row of the right camera equations differs from
      // the left camera equations.  Writing the four unique
      // simultaneous equations together gives:
      //
      // @verbatim
      //       |u    |   |f_x,   0,  c_u,      0|   |x|
      //   a * |v    | = |0,   f_y,  c_v,      0| * |y|
      //       |u - d|   |f_x,   0, c_u', -b*f_x|   |z|
      //       |1    |   |0,     0,    1,      0|   |1|
      // @endverbatim
      //
      // which is easily rearranged to:
      // 
      // @verbatim
      //       |u|   |f_x,   0,        c_u,     0|   |x|
      //   a * |v| = |0,   f_y,        c_v,     0| * |y|
      //       |d|   |0,     0, c_u - c_u', b*f_x|   |z|
      //       |1|   |0,     0,          1,     0|   |1|
      // @endverbatim
      // 
      // If we invert the 4x4 matrix in this last equation using the
      // cofactor method, and discard a scale factor of (-f_x /
      // determinant) because it simply changes the projective scale
      // factor (which we've called "a" above, but call "b" below), we get:
      //
      // @verbatim
      //       |x|   |f_y * b,       0,   0,     -f_y * c_u * b|   |u|
      //   b * |y| = |0,       f_x * b,   0,     -f_x * c_v * b| * |v|
      //       |z|   |0,             0,   0,      f_x * f_y * b|   |d|
      //       |1|   |0,             0, f_y, f_y * (c_u' - c_u)|   |1|
      // @endverbatim

      FloatType fX = intrinsics0.getFocalLengthX();
      FloatType fY = intrinsics0.getFocalLengthY();
      FloatType cU = intrinsics0.getCenterU();
      FloatType cV = intrinsics0.getCenterV();
      FloatType cUPrime = intrinsics1.getCenterU();
      FloatType fXTimesB = fX * baseline;
      FloatType fYTimesB = fY * baseline;
      return brick::numeric::Transform3D<FloatType>(
        fYTimesB, 0.0, 0.0, -fYTimesB * cU,
        0.0, fXTimesB, 0.0, -fXTimesB * cV,
        0.0, 0.0, 0.0, fXTimesB * fY,
        0.0, 0.0, fY, fY * (cUPrime - cU));
    }

    
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
      FloatType rectifiedFocalLength =
        (intrinsics0.getFocalLength() + intrinsics1.getFocalLength()) / 2.0;
      stereoRectify(intrinsics0, intrinsics1, camera0Tworld, camera1Tworld,
                    rectifiedFocalLength,
                    rectifiedIntrinsics0, rectifiedIntrinsics1,
                    rcamera0Tworld, rcamera1Tworld,
                    image0Trimage0, image1Trimage1);
    }
    
    
    template <class FloatType>
    void
    stereoRectify(CameraIntrinsicsPinhole<FloatType> const& intrinsics0,
                  CameraIntrinsicsPinhole<FloatType> const& intrinsics1,
                  numeric::Transform3D<FloatType> const& camera0Tworld,
                  numeric::Transform3D<FloatType> const& camera1Tworld,
                  FloatType const& rectifiedFocalLength,
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
      FloatType focalLength = rectifiedFocalLength;
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
