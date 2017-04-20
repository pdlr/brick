/**
***************************************************************************
* @file brick/computerVision/stereoRectify.hh
*
* Header file declaring code to compute image rectification parameters
* for stereo image pairs.
*
* Copyright (C) 2009,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_STEREORECTIFY_HH
#define BRICK_COMPUTERVISION_STEREORECTIFY_HH

#include <brick/computerVision/cameraIntrinsicsPinhole.hh>
#include <brick/numeric/transform2D.hh>
#include <brick/numeric/transform3D.hh>

namespace brick {

  namespace computerVision {

#if 0  /* Still in progress. */

    /** 
     * Given rectified intrinsics and extrinsics for a stereo camera,
     * calculate the reprojection matrix Q that transforms image
     * coordinates and disparities into 3D coordinates. 
     *
     * @param intrinsics0 This argument specifies the camera
     * intrinsics corresponding to rectified images from the first
     * (usually left) camera.
     * 
     * @param intrinsics1 This argument specifies the camera
     * intrinsics corresponding to rectified images from the second
     * (usually right) camera.  The value of this argument is normally
     * identical to intrinsics0, and is permitted to differ from
     * intrinsics0 only in the value of the U coordinate projection
     * center.
     * 
     * @param camera0Tworld This argument specifies the coordinate
     * transformation taking world coordinates and returning the
     * corresponding coordinates in the camera coordinate system
     * corresponding to intrinsics0.
     * 
     * @param camera1Tworld This argument specifies the coordinate
     * transformation taking world coordinates and returning the
     * corresponding coordinates in the camera coordinate system
     * corresponding to intrinsics1.  The rotation component of this
     * tranform must be identical to the rotation component of
     * camera0Tworld, although the translation components may differ.
     * 
     * @return The return value is the calculated reprojection matrix.
     */
    template <class FloatType>
    brick::numeric::Transform3D<FloatType>
    calculateReprojectionMatrix(
      CameraIntrinsicsPinhole<FloatType> const& intrinsics0,
      CameraIntrinsicsPinhole<FloatType> const& intrinsics1,
      numeric::Transform3D<FloatType> const& camera0Tworld,
      numeric::Transform3D<FloatType> const& camera1Tworld,
      FloatType epsilon);


    /** 
     * Given rectified intrinsics and a stereo baseline for a camera pair,
     * calculate the reprojection matrix Q that transforms image
     * coordinates and disparities into 3D coordinates.  That is:
     *
     * @code
     *    |a*x|       |u|
     *    |a*y| = Q * |v|
     *    |a*z|       |d|
     *    | a |       |1|
     * @endcode
     *
     * where 'a' is a projective scale factor, and 'd' is along-row
     * disparity between the two images.  This routine assumes that
     * the cameras are displaced in a direction parallel to the image
     * rows, as is normally the case for conventional stereo pairs.
     *
     * Note that it is permissible to have non-square image pixels.
     * 
     * @param intrinsics0 This argument specifies the camera
     * intrinsics corresponding to rectified images from the first
     * (usually left) camera.
     * 
     * @param intrinsics1 This argument specifies the camera
     * intrinsics corresponding to rectified images from the second
     * (usually right) camera.  The value of this argument is normally
     * identical to intrinsics0, and is permitted to differ from
     * intrinsics0 only in the value of the U coordinate projection
     * center.
     * 
     * @param This argument specifies the offset between left and
     * right cameras, which is assumed to be exactly in the direction
     * of the camera X axis.
     * 
     * @return The return value is the calculated reprojection matrix.
     */
    template <class FloatType>
    brick::numeric::Transform3D<FloatType>
    calculateReprojectionMatrix(
      CameraIntrinsicsPinhole<FloatType> const& intrinsics0,
      CameraIntrinsicsPinhole<FloatType> const& intrinsics1,
      FloatType baseline);
    
#endif /* #if 0 */

    
    /**
     * This function implements the stereo rectification algorithm of
     * Fusello, Trucco, and Verri [1].  Our convention is that, from
     * the perspective of rectified camera0, rectified camera1 is
     * located to the right (in the positive X direction).
     *
     * Warning: this routine assumes, but does not check, that the
     * camera?Tworld arguments are euclidean rigid body
     * transformations.  That is, the routine assumes that the 3x3
     * upper left block of each camera?Tworld argument is orthonormal,
     * and the last row is [0, 0, 0, 1].
     *
     * [1] A. Fusiello, E. Trucco, and A. Verri. A compact algorithm for
     * rectification of stereo pairs. Machine Vision and Applications,
     * 2000.
     * 
     * @param intrinsics0 This argument represents the intrinsic
     * calibration parameters of the left camera of the stereo pair.
     * 
     * @param intrinsics1 This argument represents the intrinsic
     * calibration parameters of the right camera of the stereo pair.
     * 
     * @param camera0Tworld This argument represents the extrinsic
     * parameters of the left camera of the stereo pair.  It specifies
     * a coordinate transform that takes coordinates in the world
     * coordinate system and returns camera0 coordinates representing
     * the same point.
     * 
     * @param camera1Tworld This argument represents the extrinsic
     * parameters of the right camera of the stereo pair.  It specifies
     * a coordinate transform that takes coordinates in the world
     * coordinate system and returns camera1 coordinates representing
     * the same point.
     * 
     * @param rectifiedIntrinsics0 This argument returns the updated
     * (rectified) intrinsics of the left camera.
     * 
     * @param rectifiedIntrinsics1 This argument returns the updated
     * (rectified) intrinsics of the right camera.
     * 
     * @param rcamera0Tworld This argument returns the rectified
     * extrinsic parameters of the left camera.  That is, it specifies
     * a coordinate transform that takes coordinates in the world
     * coordinate system and returns coordinates in the coordinate
     * system of the rectified left camera that represent the same
     * point.
     * 
     * @param rcamera1Tworld This argument returns the rectified
     * extrinsic parameters of the right camera.  That is, it
     * specifies a coordinate transform that takes coordinates in the
     * world coordinate system and returns coordinates in the
     * coordinate system of the rectified right camera that represent
     * the same point.
     * 
     * @param image0Trimage0 This argument returns a homography that
     * takes pixel coordinates in the unrectified left image, and
     * returns pixel coordinates in rectified left image that refer to
     * the same point in the image.
     * 
     * @param image1Trimage1 This argument returns a homography that
     * takes pixel coordinates in the unrectified right image, and
     * returns pixel coordinates in rectified right image that refer to
     * the same point in the image.
     */
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
                  numeric::Transform2D<FloatType>& image1Trimage1);


    /**
     * This function is identical to the 10-argument version of
     * stereoRectify(), except that it allows the user to explicitly
     * set the focal length of the rectified camera.  Note that the
     * focal length is specified in physical units, rather than in
     * pixels.  If you want to specify the focal length in units of
     * pixel size, you can do so by ensuring that input arguments
     * intrinsics0 and intrinsics1 have their pixel sizes set to 1.0.
     *
     * Warning: this routine assumes, but does not check, that the
     * camera?Tworld arguments are euclidean rigid body
     * transformations.  That is, the routine assumes that the 3x3
     * upper left block of each camera?Tworld argument is orthonormal,
     * and the last row is [0, 0, 0, 1].
     *
     * @param intrinsics0 This argument represents the intrinsic
     * calibration parameters of the left camera of the stereo pair.
     * 
     * @param intrinsics1 This argument represents the intrinsic
     * calibration parameters of the right camera of the stereo pair.
     * 
     * @param camera0Tworld This argument represents the extrinsic
     * parameters of the left camera of the stereo pair.  It specifies
     * a coordinate transform that takes coordinates in the world
     * coordinate system and returns camera0 coordinates representing
     * the same point.
     * 
     * @param camera1Tworld This argument represents the extrinsic
     * parameters of the right camera of the stereo pair.  It specifies
     * a coordinate transform that takes coordinates in the world
     * coordinate system and returns camera1 coordinates representing
     * the same point.
     *
     * @param rectifiedFocalLength The focal length, in physical
     * units, that is desired for the (virtual) rectified cameras that
     * will be returned.  If you don't know this, you can simply call
     * the 10-argument version of stereoRectify(), which sets the
     * returned focal length to the average of the focal lengths of
     * intrinsics0 and intrinsics1.
     * 
     * @param rectifiedIntrinsics0 This argument returns the updated
     * (rectified) intrinsics of the left camera.
     * 
     * @param rectifiedIntrinsics1 This argument returns the updated
     * (rectified) intrinsics of the right camera.
     * 
     * @param rcamera0Tworld This argument returns the rectified
     * extrinsic parameters of the left camera.  That is, it specifies
     * a coordinate transform that takes coordinates in the world
     * coordinate system and returns coordinates in the coordinate
     * system of the rectified left camera that represent the same
     * point.
     * 
     * @param rcamera1Tworld This argument returns the rectified
     * extrinsic parameters of the right camera.  That is, it
     * specifies a coordinate transform that takes coordinates in the
     * world coordinate system and returns coordinates in the
     * coordinate system of the rectified right camera that represent
     * the same point.
     * 
     * @param image0Trimage0 This argument returns a homography that
     * takes pixel coordinates in the unrectified left image, and
     * returns pixel coordinates in rectified left image that refer to
     * the same point in the image.
     * 
     * @param image1Trimage1 This argument returns a homography that
     * takes pixel coordinates in the unrectified right image, and
     * returns pixel coordinates in rectified right image that refer to
     * the same point in the image.
     */
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
                  numeric::Transform2D<FloatType>& image1Trimage1);
    

  } // namespace computerVision

} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/stereoRectify_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_STEREORECTIFY_HH */
