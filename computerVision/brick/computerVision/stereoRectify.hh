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
            

  } // namespace computerVision

} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/stereoRectify_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_STEREORECTIFY_HH */
