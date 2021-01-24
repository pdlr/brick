/**
***************************************************************************
* @file brick/computerVision/sobel.hh
*
* Header file declaring routines to compute sobel edge images.
*
* Copyright (C) 2006,2012 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_SOBEL_HH
#define BRICK_COMPUTERVISION_SOBEL_HH

#include <brick/computerVision/image.hh>

namespace brick {

  namespace computerVision {


    /**
     * This function applies the sobel edge operator in the X
     * direction.  Specifically, convolves the image with the 3x3 kernel
     *
     * @code
     *   [[-1, 0, 1],
     *    [-2, 0, 2],
     *    [-1, 0, 1]]
     * @endcode
     *
     * and then optionally rescales the result to be proportional to
     * the image gradient with scale factor 1.0.
     *
     * @param inputImage This argument is the image to be convolved.
     *
     * @param normalizeResult This argument specifies whether or not
     * to divide the resulting pixel values by 8.  Setting this
     * argument to true will result in lost precision on integer
     * types.  This feature is currently not implemented, so please
     * leave normalizeResult at its default value of false.
     *
     * @return The return value is the result of the convolution.
     */
    template <ImageFormat FORMAT>
    Image<FORMAT>
    applySobelX(const Image<FORMAT>& inputImage, bool normalizeResult=false);


    /**
     * This function applies the sobel edge operator in the Y
     * direction.  Specifically, convolves the image with the 3x3 kernel
     *
     * @code
     *   [[-1, -2, -1],
     *    [ 0,  0,  0],
     *    [ 1,  2,  1]]
     * @endcode
     *
     * and then optionally rescales the result to be proportional to
     * the image gradient with scale factor 1.0.
     *
     * @param inputImage This argument is the image to be convolved.
     *
     * @param normalizeResult This argument specifies whether or not
     * to divide the resulting pixel values by 8.  Setting this
     * argument to true will result in lost precision on integer
     * types.  This feature is currently not implemented, so please
     * leave normalizeResult at its default value of false.
     *
     * @return The return value is the result of the convolution.
     */
    template <ImageFormat FORMAT>
    Image<FORMAT>
    applySobelY(const Image<FORMAT>& inputImage, bool normalizeResult=false);

  } // namespace computerVision

} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/sobel_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_KERNEL_HH */
