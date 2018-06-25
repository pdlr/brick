/**
***************************************************************************
* @file brick/computerVision/canny.hh
*
* Header file declaring a function template to compute canny edge images.
*
* Copyright (C) 2006,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_CANNY_HH
#define BRICK_COMPUTERVISION_CANNY_HH

#include <brick/computerVision/image.hh>

namespace brick {

  namespace computerVision {


    /**
     * This function applies the canny edge detector to the input image.
     *
     * @param inputImage This argument is the image to be edge-detected.
     *
     * @param gaussianSize This argument specifies the size, in
     * pixels, of the blurring filter to be applied to the image prior
     * to edge detection.  The sigma of the blurring filter will be
     * set to 1/6 of this size.  This argument must be zero or an odd
     * number. Setting this argument to zero disables the pre-filter.
     *
     * @param upperThreshold This argument specifies the gradient
     * magnitude necessary for a pixel to be considered a "seed" point
     * from which to grow a new edge.  If this argument is less than
     * or equal to 0.0, it will be calculated automatically using the
     * value of argument autoUpperThresholdFactor.
     *
     * @param lowerThreshold This argument specifies the gradient
     * magnitude necessary for a pixel to be considered as part of an
     * existing edge.  If this argument is less than or equal to 0.0,
     * it will be calculated automatically using the value of argument
     * autoUpperThresholdFactor.
     *
     * @param autoUpperThresholdFactor If upperThreshold is less than
     * or equal to zero, then this argument is used to set calculate
     * the threshold automatically.  Smaller values (and increasingly
     * large negative values) make it easier to start an edge.  If
     * argument upperThreshold is greater than 0.0, then this argument
     * is ignored.
     *
     * @param autoLowerThresholdFactor If lowerThreshold is less than
     * or equal to zero, then this argument is used to calculate the
     * threshold automatically.  Its value should be positive and less
     * than or equal to autoUpperThresholdFactor.  Larger values make
     * edges shorter.  Smaller values (and increasingly large negative
     * values) make edges tend to stretch out longer.  If argument
     * lowerThreshold is greater than 0.0, then this argument is
     * ignored.
     *
     * @return The return value is a binary image in which all edge
     * pixels are true, and all non-edge pixels are false.
     */
    template <class FloatType, ImageFormat FORMAT>
    Image<GRAY1>
    applyCanny(const Image<FORMAT>& inputImage,
               unsigned int gaussianSize = 5,
               FloatType upperThreshold = 0.0,
               FloatType lowerThreshold = 0.0,
               FloatType autoUpperThresholdFactor = 3.0,
               FloatType autoLowerThresholdFactor = 0.0);

    /**
     * This function applies the canny edge detector to the input image.
     *
     * @param inputImage This argument is the image to be edge-detected.
     *
     * @param gaussianSize This argument specifies the size, in
     * pixels, of the blurring filter to be applied to the image prior
     * to edge detection.  The sigma of the blurring filter will be
     * set to 1/6 of this size.  This argument must be zero or an odd
     * number. Setting this argument to zero disables the pre-filter.
     *
     * @param upperThreshold This argument specifies the gradient
     * magnitude necessary for a pixel to be considered a "seed" point
     * from which to grow a new edge.  If this argument is less than
     * or equal to 0.0, it will be calculated automatically using the
     * value of argument autoUpperThresholdFactor.
     *
     * @param lowerThreshold This argument specifies the gradient
     * magnitude necessary for a pixel to be considered as part of an
     * existing edge.  If this argument is less than or equal to 0.0,
     * it will be calculated automatically using the value of argument
     * autoUpperThresholdFactor.
     *
     * @param autoUpperThresholdFactor If upperThreshold is less than
     * or equal to zero, then this argument is used to set calculate
     * the threshold automatically.  Smaller values (and increasingly
     * large negative values) make it easier to start an edge.  If
     * argument upperThreshold is greater than 0.0, then this argument
     * is ignored.
     *
     * @param autoLowerThresholdFactor If lowerThreshold is less than
     * or equal to zero, then this argument is used to calculate the
     * threshold automatically.  Its value should be positive and less
     * than or equal to autoUpperThresholdFactor.  Larger values make
     * edges shorter.  Smaller values (and increasingly large negative
     * values) make edges tend to stretch out longer.  If argument
     * lowerThreshold is greater than 0.0, then this argument is
     * ignored.
     *
     * @return The return value is a binary image in which all edge
     * pixels are true, and all non-edge pixels are false.
     */
    template <class FloatType, ImageFormat FORMAT>
    Image<GRAY1>
    applyCanny(const Image<FORMAT>& inputImage,
               brick::numeric::Array2D<FloatType>& gradientX,
               brick::numeric::Array2D<FloatType>& gradientY,
               unsigned int gaussianSize = 5,
               FloatType upperThreshold = 0.0,
               FloatType lowerThreshold = 0.0,
               FloatType autoUpperThresholdFactor = 3.0,
               FloatType autoLowerThresholdFactor = 0.0);

  } // namespace computerVision

} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/canny_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_CANNY_HH */
