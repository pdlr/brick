/**
***************************************************************************
* @file brick/computerVision/histogramEqualize.hh
*
* Header file declaring histogram equalization routines.
*
* Copyright (C) 2005,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_HHISTOGRAMEQUALIZE_H
#define BRICK_COMPUTERVISION_HHISTOGRAMEQUALIZE_H

#include <brick/common/types.hh>
#include <brick/computerVision/image.hh>


namespace brick {

  namespace computerVision {
    
    /** 
     * This function computes the histogram of an image.  That is, it
     * counts the number of pixels with each possible value and returns
     * a 1D array of counts.
     * 
     * @param inputImage This argument is the image to be histogrammed.
     * 
     * @return The return value is a 1D array in which the first element
     * indicates the number of pixels having the value 0, the second
     * element indicates the number of pixels having the value 1, and so
     * forth.
     */
    brick::numeric::Array1D<brick::common::UInt32>
    getHistogram(const Image<GRAY8>& inputImage);
  

    /** 
     * This function remaps the pixel values of the input image in such
     * a way that output pixel value increases monotonically with input
     * pixel value, and the histogram of the output image is nearly
     * flat.
     * 
     * @param inputImage This argument is the image to be equalized.
     * 
     * @return The return value is the histogram equalized image.
     */
    Image<GRAY8>
    histogramEqualize(const Image<GRAY8>& inputImage);

  } // namespace computerVision

} // namespace brick


#endif /* #ifndef BRICK_COMPUTERVISION_HHISTOGRAMEQUALIZE_H */
