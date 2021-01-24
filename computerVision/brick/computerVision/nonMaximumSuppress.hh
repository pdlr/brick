/**
***************************************************************************
* @file brick/computerVision/nonMaximumSuppress.hh
*
* Header file declaring nonMaximumSuppress() function.
*
* Copyright (C) 2006,2012 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_NONMAXIMUMSUPPRESS_HH
#define BRICK_COMPUTERVISION_NONMAXIMUMSUPPRESS_HH

#include <brick/computerVision/image.hh>

namespace brick {

  namespace computerVision {

    /**
     * This function zeros any pixels of the input image that are not
     * plausible edges.
     */
    template <class FloatType, ImageFormat FORMAT>
    Image<FORMAT>
    nonMaximumSuppress(const Image<FORMAT>& inputImage,
                       const brick::numeric::Array2D<FloatType>& gradX,
                       const brick::numeric::Array2D<FloatType>& gradY);


  } // namespace computerVision

} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/nonMaximumSuppress_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_KERNEL_HH */
