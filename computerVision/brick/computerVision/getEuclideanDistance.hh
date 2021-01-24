/**
***************************************************************************
* @file brick/computerVision/getEuclideanDistance.hh
*
* Header file declaring the getEuclideanDistance() function template.
*
* Copyright (C) 2006,2012 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_EUCLIDEANDISTANCE_HH
#define BRICK_COMPUTERVISION_EUCLIDEANDISTANCE_HH

#include <list>
#include <brick/computerVision/imageFormat.hh>
#include <brick/computerVision/imageIO.hh>
#include <brick/computerVision/image.hh>

namespace brick {

  namespace computerVision {

    template<class FloatType, ImageFormat FORMAT>
    brick::numeric::Array2D<FloatType>
    getEuclideanDistance(const Image<FORMAT>& inputImage,
                         size_t maxNumberOfPasses=10);


    template<class FloatType, ImageFormat FORMAT>
    brick::numeric::Array2D<FloatType>
    getEuclideanDistance(const Image<FORMAT>& inputImage,
                         size_t maxNumberOfPasses,
                         size_t& numberOfPassesUsed);

  } // namespace computerVision

} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/getEuclideanDistance_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_EUCLIDEANDISTANCE_HH */
