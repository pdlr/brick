/**
***************************************************************************
* @file brick/computerVision/dilate.hh
*
* Header file declaring the dilate() function template.
*
* Copyright (C) 2006,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_DILATE_HH
#define BRICK_COMPUTERVISION_DILATE_HH

#include <brick/computerVision/image.hh>

namespace brick {

  namespace computerVision {
    
    template<ImageFormat FORMAT>
    Image<FORMAT>
    dilate(const Image<FORMAT>& inputImage);

    
    template<ImageFormat FORMAT>
    Image<FORMAT>
    dilateUsingBoxIntegrator(const Image<FORMAT>& inputImage,
                             unsigned int windowWidth = 3,
                             unsigned int windowHeight = 3);
    
  } // namespace computerVision
    
} // namespace brick

// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/dilate_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_DILATE_HH */
