/**
***************************************************************************
* @file brick/computerVision/connectedComponents.hh
*
* Header file declaring the connectedComponents() function template.
*
* Copyright (C) 2006,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_CONNECTEDCOMPONENTS_HH
#define BRICK_COMPUTERVISION_CONNECTEDCOMPONENTS_HH

#include <list>
#include <brick/computerVision/imageFormat.hh>
#include <brick/computerVision/image.hh>

namespace brick {

  namespace computerVision {
  
    /**
     * This function does connected components analysis on a previously
     * segmented image.
     * 
     * @param inputImage This argument is the segmented image.  Pixels
     * which are not part of a blob must have value of 0.  All other
     * values are considered to be blob pixels.  All non-zero pixels are
     * considered to be part of the same class.  That is, adjacent
     * pixels with different non-zero values will be considered to be
     * part of the same blob.
     * 
     * @return The return value is an image of labels in which each
     * pixel describes the the corresponding pixel of the input image.
     * Blobs in the input image will be labeled 0, 1, 2, etc. in the
     * output image.  Note that the assignment of labels to blobs is
     * unspecified: it is _not_ true that the largest blob gets the
     * lowest label.
     */
    template<ImageFormat FORMAT_OUT, ImageFormat FORMAT_IN>
    Image<FORMAT_OUT>
    connectedComponents(const Image<FORMAT_IN>& inputImage);
  

    /**
     * This function is just like connectedComponents(const Image&),
     * except that it also returns (by reference) the number of
     * components in the image.
     * 
     * @param inputImage This argument is the segmented image.  Pixels
     * which are not part of a blob must have value of 0.  All other
     * values are considered to be blob pixels.  All non-zero pixels are
     * considered to be part of the same class.  That is, adjacent
     * pixels with different non-zero values will be considered to be
     * part of the same blob.
     *
     * @param numberOfComponents This argument returns by reference
     * how many distinct components were identified in the image, not
     * counting the background.
     * 
     * @return The return value is an image of labels in which each
     * pixel describes the the corresponding pixel of the input image.
     * Blobs in the input image will be labeled 0, 1, 2, etc. in the
     * output image.  Note that the assignment of labels to blobs is
     * unspecified: it is _not_ true that the largest blob gets the
     * lowest label.
     */
    template<ImageFormat FORMAT_OUT, ImageFormat FORMAT_IN>
    Image<FORMAT_OUT>
    connectedComponents(const Image<FORMAT_IN>& inputImage,
                        unsigned int& numberOfComponents);

  } // namespace computerVision
    
} // namespace brick

// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/connectedComponents_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_CONNECTEDCOMPONENTS_HH */
