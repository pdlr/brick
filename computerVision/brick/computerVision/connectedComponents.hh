/**
***************************************************************************
* @file brick/computerVision/connectedComponents.hh
*
* Header file declaring the connectedComponents() function template.
*
* Copyright (C) 2006,2012 David LaRose, dlr@davidlarose.com
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
     ** This struct controls the operation of the
     ** connectedComponents() function template.
     **/
    struct ConnectedComponentsConfig {
      enum Mode {
        /// The input image is considered to be binary, with only
        /// nonzero pixels being grouped into components.  Zero pixels
        /// are considered to be background and are ignored.  Nonzero
        /// pixels of any value can be grouped into a component
        /// together.
        FOREGROUND_BACKGROUND,

        /// Pixels that have the same value are grouped together in
        /// connected components.  Zero pixels can be grouped into
        /// components (with other zero pixels).  All pixels will be
        /// assigned to one component or another, i.e., there is no
        /// explicit background.
        SAME_COLOR
      };

      Mode mode = FOREGROUND_BACKGROUND;
    };


    /**
     ** This functor simply returns true when two pixels have the same
     ** value.  By default, this is what connectedComponents uses to
     ** decide if pixels have the same value.  Users can supply their
     ** own functor (below) in order to change how pixels are
     ** compared.
     **/
    template<class PixelType>
    struct PixelEqualityComparator {
      bool
      operator()(PixelType const& value0, PixelType const& value1) const {
        return value0 == value1;
      }
    };


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
    template<ImageFormat FORMAT_OUT, ImageFormat FORMAT_IN,
             class Comparator = PixelEqualityComparator<
               typename ImageFormatTraits<FORMAT_IN>::PixelType>>
    Image<FORMAT_OUT>
    connectedComponents(const Image<FORMAT_IN>& inputImage,
                        ConnectedComponentsConfig const& config
                        = ConnectedComponentsConfig(),
                        Comparator comparator = Comparator());


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
    template<ImageFormat FORMAT_OUT, ImageFormat FORMAT_IN,
             class Comparator = PixelEqualityComparator<
               typename ImageFormatTraits<FORMAT_IN>::PixelType>>
    Image<FORMAT_OUT>
    connectedComponents(const Image<FORMAT_IN>& inputImage,
                        unsigned int& numberOfComponents,
                        ConnectedComponentsConfig const& config
                        = ConnectedComponentsConfig(),
                        Comparator comparator = Comparator());

  } // namespace computerVision

} // namespace brick

// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/connectedComponents_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_CONNECTEDCOMPONENTS_HH */
