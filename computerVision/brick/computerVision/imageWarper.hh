/**
***************************************************************************
* @file brick/computerVision/imageWarper.hh
*
* Header file declaring ImageWarper class.
*
* Copyright (C) 2009,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_IMAGEWARPER_HH
#define BRICK_COMPUTERVISION_IMAGEWARPER_HH

#include <brick/computerVision/image.hh>
#include <brick/numeric/array2D.hh>

namespace brick {

  namespace computerVision {
    
    /**
     ** This class does lookup-table-based warping of images.  Image
     ** resampling is done using bilinear interpolation.
     ** 
     ** Template argument NumericType specifies the type (usually
     ** float or double) used to do bilinear interpolation.
     **
     ** Template argument TransformFunctor specifies the type of a
     ** functor that takes 2D pixel coordinates in the output image
     ** and returns the 2D pixel coordinates in the input image from
     ** which the output pixel should take it's color.  It must take a
     ** single Vector2D<NumericType> instance as its argument, and
     ** return a Vector2D<NumericType> instance.
     **/
    template <class NumericType, class TransformFunctor>
    class ImageWarper
    {
    public:

      /* ******** Public member functions ******** */

      /** 
       * Default constructor makes a non-functioning ImageWarper instance.
       */
      ImageWarper();

    
      /** 
       * This constructor fills in the ImageWarper state so that it's
       * ready to transform images.  Because it builds the lookup
       * table during construction, this constructor is likely to be
       * slow.
       * 
       * @param inputRows This argument specifies the height, in
       * pixels, of the images that will be warped by the ImageWarper.
       * 
       * @param inputColumns This argument specifies the width, in
       * pixels, of the images that will be warped by the ImageWarper.
       * 
       * @param outputRows This argument specifies the height, in
       * pixels, of the result images (after warping).
       * 
       * @param outputColumns This argument specifies the width, in
       * pixels, of the result images (after warping).
       * 
       * @param transformer This argument is a functor that defines
       * the warp.  Please see the documentation for class ImageWarper
       * for more information.
       */
      ImageWarper(size_t inputRows, size_t inputColumns,
                  size_t outputRows, size_t outputColumns,
                  TransformFunctor transformer);

    
      /**
       * Destroys the ImageWarper instance and deletes the internal data
       * store.
       */
      virtual
      ~ImageWarper();


      /** 
       * Warps a single image using the pre-computed lookup table.
       * 
       * @param inputImage This argument is the image to be warped.
       * 
       * @param defaultValue This argument specifies what pixel value
       * to use for pixels in the output image that map to input-image
       * pixels that lie outside the boundaries of the input image.
       * 
       * @return The return value is the warped output image.
       */
      template <ImageFormat InputFormat, ImageFormat OutputFormat>
      Image<OutputFormat>
      warpImage(Image<InputFormat> const& inputImage,
                typename Image<OutputFormat>::PixelType defaultValue) const;

    private:

      struct SampleInfo {
        NumericType c00;
        NumericType c01;
        NumericType c10;
        NumericType c11;
        size_t index00;
        bool isInBounds;
      };

      size_t m_inputColumns;
      size_t m_inputRows;
      brick::numeric::Array2D<SampleInfo> m_lookupTable;

    };

  } // namespace computerVision

} // namespace brick

// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/imageWarper_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_IMAGEWARPER_HH */
