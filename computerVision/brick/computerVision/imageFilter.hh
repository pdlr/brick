/**
***************************************************************************
* @file brick/computerVision/imageFilter.hh
*
* Header file declaring functions for performing image filtering.
*
* Copyright (C) 2006, 2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_IMAGEFILTER_HH
#define BRICK_COMPUTERVISION_IMAGEFILTER_HH

#include <brick/computerVision/image.hh>
#include <brick/computerVision/kernel.hh>
#include <brick/numeric/convolutionStrategy.hh>

namespace brick {

  namespace computerVision {

    /** 
     * This enum is used to indicate how filtering routines should
     * handle the edges of images where the filter kernel does not
     * completely overlap the image.  Currently filter() only supports
     * BRICK_CONVOLVE_PAD_RESULT.
     */
    using brick::numeric::ConvolutionStrategy;
    using brick::numeric::BRICK_CONVOLVE_PAD_RESULT;
    
    /** 
     * This function filters an image with the given kernel.  Argument
     * "convolutionStrategy" indicates what to do with the edges of the
     * filtered image (where the filter kernel only partially overlaps
     * the input image).  Currently filter() only supports
     * BRICK_CONVOLVE_PAD_RESULT.
     * 
     * @param kernel This argument is the Kernel instance with which
     * to filter.
     * 
     * @param image This argument is the Image to be filtered.
     *
     * @param fillValue This argument specifies the value with which
     * image edges should be padded.
     * 
     * @param convolutionStrategy This argument specifies how to handle the
     * edges of the image.  Please see the dlrNumeric documentation
     * for more information.
     * 
     * @return The return value is a filtered copy of image.
     */
    template<ImageFormat OutputFormat,
	     ImageFormat IntermediateFormat,
             ImageFormat ImageFormat,
             class KernelType>
    Image<OutputFormat> filter2D(
      const Kernel<KernelType>& kernel,
      const Image<ImageFormat>& image,
      const typename ImageFormatTraits<OutputFormat>::PixelType fillValue
      = typename ImageFormatTraits<OutputFormat>::PixelType(),
      ConvolutionStrategy convolutionStrategy = BRICK_CONVOLVE_PAD_RESULT);
    

    /** 
     * This function filters an image with the given kernel, placing
     * the result into a pre-constructed Image instance.  Argument
     * "convolutionStrategy" indicates what to do with the edges of the
     * filtered image (where the filter kernel only partially overlaps
     * the input image).  Currently filter() only supports
     * BRICK_CONVOLVE_PAD_RESULT.
     * 
     * @param outputImage This argument is used to return the result.
     * The associated memory is not reallocated unless outputImage has
     * a different number of rows and/or columns than image.
     *
     * @param kernel This argument is the Kernel instance with which
     * to filter.
     * 
     * @param image This argument is the Image to be filtered.
     * 
     * @param fillValue This argument specifies the value with which
     * image edges should be padded.
     * 
     * @param convolutionStrategy This argument specifies how to handle the
     * edges of the image.  Please see the dlrNumeric documentation
     * for more information.
     */
    template<ImageFormat OutputFormat,
	     ImageFormat IntermediateFormat,
             ImageFormat ImageFormat,
             class KernelType>
    void
    filter2D(
      Image<OutputFormat>& outputImage,
      const Kernel<KernelType>& kernel,
      const Image<ImageFormat>& image,
      const typename ImageFormatTraits<OutputFormat>::PixelType fillValue,
      ConvolutionStrategy convolutionStrategy = BRICK_CONVOLVE_PAD_RESULT);


  } // namespace computerVision
  
} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/imageFilter_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_IMAGEFILTER_HH */
