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
     * Note that the convolution is performed internal to the function
     * using pixel type determined by OutputFormat.  Using a
     * low-precision output format, such as GRAY8, will almost
     * certainly result in pixel overflow.  Using an integer output
     * format will introduce quantification errors in the convolution.
     * If you're low-pass filtering an integer-valued image, consider
     * using filterBinomial2D() instead.
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
     * Note that the convolution is performed internal to the function
     * using pixel type determined by OutputFormat.  Using a
     * low-precision output format, such as GRAY8, will almost
     * certainly result in pixel overflow.  Using an integer output
     * format will introduce quantification errors in the convolution.
     * If you're low-pass filtering an integer-valued image, consider
     * using filterBinomial2D() instead.
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
             ImageFormat ImageFormat,
             class KernelType>
    void
    filter2D(
      Image<OutputFormat>& outputImage,
      const Kernel<KernelType>& kernel,
      const Image<ImageFormat>& image,
      const typename ImageFormatTraits<OutputFormat>::PixelType fillValue
      = typename ImageFormatTraits<OutputFormat>::PixelType(),
      ConvolutionStrategy convolutionStrategy = BRICK_CONVOLVE_PAD_RESULT);


    /** 
     * This function low-pass filters each column of an image with
     * integer-valued pixels using a binomial approximation to a
     * Gaussian kernel.  The result is returned through the
     * pre-allocated outputImage argument.  It is an error if the size
     * of the output image differs from the input image.  This is
     * typically much faster than filtering using a floating point
     * kernel.  The returned image is the same size as the input
     * image, and regions of the returned image for which the
     * convolution result is not valid (because the convolution kernel
     * extends off the edge of the image) will be initialized to
     * fillValue.
     * 
     * @param outputImage The result of the filtering will be returned
     * through this argument, which must have the same number of rows
     * and columns as inputImage.
     *
     * @param inputImage This argument is the Image to be filtered.
     *
     * @param sigma This argument specifies the desired filter
     * sigma.  It will be rounded to the nearest supported sigma value
     * (currently, 0.707, 1.0, 1.22, corresponding to filter sizes of
     * 3, 5, and 7 pixels, respectively) Specifying a sigma value that
     * doesn't round easily to one of these is an error.
     * 
     * @param fillValue This argument specifies the value with which
     * image edges should be padded.
     *
     * @param finalShift This argument gives you control over the
     * normalization of the filter kernel, but is slightly complicated
     * to use.  Leaving it at its default value scales output pixels
     * so that the effective filter kernel integrates to one.  Setting
     * it to zero disables normalization, and preserves full precision
     * in the result. Setting it to 1 divides the unnormalized pixels
     * by two.  Setting it to 2 divides by four, and so on.  Leaving
     * it at -1 is the same as specifying 2 for the 3 element binomial
     * filter, 4 for the 5 element filter, and 6 for the 7 element
     * filter.
     */
    template<ImageFormat OutputFormat, ImageFormat ImageFormat>
    void
    filterColumnsBinomial(
      Image<OutputFormat>& outputImage,
      const Image<ImageFormat>& inputImage,
      double sigma,
      typename ImageFormatTraits<OutputFormat>::PixelType const fillValue
      = typename ImageFormatTraits<OutputFormat>::PixelType(),
      int finalShift = -1);


    /** 
     * This function low-pass filters each row of an image with
     * integer-valued pixels using a binomial approximation to a
     * Gaussian kernel.  The result is returned through the
     * pre-allocated outputImage argument.  It is an error if the size
     * of the output image differs from the input image.  This is
     * typically much faster than filtering using a floating point
     * kernel.  The returned image is the same size as the input
     * image, and regions of the returned image for which the
     * convolution result is not valid (because the convolution kernel
     * extends off the edge of the image) will be initialized to
     * fillValue.
     * 
     * @param outputImage The result of the filtering will be returned
     * through this argument, which must have the same number of rows
     * and columns as inputImage.
     *
     * @param inputImage This argument is the Image to be filtered.
     *
     * @param rowSigma This argument specifies the desired filter
     * sigma.  It will be rounded to the nearest supported sigma value
     * (currently, 0.707, 1.0, 1.22, corresponding to filter sizes of
     * 3, 5, and 7 pixels, respectively) Specifying a sigma value that
     * doesn't round easily to one of these is an error.
     * 
     * @param columnSigma This argument specifies the desired filter
     * sigma.  It will be rounded to the nearest supported sigma value
     * (currently, 0.5 pixels, 1.0 pixels, 1.5 pixels).  Specifying a
     * sigma value more than 0.5 pixels away from the closest
     * supported value is an error.
     * 
     * @param fillValue This argument specifies the value with which
     * image edges should be padded.
     *
     * @param finalShift This argument gives you control over the
     * normalization of the filter kernel, but is slightly complicated
     * to use.  Leaving it at its default value scales output pixels
     * so that the effective filter kernel integrates to one.  Setting
     * it to zero disables normalization, and preserves full precision
     * in the result. Setting it to 1 divides the unnormalized pixels
     * by two.  Setting it to 2 divides by four, and so on.  Leaving
     * it at -1 is the same as specifying 2 for the 3 element binomial
     * filter, 4 for the 5 element filter, and 6 for the 7 element
     * filter.
     */
    template<ImageFormat OutputFormat, ImageFormat ImageFormat>
    void
    filterRowsBinomial(
      Image<OutputFormat>& outputImage,
      const Image<ImageFormat>& inputImage,
      double sigma,
      typename ImageFormatTraits<OutputFormat>::PixelType const fillValue
      = typename ImageFormatTraits<OutputFormat>::PixelType(),
      int finalShift = -1);

  } // namespace computerVision
  
} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/imageFilter_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_IMAGEFILTER_HH */
