/**
***************************************************************************
* @file brick/computerVision/histogramEqualize.cc
*
* Source file defining histogram equalization routines.
*
* Copyright (C) 2005,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#include <limits>
#include <numeric>
#include <sstream>
#include <brick/computerVision/histogramEqualize.hh>

namespace brick {

  namespace computerVision {
    
    // This function computes the histogram of an image.
    numeric::Array1D<unsigned int>
    getHistogram(const Image<GRAY8>& inputImage)
    {
      if(inputImage.size() > std::numeric_limits<unsigned int>::max()) {
        std::ostringstream message;
        message << "Currently, we can only equalize images with "
                << std::numeric_limits<int>::max() << " or fewer pixels.";
        BRICK_THROW(common::ValueException, "histogramEqualize()", message.str().c_str());
      }

      numeric::Array1D<unsigned int> histogram(
        std::numeric_limits<Image<GRAY8>::PixelType>::max());
      histogram = 0;
      for(size_t pixelIndex = 0; pixelIndex < inputImage.size(); ++pixelIndex) {
        ++(histogram[inputImage[pixelIndex]]);
      }
      return histogram;
    }
  
  
    // This function remaps the pixel values of the input image in such
    // a way that output pixel value increases monotonically with input
    // pixel value, and the histogram of the output image is nearly
    // flat.
    Image<GRAY8>
    histogramEqualize(const Image<GRAY8>& inputImage)
    {
      // Compute the histogram and CDF.
      numeric::Array1D<unsigned int> histogram = getHistogram(inputImage);
      numeric::Array1D<unsigned int> cdf(histogram.size());
      std::partial_sum(histogram.begin(), histogram.end(), cdf.begin(),
                       std::plus<unsigned int>());

      // Rescale the image according to the CDF.
      Image<GRAY8> outputImage(inputImage.rows(), inputImage.columns());
      double scaleFactor = 256.0 / (inputImage.size() + 1);
      for(size_t pixelIndex = 0; pixelIndex < inputImage.size(); ++pixelIndex) {
        outputImage[pixelIndex] =
          static_cast<common::Int8>(scaleFactor * cdf[inputImage[pixelIndex]]);
      }
      return outputImage;
    }

  } // namespace computerVision    

} // namespace brick
