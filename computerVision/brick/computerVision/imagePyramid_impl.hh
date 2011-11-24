/**
***************************************************************************
* @file brick/computerVision/imagePyramid_impl.hh
*
* Header file defining a class template for constructing scale-space
* image pyramids.
*
* Copyright (C) 2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_IMAGEPYRAMID_IMPL_HH
#define BRICK_COMPUTERVISION_IMAGEPYRAMID_IMPL_HH

// This file is included by imagePyramid.hh, and should not be directly included
// by user code, so no need to include imagePyramid.hh here.
// 
// #include <brick/computerVision/imagePyramid.hh>

#include <brick/computerVision/imageFilter.hh>
#include <brick/computerVision/kernels.hh>
#include <brick/computerVision/utilities.hh>

namespace brick {

  namespace computerVision {

    template <ImageFormat Format, ImageFormat InternalFormat, class KernelType>
    ImagePyramid<Format, InternalFormat, KernelType>::
    ImagePyramid(Image<Format> const& inputImage,
                 double scaleFactorPerLevel,
                 unsigned int levels,
                 bool isBandPass)
      : m_pyramid()
    {
      if(scaleFactorPerLevel < 1.0) {
        std::ostringstream message;
        message << "Argument scaleFactorPerLevel must be greater than "
                << "or equal to 1.0, but has value " << scaleFactorPerLevel
                << ".";
        BRICK_THROW(brick::common::IndexException,
                    "ImagePyramid::ImagePyramid()", message.str().c_str());
      }

      // This constructor doesn't allow a user defined filter
      // function, so we supply an automatically generated Gaussian
      // low-pass filter.  We want the filter to cut off anything
      // higher than Nyquist for the resampled image, prior to the
      // resampling, so that we don't introduce aliasing artifacts.
      // The wavelength of the Nyquist frequency in the resampled
      // image is 2 pixels, which translates to a wavelength of 2.0 *
      // scaleFactorPerLevel in the unresampled image.
      double cutoffWavelength = 2.0 * scaleFactorPerLevel;

      // We're going to use a Gaussian low-pass filter.  We'd like to
      // attenuate frequencies above out cutoff by at least a factor
      // of 2.  The filter will have the form
      //
      //   f(x) = k_0 * exp(-(x^2)/(2.0 * sigma^2)).
      //
      // In this case, we define x to be in units of pixels.  The
      // Fourier transform of this filter is
      //
      //   F(w) = k_2 * exp(-(w^2) * (sigma^2)/2.0).
      //
      // Setting w to 2*pi / cutoffWavelength, and F(w) = 0.5 * k_2, we have
      //
      //   0.5 = exp(-2.0 * pi^2 * sigma^2 / cutoffWavelength^2),
      //
      // which gives us
      //
      //   ln(0.5) = -2.0 * pi^2 * sigma^2 / cutoffWavelength^2
      //
      //   sigma = (0.8326 / pi) * cutoffWavelength = 0.265 * cutoffWavelength.
      double sigma = 0.265 * cutoffWavelength;
      Kernel<KernelType> filterKernel = getGaussianKernel<KernelType>(
        sigma, sigma);

      // If levels is set to zero, 
      if(0 == levels) {
        // Filter kernel may extend up to 6 sigma in any direction.
        double minimumImageSize = 12 * sigma;

        // How many pyramid levels?  Well, enough that we get close to
        // minimumImageSize, but no smaller.  This implies that
        // minimumImageSize * scaleFactorPerLevel^numberOfLevels is
        // less than input image size, but minimumImageSize *
        // scaleFactorPerLevel^(numberOfLevels + 1) is greater than
        // input image size.  That is, numberOfLevels is the largest
        // integer less than variable nn in the following equation:
        //
        //   minImageSize * scaleFactor^nn = inImageSize
        //
        //   scaleFactor^nn = inImageSize / minImageSize
        //
        //   nn = ln(inImageSize / minImageSize) / ln(scaleFactor)
        unsigned int limitingDimension =
          std::min(inputImage.rows(), inputImage.columns());
        unsigned int targetSizeRatio = limitingDimension / minimumImageSize;
        levels = static_cast<unsigned int>(
          std::log(targetSizeRatio) / std::log(scaleFactorPerLevel));
        levels = std::max(levels, static_cast<unsigned int>(0));
      }

      Image<Format> currentImage = inputImage.copy();

      // Start off the pyramid.
      m_pyramid.push_back(currentImage);
      --levels;
      
      while(0 != levels) {
        Image<Format> filteredImage = filter2D<Format, InternalFormat>(
          filterKernel, currentImage);
        if(isBandPass) {
          m_pyramid[m_pyramid.size() - 1] -= filteredImage;
        }
        currentImage = this->subsampleImage(filteredImage, scaleFactorPerLevel);
        m_pyramid.push_back(currentImage);
        --levels;
      }
    }

      
    template <ImageFormat Format, ImageFormat InternalFormat, class KernelType>
    Image<Format>&
    ImagePyramid<Format, InternalFormat, KernelType>::
    getLevel(unsigned int levelIndex)
    {
      if(levelIndex >= m_pyramid.size()) {
        std::ostringstream message;
        message << "Argument levelIndex (" << levelIndex << ") "
                << "is invalid for a " << m_pyramid.size()
                << " level image pyramid.";
        BRICK_THROW(brick::common::IndexException,
                    "ImagePyramid::getLevel()", message.str().c_str());
      }
      return m_pyramid[levelIndex];
    }


    template <ImageFormat Format, ImageFormat InternalFormat, class KernelType>
    unsigned int
    ImagePyramid<Format, InternalFormat, KernelType>::
    getNumberOfLevels()
    {
      return m_pyramid.size();
    }
    

    // ============== Private member functions below this line ==============

    template <ImageFormat Format, ImageFormat InternalFormat, class KernelType>
    bool
    ImagePyramid<Format, InternalFormat, KernelType>::
    isIntegral(double scaleFactor, int& integerScaleFactor)
    {
      // If scaleFactor is within 0.01% of being integral, call it
      // close enough.
      integerScaleFactor = static_cast<int>(scaleFactor + 0.5);
      if(static_cast<double>(integerScaleFactor) == scaleFactor) {
        return true;
      }
      return false;
    }


    template <ImageFormat Format, ImageFormat InternalFormat, class KernelType>
    Image<Format>
    ImagePyramid<Format, InternalFormat, KernelType>::
    subsampleImage(Image<Format> const& inputImage,
                   double scaleFactor)
    {
      int integralScaleFactor = 0;
      if(isIntegral(scaleFactor, integralScaleFactor)) {
        return subsample(inputImage, integralScaleFactor, integralScaleFactor);
      }
      BRICK_THROW(brick::common::NotImplementedException,
                  "ImagePyramid::subsampleImage()",
                  "Non-integral scale factors are not yet supported.");
      // return resampleImageGeneralPosition(inputImage, scaleFactor);
      return subsample(inputImage, integralScaleFactor, integralScaleFactor);
    }


  } // namespace computerVision
  
} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_IMAGEPYRAMID_IMPL_HH */
