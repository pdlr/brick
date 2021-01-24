/**
***************************************************************************
* @file brick/computerVision/thresholderSauvola_impl.hh
*
* Header file defining inline and template functions declared in
* thresholderSauvola.hh.
*
* Copyright (C) 2017 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_THRESHOLDERSAUVOLA_IMPL_HH
#define BRICK_COMPUTERVISION_THRESHOLDERSAUVOLA_IMPL_HH

// This file is included by thresholderSauvola.hh, and should not be
// directly included by user code, so no need to include
// thresholderSauvola.hh here.
//
// #include <brick/computerVision/thresholderSauvola.hh>

namespace brick {

  namespace computerVision {

  } // namespace computerVision

} // namespace brick


/* ============ Definitions of inline & template functions ============ */


#include <cmath>
#include <limits>

namespace brick {

  namespace computerVision {

    // Constructor.
    template <ImageFormat Format, class Config>
    ThresholderSauvola<Format, Config>::
    ThresholderSauvola(
      uint32_t windowRadius,
      ThresholderSauvola<Format, Config>::FloatType kappa)

      : m_kappa(kappa),
        m_windowRadius(windowRadius),
        m_inputImage(),
        m_sumIntegrator(),
        m_squaredSumIntegrator()
    {
      // Empty.
    }



    // Computes a thresholded image based on the input image
    // previously set by member function setImage().
    template <ImageFormat Format, class Config>
    Image<GRAY8>
    ThresholderSauvola<Format, Config>::
    computeBinaryImage()
    {
      // Precompute some values related to the size of the ROI.
      uint32_t const windowSize = this->getWindowSize();
      FloatType const windowArea =
        static_cast<FloatType>(windowSize * windowSize);

      // This is the largest possible value for standard deviation of
      // pixel intensity.  Generally 1/2 of the maximum pixel value.
      FloatType const maxStdDev = Config::getMaxStdDev();

      // We already know (because of checks in setImage() and
      // setWindowRadius()) that the image size is larger than
      // windowSize x windowSize.  This saves us from some error
      // checking below.
      uint32_t const totalRows = this->m_inputImage.getRows();
      uint32_t const totalColumns = this->m_inputImage.getColumns();

      // Allocate space for our return value;
      Image<GRAY8> outputImage(totalRows, totalColumns);

      // Process each pixel in turn.
      for(int32_t rr = 0; rr < static_cast<int32_t>(totalRows); ++rr) {

        // For now, tolerate the inefficiency of adjusting the window
        // each time (so it fits entirely within the image).  Later,
        // we'll add specialized loops for the borders of the image so
        // we don't have to do this check.
        int32_t roiBeginRow = rr - static_cast<int32_t>(this->m_windowRadius);
        int32_t roiEndRow = roiBeginRow + windowSize;
        this->adjustWindowCoordinates(
          roiBeginRow, roiEndRow, 0,
          static_cast<int32_t>(totalRows),
          static_cast<int32_t>(windowSize));

        for(int32_t cc = 0; cc < static_cast<int32_t>(totalColumns); ++cc) {

          // For now, tolerate the inefficiency of adjusting the window
          // each time (so it fits entirely within the image).  Later,
          // we'll add specialized loops for the borders of the image so
          // we don't have to do this check.
          int32_t roiBeginColumn =
            cc - static_cast<int32_t>(this->m_windowRadius);
          int32_t roiEndColumn = roiBeginColumn + windowSize;
          this->adjustWindowCoordinates(
            roiBeginColumn, roiEndColumn, 0,
            static_cast<int32_t>(totalColumns),
            static_cast<int32_t>(windowSize));

          // The algorithm needs the mean pixel value in the window.
          SumType pixelSum = this->m_sumIntegrator.getIntegral(
            brick::numeric::Index2D(roiBeginRow, roiBeginColumn),
            brick::numeric::Index2D(roiEndRow, roiEndColumn));
          FloatType localMean = static_cast<FloatType>(pixelSum) / windowArea;

          // We'll use the unbiased estimator for variance,
          //
          // @code
          //    var = (sum_i(x_i * x_i) - N * mu * mu) / (N - 1)
          // @endcode
          //
          // where mu is the sample mean, x_i are individual samples,
          // N is the nmber of samples, and sum_i denotes summation
          // over all of the samples.  Please see
          // brick::numeric::getMeanAndVariance() for a derivation of
          // this estimator.
          FloatType localVariance = static_cast<FloatType>(
            this->m_squaredSumIntegrator.getIntegral(
              brick::numeric::Index2D(roiBeginRow, roiBeginColumn),
              brick::numeric::Index2D(roiEndRow, roiEndColumn)));
          localVariance -= (localMean * localMean * windowArea);
          localVariance /= static_cast<FloatType>(windowArea - 1);

          // For now, we need a call to sqrt here.  Later this can be
          // optimized out for small integer pixel types.
          FloatType localStdDev = brick::common::squareRoot(localVariance);

          // Here's Sauvola's threshold rule.
          FloatType normalizedStdDev = localStdDev / maxStdDev;
          FloatType multiplier = normalizedStdDev - FloatType(1);
          FloatType scaleFactor = FloatType(1) + this->m_kappa * multiplier;
          FloatType threshold = localMean * scaleFactor;

          if(static_cast<FloatType>(this->m_inputImage(rr, cc)) > threshold) {
            outputImage(rr, cc) = Config::getWhiteValue();
          } else {
            outputImage(rr, cc) = Config::getBlackValue();
          }
        }
      }

      return outputImage;
    }


    // Sets the image to be thresholded, and does preprocessing so
    // that subsequent calls to member function computeBinaryImage()
    // can execute quickly.
    template <ImageFormat Format, class Config>
    void
    ThresholderSauvola<Format, Config>::
    setImage(Image<Format> const& inputImage)
    {
      if(inputImage.rows() < this->getWindowSize()
         || inputImage.columns() < this->getWindowSize()) {
        std::ostringstream message;
        message << "Input image size (" << inputImage.rows()
                << ", " << inputImage.columns() << ") is not large enough to "
                << "accommodate window size of (" << this->getWindowSize()
                << ", " << this->getWindowSize() << ").";
        BRICK_THROW(brick::common::ValueException,
                    "ThresholderSauvola::setImage()",
                    message.str().c_str());
      }

      m_inputImage = inputImage;
      m_sumIntegrator.setArray(inputImage);
      m_squaredSumIntegrator.setArray(
        inputImage, [](PixelType const& xx) {return xx * xx;});
    }

  } // namespace computerVision

} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_THRESHOLDERSAUVOLA_IMPL_HH */
