/**
***************************************************************************
* @file brick/computerVision/canny_impl.hh
*
* Header file defining inline and template functions declared in canny.hh.
*
* Copyright (C) 2006,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_CANNY_IMPL_HH
#define BRICK_COMPUTERVISION_CANNY_IMPL_HH

// This file is included by canny.hh, and should not be directly included
// by user code, so no need to include canny.hh here.
//
// #include <brick/computerVision/canny.hh>

#include <limits>
#include <list>
#include <brick/computerVision/imageFilter.hh>
#include <brick/computerVision/kernels.hh>
#include <brick/computerVision/nonMaximumSuppress.hh>
#include <brick/computerVision/sobel.hh>
#include <brick/computerVision/utilities.hh>
#include <brick/numeric/index2D.hh>
#include <brick/numeric/utilities.hh>

namespace brick {

  namespace computerVision {

    /// @cond privateCode
    namespace privateCode {

      template <class FloatType, ImageFormat Format>
      Image<GRAY1>
      traceEdges(const Image<Format>& gradImage, FloatType threshold)
      {
        Image<GRAY1> edgeImage(gradImage.rows(), gradImage.columns());
        edgeImage = false;

        // Make a pass through the image pushing all the seed points
        // onto our list of edges to trace.
        std::list<brick::numeric::Index2D> seedList;
        size_t index0 = 0;
        for(size_t row = 0; row < gradImage.rows(); ++row) {
          for(size_t column = 0; column < gradImage.columns(); ++column) {
            if(gradImage[index0] > threshold) {
              edgeImage[index0] = true;
              seedList.push_front(brick::numeric::Index2D(static_cast<int>(row),
                                          static_cast<int>(column)));
            }
            ++index0;
          }
        }

        // Work our way through all edges that must be traced.
        size_t lastRow = gradImage.rows() - 1;
        size_t lastColumn = gradImage.columns() - 1;
        size_t columns = gradImage.columns();
        while(seedList.size() != 0) {

          // Get the next edge to trace.
          brick::numeric::Index2D seed = *seedList.begin();
          seedList.pop_front();

          // Don't trace edges on the very outside border of the image.
          size_t row = seed.getRow();
          size_t column = seed.getColumn();
          if(row == 0 || column == 0
             || row == lastRow || column == lastColumn) {
            continue;
          }

          // Use single indexing on the assumption that it's faster.
          index0 = row * columns + column;

          // Inspect each neighbor in turn, marking as appropriate,
          // and pushing any newly marked neighbor pixels onto the
          // list so that they will themselves be traced.
          size_t neighborIndex = index0 - columns - 1;
          if(gradImage(neighborIndex) != 0.0
             && edgeImage(neighborIndex) == false) {
            edgeImage(neighborIndex) = true;
            seedList.push_front(brick::numeric::Index2D(static_cast<int>(row) - 1,
                                        static_cast<int>(column) - 1));
          }
          ++neighborIndex;
          if(gradImage(neighborIndex) != 0.0
             && edgeImage(neighborIndex) == false) {
            edgeImage(neighborIndex) = true;
            seedList.push_front(brick::numeric::Index2D(static_cast<int>(row) - 1,
                                        static_cast<int>(column)));
          }
          ++neighborIndex;
          if(gradImage(neighborIndex) != 0.0
             && edgeImage(neighborIndex) == false) {
            edgeImage(neighborIndex) = true;
            seedList.push_front(brick::numeric::Index2D(static_cast<int>(row) - 1,
                                        static_cast<int>(column) + 1));
          }
          neighborIndex = index0 - 1;
          if(gradImage(neighborIndex) != 0.0
             && edgeImage(neighborIndex) == false) {
            edgeImage(neighborIndex) = true;
            seedList.push_front(brick::numeric::Index2D(static_cast<int>(row),
                                        static_cast<int>(column) - 1));
          }
          neighborIndex = index0 + 1;
          if(gradImage(neighborIndex) != 0.0
             && edgeImage(neighborIndex) == false) {
            edgeImage(neighborIndex) = true;
            seedList.push_front(brick::numeric::Index2D(static_cast<int>(row),
                                        static_cast<int>(column) + 1));
          }
          neighborIndex = index0 + columns - 1;
          if(gradImage(neighborIndex) != 0.0
             && edgeImage(neighborIndex) == false) {
            edgeImage(neighborIndex) = true;
            seedList.push_front(brick::numeric::Index2D(static_cast<int>(row) + 1,
                                        static_cast<int>(column) - 1));
          }
          ++neighborIndex;
          if(gradImage(neighborIndex) != 0.0
             && edgeImage(neighborIndex) == false) {
            edgeImage(neighborIndex) = true;
            seedList.push_front(brick::numeric::Index2D(static_cast<int>(row) + 1,
                                        static_cast<int>(column)));
          }
          ++neighborIndex;
          if(gradImage(neighborIndex) != 0.0
             && edgeImage(neighborIndex) == false) {
            edgeImage(neighborIndex) = true;
            seedList.push_front(brick::numeric::Index2D(static_cast<int>(row) + 1,
                                        static_cast<int>(column) + 1));
          }
        }
        return edgeImage;
      }

    } // namespace privateCode
    /// @endcond


    // This function applies the canny edge operator.
    template <class FloatType, ImageFormat FORMAT>
    Image<GRAY1>
    applyCanny(const Image<FORMAT>& inputImage,
               unsigned int gaussianSize,
               FloatType upperThreshold,
               FloatType lowerThreshold,
               FloatType autoUpperThresholdFactor,
               FloatType autoLowerThresholdFactor)
    {
      brick::numeric::Array2D<FloatType> gradientX;
      brick::numeric::Array2D<FloatType> gradientY;
      return applyCanny(inputImage, gradientX, gradientY, gaussianSize,
                        upperThreshold, lowerThreshold,
                        autoUpperThresholdFactor, autoLowerThresholdFactor);
    }


    // This function applies the canny edge operator.
    template <class FloatType, ImageFormat FORMAT>
    Image<GRAY1>
    applyCanny(const Image<FORMAT>& inputImage,
               brick::numeric::Array2D<FloatType>& gradientX,
               brick::numeric::Array2D<FloatType>& gradientY,
               unsigned int gaussianSize,
               FloatType upperThreshold,
               FloatType lowerThreshold,
               FloatType autoUpperThresholdFactor,
               FloatType autoLowerThresholdFactor)
    {
      // Argument checking.
      if(inputImage.rows() < gaussianSize + 3
         || inputImage.columns() < gaussianSize + 3) {
        BRICK_THROW(brick::common::ValueException, "applyCanny()",
                  "Argument inputImage has insufficient size, or argument "
                  "gaussianSize is too large.");
      }
      if(lowerThreshold > upperThreshold) {
        BRICK_THROW(brick::common::ValueException, "applyCanny()",
                  "Argument lowerThreshold must be less than or equal to "
                  "Arguments upperThreshold.");
      }
      autoLowerThresholdFactor =
        std::min(autoLowerThresholdFactor, autoUpperThresholdFactor);

      // We use non-normalized convolution kernels for the gradient,
      // which means that our user-specified thresholds don't match
      // the magnitude of our gradients.  Each gradient component is
      // 8 times as large as it should be, so the magnitude of the
      // gradient is 8*sqrt(2) times as large.  We solve this by
      // scaling the thresholds here.
      FloatType scaleFactor = static_cast<FloatType>(std::sqrt(2.0) * 8.0);
      lowerThreshold *= scaleFactor;
      upperThreshold *= scaleFactor;

      // Step 1: Blur with a gaussian kernel to reduce noise.
      Image<ImageFormatIdentifierGray<FloatType>::Format> blurredImage;
      if(gaussianSize == 0) {
        blurredImage = convertColorspace<
          ImageFormatIdentifierGray<FloatType>::Format>(inputImage);
      } else {
        Kernel<FloatType> gaussian =
          getGaussianKernelBySize<FloatType>(gaussianSize, gaussianSize);
        blurredImage =
          filter2D<ImageFormatIdentifierGray<FloatType>::Format, FORMAT,
                   FloatType>(
                     gaussian, inputImage, 0.0);
      }

      // Step 2: Compute derivatives of the blurred image, and discard
      // any which are less than the lower threshold.
      Image<ImageFormatIdentifierGray<FloatType>::Format> gradX =
        applySobelX(blurredImage);
      Image<ImageFormatIdentifierGray<FloatType>::Format> gradY =
        applySobelY(blurredImage);
      Image<ImageFormatIdentifierGray<FloatType>::Format> gradMagnitude(
          gradX.rows(), gradX.columns());

      // Communicate gradients to the calling context.
      gradientX = gradX;
      gradientY = gradY;

      // Continue with Canny algorithm.
      if(lowerThreshold > 0.0 && upperThreshold > 0.0) {
        // Discard values less than the lower threshold.
        for(size_t index0 = 0; index0 < gradX.size(); ++index0) {
          FloatType tmpVal = brick::common::squareRoot(
            gradX[index0] * gradX[index0] + gradY[index0] * gradY[index0]);
          gradMagnitude[index0] = (tmpVal > lowerThreshold) ? tmpVal : 0.0;
        }
      } else {
        // Temporarily retain all gradient values.
        for(size_t index0 = 0; index0 < gradX.size(); ++index0) {
          gradMagnitude[index0] = brick::common::squareRoot(
            gradX[index0] * gradX[index0] + gradY[index0] * gradY[index0]);
        }

        // Pick edge thresholds.
        size_t startRow = (gaussianSize + 1) / 2;
        size_t endRow = gradMagnitude.rows() - startRow;
        size_t startColumn = startRow;
        size_t endColumn = gradMagnitude.columns() - startColumn;

        // Paranoid check should never fail.
        if((startRow >= endRow) || (startColumn >= endColumn)) {
          BRICK_THROW(brick::common::ValueException, "applyCanny()",
                    "Filter kernel is too large for image.");
        }

        // Hack(xxx): Zero out borders of images to avoid spurious edges.
        for(size_t row = 0; row < startRow; ++row) {
          for(size_t column = 0; column < gradMagnitude.columns(); ++column) {
            gradMagnitude(row, column) = 0.0;
          }
        }
        for(size_t row = endRow; row < gradMagnitude.rows(); ++row) {
          for(size_t column = 0; column < gradMagnitude.columns(); ++column) {
            gradMagnitude(row, column) = 0.0;
          }
        }
        for(size_t row = 0; row < gradMagnitude.rows(); ++row) {
          for(size_t column = 0; column < startColumn; ++column) {
            gradMagnitude(row, column) = 0.0;
          }
          for(size_t column = endColumn; column < gradMagnitude.columns();
              ++column) {
            gradMagnitude(row, column) = 0.0;
          }
        }

        // Compute mean and variance of gradient values.
        size_t numberOfPixels =
          (endRow - startRow) * (endColumn - startColumn);
        FloatType sumOfGradient = 0.0;
        FloatType sumOfGradientSquared = 0.0;
        for(size_t row = startRow; row < endRow; ++row) {
          FloatType subSum = 0.0;
          FloatType subSumOfSquares = 0.0;
          for(size_t column = startColumn; column < endColumn; ++column) {
            FloatType testValue = gradMagnitude(row, column);
            // Changing how we compute threshold...
            //
            // if(testValue < minGrad) {minGrad = testValue;}
            // if(testValue > maxGrad) {maxGrad = testValue;}
            subSum += testValue;
            subSumOfSquares += testValue * testValue;
          }
          sumOfGradient += subSum;
          sumOfGradientSquared += subSumOfSquares;
        }
        FloatType gradientMean = sumOfGradient / numberOfPixels;
        FloatType gradientVariance =
          sumOfGradientSquared / numberOfPixels - gradientMean * gradientMean;
        FloatType gradientSigma = brick::common::squareRoot(gradientVariance);

        if(upperThreshold <= 0.0) {
          upperThreshold =
            gradientMean + autoUpperThresholdFactor * gradientSigma;
          upperThreshold = std::max(upperThreshold, 0.0);
        }
        if(lowerThreshold <= 0.0) {
          lowerThreshold =
            gradientMean + autoLowerThresholdFactor * gradientSigma;
          lowerThreshold = std::min(lowerThreshold, upperThreshold);
          lowerThreshold = std::max(lowerThreshold, 0.0);
        }

        // Now zero out gradients that for sure can never be edges.
        for(size_t index0 = 0; index0 < gradX.size(); ++index0) {
          FloatType tmpVal = gradMagnitude[index0];
          gradMagnitude[index0] = (tmpVal > lowerThreshold) ? tmpVal : 0.0;
        }
      }

      // Step 3: Non-maximum suppression.
      Image<ImageFormatIdentifierGray<FloatType>::Format> edgeCandidates =
        nonMaximumSuppress(gradMagnitude, gradX, gradY);

      // Step 4: Threshold with hysteresis.
      Image<GRAY1> edgeImage =
        privateCode::traceEdges(edgeCandidates, upperThreshold);
      return edgeImage;
    }

  } // namespace computerVision

} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_CANNY_IMPL_HH */
