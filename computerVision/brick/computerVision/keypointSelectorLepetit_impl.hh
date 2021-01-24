/**
***************************************************************************
* @file brick/computerVision/keypointSelectorLepetit_impl.hh
*
* Header file defining a class template for selecting stable keypoints
* from an image.
*
* Copyright (C) 2011 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_KEYPOINTSELECTORLEPETIT_IMPL_HH
#define BRICK_COMPUTERVISION_KEYPOINTSELECTORLEPETIT_IMPL_HH

// This file is included by keypointSelectorLepetit.hh, and should not be directly included
// by user code, so no need to include keypointSelectorLepetit.hh here.
//
// #include <brick/computerVision/keypointSelectorLepetit.hh>

#include <algorithm>
#include <brick/common/mathFunctions.hh>
#include <brick/computerVision/imagePyramid.hh>


namespace brick {

  namespace computerVision {

    template <ImageFormat Format>
    void
    KeypointSelectorLepetit::
    setImage(Image<Format> const& image)
    {
      // Nearly all of the following computations are easier with a
      // floating point representation.  Speed of development is the
      // highest priority right now, so we take the hit of going to
      // floating point.
      Image<GRAY_FLOAT32> floatImage = convertColorspace<GRAY_FLOAT32>(image);

      // Lepetit's paper suggests that they've simply hard-coded some
      // constants in their keypoint selection algorithm.  Here we try
      // to estimate appropriate values for those constants
      // automatically.
      this->estimateThresholds(
        floatImage, m_thresholdPixelSimilarity, m_thresholdLaplacianMagnitude);

      // Create an image pyramid, one image per octave, computing only
      // the low-pass filtered images (no need to generate
      // difference-of-Gaussian images).  For now, we'll just run on
      // the first 3 octaves, so we only need to generate three
      // pyramid levels.
      ImagePyramid<GRAY_FLOAT32, GRAY_FLOAT32, common::Float32> pyramid(
        floatImage, 2.0, 3, false);

      // Iterate over all pyramid levels looking for keypoints.
      for(unsigned int level = 0; level < pyramid.getNumberOfLevels();
          ++level) {
        Image<GRAY_FLOAT32> currentLevel = pyramid.getLevel(level);

        // We will iterate over all pixels that have a complete set of
        // valid neighbors.  This means we need to stay at least one
        // pixel away from the borders of the image, where the
        // filtered image data is no longer valid.
        unsigned int startRow = pyramid.getBorderSizeTopBottom() + 1;
        unsigned int stopRow = (currentLevel.rows()
                                - pyramid.getBorderSizeTopBottom() - 1);
        unsigned int startColumn = pyramid.getBorderSizeLeftRight() + 1;
        unsigned int stopColumn = (currentLevel.columns()
                                   - pyramid.getBorderSizeLeftRight() - 1);
        for(unsigned int row = startRow; row < stopRow; ++row) {
          for(unsigned int column = startColumn; column < stopColumn;
              ++column) {
            if(this->testPixel(
                 row, column, currentLevel, m_thresholdPixelSimilarity,
                 m_thresholdLaplacianMagnitude)) {
              brick::numeric::Index2D localKeypoint(row, column);
              brick::numeric::Index2D keyPoint =
                pyramid.convertImageCoordinates(localKeypoint, level, 0);
              m_locationVector.push_back(keyPoint);
              m_levelVector.push_back(level);
            }
          }
        }
      }
    }


    void
    KeypointSelectorLepetit::
    estimateThresholds(Image<GRAY_FLOAT32> const& image,
                       float& pixelSimilarityThreshold,
                       float& laplacianMagnitudeThreshold)
    {
      const unsigned int expectedKeypointsPerImage = 300;
      const unsigned int sparsity = 10;
      const unsigned int pixelMeasurementRadius = 3;

      // Sparsely sample over entire image to accumulate statistics
      // about what the normal pixelSimilarity and laplacianMagnitude
      // are.  Sampling doesn't go all the way to edges because the
      // stats require a 7x7 window to compute.
      std::vector<float> pixelSimilarityVector;
      std::vector<float> laplacianMagnitudeVector;
      for(unsigned int row = pixelMeasurementRadius;
          row < image.rows() - pixelMeasurementRadius;
          row += sparsity) {
        for(unsigned int column = pixelMeasurementRadius;
            column < image.columns() - pixelMeasurementRadius;
            column += sparsity) {
          float pixelSimilarity = 0.0;
          float laplacianMagnitude = 0.0;
          this->measurePixelThresholds(
            row, column, image, pixelSimilarity, laplacianMagnitude);
          pixelSimilarityVector.push_back(pixelSimilarity);
          laplacianMagnitudeVector.push_back(laplacianMagnitude);
        }
      }

      // Figure out how many keypoints we can reasonably expect to see
      // in our sparse vector of samples.
      unsigned int numPixels = image.size();
      unsigned int numSamples = pixelSimilarityVector.size();
      double proportionSampled = double(numSamples) / numPixels;
      unsigned int numKeypoints = proportionSampled * expectedKeypointsPerImage;
      numKeypoints = std::min(numKeypoints, numSamples);
      numKeypoints = std::max(numKeypoints, static_cast<unsigned int>(1));
      unsigned int numInliers = numSamples - numKeypoints;

      // Sort our two vectors of samples.  This lets us more easily
      // find threshold values that separate the keypoints from the
      // non-keypoints.
      std::sort(pixelSimilarityVector.begin(),
                pixelSimilarityVector.end());
      std::sort(laplacianMagnitudeVector.begin(),
                laplacianMagnitudeVector.end());

      // Discard the largest values of pixelSimilarity and
      // laplacianMagnitude, as those probably come from actual
      // keypoints that we accidentally sampled.  Assume that the
      // remaining values represent "normal" for non-keypoint pixels.
#if 1
      pixelSimilarityThreshold = pixelSimilarityVector[numInliers];
      laplacianMagnitudeThreshold = laplacianMagnitudeVector[numInliers];
#else
      // Compute statistics over all of the "inliers" so we can
      // build a model of what normal pixels look like.
      float pixelSimilaritySum = 0.0;
      float laplacianMagnitudeSum = 0.0;
      float pixelSimilaritySumOfSquares = 0.0;
      float laplacianMagnitudeSumOfSquares = 0.0;
      for(unsigned int ii = 0; ii < numInliers; ++ii) {
        float pixelSimilarity = pixelSimilarityVector[ii];
        float laplacianMagnitude = laplacianMagnitudeVector[ii];
        pixelSimilaritySum += pixelSimilarity;
        pixelSimilaritySumOfSquares += pixelSimilarity * pixelSimilarity;
        laplacianMagnitudeSum += laplacianMagnitude;
        laplacianMagnitudeSumOfSquares +=
          laplacianMagnitude * laplacianMagnitude;
      }
      float pixelSimilarityMean = pixelSimilaritySum / numInliers;
      float laplacianMagnitudeMean = laplacianMagnitudeSum / numInliers;
      float pixelSimilarityVariance = (
        pixelSimilaritySumOfSquares / numInliers
        - pixelSimilarityMean * pixelSimilarityMean);
      float laplacianMagnitudeVariance = (
        laplacianMagnitudeSumOfSquares / numInliers
        - laplacianMagnitudeMean * laplacianMagnitudeMean);

      // Finally, set thresholds at mean + 3*sigma, so that (assuming
      // Gaussian distribution), we discard most of the inliers.
      pixelSimilarityThreshold =
        pixelSimilarityMean + 3.0 * std::sqrt(pixelSimilarityVariance);
      laplacianMagnitudeThreshold =
        laplacianMagnitudeMean + 3.0 * std::sqrt(laplacianMagnitudeVariance);
#endif
    }


    void
    KeypointSelectorLepetit::
    measurePixelThresholds(unsigned int row, unsigned int column,
                           Image<GRAY_FLOAT32> const& image,
                           float& pixelSimilarity,
                           float& laplacianMagnitude)
    {
      float testValue = image(row, column);

      // Find the value of pixelSimilarity at which this
      // pixel makes the transition from failing to passing.  That is,
      // find the highest threshold at which this pixel still passes the
      // test.  For more information, see member function testPixel().
      float similarityValues[8];
      similarityValues[0] = std::max(
        common::absoluteValue(testValue - image(row, column - 3)),
        common::absoluteValue(testValue - image(row, column + 3)));
      similarityValues[1] = std::max(
        common::absoluteValue(testValue - image(row - 3, column)),
        common::absoluteValue(testValue - image(row + 3, column)));
      similarityValues[2] = std::max(
        common::absoluteValue(testValue - image(row + 2, column + 2)),
        common::absoluteValue(testValue - image(row - 2, column - 2)));
      similarityValues[3] = std::max(
        common::absoluteValue(testValue - image(row - 2, column + 2)),
        common::absoluteValue(testValue - image(row + 2, column - 2)));
      similarityValues[4] = std::max(
        common::absoluteValue(testValue - image(row + 1, column + 3)),
        common::absoluteValue(testValue - image(row - 1, column - 3)));
      similarityValues[5] = std::max(
        common::absoluteValue(testValue - image(row - 3, column + 1)),
        common::absoluteValue(testValue - image(row + 3, column - 1)));
      similarityValues[6] = std::max(
        common::absoluteValue(testValue - image(row + 3, column + 1)),
        common::absoluteValue(testValue - image(row - 3, column - 1)));
      similarityValues[7] = std::max(
        common::absoluteValue(testValue - image(row - 1, column + 3)),
        common::absoluteValue(testValue - image(row + 1, column -3)));

      // Select the smallest of the measured max values, 'cause that's
      // the one that would get this pixel thrown out if the threshold
      // were too high.
      pixelSimilarity =
        *(std::min_element(&(similarityValues[0]), &(similarityValues[0]) + 8));

      // Find the Laplacian-of-Gaussian approximation for this pixel.
      // This is the largest value of laplacianMagnitudeThreshold at
      // which this pixel would pass the test.
      float averageBorderValue = (
        image(row - 3, column - 1) + image(row - 3, column)
        + image(row - 3, column + 1)
        + image(row - 2, column - 2) + image(row - 2, column + 2)
        + image(row - 1, column - 3) + image(row - 1, column + 3)
        + image(row, column - 3) + image(row, column + 3)
        + image(row + 1, column - 3) + image(row + 1, column + 3)
        + image(row + 2, column - 2) + image(row + 2, column + 2)
        + image(row + 3, column - 1) + image(row + 3, column)
        + image(row + 3, column + 1)) / 16.0;
      laplacianMagnitude = brick::common::absoluteValue(
        testValue - averageBorderValue);
    }


    bool
    KeypointSelectorLepetit::
    testPixel(unsigned int row, unsigned int column,
              Image<GRAY_FLOAT32> const& image,
              float pixelSimilarityThreshold,
              float laplacianMagnitudeThreshold)
    {
      float testValue = image(row, column);

      // Counting on short-circuit evaluation to save us time here.
      // This giant if clause passes wnen the current pixel's gray
      // level is "...close to those of any two diametrically opposed
      // pixels on ... [the discretized circle surrounding the current
      // pixel, and is] ... therefore unlikely to be stable ..."
      //
      // Following Rosten's FAST feature detector (ICCV '05), we use a
      // Bresenham circle of radius 3 for the "discretized circle" of
      // Lepetit's description.
      if(((common::absoluteValue(testValue - image(row, column - 3))
           < pixelSimilarityThreshold)
          && (common::absoluteValue(testValue - image(row, column + 3))
              < pixelSimilarityThreshold))
         || ((common::absoluteValue(testValue - image(row - 3, column))
              < pixelSimilarityThreshold)
             && (common::absoluteValue(testValue - image(row + 3, column))
                 < pixelSimilarityThreshold))
         || ((common::absoluteValue(testValue - image(row + 2, column + 2))
              < pixelSimilarityThreshold)
             && (common::absoluteValue(testValue - image(row - 2, column - 2))
                 < pixelSimilarityThreshold))
         || ((common::absoluteValue(testValue - image(row - 2, column + 2))
              < pixelSimilarityThreshold)
             && (common::absoluteValue(testValue - image(row + 2, column - 2))
                 < pixelSimilarityThreshold))
         || ((common::absoluteValue(testValue - image(row + 1, column + 3))
              < pixelSimilarityThreshold)
             && (common::absoluteValue(testValue - image(row - 1, column - 3))
                 < pixelSimilarityThreshold))
         || ((common::absoluteValue(testValue - image(row - 3, column + 1))
              < pixelSimilarityThreshold)
             && (common::absoluteValue(testValue - image(row + 3, column - 1))
                 < pixelSimilarityThreshold))
         || ((common::absoluteValue(testValue - image(row + 3, column + 1))
              < pixelSimilarityThreshold)
             && (common::absoluteValue(testValue - image(row - 3, column - 1))
                 < pixelSimilarityThreshold))
         || ((common::absoluteValue(testValue - image(row - 1, column + 3))
              < pixelSimilarityThreshold)
             && (common::absoluteValue(testValue - image(row + 1, column -3))
                 < pixelSimilarityThreshold))) {
        // Unstable keypoint.  Discard it!
        return false;
      }

      // Candidate keypoint passed the "likely to be stable" test.
      // Now let's see if it's at a location with sufficiently high
      // approximate LoG value.
      float averageBorderValue = (
        image(row - 3, column - 1) + image(row - 3, column)
        + image(row - 3, column + 1)
        + image(row - 2, column - 2) + image(row - 2, column + 2)
        + image(row - 1, column - 3) + image(row - 1, column + 3)
        + image(row, column - 3) + image(row, column + 3)
        + image(row + 1, column - 3) + image(row + 1, column + 3)
        + image(row + 2, column - 2) + image(row + 2, column + 2)
        + image(row + 3, column - 1) + image(row + 3, column)
        + image(row + 3, column + 1)) / 16.0;
      float logApproximation = brick::common::absoluteValue(
        testValue - averageBorderValue);
      if(logApproximation <= laplacianMagnitudeThreshold) {
        // Nope!  Laplacian of Gaussian approximation is too small.
        // Discard this keypoint candidate.
        return false;
      }

      // All tests passed!  Return true to indicate that a new
      // keypoint should be added at (row, column).
      return true;
    }


    // ============== Private member functions below this line ==============

  } // namespace computerVision

} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_KEYPOINTSELECTORLEPETIT_IMPL_HH */
