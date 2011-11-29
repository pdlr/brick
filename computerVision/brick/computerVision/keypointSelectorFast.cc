/**
***************************************************************************
* @file brick/computerVision/keypointSelectorFast.cc
*
* Source file defining a class for selecting stable keypoints from an
* image.
*
* Copyright (C) 2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#include <brick/computerVision/keypointSelectorFast.hh>

namespace brick {

  namespace computerVision {

    KeypointSelectorFast::
    KeypointSelectorFast()
      : m_keypointVector(),
        m_threshold(0)
    {
      // Empty.
    }


    void
    KeypointSelectorFast::
    estimateThreshold(Image<GRAY8> const& inImage,
                      unsigned int startRow,
                      unsigned int startColumn,
                      unsigned int stopRow,
                      unsigned int stopColumn,
                      unsigned int expectedKeypointsPerImage)
    {
      const unsigned int sparsity = 10;
      const unsigned int pixelMeasurementRadius = 3;

      // Make sure the passed-in image bounds are legal.
      this->checkAndRepairRegionOfInterest(
        inImage, pixelMeasurementRadius, startRow, startColumn,
        stopRow, stopColumn);
      
      // Sparsely sample over entire image to accumulate statistics
      // about what the normal pixel intensity difference is.
      std::vector<common::Int16> intensityDifferenceVector;
      for(unsigned int row = startRow; row < stopRow; row += sparsity) {
        for(unsigned int column = startColumn; column < stopColumn;
            column += sparsity) {
          intensityDifferenceVector.push_back(
            this->measurePixelThreshold(inImage, row, column));
        }
      }

      // Figure out how many keypoints we can reasonably expect to see
      // in our sparse vector of samples.  
      unsigned int numPixels =
        (stopRow - startRow) * (stopColumn - startColumn);
      unsigned int numSamples = intensityDifferenceVector.size();
      double proportionSampled = double(numSamples) / numPixels;
      unsigned int numKeypoints = proportionSampled * expectedKeypointsPerImage;
      numKeypoints = std::min(numKeypoints, numSamples);
      numKeypoints = std::max(numKeypoints, static_cast<unsigned int>(0));

      // We'll call points that are _not_ keypoints "inliers."  Find
      // out how many of them we can expect.
      unsigned int numInliers = numSamples - numKeypoints;

      // Sort our vector of samples so that we can easily find the
      // threshold value that would separate the keypoints from the
      // non-keypoints (inliers).
      if(numInliers > 0) {
        std::sort(intensityDifferenceVector.begin(),
                  intensityDifferenceVector.end());
        m_threshold = intensityDifferenceVector[numInliers];
      } else {
        m_threshold = 0;
      }
    }

    
    std::vector<KeypointFast>
    KeypointSelectorFast::
    getKeypoints()
    {
      return m_keypointVector;
    }


    brick::common::Int16
    KeypointSelectorFast::
    getThreshold()
    {
      return m_threshold;
    }


    void
    KeypointSelectorFast::
    setImage(Image<GRAY8> const& inImage,
             unsigned int startRow,
             unsigned int startColumn,
             unsigned int stopRow,
             unsigned int stopColumn)
    {
      const unsigned int pixelMeasurementRadius = 3;

      // Make sure the passed-in image bounds are legal.
      this->checkAndRepairRegionOfInterest(
        inImage, pixelMeasurementRadius, startRow, startColumn,
        stopRow, stopColumn);
      
      if(m_threshold == 0) {
        // Rosten's paper suggests that they've simply hard-coded a
        // constant threshold in their keypoint selection algorithm.
        // Here we try to estimate appropriate values for those
        // constants automatically.
        this->estimateThreshold(inImage, startRow, startColumn,
                                stopRow, stopColumn);
      }

      // Test every pixel!
      KeypointFast keypoint;
      for(unsigned int row = startRow; row < stopRow; ++row) {
        for(unsigned int column = startColumn; column < stopColumn;
            ++column) {
          if(this->testPixel(inImage, row, column, m_threshold, keypoint)) {
            m_keypointVector.push_back(keypoint);
          }
        }
      }
    }


    void
    KeypointSelectorFast::
    checkAndRepairRegionOfInterest(Image<GRAY8> const& inImage,
                                   unsigned int pixelMeasurementRadius,
                                   unsigned int& startRow,
                                   unsigned int& startColumn,
                                   unsigned int& stopRow,
                                   unsigned int& stopColumn)
    {
      // The FAST feature detector uses a Bresenham circle of radius
      // 3, so this means we need to stay at least 3 pixels away from
      // the borders of the image.
      startRow = std::max(
        startRow, static_cast<unsigned int>(pixelMeasurementRadius));
      stopRow = std::min(
        stopRow, inImage.rows() - pixelMeasurementRadius);
      startColumn = std::max(
        startColumn, static_cast<unsigned int>(pixelMeasurementRadius));
      stopColumn = std::min(
        stopColumn, inImage.columns() - pixelMeasurementRadius);
    }


    common::Int16
    KeypointSelectorFast::
    measurePixelThreshold(Image<GRAY8> const& image,
                          unsigned int row, unsigned int column)
    {
      common::Int16 testValue = static_cast<common::Int16>(image(row, column));

      // Find the value of threshold at which this pixel makes the
      // transition from failing to passing.  Considering differences
      // between the 16 "neighbor" pixels and the center pixel
      // (testValue), we're looking for the largest threshold value
      // that's still exceeded by at least 12 consecutive neighbors,
      // with the added constraint that each of the 12 consecutive
      // neighbors must exceed the threshold in the same direction
      // (i.e., all 12 must be brighter than testValue, or all 12 must
      // be dimmer than testValue).  To do this, we need to find the
      // minimum differences in each direction among each neighborhood
      // of size 12.  If this doesn't make sense to you, you may
      // benefit from reading member function testPixel() before
      // continuing.

      // We'll accumulate min differences into this table so that we
      // don't have to keep recomputing them.  We've opted to
      // hard-code the dimensions here because they're integral to the
      // algorithm, and they make the following code a little clearer.
      common::Int16 minAccumulator[5][16];
      common::Int16 maxAccumulator[5][16];

      // First check all of the 16 neighbors and put differences into
      // the first row of min table.  Could probably save some cache
      // misses by rearranging order of access here.  Instead, we
      // start at the top of the circle and simply walk around it
      // clockwise.
      minAccumulator[0][0] = testValue - image(row - 3, column);
      minAccumulator[0][1] = testValue - image(row - 3, column + 1);
      minAccumulator[0][2] = testValue - image(row - 2, column + 2);
      minAccumulator[0][3] = testValue - image(row - 1, column + 3);
      minAccumulator[0][4] = testValue - image(row, column + 3);
      minAccumulator[0][5] = testValue - image(row + 1, column + 3);
      minAccumulator[0][6] = testValue - image(row + 2, column + 2);
      minAccumulator[0][7] = testValue - image(row + 3, column + 1);
      minAccumulator[0][8] = testValue - image(row + 3, column);
      minAccumulator[0][9] = testValue - image(row + 3, column - 1);
      minAccumulator[0][10] = testValue - image(row + 2, column - 2);
      minAccumulator[0][11] = testValue - image(row + 1, column - 3);
      minAccumulator[0][12] = testValue - image(row, column - 3);
      minAccumulator[0][13] = testValue - image(row - 1, column - 3);
      minAccumulator[0][14] = testValue - image(row - 2, column - 2);
      minAccumulator[0][15] = testValue - image(row - 3, column - 1);

      // We'll use minAccumulator to examine differences in which
      // testValue is brighter than the neighbor pixel (we'll be
      // accumulating the _minimum_ of the positive differences).
      // We'll use maxAccumulator to examine differences in which
      // testValue is dimmer than the neighbor pixel (we'll be
      // accumulating the least-negative, or _maximum_ negative
      // differences).  It makes the code below a little confusing,
      // but we're going to do this last step by negating the values
      // in maxAccumulator, and then finding their minimum.  It's
      // exactly equivalent...
      for(unsigned int ii = 0; ii < 16; ++ii) {
        maxAccumulator[0][ii] = -(minAccumulator[0][ii]);
        minAccumulator[0][ii] = std::max(minAccumulator[0][ii],
                                         common::Int16(0));
        maxAccumulator[0][ii] = std::max(maxAccumulator[0][ii],
                                         common::Int16(0));
      }
      
      // Second, third, and fourth rows get minimum and maximum values
      // among neighborhoods of size two, four, and eight,
      // respectively.  Note that the use of std::min() instead of
      // std::max() below is intentional (and correct).
      unsigned int step = 1;
      for(unsigned int rr = 1; rr < 4; ++rr) {
        for(unsigned int ii = 0; ii < 16; ++ii) {
          unsigned int jj = (ii + step) % 16;
          minAccumulator[rr][ii] = std::min(minAccumulator[rr - 1][ii],
                                            minAccumulator[rr - 1][jj]);
          maxAccumulator[rr][ii] = std::min(maxAccumulator[rr - 1][ii],
                                            maxAccumulator[rr - 1][jj]);
        }
        step *= 2;
      }

      // Fifth row gets minimum/maximum values among neighborhoods of
      // size 12.
      for(unsigned int ii = 0; ii < 16; ++ii) {
        unsigned int jj = (ii + step) % 16;
        minAccumulator[4][ii] = std::min(minAccumulator[3][ii],
                                         minAccumulator[2][jj]);
        maxAccumulator[4][ii] = std::min(maxAccumulator[3][ii],
                                         maxAccumulator[2][jj]);
      }

      // If threshold were less than the largest of these
      // min-of-12-element-neighborhood values, then this pixel would
      // pass.b
      common::Int16 thresholdFromMin = *std::max_element(
        &(minAccumulator[4][0]), &(minAccumulator[4][0]) + 16);
      common::Int16 thresholdFromMax = *std::max_element(
        &(maxAccumulator[4][0]), &(maxAccumulator[4][0]) + 16);
      common::Int16 threshld = std::max(thresholdFromMin, thresholdFromMax) - 1;
      return std::max(threshld, common::Int16(0));
    }


    bool
    KeypointSelectorFast::
    testPixel(Image<GRAY8> const& image,
              unsigned int row, unsigned int column,
              const common::Int16 threshold,
              KeypointFast& keypoint)
    {
      common::Int16 testValue = image(row, column);

      // From the paper: "A feature is detected at p if the
      // intensities of at least 12 contiguous pixels [on the
      // Bresenham circle of radius 3 surrounding p] are all above or
      // all below the intensity of p by some threshold, t."  We start
      // by evaluating just the four points of the compass, as these
      // let us eliminate the majority of pixels very quickly.
      common::Int16 difference0 = testValue - image(row - 3, column);
      common::Int16 difference1 = testValue - image(row, column + 3);
      common::Int16 difference3 = testValue - image(row, column - 3);
      common::Int16 difference2 = testValue - image(row + 3, column);
      
      // Much more optimizing to do here.  For now, we'll go ahead with a
      // boneheaded approach to get things running.
      int positiveCount = (int(difference0 > threshold)
                           + int(difference1 > threshold)
                           + int(difference2 > threshold)
                           + int(difference3 > threshold));
      if(positiveCount >= 3) {
        return testPixelDetails(image, row, column, testValue, threshold,
                                keypoint, true);
      }

      common::Int16 minusThreshold = -threshold;
      int negativeCount = (int(difference0 < minusThreshold)
                           + int(difference1 < minusThreshold)
                           + int(difference2 < minusThreshold)
                           + int(difference3 < minusThreshold));
      if(negativeCount >= 3) {
        return testPixelDetails(image, row, column, testValue, threshold,
                                keypoint, false);
      }

      // Preliminary test failed.  No need to do the detailed test.
      return false;
    }


    bool
    KeypointSelectorFast::
    testPixelDetails(Image<GRAY8> const& image,
                     unsigned int row, unsigned int column,
                     const common::Int16 testValue,
                     const common::Int16 threshold,
                     KeypointFast& keypoint,
                     bool isPositive)
    {
      // Slow, but clear, implementation for now.
      keypoint.isPositive = isPositive;
      keypoint.row = row;
      keypoint.column = column;
      keypoint.featureVector[0] = image(row - 3, column);
      keypoint.featureVector[1] = image(row - 3, column + 1);
      keypoint.featureVector[2] = image(row - 2, column + 2);
      keypoint.featureVector[3] = image(row - 1, column + 3);
      keypoint.featureVector[4] = image(row, column + 3);
      keypoint.featureVector[5] = image(row + 1, column + 3);
      keypoint.featureVector[6] = image(row + 2, column + 2);
      keypoint.featureVector[7] = image(row + 3, column + 1);
      keypoint.featureVector[8] = image(row + 3, column);
      keypoint.featureVector[9] = image(row + 3, column - 1);
      keypoint.featureVector[10] = image(row + 2, column - 2);
      keypoint.featureVector[11] = image(row + 1, column - 3);
      keypoint.featureVector[12] = image(row, column - 3);
      keypoint.featureVector[13] = image(row - 1, column - 3);
      keypoint.featureVector[14] = image(row - 2, column - 2);
      keypoint.featureVector[15] = image(row - 3, column - 1);

      // Only pass if 12 contiguous neighbors exceed threshold.
      unsigned int passCount = 0;
      if(isPositive) {
        for(unsigned int ii = 0; ii < 16; ++ii) {
          passCount = (((testValue - keypoint.featureVector[ii]) > threshold)
                       ? passCount + 1 : 0);
          if(passCount == 12) {
            return true;
          }
        }
      } else {
        for(unsigned int ii = 0; ii < 16; ++ii) {
          passCount = (((keypoint.featureVector[ii] - testValue) > threshold)
                       ? passCount + 1 : 0);
          if(passCount == 12) {
            return true;
          }
        }
      }
      return false;
    }
    
  } // namespace brick

} // namespace computerVision

    
