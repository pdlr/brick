/**
***************************************************************************
* @file brick/computerVision/keypointSelectorBullseye.cc
*
* Source file defining a class for selecting stable keypoints from an
* image.
*
* Copyright (C) 2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#include <brick/computerVision/keypointSelectorBullseye.hh>

namespace brick {

  namespace computerVision {

    KeypointSelectorBullseye::
    KeypointSelectorBullseye()
      : m_keypointVector(),
        m_threshold(0)
    {
      // Empty.
    }


    std::vector<KeypointBullseye>
    KeypointSelectorBullseye::
    getKeypoints() const
    {
      return m_keypointVector;
    }


    brick::common::Int16
    KeypointSelectorBullseye::
    getThreshold() const
    {
      return m_threshold;
    }


    void
    KeypointSelectorBullseye::
    setImage(Image<GRAY8> const& inImage,
             unsigned int startRow,
             unsigned int startColumn,
             unsigned int stopRow,
             unsigned int stopColumn)
    {
      const unsigned int pixelMeasurementRadius = 8;

      // Discard last image's keypoints.
      m_keypointVector.clear();
      
      // Make sure the passed-in image bounds are legal.
      this->checkAndRepairRegionOfInterest(
        inImage, pixelMeasurementRadius, startRow, startColumn,
        stopRow, stopColumn);

      // We're going to prune most of the image pixels using a
      // threshold based on local symmetry.  Things that aren't
      // symmetrical aren't bullseyes.  Here we estimate what a normal
      // amount of symmetry is for a non-bullseye pixel, so that we
      // can set the threshold higher than that.
      this->m_symmetryThreshold = this->estimateSymmetryThreshold();
      
      // Test every pixel!
      KeypointBullseye keypoint;
      for(unsigned int row = startRow; row < stopRow; ++row) {
        for(unsigned int column = startColumn; column < stopColumn;
            ++column) {
          FloatType symmetry = this->evaluateSymmetry(
            inImage, row, column, keypoint);
          if(symmetry > symmetryThreshold) {
            FloatType bullseyeMetric = this->evaluateBullseyeMetric(
              inImage, row, column, keypoint);
            if(bullseyeMetric > m_targetList[0].bullseyeMetric) {
              this->sorted_insert(keypoint, m_targetList);
            }
          }
        }
      }
    }


    void
    KeypointSelectorBullseye::
    checkAndRepairRegionOfInterest(Image<GRAY8> const& inImage,
                                   unsigned int pixelMeasurementRadius,
                                   unsigned int& startRow,
                                   unsigned int& startColumn,
                                   unsigned int& stopRow,
                                   unsigned int& stopColumn) const
    {
      startRow = std::max(
        startRow, static_cast<unsigned int>(pixelMeasurementRadius));
      stopRow = std::min(
        stopRow,
        static_cast<unsigned int>(inImage.rows() - pixelMeasurementRadius));
      startColumn = std::max(
        startColumn, static_cast<unsigned int>(pixelMeasurementRadius));
      stopColumn = std::min(
        stopColumn,
        static_cast<unsigned int>(inImage.columns() - pixelMeasurementRadius));

      // Of course, all of the above will be broken if there aren't
      // enough rows or columns in the image.
      if(inImage.rows() <= pixelMeasurementRadius) {
        startRow = 0;
        stopRow = 0;
      }
      if(inImage.columns() <= pixelMeasurementRadius) {
        startColumn = 0;
        stopColumn = 0;
      }
    }


    estimateScale(Image<GRAY8> const& image,
                  unsigned int row, unsigned int column,
                  KeypointBullseye& keypoint) const
    {
      // Record four "spokes" of image data in each major direction,
      // while computing the range of pixel values in this
      // neighborhood.
      UnsignedInt8 minimumValue = std::numeric_limits<UnsignedInt8>::max();
      UnsignedInt8 maximumValue = std::numeric_limits<UnsignedInt8>::min();
      for(unsigned int rr = 0; rr < radius; ++rr) {
        keypoint.leftSpoke[rr] = image(row, column - rr);
        keypoint.rightSpoke[rr] = image(row, column + rr);
        minimumValue = std::minimum(
          std::minimum(keypoint.leftSpoke[rr], keypoint.rightSpoke[rr]),
          minimumValue);
        maximumValue = std::maximum(
          std::maximum(keypoint.leftSpoke[rr], keypoint.rightSpoke[rr]),
          maximumValue);
      }

      for(unsigned int rr = 0; rr < radius; ++rr) {
        keypoint.topSpoke[rr] = image(row - rr, column);
        keypoint.bottomSpoke[rr] = image(row + rr, column);
        minimumValue = std::minimum(
          std::minimum(keypoint.topSpoke[rr], keypoint.bottomSpoke[rr]),
          minimumValue);
        maximumValue = std::maximum(
          std::maximum(keypoint.topSpoke[rr], keypoint.bottomSpoke[rr]),
          maximumValue);
      }

      // Compute the midpoint between min and max pixel value.  This
      // is likely to be somewhere in between the black and white
      // values of the bullseye.
      UnsignedInt8 threshold = ((static_cast<UnsignedInt16>(maximumValue)
                                 - static_cast<UnsignedInt16>(minimumValue))
                                >> 1);

      // Now generate a "scale" for the two major directions by
      // counting how far you have to go in each direction to find a
      // total of three transitions between dark and light.  These
      // distances may not be the same because we might be viewing the
      // bullseye obliquely.
      unsigned int count = 0;
      for(unsigned int rr = 1; rr < radius; ++rr) {
        if(((keypoint.topSpoke[rr - 1] < threshold)
            && (keypoint.topSpoke[rr] >= threshold))
           || ((keypoint.topSpoke[rr - 1] >= threshold)
               && (keypoint.topSpoke[rr] < threshold))) {
          ++count;
          if(count >= 3) {
            keypoint.verticalScale = rr;
          }
        }
      }

      count = 0;
      for(unsigned int rr = 1; rr < radius; ++rr) {
        if(((keypoint.leftSpoke[rr - 1] < threshold)
            && (keypoint.leftSpoke[rr] >= threshold))
           || ((keypoint.leftSpoke[rr - 1] >= threshold)
               && (keypoint.leftSpoke[rr] < threshold))) {
          ++count;
          if(count >= 3) {
            keypoint.horizontalScale = rr;
          }
        }
      }
    }


    FloatType
    KeypointSelectorBullseye::
    evaluateSymmetry(Image<GRAY8> const& image,
                     unsigned int row, unsigned int column,
                     KeypointBullseye& keypoint) const
    {
      UnsignedInt32 pixelSum = 0;
      UnsignedInt32 pixelSquaredSum = 0;
      UnsignedInt32 asymmetrySum = 0;
      for(unsigned int rr = 0; rr < radius; ++rr) {
        this->accumulateAsymmetrySums(
          image(row - rr, column - rr), image(row + rr, column + rr),
          pixelSum, pixelSquaredSum, asymmetrySum);
        this->accumulateAsymmetrySums(
          image(row - rr, column), image(row + rr, column),
          pixelSum, pixelSquaredSum, asymmetrySum);
        this->accumulateAsymmetrySums(
          image(row - rr, column + rr), image(row + rr, column - rr),
          pixelSum, pixelSquaredSum, asymmetrySum);
        this->accumulateAsymmetrySums(
          image(row, column - rr), image(row, column + rr),
          pixelSum, pixelSquaredSum, asymmetrySum);
      }

      FloatType count = 8 * radius;
      FloatType pixelMean = pixelSum / count;
      FloatType pixelVariance = pixelSquaredSum / count - pixelMean * pixelMean;
      if(pixelVariance <= 0.0) {
        return 0.0;
      }
      FloatType asymmetry = (asymmetrySum / (count >> 1)) / pixelVariance;
      return asymmetry;
    }
      
  } // namespace brick

} // namespace computerVision
