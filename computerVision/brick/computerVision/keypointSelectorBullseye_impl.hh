/**
***************************************************************************
* @file brick/computerVision/keypointSelectorBullseye_impl.hh
*
* Header file defining a class template for selecting stable keypoints
* from an image.
*
* Copyright (C) 2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_KEYPOINTSELECTORBULLSEYE_IMPL_HH
#define BRICK_COMPUTERVISION_KEYPOINTSELECTORBULLSEYE_IMPL_HH

// This file is included by keypointSelectorBullseye.hh, and should
// not be directly included by user code, so no need to include
// keypointSelectorBullseye.hh here.
// 
// #include <brick/computerVision/keypointSelectorBullseye.hh>

#include <brick/common/mathFunctions.hh>

// xxx
#include <brick/utilities/imageIO.hh>


namespace brick {

  namespace computerVision {

    template <class FloatType>
    KeypointSelectorBullseye<FloatType>::
    KeypointSelectorBullseye(brick::common::UnsignedInt32 maxNumberOfBullseyes,
                             brick::common::UnsignedInt32 maxRadius,
                             brick::common::UnsignedInt32 minRadius)
      : m_keypointVector(),
        m_maxNumberOfBullseyes(maxNumberOfBullseyes),
        m_maxRadius(maxRadius),
        m_minRadius(minRadius)
    {
      // Empty.
    }


    template <class FloatType>
    std::vector< KeypointBullseye<brick::common::Int32> >
    KeypointSelectorBullseye<FloatType>::
    getKeypoints() const
    {
      return m_keypointVector;
    }


    // Process an image to find keypoints.
    template <class FloatType>
    void
    KeypointSelectorBullseye<FloatType>::
    setImage(Image<GRAY8> const& inImage)
    {
      this->setImage(inImage, 0, 0, inImage.rows(), inImage.columns());
    }

    
    template <class FloatType>
    void
    KeypointSelectorBullseye<FloatType>::
    setImage(Image<GRAY8> const& inImage,
             unsigned int startRow,
             unsigned int startColumn,
             unsigned int stopRow,
             unsigned int stopColumn)
    {
      // Discard last image's keypoints.
      m_keypointVector.clear();
      
      // Make sure the passed-in image bounds are legal.
      this->checkAndRepairRegionOfInterest(
        inImage, m_maxRadius, startRow, startColumn, stopRow, stopColumn);

      // We're going to prune most of the image pixels using a
      // threshold based on local symmetry.  Things that aren't
      // symmetrical aren't bullseyes.  Here we estimate what a normal
      // amount of symmetry is for a non-bullseye pixel, so that we
      // can set the threshold higher than that.
      //
      // Figure out how many pixels to sample when estimating.
      unsigned int numberOfPixelsToSample =
        (stopRow - startRow) * (stopColumn - startColumn) / 1000;
      numberOfPixelsToSample = std::max(numberOfPixelsToSample,
                                        static_cast<unsigned int>(100));

      // Do the sampling and estimate the threshold.
      FloatType symmetryThreshold = this->estimateSymmetryThreshold(
        inImage, m_minRadius, startRow, startColumn, stopRow, stopColumn,
        numberOfPixelsToSample);
      
      // Test every pixel!
      // xxx
      Image<GRAY_FLOAT64> outImage(inImage.rows(), inImage.columns());
      outImage = 0.0;
      
      for(unsigned int row = startRow; row < stopRow; ++row) {
        for(unsigned int column = startColumn; column < stopColumn;
            ++column) {
          FloatType symmetry = this->evaluateSymmetry(
            inImage, m_minRadius, row, column);
          if(symmetry < symmetryThreshold) {
            KeypointBullseye<brick::common::Int32> keypoint(
              row, column, symmetry);
            // FloatType bullseyeMetric = this->evaluateBullseyeMetric(
            //   inImage, row, column, keypoint);
            this->sortedInsert(keypoint, m_keypointVector,
                               m_maxNumberOfBullseyes);
            outImage(row, column) = symmetry;
          }
        }
      }

      // xxx
      brick::utilities::writePGM("out.pgm", outImage.data(),
                                 outImage.rows(), outImage.columns(),
                                 true, true, 16);
    }


    // ============== Private member functions below this line ==============


    template <class FloatType>
    void
    KeypointSelectorBullseye<FloatType>::
    checkAndRepairRegionOfInterest(Image<GRAY8> const& inImage,
                                   unsigned int radius,
                                   unsigned int& startRow,
                                   unsigned int& startColumn,
                                   unsigned int& stopRow,
                                   unsigned int& stopColumn) const
    {
      startRow = std::max(
        startRow, static_cast<unsigned int>(radius));
      stopRow = std::min(
        stopRow,
        static_cast<unsigned int>(inImage.rows() - radius));
      startColumn = std::max(
        startColumn, static_cast<unsigned int>(radius));
      stopColumn = std::min(
        stopColumn,
        static_cast<unsigned int>(inImage.columns() - radius));

      // Of course, all of the above will be broken if there aren't
      // enough rows or columns in the image.
      if(inImage.rows() <= radius) {
        startRow = 0;
        stopRow = 0;
      }
      if(inImage.columns() <= radius) {
        startColumn = 0;
        stopColumn = 0;
      }
    }


    template <class FloatType>
    void
    KeypointSelectorBullseye<FloatType>::
    accumulateAsymmetrySums(brick::common::Int32 pixel0,
                            brick::common::Int32 pixel1,
                            brick::common::UnsignedInt32& pixelSum,
                            brick::common::UnsignedInt32& pixelSquaredSum,
                            brick::common::UnsignedInt32& asymmetrySum) const
    {
      pixelSum += pixel0 + pixel1;
      pixelSquaredSum += pixel0 * pixel0 + pixel1 * pixel1;
      brick::common::Int32 difference = pixel1 - pixel0;
      asymmetrySum += difference * difference;
    }
    

    template <class FloatType>
    void
    KeypointSelectorBullseye<FloatType>::
    estimateScale(Image<GRAY8> const& image,
                  unsigned int radius,
                  unsigned int row, unsigned int column,
                  KeypointBullseye<brick::common::Int32>& keypoint) const
    {
      // Record four "spokes" of image data in each major direction,
      // while computing the range of pixel values in this
      // neighborhood.
      brick::common::UnsignedInt8 minimumValue =
        std::numeric_limits<brick::common::UnsignedInt8>::max();
      brick::common::UnsignedInt8 maximumValue =
        std::numeric_limits<brick::common::UnsignedInt8>::min();
      for(unsigned int rr = 0; rr < radius; ++rr) {
        keypoint.leftSpoke[rr] = image(row, column - rr);
        keypoint.rightSpoke[rr] = image(row, column + rr);
        minimumValue = std::min(
          std::min(keypoint.leftSpoke[rr], keypoint.rightSpoke[rr]),
          minimumValue);
        maximumValue = std::max(
          std::max(keypoint.leftSpoke[rr], keypoint.rightSpoke[rr]),
          maximumValue);
      }

      for(unsigned int rr = 0; rr < radius; ++rr) {
        keypoint.topSpoke[rr] = image(row - rr, column);
        keypoint.bottomSpoke[rr] = image(row + rr, column);
        minimumValue = std::min(
          std::min(keypoint.topSpoke[rr], keypoint.bottomSpoke[rr]),
          minimumValue);
        maximumValue = std::max(
          std::max(keypoint.topSpoke[rr], keypoint.bottomSpoke[rr]),
          maximumValue);
      }

      // Compute the midpoint between min and max pixel value.  This
      // is likely to be somewhere in between the black and white
      // values of the bullseye.
      brick::common::UnsignedInt8 threshold =
        ((static_cast<brick::common::UnsignedInt16>(maximumValue)
          - static_cast<brick::common::UnsignedInt16>(minimumValue)) >> 1);

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
          if(count >= this->m_numberOfColorTransitions) {
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
          if(count >= this->m_numberOfColorTransitions) {
            keypoint.horizontalScale = rr;
          }
        }
      }
    }


    // Figure out what "normal" is for the symmetry measure, and
    // pick a threshold that's low enough.
    template <class FloatType>
    FloatType
    KeypointSelectorBullseye<FloatType>::
    estimateSymmetryThreshold(Image<GRAY8> const& inImage,
                              unsigned int radius,
                              unsigned int startRow,
                              unsigned int startColumn,
                              unsigned int stopRow,
                              unsigned int stopColumn,
                              unsigned int numberOfSamples) const
    {
      // Figure out how much to subsample rows and columns when
      // computing "normal" for the symmetry measure.
      unsigned int numberOfPixels = 
        (stopRow - startRow) * (stopColumn - startColumn);
      FloatType decimationFactor = (static_cast<FloatType>(numberOfPixels)
                                    / static_cast<FloatType>(numberOfSamples));

      // Of course, if we subsample in both rows and columns, we need
      // the square root of the number we just calculated.
      FloatType axisFactor = brick::common::squareRoot(decimationFactor);
      unsigned int step =
        static_cast<unsigned int>(axisFactor + FloatType(0.5));

      // Now inspect approximately numberOfSamples pixels:
      unsigned int count = 0;
      FloatType symmetrySum = 0;
      FloatType symmetrySquaredSum = 0;
      
      for(unsigned int row = startRow; row < stopRow; row += step) {
        for(unsigned int column = startColumn; column < stopColumn;
            column += step) {
          FloatType symmetry = this->evaluateSymmetry(
            inImage, radius, row, column);
          symmetrySum += symmetry;
          symmetrySquaredSum += symmetry * symmetry;
          ++count;
        }
      }

      // There is almost no chance that these samples made a Gaussian,
      // but it's easy to pretend, so for now, we just take the
      // 1-sigma bound.

      FloatType meanSymmetry = symmetrySum / count;
      FloatType varianceSymmetry =
        (symmetrySquaredSum / count) - (meanSymmetry * meanSymmetry);

      // Quick check to avoid discarding all possible targets if the
      // numbers go wrong.
      if(meanSymmetry - varianceSymmetry < 0.0) {
        return meanSymmetry;
      }

      // All done.  Return our fake 1-sigma threshold.
      return meanSymmetry - varianceSymmetry;
      
    }


    template <class FloatType>
    FloatType
    KeypointSelectorBullseye<FloatType>::
    evaluateSymmetry(Image<GRAY8> const& image,
                     unsigned int radius,
                     unsigned int row, unsigned int column) const
    {
      brick::common::UnsignedInt32 pixelSum = 0;
      brick::common::UnsignedInt32 pixelSquaredSum = 0;
      brick::common::UnsignedInt32 asymmetrySum = 0;
      for(unsigned int rr = 1; rr < radius; ++rr) {
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

      FloatType count = 8 * (radius - 1);
      FloatType pixelMean = pixelSum / count;
      FloatType pixelVariance = pixelSquaredSum / count - pixelMean * pixelMean;
      if(pixelVariance <= 0.0) {
        return std::numeric_limits<FloatType>::max();
      }
      FloatType asymmetry = (asymmetrySum / (count / 2.0)) / pixelVariance;
      return asymmetry;
    }


    template <class FloatType>
    void
    KeypointSelectorBullseye<FloatType>::
    sortedInsert(
      KeypointBullseye<brick::common::Int32> const& keypoint,
      std::vector< KeypointBullseye<brick::common::Int32> >& keypointVector,
      unsigned int maxNumberOfBullseyes)
    {
      // Special case: if keypointVector is empty, just add the new point.
      if(keypointVector.empty()) {
        keypointVector.push_back(keypoint);
        return;
      }

      // Special case: if the new point doesn't make the grade,
      // discard it... we're done with it.
      unsigned int vectorSize = keypointVector.size();
      if((vectorSize >= maxNumberOfBullseyes)
         && (keypoint.value >= keypointVector[vectorSize - 1].value)) {
        return;
      }

      // At this point, we know that the new point is going to be
      // added to the vector, so make room for it by discarding the
      // worst point we've seen so far (unless keypointVector isn't
      // full yet; in that case, no need to throw away any keypoint.
      while(vectorSize >= maxNumberOfBullseyes) {
        keypointVector.pop_back();
        --vectorSize;
      }

      // Now add the new point, and bubble-sort it into its rightful
      // place.  Inserting a point in the middle of a vector is O(n)
      // anyway, so the bubble sort isn't too much more of a penalty.
      // TBD make this stanza clearer.
      keypointVector.push_back(keypoint);
      ++vectorSize;
      unsigned int ii = vectorSize - 1;
      while(ii != 0) {
        if(keypointVector[ii].value < keypointVector[ii - 1].value) {
          break;
        }
        std::swap(keypointVector[ii], keypointVector[ii - 1]);
        --ii;
      }
    }
    

  } // namespace computerVision
  
} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_KEYPOINTSELECTORBULLSEYE_IMPL_HH */
