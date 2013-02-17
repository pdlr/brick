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
#include <brick/computerVision/canny.hh>

// Debugging code.
// #include <brick/utilities/imageIO.hh>

namespace brick {

  namespace computerVision {

    template <class FloatType>
    KeypointSelectorBullseye<FloatType>::
    KeypointSelectorBullseye(brick::common::UnsignedInt32 maxNumberOfBullseyes,
                             brick::common::UnsignedInt32 maxRadius,
                             brick::common::UnsignedInt32 minRadius)
      : m_bullseyePoints(),
        m_bullseyeEdgeCounts(3), // TBD: make this user set-able.
        m_edgePositions(3), // TBD: make this user set-able.
        m_keypointVector(),
        m_maxNumberOfBullseyes(maxNumberOfBullseyes),
        m_numberOfTransitions(3), // TBD: make this user set-able.
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


#if 0
    template <class FloatType>
    void drawBullseye(brick::geometry::Bullseye2D<FloatType> const& bullseye,
                      Image<GRAY8>& image,
                      brick::common::UInt8 color)
    {
      brick::numeric::Vector2D<FloatType> origin = bullseye.getOrigin();
      image(origin.y() + 0.5, origin.x() + 0.5) = color;

      for(unsigned int ii = 0; ii < 3; ++ii) {
        for(FloatType angle = 0.0; angle < 6.28; angle += (3.14  / 180.)) {
          brick::numeric::Vector2D<FloatType> pp =
            (origin + std::cos(angle) * bullseye.getSemimajorAxis(ii)
             + std::sin(angle) * bullseye.getSemiminorAxis(ii));
          image(pp.y() + 0.5, pp.x() + 0.5) = color;
        }
      }
    }
#endif
    
    
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
      // threshold based on local asymmetry.  Things that aren't
      // symmetrical aren't bullseyes.  Here we estimate what a normal
      // amount of asymmetry is for a non-bullseye pixel, so that we
      // can set the threshold higher than that.
      //
      // Figure out how many pixels to sample when estimating.
      unsigned int numberOfPixelsToSample =
        (stopRow - startRow) * (stopColumn - startColumn) / 1000;
      numberOfPixelsToSample = std::max(numberOfPixelsToSample,
                                        static_cast<unsigned int>(100));

      // Do the sampling and estimate the threshold.
      FloatType asymmetryThreshold = this->estimateAsymmetryThreshold(
        inImage, m_minRadius, startRow, startColumn, stopRow, stopColumn,
        numberOfPixelsToSample);

      // If a pixel has sufficiently good symmetry, it will be tested
      // with a more expensive bullseye algorithm that needs to know
      // which pixels are edges.  Compute an edge image here.  For now
      // we use the expensive Canny algorithm.
      brick::numeric::Array2D<FloatType> gradientX;
      brick::numeric::Array2D<FloatType> gradientY;
      Image<GRAY1> edgeImage = applyCanny<FloatType>(
        inImage, gradientX, gradientY);

      // Test every pixel!
      unsigned int totalPixels = 0;
      unsigned int testedPixels = 0;
      for(unsigned int row = startRow; row < stopRow; ++row) {
        for(unsigned int column = startColumn; column < stopColumn;
            ++column) {
          FloatType asymmetry = 0;
          if(this->evaluateAsymmetry(
               inImage, m_minRadius, row, column, asymmetry)) {
            if(asymmetry < asymmetryThreshold) {
              KeypointBullseye<brick::common::Int32> keypoint(
                row, column, asymmetry);
              this->evaluateBullseyeMetric(keypoint, edgeImage,
                                           gradientX, gradientY, m_maxRadius);
              this->sortedInsert(keypoint, m_keypointVector,
                                 m_maxNumberOfBullseyes);
              ++testedPixels;
            }
          }
          ++totalPixels;
        }
      }

      std::cout << "Tested " << testedPixels
                << " (" << (100.0 * testedPixels) / totalPixels << "%) of "
                << totalPixels << " pixels." << std::endl;
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
    bool
    KeypointSelectorBullseye<FloatType>::
    estimateBullseye(
      brick::geometry::Bullseye2D<FloatType>& bullseye,
      std::vector< std::vector< brick::numeric::Vector2D<FloatType> > > const&
        edgePositions,
      unsigned int numberOfTransitions)
    {
      // A common failure is to not find any points on the outside
      // ring of the bullseye.  This makes sense: we search from the
      // center, so the outside ring is the one that gets found last.
      // If there are no outside-ring points, then we can't go
      // forward.
      if(edgePositions[numberOfTransitions - 1].size() == 0) {
        return false;
      }

      // m_bullseyePoints is really just here to match the
      // Bullseye2D::estimate() interface.  By making it a member of
      // *this, we avoid reallocating every time.  Probably this is
      // false economy, but we'll profile shortly.  Copy all of the
      // edge points into it.
      m_bullseyePoints.clear();
      for(unsigned int ii = 0; ii < m_numberOfTransitions; ++ii) {
        std::copy(m_edgePositions[ii].begin(), m_edgePositions[ii].end(),
                  std::back_inserter(m_bullseyePoints));
        m_bullseyeEdgeCounts[ii] = m_edgePositions[ii].size();
      }
      
      // We require numberOfTransitions + 2 points because that's what
      // Bullseye2D::estimate() needs.
      // unsigned int const numberRequired = numberOfTransitions + 2;
      unsigned int const numberRequired = numberOfTransitions + 5;
      if(m_bullseyePoints.size() < numberRequired) {
        return false;
      }
      
      // See if the edges look like a bullseye.
      brick::numeric::Array1D<FloatType> residuals(m_bullseyePoints.size());
      try {
        bullseye.estimate(
          m_bullseyePoints.begin(), m_bullseyePoints.end(),
          m_bullseyeEdgeCounts.begin(), m_bullseyeEdgeCounts.end(),
          residuals.begin());
      } catch(brick::common::ValueException) {
        // Input points weren't good enough to define a bullseye.
        return false;
      }

      // Here we do some poor-man's robust statistics.  Assuming that
      // the preponderance of the input points lie on the bullseye,
      // discard the worst 25% of points and hope that the rest are
      // inliers.  Hopefully this gets rid of the occasional bad input
      // point.
      unsigned int numberToRetain = (m_bullseyePoints.size() * 0.75) + 0.5;
      numberToRetain = std::max(numberToRetain, numberRequired);
      if(numberToRetain >= m_bullseyePoints.size()) {
        return true;
      }

      // Now we know how many points to retain, and we know we want
      // the N elements with the smallest residuals.  Find the maximum
      // acceptable residual.  Note that nth_element scrambles
      // absResiduals, so we unscrable afterwords.
      brick::numeric::Array1D<FloatType> absResiduals =
        brick::numeric::abs(residuals);
      std::nth_element(absResiduals.begin(),
                       absResiduals.begin() + numberToRetain - 1,
                       absResiduals.end());
      FloatType maximumAcceptableResidual = absResiduals[numberToRetain - 1];
      absResiduals = brick::numeric::abs(residuals);

      // In some cases, we have more than one point with
      // maximumAcceptableResidual, meaning we can't easily retain
      // exactly numberToRetain samples.  We adjust here to avoid
      // problems later.
      numberToRetain = std::count_if(
        absResiduals.begin(), absResiduals.end(),
        std::bind2nd(std::less_equal<FloatType>(), maximumAcceptableResidual));
      if((numberToRetain >= m_bullseyePoints.size())
         || (numberToRetain < numberRequired)) {
        return true;
      }
      
      // Discard the points with the worst fit.
      std::vector< brick::numeric::Vector2D<FloatType> >
        inliers(numberToRetain);
      unsigned int currentRing = 0;
      unsigned int pointsThisRing = 0;
      unsigned int outputIndex = 0;
      for(unsigned int ii = 0; ii < m_bullseyePoints.size(); ++ii) {
        if(pointsThisRing >= m_bullseyeEdgeCounts[currentRing]) {
          ++currentRing;
          pointsThisRing = 0;
        }
        if(absResiduals[ii] > maximumAcceptableResidual) {
          // Found an outlier.  Update bookkeeping and skip it.
          if(--m_bullseyeEdgeCounts[currentRing] == 0) {
            // We require at least one point in each ring, so we're
            // done.  Fortunately, we still have the bullseye estimate
            // from our non-robust attempt, so return true to indicate
            // that the calling context should use that estimate.
            return true;
          }
          continue;
        }

        // Sanity check to make sure register precision issues don't
        // make us to have more than numberToRetain inliers.
        if(outputIndex >= numberToRetain) {
          BRICK_THROW(brick::common::LogicException,
                      "KeypointSelectorBullseye::estimateBullseye()",
                      "Found too many inliers!");
        }
        
        // Looks like this point is an inlier.  Copy it.
        inliers[outputIndex] = m_bullseyePoints[ii];
        ++outputIndex;
        ++pointsThisRing;
      }

      // Sanity check to make sure register precision issues don't
      // make us have fewer than numberToRetain inliers.
      if(outputIndex < numberToRetain) {
        BRICK_THROW(brick::common::LogicException,
                    "KeypointSelectorBullseye::estimateBullseye()",
                    "Found too few inliers!");
      }

      // Now that we've pruned our set of input points.  Redo the estimation.
      bullseye.estimate(
        inliers.begin(), inliers.end(),
        m_bullseyeEdgeCounts.begin(), m_bullseyeEdgeCounts.end());

      // All done.
      return true;
    }

    
#if 0
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
            break;
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
            break;
          }
        }
      }
    }
#endif /* #if 0 */
    

    // Figure out what "normal" is for the asymmetry measure, and
    // pick a threshold that's low enough.
    template <class FloatType>
    FloatType
    KeypointSelectorBullseye<FloatType>::
    estimateAsymmetryThreshold(Image<GRAY8> const& inImage,
                              unsigned int radius,
                              unsigned int startRow,
                              unsigned int startColumn,
                              unsigned int stopRow,
                              unsigned int stopColumn,
                              unsigned int numberOfSamples) const
    {
      // Figure out how much to subsample rows and columns when
      // computing "normal" for the asymmetry measure.
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
      FloatType asymmetrySum = 0;
      FloatType asymmetrySquaredSum = 0;
      
      for(unsigned int row = startRow; row < stopRow; row += step) {
        for(unsigned int column = startColumn; column < stopColumn;
            column += step) {
          FloatType asymmetry = 0;
          if(this->evaluateAsymmetry(inImage, radius, row, column, asymmetry)) {
            asymmetrySum += asymmetry;
            asymmetrySquaredSum += asymmetry * asymmetry;
            ++count;
          }
        }
      }

      // Protect against cases where we didn't find any valid patches
      // (perhaps a blank image).
      if(count == 0) {
        BRICK_THROW(brick::common::ValueException,
                    "KeypointSelectorBullseye::evaluateAsymmetryThreshold()",
                    "Found no valid patches to sample.");
      }
      
      // There is almost no chance that these samples made a Gaussian,
      // but it's easy to pretend, so for now, we just take the
      // 2-sigma bound.
      
      FloatType meanAsymmetry = asymmetrySum / count;
      FloatType varianceAsymmetry =
        (asymmetrySquaredSum / count) - (meanAsymmetry * meanAsymmetry);

      // Quick check to avoid discarding all possible targets if the
      // numbers go wrong.
      if(meanAsymmetry - 2 * varianceAsymmetry < 0.0) {
        std::cout << "KeypointSelectorBullseye::estimateAsymmetryThreshold() -- "
                  << "Warning: asymmetry model is broken." << std::endl;
        return meanAsymmetry;
      }

      // All done.  Return our fake 1-sigma threshold.
      std::cout << meanAsymmetry << ", " << varianceAsymmetry << ", "
                << meanAsymmetry - 2 * varianceAsymmetry << std::endl;
      return meanAsymmetry - 2 * varianceAsymmetry;
      
    }

    
    template <class FloatType>
    void
    KeypointSelectorBullseye<FloatType>::
    evaluateBullseyeMetric(
        KeypointBullseye<brick::common::Int32>& keypoint,
        Image<GRAY1> const& edgeImage,
        brick::numeric::Array2D<FloatType> const& gradientX,
        brick::numeric::Array2D<FloatType> const& gradientY,
        // unsigned int minRadius,
        unsigned int maxRadius)
    {
      // Make sure there's no cruft still left in our pre-allocated
      // buffers.
      for(unsigned int ii = 0; ii < m_numberOfTransitions; ++ii) {
        m_edgePositions[ii].clear();
      }
      
      // Find nearby edges along several major directions.  In each
      // direction, put the first-encountered edge in
      // edgePositions[0], the second in edgePositions[1], and so on.

      // Look left.
      unsigned int edgeCount = 0;
      for(unsigned int ii = 1; ii < maxRadius; ++ii) {
        if(testAndRecordEdges(
             edgeImage, keypoint.row, keypoint.column - ii, m_edgePositions,
             edgeCount, m_numberOfTransitions)) {
          break;
        }
      }

      // Look right.
      edgeCount = 0;
      for(unsigned int ii = 1; ii < maxRadius; ++ii) {
        if(testAndRecordEdges(
             edgeImage, keypoint.row, keypoint.column + ii, m_edgePositions,
             edgeCount, m_numberOfTransitions)) {
          break;
        }
      }

      // Look up.
      edgeCount = 0;
      for(unsigned int ii = 1; ii < maxRadius; ++ii) {
        if(testAndRecordEdges(
             edgeImage, keypoint.row - ii, keypoint.column, m_edgePositions,
             edgeCount, m_numberOfTransitions)) {
          break;
        }
      }

      // Look up and to the left.
      edgeCount = 0;
      for(unsigned int ii = 1; ii < maxRadius; ++ii) {
        if(testAndRecordEdgesDiagonal(
             edgeImage, keypoint.row - ii, keypoint.column - ii, -1, -1,
             m_edgePositions, edgeCount, m_numberOfTransitions)) {
          break;
        }
      }

      // Look up and to the right.
      edgeCount = 0;
      for(unsigned int ii = 1; ii < maxRadius; ++ii) {
        if(testAndRecordEdgesDiagonal(
             edgeImage, keypoint.row - ii, keypoint.column + ii, -1, 1,
             m_edgePositions, edgeCount, m_numberOfTransitions)) {
          break;
        }
      }
      
      // Look down.
      edgeCount = 0;
      for(unsigned int ii = 1; ii < maxRadius; ++ii) {
        if(testAndRecordEdges(
             edgeImage, keypoint.row + ii, keypoint.column, m_edgePositions,
             edgeCount, m_numberOfTransitions)) {
          break;
        }
      }


      // Look down and to the left.
      edgeCount = 0;
      for(unsigned int ii = 1; ii < maxRadius; ++ii) {
        if(testAndRecordEdgesDiagonal(
             edgeImage, keypoint.row + ii, keypoint.column - ii, 1, -1,
             m_edgePositions, edgeCount, m_numberOfTransitions)) {
          break;
        }
      }

      // Look down and to the right.
      edgeCount = 0;
      for(unsigned int ii = 1; ii < maxRadius; ++ii) {
        if(testAndRecordEdgesDiagonal(
             edgeImage, keypoint.row + ii, keypoint.column + ii, 1, 1,
             m_edgePositions, edgeCount, m_numberOfTransitions)) {
          break;
        }
      }
      
      // If we have a full set of edge points, find the best-fit bullseye.
      keypoint.bullseyeMetric = std::numeric_limits<FloatType>::max();
      brick::geometry::Bullseye2D<FloatType> bullseye;
      if(this->estimateBullseye(
           bullseye, m_edgePositions, m_numberOfTransitions)) {
        FloatType goodness = 0;
        if(this->validateBullseye(
             bullseye,
             // inImage,
             edgeImage, gradientX, gradientY,
             keypoint.row, keypoint.column,
             maxRadius, goodness)) {

          // OK, this bullseye passed all the tests, remember it.
          keypoint.bullseyeMetric = 1.0 / goodness;
          keypoint.bullseye = bullseye;
        }
      }
    }

#undef BRICK_CV_TESTEVALUATE_BREAK

    template <class FloatType>
    bool
    KeypointSelectorBullseye<FloatType>::
    evaluateAsymmetry(Image<GRAY8> const& image,
                     unsigned int radius,
                     unsigned int row, unsigned int column,
                     FloatType& asymmetry) const
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

      // Avoid numerical issues on blank regions of the image, which
      // are never bullseyes.
      if(pixelVariance < 1.0) {
        return false;
      }
      asymmetry = (asymmetrySum / (count / 2.0)) / pixelVariance;
      return true;
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
         && (keypoint.bullseyeMetric >=
             keypointVector[vectorSize - 1].bullseyeMetric)) {
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
        if(keypointVector[ii].bullseyeMetric
           >= keypointVector[ii - 1].bullseyeMetric) {
          break;
        }
        std::swap(keypointVector[ii], keypointVector[ii - 1]);
        --ii;
      }
    }


    template <class FloatType>
    bool
    KeypointSelectorBullseye<FloatType>::
    validateBullseye(brick::geometry::Bullseye2D<FloatType> const& bullseye,
                     // Image<GRAY8> const& inImage,
                     Image<GRAY1> const& edgeImage,
                     brick::numeric::Array2D<FloatType> const& gradientX,
                     brick::numeric::Array2D<FloatType> const& gradientY,
                     unsigned int row,
                     unsigned int column,
                     unsigned int maxRadius,
                     FloatType& goodness)
    {
      // If the pixel under consideration isn't at the center of
      // the bullseye, then this isn't the right pixel.
      FloatType differenceInX = bullseye.getOrigin().x() - column;
      FloatType differenceInY = bullseye.getOrigin().y() - row;
      if((brick::common::absoluteValue(differenceInX) > 1)
         || (brick::common::absoluteValue(differenceInY) > 1)) {
        return false;
      }

      // Only proceed if bullseye is smaller than maxRadius.
      brick::numeric::Vector2D<FloatType> semimajorAxis =
        bullseye.getSemimajorAxis(m_numberOfTransitions - 1);
      FloatType bullseyeRadiusSquared =
        brick::numeric::dot<FloatType>(semimajorAxis, semimajorAxis);
      if(bullseyeRadiusSquared > (maxRadius * maxRadius)) {
        return false;
      }

      // We expect the bullseye to have its rings spaced roughly evenly.
      std::vector<FloatType> ringWidths(m_numberOfTransitions - 1);
      FloatType averageRingWidth = 0.0;
      for(unsigned int ii = 1; ii < m_numberOfTransitions; ++ii) {
        ringWidths[ii - 1] = brick::numeric::magnitude<FloatType>(
          bullseye.getSemimajorAxis(ii) - bullseye.getSemimajorAxis(ii - 1));
        averageRingWidth += ringWidths[ii - 1];
      }
      averageRingWidth /= (m_numberOfTransitions - 1);
      for(unsigned int ii = 0; ii < ringWidths.size(); ++ii) {
        FloatType difference =
          brick::common::absoluteValue(ringWidths[ii] - averageRingWidth);
        if((difference / averageRingWidth) > 0.2) {
          return false;
        }
      }

      // Below, we'll compute a measure of how well the bullseye
      // explains the image.  To do this, we'll need a few
      // pre-computed quantities.
      std::vector< brick::numeric::Vector2D<FloatType> > majors(
        m_numberOfTransitions);
      std::vector< brick::numeric::Vector2D<FloatType> > minors(
        m_numberOfTransitions);
      std::vector<FloatType> majorMagnitudes(m_numberOfTransitions);
      std::vector<FloatType> minorMagnitudes(m_numberOfTransitions);
      for(unsigned int ii = 0; ii < m_numberOfTransitions; ++ii) {
        majors[ii] = bullseye.getSemimajorAxis(ii);
        minors[ii] = bullseye.getSemiminorAxis(ii);
        majorMagnitudes[ii] = brick::numeric::magnitude<FloatType>(majors[ii]);
        minorMagnitudes[ii] = brick::numeric::magnitude<FloatType>(minors[ii]);
      }
      
      // We compute a measure of how well the bullseye explains the
      // image based on how well the edges in the region match with
      // the rings of the bullseye.
      unsigned int kk = m_numberOfTransitions - 1;
      unsigned int xRadius =
        brick::common::absoluteValue(bullseye.getSemimajorAxis(kk).x())
        + brick::common::absoluteValue(bullseye.getSemiminorAxis(kk).x());
      unsigned int yRadius = 
        brick::common::absoluteValue(bullseye.getSemimajorAxis(kk).y())
        + brick::common::absoluteValue(bullseye.getSemiminorAxis(kk).y());
      xRadius = std::min(xRadius, maxRadius);
      yRadius = std::min(yRadius, maxRadius);
      unsigned int minRow = row - yRadius;
      unsigned int maxRow = row + yRadius;
      unsigned int minColumn = column - xRadius;
      unsigned int maxColumn = column + xRadius;
      unsigned int inBoundsCount = 0;
      FloatType onRingCount = 0.0;  // FloatType [sic]
      for(unsigned int rr = minRow; rr < maxRow; ++rr) {
        for(unsigned int cc = minColumn; cc < maxColumn; ++cc) {
          if(edgeImage(rr, cc)) {
            // There's an edge here!  Figure out if it lies on one of
            // the rings of the bullseye.
            brick::numeric::Vector2D<FloatType> vectorToEdge =
              brick::numeric::Vector2D<FloatType>(cc, rr)
              - bullseye.getOrigin();

            // Check the edge against each of the rings of the
            // bullseye to see if it's a fit.
            for(unsigned int ii = 0; ii < m_numberOfTransitions; ++ii) {

              // Points on the ellipse can be parameterized cos(theta)
              // * majorAxis + sin(theta) * minorAxis.  Major and
              // minor axes are orthogonal.  This means that the dot
              // product of a point on the ellipse with the major axis
              // gives cos(theat) * magnitudeSquared(majorAxis).  Here
              // we divide by magnitudeSquared(majorAxis) to get a
              // distance metric that -- for points on the ellipse --
              // is proportional to cos(theta).  Similarly for minor
              // axis & sin(theta).
              FloatType distance0 = (
                brick::numeric::dot<FloatType>(vectorToEdge, majors[ii])
                / (majorMagnitudes[ii] * majorMagnitudes[ii]));
              FloatType distance1 = (
                brick::numeric::dot<FloatType>(vectorToEdge, minors[ii])
                / (minorMagnitudes[ii] * minorMagnitudes[ii]));

              // For points on the ellipse, sin^2 + cos^2 = 1.  For
              // points inside the ellipse this number will be less
              // than 1, for points outside it will be greater than 1.
              FloatType normalizedRadius = 
                brick::common::squareRoot(distance0 * distance0
                                          + distance1 * distance1);

              // We'll deem edge points to be on the ellipse if
              // they're within within 10% of major axis length,
              // assuming that's more than 1 pixel.  If major axis is
              // less than 10 pixels long, then force the threshold to
              // be 1 pixel.
              FloatType tolerance =
                (std::max(FloatType(1.0), FloatType(0.1 * majorMagnitudes[ii]))
                 / majorMagnitudes[ii]);
              if(normalizedRadius < (1.0 - tolerance)) {
                // Edge is inside of this ring, but too far inside to
                // lie on the ring itself.  In any event, it's within
                // the perimeter of the bullseye.
                ++inBoundsCount;
                break;
              }
              if(normalizedRadius < (1.0 + tolerance)) {
                // Looks like this edge lies on the ring, so it's
                // certainly within the perimeter of the bullseye.
                ++inBoundsCount;

                // Before saying this edge is on the ring, make sure
                // its orientation is consistent with that of the
                // ring.  Remember that for points on the ring,
                // distance0 is equal to cos(theta) and distance1 is
                // equal to sin(theta).  Here we just differentiate by
                // theta to find the edge direction.
                brick::numeric::Vector2D<FloatType> nominalEdgeDirection =
                  distance1 * majors[ii] - distance0 * minors[ii];
                nominalEdgeDirection /= brick::numeric::magnitude<FloatType>(
                  nominalEdgeDirection);

                // Actual edge direction can be inferred from image gradients.
                brick::numeric::Vector2D<FloatType> actualEdgeDirection(
                  -gradientY(rr, cc), gradientX(rr, cc));
                actualEdgeDirection /= brick::numeric::magnitude<FloatType>(
                  actualEdgeDirection);

                // Only count this edge pixel to the extent that it's
                // consistent with the expected orientation.
                onRingCount += brick::common::absoluteValue(
                  brick::numeric::dot<FloatType>(
                    actualEdgeDirection, nominalEdgeDirection));
                
                break;
              }
            }
          }
        }
      }

      // Now we know how many edge pixels are on/off the ellipse,
      // compare that with how many we expected.
      FloatType edgeLengths = 0;
      for(unsigned int ii = 0; ii < m_numberOfTransitions; ++ii) {
        FloatType majorLength =
          brick::numeric::magnitude<FloatType>(bullseye.getSemimajorAxis(ii));
        FloatType minorLength =
          brick::numeric::magnitude<FloatType>(bullseye.getSemiminorAxis(ii));
        FloatType approxCircumference =
          FloatType(brick::common::constants::pi) *(
            FloatType(3) * (majorLength + minorLength)
            - brick::common::squareRoot((3 * majorLength + minorLength)
                                        * (majorLength + 3 * minorLength)));
        edgeLengths += approxCircumference;
      }
      goodness = 0.0;
      if(2 * onRingCount > inBoundsCount) {
        goodness = (2 * onRingCount - inBoundsCount) / edgeLengths;
      }
      return true;
    }

#if 0
// Saved code for comparing adjacent bullseye rings to see how much
// they differ in color.  Currently not being used because small
// errors in bullseye estimation (due to quantized pixel position)
// mean that it's hard to stay within the correct ring when sampling.
    {
      // Sample the rings of the bullseye.
      brick::numeric::Vector2D<FloatType> origin = bullseye.getOrigin();
      brick::numeric::Vector2D<FloatType> ring1Major = (
        (bullseye.getSemimajorAxis(0) + bullseye.getSemimajorAxis(1))
        * FloatType(0.5));
      brick::numeric::Vector2D<FloatType> ring1Minor = (
        (bullseye.getSemiminorAxis(0) + bullseye.getSemiminorAxis(1))
        * FloatType(0.5));
      brick::numeric::Vector2D<FloatType> ring2Major = (
        (bullseye.getSemimajorAxis(1) + bullseye.getSemimajorAxis(2))
        * FloatType(0.5));
      brick::numeric::Vector2D<FloatType> ring2Minor = (
        (bullseye.getSemiminorAxis(1) + bullseye.getSemiminorAxis(2))
        * FloatType(0.5));
      FloatType ring1Sum = 0.0;
      FloatType ring1SquaredSum = 0.0;
      FloatType ring2Sum = 0.0;
      FloatType ring2SquaredSum = 0.0;
      unsigned int count = 0;

      for(FloatType theta = 0.0; theta < brick::common::constants::twoPi;
          theta += (brick::common::constants::pi / 11)) {
        FloatType cosineTheta = brick::common::cosine(theta);
        FloatType sineTheta = brick::common::sine(theta);
        brick::numeric::Vector2D<FloatType> ring1Position = (
          origin + cosineTheta * ring1Major + sineTheta * ring1Minor);
        brick::numeric::Vector2D<FloatType> ring2Position = (
          origin + cosineTheta * ring2Major + sineTheta * ring2Minor);
        brick::common::UInt8 ring1Value = inImage(
          brick::common::UInt32(ring1Position.y() + 0.5),
          brick::common::UInt32(ring1Position.x() + 0.5));
        brick::common::UInt8 ring2Value = inImage(
          brick::common::UInt32(ring2Position.y() + 0.5),
          brick::common::UInt32(ring2Position.x() + 0.5));
        ring1Sum += ring1Value;
        ring2Sum += ring2Value;
        ring1SquaredSum += ring1Value * ring1Value;
        ring2SquaredSum += ring2Value * ring2Value;
        ++count;

      FloatType ring1Mean = ring1Sum / count;
      FloatType ring2Mean = ring2Sum / count;
      FloatType ring1Variance = ring1SquaredSum / count - ring1Mean * ring1Mean;
      FloatType ring2Variance = ring2SquaredSum / count - ring2Mean * ring2Mean;

      FloatType maximumVariance = std::max(ring1Variance, ring2Variance);
      goodness = (brick::common::absoluteValue(ring1Mean - ring2Mean)
                  / maximumVariance);
      return true;
    }
#endif

  } // namespace computerVision
  
} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_KEYPOINTSELECTORBULLSEYE_IMPL_HH */
