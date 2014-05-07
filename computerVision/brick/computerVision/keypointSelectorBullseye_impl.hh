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
#include <brick/computerVision/connectedComponents.hh>

// Debugging code.
// #include <brick/utilities/imageIO.hh>

namespace brick {

  namespace computerVision {

    template <class FloatType>
    KeypointSelectorBullseye<FloatType>::
    KeypointSelectorBullseye(brick::common::UInt32 maxNumberOfBullseyes,
                             brick::common::UInt32 maxRadius,
                             brick::common::UInt32 minRadius,
                             brick::common::UInt32 numberOfTransitions,
                             bool isGeneralPositionRequired)
      : m_bullseyePoints(),
        m_bullseyeEdgeCounts(numberOfTransitions),
        m_edgePositions(numberOfTransitions),
        m_isGeneralPositionRequired(isGeneralPositionRequired),
        m_keypointVector(),
        m_keypointGPVector(),
        m_maxNumberOfBullseyes(maxNumberOfBullseyes),
        m_numberOfTransitions(numberOfTransitions),
        m_maxRadius(maxRadius),
        m_minRadius(minRadius),
        m_minDynamicRange(10)
    {
      // Don't crash if the user confuses the min & max arguments.
      if(m_minRadius > m_maxRadius) {
        std::swap(m_minRadius, m_maxRadius);
      }

      // Sanity check number of transition so we have something valid
      // even if the user passes in a zero.
      if(0 == numberOfTransitions) {
        m_numberOfTransitions = 1;
        m_bullseyeEdgeCounts.resize(1);
        m_edgePositions.resize(1);
      }
    }


    template <class FloatType>
    KeypointBullseye<FloatType, FloatType>
    KeypointSelectorBullseye<FloatType>::
    fineTuneKeypoint(
      KeypointBullseye<brick::common::Int32, FloatType> const& inputKeypoint,
      Image<GRAY8> const& inImage)
    {
      // Fill in as much of the return value as we can up-front.
      KeypointBullseye<FloatType, FloatType> result(0, 0);
      result.asymmetry = inputKeypoint.asymmetry;
      result.bullseyeMetric = inputKeypoint.bullseyeMetric;
      result.darkColor = inputKeypoint.darkColor;
      result.lightColor = inputKeypoint.lightColor;
      result.seedPoints = inputKeypoint.seedPoints;
      result.bullseye = inputKeypoint.bullseye;
    
      // Figure out how big a region we need to look at to be sure we
      // get the entire center of the bullseye.
      brick::numeric::Vector2D<FloatType> semimajorAxis;
      if(inputKeypoint.bullseye.getNumberOfRings() > 1) {
        semimajorAxis = inputKeypoint.bullseye.getSemimajorAxis(
          inputKeypoint.bullseye.getNumberOfRings() - 1);
      } else {
        semimajorAxis = inputKeypoint.bullseye.getSemimajorAxis(0);
        semimajorAxis *= FloatType(2);
      }
      int radius = static_cast<int>(
        brick::numeric::magnitude<FloatType>(semimajorAxis) + 0.5);

      // Add a pixel safety margin to prevent the edge of the center
      // of the bullseye from touching the edge of our reginon of
      // interest.
      ++radius;

      // Pick a region of the image that entirely contains the center
      // of the bullseye.
      brick::common::UInt32 startRow = std::max(
        inputKeypoint.row - radius, 0);
      brick::common::UInt32 stopRow = std::min(
        inputKeypoint.row + radius, static_cast<int>(inImage.rows()));
      brick::common::UInt32 startColumn = std::max(
        inputKeypoint.column - radius, 0);
      brick::common::UInt32 stopColumn = std::min(
        inputKeypoint.column + radius, static_cast<int>(inImage.columns()));

      // Create a binary image in which the center of the bullseye is true.
      Image<GRAY1> binaryImage(stopRow - startRow, stopColumn - startColumn);
      brick::common::UInt8 threshold =
        inputKeypoint.darkColor / 2 + inputKeypoint.lightColor / 2;
      for(brick::common::UInt32 rr = 0; rr < binaryImage.rows(); ++rr) {
        for(brick::common::UInt32 cc = 0; cc < binaryImage.columns(); ++cc) {
          bool pixelValue = (inImage(startRow + rr, startColumn + cc)
                             <= threshold);
          binaryImage(rr, cc) = pixelValue;
        }
      }

      // Run connected components on the binarized image so that we
      // can select only the center of the bullseye.
      Image<GRAY8> ccImage = connectedComponents<GRAY8>(binaryImage);

      // Find out which component contains the center of the bullseye.
      brick::common::UInt8 componentNumber = ccImage(
        inputKeypoint.row - startRow, inputKeypoint.column - startColumn);
      if(componentNumber == 0) {
#if 1
        // Careful checking in member function isPlausibleBullseye()
        // should prevent this from ever happening when called from
        // member function setImage().  The exception might be thrown,
        // however, if the user calls fineTuneBullseye() directly.
        BRICK_THROW(brick::common::ValueException,
                    "KeypointSelectorBullseye::fineTuneKeypoint()",
                    "Center of bullseye is the wrong color.");
#else
        // Coloration is wrong, so we can't find the center more
        // precisely.  Punt and just return the integer coords.  This
        // bullseye is probably an outlier anyway...
        result.row = inputKeypoint.row;
        result.column = inputKeypoint.column;
        return result;
#endif
      }
      
      // Compute the centroid of the center of the bullseye.  We have
      // to be careful here, as our coordinate system convention puts
      // 0, 0 at the upper left corner of the the upper-left pixel.  A
      // blob centered on 0, 0 (for example, a blob made up of only
      // the upper-left pixel) would have its centroid in the center
      // of that pixel, at general position coordinates 0.5, 0.5.  For
      // this reason, we add 0.5 to row and column at the end of the
      // centroid calculation.
      brick::common::UInt32 count = 0;
      for(brick::common::UInt32 rr = 0; rr < binaryImage.rows(); ++rr) {
        for(brick::common::UInt32 cc = 0; cc < binaryImage.columns(); ++cc) {
          if(ccImage(rr, cc) == componentNumber) {
            result.row += rr;
            result.column += cc;
            ++count;
          }
        }
      }
      result.row /= count;
      result.column /= count;
      result.row += 0.5;
      result.column += 0.5;
  
      // Translate back into image coordinates.
      result.row += startRow;
      result.column += startColumn;

      // Keypoint position is redundantly encoded.  Update the
      // bullseye instance to reflect our new location.
      result.bullseye.setOrigin(
        brick::numeric::Vector2D<FloatType>(result.column, result.row));
  
      // Note: no checking that the subpixel position is near the
      // integer position.
      return result;
    }

    
    template <class FloatType>
    std::vector< KeypointBullseye<brick::common::Int32, FloatType> >
    KeypointSelectorBullseye<FloatType>::
    getKeypoints() const
    {
      return m_keypointVector;
    }


    template <class FloatType>
    std::vector< KeypointBullseye<FloatType, FloatType> >
    KeypointSelectorBullseye<FloatType>::
    getKeypointsGeneralPosition() const
    {
      return m_keypointGPVector;
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

      for(brick::common::UInt32 ii = 0; ii < 3; ++ii) {
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
             brick::common::UInt32 startRow,
             brick::common::UInt32 startColumn,
             brick::common::UInt32 stopRow,
             brick::common::UInt32 stopColumn)
    {
      // Discard last image's keypoints.
      m_keypointVector.clear();
      m_keypointGPVector.clear();
      
      // Make sure the passed-in image bounds are legal.
      this->checkAndRepairRegionOfInterest(
        inImage.rows(), inImage.columns(), m_maxRadius,
        startRow, startColumn, stopRow, stopColumn);

      // We're going to prune most of the image pixels using a
      // threshold based on local asymmetry.  Things that aren't
      // symmetrical aren't bullseyes.  Here we estimate what a normal
      // amount of asymmetry is for a non-bullseye pixel, so that we
      // can set the threshold higher than that.
      //
      // Figure out how many pixels to sample when estimating.
      brick::common::UInt32 numberOfPixelsToSample =
        (stopRow - startRow) * (stopColumn - startColumn) / 1000;
      numberOfPixelsToSample = std::max(numberOfPixelsToSample,
                                        static_cast<brick::common::UInt32>(100));

      // Do the sampling and estimate the threshold.
      FloatType asymmetryThreshold = this->estimateAsymmetryThreshold(
        inImage, m_minRadius, m_maxRadius, startRow, startColumn,
        stopRow, stopColumn, numberOfPixelsToSample);

      // If a pixel has sufficiently good symmetry, it will be tested
      // with a more expensive bullseye algorithm that needs to know
      // which pixels are edges.  Compute an edge image here.  For now
      // we use the expensive Canny algorithm.
      brick::numeric::Array2D<FloatType> gradientX;
      brick::numeric::Array2D<FloatType> gradientY;
      Image<GRAY1> edgeImage = applyCanny<FloatType>(
        inImage, gradientX, gradientY);

      // Test every pixel!
      brick::common::UInt32 totalPixels = 0;
      brick::common::UInt32 testedPixels = 0;
      for(brick::common::UInt32 row = startRow; row < stopRow; ++row) {
        for(brick::common::UInt32 column = startColumn; column < stopColumn;
            ++column, ++totalPixels) {
          // Create a candidate keypoint.
          KeypointBullseye<brick::common::Int32, FloatType> keypoint(
            row, column);

          // Member function evaluateBullseyeMetric() is too expensive
          // to run at every pixel.  Make absolutely sure this could
          // be a bullseye before proceeding.
          if(!this->isPlausibleBullseye(
               keypoint, inImage, m_minRadius, m_maxRadius,
               asymmetryThreshold)) {
            continue;
          }

          // All prescreening passes.  Go ahead with the expensive
          // bullseye evaluation.
          this->evaluateBullseyeMetric(keypoint, edgeImage,
                                       gradientX, gradientY,
                                       m_minRadius, m_maxRadius);
          this->sortedInsert(keypoint, m_keypointVector,
                             m_maxNumberOfBullseyes);
          ++testedPixels;
        }
      }

      // std::cout << "Tested " << testedPixels
      //           << " (" << (100.0 * testedPixels) / totalPixels << "%) of "
      //           << totalPixels << " pixels." << std::endl;

      // Get general position estimates for our
      if(m_isGeneralPositionRequired) {
        m_keypointGPVector.resize(m_keypointVector.size());
        for(brick::common::UInt32 ii = 0; ii < m_keypointVector.size(); ++ii) {
          m_keypointGPVector[ii] = this->fineTuneKeypoint(
            m_keypointVector[ii], inImage);
        }
      }
    }


    // ============== Private member functions below this line ==============


    template <class FloatType>
    void
    KeypointSelectorBullseye<FloatType>::
    checkAndRepairRegionOfInterest(brick::common::UInt32 rows,
                                   brick::common::UInt32 columns,
                                   brick::common::UInt32 radius,
                                   brick::common::UInt32& startRow,
                                   brick::common::UInt32& startColumn,
                                   brick::common::UInt32& stopRow,
                                   brick::common::UInt32& stopColumn) const
    {
      startRow = std::max(
        startRow, static_cast<brick::common::UInt32>(radius));
      stopRow = std::min(
        stopRow,
        static_cast<brick::common::UInt32>(rows - radius));
      startColumn = std::max(
        startColumn, static_cast<brick::common::UInt32>(radius));
      stopColumn = std::min(
        stopColumn,
        static_cast<brick::common::UInt32>(columns - radius));

      // Of course, all of the above will be broken if there aren't
      // enough rows or columns in the image.
      if(rows <= radius) {
        startRow = 0;
        stopRow = 0;
      }
      if(columns <= radius) {
        startColumn = 0;
        stopColumn = 0;
      }
    }


    template <class FloatType>
    void
    KeypointSelectorBullseye<FloatType>::
    accumulateAsymmetrySums(brick::common::Int32 pixel0,
                            brick::common::Int32 pixel1,
                            brick::common::UInt32& pixelSum,
                            brick::common::UInt32& pixelSquaredSum,
                            brick::common::UInt32& asymmetrySum) const
    {
      pixelSum += pixel0 + pixel1;
      pixelSquaredSum += pixel0 * pixel0 + pixel1 * pixel1;
      brick::common::Int32 difference = pixel1 - pixel0;
      asymmetrySum += difference * difference;
    }


    template <class FloatType>
    bool
    KeypointSelectorBullseye<FloatType>::
    countTransitions(std::vector<brick::common::UInt8> const& spoke,
                     brick::common::UInt32 numberOfTransitions,
                     brick::common::UInt8 minDynamicRange,
                     brick::common::UInt8& darkColor,
                     brick::common::UInt8& lightColor,
                     brick::common::UInt32 minRadius,
                     brick::common::UInt32& actualRadius) const
    {
      brick::common::UInt32 const numberOfColorOutliers = 1;
      brick::common::Int16 const pixelMin =
        std::numeric_limits<brick::common::UInt8>::min();
      brick::common::Int16 const pixelMax =
        std::numeric_limits<brick::common::UInt8>::max();

      // Local variables for tracking the color of the bullseye,
      // ignoring color outliers, and for maintaining a threshold
      // (with hysteresis) between light and dark.
      brick::common::UInt8 darkColorBuffer[numberOfColorOutliers + 1];
      brick::common::UInt8 lightColorBuffer[numberOfColorOutliers + 1];
      brick::common::Int16 threshold = 0;
      bool isDark = true;
      bool isColorUpdated = true;

      // Maintain a count of how many color transitions we've seen, so
      // we'll know when to claim success.
      brick::common::UInt32 transitionCount = 0;
      
      // Initialize our accumulation buffers to best guess light and
      // dark colors so far.
      for(brick::common::UInt32 ii = 0; ii <= numberOfColorOutliers; ++ii) {
        darkColorBuffer[ii] = darkColor;
        lightColorBuffer[ii] = lightColor;
      }
      
      // Look at each value in the input array.
      brick::common::UInt32 ii = 0;
      for(; ii < spoke.size(); ++ii) {
        brick::common::UInt8 currentColor = spoke[ii];
        
        // Update our tracking of light and dark.
        if(currentColor < darkColorBuffer[numberOfColorOutliers]) {
          darkColorBuffer[numberOfColorOutliers] = currentColor;
          for(brick::common::UInt32 jj = numberOfColorOutliers; jj > 0; --jj) {
            if(darkColorBuffer[jj] < darkColorBuffer[jj - 1]) {
              std::swap(darkColorBuffer[jj], darkColorBuffer[jj - 1]);
            } else {
              break;
            }
          }
          isColorUpdated = true;
        }
        if(currentColor > lightColorBuffer[numberOfColorOutliers]) {
          lightColorBuffer[numberOfColorOutliers] = currentColor;
          for(brick::common::UInt32 jj = numberOfColorOutliers; jj > 0; --jj) {
            if(lightColorBuffer[jj] > lightColorBuffer[jj - 1]) {
              std::swap(lightColorBuffer[jj], lightColorBuffer[jj - 1]);
            } else {
              break;
            }
          }
          isColorUpdated = true;
        }

        // If we've changed our idea of what black or white is, then
        // we need to adjust our idea of where the transition between
        // black and white happens.
        if(isColorUpdated) {
          threshold = (darkColorBuffer[numberOfColorOutliers] / 2
                       + lightColorBuffer[numberOfColorOutliers] / 2);

          // Add some hysteresis.
          if(isDark) {
            threshold += minDynamicRange / 2;
          } else {
            threshold -= minDynamicRange / 2;
          }

          // And make sure our threshold isn't out of bounds.
          threshold = std::min(threshold, pixelMax);
          threshold = std::max(threshold, pixelMin);
          isColorUpdated = false;
        }

        // Now see if the next pixel constitutes a transition from
        // light to dark or dark to light, and 
        if(isDark) {
          if(currentColor > threshold) {
            isDark = false;
            isColorUpdated = true; // Force recomputation of threshold.
            ++transitionCount;
          }
        } else {
          if(currentColor < threshold) {
            isDark = true;
            isColorUpdated = true; // Force recomputation of threshold.
            ++transitionCount;
          }
        }

        // If we've found enough transitions, then success!
        if(transitionCount >= numberOfTransitions) {
          darkColor = darkColorBuffer[numberOfColorOutliers];
          lightColor = lightColorBuffer[numberOfColorOutliers];
          actualRadius = ii;
          
          // If the bullseye is too small, indicate this to the
          // calling context.
          if(ii < minRadius) {
            return false;
          }
          return true;
        }
      }

      // If we get here, then we didn't find the requisite number of
      // transitions.
      darkColor = darkColorBuffer[numberOfColorOutliers];
      lightColor = lightColorBuffer[numberOfColorOutliers];
      actualRadius = ii;
      return false;
    }
    

    template <class FloatType>
    bool
    KeypointSelectorBullseye<FloatType>::
    estimateBullseye(
      brick::geometry::Bullseye2D<FloatType>& bullseye,
      std::vector< std::vector< brick::numeric::Vector2D<FloatType> > > const&
        edgePositions,
      brick::common::UInt32 numberOfTransitions) const
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
      const_cast<KeypointSelectorBullseye<FloatType>*>(this)->
        m_bullseyePoints.clear();
      for(brick::common::UInt32 ii = 0; ii < m_numberOfTransitions; ++ii) {
        std::copy(m_edgePositions[ii].begin(), m_edgePositions[ii].end(),
                  std::back_inserter(
                    const_cast<KeypointSelectorBullseye<FloatType>*>(this)->
                    m_bullseyePoints));
        const_cast<KeypointSelectorBullseye<FloatType>*>(this)->
          m_bullseyeEdgeCounts[ii] = m_edgePositions[ii].size();
      }
      
      // We require numberOfTransitions + 2 points because that's what
      // Bullseye2D::estimate() needs.
      // brick::common::UInt32 const numberRequired = numberOfTransitions + 2;
      brick::common::UInt32 const numberRequired = numberOfTransitions + 5;
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
      brick::common::UInt32 numberToRetain = (m_bullseyePoints.size() * 0.75) + 0.5;
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
      brick::common::UInt32 currentRing = 0;
      brick::common::UInt32 pointsThisRing = 0;
      brick::common::UInt32 outputIndex = 0;
      for(brick::common::UInt32 ii = 0; ii < m_bullseyePoints.size(); ++ii) {
        if(pointsThisRing >= m_bullseyeEdgeCounts[currentRing]) {
          ++currentRing;
          pointsThisRing = 0;
        }
        if(absResiduals[ii] > maximumAcceptableResidual) {
          // Found an outlier.  Update bookkeeping and skip it.
          // Again, we fight with our decision to make
          // m_bullseyeEdgeCounts be a class member.
          if(
            --(const_cast<KeypointSelectorBullseye<FloatType>*>(this)->
               m_bullseyeEdgeCounts[currentRing]) == 0
            ) {
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
      try {
        bullseye.estimate(
          inliers.begin(), inliers.end(),
          m_bullseyeEdgeCounts.begin(), m_bullseyeEdgeCounts.end());
      } catch(brick::common::ValueException) {
        // If this call throws, we'll just return the un-updated
        // bullseye.
      }

      // All done.
      return true;
    }

    
#if 0
    template <class FloatType>
    void
    KeypointSelectorBullseye<FloatType>::
    estimateScale(Image<GRAY8> const& image,
                  brick::common::UInt32 radius,
                  brick::common::UInt32 row, brick::common::UInt32 column,
                  KeypointBullseye<brick::common::Int32, FloatType>& keypoint)
      const
    {
      // Record four "spokes" of image data in each major direction,
      // while computing the range of pixel values in this
      // neighborhood.
      brick::common::UInt8 minimumValue =
        std::numeric_limits<brick::common::UInt8>::max();
      brick::common::UInt8 maximumValue =
        std::numeric_limits<brick::common::UInt8>::min();
      for(brick::common::UInt32 rr = 0; rr < radius; ++rr) {
        keypoint.leftSpoke[rr] = image(row, column - rr);
        keypoint.rightSpoke[rr] = image(row, column + rr);
        minimumValue = std::min(
          std::min(keypoint.leftSpoke[rr], keypoint.rightSpoke[rr]),
          minimumValue);
        maximumValue = std::max(
          std::max(keypoint.leftSpoke[rr], keypoint.rightSpoke[rr]),
          maximumValue);
      }

      for(brick::common::UInt32 rr = 0; rr < radius; ++rr) {
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
      brick::common::UInt8 threshold =
        ((static_cast<brick::common::UInt16>(maximumValue)
          - static_cast<brick::common::UInt16>(minimumValue)) >> 1);

      // Now generate a "scale" for the two major directions by
      // counting how far you have to go in each direction to find a
      // total of three transitions between dark and light.  These
      // distances may not be the same because we might be viewing the
      // bullseye obliquely.
      brick::common::UInt32 count = 0;
      for(brick::common::UInt32 rr = 1; rr < radius; ++rr) {
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
      for(brick::common::UInt32 rr = 1; rr < radius; ++rr) {
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
                              brick::common::UInt32 minRadius,
                              brick::common::UInt32 maxRadius,
                              brick::common::UInt32 startRow,
                              brick::common::UInt32 startColumn,
                              brick::common::UInt32 stopRow,
                              brick::common::UInt32 stopColumn,
                              brick::common::UInt32 numberOfSamples) const
    {
      // Figure out how much to subsample rows and columns when
      // computing "normal" for the asymmetry measure.
      brick::common::UInt32 numberOfPixels = 
        (stopRow - startRow) * (stopColumn - startColumn);
      FloatType decimationFactor = (static_cast<FloatType>(numberOfPixels)
                                    / static_cast<FloatType>(numberOfSamples));

      // Of course, if we subsample in both rows and columns, we need
      // the square root of the number we just calculated.
      FloatType axisFactor = brick::common::squareRoot(decimationFactor);
      brick::common::UInt32 step =
        static_cast<brick::common::UInt32>(axisFactor + FloatType(0.5));

      // Now inspect approximately numberOfSamples pixels:
      brick::common::UInt32 count = 0;
      FloatType asymmetrySum = 0;
      FloatType asymmetrySquaredSum = 0;

      for(brick::common::UInt32 row = startRow; row < stopRow; row += step) {
        for(brick::common::UInt32 column = startColumn; column < stopColumn;
            column += step) {
          KeypointBullseye<brick::common::Int32, FloatType> keypoint(
            row, column);
          if(this->isPlausibleBullseye(
               keypoint, inImage, minRadius, maxRadius, 0.0, true)) {
            asymmetrySum += keypoint.asymmetry;
            asymmetrySquaredSum += keypoint.asymmetry * keypoint.asymmetry;
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
      // std::cout << meanAsymmetry << ", " << varianceAsymmetry << ", "
      //           << meanAsymmetry - 2 * varianceAsymmetry << std::endl;
      return meanAsymmetry - 2 * varianceAsymmetry;
      
    }

    
    template <class FloatType>
    void
    KeypointSelectorBullseye<FloatType>::
    evaluateBullseyeMetric(
      KeypointBullseye<brick::common::Int32, FloatType>& keypoint,
      Image<GRAY1> const& edgeImage,
      brick::numeric::Array2D<FloatType> const& gradientX,
      brick::numeric::Array2D<FloatType> const& gradientY,
      brick::common::UInt32 minRadius,
      brick::common::UInt32 maxRadius) const
    {
      // Make sure there's no cruft still left in our pre-allocated
      // buffers.
      for(brick::common::UInt32 ii = 0; ii < m_numberOfTransitions; ++ii) {
        // This pre-allocated working buffer idea sucks because it
        // violates our const-ness promise.
        const_cast<KeypointSelectorBullseye<FloatType>*>(this)->
          m_edgePositions[ii].clear();
      }

      // Clean the incoming keypoint struct.
      keypoint.seedPoints.clear();

      // Find nearby edges along several major directions.  In each
      // direction, put the first-encountered edge in
      // edgePositions[0], the second in edgePositions[1], and so on.

      // Look left.
      brick::common::UInt32 edgeCount = 0;
      for(brick::common::UInt32 ii = 1; ii < maxRadius; ++ii) {
        if(testAndRecordEdges(
             edgeImage, keypoint.row, keypoint.column - ii,
             (const_cast<KeypointSelectorBullseye<FloatType>*>(this)->
              m_edgePositions),
             edgeCount, m_numberOfTransitions)) {
          break;
        }
      }

      // Look right.
      edgeCount = 0;
      for(brick::common::UInt32 ii = 1; ii < maxRadius; ++ii) {
        if(testAndRecordEdges(
             edgeImage, keypoint.row, keypoint.column + ii,
             (const_cast<KeypointSelectorBullseye<FloatType>*>(this)->
              m_edgePositions),
             edgeCount, m_numberOfTransitions)) {
          break;
        }
      }

      // Look up.
      edgeCount = 0;
      for(brick::common::UInt32 ii = 1; ii < maxRadius; ++ii) {
        if(testAndRecordEdges(
             edgeImage, keypoint.row - ii, keypoint.column,
             (const_cast<KeypointSelectorBullseye<FloatType>*>(this)->
              m_edgePositions),
             edgeCount, m_numberOfTransitions)) {
          break;
        }
      }

      // Look up and to the left.
      edgeCount = 0;
      for(brick::common::UInt32 ii = 1; ii < maxRadius; ++ii) {
        if(testAndRecordEdgesDiagonal(
             edgeImage, keypoint.row - ii, keypoint.column - ii, -1, -1,
             (const_cast<KeypointSelectorBullseye<FloatType>*>(this)->
              m_edgePositions),
             edgeCount, m_numberOfTransitions)) {
          break;
        }
      }

      // Look up and to the right.
      edgeCount = 0;
      for(brick::common::UInt32 ii = 1; ii < maxRadius; ++ii) {
        if(testAndRecordEdgesDiagonal(
             edgeImage, keypoint.row - ii, keypoint.column + ii, -1, 1,
             (const_cast<KeypointSelectorBullseye<FloatType>*>(this)->
              m_edgePositions), edgeCount, m_numberOfTransitions)) {
          break;
        }
      }
      
      // Look down.
      edgeCount = 0;
      for(brick::common::UInt32 ii = 1; ii < maxRadius; ++ii) {
        if(testAndRecordEdges(
             edgeImage, keypoint.row + ii, keypoint.column,
             (const_cast<KeypointSelectorBullseye<FloatType>*>(this)->
              m_edgePositions),
             edgeCount, m_numberOfTransitions)) {
          break;
        }
      }


      // Look down and to the left.
      edgeCount = 0;
      for(brick::common::UInt32 ii = 1; ii < maxRadius; ++ii) {
        if(testAndRecordEdgesDiagonal(
             edgeImage, keypoint.row + ii, keypoint.column - ii, 1, -1,
             (const_cast<KeypointSelectorBullseye<FloatType>*>(this)->
              m_edgePositions), edgeCount, m_numberOfTransitions)) {
          break;
        }
      }

      // Look down and to the right.
      edgeCount = 0;
      for(brick::common::UInt32 ii = 1; ii < maxRadius; ++ii) {
        if(testAndRecordEdgesDiagonal(
             edgeImage, keypoint.row + ii, keypoint.column + ii, 1, 1,
             (const_cast<KeypointSelectorBullseye<FloatType>*>(this)->
              m_edgePositions), edgeCount, m_numberOfTransitions)) {
          break;
        }
      }
      
      // If we have a full set of edge points, find the best-fit bullseye.
      keypoint.bullseyeMetric = -1.0;
      brick::geometry::Bullseye2D<FloatType> bullseye;
      if(this->estimateBullseye(
           bullseye, m_edgePositions, m_numberOfTransitions)) {
        FloatType bullseyeMetric = -1.0;
        if(this->validateBullseye(
             bullseye,
             // inImage,
             edgeImage, gradientX, gradientY,
             keypoint.row, keypoint.column,
             minRadius, maxRadius, bullseyeMetric)) {

          // OK, this bullseye passed all the tests, remember it.
          keypoint.bullseyeMetric = bullseyeMetric;
          keypoint.bullseye = bullseye;
          for(brick::common::UInt32 ii = 0; ii < m_numberOfTransitions; ++ii) {
            std::copy(m_edgePositions[ii].begin(), m_edgePositions[ii].end(),
                      std::back_inserter(keypoint.seedPoints));
          }
        }
      }
    }

#undef BRICK_CV_TESTEVALUATE_BREAK

    template <class FloatType>
    bool
    KeypointSelectorBullseye<FloatType>::
    evaluateAsymmetry(Image<GRAY8> const& image,
                     brick::common::UInt32 radius,
                     brick::common::UInt32 row, brick::common::UInt32 column,
                     FloatType& asymmetry) const
    {
      brick::common::UInt32 pixelSum = 0;
      brick::common::UInt32 pixelSquaredSum = 0;
      brick::common::UInt32 asymmetrySum = 0;
      for(brick::common::UInt32 rr = 1; rr < radius; ++rr) {
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
    bool
    KeypointSelectorBullseye<FloatType>::
    isPlausibleBullseye(
      KeypointBullseye<brick::common::Int32, FloatType>& keypoint,
      Image<GRAY8> const& inImage,
      brick::common::UInt32 minRadius,
      brick::common::UInt32 maxRadius,
      FloatType asymmetryThreshold,
      bool forceAsymmetry) const
    {
      // Figure out what light and dark mean in the neighborhood of
      // this potential bullseye.  Meanwhile, make sure that in each
      // of the major directions the image has a dark-light-dark-light
      // pattern, as if it were a bullseye.  If more than one of the
      // major (N,S,E,W) directions doesn't have this pattern, then
      // this is not a bullseye.
      std::vector<brick::common::UInt8> spoke(maxRadius);
      brick::common::UInt32 failedSpokeCount = 0;
      keypoint.darkColor = inImage(keypoint.row, keypoint.column);
      keypoint.lightColor = keypoint.darkColor;

      // While doing this, we'll also estimate how big the
      // hypothetical bullseey is.  This will be useful below, when
      // we'll want to check the asymmetry of the hypothetical
      // bullseye.  Remember the radius in each direction using this
      // array.
      brick::common::UInt32 actualRadii[4];
      
      // Look right.
      for(brick::common::UInt32 ii = 0; ii < maxRadius; ++ii) {
        spoke[ii] = inImage(keypoint.row, keypoint.column + ii);
      }

      // In each paragraph below here, we return false if 2 spokes
      // don't have the required transitions, unless forceAsymmetry is
      // specified, in which case we persever 'til the bitter end.
      if(!this->countTransitions(
           spoke, m_numberOfTransitions, m_minDynamicRange,
           keypoint.darkColor, keypoint.lightColor, minRadius,
           actualRadii[0])) {
        if(++failedSpokeCount >= 2 && !forceAsymmetry) {
          return false;
        }
      }

      // Look left.
      for(brick::common::UInt32 ii = 0; ii < maxRadius; ++ii) {
        spoke[ii] = inImage(keypoint.row, keypoint.column - ii);
      }
      if(!this->countTransitions(
           spoke, m_numberOfTransitions, m_minDynamicRange,
           keypoint.darkColor, keypoint.lightColor, minRadius,
           actualRadii[1])) {
        if(++failedSpokeCount >= 2 && !forceAsymmetry) {
          return false;
        }
      }

      // Look down.
      for(brick::common::UInt32 ii = 0; ii < maxRadius; ++ii) {
        spoke[ii] = inImage(keypoint.row + ii, keypoint.column);
      }
      if(!this->countTransitions(
           spoke, m_numberOfTransitions, m_minDynamicRange,
           keypoint.darkColor, keypoint.lightColor, minRadius,
           actualRadii[2])) {
        if(++failedSpokeCount >= 2 && !forceAsymmetry) {
          return false;
        }
      }

      // Look up.
      for(brick::common::UInt32 ii = 0; ii < maxRadius; ++ii) {
        spoke[ii] = inImage(keypoint.row - ii, keypoint.column);
      }
      if(!this->countTransitions(
           spoke, m_numberOfTransitions, m_minDynamicRange,
           keypoint.darkColor, keypoint.lightColor, minRadius,
           actualRadii[3])) {
        if(++failedSpokeCount >= 2 && !forceAsymmetry) {
          return false;
        }
      }

      // For now, we require bullseyes to be dark in the middle.  Note
      // that this assumption also comes up in countTransitions(),
      // above.
      brick::common::UInt8 threshold =
        keypoint.darkColor / 2 + keypoint.lightColor / 2;
      if(inImage(keypoint.row, keypoint.column) > threshold
         && !forceAsymmetry) {
        return false;
      }

      // We now have 4 different estimates of the size of the
      // hypothetical bullseye. Pick one of the non-extreme ones.
      std::sort(&actualRadii[0], (&actualRadii[0]) + 4);
      brick::common::UInt32 actualRadius = actualRadii[1];

      // Now that we know the approximate size of the hypothetical
      // bullseye, we can see if it's sufficiently symmetrical to be a
      // good candidate.
      if(!this->evaluateAsymmetry(
           inImage, actualRadius, keypoint.row, keypoint.column,
           keypoint.asymmetry)) {
        return false;
      }

      // If we're forcing asymmetry computation, we want to return
      // true here, indicating that the calculation was successful. If
      // we're not forcing, then we want to return true only if
      // asymmetry is below threshold.
      if((keypoint.asymmetry > asymmetryThreshold) && !forceAsymmetry) {
        return false;
      }

      // Failed to discard this candidate.  Return true so it gets
      // subjected to more tests.
      return true;
    }


    template <class FloatType>
    void
    KeypointSelectorBullseye<FloatType>::
    sortedInsert(
      KeypointBullseye<brick::common::Int32, FloatType> const& keypoint,
      std::vector< KeypointBullseye<brick::common::Int32, FloatType> >& keypointVector,
      brick::common::UInt32 maxNumberOfBullseyes)
    {
      // Special case: if keypointVector is empty, just add the new point.
      if(keypointVector.empty()) {
        keypointVector.push_back(keypoint);
        return;
      }

      // Special case: if the new point doesn't make the grade,
      // discard it... we're done with it.
      brick::common::UInt32 vectorSize = keypointVector.size();
      if((vectorSize >= maxNumberOfBullseyes)
         && (keypoint.bullseyeMetric <=
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
      brick::common::UInt32 ii = vectorSize - 1;
      while(ii != 0) {
        if(keypointVector[ii].bullseyeMetric
           <= keypointVector[ii - 1].bullseyeMetric) {
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
                     brick::common::UInt32 row,
                     brick::common::UInt32 column,
                     brick::common::UInt32 minRadius,
                     brick::common::UInt32 maxRadius,
                     FloatType& bullseyeMetric) const
    {
      // If the pixel under consideration isn't at the center of
      // the bullseye, then this isn't the right pixel.
      FloatType differenceInX = bullseye.getOrigin().x() - column;
      FloatType differenceInY = bullseye.getOrigin().y() - row;
      if((brick::common::absoluteValue(differenceInX) > 1)
         || (brick::common::absoluteValue(differenceInY) > 1)) {
        return false;
      }

      // Only proceed if bullseye is smaller than maxRadius and larger
      // than minRadius.
      brick::numeric::Vector2D<FloatType> semimajorAxis =
        bullseye.getSemimajorAxis(m_numberOfTransitions - 1);
      FloatType bullseyeRadiusSquared =
        brick::numeric::dot<FloatType>(semimajorAxis, semimajorAxis);
      if(bullseyeRadiusSquared >= (maxRadius * maxRadius)) {
        return false;
      }
      if(bullseyeRadiusSquared < (minRadius * minRadius)) {
        return false;
      }

      // We expect the bullseye to have its rings spaced roughly evenly.
      std::vector<FloatType> ringWidths(m_numberOfTransitions - 1);
      FloatType averageRingWidth = 0.0;
      for(brick::common::UInt32 ii = 1; ii < m_numberOfTransitions; ++ii) {
        ringWidths[ii - 1] = brick::numeric::magnitude<FloatType>(
          bullseye.getSemimajorAxis(ii) - bullseye.getSemimajorAxis(ii - 1));
        averageRingWidth += ringWidths[ii - 1];
      }
      averageRingWidth /= (m_numberOfTransitions - 1);
      for(brick::common::UInt32 ii = 0; ii < ringWidths.size(); ++ii) {
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
      for(brick::common::UInt32 ii = 0; ii < m_numberOfTransitions; ++ii) {
        majors[ii] = bullseye.getSemimajorAxis(ii);
        minors[ii] = bullseye.getSemiminorAxis(ii);
        majorMagnitudes[ii] = brick::numeric::magnitude<FloatType>(majors[ii]);
        minorMagnitudes[ii] = brick::numeric::magnitude<FloatType>(minors[ii]);
      }
      
      // We compute a measure of how well the bullseye explains the
      // image based on how well the edges in the region match with
      // the rings of the bullseye.
      brick::common::UInt32 kk = m_numberOfTransitions - 1;
      brick::common::UInt32 xRadius =
        brick::common::absoluteValue(bullseye.getSemimajorAxis(kk).x())
        + brick::common::absoluteValue(bullseye.getSemiminorAxis(kk).x());
      brick::common::UInt32 yRadius = 
        brick::common::absoluteValue(bullseye.getSemimajorAxis(kk).y())
        + brick::common::absoluteValue(bullseye.getSemiminorAxis(kk).y());
      xRadius = std::min(xRadius, maxRadius);
      yRadius = std::min(yRadius, maxRadius);
      brick::common::UInt32 minRow = row - yRadius;
      brick::common::UInt32 maxRow = row + yRadius;
      brick::common::UInt32 minColumn = column - xRadius;
      brick::common::UInt32 maxColumn = column + xRadius;
      brick::common::UInt32 inBoundsCount = 0;
      FloatType onRingCount = 0.0;  // FloatType [sic]
      for(brick::common::UInt32 rr = minRow; rr < maxRow; ++rr) {
        for(brick::common::UInt32 cc = minColumn; cc < maxColumn; ++cc) {
          if(edgeImage(rr, cc)) {
            // There's an edge here!  Figure out if it lies on one of
            // the rings of the bullseye.
            brick::numeric::Vector2D<FloatType> vectorToEdge =
              brick::numeric::Vector2D<FloatType>(cc, rr)
              - bullseye.getOrigin();

            // Check the edge against each of the rings of the
            // bullseye to see if it's a fit.
            for(brick::common::UInt32 ii = 0; ii < m_numberOfTransitions; ++ii) {

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
      for(brick::common::UInt32 ii = 0; ii < m_numberOfTransitions; ++ii) {
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
      bullseyeMetric = -1.0;
      if(2 * onRingCount > inBoundsCount) {
        // This number goes from somewhere south of -1 (when
        // onRingCount is zero, and inBoundsCount > edgeLengths), to a
        // max of about 1 (when onRingCount approaches edgeLengths,
        // and inBoundsCount == onRingCount).
        bullseyeMetric = (2 * onRingCount - inBoundsCount) / edgeLengths;
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
      brick::common::UInt32 count = 0;

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
