/**
***************************************************************************
* @file brick/computerVision/keypointMatcherFast.cc
*
* Source file defining a class template for matching keypoints using
* Rosten's 2005 FAST keypoint matching algorithm.
*
* Copyright (C) 2011 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#include <numeric>
#include <brick/common/constants.hh>
#include <brick/common/mathFunctions.hh>
#include <brick/computerVision/keypointMatcherFast.hh>

namespace brick {

  namespace computerVision {

    // Default constructor.
    KeypointMatcherFast::
    KeypointMatcherFast(double expectedRotation)
      : m_expectedRotation(expectedRotation),
        m_keypointMapNegative(),
        m_keypointMapPositive()
    {
      // Empty.
    }


    // Destructor.
    KeypointMatcherFast::
    ~KeypointMatcherFast()
    {
      // Empty.
    }


    bool
    KeypointMatcherFast::
    matchKeypoint(KeypointFast const& query, KeypointFast& bestMatch) const
    {
      if(query.isPositive) {
        return this->matchKeypoint(query, bestMatch,
                                   this->m_keypointMapPositive);
      }
      return this->matchKeypoint(query, bestMatch,
                                 this->m_keypointMapNegative);
    }


    double
    KeypointMatcherFast::
    computeFeatureVectorMean(KeypointFast const& keypoint) const
    {
      return std::accumulate(
        &(keypoint.featureVector[0]),
        &(keypoint.featureVector[0]) + KeypointFast::numberOfFeatures,
        double(0.0)) / KeypointFast::numberOfFeatures;
    }


    double
    KeypointMatcherFast::
    computeSSD(KeypointFast const& keypoint0, KeypointFast const& keypoint1)
      const
    {
      double ssd = 0.0;
      for(unsigned int ii = 0; ii < KeypointFast::numberOfFeatures; ++ii) {
        double difference = (keypoint1.featureVector[ii]
                             - keypoint0.featureVector[ii]);
        ssd += difference * difference;
      }
      return ssd;
    }


    double
    KeypointMatcherFast::
    computeSSDRotationInvariant(KeypointFast const& keypoint0,
                                KeypointFast const& keypoint1,
                                double expectedRotation) const
    {
      // We add 1.0 here, rather than adding 0.5, as in the standard
      // float-to-int conversion, because we want to round up.  This
      // is just like calling ceil().
      const unsigned int angleInSamples = std::min(
        static_cast<unsigned int>(
          std::fabs(expectedRotation) / common::constants::twoPi
          * KeypointFast::numberOfFeatures) + 1.0,
        KeypointFast::numberOfFeatures / 2.0);

      // Inefficient for now...  Compute SSD at each rotation.
      double minimumSsd = std::numeric_limits<double>::max();

      // For "negative" rotations.
      for(unsigned int ii = angleInSamples; ii > 0; --ii) {
        double ssd0 = 0.0;

        // Account for "wrap" at the beginning of the feature vector,
        // where the first few features of keypoint0 match up with the
        // last few features of keypoint1.
        for(unsigned int jj = 0; jj < ii; ++jj) {
          unsigned int wrappedIndex = KeypointFast::numberOfFeatures - ii;
          double difference = (keypoint1.featureVector[wrappedIndex + jj]
                               - keypoint0.featureVector[jj]);
          ssd0 += difference * difference;
        }

        // Account for the rest of the feature vector (the "unwrapped" part).
        for(unsigned int jj = ii; jj < KeypointFast::numberOfFeatures; ++jj) {
          double difference = (keypoint1.featureVector[jj - ii]
                               - keypoint0.featureVector[jj]);
          ssd0 += difference * difference;
        }

        if(ssd0 < minimumSsd) {
          minimumSsd = ssd0;
        }
      }

      // For zero rotation.
      {
        double ssd0 = 0.0;
        for(unsigned int jj = 0; jj < KeypointFast::numberOfFeatures; ++jj) {
          double difference = (keypoint1.featureVector[jj]
                               - keypoint0.featureVector[jj]);
          ssd0 += difference * difference;
        }
        if(ssd0 < minimumSsd) {
          minimumSsd = ssd0;
        }
      }

      // For positive rotation.
      for(unsigned int ii = angleInSamples; ii > 0; --ii) {

        double ssd0 = 0.0;

        // Account for "wrap" at the beginning of the feature vector,
        // where the first few features of keypoint1 match up with the
        // last few features of keypoint0.
        for(unsigned int jj = 0; jj < ii; ++jj) {
          unsigned int wrappedIndex = KeypointFast::numberOfFeatures - ii;
          double difference = (keypoint0.featureVector[wrappedIndex + jj]
                               - keypoint1.featureVector[jj]);
          ssd0 += difference * difference;
        }

        // Account for the rest of the feature vector (the "unwrapped" part).
        for(unsigned int jj = ii; jj < KeypointFast::numberOfFeatures; ++jj) {
          double difference = (keypoint0.featureVector[jj - ii]
                               - keypoint1.featureVector[jj]);
          ssd0 += difference * difference;
        }

        if(ssd0 < minimumSsd) {
          minimumSsd = ssd0;
        }
      }

      return minimumSsd;
    }


    bool
    KeypointMatcherFast::
    matchKeypoint(KeypointFast const& query, KeypointFast& bestMatch,
                  std::map<double, KeypointFast> const& keypointMap) const
    {
      // Sanity check.
      if(keypointMap.empty()) {
        return false;
      }

      // Start by finding the keypoint whose feature vector mean is
      // closest to that of the query point.  This is a good starting
      // point for a linear search.
      typedef std::map<double, KeypointFast>::const_iterator MapIterator;
      double featureVectorMean = this->computeFeatureVectorMean(query);
      MapIterator startIter = keypointMap.lower_bound(featureVectorMean);

      // Now search forward until we know for sure we're not going to
      // find a better match.  We'll know we've gone far enough when
      // the difference in feature vector means is big enough to
      // guarantee that the new SSD between feature vectors is larger
      // than our best so far.
      common::Int64 bestSSDSoFar = std::numeric_limits<common::Int64>::max();
      MapIterator currentIter = startIter;
      while(currentIter != keypointMap.end()) {
        double differenceInMeans = common::absoluteValue(
          featureVectorMean - currentIter->first);
        double boundOnNewSSD = (differenceInMeans * differenceInMeans
                                * KeypointFast::numberOfFeatures);
        if(boundOnNewSSD >= bestSSDSoFar) {
          break;
        }
        double newSSD = this->computeSSDRotationInvariant(
          query, currentIter->second, m_expectedRotation);
        if(newSSD < bestSSDSoFar) {
          bestSSDSoFar = newSSD;
          bestMatch = currentIter->second;
        }
        ++currentIter;
      }

      // Do the search again, but this time go backward, toward the
      // beginning of the map.
      currentIter = startIter;
      while(currentIter != keypointMap.begin()) {
        --currentIter;
        double differenceInMeans = common::absoluteValue(
          featureVectorMean - currentIter->first);
        double boundOnNewSSD = (differenceInMeans * differenceInMeans
                                * KeypointFast::numberOfFeatures);
        if(boundOnNewSSD >= bestSSDSoFar) {
          break;
        }
        double newSSD = this->computeSSDRotationInvariant(
          query, currentIter->second, m_expectedRotation);
        if(newSSD < bestSSDSoFar) {
          bestSSDSoFar = newSSD;
          bestMatch = currentIter->second;
        }
      }
      return true;
    }

  } // namespace computerVision

} // namespace brick
