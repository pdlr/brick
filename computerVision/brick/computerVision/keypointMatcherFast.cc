/**
***************************************************************************
* @file brick/computerVision/keypointMatcherFast.cc
*
* Source file defining a class template for matching keypoints using
* Rosten's 2005 FAST keypoint matching algorithm.
*
* Copyright (C) 2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#include <numeric>
#include <brick/common/mathFunctions.hh>
#include <brick/computerVision/keypointMatcherFast.hh>

namespace brick {

  namespace computerVision {

    // Default constructor.
    KeypointMatcherFast::
    KeypointMatcherFast()
      : m_keypointMapNegative(),
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


    bool
    KeypointMatcherFast::
    matchKeypoint(KeypointFast const& query, KeypointFast& bestMatch,
                  std::map<double, KeypointFast> const& keypointMap) const
    {
      // Sanity check.
      if(keypointMap.empty()) {
        return false;
      }
      
      // Start by finding the keypoint who's feature vector mean is
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
        double newSSD = this->computeSSD(query, currentIter->second);
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
        double newSSD = this->computeSSD(query, currentIter->second);
        if(newSSD < bestSSDSoFar) {
          bestSSDSoFar = newSSD;
          bestMatch = currentIter->second;
        }
      } 
      return true;
    }

  } // namespace computerVision
  
} // namespace brick
