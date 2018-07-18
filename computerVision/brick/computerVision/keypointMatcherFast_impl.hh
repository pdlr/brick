/**
***************************************************************************
* @file brick/computerVision/keypointMatcherFast_impl.hh
*
* Header file defining a class template for selecting stable keypoints
* from an image.
*
* Copyright (C) 2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_KEYPOINTMATCHERFAST_IMPL_HH
#define BRICK_COMPUTERVISION_KEYPOINTMATCHERFAST_IMPL_HH

// This file is included by keypointMatcherFast.hh, and should not be directly included
// by user code, so no need to include keypointMatcherFast.hh here.
//
// #include <brick/computerVision/keypointMatcherFast.hh>

namespace brick {

  namespace computerVision {

    template <class Iter>
    void
    KeypointMatcherFast::
    setKeypoints(Iter sequenceBegin, Iter sequenceEnd)
    {
      m_keypointMapNegative.clear();
      m_keypointMapPositive.clear();
      while(sequenceBegin != sequenceEnd) {
        double featureVectorMean = this->computeFeatureVectorMean(
          *sequenceBegin);
        if(sequenceBegin->isPositive) {
          m_keypointMapPositive.insert(
            std::make_pair(featureVectorMean, *sequenceBegin));
        } else {
          m_keypointMapNegative.insert(
            std::make_pair(featureVectorMean, *sequenceBegin));
        }
        ++sequenceBegin;
      }
    }

    // ============== Private member functions below this line ==============

  } // namespace computerVision

} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_KEYPOINTMATCHERFAST_IMPL_HH */
