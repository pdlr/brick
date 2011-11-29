/**
***************************************************************************
* @file brick/computerVision/keypointMatcherFast.hh
*
* Header file declaring a class template for matching keypoints using
* Rosten's 2005 FAST keypoint matching algorithm.
*
* Copyright (C) 2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_KEYPOINTMATCHERFAST_HH
#define BRICK_COMPUTERVISION_KEYPOINTMATCHERFAST_HH

namespace brick {

  namespace computerVision {

    /**
     ** This class implements Rosten's "FAST" keypoint recognition
     ** algorithm, as as described in [1], with the exception that our
     ** current implementation does not do any transformation on the
     ** feature vector to concentrate energy in the first few
     ** elements, and does not short-circuit SSD computations to speed
     ** up the search.
     **
     ** [1] E. Rosten, and T. Drummond, "Fusing Points and Lines for
     ** High Performance Tracking," International Conference on
     ** Computer Vision, 2005.
     **/
    class KeypointMatcherFast : public CameraIntrinsics {
    public:

      /** 
       * Default constructor.
       */
      KeypointMatcherFast();

      
      /** 
       * Destructor.
       */
      virtual
      ~KeypointMatcherFast() {}


      bool
      matchKeypoint(KeypointFast const& query, KeypointFast& bestMatch);


      template <class Iter>
      void
      setKeypoints(Iter sequenceBegin, Iter sequenceEnd);


    protected:

      double
      computeFeatureVectorMean(KeypointFast const& keypoint);

      double
      computeSSD(KeypointFast const& keypoint0, KeypointFast const& keypoint1);

      
      std::map<double, KeypointFast> m_keypointMap;
    };
    
  } // namespace computerVision
  
} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/keypointSelectorFast_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_KEYPOINTMATCHERFAST_HH */
