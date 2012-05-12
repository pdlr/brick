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

#include <map>
#include <brick/computerVision/keypointSelectorFast.hh>

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
    class KeypointMatcherFast {
    public:

      /** 
       * Default constructor.
       * 
       * @param expectedRotation Set this argument to an angle, in
       * radians, over which you'd like to search (something like 0.5
       * is reasonable).  This is not part of Rosten's algorithm, but
       * we find it useful.  Set to 0.0 if you want traditional,
       * non-rotation-invariant matching (which is faster).  Note that
       * the search is over both positive an negative rotations, so
       * setting expectedRotation to 0.5 will consider rotations of up
       * to 0.5 radians in either direction.
       */
      explicit
      KeypointMatcherFast(double expectedRotation = 0.0);

      
      /** 
       * Destructor.
       */
      virtual
      ~KeypointMatcherFast();


      /** 
       * Search the set of stored keypoints and find the one most
       * similar to the input "query" keypoint.  See member function
       * setKeypoints().
       * 
       * @param query This argument is the keypoint to be matched.
       * 
       * @param bestMatch If a match is found, it is returned via this
       * argument.
       * 
       * @return The return value is true if a match was found, false
       * othrwise.
       */
      bool
      matchKeypoint(KeypointFast const& query, KeypointFast& bestMatch) const;


      /** 
       * Specify the set of keypoints from which to draw matches when
       * member function matchKeypoint() is subsequently called.  All
       * previously set points will be discarded prior to adding the
       * new set.
       * 
       * @param sequenceBegin This argument is the beginning of a
       * sequence of keypoints to be added.
       * 
       * @param sequenceEnd This argument is the end of a sequence of
       * keypoints to be added.
       */
      template <class Iter>
      void
      setKeypoints(Iter sequenceBegin, Iter sequenceEnd);


    protected:

      // Get the mean of the feature vector associated with keypoint.
      double
      computeFeatureVectorMean(KeypointFast const& keypoint) const;


      // Compute the sum of squared difference between the feature
      // vectors of the two input keypoints.
      double
      computeSSD(KeypointFast const& keypoint0, KeypointFast const& keypoint1)
        const;

      // Compute the sum of squared difference between the feature
      // vectors of the two input keypoints.  Considers rotations up
      // to +/- expectedRotation.  This is not part of Rosten's
      // original algorithm.
      double
      computeSSDRotationInvariant(KeypointFast const& keypoint0,
                                  KeypointFast const& keypoint1,
                                  double expectedRotation = 0.5) const;

      // Search the specified map of stored keypoints and find the one
      // most similar to the input "query" keypoint.
      bool
      matchKeypoint(KeypointFast const& query, KeypointFast& bestMatch,
                    std::map<double, KeypointFast> const& keypointMap) const;

      double m_expectedRotation;
      std::map<double, KeypointFast> m_keypointMapNegative;
      std::map<double, KeypointFast> m_keypointMapPositive;
    };
    
  } // namespace computerVision
  
} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/keypointMatcherFast_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_KEYPOINTMATCHERFAST_HH */
