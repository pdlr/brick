/**
***************************************************************************
* @file keypointMatcherFastTest.cpp
*
* Source file defining tests for the KeypointMatcherFast class.
*
* Copyright (C) 2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/computerVision/keypointMatcherFast.hh>

#include <brick/test/testFixture.hh>

namespace brick {

  namespace computerVision {

    class KeypointMatcherFastTest
      : public brick::test::TestFixture<KeypointMatcherFastTest> {

    public:

      KeypointMatcherFastTest();
      ~KeypointMatcherFastTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testKeypointMatcherFast();
      void testKeypointMatcherFastRotationInvariant();
      void testKeypointMatcherFastRotationInvariant2();
        
    private:

      // Generate test input.
      void
      generateKeypointVectors(std::vector<KeypointFast>& keypoints,
                              std::vector<KeypointFast>& queryPoints);
      
      double m_defaultTolerance;
      
    }; // class KeypointMatcherFastTest


    /* ============== Member Function Definititions ============== */

    KeypointMatcherFastTest::
    KeypointMatcherFastTest()
      : TestFixture<KeypointMatcherFastTest>("KeypointMatcherFastTest"),
        m_defaultTolerance(1.0E-8)
    {
      BRICK_TEST_REGISTER_MEMBER(testKeypointMatcherFast);
      BRICK_TEST_REGISTER_MEMBER(testKeypointMatcherFastRotationInvariant);
      BRICK_TEST_REGISTER_MEMBER(testKeypointMatcherFastRotationInvariant2);
    }


    void
    KeypointMatcherFastTest::
    testKeypointMatcherFast()
    {
      // Carefully crafted array of keypoints.
      std::vector<KeypointFast> keypoints;
      std::vector<KeypointFast> queryPoints;
      this->generateKeypointVectors(keypoints, queryPoints);
      unsigned int numberOfKeypoints = keypoints.size();
      
      // Create and sanity-check a matcher to be tested.
      KeypointMatcherFast matcher;
      KeypointFast matchingPoint;
      BRICK_TEST_ASSERT(matcher.matchKeypoint(queryPoints[0], matchingPoint)
                        == false);

      // Now let the keypoint matcher organize the input points.
      std::random_shuffle(&(keypoints[0]),
                          &(keypoints[0]) + numberOfKeypoints);
      matcher.setKeypoints(&(keypoints[0]),
                           &(keypoints[0]) + numberOfKeypoints);

      // Make sure all keypoints are correctly matched.
      for(unsigned int ii = 0; ii < numberOfKeypoints; ++ii) {
        bool returnValue = matcher.matchKeypoint(
          queryPoints[ii], matchingPoint);
        BRICK_TEST_ASSERT(returnValue == true);
        BRICK_TEST_ASSERT(matchingPoint.row == queryPoints[ii].row);
        BRICK_TEST_ASSERT(matchingPoint.column == queryPoints[ii].column);
        BRICK_TEST_ASSERT(matchingPoint.isPositive
                          == queryPoints[ii].isPositive);
      }
      
    }


    void
    KeypointMatcherFastTest::
    testKeypointMatcherFastRotationInvariant()
    {
      // Carefully crafted array of keypoints.
      std::vector<KeypointFast> keypoints;
      std::vector<KeypointFast> queryPoints;
      this->generateKeypointVectors(keypoints, queryPoints);
      unsigned int numberOfKeypoints = keypoints.size();

      // Tweak queryPoinr 2 so that it matches keypoint 1 better, but only
      // if rotation invariance is up to 3 elements.
      queryPoints[2].featureVector[8] = 20;

      // Create and sanity-check a matcher to be tested.  Allow
      // rotation up to 2 elements.  2/16. = 0.125 rotations ~= 0.78
      // radians.  Note that matcher rounds up to the nearest number
      // of pixels, so if we set this to 0.79, rotation up to 3 pixels
      // will be accepted.
      KeypointMatcherFast matcher(0.78);
      
      KeypointFast matchingPoint;
      BRICK_TEST_ASSERT(matcher.matchKeypoint(queryPoints[0], matchingPoint)
                        == false);

      // Now let the keypoint matcher organize the input points.
      std::random_shuffle(&(keypoints[0]),
                          &(keypoints[0]) + numberOfKeypoints);
      matcher.setKeypoints(&(keypoints[0]),
                           &(keypoints[0]) + numberOfKeypoints);

      // Make sure all keypoints are correctly matched.
      for(unsigned int ii = 0; ii < numberOfKeypoints; ++ii) {
        bool returnValue = matcher.matchKeypoint(
          queryPoints[ii], matchingPoint);
        BRICK_TEST_ASSERT(returnValue == true);
        BRICK_TEST_ASSERT(matchingPoint.isPositive
                          == queryPoints[ii].isPositive);
        if(ii == 2) {
          BRICK_TEST_ASSERT(matchingPoint.row == 3);
          BRICK_TEST_ASSERT(matchingPoint.column == 3);
        } else {
          BRICK_TEST_ASSERT(matchingPoint.row == queryPoints[ii].row);
          BRICK_TEST_ASSERT(matchingPoint.column == queryPoints[ii].column);
        }
      }
      
    }


    void
    KeypointMatcherFastTest::
    testKeypointMatcherFastRotationInvariant2()
    {
      // Carefully crafted array of keypoints.
      std::vector<KeypointFast> keypoints;
      std::vector<KeypointFast> queryPoints;
      this->generateKeypointVectors(keypoints, queryPoints);
      unsigned int numberOfKeypoints = keypoints.size();

      
      // Tweak queryPoinr 2 so that it matches keypoint 1 better, but only
      // if rotation invariance is up to 3 elements.
      queryPoints[2].featureVector[8] = 20;

      // Create and sanity-check a matcher to be tested.  Allow
      // rotation up to 3 elements.  3/16. = 0.185 rotations ~= 1.15
      // radians.
      KeypointMatcherFast matcher(1.15);
      
      KeypointFast matchingPoint;
      BRICK_TEST_ASSERT(matcher.matchKeypoint(queryPoints[0], matchingPoint)
                        == false);

      // Now let the keypoint matcher organize the input points.
      std::random_shuffle(&(keypoints[0]),
                          &(keypoints[0]) + numberOfKeypoints);
      matcher.setKeypoints(&(keypoints[0]),
                           &(keypoints[0]) + numberOfKeypoints);

      // Make sure all keypoints are correctly matched.
      for(unsigned int ii = 0; ii < numberOfKeypoints; ++ii) {
        bool returnValue = matcher.matchKeypoint(
          queryPoints[ii], matchingPoint);
        BRICK_TEST_ASSERT(returnValue == true);
        BRICK_TEST_ASSERT(matchingPoint.isPositive
                          == queryPoints[ii].isPositive);
        if(ii == 2) {
          BRICK_TEST_ASSERT(matchingPoint.row == 1);
          BRICK_TEST_ASSERT(matchingPoint.column == 1);
        } else {
          BRICK_TEST_ASSERT(matchingPoint.row == queryPoints[ii].row);
          BRICK_TEST_ASSERT(matchingPoint.column == queryPoints[ii].column);
        }
      }
    }
    

    // Generate test input.
    void
    KeypointMatcherFastTest::
    generateKeypointVectors(std::vector<KeypointFast>& keypoints,
                            std::vector<KeypointFast>& queryPoints)
    {
      unsigned int const numberOfKeypoints = 20;
      keypoints.resize(20);
      queryPoints.resize(20);

      // Coordinates don't really matter here.  We'll just use them to
      // keep track of which keypoint is which.  Make a set of
      // keypoints, each with a very simple feature vector, such as
      // [0, 0, 0, ...], [1, 1, 1, ...], etc.
      for(unsigned int ii = 0; ii < numberOfKeypoints / 2; ++ii) {
        keypoints[ii].row = ii;
        keypoints[ii].column = ii;
        for(unsigned int jj = 0; jj < 16; ++jj) {
          keypoints[ii].featureVector[jj] =
            static_cast<common::UnsignedInt8>(ii);
        }
        keypoints[ii].isPositive = true;
      }
      
      // Tweak a couple of feature points so their means are not a
      // good indicator of how well they match.
      keypoints[1].featureVector[5] = 20;
      keypoints[10].featureVector[6] = 0;

      // Now make a second set of keypoints, very similar, but with
      // negative center, rather than positive.
      std::copy(&(keypoints[0]), &(keypoints[0]) + numberOfKeypoints / 2,
                &(keypoints[0]) + numberOfKeypoints / 2);
      for(unsigned int ii = numberOfKeypoints / 2; ii < numberOfKeypoints;
          ++ii) {
        keypoints[ii].isPositive = false;
      }
      
      // Make a set of keypoints for which we'll find matches among
      // the original set.  We'll start by simply copying the points
      // we just made, and then perturb the feature vectors so that
      // they don't match perfectly.
      std::copy(&(keypoints[0]), &(keypoints[0]) + numberOfKeypoints,
                &(queryPoints[0]));
      for(unsigned int jj = 0; jj < 16; jj += 4) {
        queryPoints[3].featureVector[jj] += 1;
        queryPoints[5].featureVector[jj] += 1;
        queryPoints[6].featureVector[jj] -= 2;
        queryPoints[13].featureVector[jj] += 1;
        queryPoints[15].featureVector[jj] += 1;
        queryPoints[16].featureVector[jj] -= 2;
      }
      for(unsigned int jj = 3; jj < 16; jj += 4) {
        queryPoints[0].featureVector[jj] += 1; 
        queryPoints[1].featureVector[jj] -= 1;
        queryPoints[10].featureVector[jj] += 2; 
        queryPoints[11].featureVector[jj] -= 1;
      }

    }
    
    
  } // namespace computerVision

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::computerVision::KeypointMatcherFastTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::computerVision::KeypointMatcherFastTest currentTest;

}

#endif
