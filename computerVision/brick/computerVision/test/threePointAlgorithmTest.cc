/**
***************************************************************************
* @file /brick/computerVision/test/threePointAlgorithmTest.cc
*
* Source file defining tests for threePointAlgorithm().
*
* Copyright (C) 2009,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/computerVision/threePointAlgorithm.hh>
#include <brick/computerVision/imageIO.hh>
#include <brick/linearAlgebra/linearAlgebra.hh>
#include <brick/numeric/vector2D.hh>
#include <brick/random/pseudoRandom.hh>
#include <brick/test/testFixture.hh>

namespace linalg = brick::linearAlgebra;
namespace num = brick::numeric;

namespace brick {

  namespace computerVision {
    
    class ThreePointAlgorithmTest
      : public brick::test::TestFixture<ThreePointAlgorithmTest> {

    public:

      ThreePointAlgorithmTest();
      ~ThreePointAlgorithmTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testThreePointAlgorithm();
      void testThreePointAlgorithmRobust();
      
    private:

      void
      getCameraIntrinsics(
        std::vector< CameraIntrinsicsPinhole<double> >& intrinsicsVector);
      
      void
      getCameraPoses(std::vector< num::Transform3D<double> >& worldTcamVector);

      void
      getTestPoints3D(std::vector< num::Vector3D<double> >& pVector);
      
      bool
      isApproximatelyEqual(const Vector3D<double>& vector0,
                           const Vector3D<double>& vector1);
      
      
      double m_defaultTolerance;
      
    }; // class ThreePointAlgorithmTest


    /* ============== Member Function Definititions ============== */

    ThreePointAlgorithmTest::
    ThreePointAlgorithmTest()
      : brick::test::TestFixture<ThreePointAlgorithmTest>("ThreePointAlgorithmTest"),
        m_defaultTolerance(1.0E-8)
    {
      BRICK_TEST_REGISTER_MEMBER(testThreePointAlgorithm);
      BRICK_TEST_REGISTER_MEMBER(testThreePointAlgorithmRobust);
    }


    void
    ThreePointAlgorithmTest::
    testThreePointAlgorithm()
    {
      std::vector< num::Transform3D<double> > worldTcamVector;
      this->getCameraPoses(worldTcamVector);

      std::vector< CameraIntrinsicsPinhole<double> > intrinsicsVector;
      this->getCameraIntrinsics(intrinsicsVector);

      std::vector< num::Vector3D<double> > worldPoints;
      this->getTestPoints3D(worldPoints);
      
      for(unsigned int poseNumber = 0;
          poseNumber < worldTcamVector.size();
          ++poseNumber) {
        for(unsigned int camNumber = 0;
            camNumber < intrinsicsVector.size();
            ++camNumber) {
          for(unsigned int pointIndex = 0;
              pointIndex < (worldPoints.size() - 3);
              pointIndex += 3) {
            num::Vector3D<double> testPoints_world[3];
            num::Vector3D<double> testPoints_camera[3];
            num::Vector2D<double> testPoints_image[3];

            // Get some world points to test with.
            std::copy(worldPoints.begin() + pointIndex,
                      worldPoints.begin() + pointIndex + 3, testPoints_world);

            // Transform into camera frame.
            num::Transform3D<double> camTworld = worldTcamVector[poseNumber].invert();
            std::transform(&(testPoints_world[0]), &(testPoints_world[0]) + 3,
                           &(testPoints_camera[0]), camTworld.getFunctor());

            // Project into image.
            CameraIntrinsicsPinhole<double> intrinsics = intrinsicsVector[camNumber];
            for(size_t ii = 0; ii < 3; ++ii) {
              testPoints_image[ii] = intrinsics.project(testPoints_camera[ii]);
            }

            double d0 = num::magnitude<double>(
              testPoints_camera[1] - testPoints_camera[2]);
            double d1 = num::magnitude<double>(
              testPoints_camera[0] - testPoints_camera[2]);
            double d2 = num::magnitude<double>(
              testPoints_camera[0] - testPoints_camera[1]);
            
            std::vector< num::Vector3D<double> > camPoints0;
            std::vector< num::Vector3D<double> > camPoints1;
            std::vector< num::Vector3D<double> > camPoints2;
            unsigned int numberOfSolutions = threePointAlgorithm(
              testPoints_world[0], testPoints_world[1], testPoints_world[2],
              testPoints_image[0], testPoints_image[1], testPoints_image[2],
              intrinsics, std::back_inserter(camPoints0),
              std::back_inserter(camPoints1),
              std::back_inserter(camPoints2));

            BRICK_TEST_ASSERT(numberOfSolutions >= 1);
            BRICK_TEST_ASSERT(camPoints0.size() == numberOfSolutions);
            BRICK_TEST_ASSERT(camPoints1.size() == numberOfSolutions);
            BRICK_TEST_ASSERT(camPoints2.size() == numberOfSolutions);

            for(size_t solutionNumber = 0;
                solutionNumber < camPoints0.size();
                ++solutionNumber) {
              num::Vector3D<double> pHat0 = camPoints0[solutionNumber];
              num::Vector3D<double> pHat1 = camPoints1[solutionNumber];
              num::Vector3D<double> pHat2 = camPoints2[solutionNumber];

              num::Vector2D<double> uHat0 = intrinsics.project(pHat0);
              num::Vector2D<double> uHat1 = intrinsics.project(pHat1);
              num::Vector2D<double> uHat2 = intrinsics.project(pHat2);

              double dHat0 = num::magnitude<double>(pHat1 - pHat2);
              double dHat1 = num::magnitude<double>(pHat0 - pHat2);
              double dHat2 = num::magnitude<double>(pHat0 - pHat1);
              
              BRICK_TEST_ASSERT(
                num::magnitude<double>(
                  uHat0 - testPoints_image[0]) < m_defaultTolerance);
              BRICK_TEST_ASSERT(
                num::magnitude<double>(
                  uHat1 - testPoints_image[1]) < m_defaultTolerance);
              BRICK_TEST_ASSERT(
                num::magnitude<double>(
                  uHat2 - testPoints_image[2]) < m_defaultTolerance);

              BRICK_TEST_ASSERT(
                approximatelyEqual(dHat0, d0, m_defaultTolerance));
              BRICK_TEST_ASSERT(
                approximatelyEqual(dHat1, d1, m_defaultTolerance));
              BRICK_TEST_ASSERT(
                approximatelyEqual(dHat2, d2, m_defaultTolerance));
            }
          }
        }
      }
    }


    void
    ThreePointAlgorithmTest::
    testThreePointAlgorithmRobust()
    {
      const unsigned int sampleSize = 7;
      const double residualTolerance = 1.0E-5;
      
      std::vector< num::Transform3D<double> > worldTcamVector;
      this->getCameraPoses(worldTcamVector);

      std::vector< CameraIntrinsicsPinhole<double> > intrinsicsVector;
      this->getCameraIntrinsics(intrinsicsVector);

      std::vector< num::Vector3D<double> > worldPoints;
      this->getTestPoints3D(worldPoints);

      // Create a pseudorandom number generator.
      brick::random::PseudoRandom pRandom(26);

      for(unsigned int poseNumber = 0;
          poseNumber < worldTcamVector.size();
          ++poseNumber) {
        for(unsigned int camNumber = 0;
            camNumber < intrinsicsVector.size();
            ++camNumber) {
          for(unsigned int pointIndex = 0;
              pointIndex < (worldPoints.size() - sampleSize);
              pointIndex += sampleSize) {

            num::Vector3D<double> testPoints_world[sampleSize];
            num::Vector3D<double> testPoints_camera[sampleSize];
            num::Vector2D<double> testPoints_image[sampleSize];

            // Get some world points to test with.
            std::copy(worldPoints.begin() + pointIndex,
                      worldPoints.begin() + pointIndex + sampleSize,
                      testPoints_world);

            // Transform into camera frame.
            num::Transform3D<double> camTworld = worldTcamVector[poseNumber].invert();
            std::transform(&(testPoints_world[0]),
                           &(testPoints_world[0]) + sampleSize,
                           testPoints_camera,
                           camTworld.getFunctor());

            // Project into image.
            CameraIntrinsicsPinhole<double> intrinsics = intrinsicsVector[camNumber];
            for(size_t ii = 0; ii < sampleSize; ++ii) {
              testPoints_image[ii] = intrinsics.project(testPoints_camera[ii]);
            }

            double score;
            num::Transform3D<double> camTworldEstimate = threePointAlgorithmRobust(
              &(testPoints_world[0]), &(testPoints_world[0]) + sampleSize,
              &(testPoints_image[0]), intrinsics, 1, 1.0, score, pRandom);

            BRICK_TEST_ASSERT(score < m_defaultTolerance);

            for(size_t pointNumber = 0;
                pointNumber < sampleSize;
                ++pointNumber) {
              num::Vector3D<double> pHat =
                camTworldEstimate * testPoints_world[pointNumber];
              double residual =
                num::magnitude<double>(pHat - testPoints_camera[pointNumber]);

              BRICK_TEST_ASSERT(residual < residualTolerance);
            }
          }
        }
      }
    }


    void
    ThreePointAlgorithmTest::
    getCameraIntrinsics(std::vector< CameraIntrinsicsPinhole<double> >& intrinsicsVector)
    {
      intrinsicsVector.clear();

      intrinsicsVector.push_back(
        CameraIntrinsicsPinhole<double>(320, 240, 0.006, 0.00006, 0.00005, 180, 100));

      intrinsicsVector.push_back(
        CameraIntrinsicsPinhole<double>(320, 240, 0.004, 0.00002, 0.00001, 150, 150));

      intrinsicsVector.push_back(
        CameraIntrinsicsPinhole<double>(320, 240, 0.006, 0.0001, 0.00012, 180, 100));
    }

    
    void
    ThreePointAlgorithmTest::
    getCameraPoses(std::vector< num::Transform3D<double> >& worldTcamVector)
    {
      worldTcamVector.clear();
      
      std::vector< num::Vector3D<double> > centerVector;
      std::vector<double> angleVector;

      centerVector.push_back(num::Vector3D<double>(-7, 4, -9));
      angleVector.push_back(0.05);

      centerVector.push_back(num::Vector3D<double>(-7, 1, -1));
      angleVector.push_back(0.75);

      centerVector.push_back(num::Vector3D<double>(-3, -2, 10));
      angleVector.push_back(-0.24);

      centerVector.push_back(num::Vector3D<double>(-9, -7, -8));
      angleVector.push_back(0.15);

      centerVector.push_back(num::Vector3D<double>(0, 3, 3));
      angleVector.push_back(-0.2);

      for(size_t ii = 0; ii < centerVector.size(); ++ii) {
        // Set up cameras and target point in world coordinates.
        double cos0 = std::cos(angleVector[ii]);
        double sin0 = std::sin(angleVector[ii]);
        worldTcamVector.push_back(
          num::Transform3D<double>(cos0, 0.0, -sin0, centerVector[ii].x(),
                           0.0, 1.0,   0.0, centerVector[ii].y(),
                           sin0, 0.0,  cos0, centerVector[ii].z(),
                           0.0, 0.0,   0.0, 1.0));
      }
    }
    

    void
    ThreePointAlgorithmTest::
    getTestPoints3D(std::vector< num::Vector3D<double> >& pVector)
    {
      // Create a series of arbitrary points.
      pVector.clear();

      // pVector.push_back(num::Vector2D<double>(2.0, 2.0, 15.0));
      // pVector.push_back(num::Vector2D<double>(1.0, 3.0, 12.0));
      // pVector.push_back(num::Vector2D<double>(0.0, 2.0, 19.0));
      // pVector.push_back(num::Vector2D<double>(-1.0, 4.0, 25.0));
      // pVector.push_back(num::Vector2D<double>(-2.0, 5.0, 12.0));

      pVector.push_back(num::Vector3D<double>(5, 8, 26));
      pVector.push_back(num::Vector3D<double>(-7, -5, 27));
      pVector.push_back(num::Vector3D<double>(1, 2, 11));
      pVector.push_back(num::Vector3D<double>(10, 1, 20));
      pVector.push_back(num::Vector3D<double>(-7, 8, 18));
      pVector.push_back(num::Vector3D<double>(5, 7, 26));
      pVector.push_back(num::Vector3D<double>(-2, 5, 13));
      pVector.push_back(num::Vector3D<double>(-8, -8, 17));
      pVector.push_back(num::Vector3D<double>(-9, 2, 14));
      pVector.push_back(num::Vector3D<double>(6, 5, 14));
      pVector.push_back(num::Vector3D<double>(3, -5, 21));
      pVector.push_back(num::Vector3D<double>(4, 8, 18));
      pVector.push_back(num::Vector3D<double>(7, -3, 29));
      pVector.push_back(num::Vector3D<double>(-9, -1, 22));
      pVector.push_back(num::Vector3D<double>(4, -9, 18));
      pVector.push_back(num::Vector3D<double>(1, -8, 28));
      pVector.push_back(num::Vector3D<double>(-3, -2, 18));
      pVector.push_back(num::Vector3D<double>(8, -6, 19));
      pVector.push_back(num::Vector3D<double>(2.5, 3, 21));
      pVector.push_back(num::Vector3D<double>(1.6, 1, 23));
      pVector.push_back(num::Vector3D<double>(0.2, 5, 23));
      pVector.push_back(num::Vector3D<double>(0, 0, 23));
      pVector.push_back(num::Vector3D<double>(0, 2, 27));
      pVector.push_back(num::Vector3D<double>(6, 4, 15));
      pVector.push_back(num::Vector3D<double>(-4, -4, 17));
      pVector.push_back(num::Vector3D<double>(-3, 3, 22));
      pVector.push_back(num::Vector3D<double>(-1, 5, 25));
      pVector.push_back(num::Vector3D<double>(6, 2, 25));
      pVector.push_back(num::Vector3D<double>(5, 1, 18));
    }
    

    bool
    ThreePointAlgorithmTest::
    isApproximatelyEqual(const Vector3D<double>& vector0,
                         const Vector3D<double>& vector1)
    {
      return
        (approximatelyEqual(vector0.x(), vector1.x(), m_defaultTolerance)
         && approximatelyEqual(vector0.y(), vector1.y(), m_defaultTolerance)
         && approximatelyEqual(vector0.z(), vector1.z(), m_defaultTolerance));
    }
    
  } // namespace computerVision

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::computerVision::ThreePointAlgorithmTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::computerVision::ThreePointAlgorithmTest currentTest;

}

#endif
