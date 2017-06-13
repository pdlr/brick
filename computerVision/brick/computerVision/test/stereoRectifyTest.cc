/**
***************************************************************************
* @file brick/computerVision/test/stereoRectifyTest.cc
*
* Source file defining tests for stereoRectify().
*
* Copyright (C) 2009-2017 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/computerVision/stereoRectify.hh>
#include <brick/numeric/rotations.hh>
#include <brick/numeric/utilities.hh>
#include <brick/test/testFixture.hh>

namespace num = brick::numeric;

namespace brick {

  namespace computerVision {
    
    class StereoRectifyTest
      : public brick::test::TestFixture<StereoRectifyTest> {

    public:

      StereoRectifyTest();
      ~StereoRectifyTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testGetReprojectionMatrix();
      void testStereoRectify();

    private:

      const double m_defaultTolerance;
      
    }; // class StereoRectifyTest


    /* ============== Member Function Definititions ============== */

    StereoRectifyTest::
    StereoRectifyTest()
      : brick::test::TestFixture<StereoRectifyTest>("StereoRectifyTest"),
        m_defaultTolerance(1.0E-10)
    {
      BRICK_TEST_REGISTER_MEMBER(testGetReprojectionMatrix);
      BRICK_TEST_REGISTER_MEMBER(testStereoRectify);
    }


    void
    StereoRectifyTest::
    testGetReprojectionMatrix()
    {
      // Set up an arbitrary stereo pair.
      CameraIntrinsicsPinhole<double> intrinsics0(
        640, 480, 300.0, 250.0, 310.0, 220.0);
      CameraIntrinsicsPinhole<double> intrinsics1(
        640, 480, 300.0, 250.0, 325.0, 220.0);
      double baseline = 0.21;  // meters.
      
      // Compute the Q matrix.
      brick::numeric::Transform3D<double> QQ = getReprojectionMatrix(
        intrinsics0, intrinsics1, baseline);

      // Pick some points against which to test.
      std::vector< brick::numeric::Vector3D<double> > testPoints;
      testPoints.push_back(brick::numeric::Vector3D<double>(0.0, 0.0, 1.0));
      testPoints.push_back(brick::numeric::Vector3D<double>(1.0, 0.0, 1.5));
      testPoints.push_back(brick::numeric::Vector3D<double>(-1.0, 1.0, 2.0));
      testPoints.push_back(brick::numeric::Vector3D<double>(0.0, 1.0, 1.5));
      testPoints.push_back(brick::numeric::Vector3D<double>(1.0, -1.0, 3.0));
      testPoints.push_back(brick::numeric::Vector3D<double>(-1.0, 1.0, 1.0));
      testPoints.push_back(brick::numeric::Vector3D<double>(0.0, 1.0, 2.0));
      testPoints.push_back(brick::numeric::Vector3D<double>(1.0, 1.0, 1.5));
      testPoints.push_back(brick::numeric::Vector3D<double>(-1.0, -1.0, 4.0));

      // Check against each point in turn.
      for(size_t ii = 0; ii < testPoints.size(); ++ii) {
        brick::numeric::Vector3D<double> testPoint_c0 = testPoints[ii];
        brick::numeric::Vector3D<double> testPoint_c1(
          testPoint_c0.x() - baseline, testPoint_c0.y(), testPoint_c0.z());

        // Project into both cameras.
        brick::numeric::Vector2D<double> p0 = intrinsics0.project(testPoint_c0);
        brick::numeric::Vector2D<double> p1 = intrinsics1.project(testPoint_c1);

        // Check that our projection assumtions are met.
        BRICK_TEST_ASSERT(
          approximatelyEqual(p0.y(), p1.y(), m_defaultTolerance));

        // Reconstruct the 3D point.
        double disparity = p0.x() - p1.x();
        Vector3D<double> uvd(p0.x(), p0.y(), disparity);
        Vector3D<double> reconstructedPoint = QQ * uvd;

        // Check that reprojected point matches expectations.
        BRICK_TEST_ASSERT(
          approximatelyEqual(reconstructedPoint.x(), testPoint_c0.x(),
                             m_defaultTolerance));
        BRICK_TEST_ASSERT(
          approximatelyEqual(reconstructedPoint.y(), testPoint_c0.y(),
                             m_defaultTolerance));
        BRICK_TEST_ASSERT(
          approximatelyEqual(reconstructedPoint.z(), testPoint_c0.z(),
                             m_defaultTolerance));
      }
    }

    
    void
    StereoRectifyTest::
    testStereoRectify()
    {
      // Set up a pair of (arbitrary) cameras.
      CameraIntrinsicsPinhole<double> intrinsics0(
        640, 480, 500.0, 1.0, 1.2, 300, 220);
      CameraIntrinsicsPinhole<double> intrinsics1(
        640, 480, 500.0, 1.1, 0.9, 350, 260);

      // Pick arbitrary, and slightly different, orientations for the
      // two cameras.
      brick::numeric::Transform3D<double> camera0Tworld =
        brick::numeric::rollPitchYawToTransform3D(
          brick::numeric::Vector3D<double>(0.1, 0.7, -0.3));
      brick::numeric::Transform3D<double> camera1Tworld =
        brick::numeric::rollPitchYawToTransform3D(
          brick::numeric::Vector3D<double>(-0.1, 0.5, -0.1));

      // Make cameras be looking coords origin, from a distance of,
      // say 5m.  Camera0 on the left.
      camera0Tworld.setValue<0, 3>(0.01);
      camera0Tworld.setValue<1, 3>(0.04);
      camera0Tworld.setValue<2, 3>(5.2);
      camera1Tworld.setValue<0, 3>(-0.01);
      camera1Tworld.setValue<1, 3>(-0.06);
      camera1Tworld.setValue<2, 3>(4.9);

      // Invert, since we'll need the inverses later.
      brick::numeric::Transform3D<double> worldTcamera0 =
        camera0Tworld.invert();
      brick::numeric::Transform3D<double> worldTcamera1 =
        camera1Tworld.invert();

      // Do the rectification.
      CameraIntrinsicsPinhole<double> rectifiedIntrinsics0;
      CameraIntrinsicsPinhole<double> rectifiedIntrinsics1;
      brick::numeric::Transform3D<double> rcamera0Tworld;
      brick::numeric::Transform3D<double> rcamera1Tworld;
      brick::numeric::Transform2D<double> image0Trimage0;
      brick::numeric::Transform2D<double> image1Trimage1;
      stereoRectify(intrinsics0, intrinsics1, camera0Tworld, camera1Tworld,
                    rectifiedIntrinsics0, rectifiedIntrinsics1,
                    rcamera0Tworld, rcamera1Tworld,
                    image0Trimage0, image1Trimage1);

      // Check that camera centers are still at the same place.
      brick::numeric::Transform3D<double> worldTrcamera0 = rcamera0Tworld.invert();
      brick::numeric::Transform3D<double> worldTrcamera1 = rcamera1Tworld.invert();
      brick::numeric::Vector3D<double> center0(
        worldTcamera0(0, 3), worldTcamera0(1, 3), worldTcamera0(2, 3));
      brick::numeric::Vector3D<double> center1(
        worldTcamera1(0, 3), worldTcamera1(1, 3), worldTcamera1(2, 3));
      brick::numeric::Vector3D<double> rcenter0(
        worldTrcamera0(0, 3), worldTrcamera0(1, 3), worldTrcamera0(2, 3));
      brick::numeric::Vector3D<double> rcenter1(
        worldTrcamera1(0, 3), worldTrcamera1(1, 3), worldTrcamera1(2, 3));

      BRICK_TEST_ASSERT(
        approximatelyEqual(center0.x(), rcenter0.x(), m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(center0.y(), rcenter0.y(), m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(center0.z(), rcenter0.z(), m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(center1.x(), rcenter1.x(), m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(center1.y(), rcenter1.y(), m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(center1.z(), rcenter1.z(), m_defaultTolerance));


      // Check that the rectified rotation matrices are the same.
      for(size_t row = 0; row < 3; ++row) {
        for(size_t column = 0; column < 3; ++column) {
          BRICK_TEST_ASSERT(
            approximatelyEqual(worldTrcamera0(row, column),
                               worldTrcamera1(row, column),
                               m_defaultTolerance));
        }
      }

      // Check that the rectified rotation matrices are orthonormal.
      brick::numeric::Array2D<double> rArray(3, 3);
      for(size_t row = 0; row < 3; ++row) {
        for(size_t column = 0; column < 3; ++column) {
          rArray(row, column) = worldTrcamera0(row, column);
        }
      }
      brick::numeric::Array2D<double> rrt = brick::numeric::matrixMultiply<double>(
        rArray, rArray.transpose());
      for(size_t row = 0; row < 3; ++row) {
        for(size_t column = 0; column < 3; ++column) {
          if(row == column) {
            BRICK_TEST_ASSERT(
              approximatelyEqual(rrt(row, column), 1.0, m_defaultTolerance));
          } else {
            BRICK_TEST_ASSERT(
              approximatelyEqual(rrt(row, column), 0.0, m_defaultTolerance));
          }
        }
      }

      // Check that camera X axes are parallel to the vector
      // connecting their rotation centers.
      brick::numeric::Vector3D<double> baseline_world = rcenter1 - rcenter0;
      baseline_world /= brick::numeric::magnitude<double>(baseline_world);
      brick::numeric::Vector3D<double> xAxis(
        rcamera1Tworld(0, 0), rcamera1Tworld(0, 1), rcamera1Tworld(0, 2));
      brick::numeric::Vector3D<double> crossProduct = brick::numeric::cross(baseline_world, xAxis);
      double crossMag = brick::numeric::magnitude<double>(crossProduct);
      double dotProduct = brick::numeric::dot<double>(baseline_world, xAxis);
      BRICK_TEST_ASSERT(approximatelyEqual(crossMag, 0.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(approximatelyEqual(dotProduct, 1.0, m_defaultTolerance));
      
      // Pick some points against which to test.
      std::vector< brick::numeric::Vector3D<double> > testPoints;
      testPoints.push_back(brick::numeric::Vector3D<double>(0.0, 0.0, 0.0));
      testPoints.push_back(brick::numeric::Vector3D<double>(1.0, 0.0, 0.0));
      testPoints.push_back(brick::numeric::Vector3D<double>(-1.0, 0.0, 0.0));
      testPoints.push_back(brick::numeric::Vector3D<double>(0.0, 1.0, 0.0));
      testPoints.push_back(brick::numeric::Vector3D<double>(0.0, -1.0, 0.0));
      testPoints.push_back(brick::numeric::Vector3D<double>(0.0, 0.0, 1.0));
      testPoints.push_back(brick::numeric::Vector3D<double>(0.0, 0.0, -1.0));
      testPoints.push_back(brick::numeric::Vector3D<double>(-10.0, 5.0, 2.0));
      testPoints.push_back(brick::numeric::Vector3D<double>(5.0, -7.0, -2.0));

      // Check against each point in turn.
      for(size_t ii = 0; ii < testPoints.size(); ++ii) {
        brick::numeric::Vector3D<double> testPoint = testPoints[ii];

        // Project into all cameras.
        brick::numeric::Vector3D<double> p_camera0 = camera0Tworld * testPoint;
        brick::numeric::Vector3D<double> p_camera1 = camera1Tworld * testPoint;
        brick::numeric::Vector3D<double> p_rcamera0 = rcamera0Tworld * testPoint;
        brick::numeric::Vector3D<double> p_rcamera1 = rcamera1Tworld * testPoint;
        brick::numeric::Vector2D<double> p_image0 = intrinsics0.project(p_camera0);
        brick::numeric::Vector2D<double> p_image1 = intrinsics1.project(p_camera1);
        brick::numeric::Vector2D<double> p_rimage0 =
          rectifiedIntrinsics0.project(p_rcamera0);
        brick::numeric::Vector2D<double> p_rimage1 =
          rectifiedIntrinsics1.project(p_rcamera1);

        // Check that rectified projections are on the same scan line.
        BRICK_TEST_ASSERT(
          approximatelyEqual(p_rimage0.y(), p_rimage1.y(), m_defaultTolerance));

        // Check that mappings from rectified to unrectified images
        // are correct.
        brick::numeric::Vector2D<double> pHat_image0 = image0Trimage0 * p_rimage0;
        brick::numeric::Vector2D<double> pHat_image1 = image1Trimage1 * p_rimage1;
        BRICK_TEST_ASSERT(approximatelyEqual(pHat_image0.x(), p_image0.x(),
                                           m_defaultTolerance));
        BRICK_TEST_ASSERT(approximatelyEqual(pHat_image0.y(), p_image0.y(),
                                           m_defaultTolerance));
        BRICK_TEST_ASSERT(approximatelyEqual(pHat_image1.x(), p_image1.x(),
                                           m_defaultTolerance));
        BRICK_TEST_ASSERT(approximatelyEqual(pHat_image1.y(), p_image1.y(),
                                           m_defaultTolerance));
      }
    }

  } // namespace computerVision

} // namespace brick


#if 1

int main(int /* argc */, char** /* argv */)
{
  brick::computerVision::StereoRectifyTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::computerVision::StereoRectifyTest currentTest;

}

#endif

