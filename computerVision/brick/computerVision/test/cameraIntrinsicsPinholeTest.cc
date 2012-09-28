/**
***************************************************************************
* @file cameraIntrinsicsPinholeTest.cc
*
* Source file defining tests for the CameraIntrinsicsPinhole class.
*
* Copyright (C) 2006-2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/common/functional.hh>
#include <brick/computerVision/cameraIntrinsicsPinhole.hh>
#include <brick/test/testFixture.hh>

using namespace brick::common;
using namespace brick::computerVision;
using namespace brick::geometry;
using namespace brick::numeric;
using namespace brick::test;


class CameraIntrinsicsPinholeTest
  : public TestFixture<CameraIntrinsicsPinholeTest> {

public:

  CameraIntrinsicsPinholeTest();
  ~CameraIntrinsicsPinholeTest() {}

  void setUp(const std::string& /* testName */) {}
  void tearDown(const std::string& /* testName */) {}

  // Tests.
  void testAccessorFunctions();
  void testConstructor__void();
  void testConstructor__args();
  void testGetProjectionMatrix();
  void testProject();
  void testReverseProject();
  
private:

  double m_defaultTolerance;
  double m_relaxedTolerance;
  
}; // class CameraIntrinsicsPinholeTest


/* ============== Member Function Definititions ============== */

CameraIntrinsicsPinholeTest::
CameraIntrinsicsPinholeTest()
  : brick::test::TestFixture<CameraIntrinsicsPinholeTest>("CameraIntrinsicsPinholeTest"),
    m_defaultTolerance(1.0E-10),
    m_relaxedTolerance(1.0E-7)
{
  BRICK_TEST_REGISTER_MEMBER(testAccessorFunctions);
  BRICK_TEST_REGISTER_MEMBER(testConstructor__void);
  BRICK_TEST_REGISTER_MEMBER(testConstructor__args);
  BRICK_TEST_REGISTER_MEMBER(testGetProjectionMatrix);
  BRICK_TEST_REGISTER_MEMBER(testProject);
  BRICK_TEST_REGISTER_MEMBER(testReverseProject);
}


void
CameraIntrinsicsPinholeTest::
testAccessorFunctions()
{
  // Arbitrary camera params.
  size_t numPixelsX = 320;
  size_t numPixelsY = 240;
  double focalLength = 0.03;
  double pixelSizeX = 0.001;
  double pixelSizeY = 0.002;
  double centerU = 100.0;
  double centerV = 125.0;
  CameraIntrinsicsPinhole<double> intrinsics(
    numPixelsX, numPixelsY, focalLength, pixelSizeX, pixelSizeY,
    centerU, centerV);

  BRICK_TEST_ASSERT(intrinsics.getNumPixelsX() == numPixelsX);
  BRICK_TEST_ASSERT(intrinsics.getNumPixelsY() == numPixelsY);
  BRICK_TEST_ASSERT(
    approximatelyEqual(intrinsics.getFocalLength(), focalLength,
                       m_defaultTolerance));
  BRICK_TEST_ASSERT(
    approximatelyEqual(intrinsics.getPixelSizeX(), pixelSizeX,
                       m_defaultTolerance));
  BRICK_TEST_ASSERT(
    approximatelyEqual(intrinsics.getPixelSizeY(), pixelSizeY,
                       m_defaultTolerance));
  BRICK_TEST_ASSERT(
    approximatelyEqual(intrinsics.getCenterU(), centerU,
                       m_defaultTolerance));
  BRICK_TEST_ASSERT(
    approximatelyEqual(intrinsics.getCenterV(), centerV,
                       m_defaultTolerance));
}


void
CameraIntrinsicsPinholeTest::
testConstructor__void()
{
  // Pass.
}


void
CameraIntrinsicsPinholeTest::
testConstructor__args()
{
  // Tested in testReverseProject().
}


void
CameraIntrinsicsPinholeTest::
testGetProjectionMatrix()
{
  // Arbitrary camera params.
  size_t numPixelsX = 320;
  size_t numPixelsY = 240;
  double focalLength = 0.03;
  double pixelSizeX = 0.001;
  double pixelSizeY = 0.002;
  double centerU = 100.0;
  double centerV = 125.0;
  CameraIntrinsicsPinhole<double> intrinsics(numPixelsX, numPixelsY,
                                             0.03, 0.001, 0.002, 100, 125);
  Array2D<double> projectionMatrix = intrinsics.getProjectionMatrix();

  for(double zCoord = 1.0; zCoord < 10.0; zCoord += 0.7) {
    for(double yCoord = -1.0; yCoord < 1.0; yCoord += 0.1) {
      for(double xCoord = -1.0; xCoord < 1.0; xCoord += 0.1) {
        Vector3D<double> cameraCoord(xCoord, yCoord, zCoord);
        Vector2D<double> pixelCoord(
          (projectionMatrix(0, 0) * cameraCoord.x()
           + projectionMatrix(0, 1) * cameraCoord.y()
           + projectionMatrix(0, 2) * cameraCoord.z()
           + projectionMatrix(0, 3)),
          (projectionMatrix(1, 0) * cameraCoord.x()
           + projectionMatrix(1, 1) * cameraCoord.y()
           + projectionMatrix(1, 2) * cameraCoord.z()
           + projectionMatrix(1, 3)),
          (projectionMatrix(2, 0) * cameraCoord.x()
           + projectionMatrix(2, 1) * cameraCoord.y()
           + projectionMatrix(2, 2) * cameraCoord.z()
           + projectionMatrix(2, 3)));
        double referenceU =
          centerU + (xCoord * focalLength) / (pixelSizeX * zCoord);
        double referenceV =
          centerV + (yCoord * focalLength) / (pixelSizeY * zCoord);
        BRICK_TEST_ASSERT(
          approximatelyEqual(pixelCoord.x(), referenceU, m_defaultTolerance));
        BRICK_TEST_ASSERT(
          approximatelyEqual(pixelCoord.y(), referenceV, m_defaultTolerance));
      }
    }
  }
}


void
CameraIntrinsicsPinholeTest::
testProject()
{
  // Arbitrary camera params.
  size_t numPixelsX = 320;
  size_t numPixelsY = 240;
  double focalLength = 0.03;
  double pixelSizeX = 0.001;
  double pixelSizeY = 0.002;
  double centerU = 100.0;
  double centerV = 125.0;
  CameraIntrinsicsPinhole<double> intrinsics(
    numPixelsX, numPixelsY, focalLength, pixelSizeX, pixelSizeY,
    centerU, centerV);

  for(double zCoord = 1.0; zCoord < 10.0; zCoord += 0.7) {
    for(double yCoord = -1.0; yCoord < 1.0; yCoord += 0.1) {
      for(double xCoord = -1.0; xCoord < 1.0; xCoord += 0.1) {
        Vector3D<double> cameraCoord(xCoord, yCoord, zCoord);
        Vector2D<double> pixelCoord = intrinsics.project(cameraCoord);
        double referenceU =
          centerU + (xCoord * focalLength) / (pixelSizeX * zCoord);
        double referenceV =
          centerV + (yCoord * focalLength) / (pixelSizeY * zCoord);
        BRICK_TEST_ASSERT(
          approximatelyEqual(pixelCoord.x(), referenceU, m_defaultTolerance));
        BRICK_TEST_ASSERT(
          approximatelyEqual(pixelCoord.y(), referenceV, m_defaultTolerance));
      }
    }
  }
}


void
CameraIntrinsicsPinholeTest::
testReverseProject()
{
  // We test by round trip against convertWorldPointToPixel(), which
  // has its own independent test.

  // Arbitrary camera params.
  size_t numPixelsX = 320;
  size_t numPixelsY = 240;
  CameraIntrinsicsPinhole<double> intrinsics(numPixelsX, numPixelsY,
                                     0.03, 0.001, 0.002, 100, 125);

  for(double vCoord = 0.0; vCoord < numPixelsY; vCoord += 1.2) {
    for(double uCoord = 0.0; uCoord < numPixelsX; uCoord += 1.2) {
      Vector2D<double> pixelCoord(uCoord, vCoord);
      Ray3D<double> ray = intrinsics.reverseProject(pixelCoord);
      Vector3D<double> pointOnRay =
        ray.getOrigin() + 12.0 * ray.getDirectionVector();
      Vector2D<double> recoveredPixelCoord = intrinsics.project(pointOnRay);
      BRICK_TEST_ASSERT(approximatelyEqual(
                        recoveredPixelCoord.x(), pixelCoord.x(),
                        m_defaultTolerance));
      BRICK_TEST_ASSERT(approximatelyEqual(
                        recoveredPixelCoord.y(), pixelCoord.y(),
                        m_defaultTolerance));
    }
  }
}


#if 0

int main(int argc, char** argv)
{
  CameraIntrinsicsPinholeTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  CameraIntrinsicsPinholeTest currentTest;

}

#endif
