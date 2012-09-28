/**
***************************************************************************
* @file calibrationToolsRobustTest.cpp
*
* Source file defining tests for the calibration utility routines.
*
* Copyright (C) 2009-2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/


#include <brick/common/functional.hh>
#include <brick/computerVision/calibrationToolsRobust.hh>
#include <brick/computerVision/cameraIntrinsicsPinhole.hh>
#include <brick/computerVision/cameraIntrinsicsPlumbBob.hh>
#include <brick/numeric/transform3D.hh>
#include <brick/numeric/quaternion.hh>
#include <brick/numeric/rotations.hh>
#include <brick/test/testFixture.hh>

using namespace brick::common;
using namespace brick::computerVision;
using namespace brick::geometry;
using namespace brick::numeric;
using namespace brick::test;


class CalibrationToolsRobustTest
  : public TestFixture<CalibrationToolsRobustTest> {

public:

  CalibrationToolsRobustTest();
  ~CalibrationToolsRobustTest() {}

  void setUp(const std::string& /* testName */) {}
  void tearDown(const std::string& /* testName */) {}

  // Tests.
  void testEstimateCameraParametersRobust();
  void testEstimateCameraParametersPinholeRobust();
  
private:

  bool
  checkIntrinsicsEqual(CameraIntrinsicsPinhole<double> const& intrinsics0,
                       CameraIntrinsicsPinhole<double> const& intrinsics1,
                       double tolerance);

  bool
  checkIntrinsicsEqual(CameraIntrinsicsPlumbBob<double> const& intrinsics0,
                       CameraIntrinsicsPlumbBob<double> const& intrinsics1,
                       double tolerance);
  

  bool
  checkTransformEqual(Transform3D<double> const& transform0,
                      Transform3D<double> const& transform1,
                      double tolerance);

  bool
  checkVectorEqual(Array1D<double> const& vector0,
                   Array1D<double> const& vector1,
                   double tolerance);
  
  
  template <class Intrinsics>
  void
  compute2DTestData(std::vector< Vector2D<double> >& points2D,
                    std::vector< Vector3D<double> >& points3D_camera,
                    Intrinsics const& intrinsics);

  void
  get3DTestData(std::vector< Vector3D<double> >& points3D_world,
                std::vector< Vector3D<double> >& points3D_camera,
                Transform3D<double>& cameraTworld);

  void
  getTestIntrinsicsPlumbBob(CameraIntrinsicsPlumbBob<double>& intrinsics);
  
  
  double m_defaultTolerance;
  double m_relaxedTolerance;
  double m_stringentTolerance;
  
}; // class CalibrationToolsRobustTest


/* ============== Member Function Definititions ============== */

CalibrationToolsRobustTest::
CalibrationToolsRobustTest()
  : brick::test::TestFixture<CalibrationToolsRobustTest>("CalibrationToolsRobustTest"),
    m_defaultTolerance(1.0E-10),
    m_relaxedTolerance(5.0E-5),
    m_stringentTolerance(1.0E-13)
{
  BRICK_TEST_REGISTER_MEMBER(testEstimateCameraParametersRobust);
  BRICK_TEST_REGISTER_MEMBER(testEstimateCameraParametersPinholeRobust);
}


void
CalibrationToolsRobustTest::
testEstimateCameraParametersRobust()
{
  CameraIntrinsicsPlumbBob<double> referenceIntrinsics;
  this->getTestIntrinsicsPlumbBob(referenceIntrinsics);
    
  std::vector< Vector3D<double> > points3D_world;
  std::vector< Vector3D<double> > points3D_camera;
  std::vector< Vector2D<double> > points2D;
  Transform3D<double> cameraTworld;
  this->get3DTestData(points3D_world, points3D_camera, cameraTworld);
  this->compute2DTestData(points2D, points3D_camera, referenceIntrinsics);

  // Corrupt the test data with an outlier.
  points3D_world[3] = points3D_world[5];

  // Important quantities for our robust search.
  //   Solution ought to be pretty good, as we haven't added any noise.
  double maxResidual = 0.01;
  //   Inlier probability is straightforward.
  double inlierProbability = (static_cast<double>(points2D.size() - 1)
                              / static_cast<double>(points2D.size()));

  //   This is a unit test... we'd rather it not fail.  Require
  //   confidence of 10,000,000 to one.
  double requiredConfidence = 0.9999999;
  
  // Now run the estimation code.
  CameraIntrinsicsPlumbBob<double> recoveredIntrinsics;
  Transform3D<double> recoveredCameraTworld;
  estimateCameraParametersRobust(
    recoveredIntrinsics, recoveredCameraTworld,
    referenceIntrinsics.getNumPixelsX(), referenceIntrinsics.getNumPixelsY(),
    points3D_world.begin(), points3D_world.end(), points2D.begin(),
    maxResidual, inlierProbability, requiredConfidence);

  BRICK_TEST_ASSERT(
    this->checkIntrinsicsEqual(
      recoveredIntrinsics, referenceIntrinsics, m_relaxedTolerance));
  BRICK_TEST_ASSERT(this->checkTransformEqual(
                    recoveredCameraTworld, cameraTworld, m_relaxedTolerance));
}


void
CalibrationToolsRobustTest::
testEstimateCameraParametersPinholeRobust()
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

  // Get some data with which to test.
  std::vector< Vector3D<double> > points3D_world;
  std::vector< Vector3D<double> > points3D_camera;
  std::vector< Vector2D<double> > points2D;
  Transform3D<double> cameraTworld;
  this->get3DTestData(points3D_world, points3D_camera, cameraTworld);
  this->compute2DTestData(points2D, points3D_camera, intrinsics);

  // Corrupt the test data with an outlier.
  points3D_camera[3] = points3D_camera[5];
  points3D_world[3] = points3D_world[5];

  // Important quantities for our robust search.
  //   Solution ought to be pretty good, as we haven't added any noise.
  double maxResidual = 0.01;
  //   Inlier probability is straightforward.
  double inlierProbability = (static_cast<double>(points2D.size() - 1)
                              / static_cast<double>(points2D.size()));

  //   This is a unit test... we'd rather it not fail.  Require
  //   confidence of 10,000,000 to one.
  double requiredConfidence = 0.9999999;

  // The "easy" case, where cameraTworld is identity.
  CameraIntrinsicsPinhole<double> recoveredIntrinsics0;
  brick::numeric::Transform3D<double> recoveredCameraTworld0;
  estimateCameraParametersPinholeRobust(
    recoveredIntrinsics0, recoveredCameraTworld0,
    intrinsics.getNumPixelsX(), intrinsics.getNumPixelsY(),
    points3D_camera.begin(), points3D_camera.end(), points2D.begin(),
    maxResidual, inlierProbability, requiredConfidence);
  BRICK_TEST_ASSERT(this->checkIntrinsicsEqual(
                    recoveredIntrinsics0, intrinsics, m_relaxedTolerance));
  Transform3D<double> identityXf;
  BRICK_TEST_ASSERT(this->checkTransformEqual(
                    recoveredCameraTworld0, identityXf, m_relaxedTolerance));

  // The "hard" case, where cameraTworld is not identity.
  CameraIntrinsicsPinhole<double> recoveredIntrinsics1;
  brick::numeric::Transform3D<double> recoveredCameraTworld1;
  estimateCameraParametersPinholeRobust(
    recoveredIntrinsics1, recoveredCameraTworld1,
    intrinsics.getNumPixelsX(), intrinsics.getNumPixelsY(),
    points3D_world.begin(), points3D_world.end(), points2D.begin(),
    maxResidual, inlierProbability, requiredConfidence);
  BRICK_TEST_ASSERT(this->checkIntrinsicsEqual(
                    recoveredIntrinsics1, intrinsics, m_relaxedTolerance));
  BRICK_TEST_ASSERT(this->checkTransformEqual(
                    recoveredCameraTworld1, cameraTworld, m_relaxedTolerance));
}


bool
CalibrationToolsRobustTest::
checkIntrinsicsEqual(CameraIntrinsicsPinhole<double> const& intrinsics0,
                     CameraIntrinsicsPinhole<double> const& intrinsics1,
                     double tolerance)
{
  if(intrinsics0.getNumPixelsX() != intrinsics1.getNumPixelsX()) {
    return false;
  }
  if(intrinsics0.getNumPixelsY() != intrinsics1.getNumPixelsY()) {
    return false;
  }
  if(!approximatelyEqual(
       intrinsics0.getKx(), intrinsics1.getKx(), tolerance)) {
    return false;
  }
  if(!approximatelyEqual(
       intrinsics0.getKy(), intrinsics1.getKy(), tolerance)) {
    return false;
  }
  if(!approximatelyEqual(
       intrinsics0.getCenterU(), intrinsics1.getCenterU(), tolerance)) {
    return false;
  }
  if(!approximatelyEqual(
       intrinsics0.getCenterV(), intrinsics1.getCenterV(), tolerance)) {
    return false;
  }
  return true;
}


bool
CalibrationToolsRobustTest::
checkIntrinsicsEqual(CameraIntrinsicsPlumbBob<double> const& intrinsics0,
                     CameraIntrinsicsPlumbBob<double> const& intrinsics1,
                     double tolerance)
{
  if(intrinsics0.getNumPixelsX() != intrinsics1.getNumPixelsX()) {
    return false;
  }
  if(intrinsics0.getNumPixelsY() != intrinsics1.getNumPixelsY()) {
    return false;
  }
  if(!approximatelyEqual(
       intrinsics0.getFocalLengthX(), intrinsics1.getFocalLengthX(),
       tolerance)) {
    return false;
  }
  if(!approximatelyEqual(
       intrinsics0.getFocalLengthY(), intrinsics1.getFocalLengthY(),
       tolerance)) {
    return false;
  }
  if(!approximatelyEqual(
       intrinsics0.getCenterU(), intrinsics1.getCenterU(), tolerance)) {
    return false;
  }
  if(!approximatelyEqual(
       intrinsics0.getCenterV(), intrinsics1.getCenterV(), tolerance)) {
    return false;
  }
  if(!approximatelyEqual(
       intrinsics0.getSkewCoefficient(), intrinsics1.getSkewCoefficient(),
       tolerance)) {
    return false;
  }
  if(!approximatelyEqual(
       intrinsics0.getRadialCoefficient0(), intrinsics1.getRadialCoefficient0(),
       tolerance)) {
    return false;
  }
  if(!approximatelyEqual(
       intrinsics0.getRadialCoefficient1(), intrinsics1.getRadialCoefficient1(),
       tolerance)) {
    return false;
  }
  if(!approximatelyEqual(
       intrinsics0.getRadialCoefficient2(), intrinsics1.getRadialCoefficient2(),
       tolerance)) {
    return false;
  }
  if(!approximatelyEqual(
       intrinsics0.getTangentialCoefficient0(),
       intrinsics1.getTangentialCoefficient0(),
       tolerance)) {
    return false;
  }
  if(!approximatelyEqual(
       intrinsics0.getTangentialCoefficient1(),
       intrinsics1.getTangentialCoefficient1(),
       tolerance)) {
    return false;
  }
  return true;
}


bool
CalibrationToolsRobustTest::
checkTransformEqual(Transform3D<double> const& transform0,
                    Transform3D<double> const& transform1,
                    double tolerance)
{
  for(unsigned int ii = 0; ii < 4; ++ii) {
    for(unsigned int jj = 0; jj < 4; ++jj) {
      if(!approximatelyEqual(
           transform0(ii, jj), transform1(ii, jj), tolerance)) {
        return false;
      }
    }
  }
  return true;
}


bool
CalibrationToolsRobustTest::
checkVectorEqual(Array1D<double> const& vector0,
                 Array1D<double> const& vector1,
                 double tolerance)
{
  if(vector0.size() != vector1.size()) {
    return false;
  }
  for(unsigned int ii = 0; ii < vector0.size(); ++ii) {
    if(!approximatelyEqual(vector0[ii], vector1[ii], tolerance)) {
      return false;
    }
  }
  return true;
}

    
template<class Intrinsics>
void
CalibrationToolsRobustTest::
compute2DTestData(std::vector< Vector2D<double> >& points2D,
                  std::vector< Vector3D<double> >& points3D_camera,
                  Intrinsics const& intrinsics)
{                  
  // Generate corresponding 2D poinst.
  points2D.resize(points3D_camera.size());
  for(unsigned int ii = 0; ii < points3D_camera.size(); ++ii) {
    points2D[ii] = intrinsics.project(points3D_camera[ii]);
  }
}


void
CalibrationToolsRobustTest::
get3DTestData(std::vector< Vector3D<double> >& points3D_world,
              std::vector< Vector3D<double> >& points3D_camera,
              Transform3D<double>& cameraTworld)
{
#if 1 /* Setting this to 0 makes for easier debugging. */
  cameraTworld = rollPitchYawToTransform3D(Vector3D<double>(0.1, -0.2, 0.3));
  cameraTworld.setValue<0, 3>(0.2);
  cameraTworld.setValue<1, 3>(-0.7);
  cameraTworld.setValue<2, 3>(-0.2);
#else
  cameraTworld.setValue(1.0, 0.0, 0.0, 0.0,
                        0.0, 1.0, 0.0, 0.0,
                        0.0, 0.0, 1.0, 0.0,
                        0.0, 0.0, 0.0, 1.0);
#endif

  // Generate 3D points.
#if 1  /* Setting this to 0 makes for easier debugging of functions that */
       /* don't require many points.                                     */
  for(double zCoord = 1.0; zCoord < 10.0; zCoord += 2.4) {
    for(double yCoord = -1.0; yCoord < 1.0; yCoord += 0.2) {
      for(double xCoord = -1.0; xCoord < 1.0; xCoord += 0.31) {
#else
  for(double zCoord = 1.0; zCoord < 10.0; zCoord += 7.0) {
    for(double yCoord = -1.0; yCoord < 1.0; yCoord += 1.0) {
      for(double xCoord = -1.0; xCoord < 1.0; xCoord += 1.0) {
#endif
        points3D_world.push_back(Vector3D<double>(xCoord, yCoord, zCoord));
        points3D_camera.push_back(cameraTworld * (*(points3D_world.end() - 1)));
      }
    }
  }
}


void
CalibrationToolsRobustTest::
getTestIntrinsicsPlumbBob(CameraIntrinsicsPlumbBob<double>& intrinsics)
{
  // Arbitrary camera params.
  size_t numPixelsX = 320;
  size_t numPixelsY = 240;
  double focalLengthX = 30.0;
  double focalLengthY = 15.0;
  double centerU = 100.0;
  double centerV = 125.0;
  double skewCoefficient = 0.001;
  double radialCoefficient0 = 0.02;
  double radialCoefficient1 = 0.0001;
  double radialCoefficient2 = 0.000007;
  double tangentialCoefficient0 = -0.01;
  double tangentialCoefficient1 = 0.005;
  CameraIntrinsicsPlumbBob<double> referenceIntrinsics(
    numPixelsX, numPixelsY, focalLengthX, focalLengthY, centerU, centerV,
    skewCoefficient, radialCoefficient0, radialCoefficient1, radialCoefficient2,
    tangentialCoefficient0, tangentialCoefficient1);
  intrinsics = referenceIntrinsics;
}


#if 0

int main(int argc, char** argv)
{
  CalibrationToolsRobustTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  CalibrationToolsRobustTest currentTest;

}

#endif
