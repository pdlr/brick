/**
***************************************************************************
* @file cameraIntrinsicsPlumbBobTest.cc
*
* Source file defining tests for the CameraIntrinsicsPlumbBob class.
*
* Copyright (C) 2006-2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/numeric/differentiableScalar.hh>
#include <brick/common/functional.hh>
#include <brick/computerVision/cameraIntrinsicsPlumbBob.hh>
#include <brick/optimization/gradientFunction.hh>
#include <brick/test/testFixture.hh>

using namespace brick::common;
using namespace brick::computerVision;
using namespace brick::geometry;
using namespace brick::numeric;
using namespace brick::optimization;
using namespace brick::test;


class CameraIntrinsicsPlumbBobTest
  : public TestFixture<CameraIntrinsicsPlumbBobTest> {

public:

  CameraIntrinsicsPlumbBobTest();
  ~CameraIntrinsicsPlumbBobTest() {}

  void setUp(const std::string& /* testName */) {}
  void tearDown(const std::string& /* testName */) {}

  // Tests.
  void testReverseProjectionObjectiveApplicationOperator();
  void testReverseProjectionObjectiveGradient();
  void testConstructor__void();
  void testConstructor__args();
  void testProject();
  void testReverseProject();
  void testReverseProjectEM();
  void testStreamOperators();
  void testReverseProjectWithJacobian();

private:

  CameraIntrinsicsPlumbBob<double>
  getIntrinsicsInstance();

  CameraIntrinsicsPlumbBob<double>
  getIntrinsicsInstanceMild();

  double m_defaultTolerance;
  double m_gradientTolerance;
  double m_reconstructionTolerance;
  double m_reverseProjectionTolerance;
  double m_reverseProjectionGradTolerance;

  const unsigned int m_numPixelsX;
  const unsigned int m_numPixelsY;
  const double m_focalLengthX;
  const double m_focalLengthY;
  const double m_centerU;
  const double m_centerV;
  const double m_skewCoefficient;
  const double m_radialCoefficient0;
  const double m_radialCoefficient1;
  const double m_radialCoefficient2;
  const double m_tangentialCoefficient0;
  const double m_tangentialCoefficient1;
}; // class CameraIntrinsicsPlumbBobTest


/* ============== Member Function Definititions ============== */

CameraIntrinsicsPlumbBobTest::
CameraIntrinsicsPlumbBobTest()
  : brick::test::TestFixture<CameraIntrinsicsPlumbBobTest>("CameraIntrinsicsPlumbBobTest"),
    m_defaultTolerance(1.0E-10),
    m_gradientTolerance(1.0E-5),
    // xxx
    m_reconstructionTolerance(1.0E-12),
    m_reverseProjectionTolerance(1.0E-4),
    m_reverseProjectionGradTolerance(1.0E-2),
    m_numPixelsX(320),
    m_numPixelsY(240),
    m_focalLengthX(30.0),
    m_focalLengthY(15.0),
    m_centerU(100.0),
    m_centerV(125.0),
    m_skewCoefficient(0.001),
    m_radialCoefficient0(0.02),
    m_radialCoefficient1(0.0001),
    m_radialCoefficient2(0.000007),
    m_tangentialCoefficient0(-0.01),
    m_tangentialCoefficient1(0.005)
{
  BRICK_TEST_REGISTER_MEMBER(testReverseProjectionObjectiveApplicationOperator);
  BRICK_TEST_REGISTER_MEMBER(testReverseProjectionObjectiveGradient);
  BRICK_TEST_REGISTER_MEMBER(testConstructor__void);
  BRICK_TEST_REGISTER_MEMBER(testConstructor__args);
  BRICK_TEST_REGISTER_MEMBER(testProject);
  BRICK_TEST_REGISTER_MEMBER(testReverseProject);
  BRICK_TEST_REGISTER_MEMBER(testReverseProjectEM);
  BRICK_TEST_REGISTER_MEMBER(testStreamOperators);
  BRICK_TEST_REGISTER_MEMBER(testReverseProjectWithJacobian);
}


void
CameraIntrinsicsPlumbBobTest::
testReverseProjectionObjectiveApplicationOperator()
{
  // Arbitrary camera params.
  CameraIntrinsicsPlumbBob<double> intrinsics = this->getIntrinsicsInstance();

  for(double yCoord = -1.0; yCoord < 1.0; yCoord += 0.1) {
    for(double xCoord = -1.0; xCoord < 1.0; xCoord += 0.1) {
      Vector3D<double> cameraCoord(xCoord, yCoord, 1.0);
      Vector2D<double> pixelCoord = intrinsics.project(cameraCoord);
      brick::computerVision::privateCode::ReverseProjectionObjective<double>
        objective(intrinsics, pixelCoord, 0.0, 0.0);

      Array1D<double> theta(2);
      double result;
      double offset;

      theta[0] = xCoord; theta[1] = yCoord;
      result = objective(theta);
      offset = objective.getOffset();
      BRICK_TEST_ASSERT(approximatelyEqual(result, offset, m_defaultTolerance));

      theta[0] = xCoord + 0.1; theta[1] = yCoord;
      result = objective(theta);
      BRICK_TEST_ASSERT(result > offset + m_defaultTolerance);

      theta[0] = xCoord - 0.1; theta[1] = yCoord;
      result = objective(theta);
      BRICK_TEST_ASSERT(result > offset + m_defaultTolerance);

      theta[0] = xCoord; theta[1] = yCoord + 0.1;
      result = objective(theta);
      BRICK_TEST_ASSERT(result > offset + m_defaultTolerance);

      theta[0] = xCoord; theta[1] = yCoord - 0.1;
      result = objective(theta);
      BRICK_TEST_ASSERT(result > offset + m_defaultTolerance);
    }
  }
}


void
CameraIntrinsicsPlumbBobTest::
testReverseProjectionObjectiveGradient()
{
  // Arbitrary camera params.
  CameraIntrinsicsPlumbBob<double> intrinsics = this->getIntrinsicsInstance();

  for(double yCoord = -1.0; yCoord < 1.0; yCoord += 0.1) {
    for(double xCoord = -1.0; xCoord < 1.0; xCoord += 0.1) {
      Vector3D<double> cameraCoord(xCoord, yCoord, 1.0);
      Vector2D<double> pixelCoord = intrinsics.project(cameraCoord);
      brick::computerVision::privateCode::ReverseProjectionObjective<double>
        objective(intrinsics, pixelCoord, 0.0, 0.0);

      Array1D<double> theta(2);
      theta[0] = xCoord; theta[1] = yCoord;
      Array1D<double> gradient = objective.gradient(theta);

      BRICK_TEST_ASSERT(gradient.size() == 2);
      BRICK_TEST_ASSERT(approximatelyEqual(gradient[0], 0.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(approximatelyEqual(gradient[1], 0.0, m_defaultTolerance));
    }
  }

  Vector2D<double> pixelCoord = intrinsics.project(
    Vector3D<double>(1.0, -1.0, 1.0));

  typedef brick::computerVision::privateCode::ReverseProjectionObjective<double>
    Objective;

  Objective objective(intrinsics, pixelCoord, 0.0, 0.0);
  GradientFunction<Objective> refObjective(objective);

  for(double yCoord = -1.0; yCoord < 1.0; yCoord += 0.1) {
    for(double xCoord = -1.0; xCoord < 1.0; xCoord += 0.1) {
      Array1D<double> theta(2);
      theta[0] = xCoord; theta[1] = yCoord;

      Array1D<double> testGradient = objective.gradient(theta);
      Array1D<double> refGradient = refObjective.gradient(theta);

      BRICK_TEST_ASSERT(testGradient.size() == 2);
      BRICK_TEST_ASSERT(refGradient.size() == 2);
      BRICK_TEST_ASSERT(
        approximatelyEqual(testGradient[0], refGradient[0],
                           m_gradientTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(testGradient[1], refGradient[1],
                           m_gradientTolerance));
    }
  }

  Objective objective_fov(intrinsics, pixelCoord, 0.62, 0.43);
  GradientFunction<Objective> refObjective_fov(objective_fov);

  for(double yCoord = -1.0; yCoord < 1.0; yCoord += 0.1) {
    for(double xCoord = -1.0; xCoord < 1.0; xCoord += 0.1) {
      Array1D<double> theta(2);
      theta[0] = xCoord; theta[1] = yCoord;

      Array1D<double> testGradient = objective_fov.gradient(theta);
      Array1D<double> refGradient = refObjective_fov.gradient(theta);

      // We scale m_gradientTolerance here because the
      // field-of-view-bounded errors get huge when you're far out of
      // bounds.  Scaling avoids failures due to numerical issues.
      double epsScale = (std::max(brick::numeric::absoluteValue(refGradient[0]),
                                  brick::numeric::absoluteValue(refGradient[1]))
                         * 10.0);
      epsScale = std::max(epsScale, 1.0);

      BRICK_TEST_ASSERT(testGradient.size() == 2);
      BRICK_TEST_ASSERT(refGradient.size() == 2);
      BRICK_TEST_ASSERT(
        approximatelyEqual(testGradient[0], refGradient[0],
                           epsScale * m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(testGradient[1], refGradient[1],
                           epsScale * m_defaultTolerance));
    }
  }
}


void
CameraIntrinsicsPlumbBobTest::
testConstructor__void()
{
  // Pass.
}


void
CameraIntrinsicsPlumbBobTest::
testConstructor__args()
{
  // Tested in testReverseProject().
}


void
CameraIntrinsicsPlumbBobTest::
testProject()
{
  // Arbitrary camera params.
  CameraIntrinsicsPlumbBob<double> intrinsics = this->getIntrinsicsInstance();

  for(double zCoord = 1.0; zCoord < 10.0; zCoord += 0.7) {
    for(double yCoord = -1.0; yCoord < 1.0; yCoord += 0.1) {
      for(double xCoord = -1.0; xCoord < 1.0; xCoord += 0.1) {
        Vector3D<double> cameraCoord(xCoord, yCoord, zCoord);
        Vector2D<double> pixelCoord = intrinsics.project(cameraCoord);

        double inverseZ = 1.0 / zCoord;
        double xNorm = xCoord * inverseZ;
        double yNorm = yCoord * inverseZ;
        double r2 = (xNorm * xNorm) + (yNorm * yNorm);
        double r4 = r2 * r2;
        double r6 = r2 * r2 * r2;

        // Radial distortion.
        double radialDistortion =
          (1.0 + m_radialCoefficient0 * r2 + m_radialCoefficient1 * r4
           + m_radialCoefficient2 * r6);
        double xDistorted0 = xNorm * radialDistortion;
        double yDistorted0 = yNorm * radialDistortion;

        // Tangential distortion.
        double a1 = 2.0 * xNorm * yNorm;
        double a2 = r2 + 2 * xNorm * xNorm;
        double a3 = r2 + 2 * yNorm * yNorm;
        double xTangential =
          m_tangentialCoefficient0 * a1 + m_tangentialCoefficient1 * a2;
        double yTangential =
          m_tangentialCoefficient0 * a3 + m_tangentialCoefficient1 * a1;
        double xDistorted1 = xDistorted0 + xTangential;
        double yDistorted1 = yDistorted0 + yTangential;

        // Skew.
        double xDistorted2 = xDistorted1 + m_skewCoefficient * yDistorted1;
        double yDistorted2 = yDistorted1;

        // Pinhole model.
        double referenceU = xDistorted2 * m_focalLengthX + m_centerU;
        double referenceV = yDistorted2 * m_focalLengthY + m_centerV;

        BRICK_TEST_ASSERT(
          approximatelyEqual(pixelCoord.x(), referenceU, m_defaultTolerance));
        BRICK_TEST_ASSERT(
          approximatelyEqual(pixelCoord.y(), referenceV, m_defaultTolerance));
      }
    }
  }
}


void
CameraIntrinsicsPlumbBobTest::
testReverseProject()
{
  // We test by round trip against convertWorldPointToPixel(), which
  // has its own independent test.

  // Arbitrary camera params.
  CameraIntrinsicsPlumbBob<double> intrinsics = this->getIntrinsicsInstance();

  for(double vCoord = 0.0; vCoord < m_numPixelsY; vCoord += 10.2) {
    for(double uCoord = 0.0; uCoord < m_numPixelsX; uCoord += 10.2) {
      Vector2D<double> pixelCoord(uCoord, vCoord);
      Ray3D<double> ray = intrinsics.reverseProject(pixelCoord);
      Vector3D<double> pointOnRay =
        ray.getOrigin() + 12.0 * ray.getDirectionVector();
      Vector2D<double> recoveredPixelCoord = intrinsics.project(pointOnRay);

      double residual = magnitude<double>(recoveredPixelCoord - pixelCoord);
      BRICK_TEST_ASSERT(residual < m_reverseProjectionTolerance);
    }
  }
}


void
CameraIntrinsicsPlumbBobTest::
testReverseProjectEM()
{
  // We test by round trip against convertWorldPointToPixel(), which
  // has its own independent test.

  // Arbitrary camera params.
  CameraIntrinsicsPlumbBob<double> intrinsics =
    this->getIntrinsicsInstanceMild();

  for(double vCoord = 0.0; vCoord < m_numPixelsY; vCoord += 10.2) {
    for(double uCoord = 0.0; uCoord < m_numPixelsX; uCoord += 10.2) {
      Vector2D<double> pixelCoord(uCoord, vCoord);
      Ray3D<double> ray = intrinsics.reverseProjectEM(
        pixelCoord, true, 1.0E-7);
      Vector3D<double> pointOnRay =
        ray.getOrigin() + 12.0 * ray.getDirectionVector();
      Vector2D<double> recoveredPixelCoord = intrinsics.project(pointOnRay);

      double residual = magnitude<double>(recoveredPixelCoord - pixelCoord);
      BRICK_TEST_ASSERT(residual < m_reverseProjectionTolerance);
    }
  }
}


void
CameraIntrinsicsPlumbBobTest::
testStreamOperators()
{
  CameraIntrinsicsPlumbBob<double> refIntrinsics = this->getIntrinsicsInstance();
  CameraIntrinsicsPlumbBob<double> testIntrinsics;

  std::ostringstream outputStream;
  outputStream << refIntrinsics;
  std::istringstream inputStream(outputStream.str());
  inputStream >> testIntrinsics;

  BRICK_TEST_ASSERT(
    testIntrinsics.getNumPixelsX() == refIntrinsics.getNumPixelsX());
  BRICK_TEST_ASSERT(
    testIntrinsics.getNumPixelsY() == refIntrinsics.getNumPixelsY());
  BRICK_TEST_ASSERT(approximatelyEqual(testIntrinsics.getFocalLengthX(),
                                     refIntrinsics.getFocalLengthX(),
                                     m_reconstructionTolerance));
  BRICK_TEST_ASSERT(approximatelyEqual(testIntrinsics.getFocalLengthY(),
                                     refIntrinsics.getFocalLengthY(),
                                     m_reconstructionTolerance));
  BRICK_TEST_ASSERT(approximatelyEqual(testIntrinsics.getCenterU(),
                                     refIntrinsics.getCenterU(),
                                     m_reconstructionTolerance));
  BRICK_TEST_ASSERT(approximatelyEqual(testIntrinsics.getCenterV(),
                                     refIntrinsics.getCenterV(),
                                     m_reconstructionTolerance));
  BRICK_TEST_ASSERT(approximatelyEqual(testIntrinsics.getSkewCoefficient(),
                                     refIntrinsics.getSkewCoefficient(),
                                     m_reconstructionTolerance));
  BRICK_TEST_ASSERT(approximatelyEqual(testIntrinsics.getRadialCoefficient0(),
                                     refIntrinsics.getRadialCoefficient0(),
                                     m_reconstructionTolerance));
  BRICK_TEST_ASSERT(approximatelyEqual(testIntrinsics.getRadialCoefficient1(),
                                     refIntrinsics.getRadialCoefficient1(),
                                     m_reconstructionTolerance));
  BRICK_TEST_ASSERT(approximatelyEqual(testIntrinsics.getRadialCoefficient2(),
                                     refIntrinsics.getRadialCoefficient2(),
                                     m_reconstructionTolerance));
  BRICK_TEST_ASSERT(approximatelyEqual(testIntrinsics.getTangentialCoefficient0(),
                                     refIntrinsics.getTangentialCoefficient0(),
                                     m_reconstructionTolerance));
  BRICK_TEST_ASSERT(approximatelyEqual(testIntrinsics.getTangentialCoefficient1(),
                                     refIntrinsics.getTangentialCoefficient1(),
                                     m_reconstructionTolerance));
}


void
CameraIntrinsicsPlumbBobTest::
testReverseProjectWithJacobian()
{
  double constexpr epsilon = 1.0E-8;
  double constexpr requiredPrecision = 10E-10;
  std::size_t constexpr maximumIterations = 150;

  // 4 pinhole parameters + 6 distortion parameters + 2 image coordinates
  std::size_t constexpr numParameters = 12;

  Vector2D<double> const uEpsilon(epsilon, 0.0);
  Vector2D<double> const vEpsilon(0.0, epsilon);

  // Arbitrary camera params.
  CameraIntrinsicsPlumbBob<double> intrinsics =
    this->getIntrinsicsInstanceMild();

  for(double vCoord = 0.0; vCoord < m_numPixelsY; vCoord += 10.2) {
    for(double uCoord = 0.0; uCoord < m_numPixelsX; uCoord += 10.2) {
      Vector2D<double> imagePoint(uCoord, vCoord);
      Vector2D<double> rectifiedPoint;
      Array2D<double> jacobian;
      bool success = reverseProjectWithJacobian(
        rectifiedPoint, jacobian,
        imagePoint, intrinsics, requiredPrecision, maximumIterations);

      BRICK_TEST_ASSERT(success);
      BRICK_TEST_ASSERT(jacobian.rows() == 2);
      BRICK_TEST_ASSERT(jacobian.columns() == numParameters);

      // OK, the reverse projection appears to have gone OK.  Get
      // some ground truth to compare with.  Bear in mind that
      // reverseProjectEM() is independently tested.
      Ray3D<double> ray = intrinsics.reverseProjectEM(
        imagePoint, true, requiredPrecision, maximumIterations);
      Vector2D<double> rectifiedPointGT(
        ray.getDirectionVector().x() / ray.getDirectionVector().z(),
        ray.getDirectionVector().y() / ray.getDirectionVector().z());

      BRICK_TEST_ASSERT(
        approximatelyEqual(rectifiedPoint.x(), rectifiedPointGT.x(),
                           m_reverseProjectionTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(rectifiedPoint.y(), rectifiedPointGT.y(),
                           m_reverseProjectionTolerance));

      // Fill out most of the ground truth jacobian.
      Array1D<double> parameterVector = intrinsics.getParameters();
      Array2D<double> jacobianGT(2, numParameters);

      for(std::size_t ii = 0; ii < parameterVector.size(); ++ii) {
        Array1D<double> parameterVectorPlus = parameterVector.copy();
        Array1D<double> parameterVectorMinus = parameterVector.copy();
        parameterVectorPlus[ii] += epsilon;
        parameterVectorMinus[ii] -= epsilon;

        intrinsics.setParameters(parameterVectorPlus);
        Ray3D<double> rayPlus = intrinsics.reverseProjectEM(
          imagePoint, true, requiredPrecision, maximumIterations);

        intrinsics.setParameters(parameterVectorMinus);
        Ray3D<double> rayMinus = intrinsics.reverseProjectEM(
          imagePoint, true, requiredPrecision, maximumIterations);

        Vector2D<double> rectifiedPointPlus(
          rayPlus.getDirectionVector().x() / rayPlus.getDirectionVector().z(),
          rayPlus.getDirectionVector().y() / rayPlus.getDirectionVector().z());
        Vector2D<double> rectifiedPointMinus(
          rayMinus.getDirectionVector().x() / rayMinus.getDirectionVector().z(),
          rayMinus.getDirectionVector().y() / rayMinus.getDirectionVector().z()
          );
        jacobianGT(0, ii) =
          (rectifiedPointPlus.x() - rectifiedPointMinus.x()) / (2.0 * epsilon);
        jacobianGT(1, ii) =
          (rectifiedPointPlus.y() - rectifiedPointMinus.y()) / (2.0 * epsilon);
      }

      // Final two columns of the ground truth jacobian.
      Ray3D<double> rayPlus = intrinsics.reverseProjectEM(
        imagePoint + uEpsilon, true, 150);
      Ray3D<double> rayMinus = intrinsics.reverseProjectEM(
        imagePoint - uEpsilon, true, 150);
      Vector2D<double> rectifiedPointPlus(
        rayPlus.getDirectionVector().x() / rayPlus.getDirectionVector().z(),
        rayPlus.getDirectionVector().y() / rayPlus.getDirectionVector().z());
      Vector2D<double> rectifiedPointMinus(
        rayMinus.getDirectionVector().x() / rayMinus.getDirectionVector().z(),
        rayMinus.getDirectionVector().y() / rayMinus.getDirectionVector().z());
      jacobianGT(0, numParameters - 2) =
        (rectifiedPointPlus.x() - rectifiedPointMinus.x()) / (2.0 * epsilon);
      jacobianGT(1, numParameters - 2) =
        (rectifiedPointPlus.y() - rectifiedPointMinus.y()) / (2.0 * epsilon);

      // Final column of the ground truth jacobian.
      rayPlus = intrinsics.reverseProjectEM(
        imagePoint + vEpsilon, true, 150);
      rayMinus = intrinsics.reverseProjectEM(
        imagePoint - vEpsilon, true, 150);
      rectifiedPointPlus.setValue(
        rayPlus.getDirectionVector().x() / rayPlus.getDirectionVector().z(),
        rayPlus.getDirectionVector().y() / rayPlus.getDirectionVector().z());
      rectifiedPointMinus.setValue(
        rayMinus.getDirectionVector().x() / rayMinus.getDirectionVector().z(),
        rayMinus.getDirectionVector().y() / rayMinus.getDirectionVector().z());
      jacobianGT(0, numParameters - 1) =
        (rectifiedPointPlus.x() - rectifiedPointMinus.x()) / (2.0 * epsilon);
      jacobianGT(1, numParameters - 1) =
        (rectifiedPointPlus.y() - rectifiedPointMinus.y()) / (2.0 * epsilon);

      BRICK_TEST_ASSERT(jacobian.rows() == jacobianGT.rows());
      BRICK_TEST_ASSERT(jacobian.columns() == jacobianGT.columns());
      for(std::size_t row = 0; row < jacobian.rows(); ++row) {
        for(std::size_t column = 0; column < jacobian.columns(); ++column) {
          double threshold = std::max(m_reverseProjectionGradTolerance,
                                      (std::fabs(jacobianGT(row, column))
                                       * m_reverseProjectionGradTolerance));;
          BRICK_TEST_ASSERT(
            approximatelyEqual(jacobian(row, column), jacobianGT(row, column),
                               threshold));
        }
      }
    }
  }
}


CameraIntrinsicsPlumbBob<double>
CameraIntrinsicsPlumbBobTest::
getIntrinsicsInstance()
{
  // Arbitrary camera params.
  CameraIntrinsicsPlumbBob<double> intrinsics(
    m_numPixelsX, m_numPixelsY, m_focalLengthX, m_focalLengthY,
    m_centerU, m_centerV, m_skewCoefficient, m_radialCoefficient0,
    m_radialCoefficient1, m_radialCoefficient2,
    m_tangentialCoefficient0, m_tangentialCoefficient1);
  return intrinsics;
}


CameraIntrinsicsPlumbBob<double>
CameraIntrinsicsPlumbBobTest::
getIntrinsicsInstanceMild()
{
  // Arbitrary camera params.
  CameraIntrinsicsPlumbBob<double> intrinsics(
    m_numPixelsX, m_numPixelsY, m_focalLengthX, m_focalLengthY,
    m_centerU, m_centerV, m_skewCoefficient, m_radialCoefficient0 / 100.0,
    m_radialCoefficient1 / 100.0, m_radialCoefficient2 / 100.0,
    m_tangentialCoefficient0 / 100.0, m_tangentialCoefficient1 / 100.0);
  return intrinsics;
}


#if 1

int main(int, char**)
{
  CameraIntrinsicsPlumbBobTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  CameraIntrinsicsPlumbBobTest currentTest;

}

#endif
