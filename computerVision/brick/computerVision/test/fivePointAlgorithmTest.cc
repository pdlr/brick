/**
***************************************************************************
* @file brick/computerVision/fivePointAlgorithmTest.cc
*
* Source file defining tests for fivePointAlgorithm().
*
* Copyright (C) 2009,2012 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

// Warning(xxx): Sometimes these test fail.  Figure out why!!

#define BRICK_FPA_VERBOSE 0

// xxx
#include <iomanip>

#include <brick/computerVision/fivePointAlgorithm.hh>
#include <brick/computerVision/imageIO.hh>
#include <brick/linearAlgebra/linearAlgebra.hh>
#include <brick/numeric/vector2D.hh>
#include <brick/random/pseudoRandom.hh>
#include <brick/test/testFixture.hh>

namespace cmn = brick::common;
namespace linalg = brick::linearAlgebra;
namespace num = brick::numeric;

// Anonymous namespace for locally defined functions.
namespace {

  cmn::Float64 trace(num::Array2D<cmn::Float64> const& matrix)
  {
    size_t numberOfDiagElements = std::min(matrix.rows(), matrix.columns());
    cmn::Float64 result = 0.0;
    for(size_t ii = 0; ii < numberOfDiagElements; ++ii) {
      result += matrix(ii, ii);
    }
    return result;
  }

} // namespace

namespace brick {

  namespace computerVision {

    class FivePointAlgorithmTest
      : public brick::test::TestFixture<FivePointAlgorithmTest> {

    public:

      FivePointAlgorithmTest();
      ~FivePointAlgorithmTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testFivePointAlgorithm();
      void testFivePointAlgorithmRobust__Iter_Iter_Iter_size_t();
      void testFivePointAlgorithmRobust__Iter_Iter_Iter_Iter_size_t();
      void testGetCameraMotionFromEssentialMatrix();
      void testTriangulateCalibratedImagePoint();

    private:

      void
      getCameraPoses(std::vector< num::Transform3D<cmn::Float64> >& worldTcam0Vector,
                     std::vector< num::Transform3D<cmn::Float64> >& worldTcam1Vector);

      void
      getCameraPoses(std::vector< num::Transform3D<cmn::Float64> >& worldTcam0Vector,
                     std::vector< num::Transform3D<cmn::Float64> >& worldTcam1Vector,
                     std::vector< num::Transform3D<cmn::Float64> >& worldTcam2Vector);

      void
      getTestPoints3D(std::vector< num::Vector3D<cmn::Float64> >& pVector,
                      bool isExtraPoints = false);

      void
      getTestPoints(std::vector< num::Vector2D<cmn::Float64> >& qVector,
                    std::vector< num::Vector2D<cmn::Float64> >& qPrimeVector);

      void
      getTestPoints(std::vector< num::Vector2D<cmn::Float64> >& qVector,
                    std::vector< num::Vector2D<cmn::Float64> >& qPrimeVector,
                    size_t transformNumber);

      void
      getTestPoints(std::vector< num::Vector2D<cmn::Float64> >& qVector,
                    std::vector< num::Vector2D<cmn::Float64> >& qPrimeVector,
                    std::vector< num::Vector2D<cmn::Float64> >& qPrimePrimeVector,
                    size_t transformNumber);

      bool
      isApproximatelyEqual(const Array1D<cmn::Float64>& array0,
                           const Array1D<cmn::Float64>& array1);

      bool
      isApproximatelyEqual(const Array2D<cmn::Float64>& array0,
                           const Array2D<cmn::Float64>& array1);

      bool
      isApproximatelyEqual(const Vector3D<cmn::Float64>& vector0,
                           const Vector3D<cmn::Float64>& vector1);


      cmn::Float64 m_defaultTolerance;

    }; // class FivePointAlgorithmTest


    /* ============== Member Function Definititions ============== */

    FivePointAlgorithmTest::
    FivePointAlgorithmTest()
      : brick::test::TestFixture<FivePointAlgorithmTest>("FivePointAlgorithmTest"),
        m_defaultTolerance(1.0E-5)
    {
      BRICK_TEST_REGISTER_MEMBER(testFivePointAlgorithm);
      BRICK_TEST_REGISTER_MEMBER(
        testFivePointAlgorithmRobust__Iter_Iter_Iter_size_t);
      BRICK_TEST_REGISTER_MEMBER(
        testFivePointAlgorithmRobust__Iter_Iter_Iter_Iter_size_t);
      BRICK_TEST_REGISTER_MEMBER(testGetCameraMotionFromEssentialMatrix);
      BRICK_TEST_REGISTER_MEMBER(testTriangulateCalibratedImagePoint);
    }


    void
    FivePointAlgorithmTest::
    testFivePointAlgorithm()
    {
      num::Array2D<cmn::Float64> EOrig;
      std::vector< num::Vector2D<cmn::Float64> > qVector;
      std::vector< num::Vector2D<cmn::Float64> > qPrimeVector;
      this->getTestPoints(qVector, qPrimeVector);

      std::vector< num::Array2D<cmn::Float64> > EVector = fivePointAlgorithm<cmn::Float64>(
        qVector.begin(), qVector.end(), qPrimeVector.begin());

      // Make sure EE has the appropriate qualities:
      for(size_t ii = 0; ii < EVector.size(); ++ii) {
        num::Array2D<cmn::Float64> EE = EVector[ii];

        // Check that EE is not full rank.
        cmn::Float64 det = linalg::determinant(EE);
        // xxx std::cout << "Det: " << det << std::endl;
        BRICK_TEST_ASSERT(approximatelyEqual(det, 0.0, m_defaultTolerance));

        // Check that linear constraints are satisfied.
        for(size_t jj = 0; jj < qVector.size(); ++jj) {
          num::Array1D<cmn::Float64> qq(3);
          qq[0] = qVector[jj].x();
          qq[1] = qVector[jj].y();
          qq[2] = 1;
          num::Array1D<cmn::Float64> qPrime(3);
          qPrime[0] = qPrimeVector[jj].x();
          qPrime[1] = qPrimeVector[jj].y();
          qPrime[2] = 1.0;
          cmn::Float64 residual = num::dot<cmn::Float64>(
            qPrime, num::matrixMultiply<cmn::Float64>(EE, qq));
          BRICK_TEST_ASSERT(
            approximatelyEqual(residual, 0.0, m_defaultTolerance));
        }

        // Check that trace constraint is satisfied.
        num::Array2D<cmn::Float64> term0 =
          2.0 * num::matrixMultiply<cmn::Float64>(num::matrixMultiply<cmn::Float64>(EE, EE.transpose()), EE);
        num::Array2D<cmn::Float64> term1 =
          trace(num::matrixMultiply<cmn::Float64>(EE, EE.transpose())) * EE;
        num::Array2D<cmn::Float64> residualArray = term0 - term1;
        for(size_t kk = 0; kk < residualArray.size(); ++kk) {
          BRICK_TEST_ASSERT(
            approximatelyEqual(residualArray[kk], 0.0, m_defaultTolerance));
        }
      }
    }


    void
    FivePointAlgorithmTest::
    testFivePointAlgorithmRobust__Iter_Iter_Iter_size_t()
    {
      brick::random::PseudoRandom pRandom;
      if(!BRICK_FPA_VERBOSE) {
        pRandom.setCurrentSeed(0);
      } else {
        std::cout << "Seed0: " << pRandom.getCurrentSeed() << std::endl;
      }

      std::vector< num::Vector2D<cmn::Float64> > qVector;
      std::vector< num::Vector2D<cmn::Float64> > qPrimeVector;
      this->getTestPoints(qVector, qPrimeVector, 0);

      cmn::Float64 score;
      num::Array2D<cmn::Float64> EE = fivePointAlgorithmRobust<cmn::Float64>(
        qVector.begin(), qVector.end(), qPrimeVector.begin(), 10, 0.8, score,
        pRandom);

      // Check that EE is not full rank.
      cmn::Float64 det = linalg::determinant(EE);
      // xxx std::cout << "Det: " << det << std::endl;
      BRICK_TEST_ASSERT(approximatelyEqual(det, 0.0, m_defaultTolerance));

      // Shouldn't have significant residuals for this simple test case.
      BRICK_TEST_ASSERT(score < 0.1);

      // Check that linear constraints are satisfied.
      for(size_t jj = 0; jj < qVector.size(); ++jj) {
        num::Array1D<cmn::Float64> qq(3);
        qq[0] = qVector[jj].x();
        qq[1] = qVector[jj].y();
        qq[2] = 1;
        num::Array1D<cmn::Float64> qPrime(3);
        qPrime[0] = qPrimeVector[jj].x();
        qPrime[1] = qPrimeVector[jj].y();
        qPrime[2] = 1.0;
        cmn::Float64 residual = dot<cmn::Float64>(qPrime, num::matrixMultiply<cmn::Float64>(EE, qq));
        BRICK_TEST_ASSERT(
          approximatelyEqual(residual, 0.0, m_defaultTolerance));
      }

      // Check that trace constraint is satisfied.
      num::Array2D<cmn::Float64> term0 =
        2.0 * num::matrixMultiply<cmn::Float64>(num::matrixMultiply<cmn::Float64>(EE, EE.transpose()), EE);
      num::Array2D<cmn::Float64> term1 =
        trace(num::matrixMultiply<cmn::Float64>(EE, EE.transpose())) * EE;
      num::Array2D<cmn::Float64> residualArray = term0 - term1;
      for(size_t kk = 0; kk < residualArray.size(); ++kk) {
        BRICK_TEST_ASSERT(
          approximatelyEqual(residualArray[kk], 0.0, m_defaultTolerance));
      }
    }


    void
    FivePointAlgorithmTest::
    testFivePointAlgorithmRobust__Iter_Iter_Iter_Iter_size_t()
    {
      brick::random::PseudoRandom pRandom;
      if(!BRICK_FPA_VERBOSE) {
        pRandom.setCurrentSeed(0);
      } else {
        std::cout << "Seed1: " << pRandom.getCurrentSeed() << std::endl;
      }

      std::vector< num::Vector2D<cmn::Float64> > qVector;
      std::vector< num::Vector2D<cmn::Float64> > qPrimeVector;
      std::vector< num::Vector2D<cmn::Float64> > qPrimePrimeVector;
      this->getTestPoints(qVector, qPrimeVector, qPrimePrimeVector, 0);

      cmn::Float64 score;
      num::Array2D<cmn::Float64> EE;
      num::Transform3D<cmn::Float64> cam0Tcam2;
      num::Transform3D<cmn::Float64> cam1Tcam2;
      fivePointAlgorithmRobust<cmn::Float64>(
        qVector.begin(), qVector.end(), qPrimeVector.begin(),
        qPrimePrimeVector.begin(), 10, 0.8, EE, cam0Tcam2, cam1Tcam2, score,
        pRandom);

      // Check that EE is not full rank.
      cmn::Float64 det = linalg::determinant(EE);
      // xxx std::cout << "Det: " << det << std::endl;
      BRICK_TEST_ASSERT(approximatelyEqual(det, 0.0, m_defaultTolerance));

      // Shouldn't have significant residuals for this simple test case.
      BRICK_TEST_ASSERT(score < 0.1);

      // Check that linear constraints are satisfied.
      for(size_t jj = 0; jj < qVector.size(); ++jj) {
        num::Array1D<cmn::Float64> qq(3);
        qq[0] = qVector[jj].x();
        qq[1] = qVector[jj].y();
        qq[2] = 1;
        num::Array1D<cmn::Float64> qPrimePrime(3);
        qPrimePrime[0] = qPrimePrimeVector[jj].x();
        qPrimePrime[1] = qPrimePrimeVector[jj].y();
        qPrimePrime[2] = 1.0;
        cmn::Float64 residual = dot<cmn::Float64>(qPrimePrime, num::matrixMultiply<cmn::Float64>(EE, qq));
        BRICK_TEST_ASSERT(
          approximatelyEqual(residual, 0.0, m_defaultTolerance));
      }

      // Check that trace constraint is satisfied.
      num::Array2D<cmn::Float64> term0 =
        2.0 * num::matrixMultiply<cmn::Float64>(num::matrixMultiply<cmn::Float64>(EE, EE.transpose()), EE);
      num::Array2D<cmn::Float64> term1 =
        trace(num::matrixMultiply<cmn::Float64>(EE, EE.transpose())) * EE;
      num::Array2D<cmn::Float64> residualArray = term0 - term1;
      for(size_t kk = 0; kk < residualArray.size(); ++kk) {
        BRICK_TEST_ASSERT(
          approximatelyEqual(residualArray[kk], 0.0, m_defaultTolerance));
      }
    }


    void
    FivePointAlgorithmTest::
    testGetCameraMotionFromEssentialMatrix()
    {
      brick::random::PseudoRandom pRandom;
      if(!BRICK_FPA_VERBOSE) {
        pRandom.setCurrentSeed(0);
      } else {
        // xxx std::cout << "Seed2: " << pRandom.getCurrentSeed() << std::endl;
      }

      std::vector< num::Transform3D<cmn::Float64> > worldTcam0Vector;
      std::vector< num::Transform3D<cmn::Float64> > worldTcam1Vector;
      this->getCameraPoses(worldTcam0Vector, worldTcam1Vector);

      for(size_t ii = 0; ii < worldTcam0Vector.size(); ++ii) {
        num::Transform3D<cmn::Float64> worldTcam0 = worldTcam0Vector[ii];
        num::Transform3D<cmn::Float64> worldTcam1 = worldTcam1Vector[ii];
        num::Transform3D<cmn::Float64> cam0Tworld = worldTcam0.invert();
        num::Transform3D<cmn::Float64> cam1Tworld = worldTcam1.invert();
        num::Transform3D<cmn::Float64> cam1Tcam0 = cam1Tworld * worldTcam0;
        num::Transform3D<cmn::Float64> cam0Tcam1 = cam0Tworld * worldTcam1;

        // Start with some randomly generated input points.
        std::vector< num::Vector2D<cmn::Float64> > qVector;
        std::vector< num::Vector2D<cmn::Float64> > qPrimeVector;
        this->getTestPoints(qVector, qPrimeVector, ii);

        // Recover essential matrix.
        cmn::Float64 score;
        num::Array2D<cmn::Float64> EE = fivePointAlgorithmRobust<cmn::Float64>(
          qVector.begin(), qVector.end(), qPrimeVector.begin(), 10, 0.8, score,
          pRandom);

        // Shouldn't have significant residuals for this simple test case.
        BRICK_TEST_ASSERT(score < 0.1);

        num::Transform3D<cmn::Float64> recoveredCam1Tcam0 =
          getCameraMotionFromEssentialMatrix(EE, qVector[0], qPrimeVector[0]);
        num::Vector3D<cmn::Float64> recoveredTranslation(
          recoveredCam1Tcam0(0, 3), recoveredCam1Tcam0(1, 3),
          recoveredCam1Tcam0(2, 3));

        num::Transform3D<cmn::Float64> recoveredCam1Rcam0 = recoveredCam1Tcam0;
        recoveredCam1Rcam0.setValue(0, 3, 0.0);
        recoveredCam1Rcam0.setValue(1, 3, 0.0);
        recoveredCam1Rcam0.setValue(2, 3, 0.0);

        num::Vector3D<cmn::Float64> translation(
          cam1Tcam0(0, 3), cam1Tcam0(1, 3), cam1Tcam0(2, 3));

        num::Transform3D<cmn::Float64> cam0Rcam1 = cam0Tcam1;
        cam0Rcam1.setValue(0, 3, 0.0);
        cam0Rcam1.setValue(1, 3, 0.0);
        cam0Rcam1.setValue(2, 3, 0.0);

        cmn::Float64 transDot = num::dot<cmn::Float64>(recoveredTranslation, translation);
        num::Vector3D<cmn::Float64> transCross = num::cross(
          recoveredTranslation, translation);
        num::Transform3D<cmn::Float64> recoveredIdentity = recoveredCam1Rcam0 * cam0Rcam1;

        BRICK_TEST_ASSERT(transDot > 0.0);
        BRICK_TEST_ASSERT(
          num::magnitude<cmn::Float64>(transCross) < m_defaultTolerance);
        BRICK_TEST_ASSERT(
          std::fabs(recoveredIdentity(0, 0) - 1.0) < m_defaultTolerance);
        BRICK_TEST_ASSERT(
          std::fabs(recoveredIdentity(1, 1) - 1.0) < m_defaultTolerance);
        BRICK_TEST_ASSERT(
          std::fabs(recoveredIdentity(2, 2) - 1.0) < m_defaultTolerance);
        BRICK_TEST_ASSERT(
          std::fabs(recoveredIdentity(0, 1) < m_defaultTolerance));
        BRICK_TEST_ASSERT(
          std::fabs(recoveredIdentity(0, 2) < m_defaultTolerance));
        BRICK_TEST_ASSERT(
          std::fabs(recoveredIdentity(1, 0) < m_defaultTolerance));
        BRICK_TEST_ASSERT(
          std::fabs(recoveredIdentity(1, 2) < m_defaultTolerance));
        BRICK_TEST_ASSERT(
          std::fabs(recoveredIdentity(2, 0) < m_defaultTolerance));
        BRICK_TEST_ASSERT(
          std::fabs(recoveredIdentity(2, 1) < m_defaultTolerance));
      }
    }


    void
    FivePointAlgorithmTest::
    testTriangulateCalibratedImagePoint()
    {
      // Start with some randomly generated input points.
      std::vector< num::Vector3D<cmn::Float64> > targetVector;
      this->getTestPoints3D(targetVector);

      std::vector< num::Transform3D<cmn::Float64> > worldTcam0Vector;
      std::vector< num::Transform3D<cmn::Float64> > worldTcam1Vector;
      this->getCameraPoses(worldTcam0Vector, worldTcam1Vector);

      for(size_t ii = 0; ii < targetVector.size(); ++ii) {
        // Set up cameras and target points in world coordinates.
        num::Transform3D<cmn::Float64> worldTcam0 = worldTcam0Vector[ii];
        num::Transform3D<cmn::Float64> worldTcam1 = worldTcam1Vector[ii];
        num::Transform3D<cmn::Float64> cam0Tworld = worldTcam0.invert();
        num::Transform3D<cmn::Float64> cam1Tworld = worldTcam1.invert();
        num::Transform3D<cmn::Float64> cam1Tcam0 = cam1Tworld * worldTcam0;
        num::Transform3D<cmn::Float64> cam0Tcam1 = cam0Tworld * worldTcam1;

        // Find target point in camera coordinates.
        num::Vector3D<cmn::Float64> target0 = cam0Tworld * targetVector[ii];
        num::Vector3D<cmn::Float64> target1 = cam1Tworld * targetVector[ii];

        // Find target point in calibrated image coordinates.
        num::Vector2D<cmn::Float64> q0(target0.x(), target0.y(), target0.z());
        num::Vector2D<cmn::Float64> q1(target1.x(), target1.y(), target1.z());

        // Triangulate.
        brick::numeric::Vector3D<cmn::Float64> recoveredTarget0 =
          triangulateCalibratedImagePoint(cam0Tcam1, q0, q1);
        brick::numeric::Vector3D<cmn::Float64> recoveredTarget1 =
          cam0Tcam1.invert() * recoveredTarget0;

        // Check result.
        BRICK_TEST_ASSERT(this->isApproximatelyEqual(recoveredTarget0, target0));
        BRICK_TEST_ASSERT(this->isApproximatelyEqual(recoveredTarget1, target1));
      }
    }


    void
    FivePointAlgorithmTest::
    getCameraPoses(std::vector< num::Transform3D<cmn::Float64> >& worldTcam0Vector,
                   std::vector< num::Transform3D<cmn::Float64> >& worldTcam1Vector)
    {
      std::vector< num::Transform3D<cmn::Float64> > worldTcam2Vector;
      this->getCameraPoses(worldTcam0Vector,worldTcam1Vector, worldTcam2Vector);
    }


    void
    FivePointAlgorithmTest::
    getCameraPoses(std::vector< num::Transform3D<cmn::Float64> >& worldTcam0Vector,
                   std::vector< num::Transform3D<cmn::Float64> >& worldTcam1Vector,
                   std::vector< num::Transform3D<cmn::Float64> >& worldTcam2Vector)
    {
      worldTcam0Vector.clear();
      worldTcam1Vector.clear();
      worldTcam2Vector.clear();

      std::vector< num::Vector3D<cmn::Float64> > center0Vector;
      std::vector< num::Vector3D<cmn::Float64> > center1Vector;
      std::vector< num::Vector3D<cmn::Float64> > center2Vector;
      std::vector<cmn::Float64> angle0Vector;
      std::vector<cmn::Float64> angle1Vector;
      std::vector<cmn::Float64> angle2Vector;

      center0Vector.push_back(num::Vector3D<cmn::Float64>(-7, 4, -9));
      center1Vector.push_back(num::Vector3D<cmn::Float64>(-9, -7, 3));
      center2Vector.push_back(num::Vector3D<cmn::Float64>(-5, 5, 0));
      angle0Vector.push_back(0.05);
      angle1Vector.push_back(-1.0);
      angle2Vector.push_back(-0.25);

      center0Vector.push_back(num::Vector3D<cmn::Float64>(-7, 1, -1));
      center1Vector.push_back(num::Vector3D<cmn::Float64>(0, 6, -7));
      center2Vector.push_back(num::Vector3D<cmn::Float64>(2, 4, -2));
      angle0Vector.push_back(0.75);
      angle1Vector.push_back(0.3);
      angle2Vector.push_back(0.1);

      center0Vector.push_back(num::Vector3D<cmn::Float64>(-3, -2, 10));
      center1Vector.push_back(num::Vector3D<cmn::Float64>(-1, -7, 9));
      center2Vector.push_back(num::Vector3D<cmn::Float64>(0, 0, -1));
      angle0Vector.push_back(-0.24);
      angle1Vector.push_back(-0.09);
      angle2Vector.push_back(0.09);

      center0Vector.push_back(num::Vector3D<cmn::Float64>(-9, -7, -8));
      center1Vector.push_back(num::Vector3D<cmn::Float64>(10, -2, 9));
      center2Vector.push_back(num::Vector3D<cmn::Float64>(4, 4, -3));
      angle0Vector.push_back(0.15);
      angle1Vector.push_back(1.1);
      angle2Vector.push_back(-0.7);

      center0Vector.push_back(num::Vector3D<cmn::Float64>(0, 3, 3));
      center1Vector.push_back(num::Vector3D<cmn::Float64>(-7, -3, -3));
      center2Vector.push_back(num::Vector3D<cmn::Float64>(-5, 3, 3));
      angle0Vector.push_back(-0.2);
      angle1Vector.push_back(0.3);
      angle2Vector.push_back(0.3);

      for(size_t ii = 0; ii < center0Vector.size(); ++ii) {
        // Set up cameras and target point in world coordinates.
        cmn::Float64 cos0 = std::cos(angle0Vector[ii]);
        cmn::Float64 sin0 = std::sin(angle0Vector[ii]);
        cmn::Float64 cos1 = std::cos(angle1Vector[ii]);
        cmn::Float64 sin1 = std::sin(angle1Vector[ii]);
        cmn::Float64 cos2 = std::cos(angle2Vector[ii]);
        cmn::Float64 sin2 = std::sin(angle2Vector[ii]);
        worldTcam0Vector.push_back(
          num::Transform3D<cmn::Float64>(cos0, 0.0, -sin0, center0Vector[ii].x(),
                           0.0, 1.0,   0.0, center0Vector[ii].y(),
                           sin0, 0.0,  cos0, center0Vector[ii].z(),
                           0.0, 0.0,   0.0, 1.0));
        worldTcam1Vector.push_back(
          num::Transform3D<cmn::Float64>(cos1, -sin1, 0.0, center1Vector[ii].x(),
                           sin1,  cos1, 0.0, center1Vector[ii].y(),
                           0.0,   0.0, 1.0, center1Vector[ii].z(),
                           0.0,   0.0, 0.0, 1.0));
        worldTcam2Vector.push_back(
          num::Transform3D<cmn::Float64>(cos2, -sin2, 0.0, center2Vector[ii].x(),
                           sin2,  cos2, 0.0, center2Vector[ii].y(),
                           0.0,   0.0, 1.0, center2Vector[ii].z(),
                           0.0,   0.0, 0.0, 1.0));
      }
    }


    void
    FivePointAlgorithmTest::
    getTestPoints(std::vector< num::Vector2D<cmn::Float64> >& qVector,
                  std::vector< num::Vector2D<cmn::Float64> >& qPrimeVector)
    {
      // Create a series of arbitrary q, qPrime pairs (where q and
      // qPrimae are stored in different vectors.
      qVector.clear();
      qVector.push_back(num::Vector2D<cmn::Float64>(2.0, 2.0));
      qVector.push_back(num::Vector2D<cmn::Float64>(1.0, 3.0));
      qVector.push_back(num::Vector2D<cmn::Float64>(0.0, 2.0));
      qVector.push_back(num::Vector2D<cmn::Float64>(-1.0, 4.0));
      qVector.push_back(num::Vector2D<cmn::Float64>(-2.0, 5.0));

      for(size_t ii = 0; ii < qVector.size(); ++ii) {
        qVector[ii] += num::Vector2D<cmn::Float64>(10.0, 20.0);
      }

      qPrimeVector.clear();
      qPrimeVector.push_back(num::Vector2D<cmn::Float64>(2.0, -2.0));
      qPrimeVector.push_back(num::Vector2D<cmn::Float64>(1.0, -3.0));
      qPrimeVector.push_back(num::Vector2D<cmn::Float64>(4.0, -4.0));
      qPrimeVector.push_back(num::Vector2D<cmn::Float64>(-1.0, 5.0));
      qPrimeVector.push_back(num::Vector2D<cmn::Float64>(3.0, 6.0));
    }


    void
    FivePointAlgorithmTest::
    getTestPoints(std::vector< num::Vector2D<cmn::Float64> >& qVector,
                  std::vector< num::Vector2D<cmn::Float64> >& qPrimeVector,
                  size_t transformNumber)
    {
      std::vector< num::Vector2D<cmn::Float64> > qPrimePrimeVector;
      this->getTestPoints(qVector, qPrimeVector, qPrimePrimeVector,
                          transformNumber);
    }


    void
    FivePointAlgorithmTest::
    getTestPoints(std::vector< num::Vector2D<cmn::Float64> >& qVector,
                  std::vector< num::Vector2D<cmn::Float64> >& qPrimeVector,
                  std::vector< num::Vector2D<cmn::Float64> >& qPrimePrimeVector,
                  size_t transformNumber)
    {
      qVector.clear();
      qPrimeVector.clear();
      qPrimePrimeVector.clear();

      std::vector< num::Vector3D<cmn::Float64> > targetVector;
      this->getTestPoints3D(targetVector, true);

      std::vector< num::Transform3D<cmn::Float64> > worldTcam0Vector;
      std::vector< num::Transform3D<cmn::Float64> > worldTcam1Vector;
      std::vector< num::Transform3D<cmn::Float64> > worldTcam2Vector;
      this->getCameraPoses(worldTcam0Vector, worldTcam1Vector,
                           worldTcam2Vector);

      num::Transform3D<cmn::Float64> worldTcam0 = worldTcam0Vector[transformNumber];
      num::Transform3D<cmn::Float64> worldTcam1 = worldTcam1Vector[transformNumber];
      num::Transform3D<cmn::Float64> worldTcam2 = worldTcam2Vector[transformNumber];
      num::Transform3D<cmn::Float64> cam0Tworld = worldTcam0.invert();
      num::Transform3D<cmn::Float64> cam1Tworld = worldTcam1.invert();
      num::Transform3D<cmn::Float64> cam2Tworld = worldTcam2.invert();
      num::Transform3D<cmn::Float64> cam1Tcam0 = cam1Tworld * worldTcam0;
      num::Transform3D<cmn::Float64> cam0Tcam1 = cam0Tworld * worldTcam1;
      num::Transform3D<cmn::Float64> cam2Tcam0 = cam2Tworld * worldTcam0;
      num::Transform3D<cmn::Float64> cam0Tcam2 = cam0Tworld * worldTcam2;
      num::Transform3D<cmn::Float64> cam2Tcam1 = cam2Tworld * worldTcam1;
      num::Transform3D<cmn::Float64> cam1Tcam2 = cam1Tworld * worldTcam2;

      for(size_t ii = 0; ii < targetVector.size(); ++ii) {
        // Find target point in camera coordinates.
        num::Vector3D<cmn::Float64> target0 = cam0Tworld * targetVector[ii];
        num::Vector3D<cmn::Float64> target1 = cam1Tworld * targetVector[ii];
        num::Vector3D<cmn::Float64> target2 = cam2Tworld * targetVector[ii];

        // Find target point in calibrated image coordinates.
        num::Vector2D<cmn::Float64> q0(target0.x(), target0.y(), target0.z());
        num::Vector2D<cmn::Float64> q1(target1.x(), target1.y(), target1.z());
        num::Vector2D<cmn::Float64> q2(target2.x(), target2.y(), target2.z());
        qVector.push_back(q0);
        qPrimeVector.push_back(q1);
        qPrimePrimeVector.push_back(q2);
      }
    }


    void
    FivePointAlgorithmTest::
    getTestPoints3D(std::vector< num::Vector3D<cmn::Float64> >& pVector, bool isExtraPoints)
    {
      // Create a series of arbitrary points.
      pVector.clear();

      // pVector.push_back(num::Vector2D<cmn::Float64>(2.0, 2.0, 15.0));
      // pVector.push_back(num::Vector2D<cmn::Float64>(1.0, 3.0, 12.0));
      // pVector.push_back(num::Vector2D<cmn::Float64>(0.0, 2.0, 19.0));
      // pVector.push_back(num::Vector2D<cmn::Float64>(-1.0, 4.0, 25.0));
      // pVector.push_back(num::Vector2D<cmn::Float64>(-2.0, 5.0, 12.0));

      pVector.push_back(num::Vector3D<cmn::Float64>(5, 8, 26));
      pVector.push_back(num::Vector3D<cmn::Float64>(-7, -5, 27));
      pVector.push_back(num::Vector3D<cmn::Float64>(1, 2, 11));
      pVector.push_back(num::Vector3D<cmn::Float64>(10, 1, 20));
      pVector.push_back(num::Vector3D<cmn::Float64>(-7, 8, 18));

      if(isExtraPoints) {
        pVector.push_back(num::Vector3D<cmn::Float64>(5, 7, 26));
        pVector.push_back(num::Vector3D<cmn::Float64>(-2, 5, 13));
        pVector.push_back(num::Vector3D<cmn::Float64>(-8, -8, 17));
        pVector.push_back(num::Vector3D<cmn::Float64>(-9, 2, 11));
        pVector.push_back(num::Vector3D<cmn::Float64>(6, 5, 14));
        pVector.push_back(num::Vector3D<cmn::Float64>(3, -5, 21));
        pVector.push_back(num::Vector3D<cmn::Float64>(4, 8, 18));
        pVector.push_back(num::Vector3D<cmn::Float64>(7, -3, 29));
        pVector.push_back(num::Vector3D<cmn::Float64>(-9, -1, 22));
        pVector.push_back(num::Vector3D<cmn::Float64>(4, -9, 18));
      }
    }


    bool
    FivePointAlgorithmTest::
    isApproximatelyEqual(const Array1D<cmn::Float64>& array0,
                         const Array1D<cmn::Float64>& array1)
    {
      if(array0.size() != array1.size()) {
        return false;
      }
      return std::equal(array0.begin(), array0.end(), array1.begin(),
                        ApproximatelyEqualFunctor<cmn::Float64>(m_defaultTolerance));
    }


    bool
    FivePointAlgorithmTest::
    isApproximatelyEqual(const Array2D<cmn::Float64>& array0,
                         const Array2D<cmn::Float64>& array1)
    {
      if(array0.rows() != array1.rows()) {
        return false;
      }
      if(array0.columns() != array1.columns()) {
        return false;
      }
      return std::equal(array0.begin(), array0.end(), array1.begin(),
                        ApproximatelyEqualFunctor<cmn::Float64>(m_defaultTolerance));
    }


    bool
    FivePointAlgorithmTest::
    isApproximatelyEqual(const Vector3D<cmn::Float64>& vector0,
                         const Vector3D<cmn::Float64>& vector1)
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
  brick::computerVision::FivePointAlgorithmTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::computerVision::FivePointAlgorithmTest currentTest;

}

#endif
