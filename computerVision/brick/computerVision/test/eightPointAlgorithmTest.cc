/**
***************************************************************************
* @file brick/computerVision/eightPointAlgorithmTest.cc
*
* Source file defining tests for eightPointAlgorithm().
*
* Copyright (C) 2008,2012 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/computerVision/test/testImages.hh>
#include <brick/computerVision/eightPointAlgorithm.hh>
#include <brick/computerVision/imageIO.hh>
#include <brick/test/testFixture.hh>

namespace cmn = brick::common;
namespace num = brick::numeric;

namespace brick {

  namespace computerVision {

    class EightPointAlgorithmTest
      : public brick::test::TestFixture<EightPointAlgorithmTest> {

    public:

      EightPointAlgorithmTest();
      ~EightPointAlgorithmTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testEightPointAlgorithm();
      void testNormalizePointSequence();

    private:

      void
      getTestPoints(std::vector< num::Vector2D<cmn::Float64> >& uVector,
                    std::vector< num::Vector2D<cmn::Float64> >& uPrimeVector,
                    num::Array2D<cmn::Float64>& FOrig);

      bool
      isApproximatelyEqual(const num::Array1D<cmn::Float64>& array0,
                           const num::Array1D<cmn::Float64>& array1);

      bool
      isApproximatelyEqual(const num::Array2D<cmn::Float64>& array0,
                           const num::Array2D<cmn::Float64>& array1);


      cmn::Float64 m_defaultTolerance;

    }; // class EightPointAlgorithmTest


    /* ============== Member Function Definititions ============== */

    EightPointAlgorithmTest::
    EightPointAlgorithmTest()
      : brick::test::TestFixture<EightPointAlgorithmTest>("EightPointAlgorithmTest"),
        m_defaultTolerance(1.0E-10)
    {
      BRICK_TEST_REGISTER_MEMBER(testEightPointAlgorithm);
      BRICK_TEST_REGISTER_MEMBER(testNormalizePointSequence);
    }


    void
    EightPointAlgorithmTest::
    testEightPointAlgorithm()
    {
      num::Array2D<cmn::Float64> FOrig;
      std::vector< num::Vector2D<cmn::Float64> > uVector;
      std::vector< num::Vector2D<cmn::Float64> > uPrimeVector;
      this->getTestPoints(uVector, uPrimeVector, FOrig);

      // Try to recover FOrig.
      num::Array2D<cmn::Float64> FEstimate = eightPointAlgorithm<cmn::Float64>(
        uVector.begin(), uVector.end(), uPrimeVector.begin());

      // Make sure FEstimate has the appropriate qualities:
      for(size_t ii = 0; ii < uVector.size(); ++ii) {
        num::Array1D<cmn::Float64> uu(3);
        uu[0] = uVector[ii].x();
        uu[1] = uVector[ii].y();
        uu[2] = 1.0;
        num::Array1D<cmn::Float64> uPrime(3);
        uPrime[0] = uPrimeVector[ii].x();
        uPrime[1] = uPrimeVector[ii].y();
        uPrime[2] = 1.0;

        num::Array1D<cmn::Float64> lineCoeffs =
          num::matrixMultiply<cmn::Float64>(FEstimate, uu);
        cmn::Float64 residual = num::dot<cmn::Float64>(uPrime, lineCoeffs);
        BRICK_TEST_ASSERT(std::fabs(residual) < m_defaultTolerance);
      }
    }


    void
    EightPointAlgorithmTest::
    testNormalizePointSequence()
    {
      // Create a series of input points, ostensibly coming from one image.
      num::Array2D<cmn::Float64> inputPoints("[[ 12.0,  22.0, 1.0],"
                                       " [ 11.0,  22.0, 1.0],"
                                       " [ 10.0,  22.0, 1.0],"
                                       " [ -1.0,  22.0, 1.0],"
                                       " [ -2.0,  22.0, 1.0],"
                                       " [ 12.0,  21.0, 1.0],"
                                       " [ 11.0,  21.0, 1.0],"
                                       " [ 10.0,  21.0, 1.0],"
                                       " [-11.0,  21.0, 1.0],"
                                       " [-12.0,  21.0, 1.0],"
                                       " [ 12.0,  20.0, 1.0],"
                                       " [ 11.0,  20.0, 1.0],"
                                       " [ 10.0,  20.0, 1.0],"
                                       " [ -1.0,  20.0, 1.0],"
                                       " [ -2.0,  20.0, 1.0],"
                                       " [ 12.0, -1.0, 1.0],"
                                       " [ 11.0, -1.0, 1.0],"
                                       " [ 10.0, -21.0, 1.0],"
                                       " [-11.0, -21.0, 1.0],"
                                       " [-12.0, -1.0, 1.0],"
                                       " [ 12.0, -2.0, 1.0],"
                                       " [ 11.0, -22.0, 1.0],"
                                       " [ 10.0, -2.0, 1.0],"
                                       " [ -1.0, -2.0, 1.0],"
                                       " [ -2.0, -22.0, 1.0]]");


      // Normalize these points.
      num::Array2D<cmn::Float64> outputPoints;
      num::Array2D<cmn::Float64> transform;
      normalizePointSequence(inputPoints, outputPoints, transform);

      // Warning(xxx): mean and covariance are not normalized. Check this
      // out later.
      //
      // // Check that the normalization worked.
      // num::Array1D<cmn::Float64> meanArray;
      // num::Array2D<cmn::Float64> covarianceArray;
      // num::getMeanAndCovariance(outputPoints, meanArray, covarianceArray);
      // BRICK_TEST_ASSERT(this->isApproximatelyEqual(
      //                   meanArray, num::zeros<cmn:Float64>(3)));
      // BRICK_TEST_ASSERT(this->isApproximatelyEqual(
      //                   covarianceArray,
      //                   num::identity<cmn::Float64>(3, 3)));

      // Check that the transform and output points match.
      BRICK_TEST_ASSERT(outputPoints.rows() == inputPoints.rows());
      BRICK_TEST_ASSERT(outputPoints.columns() == inputPoints.columns());
      BRICK_TEST_ASSERT(transform.rows() == 3);
      BRICK_TEST_ASSERT(transform.columns() == 3);
      for(size_t ii = 0; ii < inputPoints.rows(); ++ii) {
        BRICK_TEST_ASSERT(
          this->isApproximatelyEqual(
            num::matrixMultiply<cmn::Float64>(
              transform, inputPoints.getRow(ii)),
            outputPoints.getRow(ii)));
      }

      // Check that the output points have the desired property.
      num::Array2D<cmn::Float64> uuT =
        num::matrixMultiply<cmn::Float64>(
          outputPoints.transpose(), outputPoints);
      num::Array2D<cmn::Float64> expectedOuterProductSum =
        (static_cast<cmn::Float64>(inputPoints.rows())
         * num::identity<cmn::Float64>(3, 3));
      BRICK_TEST_ASSERT(this->isApproximatelyEqual(uuT, expectedOuterProductSum));
    }


    void
    EightPointAlgorithmTest::
    getTestPoints(std::vector< num::Vector2D<cmn::Float64> >& uVector,
                  std::vector< num::Vector2D<cmn::Float64> >& uPrimeVector,
                  num::Array2D<cmn::Float64>& FOrig)
    {
      // Create a rank-2 homogeneous 2D transform to serve as our
      // fundamental matrix.
      FOrig = num::Array2D<cmn::Float64>("[[1.0, -1.0, 3.0],"
                                   " [1.0, 0.0, -1.0],"
                                   " [1.0, -1.0, 3.0]]");

      // Create a series of "u" points, ostensibly coming from one image.
      uVector.clear();
      uVector.push_back(num::Vector2D<cmn::Float64>(2.0, 2.0));
      uVector.push_back(num::Vector2D<cmn::Float64>(1.0, 2.0));
      uVector.push_back(num::Vector2D<cmn::Float64>(0.0, 2.0));
      uVector.push_back(num::Vector2D<cmn::Float64>(-1.0, 2.0));
      uVector.push_back(num::Vector2D<cmn::Float64>(-2.0, 2.0));
      uVector.push_back(num::Vector2D<cmn::Float64>(2.0, 1.0));
      uVector.push_back(num::Vector2D<cmn::Float64>(1.0, 1.0));
      uVector.push_back(num::Vector2D<cmn::Float64>(0.0, 1.0));
      uVector.push_back(num::Vector2D<cmn::Float64>(-1.0, 1.0));
      uVector.push_back(num::Vector2D<cmn::Float64>(-2.0, 1.0));
      uVector.push_back(num::Vector2D<cmn::Float64>(2.0, 0.0));
      uVector.push_back(num::Vector2D<cmn::Float64>(1.0, 0.0));
      uVector.push_back(num::Vector2D<cmn::Float64>(0.0, 0.0));
      uVector.push_back(num::Vector2D<cmn::Float64>(-1.0, 0.0));
      uVector.push_back(num::Vector2D<cmn::Float64>(-2.0, 0.0));
      uVector.push_back(num::Vector2D<cmn::Float64>(2.0, -1.0));
      uVector.push_back(num::Vector2D<cmn::Float64>(1.0, -1.0));
      uVector.push_back(num::Vector2D<cmn::Float64>(0.0, -1.0));
      uVector.push_back(num::Vector2D<cmn::Float64>(-1.0, -1.0));
      uVector.push_back(num::Vector2D<cmn::Float64>(-2.0, -1.0));
      uVector.push_back(num::Vector2D<cmn::Float64>(2.0, -2.0));
      uVector.push_back(num::Vector2D<cmn::Float64>(1.0, -2.0));
      uVector.push_back(num::Vector2D<cmn::Float64>(0.0, -2.0));
      uVector.push_back(num::Vector2D<cmn::Float64>(-1.0, -2.0));
      uVector.push_back(num::Vector2D<cmn::Float64>(-2.0, -2.0));

      for(size_t ii = 0; ii < uVector.size(); ++ii) {
        uVector[ii] += num::Vector2D<cmn::Float64>(10.0, 20.0);
      }

      // Create a series of matching "u'" points, ostensibly coming
      // from the other image.
      uPrimeVector.resize(uVector.size());
      for(size_t ii = 0; ii < uVector.size(); ++ii) {
        num::Array1D<cmn::Float64> inputPoint(3);
        inputPoint[0] = uVector[ii].x();
        inputPoint[1] = uVector[ii].y();
        inputPoint[2] = 1.0;
        num::Array1D<cmn::Float64> lineCoeffs =
          num::matrixMultiply<cmn::Float64>(FOrig, inputPoint);

        // By definition of Fundamental Matrix, dot_(uPrime,
        // lineCoeffs) must equal 0, where dot_() is a dot product
        // that includes the third element of the 2D homogeneous
        // coords.
        //
        // if
        //
        //   lineCoords == [c_0, c_1, c_2],
        //
        // and
        //
        //   uPrime == [alpha * k_0, (1.0 - alpha) * k_0, 1]
        //
        // then
        //
        //   c_0 * alpha * k_0 + c_1 * (1.0 - alpha) * k_0 + c_2 = 0
        //
        //   k_0 * ((c_0 - c_1) * alpha + c_1) + c_2 = 0
        //
        //   k_0 = -c_2 / ((c_0 - c_1) * alpha + c_1)

        // Choose 0 < alpha < 1, with alpha not being a "clean"
        // fraction.  This makes us less likely to stumble onto an
        // unlucky combination of u and alpha that makes the solution
        // for k_0 singular.
        cmn::Float64 alpha =
          static_cast<cmn::Float64>((ii % (uVector.size() / 2)) + 1.351)
          / (uVector.size() / 2 + 2);

        // Solve for k_0 and set the value of  uPrime.
        cmn::Float64 k_0 = (-1.0 * lineCoeffs[2]
                      / ((lineCoeffs[0] - lineCoeffs[1]) * alpha
                         + lineCoeffs[1]));
        uPrimeVector[ii].setValue(alpha * k_0, (1.0 - alpha) * k_0);

        // Sanity check.
        cmn::Float64 dotProduct = (uPrimeVector[ii].x() * lineCoeffs[0]
                                   + uPrimeVector[ii].y() * lineCoeffs[1]
                                   + lineCoeffs[2]);
        if(dotProduct >= 1.0E-8) {
          // Should never get here.
          BRICK_THROW(brick::common::LogicException,
                      "EightPointAlgorithmTest::getTestPoints()",
                      "Poorly computed corresponding point!");
        }
      }
    }


    bool
    EightPointAlgorithmTest::
    isApproximatelyEqual(const num::Array1D<cmn::Float64>& array0,
                         const num::Array1D<cmn::Float64>& array1)
    {
      if(array0.size() != array1.size()) {
        return false;
      }
      return std::equal(array0.begin(), array0.end(), array1.begin(),
                        ApproximatelyEqualFunctor<cmn::Float64>(1.0E-10));
    }


    bool
    EightPointAlgorithmTest::
    isApproximatelyEqual(const num::Array2D<cmn::Float64>& array0,
                         const num::Array2D<cmn::Float64>& array1)
    {
      if(array0.rows() != array1.rows()) {
        return false;
      }
      if(array0.columns() != array1.columns()) {
        return false;
      }
      return std::equal(array0.begin(), array0.end(), array1.begin(),
                        ApproximatelyEqualFunctor<cmn::Float64>(1.0E-10));
    }

  } // namespace computerVision

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::computerVision::EightPointAlgorithmTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::computerVision::EightPointAlgorithmTest currentTest;

}

#endif
