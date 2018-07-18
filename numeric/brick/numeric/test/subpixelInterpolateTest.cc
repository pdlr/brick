/**
***************************************************************************
* @file brick/numeric/test/subpixelInterpolateTest.cc
*
* Source file defining tests for the subpixelInterpolate() function
* template.
*
* Copyright (C) 2009,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/common/functional.hh>
#include <brick/numeric/array1D.hh>
#include <brick/numeric/array2D.hh>
#include <brick/numeric/utilities.hh>
#include <brick/numeric/subpixelInterpolate.hh>
#include <brick/test/testFixture.hh>


namespace brick {

  namespace numeric {

    class SubpixelInterpolateTest
      : public brick::test::TestFixture<SubpixelInterpolateTest> {

    public:

      SubpixelInterpolateTest();
      ~SubpixelInterpolateTest() {}

      void setUp(const std::string&) {}
      void tearDown(const std::string&) {}

      // Tests.
      void testGetQuadraticCoefficients();
      void testSubpixelInterpolate__TypeX3_etc();
      void testSubpixelInterpolate__TypeX9_etc();

    private:

#if 0
      // Here's a reference implementation of
      // getQuadraticCoefficients3x3().  Unfortunately, it needs
      // linearAlgebra::pseudoInverse(), so if you want to compile and
      // run it, you'll have to go through a contortion.
      void
      privateGetQuadraticCoefficients3x3(
        double value00, double value01, double value02,
        double value10, double value11, double value12,
        double value20, double value21, double value22,
        double& k0, double& k1, double& k2, double& k3, double& k4, double& k5);
#endif /* if 0 */

      double m_defaultTolerance;

    }; // class SubpixelInterpolateTest


    /* ============== Member Function Definititions ============== */

    SubpixelInterpolateTest::
    SubpixelInterpolateTest()
      : brick::test::TestFixture<SubpixelInterpolateTest>("SubpixelInterpolateTest"),
        m_defaultTolerance(1.0E-10)
    {
      BRICK_TEST_REGISTER_MEMBER(testGetQuadraticCoefficients);
      BRICK_TEST_REGISTER_MEMBER(testSubpixelInterpolate__TypeX3_etc);
      BRICK_TEST_REGISTER_MEMBER(testSubpixelInterpolate__TypeX9_etc);
    }


    void
    SubpixelInterpolateTest::
    testGetQuadraticCoefficients()
    {
      // Make an arbitrary quadratic and fill in a 3x3 array of values.
      Array2D<double> AA("[[2.0, 3.0],"
                         " [3.0, -1.5]]");
      Array1D<double> bb("[2.7, -2.3]");
      double cc = 5.0;

      Array2D<double> myImage(3, 3);
      for(int xx = 0; xx < static_cast<int>(myImage.columns()); ++xx) {
        for(int yy = 0; yy < static_cast<int>(myImage.rows()); ++yy) {
          Array1D<double> xyVec(2);
          xyVec[0] = static_cast<double>(xx - 1);
          xyVec[1] = static_cast<double>(yy - 1);

          myImage(yy, xx) =
            (dot<double>(xyVec, matrixMultiply<double>(AA, xyVec))
             + dot<double>(xyVec, bb) + cc);
        }
      }

      // Try to recover the coefficients of the quadratic.
      double k0, k1, k2, k3, k4, k5;
      getQuadraticCoefficients3x3(
        myImage(0, 0), myImage(0, 1), myImage(0, 2),
        myImage(1, 0), myImage(1, 1), myImage(1, 2),
        myImage(2, 0), myImage(2, 1), myImage(2, 2),
        k0, k1, k2, k3, k4, k5);

      // Check that the estimates match the real coefficients.
      BRICK_TEST_ASSERT(approximatelyEqual(k0, AA(0, 0), m_defaultTolerance));
      BRICK_TEST_ASSERT(approximatelyEqual(k1, AA(0, 1), m_defaultTolerance));
      BRICK_TEST_ASSERT(approximatelyEqual(k2, AA(1, 1), m_defaultTolerance));
      BRICK_TEST_ASSERT(approximatelyEqual(k3, bb(0), m_defaultTolerance));
      BRICK_TEST_ASSERT(approximatelyEqual(k4, bb(1), m_defaultTolerance));
      BRICK_TEST_ASSERT(approximatelyEqual(k5, cc, m_defaultTolerance));
    }



    void
    SubpixelInterpolateTest::
    testSubpixelInterpolate__TypeX3_etc()
    {
      // Set up a quadratic.
      //
      //   f(p) = 2*p^2 - p - 2
      //
      // Solving manually gives extremum at 0.25.
      const double targetPosition = 0.25;
      const double minValue =
        2.0 * targetPosition * targetPosition - targetPosition - 2.0;

      Array1D<double> myImage(10);
      for(int position = 0; position < static_cast<int>(myImage.size());
          ++position) {
        myImage(position) = 2.0 * position * position - position - 2.0;
      }

      // Try to recover the min value and location in several places.
      for(int p0 = 1; p0 < static_cast<int>(myImage.size()) - 1; ++p0) {
        double extremumPosition;
        double extremeValue;
        bool returnValue = subpixelInterpolate<double>(
          static_cast<double>(p0), myImage(p0-1), myImage(p0), myImage(p0+1),
          extremumPosition, extremeValue);

        BRICK_TEST_ASSERT(returnValue == true);
        BRICK_TEST_ASSERT(
          approximatelyEqual(
            extremumPosition, targetPosition, m_defaultTolerance));
        BRICK_TEST_ASSERT(
          approximatelyEqual(extremeValue, minValue, m_defaultTolerance));
      }
    }


    void
    SubpixelInterpolateTest::
    testSubpixelInterpolate__TypeX9_etc()
    {
#if 0
      // Make an array that contains nothing but a quadratic with a
      // minimum of -5.0 at row 3.4, column 7.2.
      const double minValue = -5.0;
      const double targetRow = 3.4;
      const double targetColumn = 7.2;

      Array2D<double> myImage(10, 10);
      for(int row = 0; row < static_cast<int>(myImage.rows()); ++row) {
        for(int column = 0; column < static_cast<int>(myImage.columns());
            ++column) {
          double rowOffset = row - targetRow;
          double columnOffset = column - targetColumn;

          double pixelValue = (minValue
                               + 2.0 * columnOffset * columnOffset
                               + 1.4 * columnOffset * rowOffset);
          myImage(row, column) = pixelValue;
        }
      }
#else
      // Now a harder example:
      //   f(r,c) = [r, c] [2,   1.5] [r] + [r, c] [-1.0] - 2.0
      //                   [1.5, 3.0] [c]          [-2.5]
      //
      // Solving manually gives extremum at (-0.1, 0.46666666666666666]").
      Array2D<double> AA("[[2.0, 1.5],"
                         " [1.5, 3.0]]");
      Array1D<double> bb("[-1.0, -2.5]");
      double cc = -2.0;

      // What we're hoping to recover.
      Array1D<double> rcMin("[-0.1, 0.46666666666666666]");
      const double targetRow = rcMin[0];
      const double targetColumn = rcMin[1];
      const double minValue =
        (dot<double>(rcMin, matrixMultiply<double>(AA, rcMin))
         + dot<double>(rcMin, bb) + cc);
#endif /* if 0 */

      Array2D<double> myImage(10, 10);
      for(int row = 0; row < static_cast<int>(myImage.rows()); ++row) {
        for(int column = 0; column < static_cast<int>(myImage.columns());
            ++column) {
          Array1D<double> rcVec(2);
          rcVec(0) = row;
          rcVec(1) = column;
          double pixelValue =
            (dot<double>(rcVec, matrixMultiply<double>(AA, rcVec))
             + dot<double>(rcVec, bb) + cc);
          myImage(row, column) = pixelValue;
        }
      }

      // Try to recover the min value and location in several places.
      for(int r0 = 1; r0 < static_cast<int>(myImage.rows()) - 1; ++r0) {
        for(int c0 = 1; c0 < static_cast<int>(myImage.columns()) - 1; ++c0) {
          double extremumRow;
          double extremumColumn;
          double extremeValue;
          bool returnValue = subpixelInterpolate<double>(
            static_cast<double>(r0), static_cast<double>(c0),
            myImage(r0-1, c0-1), myImage(r0-1, c0), myImage(r0-1, c0+1),
            myImage(r0, c0-1), myImage(r0, c0), myImage(r0, c0+1),
            myImage(r0+1, c0-1), myImage(r0+1, c0), myImage(r0+1, c0+1),
            extremumRow, extremumColumn, extremeValue);

          BRICK_TEST_ASSERT(returnValue == true);
          BRICK_TEST_ASSERT(
            approximatelyEqual(
              extremumRow, targetRow, m_defaultTolerance));
          BRICK_TEST_ASSERT(
            approximatelyEqual(
              extremumColumn, targetColumn, m_defaultTolerance));
          BRICK_TEST_ASSERT(
            approximatelyEqual(extremeValue, minValue, m_defaultTolerance));
        }
      }
    }


#if 0
    // Here's a reference implementation of
    // getQuadraticCoefficients3x3().  Unfortunately, it needs
    // linearAlgebra::pseudoInverse(), so if you want to compile and
    // run it, you'll have to go through a contortion.
    void
    SubpixelInterpolateTest::
    privateGetQuadraticCoefficients3x3(
      double value00, double value01, double value02,
      double value10, double value11, double value12,
      double value20, double value21, double value22,
      double& k0, double& k1, double& k2, double& k3, double& k4, double& k5)
    {
      Array2D<double> AMatrix(9, 6);
      size_t ii = 0;
      for(int yy = -1; yy < 2; ++yy) {
        for(int xx = -1; xx < 2; ++xx) {
          AMatrix(ii, 0) = xx * xx;
          AMatrix(ii, 1) = 2.0 * xx * yy;
          AMatrix(ii, 2) = yy * yy;
          AMatrix(ii, 3) = xx;
          AMatrix(ii, 4) = yy;
          AMatrix(ii, 5) = 1.0;
          ++ii;
        }
      }

      Array2D<double> APseudoInv = linearAlgebra::pseudoinverse(AMatrix);
      Array1D<double> bVector(9);
      bVector(0) = value00;
      bVector(1) = value01;
      bVector(2) = value02;
      bVector(3) = value10;
      bVector(4) = value11;
      bVector(5) = value12;
      bVector(6) = value20;
      bVector(7) = value21;
      bVector(8) = value22;

      Array1D<double> xVector = matrixMultiply<double>(APseudoInv, bVector);
      k0 = xVector(0);
      k1 = xVector(1);
      k2 = xVector(2);
      k3 = xVector(3);
      k4 = xVector(4);
      k5 = xVector(5);
    }
#endif /* if 0 */

  } // namespace numeric

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::numeric::SubpixelInterpolateTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::numeric::SubpixelInterpolateTest currentTest;

}

#endif
