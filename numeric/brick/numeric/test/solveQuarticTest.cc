/**
***************************************************************************
* @file brick/numeric/test/solveQuarticTest.cpp
*
* Source file tests for solveQuartic functions.
*
* Copyright (C) 2005 - 2009, 2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <stdlib.h>
#include <brick/common/functional.hh>
#include <brick/numeric/solveQuartic.hh>
#include <brick/test/testFixture.hh>

namespace brick {

  namespace numeric {

    class SolveQuarticTest : public brick::test::TestFixture<SolveQuarticTest> {

    public:

      SolveQuarticTest();
      ~SolveQuarticTest() {}

      // void setUp(const std::string& /* testName */) {}
      // void tearDown(const std::string& /* testName */) {}
      void setUp(const std::string&) {}
      void tearDown(const std::string&) {}

      void testSolveQuartic();

    private:

      double m_defaultTolerance;
    
    }; // class SolveQuarticTest


    /* ============== Member Function Definititions ============== */

    SolveQuarticTest::
    SolveQuarticTest()
      : brick::test::TestFixture<SolveQuarticTest>("SolveQuarticTest"),
        m_defaultTolerance(1.0E-12)
    {
      // Register all tests.
      BRICK_TEST_REGISTER_MEMBER(testSolveQuartic);
    }


    void
    SolveQuarticTest::
    testSolveQuartic()
    {
      // Test for (2x + 4)(-x + 7)(3x - 5)(x - 3) = 0
      //   : (-2x^2 + 14x - 4x + 28)(3x - 5)(x - 3) = 0
      //   : (-2x^2 + 10x + 28)(3x - 5)(x - 3) = 0
      //   : (-6x^3 + 40x^2 + 34x - 140)(x - 3) = 0
      //   : -6x^4 + 58x^3 - 86x^2 - 242x + 420 = 0
      //   : x^4 - (29/3)x^3 + (43/3)x^2 + (121/3)*x - 70 = 0
      double c0 = -29.0 / 3.0;
      double c1 = 43.0 / 3.0;
      double c2 = 121.0 / 3.0;
      double c3 = -70.0;

      brick::common::ComplexNumber<double> root0;
      brick::common::ComplexNumber<double> root1;
      brick::common::ComplexNumber<double> root2;
      brick::common::ComplexNumber<double> root3;
      solveQuartic(c0, c1, c2, c3, root0, root1, root2, root3);

      if(root0.getRealPart() > root1.getRealPart()) {
        std::swap(root0, root1);
      }
      if(root0.getRealPart() > root2.getRealPart()) {
        std::swap(root0, root2);
      }
      if(root0.getRealPart() > root3.getRealPart()) {
        std::swap(root0, root3);
      }
      if(root1.getRealPart() > root2.getRealPart()) {
        std::swap(root1, root2);
      }
      if(root1.getRealPart() > root3.getRealPart()) {
        std::swap(root1, root3);
      }
      if(root2.getRealPart() > root3.getRealPart()) {
        std::swap(root2, root3);
      }
      BRICK_TEST_ASSERT(
        approximatelyEqual(root0.getRealPart(), -2.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root0.getImaginaryPart(), 0.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root1.getRealPart(), 5.0 / 3.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root1.getImaginaryPart(), 0.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root2.getRealPart(), 3.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root2.getImaginaryPart(), 0.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root3.getRealPart(), 7.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root3.getImaginaryPart(), 0.0, m_defaultTolerance));

      
      // Test for (2x + 4i)(4x - 8i)(3x + 9)(x - 3) = 0
      //   : (8x^2 + 32)(3x + 9)(x - 3) = 0
      //   : (24x^3 + 72x^2 + 96x + 288)(x - 3) = 0
      //   : 24x^4 - 120x^2 - 864 = 0
      //   : x^4 - 5x^2 - 36 = 0
      c0 = 0.0;
      c1 = -5.0;
      c2 = 0.0;
      c3 = -36;
      solveQuartic(c0, c1, c2, c3, root0, root1, root2, root3);

      double r0Tag = 10.0 * root0.getRealPart() + root0.getImaginaryPart();
      double r1Tag = 10.0 * root1.getRealPart() + root1.getImaginaryPart();
      double r2Tag = 10.0 * root2.getRealPart() + root2.getImaginaryPart();
      double r3Tag = 10.0 * root3.getRealPart() + root3.getImaginaryPart();
      if(r0Tag > r1Tag) {std::swap(root0, root1); std::swap(r0Tag, r1Tag);}
      if(r0Tag > r2Tag) {std::swap(root0, root2); std::swap(r0Tag, r2Tag);}
      if(r0Tag > r3Tag) {std::swap(root0, root3); std::swap(r0Tag, r3Tag);}
      if(r1Tag > r2Tag) {std::swap(root1, root2); std::swap(r1Tag, r2Tag);}
      if(r1Tag > r3Tag) {std::swap(root1, root3); std::swap(r1Tag, r3Tag);}
      if(r2Tag > r3Tag) {std::swap(root2, root3); std::swap(r2Tag, r3Tag);}
      
      BRICK_TEST_ASSERT(
        approximatelyEqual(root0.getRealPart(), -3.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root0.getImaginaryPart(), 0.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root1.getRealPart(), 0.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root1.getImaginaryPart(), -2.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root2.getRealPart(), 0.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root2.getImaginaryPart(), 2.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root3.getRealPart(), 3.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root3.getImaginaryPart(), 0.0, m_defaultTolerance));

      // Test a bunch of other quartic equations.
      for(c0 = -3.0; c0 < 3.5; c0 += 1.0) {
        for(c1 = -3.0; c1 < 3.5; c1 += 1.0) {
          for(c2 = -3.0; c2 < 3.5; c2 += 1.0) {
            for(c3 = -3.0; c3 < 3.5; c3 += 1.0) {
              solveQuartic(c0, c1, c2, c3, root0, root1, root2, root3);

              brick::common::ComplexNumber<double> result0 =
                (root0 * root0 * root0 * root0 + c0 * root0 * root0 * root0
                 + c1 * root0 * root0 + c2 * root0 + c3);
              brick::common::ComplexNumber<double> result1 =
                (root1 * root1 * root1 * root1 + c0 * root1 * root1 * root1
                 + c1 * root1 * root1 + c2 * root1 + c3);
              brick::common::ComplexNumber<double> result2 =
                (root2 * root2 * root2 * root2 + c0 * root2 * root2 * root2
                 + c1 * root2 * root2 + c2 * root2 + c3);
              brick::common::ComplexNumber<double> result3 =
                (root3 * root3 * root3 * root3 + c0 * root3 * root3 * root3
                 + c1 * root3 * root3 + c2 * root3 + c3);

              BRICK_TEST_ASSERT(
                approximatelyEqual(
                  result0.getRealPart(), 0.0, m_defaultTolerance));
              BRICK_TEST_ASSERT(
                approximatelyEqual(
                  result0.getImaginaryPart(), 0.0, m_defaultTolerance));
              BRICK_TEST_ASSERT(
                approximatelyEqual(
                  result1.getRealPart(), 0.0, m_defaultTolerance));
              BRICK_TEST_ASSERT(
                approximatelyEqual(
                  result1.getImaginaryPart(), 0.0, m_defaultTolerance));
              BRICK_TEST_ASSERT(
                approximatelyEqual(
                  result2.getRealPart(), 0.0, m_defaultTolerance));
              BRICK_TEST_ASSERT(
                approximatelyEqual(
                  result2.getImaginaryPart(), 0.0, m_defaultTolerance));
              BRICK_TEST_ASSERT(
                approximatelyEqual(
                  result3.getRealPart(), 0.0, m_defaultTolerance));
              BRICK_TEST_ASSERT(
                approximatelyEqual(
                  result3.getImaginaryPart(), 0.0, m_defaultTolerance));
            }
          }
        }
      }
    }

  } // namespace numeric
  
} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::numeric::SolveQuarticTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::numeric::SolveQuarticTest currentTest;

}

#endif
