/**
***************************************************************************
* @file solveCubicTest.cpp
*
* Source file tests for solveCubic functions.
*
* Copyright (C) 2005 - 2009, 2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <stdlib.h>
#include <brick/common/functional.hh>
#include <brick/numeric/solveCubic.hh>
#include <brick/test/testFixture.hh>

namespace brick {

  namespace numeric {

    class SolveCubicTest : public brick::test::TestFixture<SolveCubicTest> {

    public:

      SolveCubicTest();
      ~SolveCubicTest() {}

      // void setUp(const std::string& /* testName */) {}
      // void tearDown(const std::string& /* testName */) {}
      void setUp(const std::string&) {}
      void tearDown(const std::string&) {}

      void testSolveCubic__Type_Type_Type_Type_Type_Type();
      void testSolveCubic__Type_Type_Type_complex_complex_complex();

    private:

      double m_defaultTolerance;
    
    }; // class SolveCubicTest


    /* ============== Member Function Definititions ============== */

    SolveCubicTest::
    SolveCubicTest()
      : TestFixture<SolveCubicTest>("SolveCubicTest"),
        m_defaultTolerance(1.0E-12)
    {
      // Register all tests.
      BRICK_TEST_REGISTER_MEMBER(
        testSolveCubic__Type_Type_Type_Type_Type_Type);
      BRICK_TEST_REGISTER_MEMBER(
        testSolveCubic__Type_Type_Type_complex_complex_complex);
    }


    void
    SolveCubicTest::
    testSolveCubic__Type_Type_Type_Type_Type_Type()
    {
      // Test for (2x + 4)(-x + 7)(3x - 5) = 0
      //   : (-2x^2 + 14x - 4x + 28)(3x - 5) = 0
      //   : (-2x^2 + 10x + 28)(3x - 5) = 0
      //   : -6x^3 + 30x^2 + 84x + 10x^2 - 50x - 140 = 0
      //   : -6x^3 + 40x^2 + 34x - 140 = 0
      //   : x^3 - (20/3)x^2 - (17/3)x + 70/3 = 0
      double c0 = -20.0 / 3.0;
      double c1 = -17.0 / 3.0;
      double c2 = 70.0 / 3.0;

      double root0;
      double root1;
      double root2;
      bool valid = solveCubic(c0, c1, c2, root0, root1, root2);

      if(root0 > root1) {
        std::swap(root0, root1);
      }
      if(root0 > root2) {
        std::swap(root0, root2);
      }
      if(root1 > root2) {
        std::swap(root1, root2);
      }
      BRICK_TEST_ASSERT(valid == true);
      BRICK_TEST_ASSERT(approximatelyEqual(root0, -2.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(approximatelyEqual(root1, 5.0 / 3.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(approximatelyEqual(root2, 7.0, m_defaultTolerance));

      
      // Test for (2x + 4i)(4x - 8i)(3x + 9) = 0
      //   : (8x^2 + 32)(3x + 9) = 0
      //   : 24x^3 + 72x^2 + 96x + 288 = 0
      //   : x^3 + 3x^2 + 4x + 12 = 0
      c0 = 3.0;
      c1 = 4.0;
      c2 = 12.0;

      root1 = 0.0;
      root2 = 0.0;
      valid = solveCubic(c0, c1, c2, root0, root1, root2);

      BRICK_TEST_ASSERT(valid == false);
      BRICK_TEST_ASSERT(approximatelyEqual(root0, -3.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(root1 == 0.0);
      BRICK_TEST_ASSERT(root2 == 0.0);
    }


    void
    SolveCubicTest::
    testSolveCubic__Type_Type_Type_complex_complex_complex()
    {
      // Test for (2x + 4)(-x + 7)(3x - 5) = 0
      //   : (-2x^2 + 14x - 4x + 28)(3x - 5) = 0
      //   : (-2x^2 + 10x + 28)(3x - 5) = 0
      //   : -6x^3 + 30x^2 + 84x + 10x^2 - 50x - 140 = 0
      //   : -6x^3 + 40x^2 + 34x - 140 = 0
      //   : x^3 - (20/3)x^2 - (17/3)x + 70/3 = 0
      double c0 = -20.0 / 3.0;
      double c1 = -17.0 / 3.0;
      double c2 = 70.0 / 3.0;

      std::complex<double> root0;
      std::complex<double> root1;
      std::complex<double> root2;
      solveCubic(c0, c1, c2, root0, root1, root2);

      if(root0.real() > root1.real()) {
        std::swap(root0, root1);
      }
      if(root0.real() > root2.real()) {
        std::swap(root0, root2);
      }
      if(root1.real() > root2.real()) {
        std::swap(root1, root2);
      }
      BRICK_TEST_ASSERT(
        approximatelyEqual(root0.real(), -2.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root0.imag(), 0.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root1.real(), 5.0 / 3.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root1.imag(), 0.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root2.real(), 7.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root2.imag(), 0.0, m_defaultTolerance));

      
      // Test for (2x + 4i)(4x - 8i)(3x + 9) = 0
      //   : (8x^2 + 32)(3x + 9) = 0
      //   : 24x^3 + 72x^2 + 96x + 288 = 0
      //   : x^3 + 3x^2 + 4x + 12 = 0
      c0 = 3.0;
      c1 = 4.0;
      c2 = 12.0;
      solveCubic(c0, c1, c2, root0, root1, root2);

      if(root1.imag() > root2.imag()) {
        std::swap(root1, root2);
      }
      BRICK_TEST_ASSERT(
        approximatelyEqual(root0.real(), -3.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root0.imag(), 0.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root1.real(), 0.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root1.imag(), -2.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root2.real(), 0.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root2.imag(), 2.0, m_defaultTolerance));

      // Test a bunch of other cubic equations.
      for(c0 = -3.0; c0 < 3.5; c0 += 1.0) {
        for(c1 = -3.0; c1 < 3.5; c1 += 1.0) {
          for(c2 = -3.0; c2 < 3.5; c2 += 1.0) {
            solveCubic(c0, c1, c2, root0, root1, root2);

            std::complex<double> result0 =
              root0 * root0 * root0 + c0 * root0 * root0 + c1 * root0 + c2;
            std::complex<double> result1 =
              root1 * root1 * root1 + c0 * root1 * root1 + c1 * root1 + c2;
            std::complex<double> result2 =
              root2 * root2 * root2 + c0 * root2 * root2 + c1 * root2 + c2;

            BRICK_TEST_ASSERT(
              approximatelyEqual(result0.real(), 0.0, m_defaultTolerance));
            BRICK_TEST_ASSERT(
              approximatelyEqual(result0.imag(), 0.0, m_defaultTolerance));
            BRICK_TEST_ASSERT(
              approximatelyEqual(result1.real(), 0.0, m_defaultTolerance));
            BRICK_TEST_ASSERT(
              approximatelyEqual(result1.imag(), 0.0, m_defaultTolerance));
            BRICK_TEST_ASSERT(
              approximatelyEqual(result2.real(), 0.0, m_defaultTolerance));
            BRICK_TEST_ASSERT(
              approximatelyEqual(result2.imag(), 0.0, m_defaultTolerance));
          }
        }
      }
    }

  } // namespace numeric
  
} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::numeric::SolveCubicTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::numeric::SolveCubicTest currentTest;

}

#endif
