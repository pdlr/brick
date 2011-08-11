/**
***************************************************************************
* @file solveQuadraticTest.cpp
*
* Source file tests for solveQuadratic function templates.
*
* Copyright (C) 2005 - 2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <stdlib.h>
#include <brick/common/functional.hh>
#include <brick/numeric/solveQuadratic.hh>
#include <brick/test/testFixture.hh>

namespace brick {

  namespace numeric {

    class SolveQuadraticTest : public test::TestFixture<SolveQuadraticTest> {

    public:

      SolveQuadraticTest();
      ~SolveQuadraticTest() {}

      // void setUp(const std::string& testName) {}
      // void tearDown(const std::string& testName) {}
      void setUp(const std::string&) {}
      void tearDown(const std::string&) {}

      void testSolveQuadratic__Type_Type_Type_Type_Type();
      void testSolveQuadratic__Type_Type_Type_complex_complex();
      void testSolveQuadratic__complex_complex_complex_complex();

    private:

      double m_defaultTolerance;
    
    }; // class SolveQuadraticTest


    /* ============== Member Function Definititions ============== */

    SolveQuadraticTest::
    SolveQuadraticTest()
      : TestFixture<SolveQuadraticTest>("SolveQuadraticTest"),
        m_defaultTolerance(1.0E-12)
    {
      // Register all tests.
      BRICK_TEST_REGISTER_MEMBER(
        testSolveQuadratic__Type_Type_Type_Type_Type);
      BRICK_TEST_REGISTER_MEMBER(
        testSolveQuadratic__Type_Type_Type_complex_complex);
      BRICK_TEST_REGISTER_MEMBER(
        testSolveQuadratic__complex_complex_complex_complex);
    }


    void
    SolveQuadraticTest::
    testSolveQuadratic__Type_Type_Type_Type_Type()
    {
      // Test for (2x + 4)(-x + 7) = 0;
      double c0 = -2.0;
      double c1 = 10.0;
      double c2 = 28.0;

      double root0;
      double root1;
      bool valid = solveQuadratic(c0, c1, c2, root0, root1);

      if(root0 > root1) {
        std::swap(root0, root1);
      }
      BRICK_TEST_ASSERT(valid == true);
      BRICK_TEST_ASSERT(approximatelyEqual(root0, -2.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(approximatelyEqual(root1, 7.0, m_defaultTolerance));


      // Test for (2x + 4)(-x - 7) = 0;
      c0 = -2.0;
      c1 = -18.0;
      c2 = -28.0;

      valid = solveQuadratic(c0, c1, c2, root0, root1);

      if(root0 > root1) {
        std::swap(root0, root1);
      }
      BRICK_TEST_ASSERT(valid == true);
      BRICK_TEST_ASSERT(approximatelyEqual(root0, -7.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(approximatelyEqual(root1, -2.0, m_defaultTolerance));

      
      // Test for (2x + 4i)(4x - 8i) = 0;
      c0 = 8.0;
      c1 = 0.0;
      c2 = 32.0;

      root0 = 0.0;
      root1 = 0.0;
      valid = solveQuadratic(c0, c1, c2, root0, root1);

      BRICK_TEST_ASSERT(valid == false);
      BRICK_TEST_ASSERT(root0 == 0.0);
      BRICK_TEST_ASSERT(root1 == 0.0);
    }


    void
    SolveQuadraticTest::
    testSolveQuadratic__Type_Type_Type_complex_complex()
    {
      // Test for (2x + 4)(-x + 7) = 0;
      double c0 = -2.0;
      double c1 = 10.0;
      double c2 = 28.0;

      std::complex<double> root0;
      std::complex<double> root1;
      solveQuadratic(c0, c1, c2, root0, root1);

      if(root0.real() > root1.real()) {
        std::swap(root0.real(), root1.real());
      }
      BRICK_TEST_ASSERT(
        approximatelyEqual(root0.real(), -2.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root0.imag(), 0.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root1.real(), 7.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root1.imag(), 0.0, m_defaultTolerance));


      // Test for (2x + 4)(-x - 7) = 0;
      c0 = -2.0;
      c1 = -18.0;
      c2 = -28.0;

      solveQuadratic(c0, c1, c2, root0, root1);

      if(root0.real() > root1.real()) {
        std::swap(root0.real(), root1.real());
      }
      BRICK_TEST_ASSERT(
        approximatelyEqual(root0.real(), -7.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root0.imag(), 0.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root1.real(), -2.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root1.imag(), 0.0, m_defaultTolerance));

      // Test for (2x + 4i)(4x - 8i) = 0;
      c0 = 8.0;
      c1 = 0.0;
      c2 = 32.0;

      solveQuadratic(c0, c1, c2, root0, root1);

      if(root0.imag() > root1.imag()) {
        std::swap(root0, root1);
      }
      BRICK_TEST_ASSERT(
        approximatelyEqual(root0.real(), 0.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root0.imag(), -2.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root1.real(), 0.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root1.imag(), 2.0, m_defaultTolerance));
    }


    void
    SolveQuadraticTest::
    testSolveQuadratic__complex_complex_complex_complex()
    {
      // Test for (2x + 4)(-x + 7) = 0;
      std::complex<double> c0(-2.0, 0.0);
      std::complex<double> c1(10.0, 0.0);
      std::complex<double> c2(28.0, 0.0);

      std::complex<double> root0;
      std::complex<double> root1;
      solveQuadratic(c1 / c0, c2 / c0, root0, root1);

      if(root0.real() > root1.real()) {
        std::swap(root0.real(), root1.real());
      }
      BRICK_TEST_ASSERT(
        approximatelyEqual(root0.real(), -2.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root0.imag(), 0.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root1.real(), 7.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root1.imag(), 0.0, m_defaultTolerance));


      // Test for (2x + 4)(-x - 7) = 0;
      c0 = std::complex<double>(-2.0, 0.0);
      c1 = std::complex<double>(-18.0, 0.0);
      c2 = std::complex<double>(-28.0, 0.0);
      solveQuadratic(c1 / c0, c2 / c0, root0, root1);

      if(root0.real() > root1.real()) {
        std::swap(root0.real(), root1.real());
      }
      BRICK_TEST_ASSERT(
        approximatelyEqual(root0.real(), -7.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root0.imag(), 0.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root1.real(), -2.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root1.imag(), 0.0, m_defaultTolerance));

      // Test for (2x + 4i)(4x - 8i) = 0;
      c0 = std::complex<double>(8.0, 0.0);
      c1 = std::complex<double>(0.0, 0.0);
      c2 = std::complex<double>(32.0, 0.0);
      solveQuadratic(c1 / c0, c2 / c0, root0, root1);

      if(root0.imag() > root1.imag()) {
        std::swap(root0, root1);
      }
      BRICK_TEST_ASSERT(
        approximatelyEqual(root0.real(), 0.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root0.imag(), -2.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root1.real(), 0.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(root1.imag(), 2.0, m_defaultTolerance));
    }
    
  } // namespace numeric
  
} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::numeric::SolveQuadraticTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::numeric::SolveQuadraticTest currentTest;

}

#endif
