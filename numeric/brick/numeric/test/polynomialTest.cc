/**
***************************************************************************
* @file brick/numeric/test/polynomialTest.cc
*
* Source file defining PolynomialTest class.
*
* Copyright (C) 2004-2005 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/common/functional.hh>
#include <brick/numeric/polynomial.hh>
#include <brick/test/testFixture.hh>


namespace brick {

  namespace numeric {

    class PolynomialTest : public brick::test::TestFixture<PolynomialTest> {

    public:

      PolynomialTest();
      ~PolynomialTest() {};

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests of member functions.
      void testConstructor__void();
      void testConstructor__Type();
      void testConstructor__Type__Type();
      void testConstructor__Type__Type__Type();
      void testConstructor__Array1D();
      void testConstructor__Iter__Iter__bool();
      void testConstructor__Polynomial();
      void testGetCoefficients();
      void testGetOrder();
      void testAssignmentOperator();
      void testApplicationOperator();
      void testTimesEqualsOperator();
      void testPlusEqualsOperator();
      void testMinusEqualsOperator();

      // Tests of non-member functions.
      void testOperatorTimes();
      void testOperatorPlus();
      void testOperatorMinus();

    private:

    }; // class PolynomialTest


    /* ============== Member Function Definititions ============== */

    PolynomialTest::
    PolynomialTest()
      : brick::test::TestFixture<PolynomialTest>("PolynomialTest")
    {
      // Register all tests.
      BRICK_TEST_REGISTER_MEMBER(testConstructor__void);
      BRICK_TEST_REGISTER_MEMBER(testConstructor__Type);
      BRICK_TEST_REGISTER_MEMBER(testConstructor__Type__Type);
      BRICK_TEST_REGISTER_MEMBER(testConstructor__Type__Type__Type);
      BRICK_TEST_REGISTER_MEMBER(testConstructor__Array1D);
      BRICK_TEST_REGISTER_MEMBER(testConstructor__Iter__Iter__bool);
      BRICK_TEST_REGISTER_MEMBER(testConstructor__Polynomial);
      BRICK_TEST_REGISTER_MEMBER(testGetCoefficients);
      BRICK_TEST_REGISTER_MEMBER(testGetOrder);
      BRICK_TEST_REGISTER_MEMBER(testAssignmentOperator);
      BRICK_TEST_REGISTER_MEMBER(testApplicationOperator);
      BRICK_TEST_REGISTER_MEMBER(testTimesEqualsOperator);
      BRICK_TEST_REGISTER_MEMBER(testPlusEqualsOperator);
      BRICK_TEST_REGISTER_MEMBER(testMinusEqualsOperator);
      BRICK_TEST_REGISTER_MEMBER(testOperatorTimes);
      BRICK_TEST_REGISTER_MEMBER(testOperatorPlus);
      BRICK_TEST_REGISTER_MEMBER(testOperatorMinus);
    }


    void
    PolynomialTest::
    testConstructor__void()
    {
      Polynomial<double> polynomial;
      for(double xValue = -10.0; xValue < 10.0; xValue += 0.1) {
        BRICK_TEST_ASSERT(approximatelyEqual(polynomial(xValue), 1.0, 1.0E-11));
      }
    }


    void
    PolynomialTest::
    testConstructor__Type()
    {
      double c0 = 2.3;
      Polynomial<double> polynomial(c0);
      for(double xValue = -10.0; xValue < 10.0; xValue += 0.1) {
        double referenceValue = c0;
        BRICK_TEST_ASSERT(
          approximatelyEqual(polynomial(xValue), referenceValue, 1.0E-11));
      }
    }


    void
    PolynomialTest::
    testConstructor__Type__Type()
    {
      double c0 = 2.3;
      double c1 = -1.5;
      Polynomial<double> polynomial(c1, c0);
      for(double xValue = -10.0; xValue < 10.0; xValue += 0.1) {
        double referenceValue = c1 * xValue + c0;
        BRICK_TEST_ASSERT(
          approximatelyEqual(polynomial(xValue), referenceValue, 1.0E-11));
      }
    }


    void
    PolynomialTest::
    testConstructor__Type__Type__Type()
    {
      double c0 = 2.3;
      double c1 = -1.5;
      double c2 = 0.2;
      Polynomial<double> polynomial(c2, c1, c0);
      for(double xValue = -10.0; xValue < 10.0; xValue += 0.1) {
        double referenceValue = c2 * xValue * xValue + c1 * xValue + c0;
        BRICK_TEST_ASSERT(
          approximatelyEqual(polynomial(xValue), referenceValue, 1.0E-11));
      }
    }


    void
    PolynomialTest::
    testConstructor__Array1D()
    {
      Array1D<double> coefficientArray(4);
      coefficientArray[0] = 2.3;
      coefficientArray[1] = -1.5;
      coefficientArray[2] = 0.2;
      coefficientArray[3] = -0.15;
      Polynomial<double> polynomial(coefficientArray);
      for(double xValue = -10.0; xValue < 10.0; xValue += 0.1) {
        double referenceValue =
          coefficientArray[3] * xValue * xValue * xValue
          + coefficientArray[2] * xValue * xValue
          + coefficientArray[1] * xValue
          + coefficientArray[0];
        BRICK_TEST_ASSERT(
          approximatelyEqual(polynomial(xValue), referenceValue, 1.0E-11));
      }
    }


    void
    PolynomialTest::
    testConstructor__Iter__Iter__bool()
    {
#if BRICK_NUMERIC_POLYNOMIAL_TEMPLATED_CONSTRUCTORS
      std::vector<double> coefficientArray(4);
      coefficientArray[0] = 2.3;
      coefficientArray[1] = -1.5;
      coefficientArray[2] = 0.2;
      coefficientArray[3] = -0.15;
      Polynomial<double> polynomial(
        coefficientArray.begin(), coefficientArray.end(), true);
      for(double xValue = -10.0; xValue < 10.0; xValue += 0.1) {
        double referenceValue =
          coefficientArray[3] * xValue * xValue * xValue
          + coefficientArray[2] * xValue * xValue
          + coefficientArray[1] * xValue
          + coefficientArray[0];
        BRICK_TEST_ASSERT(
          approximatelyEqual(polynomial(xValue), referenceValue, 1.0E-11));
      }
#endif /* BRICK_NUMERIC_POLYNOMIAL_TEMPLATED_CONSTRUCTORS */
    }


    void
    PolynomialTest::
    testConstructor__Polynomial()
    {
      double c0 = 2.3;
      double c1 = -1.5;
      double c2 = 0.2;
      Polynomial<double> polynomial0(c2, c1, c0);
      Polynomial<double> polynomial1(polynomial0);

      // Change polynomial0;
      Polynomial<double> polynomial2(c0, c2);
      polynomial0 = polynomial2;

      // Make sure polynomial1 still has the right behavior.
      for(double xValue = -10.0; xValue < 10.0; xValue += 0.1) {
        double referenceValue = c2 * xValue * xValue + c1 * xValue + c0;
        BRICK_TEST_ASSERT(
          approximatelyEqual(polynomial1(xValue), referenceValue, 1.0E-11));
      }
    }


    void
    PolynomialTest::
    testGetCoefficients()
    {
      Array1D<double> coefficientArray(4);
      coefficientArray[0] = 2.3;
      coefficientArray[1] = -1.5;
      coefficientArray[2] = 0.2;
      coefficientArray[3] = -0.15;
      Polynomial<double> polynomial(coefficientArray);
      Array1D<double> coefficientArray2 = polynomial.getCoefficientArray();
      BRICK_TEST_ASSERT(coefficientArray2.size() == coefficientArray.size());
      BRICK_TEST_ASSERT(
        std::equal(coefficientArray.begin(), coefficientArray.end(),
                   coefficientArray2.begin()));
    }


    void
    PolynomialTest::
    testGetOrder()
    {
      BRICK_TEST_ASSERT(Polynomial<double>(1.0).getOrder() == 0);
      BRICK_TEST_ASSERT(Polynomial<double>(1.0, 1.0).getOrder() == 1);
      BRICK_TEST_ASSERT(Polynomial<double>(1.0, 1.0, 1.0).getOrder() == 2);
    }


    void
    PolynomialTest::
    testAssignmentOperator()
    {
      double c0 = 2.3;
      double c1 = -1.5;
      double c2 = 0.2;
      Polynomial<double> polynomial0(c2, c1, c0);
      Polynomial<double> polynomial1;
      polynomial1 = polynomial0;

      // Change polynomial0;
      Polynomial<double> polynomial2(c0, c2);
      polynomial0 = polynomial2;

      // Make sure polynomial1 still has the right behavior.
      for(double xValue = -10.0; xValue < 10.0; xValue += 0.1) {
        double referenceValue = c2 * xValue * xValue + c1 * xValue + c0;
        BRICK_TEST_ASSERT(
          approximatelyEqual(polynomial1(xValue), referenceValue, 1.0E-11));
      }
    }


    void
    PolynomialTest::
    testApplicationOperator()
    {
      // Already handled in constructor tests.
    }


    void
    PolynomialTest::
    testTimesEqualsOperator()
    {
      Polynomial<double> polynomial0(2.0, 3.0, -1.0);
      Polynomial<double> polynomial1(4.0, 5.0);
      polynomial0 *= polynomial1;
      Array1D<double> coefficientArray = polynomial0.getCoefficientArray();
      BRICK_TEST_ASSERT(polynomial0.getOrder() == 3);
      BRICK_TEST_ASSERT(coefficientArray.size() == 4);
      BRICK_TEST_ASSERT(approximatelyEqual(coefficientArray[0], -5.0, 1.0E-10));
      BRICK_TEST_ASSERT(approximatelyEqual(coefficientArray[1], 11.0, 1.0E-10));
      BRICK_TEST_ASSERT(approximatelyEqual(coefficientArray[2], 22.0, 1.0E-10));
      BRICK_TEST_ASSERT(approximatelyEqual(coefficientArray[3], 8.0, 1.0E-10));
    }


    void
    PolynomialTest::
    testPlusEqualsOperator()
    {
      Polynomial<double> polynomial0(2.0, 3.0, -1.0);
      Polynomial<double> polynomial1(4.0, 5.0);
      polynomial0 += polynomial1;
      Array1D<double> coefficientArray = polynomial0.getCoefficientArray();
      BRICK_TEST_ASSERT(polynomial0.getOrder() == 2);
      BRICK_TEST_ASSERT(coefficientArray.size() == 3);
      BRICK_TEST_ASSERT(approximatelyEqual(coefficientArray[0], 4.0, 1.0E-10));
      BRICK_TEST_ASSERT(approximatelyEqual(coefficientArray[1], 7.0, 1.0E-10));
      BRICK_TEST_ASSERT(approximatelyEqual(coefficientArray[2], 2.0, 1.0E-10));
    }


    void
    PolynomialTest::
    testMinusEqualsOperator()
    {
      Polynomial<double> polynomial0(2.0, 3.0, -1.0);
      Polynomial<double> polynomial1(4.0, 5.0);
      polynomial0 -= polynomial1;
      Array1D<double> coefficientArray = polynomial0.getCoefficientArray();
      BRICK_TEST_ASSERT(polynomial0.getOrder() == 2);
      BRICK_TEST_ASSERT(coefficientArray.size() == 3);
      BRICK_TEST_ASSERT(approximatelyEqual(coefficientArray[0], -6.0, 1.0E-10));
      BRICK_TEST_ASSERT(approximatelyEqual(coefficientArray[1], -1.0, 1.0E-10));
      BRICK_TEST_ASSERT(approximatelyEqual(coefficientArray[2], 2.0, 1.0E-10));
    }


    void
    PolynomialTest::
    testOperatorTimes()
    {
      Polynomial<double> polynomial0(2.0, 3.0, -1.0);
      Polynomial<double> polynomial1(4.0, 5.0);
      Polynomial<double> polynomial2 = polynomial0 * polynomial1;

      Array1D<double> coefficientArray = polynomial2.getCoefficientArray();
      BRICK_TEST_ASSERT(polynomial2.getOrder() == 3);
      BRICK_TEST_ASSERT(coefficientArray.size() == 4);
      BRICK_TEST_ASSERT(approximatelyEqual(coefficientArray[0], -5.0, 1.0E-10));
      BRICK_TEST_ASSERT(approximatelyEqual(coefficientArray[1], 11.0, 1.0E-10));
      BRICK_TEST_ASSERT(approximatelyEqual(coefficientArray[2], 22.0, 1.0E-10));
      BRICK_TEST_ASSERT(approximatelyEqual(coefficientArray[3], 8.0, 1.0E-10));

      for(double xValue = -10.0; xValue < 10.0; xValue += 0.1) {
        double referenceValue = polynomial0(xValue) * polynomial1(xValue);
        BRICK_TEST_ASSERT(
          approximatelyEqual(polynomial2(xValue), referenceValue, 1.0E-11));
      }
    }


    void
    PolynomialTest::
    testOperatorPlus()
    {
      Polynomial<double> polynomial0(2.0, 3.0, -1.0);
      Polynomial<double> polynomial1(4.0, 5.0);
      Polynomial<double> polynomial2 = polynomial0 + polynomial1;

      Array1D<double> coefficientArray = polynomial2.getCoefficientArray();
      BRICK_TEST_ASSERT(polynomial2.getOrder() == 2);
      BRICK_TEST_ASSERT(coefficientArray.size() == 3);
      BRICK_TEST_ASSERT(approximatelyEqual(coefficientArray[0], 4.0, 1.0E-10));
      BRICK_TEST_ASSERT(approximatelyEqual(coefficientArray[1], 7.0, 1.0E-10));
      BRICK_TEST_ASSERT(approximatelyEqual(coefficientArray[2], 2.0, 1.0E-10));

      for(double xValue = -10.0; xValue < 10.0; xValue += 0.1) {
        double referenceValue = polynomial0(xValue) + polynomial1(xValue);
        BRICK_TEST_ASSERT(
          approximatelyEqual(polynomial2(xValue), referenceValue, 1.0E-11));
      }
    }


    void
    PolynomialTest::
    testOperatorMinus()
    {
      Polynomial<double> polynomial0(2.0, 3.0, -1.0);
      Polynomial<double> polynomial1(4.0, 5.0);
      Polynomial<double> polynomial2 = polynomial0 - polynomial1;

      Array1D<double> coefficientArray = polynomial2.getCoefficientArray();
      BRICK_TEST_ASSERT(polynomial2.getOrder() == 2);
      BRICK_TEST_ASSERT(coefficientArray.size() == 3);
      BRICK_TEST_ASSERT(approximatelyEqual(coefficientArray[0], -6.0, 1.0E-10));
      BRICK_TEST_ASSERT(approximatelyEqual(coefficientArray[1], -1.0, 1.0E-10));
      BRICK_TEST_ASSERT(approximatelyEqual(coefficientArray[2], 2.0, 1.0E-10));

      for(double xValue = -10.0; xValue < 10.0; xValue += 0.1) {
        double referenceValue = polynomial0(xValue) - polynomial1(xValue);
        BRICK_TEST_ASSERT(
          approximatelyEqual(polynomial2(xValue), referenceValue, 1.0E-11));
      }
    }

  } //  namespace numeric

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::numeric::PolynomialTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::numeric::PolynomialTest currentTest;

}

#endif
