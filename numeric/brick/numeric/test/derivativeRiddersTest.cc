/**
***************************************************************************
* @file brick/numeric/test/derivativeRiddersTest.cpp
*
* Source file defining tests for numerical differentiation.
*
* Copyright (C) 2009 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <iostream>
#include <brick/numeric/derivativeRidders.hh>
#include <brick/test/testFixture.hh>


namespace brick {

  namespace numeric {

    class DerivativeRiddersTest
      : public test::TestFixture<DerivativeRiddersTest> {

    public:

      DerivativeRiddersTest();
      ~DerivativeRiddersTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testDerivativeRidders();
    
    private:

      double m_defaultTolerance;
    
    }; // class DerivativeRiddersTest


    // This functor implements the function
    // @code
    //   F(x) = exp(x) / (sin(x) - x^2),
    // @endcode
    // which Ridders used as an example in his paper.
    struct RiddersFunctor0
      : public std::unary_function<double, double>
    {
      double operator()(double xx) {
        return std::exp(xx) / (std::sin(xx) - xx * xx);
      }

      double derivative(double xx) {
        double numerator = std::exp(xx);
        double denominator = std::sin(xx) - xx * xx;
        double numeratorPrime = xx * numerator;
        double denominatorPrime = std::cos(xx) - 2 * xx;
        return (numeratorPrime / denominator
                - numerator * denominatorPrime / (denominator * denominator));
      }
    };

    
    /* ============== Member Function Definititions ============== */

    DerivativeRiddersTest::
    DerivativeRiddersTest()
      : test::TestFixture<DerivativeRiddersTest>("DerivativeRiddersTest"),
        m_defaultTolerance(1.0E-10)
    {
      BRICK_TEST_REGISTER_MEMBER(testDerivativeRidders);
    }


    void
    DerivativeRiddersTest::
    testDerivativeRidders()
    {
      // To test the this class, we use the example from Ridders' paper.
      RiddersFunctor0 functor;
      // DerivativeRidders<RiddersFunctor0> derivativeRidders(functor);

      // Duplicate the parameters used in Ridders' paper.
      DerivativeRidders<RiddersFunctor0> derivativeRidders(
        functor, 0.01, 2.0, 5);

      double derivative = functor.derivative(1.0);

      double errorValue;
      double derivativeEst = derivativeRidders.estimateDerivative(
        1.0, errorValue);

      // std::cout << "Error is "
      //           << common::absoluteValue(derivative - derivativeEst)
      //           << " (+/- " << errorValue << ")" << std::endl;
      
      BRICK_TEST_ASSERT(approximatelyEqual(derivative, derivativeEst, 1.0e-9));
      BRICK_TEST_ASSERT(errorValue > 1.0E-14);
      BRICK_TEST_ASSERT(errorValue < 1.0E-7);
    }

  } // namespace numeric
  
} // namespace brick

#if 0

int main(int argc, char** argv)
{
  brick::DerivativeRiddersTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::numeric::DerivativeRiddersTest currentTest;
  
}

#endif
