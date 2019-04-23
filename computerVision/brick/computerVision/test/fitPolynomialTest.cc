/**
***************************************************************************
* @file brick/computerVision/fitPolynomialTest.cc
*
* Source file defining tests for fitPolynomial().
*
* Copyright (C) 2019 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/


#include <brick/computerVision/fitPolynomial.hh>
#include <brick/test/functors.hh>
#include <brick/test/testFixture.hh>

// Anonymous namespace for locally defined functions.
namespace {

} // namespace


namespace brick {

  namespace computerVision {

    class FitPolynomialTest
      : public brick::test::TestFixture<FitPolynomialTest> {

    public:

      FitPolynomialTest();
      ~FitPolynomialTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testFitPolynomial_order0();
      void testFitPolynomial_order1();
      void testFitPolynomial_order2();
      void testFitPolynomial_order3();
      void testFitPolynomial_order4();
      void testFitPolynomial_order5();

    private:

      template<class FloatType>
      void
      testFitPolynomial_arbitraryOrder(
        std::vector<FloatType> const& coefficients);
      
      brick::common::Float64 m_defaultTolerance;

    }; // class FitPolynomialTest


    /* ============== Member Function Definititions ============== */

    FitPolynomialTest::
    FitPolynomialTest()
      : brick::test::TestFixture<FitPolynomialTest>("FitPolynomialTest"),
        m_defaultTolerance(1.0E-5)
    {
      BRICK_TEST_REGISTER_MEMBER(testFitPolynomial_order0);
      BRICK_TEST_REGISTER_MEMBER(testFitPolynomial_order1);
      BRICK_TEST_REGISTER_MEMBER(testFitPolynomial_order2);
      BRICK_TEST_REGISTER_MEMBER(testFitPolynomial_order3);
      BRICK_TEST_REGISTER_MEMBER(testFitPolynomial_order4);
      BRICK_TEST_REGISTER_MEMBER(testFitPolynomial_order5);
    }


    void
    FitPolynomialTest::
    testFitPolynomial_order0()
    {
      using brick::common::Float64;
      using brick::common::Int32;
      using brick::numeric::Array1D;
      using brick::numeric::Polynomial;
      using brick::test::approximatelyEqual;
      
      // Polynomial of order 0 should just compute the mean of the Y
      // coordinates.
      std::vector<Float64> xCoordinates;
      std::vector<Float64> yCoordinates;
      Float64 constexpr meanValue{2.0};
      Int32 constexpr minX{-3};
      Int32 constexpr maxX{3};
      for(Int32 ii = minX; ii <= maxX; ++ii) {
        xCoordinates.emplace_back(Float64(ii));
        yCoordinates.emplace_back(Float64(ii + meanValue));
      }

      Polynomial<Float64> polynomial = fitPolynomial<Float64>(
        xCoordinates.begin(), xCoordinates.end(), yCoordinates.begin(), 0);

      BRICK_TEST_ASSERT(polynomial.getOrder() == 0);

      Array1D<Float64> coefficients = polynomial.getCoefficientArray();
      BRICK_TEST_ASSERT(coefficients.size() == 1);
      BRICK_TEST_ASSERT(approximatelyEqual(coefficients[0], meanValue,
                                           this->m_defaultTolerance));

      for(Int32 ii = minX; ii <= maxX; ++ii) {
        BRICK_TEST_ASSERT(approximatelyEqual(polynomial(xCoordinates[ii]),
                                             meanValue,
                                             this->m_defaultTolerance));
      }
    }


    void
    FitPolynomialTest::
    testFitPolynomial_order1()
    {
      using brick::common::Float64;
      std::vector<Float64> const coefficients {-2.0, 3.0};
      testFitPolynomial_arbitraryOrder(coefficients);
    }      


    void
    FitPolynomialTest::
    testFitPolynomial_order2()
    {
      using brick::common::Float64;
      std::vector<Float64> const coefficients {1.5, -2.0, 2.0};
      testFitPolynomial_arbitraryOrder(coefficients);
    }


    void
    FitPolynomialTest::
    testFitPolynomial_order3()
    {
      using brick::common::Float64;
      std::vector<Float64> const coefficients {9.0, 1.0, -2.0, 2.0};
      testFitPolynomial_arbitraryOrder(coefficients);
    }


    void
    FitPolynomialTest::
    testFitPolynomial_order4()
    {
      using brick::common::Float64;
      std::vector<Float64> const coefficients {9.0, 1.0, -2.0, 2.0, 1.0};
      testFitPolynomial_arbitraryOrder(coefficients);
    }

    void
    FitPolynomialTest::
    testFitPolynomial_order5()
    {
      using brick::common::Float64;
      std::vector<Float64> const coefficients {9.0, 1.0, -2.0, 2.0, 1.0, -2.0};
      testFitPolynomial_arbitraryOrder(coefficients);
    }


    template<class FloatType>
    void
    FitPolynomialTest::
    testFitPolynomial_arbitraryOrder(std::vector<FloatType> const& coefficients)
    {
      using brick::common::Int32;
      using brick::numeric::Array1D;
      using brick::numeric::Polynomial;
      using brick::test::approximatelyEqual;
      
      std::vector<FloatType> xCoordinates;
      std::vector<FloatType> yCoordinates;
      Int32 constexpr minX{-3};
      Int32 constexpr maxX{3};
      for(Int32 xx = minX; xx <= maxX; ++xx) {
        FloatType xValue(xx);
        xCoordinates.push_back(xValue);

        FloatType yValue {0};
        FloatType xToTheN {1.0};
        for(Int32 jj = 0; jj < static_cast<Int32>(coefficients.size()); ++jj) {
          yValue += coefficients[jj] * xToTheN;
          xToTheN *= xValue;
        }
        yCoordinates.push_back(yValue);
      }

      Polynomial<FloatType> polynomial = fitPolynomial<FloatType>(
        xCoordinates.begin(), xCoordinates.end(), yCoordinates.begin(),
        coefficients.size() - 1);
      
      BRICK_TEST_ASSERT(polynomial.getOrder() == coefficients.size() - 1); 
      for(Int32 ii = 0; ii < static_cast<int>(xCoordinates.size()); ++ii) {
        BRICK_TEST_ASSERT(approximatelyEqual(polynomial(xCoordinates[ii]),
                                             yCoordinates[ii],
                                             this->m_defaultTolerance));
      }
    }
  
  } // namespace computerVision

} // namespace brick


#if 1

int main(int /* argc */, char** /* argv */)
{
  brick::computerVision::FitPolynomialTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::computerVision::FitPolynomialTest currentTest;

}

#endif
