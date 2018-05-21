/**
***************************************************************************
* @file brick/numeric/test/sampledFunctionsTest.cpp
* 
* Source file defining sampledFunctionsTest class.
*
* Copyright (C) 2006-2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/numeric/numericTraits.hh>
#include <brick/numeric/sampledFunctions.hh>
#include <brick/numeric/utilities.hh>

#include <brick/common/functional.hh>
#include <brick/test/testFixture.hh>

namespace brick {

  namespace numeric {

    class SampledFunctionsTest
      : public brick::test::TestFixture<SampledFunctionsTest> {

    public:

      SampledFunctionsTest();
      ~SampledFunctionsTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      void testGetGaussian1D();
      void testGetHammingWindow1D();

    private:

      double m_defaultTolerance;
    
    }; // class SampledFunctionsTest


    /* ============== Member Function Definititions ============== */

    SampledFunctionsTest::
    SampledFunctionsTest()
      : brick::test::TestFixture<SampledFunctionsTest>("SampledFunctionsTest"),
        m_defaultTolerance(1.0E-10)
    {
      // Register all tests.
      BRICK_TEST_REGISTER_MEMBER(testGetGaussian1D);
      BRICK_TEST_REGISTER_MEMBER(testGetHammingWindow1D);
    }


    void
    SampledFunctionsTest::
    testGetGaussian1D()
    {
      double const sigma = 1.1;
      size_t kernelSize = 7;
      Array1D<double> kernel = getGaussian1D<double>(sigma);
      BRICK_TEST_ASSERT(kernel.size() == kernelSize);
      for(size_t ii = 0; ii < kernelSize; ++ii) {
        double x = ii - ((kernelSize - 1) / 2.0);
        double exponent = -(x * x) / (2.0 * sigma * sigma);
        double referenceValue =
          (1.0 / (std::sqrt(2 * brick::common::constants::pi) * sigma))
          * std::exp(exponent);
        BRICK_TEST_ASSERT(
          approximatelyEqual(kernel[ii], referenceValue, m_defaultTolerance));
      }

      kernelSize = 9;
      kernel = getGaussian1D<double>(sigma, kernelSize);
      BRICK_TEST_ASSERT(kernel.size() == kernelSize);
      double sum = 0.0;
      for(size_t ii = 0; ii < kernelSize; ++ii) {
        double x = ii - ((kernelSize - 1) / 2.0);
        double exponent = -(x * x) / (2.0 * sigma * sigma);
        double referenceValue =
          (1.0 / (std::sqrt(2 * brick::common::constants::pi) * sigma))
          * std::exp(exponent);
        BRICK_TEST_ASSERT(
          approximatelyEqual(kernel[ii], referenceValue, m_defaultTolerance));
        sum += kernel[ii];
      }

      kernel = getGaussian1D<double>(sigma, kernelSize, true);
      BRICK_TEST_ASSERT(kernel.size() == kernelSize);
      for(size_t ii = 0; ii < kernelSize; ++ii) {
        double x = ii - ((kernelSize - 1) / 2.0);
        double exponent = -(x * x) / (2.0 * sigma * sigma);
        double referenceValue =
          ((1.0 / (std::sqrt(2 * brick::common::constants::pi) * sigma))
           * std::exp(exponent) / sum);
        BRICK_TEST_ASSERT(
          approximatelyEqual(kernel[ii], referenceValue, m_defaultTolerance));
      }
    }


    void
    SampledFunctionsTest::
    testGetHammingWindow1D()
    {
      constexpr std::size_t windowSize = 10;
      Array1D<double> hamming = getHammingWindow1D<double>(windowSize);

      // Minimal test just makes sure the window is unimodal in the
      // center.
      std::size_t ii = 1;
      while(ii < (hamming.size() / 2)) {
        BRICK_TEST_ASSERT(hamming[ii] > hamming[ii - 1]);
        ++ii;
      }

      // Even windows will have two identical values in the middle.
      if(windowSize % 2 == 0) {
        double difference = hamming[ii] - hamming[ii - 1];
        BRICK_TEST_ASSERT(absoluteValue(difference)
                          < NumericTraits<double>::epsilon());
        ++ii;
      }

      // Now down the slope to the other side.
      while(ii < hamming.size()) {
        BRICK_TEST_ASSERT(hamming[ii] < hamming[ii - 1]);
        ++ii;
      }
    }

  } // namespace numeric
  
} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::numeric::SampledFunctionsTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::numeric::SampledFunctionsTest currentTest;

}

#endif
