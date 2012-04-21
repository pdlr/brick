/**
***************************************************************************
* @file brick/computerVision/nChooseKSampleSelectorTest.cc
*
* Source file defining tests for NChooseKSampleSelector.
*
* Copyright (C) 2008,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/common/exception.hh>
#include <brick/computerVision/nChooseKSampleSelector.hh>
#include <brick/numeric/array2D.hh>
#include <brick/test/testFixture.hh>

namespace brick {

  namespace computerVision {
    
    class NChooseKSampleSelectorTest
      : public brick::test::TestFixture<NChooseKSampleSelectorTest> {

    public:

      NChooseKSampleSelectorTest();
      ~NChooseKSampleSelectorTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testGetNumberOfSamples();
      void testGetPool();
      void testGetPoolSize();
      void testGetSample();
      
    private:

      std::vector<int>
      getTestVector(size_t numberOfElements);

    }; // class NChooseKSampleSelectorTest


    /* ============== Member Function Definititions ============== */

    NChooseKSampleSelectorTest::
    NChooseKSampleSelectorTest()
      : TestFixture<NChooseKSampleSelectorTest>(
        "NChooseKSampleSelectorTest")
    {
      BRICK_TEST_REGISTER_MEMBER(testGetNumberOfSamples);
      BRICK_TEST_REGISTER_MEMBER(testGetPool);
      BRICK_TEST_REGISTER_MEMBER(testGetPoolSize);
      BRICK_TEST_REGISTER_MEMBER(testGetSample);
    }


    void
    NChooseKSampleSelectorTest::
    testGetNumberOfSamples()
    {
      // Test against numbers we worked out by hand.
      std::vector<int> testVector = this->getTestVector(5);
      NChooseKSampleSelector<int> selector0(
        3, testVector.begin(), testVector.end());
      NChooseKSampleSelector<int> selector1(
        4, testVector.begin(), testVector.end());
      NChooseKSampleSelector<int> selector2(
        5, testVector.begin(), testVector.end());
      
      BRICK_TEST_ASSERT(selector0.getNumberOfSamples() == 10);
      BRICK_TEST_ASSERT(selector1.getNumberOfSamples() == 5);
      BRICK_TEST_ASSERT(selector2.getNumberOfSamples() == 1);
    }

    
    void
    NChooseKSampleSelectorTest::
    testGetPool()
    {
      for(size_t poolSize = 3; poolSize < 8; ++poolSize) {
        std::vector<int> testVector = this->getTestVector(poolSize);
        NChooseKSampleSelector<int> selector(
          3, testVector.begin(), testVector.end());
        testVector.clear();
      

        std::vector<int> referenceVector = this->getTestVector(poolSize);
        NChooseKSampleSelector<int>::SampleSequenceType pool =
          selector.getPool();

        BRICK_TEST_ASSERT(
          static_cast<size_t>(pool.second - pool.first) == poolSize);
        BRICK_TEST_ASSERT(
          std::equal(pool.first, pool.second, referenceVector.begin()));
      }
    }
    

    void
    NChooseKSampleSelectorTest::
    testGetPoolSize()
    {
      for(size_t poolSize = 3; poolSize < 8; ++poolSize) {
        std::vector<int> testVector = this->getTestVector(poolSize);
        NChooseKSampleSelector<int> selector(
          3, testVector.begin(), testVector.end());
        testVector.clear();

        BRICK_TEST_ASSERT(selector.getPoolSize() == poolSize);
      }
    }
    

    void
    NChooseKSampleSelectorTest::
    testGetSample()
    {
      // Test against numbers we worked out by hand.
      std::vector<int> testVector = this->getTestVector(5);
      NChooseKSampleSelector<int> selector0(
        3, testVector.begin(), testVector.end());
      NChooseKSampleSelector<int> selector1(
        4, testVector.begin(), testVector.end());
      NChooseKSampleSelector<int> selector2(
        5, testVector.begin(), testVector.end());

      numeric::Array2D<size_t> referenceArray0("[[0, 1, 2], "
                                               " [0, 1, 3], "
                                               " [0, 1, 4], "
                                               " [0, 2, 3], "
                                               " [0, 2, 4], "
                                               " [0, 3, 4], "
                                               " [1, 2, 3], "
                                               " [1, 2, 4], "
                                               " [1, 3, 4], "
                                               " [2, 3, 4]]");
      numeric::Array2D<size_t> referenceArray1("[[0, 1, 2, 3], "
                                               " [0, 1, 2, 4], "
                                               " [0, 1, 3, 4], "
                                               " [0, 2, 3, 4], "
                                               " [1, 2, 3, 4]]");
      numeric::Array2D<size_t> referenceArray2("[[0, 1, 2, 3, 4]]");

      // Increment by twos for this first test loop, so that we
      // exercise private member function statelessGetSample().
      for(size_t ii = 0; ii < referenceArray0.rows(); ii += 2) {
        NChooseKSampleSelector<int>::SampleSequenceType sample =
          selector0.getSample(ii);
        BRICK_TEST_ASSERT(
          static_cast<size_t>(sample.second - sample.first)
          == referenceArray0.columns());
        BRICK_TEST_ASSERT(
          std::equal(sample.first, sample.second,
                     referenceArray0.getRow(ii).begin()));
      }
      BRICK_TEST_ASSERT_EXCEPTION(
        common::IndexException, selector0.getSample(referenceArray0.rows()));

      for(size_t ii = 0; ii < referenceArray1.rows(); ++ii) {
        NChooseKSampleSelector<int>::SampleSequenceType sample =
          selector1.getSample(ii);
        BRICK_TEST_ASSERT(
          static_cast<size_t>(sample.second - sample.first)
          == referenceArray1.columns());
        BRICK_TEST_ASSERT(
          std::equal(sample.first, sample.second,
                     referenceArray1.getRow(ii).begin()));
      }
      BRICK_TEST_ASSERT_EXCEPTION(
        common::IndexException, selector1.getSample(referenceArray1.rows()));

      for(size_t ii = 0; ii < referenceArray2.rows(); ++ii) {
        NChooseKSampleSelector<int>::SampleSequenceType sample =
          selector2.getSample(ii);
        BRICK_TEST_ASSERT(
          static_cast<size_t>(sample.second - sample.first)
          == referenceArray2.columns());
        BRICK_TEST_ASSERT(
          std::equal(sample.first, sample.second,
                     referenceArray2.getRow(ii).begin()));
      }
      BRICK_TEST_ASSERT_EXCEPTION(
        common::IndexException, selector2.getSample(referenceArray2.rows()));
    }


    std::vector<int>
    NChooseKSampleSelectorTest::
    getTestVector(size_t numberOfElements)
    {
      std::vector<int> result(numberOfElements);
      for(size_t ii = 0; ii < numberOfElements; ++ii) {
        result[ii] = ii;
      }
      return result;
    }
      
    
  } // namespace computerVision

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::computerVision::NChooseKSampleSelectorTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::computerVision::NChooseKSampleSelectorTest currentTest;

}

#endif
