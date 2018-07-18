/**
***************************************************************************
* @file brick/numeric/test/normalizedCorrelatorTest.cc
*
* Source file defining NormalizedCorrelatorTest class.
*
* Copyright (C) 2004-2005,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/numeric/normalizedCorrelator.hh>
#include <brick/numeric/utilities.hh>
#include <brick/test/testFixture.hh>


namespace brick {

  namespace numeric {

    class NormalizedCorrelatorTest
      : public brick::test::TestFixture<NormalizedCorrelatorTest> {

    public:

      NormalizedCorrelatorTest();
      ~NormalizedCorrelatorTest() {};

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests of member functions.
      void testConstructor__void();
      void testConstructor__IterType0__IterType0__IterType1();
      void testAddSamples();
      void testGetCount();
      void testGetNormalizedCorrelation();
      void testRemoveSamples();

    private:

      double m_defaultTolerance;

    }; // class NormalizedCorrelatorTest


    /* ============== Member Function Definititions ============== */

    NormalizedCorrelatorTest::
    NormalizedCorrelatorTest()
      : brick::test::TestFixture<NormalizedCorrelatorTest>("NormalizedCorrelatorTest"),
        m_defaultTolerance(1.0E-12)
    {
      // Register all tests.
      BRICK_TEST_REGISTER_MEMBER(testConstructor__void);
      BRICK_TEST_REGISTER_MEMBER(
        testConstructor__IterType0__IterType0__IterType1);
      BRICK_TEST_REGISTER_MEMBER(testAddSamples);
      BRICK_TEST_REGISTER_MEMBER(testGetCount);
      BRICK_TEST_REGISTER_MEMBER(testGetNormalizedCorrelation);
      BRICK_TEST_REGISTER_MEMBER(testRemoveSamples);
    }


    void
    NormalizedCorrelatorTest::
    testConstructor__void()
    {
      NormalizedCorrelator<double> normalizedCorrelator;
      BRICK_TEST_ASSERT(normalizedCorrelator.getCount() == 0);
      BRICK_TEST_ASSERT(normalizedCorrelator.getNormalizedCorrelation() == 1.0);
    }


    void
    NormalizedCorrelatorTest::
    testConstructor__IterType0__IterType0__IterType1()
    {
      Array1D<double> inputArray0("[1, 2, 0, 4, 0, 5, -1, 0]");
      Array1D<double> inputArray1("[7, 9, 2, 4, 0, 3,  6, 7]");
      NormalizedCorrelator<double> normalizedCorrelator(
        inputArray0.begin(), inputArray0.end(), inputArray1.begin());
      double testValue = normalizedCorrelator.getNormalizedCorrelation();
      double referenceValue = normalizedCorrelation<double>(
        inputArray0, inputArray1);
      BRICK_TEST_ASSERT(normalizedCorrelator.getCount() == inputArray0.size());
      BRICK_TEST_ASSERT(
        approximatelyEqual(testValue, referenceValue, m_defaultTolerance));
    }


    void
    NormalizedCorrelatorTest::
    testAddSamples()
    {
      Array1D<double> inputArray0a("[1, 2, 0, 4, 0, 5, -1, 0]");
      Array1D<double> inputArray1a("[7, 9, 2, 4, 0, 3,  6, 7]");
      Array1D<double> inputArray0b("[3, -1, -1, 11, 18]");
      Array1D<double> inputArray1b("[6, 10, -2, -9, 12]");
      Array1D<double> inputArray0(inputArray0a.size() + inputArray0b.size());
      Array1D<double> inputArray1(inputArray1a.size() + inputArray1b.size());
      std::copy(inputArray0a.begin(), inputArray0a.end(), inputArray0.begin());
      std::copy(inputArray0b.begin(), inputArray0b.end(),
                inputArray0.begin() + inputArray0a.size());
      std::copy(inputArray1a.begin(), inputArray1a.end(), inputArray1.begin());
      std::copy(inputArray1b.begin(), inputArray1b.end(),
                inputArray1.begin() + inputArray1a.size());


      NormalizedCorrelator<double> normalizedCorrelator(
        inputArray0a.begin(), inputArray0a.end(), inputArray1a.begin());
      normalizedCorrelator.addSamples(
        inputArray0b.begin(), inputArray0b.end(), inputArray1b.begin());

      double testValue = normalizedCorrelator.getNormalizedCorrelation();
      double referenceValue = normalizedCorrelation<double>(
        inputArray0, inputArray1);
      BRICK_TEST_ASSERT(normalizedCorrelator.getCount() == inputArray0.size());
      BRICK_TEST_ASSERT(
        approximatelyEqual(testValue, referenceValue, m_defaultTolerance));
    }


    void
    NormalizedCorrelatorTest::
    testGetCount()
    {
      // Already tested by constructor tests.
    }


    void
    NormalizedCorrelatorTest::
    testGetNormalizedCorrelation()
    {
      // Already tested by constructor tests.
    }


    void
    NormalizedCorrelatorTest::
    testRemoveSamples()
    {
      Array1D<double> inputArray0a("[1, 2, 0, 4, 0, 5, -1, 0]");
      Array1D<double> inputArray1a("[7, 9, 2, 4, 0, 3,  6, 7]");
      Array1D<double> inputArray0b("[3, -1, -1, 11, 18]");
      Array1D<double> inputArray1b("[6, 10, -2, -9, 12]");
      Array1D<double> inputArray0(inputArray0a.size() + inputArray0b.size());
      Array1D<double> inputArray1(inputArray1a.size() + inputArray1b.size());
      std::copy(inputArray0a.begin(), inputArray0a.end(), inputArray0.begin());
      std::copy(inputArray0b.begin(), inputArray0b.end(),
                inputArray0.begin() + inputArray0a.size());
      std::copy(inputArray1a.begin(), inputArray1a.end(), inputArray1.begin());
      std::copy(inputArray1b.begin(), inputArray1b.end(),
                inputArray1.begin() + inputArray1a.size());


      NormalizedCorrelator<double> normalizedCorrelator(
        inputArray0.begin(), inputArray0.end(), inputArray1.begin());
      normalizedCorrelator.removeSamples(
        inputArray0a.begin(), inputArray0a.end(), inputArray1a.begin());

      double testValue = normalizedCorrelator.getNormalizedCorrelation();
      double referenceValue = normalizedCorrelation<double>(inputArray0b, inputArray1b);
      BRICK_TEST_ASSERT(normalizedCorrelator.getCount() == inputArray0b.size());
      BRICK_TEST_ASSERT(
        approximatelyEqual(testValue, referenceValue, m_defaultTolerance));
    }

  } // namespace numeric

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::numeric::NormalizedCorrelatorTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::numeric::NormalizedCorrelatorTest currentTest;

}

#endif
