/**
***************************************************************************
* @file brick/numeric/test/minRecorderTest.cpp
*
* Source file defining minRecorderTest class.
*
* Copyright (C) 2005, 2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/numeric/minRecorder.hh>

#include <brick/common/functional.hh>
#include <brick/test/testFixture.hh>

namespace brick {

  namespace numeric {

    class MinRecorderTest : public brick::test::TestFixture<MinRecorderTest> {

    public:

      MinRecorderTest();
      ~MinRecorderTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      void testConstructor();
      void testConstructor__Type__Payload();
      void testTest();
      void testGetMin();
      void testGetMinimum();
      void testGetPayload();
      void testReset();

    private:

    }; // class MinRecorderTest


    /* ============== Member Function Definititions ============== */

    MinRecorderTest::
    MinRecorderTest()
      : brick::test::TestFixture<MinRecorderTest>("MinRecorderTest")
    {
      // Register all tests.
      BRICK_TEST_REGISTER_MEMBER(testConstructor);
      BRICK_TEST_REGISTER_MEMBER(testConstructor__Type__Payload);
      BRICK_TEST_REGISTER_MEMBER(testTest);
      BRICK_TEST_REGISTER_MEMBER(testGetMin);
      BRICK_TEST_REGISTER_MEMBER(testGetMinimum);
      BRICK_TEST_REGISTER_MEMBER(testGetPayload);
      BRICK_TEST_REGISTER_MEMBER(testReset);
    }


    void
    MinRecorderTest::
    testConstructor()
    {
      MinRecorder<int, unsigned char> minRecorder0;
      BRICK_TEST_ASSERT(minRecorder0.getMinimum() ==
                        std::numeric_limits<int>::max());
      BRICK_TEST_ASSERT(minRecorder0.getPayload() ==
                        static_cast<unsigned char>(0));

      MinRecorder<unsigned int, int> minRecorder1;
      BRICK_TEST_ASSERT(minRecorder1.getMinimum() ==
                        std::numeric_limits<unsigned int>::max());
      BRICK_TEST_ASSERT(minRecorder1.getPayload() == 0);
    }


    void
    MinRecorderTest::
    testConstructor__Type__Payload()
    {
      for(int ii = -10; ii < 10; ++ii) {
        for(unsigned int jj = 0; jj < 20; ++jj) {
          MinRecorder<int, unsigned int> minRecorder0(ii, jj);
          BRICK_TEST_ASSERT(minRecorder0.getMinimum() == ii);
          BRICK_TEST_ASSERT(minRecorder0.getPayload() == jj);
        }
      }
    }


    void
    MinRecorderTest::
    testTest()
    {
      for(int ii = -10; ii < 10; ++ii) {
        for(unsigned int jj = 0; jj < 20; ++jj) {
          int const kBound = 5;
          MinRecorder<int, unsigned int> minRecorder0(ii, jj);
          int kk;
          for(kk = kBound; kk >= -kBound; --kk) {
            unsigned int payload = static_cast<unsigned int>(kk + kBound);
            minRecorder0.test(kk, payload);
            if(kk < ii) {
              BRICK_TEST_ASSERT(minRecorder0.getMin() == kk);
              BRICK_TEST_ASSERT(minRecorder0.getMinimum() == kk);
              BRICK_TEST_ASSERT(minRecorder0.getPayload() == payload);
            } else {
              BRICK_TEST_ASSERT(minRecorder0.getMin() == ii);
              BRICK_TEST_ASSERT(minRecorder0.getMinimum() == ii);
              BRICK_TEST_ASSERT(minRecorder0.getPayload() == jj);
            }
          }
          for(kk = -kBound; kk <= kBound; ++kk) {
            unsigned int payload = static_cast<unsigned int>(kk + 2 * kBound);
            minRecorder0.test(kk, payload);
            if(-kBound < ii) {
              BRICK_TEST_ASSERT(minRecorder0.getMin() == -kBound);
              BRICK_TEST_ASSERT(minRecorder0.getMinimum() == -kBound);
              BRICK_TEST_ASSERT(minRecorder0.getPayload() == 0);
            } else {
              BRICK_TEST_ASSERT(minRecorder0.getMin() == ii);
              BRICK_TEST_ASSERT(minRecorder0.getMinimum() == ii);
              BRICK_TEST_ASSERT(minRecorder0.getPayload() == jj);
            }
          }
        }
      }
    }


    void
    MinRecorderTest::
    testGetMin()
    {
      // Already tested in testTest().
    }


    void
    MinRecorderTest::
    testGetMinimum()
    {
      // Already tested in testTest().
    }


    void
    MinRecorderTest::
    testGetPayload()
    {
      // Already tested in testTest().
    }


    void
    MinRecorderTest::
    testReset()
    {
      for(unsigned int ii = 0; ii < 20; ++ii) {
        for(unsigned int jj = 0; jj < 20; ++jj) {
          MinRecorder<int, unsigned int> minRecorder0(int(ii), jj);
          minRecorder0.reset();
          BRICK_TEST_ASSERT(minRecorder0.getMinimum() ==
                            std::numeric_limits<int>::max());
          BRICK_TEST_ASSERT(minRecorder0.getPayload() ==
                            static_cast<unsigned int>(0));

          MinRecorder<unsigned int, int> minRecorder1(ii, int(jj));
          minRecorder1.reset();
          BRICK_TEST_ASSERT(minRecorder1.getMinimum() ==
                            std::numeric_limits<unsigned int>::max());
          BRICK_TEST_ASSERT(minRecorder1.getPayload() == 0);
        }
      }
    }

  } // namespace numeric

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::numeric::MinRecorderTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::numeric::MinRecorderTest currentTest;

}

#endif
