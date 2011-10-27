/**
***************************************************************************
* @file maxRecorderTest.cpp
* 
* Source file defining maxRecorderTest class.
*
* Copyright (C) 2005, 2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/numeric/maxRecorder.hh>

#include <brick/common/functional.hh>
#include <brick/test/testFixture.hh>

namespace brick {

  namespace numeric {

    class MaxRecorderTest : public brick::test::TestFixture<MaxRecorderTest> {

    public:

      MaxRecorderTest();
      ~MaxRecorderTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      void testConstructor();
      void testConstructor__Type__Payload();
      void testTest();
      void testGetMax();
      void testGetMaximum();
      void testGetPayload();
      void testReset();

    private:
    
    }; // class MaxRecorderTest


    /* ============== Member Function Definititions ============== */

    MaxRecorderTest::
    MaxRecorderTest()
      : brick::test::TestFixture<MaxRecorderTest>("MaxRecorderTest")
    {
      // Register all tests.
      BRICK_TEST_REGISTER_MEMBER(testConstructor);
      BRICK_TEST_REGISTER_MEMBER(testConstructor__Type__Payload);
      BRICK_TEST_REGISTER_MEMBER(testTest);
      BRICK_TEST_REGISTER_MEMBER(testGetMax);
      BRICK_TEST_REGISTER_MEMBER(testGetMaximum);
      BRICK_TEST_REGISTER_MEMBER(testGetPayload);
      BRICK_TEST_REGISTER_MEMBER(testReset);
    }


    void
    MaxRecorderTest::
    testConstructor()
    {
      MaxRecorder<int, unsigned char> maxRecorder0;
      BRICK_TEST_ASSERT(maxRecorder0.getMaximum()
                        == -std::numeric_limits<int>::max());
      BRICK_TEST_ASSERT(maxRecorder0.getPayload()
                        == static_cast<unsigned char>(0));

      MaxRecorder<unsigned int, int> maxRecorder1;
      BRICK_TEST_ASSERT(maxRecorder1.getMaximum() == 0);
      BRICK_TEST_ASSERT(maxRecorder1.getPayload() == 0);      
    }

  
    void
    MaxRecorderTest::
    testConstructor__Type__Payload()
    {
      for(int ii = -10; ii < 10; ++ii) {
        for(unsigned int jj = 0; jj < 20; ++jj) {
          MaxRecorder<int, unsigned int> maxRecorder0(ii, jj);
          BRICK_TEST_ASSERT(maxRecorder0.getMaximum() == ii);
          BRICK_TEST_ASSERT(maxRecorder0.getPayload() == jj);
        }
      }
    }


    void
    MaxRecorderTest::
    testTest()
    {
      for(int ii = -10; ii < 10; ++ii) {
        for(unsigned int jj = 0; jj < 20; ++jj) {
          int const kBound = 5;
          MaxRecorder<int, unsigned int> maxRecorder0(ii, jj);
          int kk;
          for(kk = -kBound; kk <= kBound; ++kk) {
            unsigned int payload = static_cast<unsigned int>(kk + kBound);
            maxRecorder0.test(kk, payload);
            if(kk > ii) {
              BRICK_TEST_ASSERT(maxRecorder0.getMax() == kk);
              BRICK_TEST_ASSERT(maxRecorder0.getMaximum() == kk);
              BRICK_TEST_ASSERT(maxRecorder0.getPayload() == payload);
            } else {
              BRICK_TEST_ASSERT(maxRecorder0.getMax() == ii);
              BRICK_TEST_ASSERT(maxRecorder0.getMaximum() == ii);
              BRICK_TEST_ASSERT(maxRecorder0.getPayload() == jj);
            }            
          }
          for(kk = kBound; kk >= -kBound; --kk) {
            unsigned int payload = static_cast<unsigned int>(kk + 2 * kBound);
            maxRecorder0.test(kk, payload);
            if(kBound > ii) {
              BRICK_TEST_ASSERT(maxRecorder0.getMax() == kBound);
              BRICK_TEST_ASSERT(maxRecorder0.getMaximum() == kBound);
              BRICK_TEST_ASSERT(maxRecorder0.getPayload() == 2 * kBound);
            } else {
              BRICK_TEST_ASSERT(maxRecorder0.getMax() == ii);
              BRICK_TEST_ASSERT(maxRecorder0.getMaximum() == ii);
              BRICK_TEST_ASSERT(maxRecorder0.getPayload() == jj);
            }            
          }
        }
      }
    }

  
    void
    MaxRecorderTest::
    testGetMax()
    {
      // Already tested in testTest().
    }

  
    void
    MaxRecorderTest::
    testGetMaximum()
    {
      // Already tested in testTest().
    }

  
    void
    MaxRecorderTest::
    testGetPayload()
    {
      // Already tested in testTest().
    }

  
    void
    MaxRecorderTest::
    testReset()
    {
      for(unsigned int ii = 0; ii < 20; ++ii) {
        for(unsigned int jj = 0; jj < 20; ++jj) {
          MaxRecorder<int, unsigned int> maxRecorder0(int(ii), jj);
          maxRecorder0.reset();
          BRICK_TEST_ASSERT(maxRecorder0.getMaximum()
                            == -std::numeric_limits<int>::max());
          BRICK_TEST_ASSERT(maxRecorder0.getPayload()
                            == static_cast<unsigned int>(0));

          MaxRecorder<unsigned int, int> maxRecorder1(ii, int(jj));
          maxRecorder1.reset();
          BRICK_TEST_ASSERT(maxRecorder1.getMaximum() == 0);
          BRICK_TEST_ASSERT(maxRecorder1.getPayload() == 0);      
        }
      }
    }

  } // namespace numeric
  
} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::numeric::MaxRecorderTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::numeric::MaxRecorderTest currentTest;

}

#endif
