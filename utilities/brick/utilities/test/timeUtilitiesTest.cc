/**
***************************************************************************
* @file brick/utilities/test/timeUtilitiesTest.cpp
*
* Source file defining TimeUtilitiesTest class.
*
* Copyright (C) 2006-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/utilities/timeUtilities.hh>
#include <brick/test/testFixture.hh>

using brick::test::TestFixture;

namespace brick {

  namespace utilities {

    class TimeUtilitiesTest
      : public TestFixture<TimeUtilitiesTest>
    {
    public:

      typedef TimeUtilitiesTest TestFixtureType;
    

      TimeUtilitiesTest();
      ~TimeUtilitiesTest() {}

      void setUp(const std::string&) {}
      void tearDown(const std::string&) {}

      void testGetCurrentTime();

    private:

    
    }; // class TimeUtilitiesTest


    /* ============== Member Function Definititions ============== */

    TimeUtilitiesTest::
    TimeUtilitiesTest()
      : TestFixture<TimeUtilitiesTest>("TimeUtilitiesTest")
    {
      // Register all tests.
      BRICK_TEST_REGISTER_MEMBER(testGetCurrentTime);
    }


    void
    TimeUtilitiesTest::
    testGetCurrentTime()
    {
      // Make sure time resolution is better than 1ms.  On a typical
      // LINUX box you can change this to a number as small as
      // 0.00001 (10 microseconds) and still pass the test.
      const double requiredClockResolution = 0.001;
      const double testDuration = 1.0;

      double previousTime = getCurrentTime();
      double endTime = previousTime + testDuration;
      double updateAccumulator = 0.0;
      unsigned int numberOfUpdates = 0;
      double currentTime = 0;
      while(currentTime < endTime) {
        currentTime = getCurrentTime();
        BRICK_TEST_ASSERT(currentTime >= previousTime);
        if(currentTime > previousTime) {
          double difference = currentTime - previousTime;
          updateAccumulator += difference;
          ++numberOfUpdates;
          previousTime = currentTime;
        }
      }
      BRICK_TEST_ASSERT(numberOfUpdates != 0);
      double averageUpdateTime = updateAccumulator / numberOfUpdates;
      BRICK_TEST_ASSERT(averageUpdateTime < requiredClockResolution);
    }

  } // namespace utilities
  
} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::utilities::TimeUtilitiesTest currentTest;
  return currentTest.run();
}

#else

namespace {

  brick::utilities::TimeUtilitiesTest currentTest;
  
}

#endif
