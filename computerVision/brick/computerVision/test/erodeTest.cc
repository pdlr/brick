/**
***************************************************************************
* @file erodeTest.cpp
*
* Source file defining tests for erode().
*
* Copyright (C) 2006 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/computerVision/test/testImages.hh>
#include <brick/computerVision/erode.hh>
#include <brick/computerVision/imageIO.hh>
#include <brick/test/testFixture.hh>
#include <brick/utilities/timeUtilities.hh>

namespace brick {

  namespace computerVision {
    
    class ErodeTest
      : public brick::test::TestFixture<ErodeTest> {

    public:

      ErodeTest();
      ~ErodeTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testErode();
      void testErodeUsingBoxIntegrator();

    private:

    }; // class ErodeTest


    /* ============== Member Function Definititions ============== */

    ErodeTest::
    ErodeTest()
      : TestFixture<ErodeTest>("ErodeTest")
    {
      BRICK_TEST_REGISTER_MEMBER(testErode);
      BRICK_TEST_REGISTER_MEMBER(testErodeUsingBoxIntegrator);
    }


    void
    ErodeTest::
    testErode()
    {
      Image<GRAY8> inputImage = readPGM8(getDilateErodeFileNamePGM0());
      Image<GRAY8> referenceImage = readPGM8(getErodedFileNamePGM0());
      Image<GRAY8> erodedImage = erode<GRAY8>(inputImage);
      BRICK_TEST_ASSERT(erodedImage.rows() == referenceImage.rows());
      BRICK_TEST_ASSERT(erodedImage.columns() == referenceImage.columns());
      for(size_t index0 = 0; index0 < inputImage.size(); ++index0) {
        BRICK_TEST_ASSERT(erodedImage[index0] == referenceImage[index0]);
      }
    }


    void
    ErodeTest::
    testErodeUsingBoxIntegrator()
    {
      Image<GRAY8> inputImage = readPGM8(getDilateErodeFileNamePGM0());
      Image<GRAY8> referenceImage = readPGM8(getErodedFileNamePGM0());
      Image<GRAY8> erodedImage = erodeUsingBoxIntegrator<GRAY8>(inputImage);
      BRICK_TEST_ASSERT(erodedImage.rows() == referenceImage.rows());
      BRICK_TEST_ASSERT(erodedImage.columns() == referenceImage.columns());
      for(size_t index0 = 0; index0 < inputImage.size(); ++index0) {
        BRICK_TEST_ASSERT(erodedImage[index0] == referenceImage[index0]);
      }

#if 0
      // Now evaluate execution speed.
      const unsigned int numIterations = 100;

      brick::utilities::Timer timer0;
      for(unsigned int iteration = 0; iteration < numIterations; ++iteration) {
        erodedImage = erode<GRAY8>(inputImage);
      }
      double elapsedTimeThreshold = timer0.reset();
      for(unsigned int iteration = 0; iteration < numIterations; ++iteration) {
        erodedImage = erodeUsingBoxIntegrator<GRAY8>(inputImage);
      }
      double elapsedTime = timer0.reset();
      try {
        BRICK_TEST_ASSERT(elapsedTime < elapsedTimeThreshold);
      } catch(...) {
        std::cout << "elapsedTime = " << elapsedTime
                  << ", elapsedTimeThreshold = " << elapsedTimeThreshold
                  << "." << std::endl;
        throw;
      }
#endif
    }

  } // namespace computerVision

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::computerVision::ErodeTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::computerVision::ErodeTest currentTest;

}

#endif

