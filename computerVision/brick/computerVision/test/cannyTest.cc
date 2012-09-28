/**
***************************************************************************
* @file cannyTest.cc
*
* Source file defining tests for routines defined in
* brick/computerVision/canny.hh.
*
* Copyright (C) 2006,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/computerVision/test/testImages.hh>
#include <brick/computerVision/canny.hh>
#include <brick/computerVision/imageIO.hh>
#include <brick/computerVision/utilities.hh>
#include <brick/test/testFixture.hh>


namespace brick {

  namespace computerVision {
    
    class CannyTest
      : public brick::test::TestFixture<CannyTest> {

    public:

      CannyTest();
      ~CannyTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testCanny();

    private:

    }; // class CannyTest


    /* ============== Member Function Definititions ============== */

    CannyTest::
    CannyTest()
      : brick::test::TestFixture<CannyTest>("CannyTest")
    {
      BRICK_TEST_REGISTER_MEMBER(testCanny);
    }


    void
    CannyTest::
    testCanny()
    {
      Image<GRAY8> inputImage0 = readPGM8(getTestImageFileNamePGM0());
      Image<GRAY8> referenceImage = readPGM8(getEdgeImageFileNamePGM0());
      Image<GRAY1> binaryImage = applyCanny<double>(inputImage0, 5, 5.0, 1.0);
      Image<GRAY8> edgeImage = convertColorspace<GRAY8>(binaryImage);

      BRICK_TEST_ASSERT(edgeImage.rows() == referenceImage.rows());
      BRICK_TEST_ASSERT(edgeImage.columns() == referenceImage.columns());
      for(size_t index0 = 0; index0 < edgeImage.size(); ++index0) {
        BRICK_TEST_ASSERT(edgeImage[index0] == referenceImage[index0]);
      }
    }

  } // namespace computerVision

} // namespace brick

#if 0

int main(int argc, char** argv)
{
  brick::computerVision::CannyTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::computerVision::CannyTest currentTest;
  
}

#endif
