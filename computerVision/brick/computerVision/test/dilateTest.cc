/**
***************************************************************************
* @file dilateTest.cpp
*
* Source file defining tests for dilate().
*
* Copyright (C) 2006 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/computerVision/test/testImages.hh>
#include <brick/computerVision/dilate.hh>
#include <brick/computerVision/erode.hh>
#include <brick/computerVision/imageIO.hh>
#include <brick/test/testFixture.hh>


namespace brick {

  namespace computerVision {
    
    class DilateTest
      : public brick::test::TestFixture<DilateTest> {

    public:

      DilateTest();
      ~DilateTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testDilate();
      void testDilate__radius();

    private:

    }; // class DilateTest


    /* ============== Member Function Definititions ============== */

    DilateTest::
    DilateTest()
      : TestFixture<DilateTest>("DilateTest")
    {
      BRICK_TEST_REGISTER_MEMBER(testDilate);
      BRICK_TEST_REGISTER_MEMBER(testDilate__radius);
    }


    void
    DilateTest::
    testDilate()
    {
      Image<GRAY8> inputImage = readPGM8(getDilateErodeFileNamePGM0());
      Image<GRAY8> referenceImage = readPGM8(getDilatedFileNamePGM0());
      Image<GRAY8> dilatedImage = dilate<GRAY8>(inputImage) * brick::common::UnsignedInt8(255);
      BRICK_TEST_ASSERT(dilatedImage.rows() == referenceImage.rows());
      BRICK_TEST_ASSERT(dilatedImage.columns() == referenceImage.columns());
      for(size_t index0 = 0; index0 < inputImage.size(); ++index0) {
        BRICK_TEST_ASSERT(dilatedImage[index0] == referenceImage[index0]);
      }
    }

    void
    DilateTest::
    testDilate__radius()
    {
      Image<GRAY8> inputImage = readPGM8(getDilateErodeFileNamePGM0());

      // Make sure 3x3 dilation matches the hardcoded routine tested
      // above.
      Image<GRAY8> referenceImage = dilate<GRAY8>(inputImage);
      Image<GRAY8> dilatedImage = dilateUsingBoxIntegrator<GRAY8>(inputImage, 3, 3);
      BRICK_TEST_ASSERT(dilatedImage.rows() == referenceImage.rows());
      BRICK_TEST_ASSERT(dilatedImage.columns() == referenceImage.columns());
      for(size_t index0 = 0; index0 < dilatedImage.size(); ++index0) {
        BRICK_TEST_ASSERT(dilatedImage[index0] == referenceImage[index0]);
      }

      // Use the hardcoded routine again to test 5x5 dilation.
      Image<GRAY8> referenceImage2 = dilate<GRAY8>(referenceImage);
      Image<GRAY8> dilatedImage2 = dilateUsingBoxIntegrator<GRAY8>(inputImage, 5, 5);
      BRICK_TEST_ASSERT(dilatedImage2.rows() == referenceImage2.rows());
      BRICK_TEST_ASSERT(dilatedImage2.columns() == referenceImage2.columns());
      for(size_t index0 = 0; index0 < dilatedImage2.size(); ++index0) {
        BRICK_TEST_ASSERT(dilatedImage2[index0] == referenceImage2[index0]);
      }
      
      
      // Actually, the 5x5 dilation is mostly white.  Test again using a
      // sparser image.
      Image<GRAY8> inputImage3 = erode<GRAY8>(inputImage);
      Image<GRAY8> referenceImage3 = dilate<GRAY8>(dilate<GRAY8>(inputImage3));
      Image<GRAY8> dilatedImage3 = dilateUsingBoxIntegrator<GRAY8>(inputImage3, 5, 5);
      BRICK_TEST_ASSERT(dilatedImage3.rows() == referenceImage3.rows());
      BRICK_TEST_ASSERT(dilatedImage3.columns() == referenceImage3.columns());
      for(size_t index0 = 0; index0 < dilatedImage3.size(); ++index0) {
        BRICK_TEST_ASSERT(dilatedImage3[index0] == referenceImage3[index0]);
      }
      
    }

  } // namespace computerVision

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::computerVision::DilateTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::computerVision::DilateTest currentTest;

}

#endif
