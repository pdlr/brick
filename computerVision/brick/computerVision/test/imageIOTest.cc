/**
***************************************************************************
* @file imageIOTest.cpp
*
* Source file defining tests for image I/O routines.
*
* Copyright (C) 2014 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/computerVision/test/testImages.hh>
#include <brick/computerVision/image.hh>
#include <brick/computerVision/imageIO.hh>
#include <brick/test/testFixture.hh>

namespace brick {

  namespace computerVision {
    
    class ImageIOTest : public brick::test::TestFixture<ImageIOTest> {

    public:

      ImageIOTest();
      ~ImageIOTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests of member functions.
      void testWritePNG8_GRAY8();
      void testWritePNG8_RGB8();

    private:

    }; // class ImageIOTest


    /* ============== Member Function Definititions ============== */

    ImageIOTest::
    ImageIOTest()
      : brick::test::TestFixture<ImageIOTest>("ImageIOTest")
    {
      BRICK_TEST_REGISTER_MEMBER(testWritePNG8_GRAY8);
      BRICK_TEST_REGISTER_MEMBER(testWritePNG8_RGB8);
    }


    void
    ImageIOTest::
    testWritePNG8_GRAY8()
    {
      Image<GRAY8> referenceImage = readPGM8(getTestImageFileNamePGM0());
      // TBD(xxx): get a real temp file name.
      std::string outputFileName = "/var/tmp/testImage.png";
      writePNG8(outputFileName, referenceImage);
      std::string commentString;
      Image<GRAY8> resultImage = readPNG8<GRAY8>(
        outputFileName, commentString);
      
      BRICK_TEST_ASSERT(std::equal(resultImage.begin(), resultImage.end(),
                                   referenceImage.begin()));
    }

 
    void
    ImageIOTest::
    testWritePNG8_RGB8()
    {
      Image<RGB8> referenceImage = readPPM8(getTestImageFileNamePPM0());
      // TBD(xxx): get a real temp file name.
      std::string outputFileName = "/var/tmp/testImage.png";
      writePNG8(outputFileName, referenceImage);
      std::string commentString;
      Image<RGB8> resultImage = readPNG8<RGB8>(
        outputFileName, commentString);
      
      BRICK_TEST_ASSERT(std::equal(resultImage.begin(), resultImage.end(),
                                   referenceImage.begin()));
    }

  } // namespace computerVision

} // namespace brick

#if 0

int main(int argc, char** argv)
{
  brick::computerVision::ImageIOTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::computerVision::ImageIOTest currentTest;

}

#endif
