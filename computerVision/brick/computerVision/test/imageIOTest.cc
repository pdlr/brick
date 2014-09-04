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
#include <brick/computerVision/utilities.hh>
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
      void testWritePNG_GRAY8();
      void testWritePNG_RGB8();
      void testWritePNG_GRAY16();
      void testWritePNG_RGB16();

    private:

    }; // class ImageIOTest


    /* ============== Member Function Definititions ============== */

    ImageIOTest::
    ImageIOTest()
      : brick::test::TestFixture<ImageIOTest>("ImageIOTest")
    {
      BRICK_TEST_REGISTER_MEMBER(testWritePNG_GRAY8);
      BRICK_TEST_REGISTER_MEMBER(testWritePNG_RGB8);
      BRICK_TEST_REGISTER_MEMBER(testWritePNG_GRAY16);
      BRICK_TEST_REGISTER_MEMBER(testWritePNG_RGB16);
    }


    void
    ImageIOTest::
    testWritePNG_GRAY8()
    {
      Image<GRAY8> referenceImage = readPGM8(getTestImageFileNamePGM0());
      // TBD(xxx): get a real temp file name.
      std::string outputFileName = "/var/tmp/testImage.png";
      writePNG(outputFileName, referenceImage);
      std::string commentString;
      Image<GRAY8> resultImage = readPNG<GRAY8>(
        outputFileName, commentString);
      
      BRICK_TEST_ASSERT(std::equal(resultImage.begin(), resultImage.end(),
                                   referenceImage.begin()));
    }

 
    void
    ImageIOTest::
    testWritePNG_RGB8()
    {
      Image<RGB8> referenceImage = readPPM8(getTestImageFileNamePPM0());
      // TBD(xxx): get a real temp file name.
      std::string outputFileName = "/var/tmp/testImage.png";
      writePNG(outputFileName, referenceImage);
      std::string commentString;
      Image<RGB8> resultImage = readPNG<RGB8>(
        outputFileName, commentString);
      
      BRICK_TEST_ASSERT(std::equal(resultImage.begin(), resultImage.end(),
                                   referenceImage.begin()));
    }


    void
    ImageIOTest::
    testWritePNG_GRAY16()
    {
      Image<GRAY8> referenceImage = readPGM8(getTestImageFileNamePGM0());
      Image<GRAY16> referenceImage16 =
        convertColorspace<GRAY16>(referenceImage);
      
      // TBD(xxx): get a real temp file name.
      std::string outputFileName = "/var/tmp/testImage.png";
      writePNG(outputFileName, referenceImage16);
      std::string commentString;
      Image<GRAY16> resultImage = readPNG<GRAY16>(
        outputFileName, commentString);
      
      BRICK_TEST_ASSERT(std::equal(resultImage.begin(), resultImage.end(),
                                   referenceImage16.begin()));
    }

 
    void
    ImageIOTest::
    testWritePNG_RGB16()
    {
      Image<RGB8> referenceImage = readPPM8(getTestImageFileNamePPM0());
      Image<RGB16> referenceImage16 =
        convertColorspace<RGB16>(referenceImage);
      
      // TBD(xxx): get a real temp file name.
      std::string outputFileName = "/var/tmp/testImage.png";
      writePNG(outputFileName, referenceImage16);
      std::string commentString;
      Image<RGB16> resultImage = readPNG<RGB16>(
        outputFileName, commentString);
      
      BRICK_TEST_ASSERT(std::equal(resultImage.begin(), resultImage.end(),
                                   referenceImage16.begin()));
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
