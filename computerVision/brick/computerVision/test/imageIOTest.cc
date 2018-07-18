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

#include <stdint.h>

#include <fstream>
#include <iostream>

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
      void testReadPGM16();

#if HAVE_LIBPNG
      void testWritePNG_GRAY8();
      void testWritePNG_RGB8();
      void testWritePNG_GRAY16();
      void testWritePNG_RGB16();
#endif
    private:

    }; // class ImageIOTest


    /* ============== Member Function Definititions ============== */

    ImageIOTest::
    ImageIOTest()
      : brick::test::TestFixture<ImageIOTest>("ImageIOTest")
    {
      BRICK_TEST_REGISTER_MEMBER(testReadPGM16);
#if HAVE_LIBPNG
      BRICK_TEST_REGISTER_MEMBER(testWritePNG_GRAY8);
      BRICK_TEST_REGISTER_MEMBER(testWritePNG_RGB8);
      BRICK_TEST_REGISTER_MEMBER(testWritePNG_GRAY16);
      BRICK_TEST_REGISTER_MEMBER(testWritePNG_RGB16);
#endif
    }


    void
    ImageIOTest::
    testReadPGM16()
    {
      std::string testImageFileName = "/var/tmp/brickTestImage.pgm";
      uint32_t const imageWidth = 3;
      uint32_t const imageHeight = 2;
      uint8_t imageData[] = {0x00, 0x00, 0x00, 0xe8, 0x15, 0x11,
                             0xfd, 0xf0, 0x22, 0x22, 0x21, 0x01};

      std::ofstream outputStream(testImageFileName.c_str());
      if(!outputStream) {
        BRICK_THROW(brick::common::IOException, "ImageIOTest::testReadPGM16()",
                    "Unable to open output stream to "
                    "/var/tmp/brickTestImage.pgm");
      }

      outputStream << "P5\n"
                   << imageWidth << " " << imageHeight << "\n"
                   << "65535\n";
      outputStream.write(reinterpret_cast<const char*>(imageData),
                         imageHeight * imageWidth * 2);
      outputStream.close();
      Image<GRAY16> inputImage = readPGM16(testImageFileName);

      BRICK_TEST_ASSERT(inputImage.rows() == imageHeight);
      BRICK_TEST_ASSERT(inputImage.columns() == imageWidth);

      uint32_t ii = 0;
      for(uint32_t rr = 0; rr < inputImage.rows(); ++rr) {
        for(uint32_t cc = 0; cc < inputImage.columns(); ++cc) {
          uint16_t pixelValue = inputImage(rr, cc);
          BRICK_TEST_ASSERT(static_cast<uint8_t>((pixelValue & 0xff00) >> 8)
                            == imageData[ii]);
          BRICK_TEST_ASSERT(static_cast<uint8_t>(pixelValue & 0x00ff)
                            == imageData[ii + 1]);
          ii += 2;
        }
      }
    }

#if HAVE_LIBPNG
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
#endif
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
