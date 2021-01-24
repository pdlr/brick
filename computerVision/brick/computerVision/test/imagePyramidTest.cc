/**
***************************************************************************
* @file imagePyramidTest.cpp
*
* Source file defining tests for the ImagePyramid class template.
*
* Copyright (C) 2011 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/computerVision/imageIO.hh>
#include <brick/computerVision/imagePyramid.hh>
#include <brick/computerVision/test/testImages.hh>
#include <brick/test/testFixture.hh>

#include <brick/utilities/timeUtilities.hh>

namespace brick {

  namespace computerVision {


    // This test is copied
    class ImagePyramidTest
      : public brick::test::TestFixture<ImagePyramidTest> {

    public:

      ImagePyramidTest();
      ~ImagePyramidTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testImagePyramid();

    private:

      double m_defaultTolerance;

    }; // class ImagePyramidTest


    /* ============== Member Function Definititions ============== */

    ImagePyramidTest::
    ImagePyramidTest()
      : brick::test::TestFixture<ImagePyramidTest>("ImagePyramidTest"),
        m_defaultTolerance(1.0E-8)
    {
      BRICK_TEST_REGISTER_MEMBER(testImagePyramid);
    }


    template<class Type0, class Type1>
    inline PixelRGB<Type0>
    subtractMe(const PixelRGB<Type0>& pixel0, const PixelRGB<Type1>& pixel1)
    {
      if(pixel1.red) {
        return PixelRGB<Type0>(pixel0);
      }
      return PixelRGB<Type0>(pixel1);
    }


    void
    ImagePyramidTest::
    testImagePyramid()
    {
      Image<GRAY8> inputImageGray = readPGM8(getTestImageFileNamePGM0());
      // Image<GRAY8> inputImageGray2 = readPGM8(getTestImageFileNamePGM0());
      // inputImageGray -= inputImageGray2;
      Image<GRAY_FLOAT32> floatImageGray = convertColorspace<GRAY_FLOAT32>(
        inputImageGray);

      double time0 = utilities::getCurrentTime();
      ImagePyramid<GRAY_FLOAT32, GRAY_FLOAT32, double> pyramidGray(
        floatImageGray);
      double time1 = utilities::getCurrentTime();
      std::cout << "ET: " << time1 - time0 << " to build a pyramid for a "
                << inputImageGray.rows() << "x" << inputImageGray.columns()
                << " gray image." << std::endl;

#if 0
      for(unsigned int ii = 0; ii < pyramidGray.getNumberOfLevels(); ++ii) {
        Image<GRAY_FLOAT32> currentImage = pyramidGray.getLevel(ii);
        std::ostringstream fileNameStream;
        fileNameStream << "outImage" << ii << ".pgm";
        writePGM(fileNameStream.str(), currentImage.data(),
                 currentImage.rows(), currentImage.columns(),
                 true);
      }
#endif

      Image<RGB8> inputImageRGB = readPPM8(getTestImageFileNamePPM0());
      // Image<RGB8> inputImageRGB2 = readPPM8(getTestImageFileNamePPM0());
      // PixelRGB8 pixel;
      // PixelRGB8 pixel2;
      // PixelRGB8 pixel3 = subtractMe<brick::common::UnsignedInt8,
      //    brick::common::UnsignedInt8>(pixel, pixel2);
      // inputImageRGB -= inputImageRGB2;
      Image<RGB_FLOAT32> floatImageRGB = convertColorspace<RGB_FLOAT32>(
        inputImageRGB);

      time0 = utilities::getCurrentTime();
      ImagePyramid<RGB_FLOAT32, RGB_FLOAT32, double> pyramidRGB(floatImageRGB);
      time1 = utilities::getCurrentTime();
      std::cout << "ET: " << time1 - time0 << " to build a pyramid for a "
                << inputImageGray.rows() << "x" << inputImageGray.columns()
                << " RGB image." << std::endl;

    }

  } // namespace computerVision

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::computerVision::ImagePyramidTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::computerVision::ImagePyramidTest currentTest;

}

#endif
