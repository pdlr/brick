/**
***************************************************************************
* @file imagePyramidBinomialTest.cpp
*
* Source file defining tests for the ImagePyramidBinomial class template.
*
* Copyright (C) 2011,2015 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/computerVision/imageIO.hh>
#include <brick/computerVision/imagePyramidBinomial.hh>
#include <brick/computerVision/test/testImages.hh>
#include <brick/test/testFixture.hh>

#include <brick/utilities/timeUtilities.hh>

namespace brick {

  namespace computerVision {

    
    // This test is copied 
    class ImagePyramidBinomialTest
      : public brick::test::TestFixture<ImagePyramidBinomialTest> {

    public:

      ImagePyramidBinomialTest();
      ~ImagePyramidBinomialTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testImagePyramidBinomial();

    private:

      double m_defaultTolerance;
      
    }; // class ImagePyramidBinomialTest


    /* ============== Member Function Definititions ============== */

    ImagePyramidBinomialTest::
    ImagePyramidBinomialTest()
      : brick::test::TestFixture<ImagePyramidBinomialTest>("ImagePyramidBinomialTest"),
        m_defaultTolerance(1.0E-8)
    {
      BRICK_TEST_REGISTER_MEMBER(testImagePyramidBinomial);
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
    ImagePyramidBinomialTest::
    testImagePyramidBinomial()
    {
      Image<GRAY8> inputImageGray = readPGM8(getTestImageFileNamePGM0());
      // Image<GRAY8> inputImageGray2 = readPGM8(getTestImageFileNamePGM0());
      // inputImageGray -= inputImageGray2;
      Image<GRAY_FLOAT32> floatImageGray = convertColorspace<GRAY_FLOAT32>(
        inputImageGray);

      double time0 = utilities::getCurrentTime();
      ImagePyramidBinomial<GRAY_FLOAT32, GRAY_FLOAT32> pyramidGray(
        floatImageGray, 0, 6, false);
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
      ImagePyramidBinomial<RGB_FLOAT32, RGB_FLOAT32> pyramidRGB(floatImageRGB);
      time1 = utilities::getCurrentTime();
      std::cout << "ET: " << time1 - time0 << " to build a pyramid for a "
                << inputImageGray.rows() << "x" << inputImageGray.columns()
                << " RGB image." << std::endl;
      
    }

  } // namespace computerVision

} // namespace brick


#if 1

int main(/* int argc, char** argv */)
{
  brick::computerVision::ImagePyramidBinomialTest currentTest;
  currentTest.testImagePyramidBinomial();
  
  // bool result = currentTest.run();
  bool result = true;
  return (result ? 0 : 1);
}

#else

namespace {

  brick::computerVision::ImagePyramidBinomialTest currentTest;

}

#endif
