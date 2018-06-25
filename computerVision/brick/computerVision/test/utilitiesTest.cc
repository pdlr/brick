/**
***************************************************************************
* @file utilitiesTest.cpp
*
* Source file defining tests for dlrComputerVision library utilities.
*
* Copyright (C) 2006 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/computerVision/test/testImages.hh>
#include <brick/computerVision/colorspaceConverter.hh>
#include <brick/computerVision/image.hh>
#include <brick/computerVision/imageIO.hh>
#include <brick/computerVision/utilities.hh>
#include <brick/numeric/transform2D.hh>
#include <brick/numeric/vector2D.hh>
#include <brick/test/testFixture.hh>


namespace brick {

  namespace computerVision {

    class UtilitiesTest : public brick::test::TestFixture<UtilitiesTest> {

    public:

      UtilitiesTest();
      ~UtilitiesTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testAssociateColorComponents();
      void testConvertColorspace();
      void testDissociateColorComponents();
      void testEstimateAffineTransform();
      void testSubsample();
      void testSupersample();
      void testToArray();

    private:

    }; // class UtilitiesTest


    /* ============== Member Function Definititions ============== */

    UtilitiesTest::
    UtilitiesTest()
      : brick::test::TestFixture<UtilitiesTest>("UtilitiesTest")
    {
      BRICK_TEST_REGISTER_MEMBER(testAssociateColorComponents);
      BRICK_TEST_REGISTER_MEMBER(testConvertColorspace);
      BRICK_TEST_REGISTER_MEMBER(testDissociateColorComponents);
      BRICK_TEST_REGISTER_MEMBER(testEstimateAffineTransform);
      BRICK_TEST_REGISTER_MEMBER(testSubsample);
      BRICK_TEST_REGISTER_MEMBER(testSupersample);
      BRICK_TEST_REGISTER_MEMBER(testToArray);
    }


    void
    UtilitiesTest::
    testAssociateColorComponents()
    {
      Image<RGB8> inputImage = readPPM8(getTestImageFileNamePPM0());
      brick::numeric::Array2D<brick::common::UnsignedInt8> dataArray =
        dissociateColorComponents(inputImage);
      Image<RGB8> imageAlias = associateColorComponents<RGB8>(dataArray);
      BRICK_TEST_ASSERT(imageAlias.rows() == inputImage.rows());
      BRICK_TEST_ASSERT(imageAlias.columns() == inputImage.columns());
      BRICK_TEST_ASSERT(imageAlias.data() == inputImage.data());

      Image<HSV_FLOAT64> hsvImage = convertColorspace<HSV_FLOAT64>(inputImage);
      brick::numeric::Array2D<double> dataArray2 =
        dissociateColorComponents(hsvImage);
      Image<HSV_FLOAT64> imageAlias2 =
        associateColorComponents<HSV_FLOAT64>(dataArray2);
      BRICK_TEST_ASSERT(imageAlias2.rows() == hsvImage.rows());
      BRICK_TEST_ASSERT(imageAlias2.columns() == hsvImage.columns());
      BRICK_TEST_ASSERT(imageAlias2.data() == hsvImage.data());
    }


    void
    UtilitiesTest::
    testConvertColorspace()
    {
      Image<RGB8> inputImage2 = readPPM8(getTestImageFileNamePPM0());
      {
        Image<GRAY8> grayImage = convertColorspace<GRAY8>(inputImage2);
        BRICK_TEST_ASSERT(grayImage.rows() == inputImage2.rows());
        BRICK_TEST_ASSERT(grayImage.columns() == inputImage2.columns());
        ColorspaceConverter<RGB8, GRAY8> converter;
        for(size_t pixelIndex = 0; pixelIndex < inputImage2.size();
            ++pixelIndex) {
          BRICK_TEST_ASSERT(
            grayImage[pixelIndex] == converter(inputImage2[pixelIndex]));
        }
      }

      {
        Image<RGBA8> rgbaImage = convertColorspace<RGBA8>(inputImage2);
        BRICK_TEST_ASSERT(rgbaImage.rows() == inputImage2.rows());
        BRICK_TEST_ASSERT(rgbaImage.columns() == inputImage2.columns());
        ColorspaceConverter<RGB8, RGBA8> converter;
        for(size_t pixelIndex = 0; pixelIndex < inputImage2.size();
            ++pixelIndex) {
          BRICK_TEST_ASSERT(
            rgbaImage[pixelIndex] == converter(inputImage2[pixelIndex]));
        }
      }
    }


    void
    UtilitiesTest::
    testDissociateColorComponents()
    {
      Image<RGB8> inputImage = readPPM8(getTestImageFileNamePPM0());
      brick::numeric::Array2D<brick::common::UnsignedInt8> dataArray =
        dissociateColorComponents(inputImage);
      BRICK_TEST_ASSERT(dataArray.rows() == inputImage.rows());
      BRICK_TEST_ASSERT(dataArray.columns() == 3 * inputImage.columns());
      size_t dataArrayIndex = 0;
      for(size_t pixelIndex = 0; pixelIndex < inputImage.size(); ++pixelIndex) {
        BRICK_TEST_ASSERT(
          inputImage[pixelIndex].red == dataArray[dataArrayIndex++]);
        BRICK_TEST_ASSERT(
          inputImage[pixelIndex].green == dataArray[dataArrayIndex++]);
        BRICK_TEST_ASSERT(
          inputImage[pixelIndex].blue == dataArray[dataArrayIndex++]);
      }

      Image<HSV_FLOAT64> hsvImage = convertColorspace<HSV_FLOAT64>(inputImage);
      brick::numeric::Array2D<double> dataArray2 =
        dissociateColorComponents(hsvImage);
      BRICK_TEST_ASSERT(dataArray2.rows() == hsvImage.rows());
      BRICK_TEST_ASSERT(dataArray2.columns() == 3 * hsvImage.columns());
      size_t dataArray2Index = 0;
      for(size_t pixelIndex = 0; pixelIndex < hsvImage.size(); ++pixelIndex) {
        BRICK_TEST_ASSERT(
          hsvImage[pixelIndex].hue == dataArray2[dataArray2Index++]);
        BRICK_TEST_ASSERT(
          hsvImage[pixelIndex].saturation == dataArray2[dataArray2Index++]);
        BRICK_TEST_ASSERT(
          hsvImage[pixelIndex].value == dataArray2[dataArray2Index++]);
      }

    }


    void
    UtilitiesTest::
    testEstimateAffineTransform()
    {
      const unsigned int numberOfPoints = 4;
      const common::Float64 defaultTolerance = 1.0E-8;

      std::vector< numeric::Vector2D<common::Float64> > fromPoints(
        numberOfPoints);
      std::vector< numeric::Vector2D<common::Float64> > toPoints(
        numberOfPoints);
      numeric::Transform2D<common::Float64> nominalTransform(2.0, -1.0, 3.0,
                                                             -4.0, 5.0, 6.0,
                                                             0.0, 0.0, 1.0);
      fromPoints[0] = numeric::Vector2D<common::Float64>(1.0, 4.0);
      fromPoints[1] = numeric::Vector2D<common::Float64>(-1.0, 2.0);
      fromPoints[2] = numeric::Vector2D<common::Float64>(-1.0, -3.0);
      fromPoints[3] = numeric::Vector2D<common::Float64>(2.0, 1.0);
      for(unsigned int ii = 0; ii < numberOfPoints; ++ii) {
        toPoints[ii] = nominalTransform * fromPoints[ii];
      }

      // Solve for thetransform, and make sure we came close to the
      // right answer.
      numeric::Transform2D<common::Float64> result =
        estimateAffineTransform<common::Float64>(
          toPoints.begin(), toPoints.end(), fromPoints.begin());
      for(unsigned int row = 0; row < 3; ++row) {
        for(unsigned int column = 0; column < 3; ++column) {
          BRICK_TEST_ASSERT(
            approximatelyEqual(
              result(row, column), nominalTransform(row, column),
              defaultTolerance));
        }
      }
    }


    void
    UtilitiesTest::
    testSubsample()
    {
      const size_t imageRows = 2048;
      const size_t imageColumns = 1962;

      Image<RGB8> inputImage(imageRows, imageColumns);
      for(size_t rowIndex = 0; rowIndex < imageRows; ++rowIndex) {
        for(size_t columnIndex = 0; columnIndex < imageColumns; ++columnIndex) {
          PixelRGB8 inputPixel(static_cast<brick::common::UnsignedInt8>(rowIndex),
                               static_cast<brick::common::UnsignedInt8>(columnIndex),
                               static_cast<brick::common::UnsignedInt8>(0));
          inputImage(rowIndex, columnIndex) = inputPixel;
        }
      }


      for(size_t rowStep = 2; rowStep < 4; ++rowStep) {
        for(size_t columnStep = 2; columnStep < 4; ++columnStep) {
          Image<RGB8> outputImage0 = subsample(inputImage, rowStep, columnStep);

          BRICK_TEST_ASSERT(
            (outputImage0.rows() - 1) * rowStep + 1
            <= inputImage.rows());
          BRICK_TEST_ASSERT(
            (outputImage0.columns() - 1) * columnStep + 1
            <= inputImage.columns());
          BRICK_TEST_ASSERT(
            (outputImage0.rows() - 1) * rowStep + 1
            > inputImage.rows() - rowStep);
          BRICK_TEST_ASSERT(
            (outputImage0.columns() - 1) * columnStep + 1
            > inputImage.columns() - columnStep);

          for(size_t rowIndex = 0; rowIndex < outputImage0.rows();
              ++rowIndex) {
            for(size_t columnIndex = 0; columnIndex < outputImage0.columns();
                ++columnIndex) {
              BRICK_TEST_ASSERT(
                outputImage0(rowIndex, columnIndex)
                == inputImage(rowIndex * rowStep, columnIndex * columnStep));
            }
          }
        }
      }
    }


    void
    UtilitiesTest::
    testSupersample()
    {
      const size_t imageRows = 100;
      const size_t imageColumns = 121;

      Image<GRAY8> inputImage0(imageRows, imageColumns);
      Image<RGB8> inputImage1(imageRows, imageColumns);
      for(size_t rowIndex = 0; rowIndex < imageRows; ++rowIndex) {
        for(size_t columnIndex = 0; columnIndex < imageColumns; ++columnIndex) {
          inputImage0(rowIndex, columnIndex) =
            static_cast<brick::common::UnsignedInt8>(rowIndex + columnIndex);

          PixelRGB8 inputPixel1(static_cast<brick::common::UnsignedInt8>(rowIndex),
                                static_cast<brick::common::UnsignedInt8>(columnIndex),
                                static_cast<brick::common::UnsignedInt8>(0));
          inputImage1(rowIndex, columnIndex) = inputPixel1;
        }
      }


      Image<GRAY8> outputImage0 =
        supersample<GRAY8, GRAY8, GRAY16>(inputImage0);
      Image<RGB8> outputImage1 =
        supersample<RGB8, RGB8, RGB16>(inputImage1);

      BRICK_TEST_ASSERT(outputImage0.rows() == inputImage0.rows() * 2 - 1);
      BRICK_TEST_ASSERT(outputImage0.columns() == inputImage0.columns() * 2 - 1);
      for(size_t rowIndex = 0; rowIndex < outputImage0.rows();
          ++rowIndex) {
        for(size_t columnIndex = 0; columnIndex < outputImage0.columns();
            ++columnIndex) {
          brick::common::UnsignedInt8 outputValue =
            static_cast<brick::common::UnsignedInt8>((rowIndex + columnIndex) / 2);
          BRICK_TEST_ASSERT(outputImage0(rowIndex, columnIndex) == outputValue);
        }
      }

      BRICK_TEST_ASSERT(outputImage1.rows() == inputImage1.rows() * 2 - 1);
      BRICK_TEST_ASSERT(outputImage1.columns() == inputImage1.columns() * 2 - 1);
      for(size_t rowIndex = 0; rowIndex < outputImage0.rows();
          ++rowIndex) {
        for(size_t columnIndex = 0; columnIndex < outputImage0.columns();
            ++columnIndex) {
          PixelRGB8 outputPixel = outputImage1(rowIndex, columnIndex);
          BRICK_TEST_ASSERT(outputPixel.red == rowIndex / 2);
          BRICK_TEST_ASSERT(outputPixel.green == columnIndex / 2);
          BRICK_TEST_ASSERT(outputPixel.blue == 0);
        }
      }
    }


    void
    UtilitiesTest::
    testToArray()
    {
      Image<RGB8> inputImage = readPPM8(getTestImageFileNamePPM0());
      brick::numeric::Array2D<double> dataArray = toArray<double>(inputImage);

      BRICK_TEST_ASSERT(dataArray.rows() == inputImage.rows());
      BRICK_TEST_ASSERT(dataArray.columns() == 3 * inputImage.columns());
      size_t dataArrayIndex = 0;
      for(size_t pixelIndex = 0; pixelIndex < inputImage.size(); ++pixelIndex) {
        BRICK_TEST_ASSERT(
          inputImage[pixelIndex].red == dataArray[dataArrayIndex++]);
        BRICK_TEST_ASSERT(
          inputImage[pixelIndex].green == dataArray[dataArrayIndex++]);
        BRICK_TEST_ASSERT(
          inputImage[pixelIndex].blue == dataArray[dataArrayIndex++]);
      }
    }

  } // namespace computerVision

} // namespace brick


#if 1

int main(int /* argc */, char** /* argv */)
{
  brick::computerVision::UtilitiesTest currentTest;
  currentTest.testSupersample();
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::computerVision::UtilitiesTest currentTest;

}

#endif
