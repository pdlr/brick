/**
***************************************************************************
* @file filterTest.cpp
*
* Source file defining tests for ColorspaceConverter classes.
*
* Copyright (C) 2006 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/computerVision/image.hh>
#include <brick/computerVision/colorspaceConverter.hh>
#include <brick/computerVision/pixelRGB.hh>
#include <brick/test/testFixture.hh>

namespace brick {

  namespace computerVision {

    class ColorspaceConverterTest
      : public brick::test::TestFixture<ColorspaceConverterTest> {

    public:

      ColorspaceConverterTest();
      ~ColorspaceConverterTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testRGB8ToGRAY8();
      void testRGB8ToBGRA8();
      void testRGB8ToRGBA8();
      void testRGB8ToHSV_FLOAT64();
      void testRGB8ToYIQ_FLOAT64();
      void testRGB_FLOAT64ToHSV_FLOAT64();
      void testRGB_FLOAT64ToYIQ_FLOAT64();
      void testBGRA8ToRGB8();
      void testRGBA8ToRGB8();
      void testHSV_FLOAT64ToRGB8();
    private:

    }; // class ColorspaceConverterTest


    /* ============== Member Function Definititions ============== */

    ColorspaceConverterTest::
    ColorspaceConverterTest()
      : brick::test::TestFixture<ColorspaceConverterTest>(
        "ColorspaceConverterTest")
    {
      BRICK_TEST_REGISTER_MEMBER(testRGB8ToGRAY8);
      BRICK_TEST_REGISTER_MEMBER(testRGB8ToBGRA8);
      BRICK_TEST_REGISTER_MEMBER(testRGB8ToRGBA8);
      BRICK_TEST_REGISTER_MEMBER(testRGB8ToHSV_FLOAT64);
      BRICK_TEST_REGISTER_MEMBER(testRGB8ToYIQ_FLOAT64);
      BRICK_TEST_REGISTER_MEMBER(testRGB_FLOAT64ToHSV_FLOAT64);
      BRICK_TEST_REGISTER_MEMBER(testRGB_FLOAT64ToYIQ_FLOAT64);
      BRICK_TEST_REGISTER_MEMBER(testRGBA8ToRGB8);
      BRICK_TEST_REGISTER_MEMBER(testRGBA8ToRGB8);
      BRICK_TEST_REGISTER_MEMBER(testHSV_FLOAT64ToRGB8);
    }


    void
    ColorspaceConverterTest::
    testRGB8ToGRAY8()
    {
      ColorspaceConverter<RGB8, GRAY8> converter;
      for(brick::common::UnsignedInt16 redValue = 0; redValue < 256; redValue += 7) {
        for(brick::common::UnsignedInt16 greenValue = 0; greenValue < 256; greenValue += 7) {
          for(brick::common::UnsignedInt16 blueValue = 0; blueValue < 256; blueValue += 7) {
            double redDbl = redValue;
            double greenDbl = greenValue;
            double blueDbl = blueValue;
            double grayDbl = 0.3 * redDbl + 0.59 * greenDbl + 0.11 * blueDbl;
            brick::common::UnsignedInt8 grayValue = static_cast<brick::common::UnsignedInt8>(grayDbl + 0.5);

            PixelRGB8 inputPixel(static_cast<brick::common::UnsignedInt8>(redValue),
								 static_cast<brick::common::UnsignedInt8>(greenValue),
								 static_cast<brick::common::UnsignedInt8>(blueValue));
            BRICK_TEST_ASSERT(converter(inputPixel) == grayValue);
          }
        }
      }
      // Special case tests follow.
      PixelRGB8 inputPixel(255, 255, 255);
      BRICK_TEST_ASSERT(converter(inputPixel) == 255);
    }


    void
    ColorspaceConverterTest::
    testRGB8ToBGRA8()
    {
      ColorspaceConverter<RGB8, BGRA8> converter;
      for(brick::common::UnsignedInt16 redValue = 0; redValue < 256; redValue += 7) {
        for(brick::common::UnsignedInt16 greenValue = 0; greenValue < 256; greenValue += 7) {
          for(brick::common::UnsignedInt16 blueValue = 0; blueValue < 256; blueValue += 7) {
            PixelRGB8 inputPixel(
              static_cast<brick::common::UnsignedInt8>(redValue),
              static_cast<brick::common::UnsignedInt8>(greenValue),
              static_cast<brick::common::UnsignedInt8>(blueValue));
            PixelBGRA8 referencePixel(
              static_cast<brick::common::UnsignedInt8>(blueValue),
              static_cast<brick::common::UnsignedInt8>(greenValue),
              static_cast<brick::common::UnsignedInt8>(redValue),
              255);
            BRICK_TEST_ASSERT(converter(inputPixel) == referencePixel);
          }
        }
      }
    }


    void
    ColorspaceConverterTest::
    testRGB8ToRGBA8()
    {
      ColorspaceConverter<RGB8, RGBA8> converter;
      for(brick::common::UnsignedInt16 redValue = 0; redValue < 256; redValue += 7) {
        for(brick::common::UnsignedInt16 greenValue = 0; greenValue < 256; greenValue += 7) {
          for(brick::common::UnsignedInt16 blueValue = 0; blueValue < 256; blueValue += 7) {
            PixelRGB8 inputPixel(
              static_cast<brick::common::UnsignedInt8>(redValue),
              static_cast<brick::common::UnsignedInt8>(greenValue),
              static_cast<brick::common::UnsignedInt8>(blueValue));
            PixelRGBA8 referencePixel(
              static_cast<brick::common::UnsignedInt8>(redValue),
              static_cast<brick::common::UnsignedInt8>(greenValue),
              static_cast<brick::common::UnsignedInt8>(blueValue),
              255);
            BRICK_TEST_ASSERT(converter(inputPixel) == referencePixel);
          }
        }
      }
    }


    void
    ColorspaceConverterTest::
    testRGB8ToHSV_FLOAT64()
    {
      std::vector<PixelRGB8> inputPixels;
      std::vector< PixelHSV<brick::common::Float64> > targetPixels;

      inputPixels.push_back(PixelRGB8(192, 128, 108));
      inputPixels.push_back(PixelRGB8(177, 192, 108));
      inputPixels.push_back(PixelRGB8(108, 192, 134));
      inputPixels.push_back(PixelRGB8(108, 167, 192));
      inputPixels.push_back(PixelRGB8(132, 108, 192));
      inputPixels.push_back(PixelRGB8(192, 108, 169));

      // HSV values computed manually using:
      //
      //   v = max(r,g,b) / 255.0.
      //   s = (0, if v == 0
      //        (max(r, g, b) - min(r, g, b)) / max, otherwise.
      //   h = {0, if max == min
      //        (60 * (g - b) / (max - min) mod 360) / 360, if max == r
      //        (60 * (b - r) / (max - min) + 120) / 360, if max == g
      //        (60 * (r - g) / (max - min) + 240) / 360, if max == b
      targetPixels.push_back(PixelHSV<brick::common::Float64>(0.0397, 0.4375, 0.7529));
      targetPixels.push_back(PixelHSV<brick::common::Float64>(0.1964, 0.4375, 0.7529));
      targetPixels.push_back(PixelHSV<brick::common::Float64>(0.3849, 0.4375, 0.7529));
      targetPixels.push_back(PixelHSV<brick::common::Float64>(0.5496, 0.4375, 0.7529));
      targetPixels.push_back(PixelHSV<brick::common::Float64>(0.7143, 0.4375, 0.7529));
      targetPixels.push_back(PixelHSV<brick::common::Float64>(0.8790, 0.4375, 0.7529));

      ColorspaceConverter<RGB8, HSV_FLOAT64> converter;

      for(size_t ii = 0; ii < inputPixels.size(); ++ii) {
	PixelHSV<brick::common::Float64> outputPixel = converter(inputPixels[ii]);
	BRICK_TEST_ASSERT(approximatelyEqual(
          outputPixel.hue, targetPixels[ii].hue, 1.0E-4));
	BRICK_TEST_ASSERT(approximatelyEqual(
	  outputPixel.saturation, targetPixels[ii].saturation, 1.0E-4));
	BRICK_TEST_ASSERT(approximatelyEqual(
	  outputPixel.value, targetPixels[ii].value, 1.0E-4));
      }
    }


    void
    ColorspaceConverterTest::
    testRGB8ToYIQ_FLOAT64()
    {
      std::vector<PixelRGB8> inputPixels;

      inputPixels.push_back(PixelRGB8(192, 128, 108));
      inputPixels.push_back(PixelRGB8(177, 192, 108));
      inputPixels.push_back(PixelRGB8(108, 192, 134));
      inputPixels.push_back(PixelRGB8(108, 167, 192));
      inputPixels.push_back(PixelRGB8(132, 108, 192));
      inputPixels.push_back(PixelRGB8(192, 108, 169));

      std::vector<PixelYIQFloat64> targetPixels(inputPixels.size());
      for(unsigned int ii = 0; ii < inputPixels.size(); ++ii) {
        double red = inputPixels[ii].red / 255.0;
        double green = inputPixels[ii].green / 255.0;
        double blue = inputPixels[ii].blue / 255.0;
        targetPixels[ii].luma =
          red * 0.299 + green * 0.587 + blue * 0.114;
        targetPixels[ii].inPhase =
          red * 0.595716 + green * -0.274453 + blue * -0.321263;
        targetPixels[ii].quadrature =
          red * 0.211456 + green * -0.522591 + blue * 0.311135;
      }

      ColorspaceConverter<RGB8, YIQ_FLOAT64> converter;

      const double tolerance = 1.0E-9;
      for(size_t ii = 0; ii < inputPixels.size(); ++ii) {
	PixelYIQ<brick::common::Float64> outputPixel = converter(inputPixels[ii]);
	BRICK_TEST_ASSERT(approximatelyEqual(
          outputPixel.luma, targetPixels[ii].luma, tolerance));
	BRICK_TEST_ASSERT(approximatelyEqual(
	  outputPixel.inPhase, targetPixels[ii].inPhase, tolerance));
	BRICK_TEST_ASSERT(approximatelyEqual(
	  outputPixel.quadrature, targetPixels[ii].quadrature, tolerance));
      }
    }


    void
    ColorspaceConverterTest::
    testRGB_FLOAT64ToHSV_FLOAT64()
    {
      std::vector<PixelRGBFloat64> inputPixels;

      inputPixels.push_back(
        PixelRGBFloat64(192 / 255.0, 128 / 255.0, 108 / 255.0));
      inputPixels.push_back(
        PixelRGBFloat64(177 / 255.0, 192 / 255.0, 108 / 255.0));
      inputPixels.push_back(
        PixelRGBFloat64(108 / 255.0, 192 / 255.0, 134 / 255.0));
      inputPixels.push_back(
        PixelRGBFloat64(108 / 255.0, 167 / 255.0, 192 / 255.0));
      inputPixels.push_back(
        PixelRGBFloat64(132 / 255.0, 108 / 255.0, 192 / 255.0));
      inputPixels.push_back(
        PixelRGBFloat64(192 / 255.0, 108 / 255.0, 169 / 255.0));

      ColorspaceConverter<RGB_FLOAT64, HSV_FLOAT64> converter0;
      ColorspaceConverter<HSV_FLOAT64, RGB_FLOAT64> converter1;

      for(unsigned int ii = 0; ii < inputPixels.size(); ++ii) {
        PixelHSV<brick::common::Float64> hsvPixel = converter0(inputPixels[ii]);
        PixelRGBFloat64 rgbPixel = converter1(hsvPixel);

        const double tolerance = 1.0E-9;
	BRICK_TEST_ASSERT(
          approximatelyEqual(rgbPixel.red, inputPixels[ii].red,
                             tolerance));
	BRICK_TEST_ASSERT(
          approximatelyEqual(rgbPixel.green, inputPixels[ii].green,
                             tolerance));
	BRICK_TEST_ASSERT(
          approximatelyEqual(rgbPixel.blue, inputPixels[ii].blue,
                             tolerance));
      }
    }


    void
    ColorspaceConverterTest::
    testRGB_FLOAT64ToYIQ_FLOAT64()
    {
      std::vector<PixelRGBFloat64> inputPixels;

      inputPixels.push_back(
        PixelRGBFloat64(192 / 255.0, 128 / 255.0, 108 / 255.0));
      inputPixels.push_back(
        PixelRGBFloat64(177 / 255.0, 192 / 255.0, 108 / 255.0));
      inputPixels.push_back(
        PixelRGBFloat64(108 / 255.0, 192 / 255.0, 134 / 255.0));
      inputPixels.push_back(
        PixelRGBFloat64(108 / 255.0, 167 / 255.0, 192 / 255.0));
      inputPixels.push_back(
        PixelRGBFloat64(132 / 255.0, 108 / 255.0, 192 / 255.0));
      inputPixels.push_back(
        PixelRGBFloat64(192 / 255.0, 108 / 255.0, 169 / 255.0));

      std::vector<PixelYIQFloat64> targetPixels(inputPixels.size());
      for(unsigned int ii = 0; ii < inputPixels.size(); ++ii) {
        double red = inputPixels[ii].red;
        double green = inputPixels[ii].green;
        double blue = inputPixels[ii].blue;
        targetPixels[ii].luma =
          red * 0.299 + green * 0.587 + blue * 0.114;
        targetPixels[ii].inPhase =
          red * 0.595716 + green * -0.274453 + blue * -0.321263;
        targetPixels[ii].quadrature =
          red * 0.211456 + green * -0.522591 + blue * 0.311135;
      }

      ColorspaceConverter<RGB_FLOAT64, YIQ_FLOAT64> converter;

      const double tolerance = 1.0E-9;
      for(size_t ii = 0; ii < inputPixels.size(); ++ii) {
	PixelYIQ<brick::common::Float64> outputPixel = converter(inputPixels[ii]);
	BRICK_TEST_ASSERT(approximatelyEqual(
          outputPixel.luma, targetPixels[ii].luma, tolerance));
	BRICK_TEST_ASSERT(approximatelyEqual(
	  outputPixel.inPhase, targetPixels[ii].inPhase, tolerance));
	BRICK_TEST_ASSERT(approximatelyEqual(
	  outputPixel.quadrature, targetPixels[ii].quadrature, tolerance));
      }
    }


    void
    ColorspaceConverterTest::
    testBGRA8ToRGB8()
    {
      ColorspaceConverter<BGRA8, RGB8> converter;
      for(brick::common::UnsignedInt16 redValue = 0; redValue < 256; redValue += 7) {
        for(brick::common::UnsignedInt16 greenValue = 0; greenValue < 256; greenValue += 7) {
          for(brick::common::UnsignedInt16 blueValue = 0; blueValue < 256; blueValue += 7) {
            for(brick::common::UnsignedInt16 alphaValue = 0; alphaValue < 256;
                alphaValue += 7) {
              PixelBGRA8 inputPixel(
                static_cast<brick::common::UnsignedInt8>(blueValue),
                static_cast<brick::common::UnsignedInt8>(greenValue),
                static_cast<brick::common::UnsignedInt8>(redValue),
                static_cast<brick::common::UnsignedInt8>(alphaValue));
              PixelRGB8 referencePixel(
                static_cast<brick::common::UnsignedInt8>(redValue),
                static_cast<brick::common::UnsignedInt8>(greenValue),
                static_cast<brick::common::UnsignedInt8>(blueValue));
              BRICK_TEST_ASSERT(converter(inputPixel) == referencePixel);
            }
          }
        }
      }
    }

    void
    ColorspaceConverterTest::
    testRGBA8ToRGB8()
    {
      ColorspaceConverter<RGBA8, RGB8> converter;
      for(brick::common::UnsignedInt16 redValue = 0; redValue < 256; redValue += 7) {
        for(brick::common::UnsignedInt16 greenValue = 0; greenValue < 256; greenValue += 7) {
          for(brick::common::UnsignedInt16 blueValue = 0; blueValue < 256; blueValue += 7) {
            for(brick::common::UnsignedInt16 alphaValue = 0; alphaValue < 256;
                alphaValue += 7) {
              PixelRGBA8 inputPixel(
                static_cast<brick::common::UnsignedInt8>(redValue),
                static_cast<brick::common::UnsignedInt8>(greenValue),
                static_cast<brick::common::UnsignedInt8>(blueValue),
                static_cast<brick::common::UnsignedInt8>(alphaValue));
              PixelRGB8 referencePixel(
                static_cast<brick::common::UnsignedInt8>(redValue),
                static_cast<brick::common::UnsignedInt8>(greenValue),
                static_cast<brick::common::UnsignedInt8>(blueValue));
              BRICK_TEST_ASSERT(converter(inputPixel) == referencePixel);
            }
          }
        }
      }
    }


    void
    ColorspaceConverterTest::
    testHSV_FLOAT64ToRGB8()
    {
      std::vector<PixelRGBFloat64> inputPixels;

      inputPixels.push_back(PixelRGB8(192, 128, 108));
      inputPixels.push_back(PixelRGB8(177, 192, 108));
      inputPixels.push_back(PixelRGB8(108, 192, 134));
      inputPixels.push_back(PixelRGB8(108, 167, 192));
      inputPixels.push_back(PixelRGB8(132, 108, 192));
      inputPixels.push_back(PixelRGB8(192, 108, 169));

      ColorspaceConverter<RGB8, HSV_FLOAT64> converter0;
      ColorspaceConverter<HSV_FLOAT64, RGB8> converter1;

      for(unsigned int ii = 0; ii < inputPixels.size(); ++ii) {
        PixelHSV<brick::common::Float64> hsvPixel = converter0(inputPixels[ii]);
        PixelRGB8 rgbPixel = converter1(hsvPixel);

	BRICK_TEST_ASSERT(rgbPixel.red == inputPixels[ii].red);
	BRICK_TEST_ASSERT(rgbPixel.green == inputPixels[ii].green);
	BRICK_TEST_ASSERT(rgbPixel.blue == inputPixels[ii].blue);
      }
    }

  } // namespace computerVision

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::computerVision::ColorspaceConverterTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::computerVision::ColorspaceConverterTest currentTest;

}

#endif
