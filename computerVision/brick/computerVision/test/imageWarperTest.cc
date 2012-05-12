/**
***************************************************************************
* @file imageWarperTest.cpp
*
* Source file defining tests for the ImageWarper class template.
*
* Copyright (C) 2006 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/computerVision/imageWarper.hh>
#include <brick/test/testFixture.hh>

namespace num = brick::numeric;

namespace brick {

  namespace computerVision {
    
    class ImageWarperTest
      : public brick::test::TestFixture<ImageWarperTest> {

    public:

      ImageWarperTest();
      ~ImageWarperTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testImageWarper();
      void testImageWarperRGB();

    private:

      template <class FloatType>
      struct ShiftWarpFunctor {
        ShiftWarpFunctor(FloatType xShift, FloatType yShift)
          : m_shift(xShift, yShift) {}
        
        num::Vector2D<FloatType>
        operator()(num::Vector2D<FloatType> const& arg) const {return arg + m_shift;}

        num::Vector2D<FloatType> m_shift;
      };

      template <class FloatType>
      struct StretchXWarpFunctor {
        StretchXWarpFunctor(FloatType stretchFactor)
          : m_factor(stretchFactor) {}
        
        num::Vector2D<FloatType>
        operator()(num::Vector2D<FloatType> const& arg) const {
          return num::Vector2D<FloatType>(arg.x() / m_factor, arg.y());
        }

        FloatType m_factor;
      };

      template <class FloatType>
      struct StretchYWarpFunctor {
        StretchYWarpFunctor(FloatType stretchFactor)
          : m_factor(stretchFactor) {}

        num::Vector2D<FloatType>
        operator()(num::Vector2D<FloatType> const& arg) const {
          return num::Vector2D<FloatType>(arg.x(), arg.y() / m_factor);
        }

        FloatType m_factor;
      };


      template <ImageFormat Format, class FloatType>
      bool
      isInBounds(num::Vector2D<FloatType> const& coordinate,
                 Image<Format> const& image);
      
        
      double m_defaultTolerance;
      float m_defaultFloatTolerance;
      
    }; // class ImageWarperTest


    /* ============== Member Function Definititions ============== */

    ImageWarperTest::
    ImageWarperTest()
      : TestFixture<ImageWarperTest>("ImageWarperTest"),
        m_defaultTolerance(1.0E-10),
        m_defaultFloatTolerance(1.0E-6)
    {
      BRICK_TEST_REGISTER_MEMBER(testImageWarper);
      BRICK_TEST_REGISTER_MEMBER(testImageWarperRGB);
    }


    void
    ImageWarperTest::
    testImageWarper()
    {
      Image<GRAY_FLOAT64> xImage(4, 3);
      Image<GRAY_FLOAT64> yImage(4, 3);

      for(size_t row = 0; row < xImage.rows(); ++row) {
        for(size_t column = 0; column < xImage.columns(); ++column) {
          xImage(row, column) = column;
          yImage(row, column) = row;
        }
      }

      size_t outputRows = 5 * xImage.rows();
      size_t outputColumns = 6 * xImage.rows();
      
      common::Float64 defaultValue = -1.0;
      common::Float64 xShift = -1.2;
      common::Float64 yShift = -2.5;
      common::Float64 scaleFactor = 4.1;
      ShiftWarpFunctor<common::Float64> shiftWarpFunctor(xShift, yShift);
      StretchXWarpFunctor<common::Float64> stretchXWarpFunctor(scaleFactor);
      StretchYWarpFunctor<common::Float64> stretchYWarpFunctor(scaleFactor);
      
      ImageWarper< common::Float64, ShiftWarpFunctor<common::Float64> > shiftWarper(
        xImage.rows(), xImage.columns(), outputRows, outputColumns,
        shiftWarpFunctor);
      ImageWarper< common::Float64, StretchXWarpFunctor<common::Float64> > xStretchWarper(
        xImage.rows(), xImage.columns(), outputRows, outputColumns,
        stretchXWarpFunctor);
      ImageWarper< common::Float64, StretchYWarpFunctor<common::Float64> > yStretchWarper(
        xImage.rows(), xImage.columns(), outputRows, outputColumns,
        stretchYWarpFunctor);

      Image<GRAY_FLOAT64> shiftedXImage =
        shiftWarper.warpImage<GRAY_FLOAT64, GRAY_FLOAT64>(
          xImage, defaultValue);
      Image<GRAY_FLOAT64> shiftedYImage =
        shiftWarper.warpImage<GRAY_FLOAT64, GRAY_FLOAT64>(
          yImage, defaultValue);
      Image<GRAY_FLOAT64> xStretchedXImage =
        xStretchWarper.warpImage<GRAY_FLOAT64, GRAY_FLOAT64>(
          xImage, defaultValue);
      Image<GRAY_FLOAT64> xStretchedYImage =
        xStretchWarper.warpImage<GRAY_FLOAT64, GRAY_FLOAT64>(
          yImage, defaultValue);
      Image<GRAY_FLOAT64> yStretchedXImage =
        yStretchWarper.warpImage<GRAY_FLOAT64, GRAY_FLOAT64>(
          xImage, defaultValue);
      Image<GRAY_FLOAT64> yStretchedYImage =
        yStretchWarper.warpImage<GRAY_FLOAT64, GRAY_FLOAT64>(
          yImage, defaultValue);

      for(size_t row = 0; row < outputRows; ++row) {
        for(size_t column = 0; column < outputColumns; ++column) {
          if(this->isInBounds(
               shiftWarpFunctor(num::Vector2D<common::Float64>(column, row)), xImage)) {
            BRICK_TEST_ASSERT(
              approximatelyEqual(shiftedXImage(row, column), column + xShift,
                                 m_defaultTolerance));
            BRICK_TEST_ASSERT(
              approximatelyEqual(shiftedYImage(row, column), row + yShift,
                                 m_defaultTolerance));
          } else {
            BRICK_TEST_ASSERT(shiftedXImage(row, column) == defaultValue);
            BRICK_TEST_ASSERT(shiftedYImage(row, column) == defaultValue);
          }
          if(this->isInBounds(
               stretchXWarpFunctor(num::Vector2D<common::Float64>(column, row)), xImage)) {
            BRICK_TEST_ASSERT(
              approximatelyEqual(xStretchedXImage(row, column),
                                 column / scaleFactor,
                                 m_defaultTolerance));
            BRICK_TEST_ASSERT(
              approximatelyEqual(xStretchedYImage(row, column),
                                 static_cast<common::Float64>(row),
                                 m_defaultTolerance));
          } else {
            BRICK_TEST_ASSERT(xStretchedXImage(row, column) == defaultValue);
            BRICK_TEST_ASSERT(xStretchedYImage(row, column) == defaultValue);
          }
          if(this->isInBounds(
               stretchYWarpFunctor(num::Vector2D<common::Float64>(column, row)), xImage)) {
            BRICK_TEST_ASSERT(
              approximatelyEqual(yStretchedXImage(row, column),
                                 static_cast<common::Float64>(column),
                                 m_defaultTolerance));
            BRICK_TEST_ASSERT(
              approximatelyEqual(yStretchedYImage(row, column),
                                 row / scaleFactor,
                                 m_defaultTolerance));
          } else {
            BRICK_TEST_ASSERT(yStretchedXImage(row, column) == defaultValue);
            BRICK_TEST_ASSERT(yStretchedYImage(row, column) == defaultValue);
          }
        }
      }
    }


    void
    ImageWarperTest::
    testImageWarperRGB()
    {
      Image<RGB_FLOAT32> xImage(4, 3);
      Image<RGB_FLOAT32> yImage(4, 3);

      for(size_t row = 0; row < xImage.rows(); ++row) {
        for(size_t column = 0; column < xImage.columns(); ++column) {
          xImage(row, column).red = column;
          xImage(row, column).green = column + 5;
          xImage(row, column).blue = 2 * column;
          yImage(row, column).red = row;
          yImage(row, column).green = row + 10;
          yImage(row, column).blue = 3 * row;
        }
      }

      size_t outputRows = 5 * xImage.rows();
      size_t outputColumns = 6 * xImage.rows();
      
      PixelRGBFloat32 defaultValue(-1.0, -1.0, -1.0);
      common::Float32 xShift = -1.2;
      common::Float32 yShift = -2.5;
      common::Float32 scaleFactor = 4.1;
      ShiftWarpFunctor<common::Float32> shiftWarpFunctor(xShift, yShift);
      StretchXWarpFunctor<common::Float32> stretchXWarpFunctor(scaleFactor);
      StretchYWarpFunctor<common::Float32> stretchYWarpFunctor(scaleFactor);
      
      ImageWarper< common::Float32, ShiftWarpFunctor<common::Float32> > shiftWarper(
        xImage.rows(), xImage.columns(), outputRows, outputColumns,
        shiftWarpFunctor);
      ImageWarper< common::Float32, StretchXWarpFunctor<common::Float32> > xStretchWarper(
        xImage.rows(), xImage.columns(), outputRows, outputColumns,
        stretchXWarpFunctor);
      ImageWarper< common::Float32, StretchYWarpFunctor<common::Float32> > yStretchWarper(
        xImage.rows(), xImage.columns(), outputRows, outputColumns,
        stretchYWarpFunctor);

      Image<RGB_FLOAT32> shiftedXImage =
        shiftWarper.warpImage<RGB_FLOAT32, RGB_FLOAT32>(
          xImage, defaultValue);
      Image<RGB_FLOAT32> shiftedYImage =
        shiftWarper.warpImage<RGB_FLOAT32, RGB_FLOAT32>(
          yImage, defaultValue);
      Image<RGB_FLOAT32> xStretchedXImage =
        xStretchWarper.warpImage<RGB_FLOAT32, RGB_FLOAT32>(
          xImage, defaultValue);
      Image<RGB_FLOAT32> xStretchedYImage =
        xStretchWarper.warpImage<RGB_FLOAT32, RGB_FLOAT32>(
          yImage, defaultValue);
      Image<RGB_FLOAT32> yStretchedXImage =
        yStretchWarper.warpImage<RGB_FLOAT32, RGB_FLOAT32>(
          xImage, defaultValue);
      Image<RGB_FLOAT32> yStretchedYImage =
        yStretchWarper.warpImage<RGB_FLOAT32, RGB_FLOAT32>(
          yImage, defaultValue);

      for(size_t row = 0; row < outputRows; ++row) {
        for(size_t column = 0; column < outputColumns; ++column) {
          if(this->isInBounds(
               shiftWarpFunctor(num::Vector2D<common::Float32>(column, row)), xImage)) {
            BRICK_TEST_ASSERT(
              approximatelyEqual(
                shiftedXImage(row, column).red,
                common::Float32(column + xShift),
                m_defaultFloatTolerance));
            BRICK_TEST_ASSERT(
              approximatelyEqual(
                shiftedYImage(row, column).red,
                common::Float32(row + yShift),
                m_defaultFloatTolerance));
          } else {
            BRICK_TEST_ASSERT(shiftedXImage(row, column) == defaultValue);
            BRICK_TEST_ASSERT(shiftedYImage(row, column) == defaultValue);
          }
          if(this->isInBounds(
               stretchXWarpFunctor(num::Vector2D<common::Float32>(column, row)), xImage)) {
            BRICK_TEST_ASSERT(
              approximatelyEqual(
                xStretchedXImage(row, column).red,
                common::Float32(column / scaleFactor),
                m_defaultFloatTolerance));
            BRICK_TEST_ASSERT(
              approximatelyEqual(
                xStretchedYImage(row, column).red,
                static_cast<common::Float32>(row),
                m_defaultFloatTolerance));
          } else {
            BRICK_TEST_ASSERT(xStretchedXImage(row, column) == defaultValue);
            BRICK_TEST_ASSERT(xStretchedYImage(row, column) == defaultValue);
          }
          if(this->isInBounds(
               stretchYWarpFunctor(num::Vector2D<common::Float32>(column, row)),
               xImage)) {
            BRICK_TEST_ASSERT(
              approximatelyEqual(
                yStretchedXImage(row, column).red,
                static_cast<common::Float32>(column),
                m_defaultFloatTolerance));
            BRICK_TEST_ASSERT(
              approximatelyEqual(
                yStretchedYImage(row, column).red,
                common::Float32(row / scaleFactor),
                m_defaultFloatTolerance));
          } else {
            BRICK_TEST_ASSERT(yStretchedXImage(row, column) == defaultValue);
            BRICK_TEST_ASSERT(yStretchedYImage(row, column) == defaultValue);
          }
        }
      }
    }

    template <ImageFormat Format, class FloatType>
    bool
    ImageWarperTest::
    isInBounds(num::Vector2D<FloatType> const& coordinate,
               Image<Format> const& image)
    {
      return ((coordinate.x() >= FloatType(0))
              && (coordinate.y() >= FloatType(0))
              && (coordinate.x() < FloatType((image.columns() - 1)))
              && (coordinate.y() < FloatType((image.rows() - 1))));
    }
    
    
  } // namespace computerVision

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::computerVision::ImageWarperTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::computerVision::ImageWarperTest currentTest;

}

#endif

