/**
***************************************************************************
* @file brick/computerVision/erode.hh
*
* Header file defining inline and template functions declared in erode.hh.
*
* Copyright (C) 2006,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_ERODE_IMPL_HH
#define BRICK_COMPUTERVISION_ERODE_IMPL_HH

// This file is included by erode.hh, and should not be directly included
// by user code, so no need to include erode.hh here.
//
// #include <brick/computerVision/erode.hh>

#include <cmath>
#include <brick/numeric/boxIntegrator2D.hh>

namespace brick {

  namespace computerVision {

    // Private functor that will support erodeUsingBoxIntegrator().
    namespace privateCode {

      template <class Type>
      struct CountingFunctor {
        int operator()(Type const& arg0) {
          return arg0 ? 1 : 0;
        }
      };

    } // namespace privateCode


    template<ImageFormat FORMAT>
    Image<FORMAT>
    erode(const Image<FORMAT>& inputImage)
    {
      typedef typename Image<FORMAT>::value_type ValueType;

      Image<FORMAT> outputImage(inputImage.rows(), inputImage.columns());

      size_t index0 = 0;
      size_t row = 0;
      size_t rowBoundary0 = 1;
      size_t rowBoundary1 = inputImage.rows() - 1;
      for(; row < rowBoundary0; ++row) {
        for(size_t column = 0; column < inputImage.columns(); ++column) {
          outputImage[index0] = ValueType(0);
          ++index0;
        }
      }


      size_t colBoundary0 = 1;
      size_t colBoundary1 = inputImage.columns() - 1;
      for(; row < rowBoundary1; ++row) {
        size_t column = 0;
        for(; column < colBoundary0; ++column) {
          outputImage[index0] = ValueType(0);
          ++index0;
        }
        for(; column < colBoundary1; ++column) {
          if(inputImage[index0 - 1]
             && inputImage[index0 + 1]
             && inputImage[index0 - inputImage.columns() - 1]
             && inputImage[index0 - inputImage.columns()]
             && inputImage[index0 - inputImage.columns() + 1]
             && inputImage[index0 + inputImage.columns() - 1]
             && inputImage[index0 + inputImage.columns()]
             && inputImage[index0 + inputImage.columns() + 1]) {
            outputImage[index0] = inputImage[index0];
          } else {
            outputImage[index0] = ValueType(0);
          }
          ++index0;
        }
        for(; column < inputImage.columns(); ++column) {
          outputImage[index0] = ValueType(0);
          ++index0;
        }
      }

      for(; row < inputImage.rows(); ++row) {
        for(size_t column = 0; column < inputImage.columns(); ++column) {
          outputImage[index0] = ValueType(0);
          ++index0;
        }
      }

      return outputImage;
    }


    // This function will almost certainly be faster for larger
    // windows, but doesn't appear to be faster than the hard-coded
    // 3x3 version in our test case.
    template<ImageFormat FORMAT>
    Image<FORMAT>
    erodeUsingBoxIntegrator(const Image<FORMAT>& inputImage,
                            unsigned int windowWidth,
                            unsigned int windowHeight)
    {
      typedef typename Image<FORMAT>::value_type ValueType;

      // The code below assumes odd window sizes.  Silently arrange
      // for that to be true.
      if(windowHeight % 2 == 0) {
        ++windowHeight;
      }
      if(windowWidth % 2 == 0) {
        ++windowWidth;
      }

      // Constants to make code below more readable.
      unsigned int const windowRadiusH = windowHeight / 2;
      unsigned int const windowRadiusW = windowWidth / 2;
      int const regionSize = windowHeight * windowWidth;

      brick::numeric::BoxIntegrator2D<ValueType, int> integrator(
        inputImage, privateCode::CountingFunctor<ValueType>());
      Image<FORMAT> outputImage(inputImage.rows(), inputImage.columns());

      size_t index0 = 0;
      size_t row = 0;
      size_t const rowBoundary0 = windowRadiusH;  // Integer division.
      size_t const rowBoundary1 = inputImage.rows() - windowRadiusH;
      for(; row < rowBoundary0; ++row) {
        for(size_t column = 0; column < inputImage.columns(); ++column) {
          outputImage[index0] = ValueType(0);
          ++index0;
        }
      }


      size_t const colBoundary0 = windowRadiusW;  // Integer division.
      size_t colBoundary1 = inputImage.columns() - windowRadiusW;
      for(; row < rowBoundary1; ++row) {
        size_t column = 0;
        for(; column < colBoundary0; ++column) {
          outputImage[index0] = ValueType(0);
          ++index0;
        }
        for(; column < colBoundary1; ++column) {
          if(integrator.getIntegral(
               brick::numeric::Index2D(row - windowRadiusH,
                                       column - windowRadiusW),
               brick::numeric::Index2D(row + windowRadiusH + 1,
                                       column + windowRadiusW + 1))
             == regionSize) {
            outputImage[index0] = inputImage[index0];
          } else {
            outputImage[index0] = ValueType(0);
          }
          ++index0;
        }
        for(; column < inputImage.columns(); ++column) {
          outputImage[index0] = ValueType(0);
          ++index0;
        }
      }

      for(; row < inputImage.rows(); ++row) {
        for(size_t column = 0; column < inputImage.columns(); ++column) {
          outputImage[index0] = ValueType(0);
          ++index0;
        }
      }

      return outputImage;
    }

  } // namespace computerVision

} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_ERODE_IMPL_HH */
