/**
***************************************************************************
* @file brick/computerVision/dilate.hh
*
* Header file declaring the dilate() function template.
*
* Copyright (C) 2006 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_DILATE_IMPL_HH
#define BRICK_COMPUTERVISION_DILATE_IMPL_HH

// This file is included by dilate.hh, and should not be directly included
// by user code, so no need to include dilate.hh here.
// 
// #include <brick/computerVision/dilate.hh>

#include <cmath>
#include <brick/computerVision/erode.hh> // For privateCode::CountingFunctor.

namespace brick {

  namespace computerVision {
  
    template<ImageFormat FORMAT>
    Image<FORMAT>
    dilate(const Image<FORMAT>& inputImage)
    {
      typedef typename Image<FORMAT>::value_type ValueType;
    
      Image<FORMAT> outputImage(inputImage.rows(), inputImage.columns());

      size_t index0 = 0;
      if(inputImage[0]
         || inputImage[1]
         || inputImage[inputImage.columns()]
         || inputImage[inputImage.columns() + 1]) {
        outputImage[index0] = ValueType(1);
      } else {
        outputImage[index0] = ValueType(0);
      }
      ++index0;

      for(size_t column = 1; column < inputImage.columns() - 1; ++column) {
        if(inputImage[index0]
           || inputImage[index0 - 1]
           || inputImage[index0 + 1]
           || inputImage[index0 + inputImage.columns() - 1]
           || inputImage[index0 + inputImage.columns()]
           || inputImage[index0 + inputImage.columns() + 1]) {
          outputImage[index0] = ValueType(1);
        } else {
          outputImage[index0] = ValueType(0);
        }
        ++index0;
      }
        
      if(inputImage[index0]
         || inputImage[index0 - 1]
         || inputImage[index0 + inputImage.columns() - 1]
         || inputImage[index0 + inputImage.columns()]) {
        outputImage[index0] = ValueType(1);
      } else {
        outputImage[index0] = ValueType(0);
      }
      ++index0;

      for(size_t row = 1; row < inputImage.rows() - 1; ++row) {
        if(inputImage[index0]
           || inputImage[index0 + 1]
           || inputImage[index0 - inputImage.columns()]
           || inputImage[index0 - inputImage.columns() + 1]
           || inputImage[index0 + inputImage.columns()]
           || inputImage[index0 + inputImage.columns() + 1]) {
          outputImage[index0] = ValueType(1);
        } else {
          outputImage[index0] = ValueType(0);
        }
        ++index0;

        for(size_t column = 1; column < inputImage.columns() - 1; ++column) {
          if(inputImage[index0]
             || inputImage[index0 - 1]
             || inputImage[index0 + 1]
             || inputImage[index0 - inputImage.columns() - 1]
             || inputImage[index0 - inputImage.columns()]
             || inputImage[index0 - inputImage.columns() + 1]
             || inputImage[index0 + inputImage.columns() - 1]
             || inputImage[index0 + inputImage.columns()]
             || inputImage[index0 + inputImage.columns() + 1]) {
            outputImage[index0] = ValueType(1);
          } else {
            outputImage[index0] = ValueType(0);
          }
          ++index0;
        }

        if(inputImage[index0]
           || inputImage[index0 - 1]
           || inputImage[index0 - inputImage.columns()]
           || inputImage[index0 - inputImage.columns() - 1]
           || inputImage[index0 + inputImage.columns()]
           || inputImage[index0 + inputImage.columns() - 1]) {
          outputImage[index0] = ValueType(1);
        } else {
          outputImage[index0] = ValueType(0);
        }
        ++index0;
      }

      if(inputImage[index0]
         || inputImage[index0 + 1]
         || inputImage[index0 - inputImage.columns()]
         || inputImage[index0 - inputImage.columns() + 1]) {
        outputImage[index0] = ValueType(1);
      } else {
        outputImage[index0] = ValueType(0);
      }
      ++index0;

      for(size_t column = 1; column < inputImage.columns() - 1; ++column) {
        if(inputImage[index0]
           || inputImage[index0 - 1]
           || inputImage[index0 + 1]
           || inputImage[index0 - inputImage.columns() - 1]
           || inputImage[index0 - inputImage.columns()]
           || inputImage[index0 - inputImage.columns() + 1]) {
          outputImage[index0] = ValueType(1);
        } else {
          outputImage[index0] = ValueType(0);
        }
        ++index0;
      }
        
      if(inputImage[index0]
         || inputImage[index0 - 1]
         || inputImage[index0 - inputImage.columns() - 1]
         || inputImage[index0 - inputImage.columns()]) {
        outputImage[index0] = ValueType(1);
      } else {
        outputImage[index0] = ValueType(0);
      }
      ++index0;

      return outputImage;
    }


    template<ImageFormat FORMAT>
    Image<FORMAT>
    dilateUsingBoxIntegrator(const Image<FORMAT>& inputImage,
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
      
      brick::numeric::BoxIntegrator2D<ValueType, int> integrator(
        inputImage, privateCode::CountingFunctor<ValueType>());
      Image<FORMAT> outputImage(inputImage.rows(), inputImage.columns());

      size_t index0 = 0;
      size_t row = 0;
      size_t const rowBoundary0 = windowRadiusH;  // Integer division.
      size_t const rowBoundary1 = inputImage.rows() - windowRadiusH;
      size_t const colBoundary0 = windowRadiusW;  // Integer division.
      size_t const colBoundary1 = inputImage.columns() - windowRadiusW;
      for(; row < rowBoundary0; ++row) {
        size_t column = 0;
        for(; column < colBoundary0; ++column) {
          if(integrator.getIntegral(
               brick::numeric::Index2D(0, 0),
               brick::numeric::Index2D(row + windowRadiusH + 1,
                                       column + windowRadiusW + 1))
             != 0) {
            outputImage[index0] = ValueType(1);
          } else {
            outputImage[index0] = ValueType(0);
          }
          ++index0;
        }        
        for(; column < colBoundary1; ++column) {
          if(integrator.getIntegral(
               brick::numeric::Index2D(0, column - windowRadiusW),
               brick::numeric::Index2D(row + windowRadiusH + 1,
                                       column + windowRadiusW + 1))
             != 0) {
            outputImage[index0] = ValueType(1);
          } else {
            outputImage[index0] = ValueType(0);
          }
          ++index0;
        }
        for(; column < inputImage.columns(); ++column) {
          if(integrator.getIntegral(
               brick::numeric::Index2D(0, column - windowRadiusW),
               brick::numeric::Index2D(row + windowRadiusH + 1,
                                       inputImage.columns()))
             != 0) {
            outputImage[index0] = ValueType(1);
          } else {
            outputImage[index0] = ValueType(0);
          }
          ++index0;
        }
      }


      for(; row < rowBoundary1; ++row) {
        size_t column = 0;
        for(; column < colBoundary0; ++column) {
          if(integrator.getIntegral(
               brick::numeric::Index2D(row - windowRadiusH, 0),
               brick::numeric::Index2D(row + windowRadiusH + 1,
                                       column + windowRadiusW + 1))
             != 0) {
            outputImage[index0] = ValueType(1);
          } else {
            outputImage[index0] = ValueType(0);
          }
          ++index0;
        }        
        for(; column < colBoundary1; ++column) {
          if(integrator.getIntegral(
               brick::numeric::Index2D(row - windowRadiusH,
                                       column - windowRadiusW),
               brick::numeric::Index2D(row + windowRadiusH + 1,
                                       column + windowRadiusW + 1))
             != 0) {
            outputImage[index0] = ValueType(1);
          } else {
            outputImage[index0] = ValueType(0);
          }
          ++index0;
        }
        for(; column < inputImage.columns(); ++column) {
          if(integrator.getIntegral(
               brick::numeric::Index2D(row - windowRadiusH,
                                       column - windowRadiusW),
               brick::numeric::Index2D(row + windowRadiusH + 1,
                                       inputImage.columns()))
             != 0) {
            outputImage[index0] = ValueType(1);
          } else {
            outputImage[index0] = ValueType(0);
          }
          ++index0;
        }
      }

      for(; row < inputImage.rows(); ++row) {
        size_t column = 0;
        for(; column < colBoundary0; ++column) {
          if(integrator.getIntegral(
               brick::numeric::Index2D(row - windowRadiusH, 0),
               brick::numeric::Index2D(inputImage.rows(),
                                       column + windowRadiusW + 1))
             != 0) {
            outputImage[index0] = ValueType(1);
          } else {
            outputImage[index0] = ValueType(0);
          }
          ++index0;
        }        
        for(; column < colBoundary1; ++column) {
          if(integrator.getIntegral(
               brick::numeric::Index2D(row - windowRadiusH,
                                       column - windowRadiusW),
               brick::numeric::Index2D(inputImage.rows(),
                                       column + windowRadiusW + 1))
             != 0) {
            outputImage[index0] = ValueType(1);
          } else {
            outputImage[index0] = ValueType(0);
          }
          ++index0;
        }
        for(; column < inputImage.columns(); ++column) {
          if(integrator.getIntegral(
               brick::numeric::Index2D(row - windowRadiusH,
                                       column - windowRadiusW),
               brick::numeric::Index2D(inputImage.rows(),
                                       inputImage.columns()))
             != 0) {
            outputImage[index0] = ValueType(1);
          } else {
            outputImage[index0] = ValueType(0);
          }
          ++index0;
        }
      }      

      return outputImage;
    }

    
  } // namespace computerVision
    
} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_DILATE_IMPL_HH */
