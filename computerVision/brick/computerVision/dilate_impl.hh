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
#include <brick/numeric/stencil2D.hh>

namespace brick {

  namespace computerVision {

    // Private functions that will support
    // dilate(Image, unsigned int, unsigned int)
    namespace privateCode {

      template<ImageFormat FORMAT>
      void
      dilateBottomBorder(const Image<FORMAT>& inputImage, size_t radius,
                         Image<FORMAT>& outputImage)
      {
        typedef typename Image<FORMAT>::value_type ValueType;

        // Some constants to help with loops below.
        const size_t columns = inputImage.columns();
        const size_t columnsMinusRadius = inputImage.columns() - radius;
        const size_t rows = inputImage.rows();
        const size_t rowsMinusRadius = inputImage.rows() - radius;

        // This variable will scan the entire border in raster order.
        size_t outputIndex = rowsMinusRadius * columns;
        
        // Bottom border is made up of last radius rows.
        for(size_t row = rowsMinusRadius; row < rows; ++row) {

          // First radius columns are a special case.
          for(size_t column = 0; column < radius; ++column) {
            outputImage[outputIndex] = ValueType(0);
            for(size_t row2 = row - radius; row2 < rows; ++row2) {
              for(size_t column2 = 0; column2 <= column + radius; ++column2) {
                if(inputImage(row2, column2)) {
                  outputImage[outputIndex] = ValueType(1);
                  break;
                }
              }
              if(outputImage[outputIndex]) {
                break;
              }
            }
            ++outputIndex;
          }

          for(size_t column = radius; column < columnsMinusRadius; ++column) {
            outputImage[outputIndex] = ValueType(0);
            for(size_t row2 = row - radius; row2 < rows; ++row2) {
              for(size_t column2 = column - radius; column2 <= column + radius;
                  ++column2) {
                if(inputImage(row2, column2)) {
                  outputImage[outputIndex] = ValueType(1);
                  break;
                }
              }
              if(outputImage[outputIndex]) {
                break;
              }
            }
            ++outputIndex;
          }


          // Last radius columns are a special case.
          for(size_t column = columnsMinusRadius; column < columns; ++column) {
            outputImage[outputIndex] = ValueType(0);
            for(size_t row2 = row - radius; row2 < rows; ++row2) {
              for(size_t column2 = column - radius; column2 < columns;
                  ++column2) {
                if(inputImage(row2, column2)) {
                  outputImage[outputIndex] = ValueType(1);
                  break;
                }
              }
              if(outputImage[outputIndex]) {
                break;
              }
            }
            ++outputIndex;
          }
        }
      }
      

      template<ImageFormat FORMAT>
      void
      dilateTopBorder(const Image<FORMAT>& inputImage, size_t radius,
                      Image<FORMAT>& outputImage)
      {
        typedef typename Image<FORMAT>::value_type ValueType;

        // Some constants to help with loops below.
        const size_t columns = inputImage.columns();
        const size_t columnsMinusRadius = inputImage.columns() - radius;
      
        // This variable will scan the entire output image in raster order.
        size_t outputIndex = 0;

        // The top border consists of the first radius rows.
        for(size_t row = 0; row < radius; ++row) {

          // First radius columns are a special case.
          for(size_t column = 0; column < radius; ++column) {
            outputImage[outputIndex] = ValueType(0);
            for(size_t row2 = 0; row2 <= row + radius; ++row2) {
              for(size_t column2 = 0; column2 <= column + radius; ++column2) {
                if(inputImage(row2, column2)) {
                  outputImage[outputIndex] = ValueType(1);
                  break;
                }
              }
              if(outputImage[outputIndex]) {
                break;
              }
            }
            ++outputIndex;
          }

          for(size_t column = radius; column < columnsMinusRadius; ++column) {
            outputImage[outputIndex] = ValueType(0);
            for(size_t row2 = 0; row2 <= row + radius; ++row2) {
              for(size_t column2 = column - radius; column2 <= column + radius;
                  ++column2) {
                if(inputImage(row2, column2)) {
                  outputImage[outputIndex] = ValueType(1);
                  break;
                }
              }
              if(outputImage[outputIndex]) {
                break;
              }
            }
            ++outputIndex;
          }
        
          // Last radius columns are a special case.
          for(size_t column = columnsMinusRadius; column < columns; ++column) {
            outputImage[outputIndex] = ValueType(0);
            for(size_t row2 = 0; row2 <= row + radius; ++row2) {
              for(size_t column2 = column - radius; column2 < columns;
                  ++column2) {
                if(inputImage(row2, column2)) {
                  outputImage[outputIndex] = ValueType(1);
                  break;
                }
              }
              if(outputImage[outputIndex]) {
                break;
              }
            }
            ++outputIndex;
          }
        }
      }

      
      template<ImageFormat FORMAT, size_t StencilSize>
      void
      sizedDilate(const Image<FORMAT>& inputImage, size_t radius,
                  Image<FORMAT>& outputImage)
      {
        typedef typename Image<FORMAT>::value_type ValueType;

        // Some constants to help with loops below.
        const size_t columns = inputImage.columns();
        const size_t columnsMinusRadius = inputImage.columns() - radius;
        // const size_t rows = inputImage.rows();
        const size_t rowsMinusRadius = inputImage.rows() - radius;

        // A stencil to avoid quadruple-looping over most of the image.
        typedef typename ImageFormatTraits<FORMAT>::PixelType PixType;
        brick::numeric::Stencil2D<const PixType, StencilSize> inputStencil(
          2 * radius + 1, 2 * radius + 1);
        inputStencil.setTarget(inputImage);
        typedef brick::numeric::StencilIterator<const PixType, StencilSize>
          InputIterator;

        // Now do the bulk of the image (up to the last radius rows).
        for(size_t row = radius; row < rowsMinusRadius; ++row) {
          size_t outputIndex = row * outputImage.columns();
          
          // First radius columns are a special case.
          for(size_t column = 0; column < radius; ++column) {
            outputImage[outputIndex] = ValueType(0);
            for(size_t row2 = row - radius; row2 <= row + radius; ++row2) {
              for(size_t column2 = 0; column2 <= column + radius; ++column2) {
                if(inputImage(row2, column2)) {
                  outputImage[outputIndex] = ValueType(1);
                  break;
                }
              }
              if(outputImage[outputIndex]) {
                break;
              }
            }
            ++outputIndex;
          }
          
          inputStencil.goTo(row - radius, 0);
          for(size_t column = radius; column < columnsMinusRadius; ++column) {
            InputIterator inputIterator = inputStencil.begin(); 
            InputIterator endIterator = inputStencil.end(); 
            outputImage[outputIndex] = ValueType(0);
            while(inputIterator != endIterator) {
              if(*inputIterator) {
                outputImage[outputIndex] = ValueType(1);
                break;
              }
              ++inputIterator;
            }
            inputStencil.advance();
            ++outputIndex;
          }


          // Last radius columns are a special case.
          for(size_t column = columnsMinusRadius; column < columns; ++column) {
            outputImage[outputIndex] = ValueType(0);
            for(size_t row2 = row - radius; row2 <= row + radius; ++row2) {
              for(size_t column2 = column - radius; column2 < columns;
                  ++column2) {
                if(inputImage(row2, column2)) {
                  outputImage[outputIndex] = ValueType(1);
                  break;
                }
              }
              if(outputImage[outputIndex]) {
                break;
              }
            }
            ++outputIndex;
          }
        }
      }

    } // namespace privateCode
    
  
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
    dilate(const Image<FORMAT>& inputImage, size_t radius)
    {
      typedef typename Image<FORMAT>::value_type ValueType;

      if(radius >= inputImage.rows() || radius >= inputImage.columns()) {
        BRICK_THROW(brick::common::ValueException, "dilate()",
                  "Argument radius must be less than both the width and "
                  "height of the input image.");
      }
      
      Image<FORMAT> outputImage(inputImage.rows(), inputImage.columns());

      privateCode::dilateTopBorder(inputImage, radius, outputImage);

      // Now do the bulk of the image (up to the last radius rows).
      if(radius < 5) {
        privateCode::sizedDilate<FORMAT, 81>(inputImage, radius, outputImage);
      } else if(radius < 10) {
        privateCode::sizedDilate<FORMAT, 361>(inputImage, radius, outputImage);
      } else if(radius < 20) {
        privateCode::sizedDilate<FORMAT, 1521>(inputImage, radius, outputImage);
      } else {
        BRICK_THROW(brick::common::NotImplementedException, "dilate()",
                  "Dilations with radius >= 20 are not currently supported.");
      }

      privateCode::dilateBottomBorder(inputImage, radius, outputImage);
      
      return outputImage;
    }
    
  } // namespace computerVision
    
} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_DILATE_IMPL_HH */
