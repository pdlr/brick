/**
***************************************************************************
* @file brick/computerVision/nonMaximumSuppress_impl.hh
*
* Header file defining inline and template functions declared in
* nonMaxmumSuppress.hh.
*
* Copyright (C) 2006,2012 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_NONMAXIMUMSUPPRESS_IMPL_HH
#define BRICK_COMPUTERVISION_NONMAXIMUMSUPPRESS_IMPL_HH

// This file is included by nonMaximumSuppress.hh, and should not be
// directly included by user code, so no need to include
// nonMaximumSuppress.hh here.
//
// #include <brick/computerVision/nonMaximumSuppress.hh>

#include <brick/common/mathFunctions.hh>

namespace brick {

  namespace computerVision {


    // This function zeros any pixels of the input image which are not
    // plausible edges.
    template <class FloatType, ImageFormat FORMAT>
    Image<FORMAT>
    nonMaximumSuppress(const Image<FORMAT>& inputImage,
                       const brick::numeric::Array2D<FloatType>& gradX,
                       const brick::numeric::Array2D<FloatType>& gradY)
    {
      // Argument checking.
      if(inputImage.rows() == 0 || inputImage.columns() == 0) {
        BRICK_THROW(brick::common::ValueException, "nonMaximumSuppress()",
                    "Argument inputImage must have non-zero size.");
      }
      if(inputImage.rows() != gradX.rows()
         || inputImage.columns() < gradX.columns()) {
        BRICK_THROW(brick::common::ValueException, "nonMaximumSuppress()",
                    "Arguments inputImage and gradX must have the same shape.");
      }
      if(inputImage.rows() != gradY.rows()
         || inputImage.columns() < gradY.columns()) {
        BRICK_THROW(brick::common::ValueException, "nonMaximumSuppress()",
                    "Arguments inputImage and gradY must have the same shape.");
      }

      // Create an output image.
      Image<FORMAT> suppressedImage(inputImage.rows(), inputImage.columns());
      suppressedImage = static_cast<typename Image<FORMAT>::PixelType>(0);

      // Test each pixel individually.
      size_t columns = inputImage.columns();
      size_t rowsMinusOne = inputImage.rows() - 1;
      size_t columnsMinusOne = inputImage.columns() - 1;
      for(size_t row = 1; row < rowsMinusOne; ++row) {
        size_t index0 = row * inputImage.columns() + 1;
        for(size_t column = 1; column < columnsMinusOne; ++column) {
          if(inputImage[index0]) {
            double gradXComponent = gradX[index0];
            double gradYComponent = gradY[index0];
            size_t neighbor0Index;
            size_t neighbor1Index;
            if(gradXComponent == 0.0) {
              neighbor0Index = index0 + columns;
              neighbor1Index = index0 - columns;
            } else if(brick::common::absoluteValue(gradXComponent)
                      >= brick::common::absoluteValue(gradYComponent)) {
              double indicator = gradYComponent / gradXComponent;
              if(indicator >= 0.5) {
                neighbor0Index = index0 + columns + 1;
                neighbor1Index = index0 - columns - 1;
              } else if(indicator < -0.5) {
                neighbor0Index = index0 + columns - 1;
                neighbor1Index = index0 - columns + 1;
              } else {
                neighbor0Index = index0 + 1;
                neighbor1Index = index0 - 1;
              }
            } else { // fabs(gradYComponent) > fabs(gradXComponent)
              double indicator = gradXComponent / gradYComponent;
              if(indicator >= 0.5) {
                neighbor0Index = index0 + columns + 1;
                neighbor1Index = index0 - columns - 1;
              } else if(indicator < -0.5) {
                neighbor0Index = index0 + columns - 1;
                neighbor1Index = index0 - columns + 1;
              } else {
                neighbor0Index = index0 + columns;
                neighbor1Index = index0 - columns;
              }
            }
            if(inputImage[index0] > inputImage[neighbor0Index]
               && inputImage[index0] > inputImage[neighbor1Index]) {
              suppressedImage[index0] = inputImage[index0];
            }
          } // if(inputImage[index0] != 0.0)
          ++index0;
        } // for(size_t column...)
      } // for(size_t row...)

      return suppressedImage;
    } // nonMaximumSuppress()


  } // namespace computerVision

} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_KERNEL_IMPL_HH */
