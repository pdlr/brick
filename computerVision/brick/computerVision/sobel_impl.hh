/**
***************************************************************************
* @file brick/computerVision/sobel_impl.hh
*
* Header file defining inline and template functions from sobel.hh.
*
* Copyright (C) 2006,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_SOBEL_IMPL_HH
#define BRICK_COMPUTERVISION_SOBEL_IMPL_HH

// This file is included by sobel.hh, and should not be directly included
// by user code, so no need to include sobel.hh here.
// 
// #include <brick/computerVision/sobel.hh>

namespace brick {

  namespace computerVision {

    // This function applies the sobel edge operator in the X
    // direction.
    template <ImageFormat FORMAT>
    Image<FORMAT>
    applySobelX(const Image<FORMAT>& inputImage, bool normalizeResult)
    {
      // Argument checking.
      if(normalizeResult == true) {
        BRICK_THROW(brick::common::NotImplementedException, "applySobelX()",
                  "Argument normalizeResult must be false for now.");
      }
      if(inputImage.rows() < 2 || inputImage.columns() < 2) {
        BRICK_THROW(brick::common::ValueException, "applySobelX()",
                  "Argument inputImage must be 2x2 or larger.");
      }

      // Prepare a space for the result.
      Image<FORMAT> gradientImage(inputImage.rows(), inputImage.columns());

      // We'll keep a running index to avoid double-indexing the images.
      size_t index0 = 0;
      const size_t rows = inputImage.rows();
      const size_t cols = inputImage.columns();
      const size_t colsMinusOne = cols - 1;

      // NOTE: In all of the code below, We assume that the image
      // continues with constant first derivative into the (fictional)
      // out-of-bounds rows and columns.  This has the somewhat
      // surprising effect of making the first row of convolution
      // results independent of the second image row.
      
      // Use a very reduced kernel for the upper-left corner.
      gradientImage[index0] = (
        8 * (inputImage[index0 + 1] - inputImage[index0]));
      ++index0;
      
      // Use a reduced kernel for the first row.
      for(size_t column = 1; column < colsMinusOne; ++column) {
        gradientImage[index0] = (
          4 * (inputImage[index0 + 1] - inputImage[index0 - 1]));
        ++index0;
      }

      // Use a very reduced kernel for the upper-right corner.
      gradientImage[index0] = (
        8 * (inputImage[index0] - inputImage[index0 - 1]));
      ++index0;
      
      // Convolve the bulk of the image.
      for(size_t row = 1; row < rows - 1; ++row) {

        // Use a reduced kernel for first pixel in the row.
        gradientImage[index0] = (
          2 * (inputImage[index0 - cols + 1] - inputImage[index0 - cols])
          + 4 * (inputImage[index0 + 1] - inputImage[index0])
          + 2 * (inputImage[index0 + cols + 1] - inputImage[index0 + cols]));
        ++index0;
        
        // Convolve the bulk of the row.
        for(size_t column = 1; column < colsMinusOne; ++column) {
          gradientImage[index0] = (
            (inputImage[index0 - cols + 1] - inputImage[index0 - cols - 1])
            + 2 * (inputImage[index0 + 1] - inputImage[index0 - 1])
            + (inputImage[index0 + cols + 1] - inputImage[index0 + cols - 1]));
          ++index0;
        }

        // Use a reduced kernel for last pixel in the row.
        gradientImage[index0] = (
          2 * (inputImage[index0 - cols] - inputImage[index0 - cols - 1])
          + 4 * (inputImage[index0] - inputImage[index0 - 1])
          + 2 * (inputImage[index0 + cols] - inputImage[index0 + cols - 1]));
        ++index0;
      }

      // Use a very reduced kernel for the lower-left corner.
      gradientImage[index0] = (
        8 * (inputImage[index0 + 1] - inputImage[index0]));
      ++index0;
      
      // Use a reduced kernel for the last row.
      for(size_t column = 1; column < colsMinusOne; ++column) {
        gradientImage[index0] = (
          4 * (inputImage[index0 + 1] - inputImage[index0 - 1]));
        ++index0;
      }

      // Use a very reduced kernel for the lower-right corner.
      gradientImage[index0] = (
        8 * (inputImage[index0] - inputImage[index0 - 1]));
      ++index0;

      return gradientImage;
    }


    // This function applies the sobel edge operator in the Y
    // direction.
    template <ImageFormat FORMAT>
    Image<FORMAT>
    applySobelY(const Image<FORMAT>& inputImage, bool normalizeResult)
    {
      // Argument checking.
      if(normalizeResult == true) {
        BRICK_THROW(brick::common::NotImplementedException, "applySobelY()",
                  "Argument normalizeResult must be false for now.");
      }
      if(inputImage.rows() < 2 || inputImage.columns() < 2) {
        BRICK_THROW(brick::common::ValueException, "applySobelY()",
                  "Argument inputImage must be 2x2 or larger.");
      }

      // Prepare a space for the result.
      Image<FORMAT> gradientImage(inputImage.rows(), inputImage.columns());

      // We'll keep a running index to avoid double-indexing the images.
      size_t index0 = 0;
      const size_t rows = inputImage.rows();
      const size_t cols = inputImage.columns();
      const size_t colsMinusOne = cols - 1;

      // NOTE: In all of the code below, We assume that the image
      // continues with constant first derivative into the (fictional)
      // out-of-bounds rows and columns.  This has the somewhat
      // surprising effect of making the first row of convolution
      // results independent of the second image row.
      
      // Use a very reduced kernel for the upper-left corner.
      gradientImage[index0] = (
        8 * (inputImage[index0 + cols] - inputImage[index0]));
      ++index0;
      
      // Use a reduced kernel for the first row.
      for(size_t column = 1; column < colsMinusOne; ++column) {
        gradientImage[index0] = (
          2 * (inputImage[index0 + cols - 1] - inputImage[index0 - 1])
          + 4 * (inputImage[index0 + cols] - inputImage[index0])
          + 2 * (inputImage[index0 + cols + 1] - inputImage[index0 + 1]));
        ++index0;
      }

      // Use a very reduced kernel for the upper-right corner.
      gradientImage[index0] = (
        8 * (inputImage[index0 + cols] - inputImage[index0]));
      ++index0;
      
      // Convolve the bulk of the image.
      for(size_t row = 1; row < rows - 1; ++row) {

        // Use a reduced kernel for first pixel in the row.
        gradientImage[index0] = (
          4 * (inputImage[index0 + cols] - inputImage[index0 - cols]));
        ++index0;
        
        // Convolve the bulk of the row.
        for(size_t column = 1; column < colsMinusOne; ++column) {
          gradientImage[index0] = (
            (inputImage[index0 + cols - 1] - inputImage[index0 - cols - 1])
            + 2 * (inputImage[index0 + cols] - inputImage[index0 - cols])
            + (inputImage[index0 + cols + 1] - inputImage[index0 - cols + 1]));
          ++index0;
        }

        // Use a reduced kernel for last pixel in the row.
        gradientImage[index0] = (
          4 * (inputImage[index0 + cols] - inputImage[index0 - cols]));
        ++index0;
      }

      // Use a very reduced kernel for the lower-left corner.
      gradientImage[index0] = (
        8 * (inputImage[index0] - inputImage[index0 - cols]));
      ++index0;
      
      // Use a reduced kernel for the last row.
      for(size_t column = 1; column < colsMinusOne; ++column) {
        gradientImage[index0] = (
          2 * (inputImage[index0 - 1] - inputImage[index0 - cols - 1])
          + 4 * (inputImage[index0] - inputImage[index0 - cols])
          + 2 * (inputImage[index0 + 1] - inputImage[index0 - cols + 1]));
        ++index0;
      }

      // Use a very reduced kernel for the lower-right corner.
      gradientImage[index0] = (
        8 * (inputImage[index0] - inputImage[index0 - cols]));
      ++index0;

      return gradientImage;
    }


  } // namespace computerVision
  
} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_KERNEL_IMPL_HH */
