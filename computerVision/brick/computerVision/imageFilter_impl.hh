/**
***************************************************************************
* @file brick/computerVision/imageFilter_impl.hh
*
* Header file defining functions for performing image filtering.
*
* Copyright (C) 2006, 2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_IMAGEFILTER_IMPL_HH
#define BRICK_COMPUTERVISION_IMAGEFILTER_IMPL_HH

// This file is included by imageFilter.hh, and should not be directly included
// by user code, so no need to include imageFilter.hh here.
// 
// #include <brick/computerVision/imageFilter.hh>

// #include <cmath>
#include <brick/numeric/convolve2D.hh>
// #include <brick/numeric/functional.hh>

namespace brick {

  namespace computerVision {

#if 0

    /// @cond privateCode
    namespace privateCode {

      template<class OUTPUT_TYPE, class INPUT_TYPE, class KernelType>
      void
      filter2DZeroPad(brick::numeric::Array2D<OUTPUT_TYPE>& outputImage,
                      const brick::numeric::Array2D<KernelType>& kernel,
                      const brick::numeric::Array2D<INPUT_TYPE>& image)
      {
        if(outputImage.rows() != image.rows()
           || outputImage.columns() != image.columns()) {
          std::ostringstream message;
          message << "Size of output image (" << outputImage.rows()
                  << ", " << outputImage.columns()
                  << ") does not match size if input image ("
                  << image.rows() << ", " << image.columns() << ").";
          BRICK_THROW(ValueException, "filter2DZeroPad()",
                    message.str().c_str());
        }

        if((kernel.rows() % 2 != 1) || (kernel.columns() % 2 != 1)) {
          std::ostringstream message;
          message << "Kernel must have an odd number of rows and columns, "
                  << "but is (" <<  kernel.rows() << ", " << kernel.columns()
                  << ").";
          BRICK_THROW(ValueException, "filter2DZeroPad()",
                    message.str().c_str());
        }

        if((kernel.rows() >= image.rows())
           || (kernel.columns() >= image.columns())
           || (kernel.size() == 0)) {
          outputImage = static_cast<OUTPUT_TYPE>(0);
          return;
        }
        
        size_t halfKernelColumns = kernel.columns() / 2;
        size_t halfKernelRows = kernel.rows() / 2;
        size_t imageColumns = image.columns();
        size_t imageRows = image.rows();
        size_t kernelColumns = kernel.columns();
        size_t kernelRows = kernel.rows();
        size_t lastAffectedColumn = image.columns() - halfKernelColumns;
        size_t lastAffectedRow = image.rows() - halfKernelRows;

        // Get a pointer to the output image data.
        typename brick::numeric::Array2D<OUTPUT_TYPE>::iterator
          resultIterator = outputImage.begin();

        // Zero untouched rows at the top of the image.
        size_t rowIndex = 0;
        while(rowIndex < halfKernelRows) {
          for(size_t columnIndex = 0; columnIndex < image.columns();
              ++columnIndex) {
            *(resultIterator++) = static_cast<OUTPUT_TYPE>(0);
          }
          ++rowIndex;          
        }

        // Now process all of the rows which are affected by the
        // convolution.
        NumericTypeConversionFunctor<KernelType, OUTPUT_TYPE>
          numericTypeConversionFunctor;
        while(rowIndex < lastAffectedRow) {
          // Zero untouched pixels at the beginning of the row.
          size_t columnIndex = 0;
          while(columnIndex < halfKernelColumns) {
            *(resultIterator++) = static_cast<OUTPUT_TYPE>(0);
            ++columnIndex;
          }

          // Construct an iterator to keep track of where the
          // upper-left corner of the kernel falls in the input image.
          typename brick::numeric::Array2D<INPUT_TYPE>::const_iterator
            imageIterator = (image.begin()
                             + (rowIndex - halfKernelRows) * imageColumns);

          // Process most of the pixels in the current row.
          while(columnIndex < lastAffectedColumn) {
              
            // Calculate one dot-product between the kernel and the
            // current bit of image.
            KernelType dotProduct = static_cast<KernelType>(0);
            typename brick::numeric::Array2D<INPUT_TYPE>::const_iterator
              imageIteratorCopy = imageIterator;
            typename brick::numeric::Array2D<KernelType>::const_iterator
              kernelIterator = kernel.begin();
            for(size_t kernelRow = 0; kernelRow < kernelRows; ++kernelRow) {
              for(size_t kernelColumn = 0; kernelColumn < kernelColumns;
                  ++kernelColumn) {
                dotProduct += *(kernelIterator++) * *(imageIteratorCopy++);
              }
              imageIteratorCopy += imageColumns - kernelColumns;
            }
            *resultIterator = numericTypeConversionFunctor(dotProduct);

            // Advance to the next pixel.
            ++imageIterator;
            ++resultIterator;
            ++columnIndex;
          }

          // Zero untouched pixels at the end of the row.
          while(columnIndex < imageColumns) {
            *(resultIterator++) = static_cast<OUTPUT_TYPE>(0);
            ++columnIndex;
          }

          // Advance to the next row.
          ++rowIndex;
        }

        // Zero untouched rows at the bottom of the image.
        while(rowIndex < imageRows) {
          for(size_t columnIndex = 0; columnIndex < imageColumns;
              ++columnIndex) {
            *(resultIterator++) = static_cast<OUTPUT_TYPE>(0);
          }
          ++rowIndex;          
        }
      }


      
      template<class OUTPUT_TYPE, class INPUT_TYPE, class KernelType>
      void
      filterColumnsZeroPad(brick::numeric::Array2D<OUTPUT_TYPE>& outputImage,
                           const brick::numeric::Array1D<KernelType>& kernel,
                           const brick::numeric::Array2D<INPUT_TYPE>& image,
                           size_t skipColumns=0)
      {
        if(outputImage.rows() != image.rows()
           || outputImage.columns() != image.columns()) {
          std::ostringstream message;
          message << "Size of output image (" << outputImage.rows()
                  << ", " << outputImage.columns()
                  << ") does not match size if input image ("
                  << image.rows() << ", " << image.columns() << ").";
          BRICK_THROW(ValueException, "filterColumnsZeroPad()",
                    message.str().c_str());
        }

        if(kernel.size() % 2 != 1) {
          std::ostringstream message;
          message << "Kernel must have an odd number of elements, "
                  << "but has size of " <<  kernel.size() << ".";
          BRICK_THROW(ValueException, "filterColumnsZeroPad()",
                    message.str().c_str());
        }

        if((kernel.size() >= image.rows())
           || (kernel.size() == 0)) {
          outputImage = static_cast<OUTPUT_TYPE>(0);
          return;
        }
        
        size_t halfKernelSize = kernel.size() / 2;
        size_t imageColumns = image.columns();
        size_t imageRows = image.rows();
        size_t kernelSize = kernel.size();
        size_t lastAffectedColumn = imageColumns - skipColumns;
        size_t lastAffectedRow = imageRows - halfKernelSize;

        // Zero untouched rows at the top of the image.
        brick::numeric::Array1D<OUTPUT_TYPE> outputRow;
        size_t rowIndex = 0;
        while(rowIndex < halfKernelSize) {
          outputRow = outputImage.row(rowIndex);
          for(size_t columnIndex = skipColumns;
              columnIndex < lastAffectedColumn;
              ++columnIndex) {
            outputRow[columnIndex] = static_cast<OUTPUT_TYPE>(0);
          }
          ++rowIndex;          
        }

        // Now process all of the rows which are affected by the
        // convolution.
        brick::numeric::Array1D<INPUT_TYPE> inputRow;
        brick::numeric::Array1D<KernelType> accumulator(image.columns());
        while(rowIndex < lastAffectedRow) {
          outputRow = outputImage.row(rowIndex);
          inputRow = image.row(rowIndex - halfKernelSize);

	  size_t columnIndex0 = 0;
          while(columnIndex0 < skipColumns) {
	    outputRow(columnIndex0) = static_cast<OUTPUT_TYPE>(0);
	    ++columnIndex0;
	  }
	  columnIndex0 = lastAffectedColumn;
          while(columnIndex0 < outputImage.columns()) {
	    outputRow(columnIndex0) = static_cast<OUTPUT_TYPE>(0);
	    ++columnIndex0;
	  }
	    
          for(size_t columnIndex1 = skipColumns;
              columnIndex1 < lastAffectedColumn;
              ++columnIndex1) {
            accumulator(columnIndex1) = inputRow(columnIndex1) * kernel[0];
          }
          for(size_t kernelRow = 1; kernelRow < kernelSize; ++kernelRow) {
            inputRow = image.row(rowIndex - halfKernelSize + kernelRow);
            for(size_t columnIndex2 = skipColumns;
                columnIndex2 < lastAffectedColumn;
                ++columnIndex2) {
              accumulator(columnIndex2) +=
                inputRow(columnIndex2) * kernel[kernelRow];
            }
          }
          for(size_t columnIndex3 = skipColumns;
              columnIndex3 < lastAffectedColumn;
              ++columnIndex3) {
            outputRow(columnIndex3) =
              static_cast<OUTPUT_TYPE>(accumulator(columnIndex3) + 0.5);
          }

          // Advance to the next row.
          ++rowIndex;
        }

        // Zero untouched rows at the bottom of the image.
        while(rowIndex < imageRows) {
          outputRow = outputImage.row(rowIndex);
          for(size_t columnIndex = skipColumns;
              columnIndex < lastAffectedColumn;
              ++columnIndex) {
            outputRow[columnIndex] = static_cast<OUTPUT_TYPE>(0);
          }
          ++rowIndex;          
        }
      }


      template<class OUTPUT_TYPE, class INPUT_TYPE, class KernelType>
      void
      filterRowsZeroPad(brick::numeric::Array2D<OUTPUT_TYPE>& outputImage,
                        const brick::numeric::Array1D<KernelType>& kernel,
                        const brick::numeric::Array2D<INPUT_TYPE>& image,
                        size_t skipRows=0)
      {
        if(outputImage.rows() != image.rows()
           || outputImage.columns() != image.columns()) {
          std::ostringstream message;
          message << "Size of output image (" << outputImage.rows()
                  << ", " << outputImage.columns()
                  << ") does not match size if input image ("
                  << image.rows() << ", " << image.columns() << ").";
          BRICK_THROW(ValueException, "filterRowsZeroPad()",
                    message.str().c_str());
        }

        if(kernel.size() % 2 != 1) {
          std::ostringstream message;
          message << "Kernel must have an odd number of elements, "
                  << "but has size of " <<  kernel.size() << ".";
          BRICK_THROW(ValueException, "filterRowsZeroPad()",
                    message.str().c_str());
        }

        if((kernel.size() >= image.columns())
           || (kernel.size() == 0)) {
          outputImage = static_cast<OUTPUT_TYPE>(0);
          return;
        }
        
        size_t halfKernelSize = kernel.size() / 2;
        size_t imageColumns = image.columns();
        size_t imageRows = image.rows();
        size_t kernelSize = kernel.size();
        size_t lastAffectedColumn = imageColumns - halfKernelSize;
        size_t lastAffectedRow = imageRows - skipRows;

	// Zero out skipped rows.
	for(size_t index0 = 0; index0 < skipRows; ++index0) {
	  std::fill(outputImage.rowBegin(index0),
		    outputImage.rowBegin(index0),
		    static_cast<OUTPUT_TYPE>(0));
	}
	
        // Get a pointer to the output image data.
        typename brick::numeric::Array2D<OUTPUT_TYPE>::iterator
          resultIterator = outputImage.begin() + skipRows * imageColumns;

        // Now process all of the rows which are affected by the
        // convolution.
        size_t rowIndex = skipRows;
        NumericTypeConversionFunctor<KernelType, OUTPUT_TYPE>
          numericTypeConversionFunctor;
        while(rowIndex < lastAffectedRow) {
          // Zero untouched pixels at the beginning of the row.
          size_t columnIndex = 0;
          while(columnIndex < halfKernelSize) {
            *(resultIterator++) = static_cast<OUTPUT_TYPE>(0);
            ++columnIndex;
          }

          // Construct an iterator to keep track of where the
          // leftmost element of the kernel falls in the input image.
          typename brick::numeric::Array2D<INPUT_TYPE>::const_iterator
            imageIterator = image.begin() + rowIndex * imageColumns;

          // Process most of the pixels in the current row.
          while(columnIndex < lastAffectedColumn) {
            // Calculate one dot-product between the kernel and the
            // current bit of image.
            KernelType dotProduct = static_cast<KernelType>(0);
            typename brick::numeric::Array2D<INPUT_TYPE>::const_iterator
              imageIteratorCopy = imageIterator;
            typename brick::numeric::Array2D<KernelType>::const_iterator
              kernelIterator = kernel.begin();
            for(size_t kernelElement = 0; kernelElement < kernelSize;
                ++kernelElement) {
              // xxx problem here.
              dotProduct += *(kernelIterator++) * *(imageIteratorCopy++);
            }
            *resultIterator = numericTypeConversionFunctor(dotProduct);

            // Advance to the next pixel.
            ++imageIterator;
            ++resultIterator;
            ++columnIndex;
          }

          // Zero untouched pixels at the end of the row.
          while(columnIndex < imageColumns) {
            *(resultIterator++) = static_cast<OUTPUT_TYPE>(0);
            ++columnIndex;
          }

          // Advance to the next row.
          ++rowIndex;
        }

	// Zero out skipped rows.
	for(size_t index1 = lastAffectedRow; index1 < outputImage.rows();
	    ++index1) {
	  std::fill(outputImage.rowBegin(index1),
		    outputImage.rowBegin(index1),
		    static_cast<OUTPUT_TYPE>(0));
	}
      }

    } // namespace privateCode
    /// @endcond

#endif    

    // This function filters an image with the given kernel.
    template<ImageFormat OutputFormat,
	     ImageFormat IntermediateFormat,
             ImageFormat ImageFormat,
             class KernelType>
    Image<OutputFormat>
    filter2D(
      const Kernel<KernelType>& kernel,
      const Image<ImageFormat>& image,
      const typename ImageFormatTraits<OutputFormat>::PixelType fillValue,
      ConvolutionStrategy convolutionStrategy)
    {
      Image<OutputFormat> returnImage(image.rows(), image.columns());
      // filter2D<IntermediateFormat>(
      filter2D<OutputFormat, IntermediateFormat, ImageFormat, KernelType>(
	returnImage, kernel, image, fillValue, convolutionStrategy);
      return returnImage;
    }
    

    // This function filters an image with the given kernel, placing
    // the result into a pre-constructed Image instance.
    template<ImageFormat OutputFormat,
	     ImageFormat IntermediateFormat,
             ImageFormat ImageFormat,
             class KernelType>
    void
    filter2D(
      Image<OutputFormat>& outputImage,
      const Kernel<KernelType>& kernel,
      const Image<ImageFormat>& image,
      const typename ImageFormatTraits<OutputFormat>::PixelType fillValue,
      ConvolutionStrategy convolutionStrategy)
    {
      if(convolutionStrategy != brick::numeric::BRICK_CONVOLVE_PAD_RESULT) {
        BRICK_THROW(brick::common::NotImplementedException, "filter2D()",
                    "Currently, BRICK_CONVOLVE_PAD_RESULT is the only "
                    "supported convolution strategy.");
      }
      typedef typename ImageFormatTraits<OutputFormat>::PixelType
	OutputPixelType;
      typedef typename ImageFormatTraits<IntermediateFormat>::PixelType
	IntermediatePixelType;
      
      if(kernel.isSeparable()) {
	brick::numeric::Array1D<KernelType> kernelRowComponent =
	  kernel.getRowComponent();
	brick::numeric::Array1D<KernelType> kernelColumnComponent =
	  kernel.getColumnComponent();
	brick::numeric::Array2D<KernelType> rowKernel(
	  1, kernelRowComponent.size(), kernelRowComponent.data());
	brick::numeric::Array2D<KernelType> columnKernel(
	  kernelColumnComponent.size(), 1, kernelColumnComponent.data());
	
	Image<IntermediateFormat> intermediateImage =
	  numeric::correlate2D<IntermediatePixelType, IntermediatePixelType>(
	    rowKernel, image,
	    numeric::BRICK_CONVOLVE_PAD_RESULT,
	    numeric::BRICK_CONVOLVE_ROI_SAME,
	    fillValue);
	outputImage =
	  numeric::correlate2D<OutputPixelType, IntermediatePixelType>(
	    columnKernel, intermediateImage,
	    numeric::BRICK_CONVOLVE_PAD_RESULT,
	    numeric::BRICK_CONVOLVE_ROI_SAME,
	    fillValue);
      } else {
        outputImage =
	  numeric::correlate2D<OutputPixelType, IntermediatePixelType>(
	  kernel.getArray2D(), image,
	  numeric::BRICK_CONVOLVE_PAD_RESULT,
	  numeric::BRICK_CONVOLVE_ROI_SAME,
	  fillValue);
      }
    }

  } // namespace computerVision
  
} // namespace brick


#endif /* #ifndef BRICK_COMPUTERVISION_IMAGEFILTER_IMPL_HH */
