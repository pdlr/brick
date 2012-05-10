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
#include <brick/computerVision/utilities.hh>
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

    namespace privateCode {

      template <ImageFormat OutputFormat, ImageFormat InputFormat>
      void
      filterRows121(
        Image<OutputFormat>& outImage,
        Image<InputFormat> const& inImage,
        typename ImageFormatTraits<OutputFormat>::PixelType const fillValue,
        int const finalShift = -2);

      template <ImageFormat OutputFormat, ImageFormat InputFormat>
      void
      filterColumns121(
        Image<OutputFormat>& outImage,
        Image<InputFormat> const& inImage,
        typename ImageFormatTraits<OutputFormat>::PixelType const fillValue,
        int const finalShift = -2);
      
      template <ImageFormat OutputFormat, ImageFormat InputFormat>
      void
      filterRows14641(
        Image<OutputFormat>& outImage,
        Image<InputFormat> const& inImage,
        typename ImageFormatTraits<OutputFormat>::PixelType const fillValue,
        int const finalShift = -4);

      template <ImageFormat OutputFormat, ImageFormat InputFormat>
      void
      filterColumns14641(
        Image<OutputFormat>& outImage,
        Image<InputFormat> const& inImage,
        typename ImageFormatTraits<OutputFormat>::PixelType const fillValue,
        int const finalShift = -4);
      
      template <ImageFormat OutputFormat, ImageFormat InputFormat>
      void
      filterColumns1_6_15_20_15_6_1(
        Image<OutputFormat>& outImage,
        Image<InputFormat> const& inImage,
        typename ImageFormatTraits<OutputFormat>::PixelType const fillValue,
        int const finalShift = -6);

      template <ImageFormat OutputFormat, ImageFormat InputFormat>
      void
      filterRows1_6_15_20_15_6_1(
        Image<OutputFormat>& outImage,
        Image<InputFormat> const& inImage,
        typename ImageFormatTraits<OutputFormat>::PixelType const fillValue,
        int const finalShift = -6);
      
    } // namespace privateCode

    
    // This function filters an image with the given kernel.
    template<ImageFormat OutputFormat,
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
      filter2D<OutputFormat, ImageFormat, KernelType>(
	returnImage, kernel, image, fillValue, convolutionStrategy);
      return returnImage;
    }
    

    // This function filters an image with the given kernel, placing
    // the result into a pre-constructed Image instance.
    template<ImageFormat OutputFormat,
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
      
      if(kernel.isSeparable()) {
	brick::numeric::Array1D<KernelType> kernelRowComponent =
	  kernel.getRowComponent();
	brick::numeric::Array1D<KernelType> kernelColumnComponent =
	  kernel.getColumnComponent();
	brick::numeric::Array2D<KernelType> rowKernel(
	  1, kernelRowComponent.size(), kernelRowComponent.data());
	brick::numeric::Array2D<KernelType> columnKernel(
	  kernelColumnComponent.size(), 1, kernelColumnComponent.data());
	
	Image<OutputFormat> intermediateImage =
	  numeric::correlate2D<OutputPixelType, OutputPixelType>(
	    rowKernel, image,
	    numeric::BRICK_CONVOLVE_PAD_RESULT,
	    numeric::BRICK_CONVOLVE_ROI_SAME,
	    fillValue);
	outputImage =
	  numeric::correlate2D<OutputPixelType, OutputPixelType>(
	    columnKernel, intermediateImage,
	    numeric::BRICK_CONVOLVE_PAD_RESULT,
	    numeric::BRICK_CONVOLVE_ROI_SAME,
	    fillValue);
      } else {
        outputImage =
	  numeric::correlate2D<OutputPixelType, OutputPixelType>(
	  kernel.getArray2D(), image,
	  numeric::BRICK_CONVOLVE_PAD_RESULT,
	  numeric::BRICK_CONVOLVE_ROI_SAME,
	  fillValue);
      }
    }


    // This function low-pass filters each column of an image with
    // integer-valued pixels using a binomial approximation to a
    // Gaussian kernel.
    template<ImageFormat OutputFormat, ImageFormat ImageFormat>
    void
    filterColumnsBinomial(
      Image<OutputFormat>& outputImage,
      const Image<ImageFormat>& inputImage,
      double sigma,
      typename ImageFormatTraits<OutputFormat>::PixelType const fillValue,
      int finalShift)
    {
      // Make sure outputImage is appropriately sized.
      if(outputImage.rows() != inputImage.rows()
         || outputImage.columns() != inputImage.columns()) {
        BRICK_THROW(brick::common::ValueException, "filterColumnsBinomial()",
                    "Arguments outputImage and inputImage must have "
                    "the same shape.");
      }

      // Make sure filter sizes are supported.
      if(sigma < 0.0 || sigma > 1.25) {
        BRICK_THROW(brick::common::NotImplementedException,
                    "filterColumnsBinomial()",
                    "Argument sigma must be between 0.0 and 1.25, inclusive");
      }

      // We will dispatch to custom routines for each kernel size.
      // First figure out which routine to call for each pass.  The
      // three element filter has variance 0.5, the 5 element filter
      // has variance 1.0, and the 7 element filter has variance 1.5.
      // This works out to filter size == 4 * variance + 1.
      unsigned int filterSize = static_cast<unsigned int>(
        4.0 * sigma * sigma + 1.5);

      // Now do the actual filtering.
      switch(filterSize) {
      case 0:
        outputImage.copy(inputImage);
        break;
      case 1:
      case 2:
      case 3:
        if(-1 == finalShift) {finalShift = 2;}
        privateCode::filterColumns121(outputImage, inputImage,
                                      fillValue, finalShift);
        break;
      case 4:
      case 5:
        if(-1 == finalShift) {finalShift = 4;}
        privateCode::filterColumns14641(outputImage, inputImage,
                                        fillValue, finalShift);
        break;
      case 6:
      case 7:
        if(-1 == finalShift) {finalShift = 6;}
        privateCode::filterColumns1_6_15_20_15_6_1(outputImage, inputImage,
                                                   fillValue, finalShift);
        break;
      default:
        BRICK_THROW(brick::common::LogicException, "filterClumnsBinomial()",
                    "Math error in computing filter kernel size.");
      }
    }


    // This function low-pass filters each row of an image with
    // integer-valued pixels using a binomial approximation to a
    // Gaussian kernel.
    template<ImageFormat OutputFormat, ImageFormat ImageFormat>
    void
    filterRowsBinomial(
      Image<OutputFormat>& outputImage,
      const Image<ImageFormat>& inputImage,
      double sigma,
      typename ImageFormatTraits<OutputFormat>::PixelType const fillValue,
      int finalShift)
    {
      // Make sure outputImage is appropriately sized.
      if(outputImage.rows() != inputImage.rows()
         || outputImage.columns() != inputImage.columns()) {
        BRICK_THROW(brick::common::ValueException, "filterRowsBinomial()",
                    "Arguments outputImage and inputImage must have "
                    "the same shape.");
      }

      // Make sure filter sizes are supported.
      if(sigma < 0.0 || sigma > 1.25) {
        BRICK_THROW(brick::common::NotImplementedException,
                    "filterColumnsBinomial()",
                    "Argument sigma must be between 0.0 and 1.25, inclusive");
      }

      // We will dispatch to custom routines for each kernel size.
      // First figure out which routine to call for each pass.  The
      // three element filter has variance 0.5, the 5 element filter
      // has variance 1.0, and the 7 element filter has variance 1.5.
      // This works out to filter size == 4 * variance + 1.
      unsigned int filterSize = static_cast<unsigned int>(
        4.0 * sigma * sigma + 1.5);

      // Now do the actual filtering.
      switch(filterSize) {
      case 0:
        outputImage.copy(inputImage);
        break;
      case 1:
      case 2:
      case 3:
        if(-1 == finalShift) {finalShift = 2;}
        privateCode::filterRows121(outputImage, inputImage,
                                   fillValue, finalShift);
        break;
      case 4:
      case 5:
        if(-1 == finalShift) {finalShift = 4;}
        privateCode::filterRows14641(outputImage, inputImage,
                                     fillValue, finalShift);
        break;
      case 6:
      case 7:
        if(-1 == finalShift) {finalShift = 6;}
        privateCode::filterRows1_6_15_20_15_6_1(outputImage, inputImage,
                                                fillValue, finalShift);
        break;
      default:
        BRICK_THROW(brick::common::LogicException, "filterClumnsBinomial()",
                    "Math error in computing filter kernel size.");
      }
    }

  } // namespace computerVision

} // namespace brick
    

// Definitions of local symbols below.

namespace brick {

  namespace computerVision {

    namespace privateCode {

      template <ImageFormat OutputFormat, ImageFormat InputFormat>
      void
      filterColumns121(
        Image<OutputFormat>& outImage,
        Image<InputFormat> const& inImage,
        typename ImageFormatTraits<OutputFormat>::PixelType const fillValue,
        int const finalShift)
      {
        typedef typename ImageFormatTraits<InputFormat>::PixelType
          InElementType;
        typedef typename ImageFormatTraits<OutputFormat>::PixelType
          OutElementType;

        // Make sure outImage is appropriately sized.
        if(outImage.rows() != inImage.rows()
           || outImage.columns() != inImage.columns()) {
          BRICK_THROW(brick::common::ValueException, "filterColumns121()",
                      "Arguments outImage and inImage must have "
                      "the same shape.");
        }

        // Get all our dimensions straight.
        unsigned int const stopRow = inImage.rows() - 1;
        unsigned int const numColumns = inImage.columns();
        unsigned int const inRowStep = inImage.getRowStep();

        // Fill pixels for which we won't have valid results.
        OutElementType* outPtr0 = const_cast<OutElementType*>(
          outImage.data(0, 0));
        for(unsigned int index = 0; index < numColumns; ++index) {
          outPtr0[index] = fillValue;
        }
        
        // Now convolve, but do it row-wise to minimize cache misses.
        for(unsigned int row = 1; row < stopRow; ++row) {
          InElementType* inPtr = const_cast<InElementType*>(
            inImage.data(row, 0));
          OutElementType* outPtr = &(outImage(row, 0));
          for(unsigned int index = 0; index < numColumns; ++index) {
            outPtr[index] = static_cast<OutElementType>(
              (inPtr[index - inRowStep]
               + (inPtr[index] << 1)
               + inPtr[index + inRowStep]) >> finalShift);
          }
        }

        // Fill pixels for which we won't have valid results.
        outPtr0 = const_cast<OutElementType*>(outImage.data(stopRow, 0));
        for(unsigned int index = 0; index < numColumns; ++index) {
          outPtr0[index] = fillValue;
        }
      }


      template <ImageFormat OutputFormat, ImageFormat InputFormat>
      void
      filterColumns14641(
        Image<OutputFormat>& outImage,
        Image<InputFormat> const& inImage,
        typename ImageFormatTraits<OutputFormat>::PixelType const fillValue,
        int const finalShift)
      {
        typedef typename ImageFormatTraits<InputFormat>::PixelType
          InElementType;
        typedef typename ImageFormatTraits<OutputFormat>::PixelType
          OutElementType;

        // Make sure outImage is appropriately sized.
        if(outImage.rows() != inImage.rows()
           || outImage.columns() != inImage.columns()) {
          BRICK_THROW(brick::common::ValueException, "filterColumns121()",
                      "Arguments outImage and inImage must have "
                      "the same shape.");
        }

        // Get all our dimensions straight.
        unsigned int const stopRow = inImage.rows() - 2;
        unsigned int const numColumns = inImage.columns();
        unsigned int const inRowStep = inImage.getRowStep();
        unsigned int const twoInRowStep = inRowStep * 2;
        unsigned int const outRowStep = outImage.getRowStep();

        // Fill pixels for which we won't have valid results.
        OutElementType* outPtr0 = const_cast<OutElementType*>(
          outImage.data(0, 0));
        for(unsigned int index = 0; index < numColumns; ++index) {
          outPtr0[index] = fillValue;
          outPtr0[index + outRowStep] = fillValue;
        }
        
        // Now convolve, but do it row-wise to minimize cache misses.
        for(unsigned int row = 2; row < stopRow; ++row) {
          InElementType* inPtr = const_cast<InElementType*>(
            inImage.data(row, 0));
          OutElementType* outPtr = &(outImage(row, 0));
          for(unsigned int index = 0; index < numColumns; ++index) {
            outPtr[index] = static_cast<OutElementType>(
              (inPtr[index - twoInRowStep]
               + (inPtr[index - inRowStep] << 2)
               + (inPtr[index] * 6)
               + (inPtr[index + inRowStep] << 2)
               + inPtr[index + twoInRowStep]) >> finalShift);
          }
        }

        // Fill pixels for which we won't have valid results.
        outPtr0 = const_cast<OutElementType*>(outImage.data(stopRow, 0));
        for(unsigned int index = 0; index < numColumns; ++index) {
          outPtr0[index] = fillValue;
          outPtr0[index + outRowStep] = fillValue;
        }
      }

      
      template <ImageFormat OutputFormat, ImageFormat InputFormat>
      void
      filterColumns1_6_15_20_15_6_1(
        Image<OutputFormat>& outImage,
        Image<InputFormat> const& inImage,
        typename ImageFormatTraits<OutputFormat>::PixelType const fillValue,
        int const finalShift)
      {
        typedef typename ImageFormatTraits<InputFormat>::PixelType
          InElementType;
        typedef typename ImageFormatTraits<OutputFormat>::PixelType
          OutElementType;

        // Make sure outImage is appropriately sized.
        if(outImage.rows() != inImage.rows()
           || outImage.columns() != inImage.columns()) {
          BRICK_THROW(brick::common::ValueException, "filterColumns121()",
                      "Arguments outImage and inImage must have "
                      "the same shape.");
        }

        // Get all our dimensions straight.
        unsigned int const stopRow = inImage.rows() - 3;
        unsigned int const numColumns = inImage.columns();
        unsigned int const inRowStep = inImage.getRowStep();
        unsigned int const twoInRowStep = inRowStep * 2;
        unsigned int const threeInRowStep = inRowStep * 3;
        unsigned int const outRowStep = outImage.getRowStep();
        unsigned int const twoOutRowStep = outRowStep * 2;
        
        // Fill pixels for which we won't have valid results.
        OutElementType* outPtr0 = const_cast<OutElementType*>(
          outImage.data(0, 0));
        for(unsigned int index = 0; index < numColumns; ++index) {
          outPtr0[index] = fillValue;
          outPtr0[index + outRowStep] = fillValue;
          outPtr0[index + twoOutRowStep] = fillValue;
        }
        
        // Now convolve, but do it row-wise to minimize cache misses.
        for(unsigned int row = 3; row < stopRow; ++row) {
          InElementType* inPtr = const_cast<InElementType*>(
            inImage.data(row, 0));
          OutElementType* outPtr = &(outImage(row, 0));
          for(unsigned int index = 0; index < numColumns; ++index) {
            outPtr[index] = static_cast<OutElementType>((
              inPtr[index - threeInRowStep] + inPtr[index + threeInRowStep]
              + (inPtr[index - twoInRowStep] + inPtr[index + twoInRowStep]) * 6
              + (inPtr[index - inRowStep] + inPtr[index + inRowStep]) * 15
              + (inPtr[index] * 20)) >> finalShift);
          }
        }

        // Fill pixels for which we won't have valid results.
        outPtr0 = const_cast<OutElementType*>(outImage.data(stopRow, 0));
        for(unsigned int index = 0; index < numColumns; ++index) {
          outPtr0[index] = fillValue;
          outPtr0[index + outRowStep] = fillValue;
          outPtr0[index + twoOutRowStep] = fillValue;
        }
      }


      template <ImageFormat OutputFormat, ImageFormat InputFormat>
      void
      filterRows121(
        Image<OutputFormat>& outImage,
        Image<InputFormat> const& inImage,
        typename ImageFormatTraits<OutputFormat>::PixelType const fillValue,
        int const finalShift)
      {
        typedef typename ImageFormatTraits<InputFormat>::PixelType
          InElementType;
        typedef typename ImageFormatTraits<OutputFormat>::PixelType
          OutElementType;

        // Make sure outImage is appropriately sized.
        if(outImage.rows() != inImage.rows()
           || outImage.columns() != inImage.columns()) {
          BRICK_THROW(brick::common::ValueException, "filterRows121()",
                      "Arguments outImage and inImage must have "
                      "the same shape.");
        }
        
        // Now convolve each row in turn.
        unsigned int stopColumn = inImage.columns() - 1;
        for(unsigned int row = 0; row < inImage.rows(); ++row) {
          InElementType* inPtr = const_cast<InElementType*>(
            inImage.data(row, 0));
          OutElementType* outPtr = outImage.data(row, 0);
          outPtr[0] = fillValue;
          for(unsigned int column = 1; column < stopColumn; ++column) {
            outPtr[column] = static_cast<OutElementType>(
              (inPtr[column - 1]
               + (inPtr[column] << 1)
               + inPtr[column + 1]) >> finalShift);
          }
          outPtr[stopColumn] = fillValue;
        }
      }

      
      template <ImageFormat OutputFormat, ImageFormat InputFormat>
      void
      filterRows14641(
        Image<OutputFormat>& outImage,
        Image<InputFormat> const& inImage,
        typename ImageFormatTraits<OutputFormat>::PixelType const fillValue,
        int const finalShift)
      {
        typedef typename ImageFormatTraits<InputFormat>::PixelType
          InElementType;
        typedef typename ImageFormatTraits<OutputFormat>::PixelType
          OutElementType;

        // Make sure outImage is appropriately sized.
        if(outImage.rows() != inImage.rows()
           || outImage.columns() != inImage.columns()) {
          BRICK_THROW(brick::common::ValueException, "filterRows14641()",
                      "Arguments outImage and inImage must have "
                      "the same shape.");
        }

        // Now convolve each row in turn.
        unsigned int stopColumn = inImage.columns() - 2;
        for(unsigned int row = 0; row < inImage.rows(); ++row) {
          InElementType* inPtr = const_cast<InElementType*>(
            inImage.data(row, 0));
          OutElementType* outPtr = outImage.data(row, 0);

          outPtr[0] = fillValue;
          outPtr[1] = fillValue;
          for(unsigned int column = 2; column < stopColumn; ++column) {
            outPtr[column] = static_cast<OutElementType>(
              (inPtr[column - 2]
               + (inPtr[column - 1] << 2)
               + (inPtr[column] * 6)
               + (inPtr[column + 1] << 2)
               + inPtr[column + 2]) >> finalShift);
          }
          outPtr[stopColumn] = fillValue;
          outPtr[stopColumn + 11] = fillValue;
        }
      }


      template <ImageFormat OutputFormat, ImageFormat InputFormat>
      void
      filterRows1_6_15_20_15_6_1(
        Image<OutputFormat>& outImage,
        Image<InputFormat> const& inImage,
        typename ImageFormatTraits<OutputFormat>::PixelType const fillValue,
        int const finalShift)
      {
        typedef typename ImageFormatTraits<InputFormat>::PixelType
          InElementType;
        typedef typename ImageFormatTraits<OutputFormat>::PixelType
          OutElementType;

        // Make sure outImage is appropriately sized.
        if(outImage.rows() != inImage.rows()
           || outImage.columns() != inImage.columns()) {
          BRICK_THROW(brick::common::ValueException, "filterRows14641()",
                      "Arguments outImage and inImage must have "
                      "the same shape.");
        }

        // Now convolve each row in turn.
        unsigned int stopColumn = inImage.columns() - 3;
        for(unsigned int row = 0; row < inImage.rows(); ++row) {
          InElementType* inPtr = const_cast<InElementType*>(
            inImage.data(row, 0));
          OutElementType* outPtr = outImage.data(row, 0);

          outPtr[0] = fillValue;
          outPtr[1] = fillValue;
          outPtr[2] = fillValue;
          for(unsigned int column = 2; column < stopColumn; ++column) {
            outPtr[column] = static_cast<OutElementType>(
              (inPtr[column - 3] + inPtr[column + 3]
               + ((inPtr[column - 2] + inPtr[column + 2]) * 6)
               + ((inPtr[column - 1] + inPtr[column + 1]) * 15)
               + (inPtr[column] * 20)) >> finalShift);
          }
          outPtr[stopColumn] = fillValue;
          outPtr[stopColumn + 1] = fillValue;
          outPtr[stopColumn + 2] = fillValue;
        }

      }
      
    } // namespace privateCode
    
  } // namespace computerVision
  
} // namespace brick


// Test routines we're keeping around in case we want to borrow from
// them later.
namespace {

  template <brick::computerVision::ImageFormat Format>
  brick::computerVision::Image<Format> 
  filterColumns121_simple(
    brick::computerVision::Image<Format> const& inImage,
    typename brick::computerVision::ImageFormatTraits<Format>::PixelType const /* fillValue */,
    int finalShift)
  {
    typedef typename brick::computerVision::ImageFormatTraits<Format>::PixelType ElementType;

    // Allocate space for output.
    brick::computerVision::Image<Format> resultImage(inImage.rows(), inImage.columns());

    // Now convolve each column in turn.
    if(finalShift < 0) {finalShift = 2;}
    for(unsigned int column = 0; column < inImage.columns(); ++column) {
      ElementType* inColumnStart = const_cast<ElementType*>(
        inImage.data(0, column));
      ElementType* outColumnStart = &(resultImage(0, column));
      ElementType* inPtr = inColumnStart;
      ElementType* outPtr = outColumnStart;
      const int numColumns = inImage.columns();
      unsigned int stopIndex = (inImage.rows() - 1) * numColumns;
      for(unsigned int index = numColumns; index < stopIndex;
          index += numColumns) {
        outPtr[index] = (inPtr[index - numColumns]
                         + (inPtr[index] << 1)
                         + inPtr[index + numColumns]) >> finalShift;
      }
    }
    return resultImage;
  }


  template <brick::computerVision::ImageFormat Format>
  brick::computerVision::Image<Format> 
  filterRows121_multipass(
    brick::computerVision::Image<Format> const& inImage,
    unsigned int numberOfPasses,
    typename brick::computerVision::ImageFormatTraits<Format>::PixelType const fillValue,
    int finalShift)
  {
    typedef typename brick::computerVision::ImageFormatTraits<Format>::PixelType ElementType;

    // Temporary storage for multi-pass filtering.  This is a
    // wasted allocation if we are only doing a single pass.
    brick::numeric::Array1D<typename brick::computerVision::ImageFormatTraits<Format>::PixelType>
      intermediateRow(inImage.columns());

    // Allocate space for output.
    brick::computerVision::Image<Format> resultImage(inImage.rows(), inImage.columns());

    // Now convolve each row in turn.
    if(finalShift < 0) {
      finalShift = 2 * numberOfPasses;
    }
    for(unsigned int row = 0; row < inImage.rows(); ++row) {

      // Build up by repeatedly convolving with [1, 2, 1].
      ElementType* inRowStart = const_cast<ElementType*>(
        inImage.data(row, 0));
      ElementType* outRowStart = &(resultImage(row, 0));
      for(unsigned int pass = 0; pass < numberOfPasses; ++pass) {
        ElementType* inPtr = inRowStart;
        ElementType* outPtr = outRowStart;
        unsigned int stopColumn = inImage.columns() - (pass + 1);
        for(unsigned int column = pass + 1; column < stopColumn; ++column) {
          outPtr[column] =
            inPtr[column - 1] + (inPtr[column] << 1) + inPtr[column + 1];
        }

        // Ping pong back and forth between the two images on
        // alternate passes.
        if(pass == 0) {
          inRowStart = outRowStart;
          outRowStart = &(intermediateRow[0]);
        } else {
          std::swap(inRowStart, outRowStart);
        }
      }

      // xxx Wrong!
      if(finalShift > 0) {
        ElementType* finalRow = &(resultImage(row, 0));
        unsigned int stopColumn = inImage.columns() - numberOfPasses;
        for(unsigned int column = 0; column < numberOfPasses; ++column) {
          finalRow[column] = fillValue;
        }
        for(unsigned int column = numberOfPasses; column < stopColumn;
            ++column) {
          finalRow[column] = inRowStart[column] >> finalShift;
        }
        for(unsigned int column = stopColumn; column < inImage.columns();
            ++column) {
          finalRow[column] = fillValue;
        }
      }
    }
    return resultImage;
  }


  template <brick::computerVision::ImageFormat Format>
  brick::computerVision::Image<Format> 
  filter_1_6_15_20_15_6_1(brick::computerVision::Image<Format> const& inImage,
                          const unsigned int stepSize,
                          const unsigned int startRow,
                          const unsigned int stopRow,
                          const unsigned int startColumn,
                          const unsigned int stopColumn)
  {
    typedef typename brick::computerVision::ImageFormatTraits<Format>::PixelType ElementType;

    // Allocate space for output.
    brick::computerVision::Image<Format> resultImage(inImage.rows(), inImage.columns());

    // Now convolve each row in turn.
    const unsigned int stepSizeX2 = 2 * stepSize;
    const unsigned int stepSizeX3 = 3 * stepSize;
    for(unsigned int row = startRow; row < stopRow; ++row) {
      ElementType const* inPtr = inImage.getData(row, 0);
      ElementType* outPtr = resultImage.getData(row, 0);
      for(unsigned int column = startColumn; column < stopColumn;
          ++column) {
        outPtr[column] = (
          inPtr[column - stepSizeX3] + inPtr[column + stepSizeX3]
          + (inPtr[column - stepSizeX2] + inPtr[column + stepSizeX2]) * 6
          + (inPtr[column - stepSize] + inPtr[column + stepSize]) * 15
          + (inPtr[column] * 20)) >> 6;
      }
    }
    return resultImage;
  }
  
} // namespace

#endif /* #ifndef BRICK_COMPUTERVISION_IMAGEFILTER_IMPL_HH */
