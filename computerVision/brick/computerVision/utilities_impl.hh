/**
***************************************************************************
* @file brick/computerVision/utilities_impl.hh
*
* Header file defining utility functions for brick::computerVision
* library.
*
* Copyright (C) 2006, 2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_UTILITIES_IMPL_HH
#define BRICK_COMPUTERVISION_UTILITIES_IMPL_HH

// This file is included by utilities.hh, and should not be directly included
// by user code, so no need to include utilities.hh here.
//
// #include <brick/computerVision/utilities.hh>

#include <brick/linearAlgebra/linearAlgebra.hh>

namespace brick {

  namespace computerVision {

    /// @cond privateCode
    namespace privateCode {

      // This function simplifies template specialization of
      // associateColorComponents, which makes things easier on the
      // compiler.
      template<ImageFormat Format, class ComponentType>
      inline bool
      associateColorComponents(
	brick::numeric::Array2D<ComponentType>& inputArray,
	Image<Format>& outputImage)
      {
	// Some typedefs for notational convenience.
	typedef typename ImageFormatTraits<Format>::PixelType PixelType;
	typedef typename Image<Format>::iterator ImageIterator;
	typedef typename brick::numeric::Array2D<ComponentType>::iterator
          ArrayIterator;

	// Argument checking.
	if(inputArray.columns()
	   % ImageFormatTraits<Format>::getNumberOfComponents()
	   != 0) {
	  BRICK_THROW(brick::common::ValueException,
                      "associateColorComponents()",
                      "Input image must have a number of columns which is "
                      "evenly divisible by the number of color components.");
	}

	// Calculate output image width.
	size_t numberOfOutputColumns = inputArray.columns()
	  / ImageFormatTraits<Format>::getNumberOfComponents();

	bool returnValue;
	if(PixelType::isContiguous()) {
	  // Build the output image.
	  outputImage = Image<Format>(
	    inputArray.rows(), numberOfOutputColumns,
	    reinterpret_cast<PixelType*>(inputArray.data()));
	  returnValue = true;
	} else {
	  if(outputImage.rows() != inputArray.rows()
	     || outputImage.columns() != numberOfOutputColumns) {
	    outputImage = Image<Format>(
	      inputArray.rows(), numberOfOutputColumns);
	  }
	  ArrayIterator arrayIterator = inputArray.begin();
	  for(ImageIterator imageIterator = outputImage.begin();
	      imageIterator != outputImage.end(); ++imageIterator) {
	    imageIterator->copyFromIterator(arrayIterator);
	    ++imageIterator;
	  }
	  returnValue = false;
	}
	return returnValue;
      }


      // Explicit specializations of
      // privateCode::associateColorComponents() are needed for those
      // image formats which use builtin types as pixels.
      template<>
      inline bool
      associateColorComponents(
	brick::numeric::Array2D<bool>& inputArray,
	Image<GRAY1>& outputImage)
      {
	outputImage = inputArray;
	return true;
      }


      template<>
      inline bool
      associateColorComponents(
        brick::numeric::Array2D<brick::common::UnsignedInt8>& inputArray,
        Image<GRAY8>& outputImage)
      {
	outputImage = inputArray;
	return true;
      }


      template<>
      inline bool
      associateColorComponents(
        brick::numeric::Array2D<brick::common::UnsignedInt16>& inputArray,
        Image<GRAY16>& outputImage)
      {
	outputImage = inputArray;
	return true;
      }


      template<>
      inline bool
      associateColorComponents(
        brick::numeric::Array2D<brick::common::UnsignedInt32>& inputArray,
        Image<GRAY32>& outputImage)
      {
	outputImage = inputArray;
	return true;
      }


      template<>
      inline bool
      associateColorComponents(
        brick::numeric::Array2D<brick::common::Int16>& inputArray,
        Image<GRAY_SIGNED16>& outputImage)
      {
	outputImage = inputArray;
	return true;
      }


      template<>
      inline bool
      associateColorComponents(
        brick::numeric::Array2D<brick::common::Int32>& inputArray,
        Image<GRAY_SIGNED32>& outputImage)
      {
	outputImage = inputArray;
	return true;
      }


      template<>
      inline bool
      associateColorComponents(
        brick::numeric::Array2D<brick::common::Float32>& inputArray,
        Image<GRAY_FLOAT32>& outputImage)
      {
	outputImage = inputArray;
	return true;
      }


      template<>
      inline bool
      associateColorComponents(
        brick::numeric::Array2D<brick::common::Float64>& inputArray,
        Image<GRAY_FLOAT64>& outputImage)
      {
	outputImage = inputArray;
	return true;
      }

      // This function returns by reference an array which either shares
      // or copies the data from the input image.
      template<ImageFormat Format, class ComponentType>
      inline bool
      dissociateColorComponents(
	Image<Format>& inputImage,
	brick::numeric::Array2D<ComponentType>& outputArray)
      {
	// Some typedefs for notational convenience.
	typedef typename ImageFormatTraits<Format>::PixelType PixelType;
	typedef typename Image<Format>::iterator ImageIterator;
	typedef typename brick::numeric::Array2D<ComponentType>::iterator ArrayIterator;

	// Calculate output array width.
	size_t numberOfOutputColumns = inputImage.columns()
	  * ImageFormatTraits<Format>::getNumberOfComponents();

	// Now build the output array.
	bool returnValue;
	if(PixelType::isContiguous()) {
	  outputArray = brick::numeric::Array2D<ComponentType>(
	    inputImage.rows(), numberOfOutputColumns,
	    reinterpret_cast<ComponentType*>(inputImage.data()));
	  returnValue = true;
	} else {
	  if(outputArray.rows() != inputImage.rows()
	     || outputArray.columns() != numberOfOutputColumns) {
	    outputArray = brick::numeric::Array2D<ComponentType>(
	      inputImage.rows(), numberOfOutputColumns);
	  }
	  ArrayIterator arrayIterator = outputArray.begin();
	  for(ImageIterator imageIterator = inputImage.begin();
	      imageIterator != inputImage.end(); ++imageIterator) {
	    imageIterator->copyToIterator(arrayIterator);
	    ++imageIterator;
	  }
	  returnValue = false;
	}
	return returnValue;
      }


      // Explicit specializations of dissociateColorComponents() are needed
      // for those image formats which use builtin types as pixels.
      template<>
      inline bool
      dissociateColorComponents(
        Image<GRAY1>& inputImage,
        brick::numeric::Array2D<bool>& outputArray)
      {
	outputArray = inputImage;
	return true;
      }


      template<>
      inline bool
      dissociateColorComponents(
        Image<GRAY8>& inputImage,
        brick::numeric::Array2D<brick::common::UnsignedInt8>& outputArray)
      {
	outputArray = inputImage;
	return true;
      }


      template<>
      inline bool
      dissociateColorComponents(
        Image<GRAY16>& inputImage,
        brick::numeric::Array2D<brick::common::UnsignedInt16>& outputArray)
      {
	outputArray = inputImage;
	return true;
      }


      template<>
      inline bool
      dissociateColorComponents(
        Image<GRAY32>& inputImage,
        brick::numeric::Array2D<brick::common::UnsignedInt32>& outputArray)
      {
	outputArray = inputImage;
	return true;
      }


      template<>
      inline bool
      dissociateColorComponents(
        Image<GRAY_SIGNED16>& inputImage,
        brick::numeric::Array2D<brick::common::Int16>& outputArray)
      {
	outputArray = inputImage;
	return true;
      }


      template<>
      inline bool
      dissociateColorComponents(
        Image<GRAY_SIGNED32>& inputImage,
        brick::numeric::Array2D<brick::common::Int32>& outputArray)
      {
	outputArray = inputImage;
	return true;
      }


      template<>
      inline bool
      dissociateColorComponents(
        Image<GRAY_FLOAT32>& inputImage,
        brick::numeric::Array2D<brick::common::Float32>& outputArray)
      {
	outputArray = inputImage;
	return true;
      }


      template<>
      inline bool
      dissociateColorComponents(
        Image<GRAY_FLOAT64>& inputImage,
        brick::numeric::Array2D<brick::common::Float64>& outputArray)
      {
	outputArray = inputImage;
	return true;
      }

    } // namespace privateCode
    /// @endcond


    // This function returns by reference an image which either shares
    // or copies the data from the input array.
    template<ImageFormat FORMAT>
    bool
    associateColorComponents(
      brick::numeric::Array2D<typename ImageFormatTraits<FORMAT>
      ::ComponentType>& inputArray,
      Image<FORMAT>& outputImage)
    {
      return privateCode::associateColorComponents(inputArray, outputImage);
    }


    // This deprecated function tries to return an Image which
    // references the same memory as the input array, but in which each
    // pixel is the aggregate of the appropriate number of elements from
    // the array.  If it is not possible to do so, then this function
    template<ImageFormat FORMAT>
    Image<FORMAT>
    associateColorComponents(
      brick::numeric::Array2D<typename ImageFormatTraits<FORMAT>
      ::ComponentType>& inputArray)
    {
      typedef typename ImageFormatTraits<FORMAT>::PixelType PixelType;
      if(!PixelType::isContiguous()) {
        BRICK_THROW(brick::common::LogicException, "associateColorComponents()",
                  "The deprecated single-argument version of "
                  "associateColorComponents() can only be use with "
                  "contiguous pixel types.");
      }

      Image<FORMAT> outputImage;
      associateColorComponents(inputArray, outputImage);
      return outputImage;
    }


    // This function takes an image in one colorspace and generates a
    // corresponding image in a second colorspace.
    template<ImageFormat OUTPUT_FORMAT, ImageFormat INPUT_FORMAT>
    Image<OUTPUT_FORMAT>
    convertColorspace(const Image<INPUT_FORMAT>& inputImage)
    {
      Image<OUTPUT_FORMAT> outputImage(
	inputImage.rows(), inputImage.columns());
      ColorspaceConverter<INPUT_FORMAT, OUTPUT_FORMAT> converter;
      std::transform(inputImage.begin(), inputImage.end(),
                     outputImage.begin(), converter);
      return outputImage;
    }


    // This function returns by reference an array which either shares
    // or copies the data from the input image.
    template<ImageFormat FORMAT>
    bool
    dissociateColorComponents(
      Image<FORMAT>& inputImage,
      brick::numeric::Array2D<typename ImageFormatTraits<FORMAT>
      ::ComponentType>& outputArray)
    {
      return privateCode::dissociateColorComponents(inputImage, outputArray);
    }


    // This deprecated function tries to return an Array2D which
    // references the same memory as the input image, but in which each
    // pixel has been "flattened" so that the returned array has a
    // separate element for each color component of each pixel.  If it
    template<ImageFormat FORMAT>
    brick::numeric::Array2D<typename ImageFormatTraits<FORMAT>::ComponentType>
    dissociateColorComponents(Image<FORMAT>& inputImage)
    {
      typedef typename ImageFormatTraits<FORMAT>::PixelType PixelType;
      if(!PixelType::isContiguous()) {
        BRICK_THROW(brick::common::LogicException, "dissociateColorComponents()",
                  "The deprecated single-argument version of "
                  "dissociateColorComponents() can only be use with "
                  "contiguous pixel types.");
      }
      brick::numeric::Array2D<typename ImageFormatTraits<FORMAT>::ComponentType> outputArray;
      dissociateColorComponents(inputImage, outputArray);
      return outputArray;
    }


    // Return the transform that, in the least-squares sense, most
    // nearly transforms each point in fromPoints to the corresponding
    // point in toPoints.
    template <class Type, class Iter0, class Iter1>
    brick::numeric::Transform2D<Type>
    estimateAffineTransform(Iter0 toPointsBegin, Iter0 toPointsEnd,
                            Iter1 fromPointsBegin)
    {
      // Affine transforms can be expressed in matrix form as:
      //
      //       | a, b, c |
      //   T = | c, d, e |
      //       | 0, 0, 1 |
      //
      // So that (using homogeneous coordinates), the point (x, y)
      // gets transformed to (u, v) like this:
      //
      //   | u |       | x |
      //   | v | = T * | y |
      //   | 1 |       | 1 |
      //
      // Rearranging this, we have:
      //
      //   | x, y, 1, 0, 0, 0 |   | a |   | u |
      //   | 0, 0, 0, x, y, 1 | * | b | = | v |
      //                          | c |
      //                          | d |
      //                          | e |
      //                          | f |
      //
      // Given more than one (x, y) to (u, v) correspondence, we
      // combine all of the simultaneous equations into one big linear
      // regression, and solve for the unknowns a through f.
      //
      //       | a |
      //   A * | b | = B
      //       | c |
      //       | d |
      //       | e |
      //       | f |
      //
      unsigned int numberOfPoints = toPointsEnd - toPointsBegin;
      if(numberOfPoints < 3) {
        BRICK_THROW(brick::common::ValueException, "estimateAffineTransform()",
                    "Input sequences must have at least three elements.");
      }

      numeric::Array2D<brick::common::Float64> AMatrix(numberOfPoints * 2, 6);
      numeric::Array1D<brick::common::Float64> BVector(numberOfPoints * 2);

      // Zero out AMatrix to start with, then populate it.
      AMatrix = brick::common::Float64(0);
      unsigned int rowNumber = 0;
      while(toPointsBegin != toPointsEnd) {
        AMatrix(rowNumber, 0) = fromPointsBegin->getX();
        AMatrix(rowNumber, 1) = fromPointsBegin->getY();
        AMatrix(rowNumber, 2) = 1;
        BVector[rowNumber] = toPointsBegin->getX();
        ++rowNumber;
        AMatrix(rowNumber, 3) = fromPointsBegin->getX();
        AMatrix(rowNumber, 4) = fromPointsBegin->getY();
        AMatrix(rowNumber, 5) = 1;
        BVector[rowNumber] = toPointsBegin->getY();
        ++rowNumber;
        ++toPointsBegin;
        ++fromPointsBegin;
      }

      // Solve for the transform.
      brick::numeric::Array1D<brick::common::Float64> result =
        brick::linearAlgebra::linearLeastSquares(AMatrix, BVector);
      return brick::numeric::Transform2D<brick::common::Float64>(
        result[0], result[1], result[2],
        result[3], result[4], result[5],
        brick::common::Float64(0), brick::common::Float64(0),
        brick::common::Float64(1));
    }


    // This function subsamples its input to create a new, smaller
    // image.
    template<ImageFormat Format>
    Image<Format>
    subsample(const Image<Format>& inputImage,
              size_t rowStep,
              size_t columnStep)
    {
      // Argument checking.
      if(inputImage.size() == 0) {
        return Image<Format>();
      }
      if(rowStep == 0) {
        rowStep = 1;
      }
      if(columnStep == 0) {
        columnStep = 1;
      }

      // Output image size is at least (1, 1).  All rows and columns
      // must map to valid positions in the input image.
      Image<Format> outputImage((inputImage.rows() - 1) / rowStep + 1,
                                (inputImage.columns() - 1) / columnStep + 1);
      for(size_t outputRowIndex = 0; outputRowIndex < outputImage.rows();
          ++outputRowIndex) {
        brick::numeric::Array1D<typename Image<Format>::PixelType> inputRow =
          inputImage.row(rowStep * outputRowIndex);
        brick::numeric::Array1D<typename Image<Format>::PixelType> outputRow =
          outputImage.row(outputRowIndex);
        size_t inputColumnIndex = 0;
        for(size_t outputColumnIndex = 0;
            outputColumnIndex < outputImage.columns();
            ++outputColumnIndex) {
          outputRow[outputColumnIndex] = inputRow[inputColumnIndex];
          inputColumnIndex += columnStep;
        }
      }
      return outputImage;
    }


    // This function interpolates its input to create a new, larger
    // image.
    template<ImageFormat InputFormat,
             ImageFormat OutputFormat,
             ImageFormat IntermediateFormat>
    Image<OutputFormat>
    supersample(const Image<InputFormat>& inputImage)
    {
      typedef typename Image<InputFormat>::PixelType InputType;
      typedef typename Image<OutputFormat>::PixelType OutputType;
      typedef typename Image<IntermediateFormat>::PixelType IntermediateType;

      Image<OutputFormat> outputImage(inputImage.rows() * 2 - 1,
                                      inputImage.columns() * 2 - 1);
      ColorspaceConverter<InputFormat, OutputFormat> converter;

      // Working variables.
      IntermediateType accumulator;
      brick::numeric::Array1D<InputType> inputRow0;
      brick::numeric::Array1D<InputType> inputRow1;
      brick::numeric::Array1D<OutputType> outputRow0;
      brick::numeric::Array1D<OutputType> outputRow1;

      // The main supersampling loop.
      for(size_t inputRowIndex = 0; inputRowIndex < (inputImage.rows() - 1);
          ++inputRowIndex) {
        inputRow0 = inputImage.row(inputRowIndex);
        inputRow1 = inputImage.row(inputRowIndex + 1);
        outputRow0 = outputImage.row(2 * inputRowIndex);
        outputRow1 = outputImage.row(2 * inputRowIndex + 1);
        size_t outputColumnIndex = 0;
        for(size_t inputColumnIndex = 0;
            inputColumnIndex < (inputImage.columns() - 1);
            ++inputColumnIndex) {
          outputRow0[outputColumnIndex] =
            converter(inputRow0[inputColumnIndex]);

          // Note(xxx): This can be made more efficient.
          accumulator = (IntermediateType(inputRow0[inputColumnIndex])
                         + IntermediateType(inputRow0[inputColumnIndex + 1]));
          accumulator /= 2;
          outputRow0[outputColumnIndex + 1] = accumulator;

          accumulator = (IntermediateType(inputRow0[inputColumnIndex])
                         + IntermediateType(inputRow1[inputColumnIndex]));
          accumulator /= 2;
          outputRow1[outputColumnIndex] = accumulator;

          accumulator = (IntermediateType(inputRow0[inputColumnIndex])
                         + IntermediateType(inputRow0[inputColumnIndex + 1])
                         + IntermediateType(inputRow1[inputColumnIndex])
                         + IntermediateType(inputRow1[inputColumnIndex + 1]));
          accumulator /= 4;
          outputRow1[outputColumnIndex + 1] = accumulator;

          outputColumnIndex += 2;
        }
        // Handle pixels at end of rows.
        outputRow0[outputColumnIndex] =
          converter(inputRow0[inputRow0.size() - 1]);

        accumulator = (IntermediateType(inputRow0[inputRow0.size() - 1])
                       + IntermediateType(inputRow1[inputRow0.size() - 1]));
        accumulator /= 2;
        outputRow1[outputColumnIndex] = accumulator;
      }

      // Handle final row.
      inputRow0 = inputImage.row(inputImage.rows() - 1);
      outputRow0 = outputImage.row(outputImage.rows() - 1);
      size_t outputColumnIndex = 0;
      for(size_t inputColumnIndex = 0;
          inputColumnIndex < (inputImage.columns() - 1);
          ++inputColumnIndex) {
        outputRow0[outputColumnIndex] =
          converter(inputRow0[inputColumnIndex]);

        accumulator = (IntermediateType(inputRow0[inputColumnIndex])
                       + IntermediateType(inputRow0[inputColumnIndex + 1]));
        accumulator /= 2;
        outputRow0[outputColumnIndex + 1] = accumulator;

        outputColumnIndex += 2;
      }
      outputRow0[outputImage.columns() - 1] =
        converter(inputRow0[inputImage.columns() - 1]);

      return outputImage;
    }


    // This function creates a new array and copies into it the pixel
    // data from the input image.
    template <class Type, ImageFormat FORMAT>
    brick::numeric::Array2D<Type>
    toArray(const Image<FORMAT>& inputImage)
    {
      // A typedef for notational convenience.
      typedef typename ImageFormatTraits<FORMAT>::ComponentType ComponentType;

      size_t numberOfOutputColumns = inputImage.columns()
        * ImageFormatTraits<FORMAT>::getNumberOfComponents();

      brick::numeric::Array2D<ComponentType> tempArray;
      dissociateColorComponents(
        const_cast< Image<FORMAT>& >(inputImage), tempArray);

      brick::numeric::Array2D<Type> returnValue(
        inputImage.rows(), numberOfOutputColumns);
      returnValue.copy(tempArray);
      return returnValue;
    }

  } // namespace computerVision

} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_UTILITIES_IMPL_HH */
