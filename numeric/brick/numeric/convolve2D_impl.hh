/**
***************************************************************************
* @file brick/numeric/convolve2D_impl.hh
*
* Header file defining 1D correlation and convolution functions.
*
* Copyright (C) 2006-2007, 2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_CONVOLVE2D_IMPL_HH
#define BRICK_NUMERIC_CONVOLVE2D_IMPL_HH

// This file is included by convolve2D.hh, and should not be directly included
// by user code, so no need to include convolve2D.hh here.
// 
// #include <brick/numeric/convolve2D.hh>

#include <algorithm> // For std::reverse_copy()
#include <numeric> // For std::partial_sum()
#include <brick/common/functional.hh>
#include <brick/numeric/numericTraits.hh>
#include <brick/numeric/stencil2D.hh>

namespace brick {

  namespace numeric {

    /// @cond privateCode    
    namespace privateCode {

      template <class AccumulatorType, class KernelType, class SignalType>
      Array2D<AccumulatorType>
      doubleAccumulateKernel(const Array2D<KernelType>& kernel,
			     SignalType fillValue)
      {
	typedef typename ArithmeticTraits<KernelType, KernelType>::SumType
          KernelSumType;
	Array2D<KernelSumType> tempKernel(
	  kernel.rows() + 1, kernel.columns() + 1);
	std::fill(tempKernel.rowBegin(0), tempKernel.rowEnd(0),
		  static_cast<KernelSumType>(0));
	for(size_t row = 0; row < kernel.rows(); ++row) {
	  tempKernel(row + 1, 0) = static_cast<KernelSumType>(0);
	  std::copy(kernel.rowBegin(row), kernel.rowEnd(row),
		    tempKernel.rowBegin(row + 1) + 1);		    
	}

	// Accumulate horizontally.
	Array2D<KernelSumType> accumulatedKernel(
	  kernel.rows() + 1, kernel.columns() + 1);
	std::partial_sum(tempKernel.begin(), tempKernel.end(),
			 accumulatedKernel.begin());

	// Accumulate vertically.
	tempKernel = accumulatedKernel.transpose();
	accumulatedKernel.reshape(
	  static_cast<int>(kernel.columns()) + 1,
	  static_cast<int>(kernel.rows()) + 1);
	std::partial_sum(tempKernel.begin(), tempKernel.end(),
			 accumulatedKernel.begin());

	Array2D<AccumulatorType> result(
	  kernel.rows() + 1, kernel.columns() + 1);
	for(size_t row = 0; row < result.rows(); ++row) {
	  for(size_t column = 0; column < result.columns(); ++column) {
	    result(row, column) =
	      static_cast<AccumulatorType>(
		accumulatedKernel(column, row)
		* static_cast<AccumulatorType>(fillValue));
	  }
	}
	return result;
      }


      template<class OutputType, class InputType>
      inline OutputType
      integrateArea(const Array2D<InputType>& accumulatedKernel,
		    size_t row0, size_t row1,
		    size_t column0, size_t column1)
      {
	return static_cast<OutputType>(
	  accumulatedKernel(row1, column1)
	  - accumulatedKernel(row0, column1)
	  - accumulatedKernel(row1, column0)
	  + accumulatedKernel(row0, column0));
      }
      
      
      template <class OutputType, class AccumulatorType,
		class KernelType, class SignalType,
		size_t StencilSize>
      void
      sizedCorrelate2DCommon(const Array2D<KernelType>& kernel,
			     const Array2D<SignalType>& signal,
			     Array2D<OutputType>& result,
			     const Index2D& corner0, const Index2D& corner1,
			     const Index2D& resultCorner0)
      {
	typedef StencilIterator<const SignalType, StencilSize> SignalIterator;
	typedef typename Array2D<KernelType>::const_iterator KernelIterator;

	if(kernel.size() == 0) {
	  BRICK_THROW(brick::common::ValueException, "sizedCorrelate2DCommon()",
		    "Argument kernel has zero size.");
	}
	if(signal.size() == 0) {
	  BRICK_THROW(brick::common::ValueException, "sizedCorrelate2DCommon()",
		    "Argument signal has zero size.");
	}

	Stencil2D<const SignalType, StencilSize> signalStencil(
	  kernel.rows(), kernel.columns());
	signalStencil.setTarget(signal);

	const size_t startRow = corner0.getRow() - kernel.rows() / 2;
	const size_t stopRow = corner1.getRow() - kernel.rows() / 2;
	const size_t startColumn = corner0.getColumn() - kernel.columns() / 2;
	const size_t stopColumn = corner1.getColumn() - kernel.columns() / 2;
	const size_t resultRow0 = resultCorner0.getRow();
	const size_t resultColumn0 = resultCorner0.getColumn();
	for(size_t row = startRow; row < stopRow; ++row) {
	  size_t resultIndex =
	    ((resultRow0 + row - startRow) * result.columns() + resultColumn0);
	  signalStencil.goTo(row, startColumn);
	  for(size_t column = startColumn; column < stopColumn; ++column) {
	    AccumulatorType dotProduct = static_cast<AccumulatorType>(0);
	    KernelIterator kernelIter = kernel.begin();
	    KernelIterator endIter = kernel.end();
	    SignalIterator signalIter = signalStencil.begin();
	    while(kernelIter != endIter) {
	      dotProduct +=
		static_cast<AccumulatorType>(
		  *kernelIter * static_cast<AccumulatorType>(*signalIter));
	      ++kernelIter;
	      ++signalIter;
	    }
	    result(resultIndex) = static_cast<OutputType>(dotProduct);
	    ++resultIndex;
	    signalStencil.advance();
	  }
	}
      }


      template <class OutputType, class AccumulatorType,
		class KernelType, class SignalType>
      void
      correlate2DCommon(const Array2D<KernelType>& kernel,
			const Array2D<SignalType>& signal,
			Array2D<OutputType>& result,
			const Index2D& corner0, const Index2D& corner1,
			const Index2D& resultCorner0)
      {
	if(kernel.size() <= 9) {
	  sizedCorrelate2DCommon<
	    OutputType, AccumulatorType, KernelType, SignalType, 9>(
	      kernel, signal, result, corner0, corner1, resultCorner0);
	} else if(kernel.size() <= 25) {
	  sizedCorrelate2DCommon<
	    OutputType, AccumulatorType, KernelType, SignalType, 25>(
	      kernel, signal, result, corner0, corner1, resultCorner0);
	} else if(kernel.size() <= 49) {
	  sizedCorrelate2DCommon<
	    OutputType, AccumulatorType, KernelType, SignalType, 49>(
	      kernel, signal, result, corner0, corner1, resultCorner0);
	} else if(kernel.size() <= 81) {
	  sizedCorrelate2DCommon<
	    OutputType, AccumulatorType, KernelType, SignalType, 81>(
	      kernel, signal, result, corner0, corner1, resultCorner0);
	} else if(kernel.size() <= 121) {
	  sizedCorrelate2DCommon<
	    OutputType, AccumulatorType, KernelType, SignalType, 121>(
	      kernel, signal, result, corner0, corner1, resultCorner0);
	} else if(kernel.size() <= 255) {
	  sizedCorrelate2DCommon<
	    OutputType, AccumulatorType, KernelType, SignalType, 255>(
	      kernel, signal, result, corner0, corner1, resultCorner0);
	} else if(kernel.size() <= 1023) {
	  sizedCorrelate2DCommon<
	    OutputType, AccumulatorType, KernelType, SignalType, 1023>(
	      kernel, signal, result, corner0, corner1, resultCorner0);
	} else if(kernel.size() <= 4095) {
	  sizedCorrelate2DCommon<
	    OutputType, AccumulatorType, KernelType, SignalType, 4095>(
	      kernel, signal, result, corner0, corner1, resultCorner0);
	} else if(kernel.size() <= 16383) {
	  sizedCorrelate2DCommon<
	    OutputType, AccumulatorType, KernelType, SignalType, 16383>(
	      kernel, signal, result, corner0, corner1, resultCorner0);
	} else if(kernel.size() <= 65535) {
	  sizedCorrelate2DCommon<
	    OutputType, AccumulatorType, KernelType, SignalType, 65535>(
	      kernel, signal, result, corner0, corner1, resultCorner0);
	} else if(kernel.size() <= 262143) {
	  sizedCorrelate2DCommon<
	    OutputType, AccumulatorType, KernelType, SignalType, 262143>(
	      kernel, signal, result, corner0, corner1, resultCorner0);
	} else if(kernel.size() <= 1048575) {
	  sizedCorrelate2DCommon<
	    OutputType, AccumulatorType, KernelType, SignalType, 1048575>(
	      kernel, signal, result, corner0, corner1, resultCorner0);
	} else {
	  BRICK_THROW(
            brick::common::NotImplementedException, "correlate2DCommon()",
            "Kernels with more than 2^24 elements are not supported.");
	}
      }

      
      template <class OutputType, class AccumulatorType,
		class KernelType, class SignalType>
      inline OutputType
      innerProduct2D(const Array2D<KernelType>& kernel,
		     const Array2D<SignalType>& signal,
		     size_t kernelRow0, size_t kernelRow1,
		     size_t kernelColumn0, size_t kernelColumn1,
		     size_t signalRow0, size_t signalColumn0)
      {
	typedef typename Array2D<KernelType>::const_iterator KernelIterator;
	typedef typename Array2D<SignalType>::const_iterator SignalIterator;

	size_t overlapWidth = kernelColumn1 - kernelColumn0;
	AccumulatorType dotProduct = static_cast<AccumulatorType>(0);
	size_t signalRow = signalRow0;
	for(size_t kernelRow = kernelRow0; kernelRow < kernelRow1;
	    ++kernelRow) {
	  KernelIterator kernelIter =
	    kernel.rowBegin(kernelRow) + kernelColumn0;
	  KernelIterator kernelEnd = kernelIter + overlapWidth;
	  SignalIterator signalIter =
	    signal.rowBegin(signalRow) + signalColumn0;
	  while(kernelIter != kernelEnd) {
	    dotProduct += static_cast<AccumulatorType>(
	      *kernelIter * static_cast<AccumulatorType>(*signalIter));
	    ++kernelIter;
	    ++signalIter;
	  }
	  ++signalRow;
	}
	return static_cast<OutputType>(dotProduct);
      }
      

      template <class OutputType, class AccumulatorType,
		class KernelType, class SignalType>
      inline OutputType
      correlate2DLeftRightWrap(const Array2D<KernelType>& kernel,
			       const Array2D<SignalType>& signal,
			       size_t kernelRow0, size_t kernelRow1,
			       size_t signalRow0,
			       size_t leftOverlap)
      {
	AccumulatorType result;
	size_t kernelStartColumn = kernel.columns() - leftOverlap;
	size_t signalRightStartColumn = signal.columns() - kernelStartColumn;
	result = innerProduct2D<
	  AccumulatorType, AccumulatorType, KernelType, SignalType>(
	    kernel, signal, kernelRow0, kernelRow1,
	    kernelStartColumn, kernel.columns(), signalRow0, 0);
	result += innerProduct2D<
	  AccumulatorType, AccumulatorType, KernelType, SignalType>(
	    kernel, signal, kernelRow0, kernelRow1,
	    0, kernelStartColumn, signalRow0, signalRightStartColumn);
	return static_cast<OutputType>(result);
      }
      

      template <class OutputType, class AccumulatorType,
		class KernelType, class SignalType>
      inline OutputType
      correlate2DLeftReflection(const Array2D<KernelType>& kernel,
				const Array2D<SignalType>& signal,
				size_t kernelRow0, size_t kernelRow1,
				size_t signalRow0,
				size_t leftOverlap)
      {
	typedef typename Array2D<KernelType>::const_iterator KernelIterator;
	typedef typename Array2D<SignalType>::const_iterator SignalIterator;

	AccumulatorType dotProduct = static_cast<AccumulatorType>(0);
	size_t signalRow = signalRow0;
	size_t reflectionSize = kernel.columns() - leftOverlap;
	for(size_t kernelRow = kernelRow0; kernelRow < kernelRow1;
	    ++kernelRow) {
	  // Compute reflected component.
	  KernelIterator kernelIter = kernel.rowBegin(kernelRow);
	  KernelIterator kernelMiddle = kernelIter + reflectionSize;
	  KernelIterator kernelEnd = kernel.rowEnd(kernelRow);
	  SignalIterator signalIter =
	    signal.rowBegin(signalRow) + (reflectionSize - 1);
	  while(kernelIter != kernelMiddle) {
	    dotProduct += static_cast<AccumulatorType>(
	      *kernelIter * static_cast<AccumulatorType>(*signalIter));
	    ++kernelIter;
	    --signalIter;
	  }

	  // Compute in-bounds component.
	  ++signalIter;
	  while(kernelIter != kernelEnd) {
	    dotProduct += static_cast<AccumulatorType>(
	      *kernelIter * static_cast<AccumulatorType>(*signalIter));
	    ++kernelIter;
	    ++signalIter;
	  }

	  ++signalRow;
	}
	return static_cast<OutputType>(dotProduct);
      }
      

      template <class OutputType, class AccumulatorType,
		class KernelType, class SignalType>
      inline OutputType
      correlate2DRightReflection(const Array2D<KernelType>& kernel,
				 const Array2D<SignalType>& signal,
				 size_t kernelRow0, size_t kernelRow1,
				 size_t signalRow0,
				 size_t rightOverlap)
      {
	typedef typename Array2D<KernelType>::const_iterator KernelIterator;
	typedef typename Array2D<SignalType>::const_iterator SignalIterator;

	AccumulatorType dotProduct = static_cast<AccumulatorType>(0);
	size_t signalRow = signalRow0;
	for(size_t kernelRow = kernelRow0; kernelRow < kernelRow1;
	    ++kernelRow) {
	  KernelIterator kernelIter = kernel.rowBegin(kernelRow);
	  KernelIterator kernelMiddle = kernelIter + rightOverlap;
	  KernelIterator kernelEnd = kernel.rowEnd(kernelRow);
	  SignalIterator signalIter = signal.rowEnd(signalRow) - rightOverlap;
	  
	  // Compute in-bounds component.
	  while(kernelIter != kernelMiddle) {
	    dotProduct += static_cast<AccumulatorType>(
	      *kernelIter * static_cast<AccumulatorType>(*signalIter));
	    ++kernelIter;
	    ++signalIter;
	  }

	  // Compute reflected component.
	  --signalIter;
	  while(kernelIter != kernelEnd) {
	    dotProduct += static_cast<AccumulatorType>(
	      *kernelIter * static_cast<AccumulatorType>(*signalIter));
	    ++kernelIter;
	    --signalIter;
	  }
	  ++signalRow;
	}
	return static_cast<OutputType>(dotProduct);
      }
      

      template <class OutputType, class AccumulatorType,
		class KernelType, class SignalType>
      Array2D<OutputType>
      correlate2DTruncateResult(const Array2D<KernelType>& kernel,
				const Array2D<SignalType>& signal)
      {
        Array2D<OutputType> result(signal.rows() - kernel.rows() + 1,
				   signal.columns() - kernel.columns() + 1);
	Index2D corner0(static_cast<int>(kernel.rows()) / 2,
			static_cast<int>(kernel.columns()) / 2);
	Index2D corner1(
	  static_cast<int>(signal.rows())
	  - static_cast<int>(kernel.rows()) / 2,
	  static_cast<int>(signal.columns())
	  - static_cast<int>(kernel.columns()) / 2);
	Index2D resultCorner0(0, 0);
        correlate2DCommon<OutputType, AccumulatorType, KernelType, SignalType>(
	  kernel, signal, result, corner0, corner1, resultCorner0);
        return result;
      }


      template <class OutputType, class AccumulatorType,
		class KernelType, class SignalType>
      Array2D<OutputType>
      correlate2DPadResult(const Array2D<KernelType>& kernel,
			   const Array2D<SignalType>& signal,
			   const Index2D& corner0,
			   const Index2D& corner1,
			   OutputType fillValue)
      {
	typedef typename Array2D<OutputType>::iterator ResultIterator;

	size_t outputRows =
	  corner1.getRow() - corner0.getRow();
	size_t outputColumns =
	  corner1.getColumn() - corner0.getColumn();
        Array2D<OutputType> result(outputRows, outputColumns);

	// Check for degenerate case.
	if(kernel.rows() > signal.rows()
	   || kernel.columns() > signal.columns()) {
	  result = fillValue;
	  return result;
	}

	// Establish significant row and column coordinates in the
	// output image.

	// Constants specified without reference to signal or result.
	const int kRowOverTwo = static_cast<int>(kernel.rows()) / 2;
	const int kColOverTwo = static_cast<int>(kernel.columns()) / 2;
	
	// Constants specified with respect to rows & columns in
	// argument signal.
	const int startRow = corner0.getRow();
	const int stopRow = corner1.getRow();
	const int transitionRow0 = kRowOverTwo;
	const int transitionRow1 =
	  static_cast<int>(signal.rows()) - transitionRow0;
	const int startColumn = corner0.getColumn();
	const int stopColumn = corner1.getColumn();
	const int transitionColumn0 = kColOverTwo;
	const int transitionColumn1 =
	  static_cast<int>(signal.columns()) - transitionColumn0;
	const int clippedTransitionRow0 =
	  brick::common::clip(transitionRow0, startRow, stopRow);
	const int clippedTransitionRow1 =
	  brick::common::clip(transitionRow1, startRow, stopRow);
	const int clippedTransitionColumn0 =
	  brick::common::clip(transitionColumn0, startColumn, stopColumn);
	const int clippedTransitionColumn1 =
	  brick::common::clip(transitionColumn1, startColumn, stopColumn);

	// Constants specified with respect to rows & columns in
	// result.
	const int resultTransitionRow0 =
	  clippedTransitionRow0 - startRow;
	const int resultTransitionColumn0 =
	  clippedTransitionColumn0 - startColumn;

	// Fill in areas along top of image.
	int row = startRow;
	int outputRow = 0;
        while(row < clippedTransitionRow0) {
	  ResultIterator rowBegin = result.rowBegin(outputRow);
	  ResultIterator rowEnd = result.rowEnd(outputRow);
	  std::fill(rowBegin, rowEnd, fillValue);
	  ++row;
	  ++outputRow;
	}

	// Fill in areas along sides.
	int leftBorderIndex = clippedTransitionColumn0 - startColumn;
	int rightBorderIndex = clippedTransitionColumn1 - startColumn;
	while(row < clippedTransitionRow1) {
	  ResultIterator rowBegin = result.rowBegin(outputRow);
	  ResultIterator rowEnd = result.rowEnd(outputRow);
	  std::fill(rowBegin, rowBegin + leftBorderIndex, fillValue);
	  std::fill(rowBegin + rightBorderIndex, rowEnd, fillValue);
	  ++row;
	  ++outputRow;
	}

	// Fill in areas along bottom of image.
        while(row < stopRow) {
	  ResultIterator rowBegin = result.rowBegin(outputRow);
	  ResultIterator rowEnd = result.rowEnd(outputRow);
	  std::fill(rowBegin, rowEnd, fillValue);
	  ++row;
	  ++outputRow;
	}
	
	Index2D signalCorner0(clippedTransitionRow0, clippedTransitionColumn0);
	Index2D signalCorner1(clippedTransitionRow1, clippedTransitionColumn1);
	Index2D resultCorner0(resultTransitionRow0, resultTransitionColumn0);
        correlate2DCommon<OutputType, AccumulatorType, KernelType, SignalType>(
	  kernel, signal, result, signalCorner0, signalCorner1, resultCorner0);
        return result;
      }
    
                   
      template <class OutputType, class AccumulatorType,
		class KernelType, class SignalType>
      Array2D<OutputType>
      correlate2DZeroPadSignal(const Array2D<KernelType>& kernel,
			       const Array2D<SignalType>& signal,
			       const Index2D& corner0,
			       const Index2D& corner1)
      {
	size_t outputRows =
	  corner1.getRow() - corner0.getRow();
	size_t outputColumns =
	  corner1.getColumn() - corner0.getColumn();
        Array2D<OutputType> result(outputRows, outputColumns);

	// Constants specified without reference to signal or result.
	const int kRowOverTwo = static_cast<int>(kernel.rows()) / 2;
	const int kColOverTwo = static_cast<int>(kernel.columns()) / 2;

	// Constants specified with respect to rows & columns in
	// argument signal.
	const int startRow = corner0.getRow();
	const int stopRow = corner1.getRow();
	const int transitionRow0 = kRowOverTwo;
	const int transitionRow1
	  = static_cast<int>(signal.rows()) - transitionRow0;
	const int startColumn = corner0.getColumn();
	const int stopColumn = corner1.getColumn();
	const int transitionColumn0 = kColOverTwo;
	const int transitionColumn1 =
	  static_cast<int>(signal.columns()) - transitionColumn0;
	const int clippedTransitionRow0 =
	  brick::common::clip(transitionRow0, startRow, stopRow);
	const int clippedTransitionRow1 =
	  brick::common::clip(transitionRow1, startRow, stopRow);
	const int clippedTransitionColumn0 =
	  brick::common::clip(transitionColumn0, startColumn, stopColumn);
	const int clippedTransitionColumn1 =
	  brick::common::clip(transitionColumn1, startColumn, stopColumn);

	// Fill in areas along top of image.
	int row = startRow;
	int outputRow = 0;
        while(row < clippedTransitionRow0) {
	  // Fill in top left corner.
	  int column = startColumn;
	  int outputColumn = 0;
	  while(column < clippedTransitionColumn0) {
	    result(outputRow, outputColumn) = innerProduct2D<
	      OutputType, AccumulatorType, KernelType, SignalType>(
		kernel, signal, transitionRow0 - row, kernel.rows(),
		transitionColumn0 - column, kernel.columns(), 0, 0);
	    ++column;
	    ++outputColumn;
	  }

	  // Fill in top edge.
	  while(column < clippedTransitionColumn1) {
	    result(outputRow, outputColumn) = innerProduct2D<
	      OutputType, AccumulatorType, KernelType, SignalType>(
		kernel, signal, transitionRow0 - row, kernel.rows(),
		0, kernel.columns(), 0, column - transitionColumn0);
	    ++column;
	    ++outputColumn;
	  }

	  // Fill in top right corner.
	  while(column < stopColumn) {
	    result(outputRow, outputColumn) = innerProduct2D<
	      OutputType, AccumulatorType, KernelType, SignalType>(
		kernel, signal, transitionRow0 - row, kernel.rows(),
		0, signal.columns() - column + kColOverTwo,
		0, column - transitionColumn0);
	    ++column;
	    ++outputColumn;
	  }

	  ++row;
	  ++outputRow;
	}

	// Fill in areas along left and right edges of image.
	while(row < clippedTransitionRow1) {
	  // Fill in left edge.
	  int column = startColumn;
	  int outputColumn = 0;
	  while(column < clippedTransitionColumn0) {
	    result(outputRow, outputColumn) = innerProduct2D<
	      OutputType, AccumulatorType, KernelType, SignalType>(
		kernel, signal, 0, kernel.rows(),
		transitionColumn0 - column, kernel.columns(),
		row - transitionRow0, 0);
	    ++column;
	    ++outputColumn;
	  }
	    
	  // Fill in right edge.
	  column = clippedTransitionColumn1;
	  outputColumn +=
	    clippedTransitionColumn1 - clippedTransitionColumn0;
	  while(column < stopColumn) {
	    result(outputRow, outputColumn) = innerProduct2D<
	      OutputType, AccumulatorType, KernelType, SignalType>(
		kernel, signal, 0, kernel.rows(),
		0, signal.columns() - column + kColOverTwo,
		row - transitionRow0, column - transitionColumn0);
	    ++column;
	    ++outputColumn;
	  }
	  ++row;
	  ++outputRow;
	}
	
	// Fill in areas along bottom of image.
        while(row < stopRow) {
	  // Fill in bottom left corner.
	  int column = startColumn;
	  int outputColumn = 0;
	  while(column < clippedTransitionColumn0) {
	    result(outputRow, outputColumn) = innerProduct2D<
	      OutputType, AccumulatorType, KernelType, SignalType>(
		kernel, signal, 0, signal.rows() - row + kRowOverTwo,
		transitionColumn0 - column, kernel.columns(),
		row - transitionRow0, 0);
	    ++column;
	    ++outputColumn;
	  }

	  // Fill in bottom edge.
	  while(column < clippedTransitionColumn1) {
	    result(outputRow, outputColumn) = innerProduct2D<
	      OutputType, AccumulatorType, KernelType, SignalType>(
		kernel, signal, 0, signal.rows() - row + kRowOverTwo,
		0, kernel.columns(),
		row - transitionRow0, column - transitionColumn0);
	    ++column;
	    ++outputColumn;
	  }

	  // Fill in bottom right corner.
	  while(column < stopColumn) {
	    result(outputRow, outputColumn) = innerProduct2D<
	      OutputType, AccumulatorType, KernelType, SignalType>(
		kernel, signal, 0, signal.rows() - row + kRowOverTwo,
		0, signal.columns() - column + kColOverTwo,
		row - transitionRow0, column - transitionColumn0);
	    ++column;
	    ++outputColumn;
	  }
	  ++row;
	  ++outputRow;
	}

	// Fill in the middle of the image.
	Index2D newCorner0(kRowOverTwo, kColOverTwo);
	Index2D newCorner1(static_cast<int>(signal.rows()) - kRowOverTwo,
			   static_cast<int>(signal.columns()) - kColOverTwo);
	Index2D resultCorner0(static_cast<int>(kernel.rows()) - 1,
			      static_cast<int>(kernel.columns()) - 1);
        correlate2DCommon<OutputType, AccumulatorType, KernelType, SignalType>(
	  kernel, signal, result, newCorner0, newCorner1, resultCorner0);
        return result;
      }


      template <class OutputType, class AccumulatorType,
		class KernelType, class SignalType>
      Array2D<OutputType>
      correlate2DPadSignal(const Array2D<KernelType>& kernel,
			   const Array2D<SignalType>& signal,
			   const Index2D& corner0,
			   const Index2D& corner1,
			   SignalType fillValue)
      {
	size_t outputRows =
	  corner1.getRow() - corner0.getRow();
	size_t outputColumns =
	  corner1.getColumn() - corner0.getColumn();
        Array2D<OutputType> result(outputRows, outputColumns);

	// Precompute numbers which will make it easy to total up the
	// portions of the kernel which overlap padded areas of the
	// signal.
	Array2D<AccumulatorType> accumulatedKernel =
	  doubleAccumulateKernel<AccumulatorType, KernelType, SignalType>(
	    kernel, fillValue);
	
	// Constants specified without reference to signal or result.
	const int kRowOverTwo = static_cast<int>(kernel.rows()) / 2;
	const int kColOverTwo = static_cast<int>(kernel.columns()) / 2;
	
	// Constants specified with respect to rows & columns in
	// argument signal.
	const int startRow = corner0.getRow();
	const int stopRow = corner1.getRow();
	const int transitionRow0 = kRowOverTwo;
	const int transitionRow1 =
	  static_cast<int>(signal.rows()) - transitionRow0;
	const int startColumn = corner0.getColumn();
	const int stopColumn = corner1.getColumn();
	const int transitionColumn0 = kColOverTwo;
	const int transitionColumn1 =
	  static_cast<int>(signal.columns()) - transitionColumn0;
	const int clippedTransitionRow0 =
	  brick::common::clip(transitionRow0, startRow, stopRow);
	const int clippedTransitionRow1 =
	  brick::common::clip(transitionRow1, startRow, stopRow);
	const int clippedTransitionColumn0 =
	  brick::common::clip(transitionColumn0, startColumn, stopColumn);
	const int clippedTransitionColumn1 =
	  brick::common::clip(transitionColumn1, startColumn, stopColumn);

	// Fill in areas along top of image.
	int row = startRow;
	int outputRow = 0;
        while(row < clippedTransitionRow0) {
	  // Fill in top left corner.
	  int column = startColumn;
	  int outputColumn = 0;
	  while(column < clippedTransitionColumn0) {
	    result(outputRow, outputColumn) = static_cast<OutputType>(
	      innerProduct2D<
	      AccumulatorType, AccumulatorType, KernelType, SignalType>(
		kernel, signal, transitionRow0 - row, kernel.rows(),
		transitionColumn0 - column, kernel.columns(), 0, 0)
	      + integrateArea<AccumulatorType>(
		accumulatedKernel, 0, transitionRow0 - row,
		0, kernel.columns())
	      + integrateArea<AccumulatorType>(
		accumulatedKernel, transitionRow0 - row, kernel.rows(),
		0, transitionColumn0 - column));
	    ++column;
	    ++outputColumn;
	  }

	  // Fill in top edge.
	  while(column < clippedTransitionColumn1) {
	    result(outputRow, outputColumn) = static_cast<OutputType>(
	      innerProduct2D<
	      AccumulatorType, AccumulatorType, KernelType, SignalType>(
		kernel, signal, transitionRow0 - row, kernel.rows(),
		0, kernel.columns(), 0, column - transitionColumn0)
	      + integrateArea<AccumulatorType>(
		accumulatedKernel, 0, transitionRow0 - row,
		0, kernel.columns()));
	    ++column;
	    ++outputColumn;
	  }

	  // Fill in top right corner.
	  while(column < stopColumn) {
	    result(outputRow, outputColumn) = static_cast<OutputType>(
	      innerProduct2D<
	      AccumulatorType, AccumulatorType, KernelType, SignalType>(
		kernel, signal, transitionRow0 - row, kernel.rows(),
		0, signal.columns() - column + kColOverTwo,
		0, column - transitionColumn0)
	      + integrateArea<AccumulatorType>(
		accumulatedKernel, 0, transitionRow0 - row,
		0, kernel.columns())
	      + integrateArea<AccumulatorType>(
		accumulatedKernel, transitionRow0 - row, kernel.rows(),
		signal.columns() - column + kColOverTwo, kernel.columns()));
	    ++column;
	    ++outputColumn;
	  }

	  ++row;
	  ++outputRow;
	}

	// Fill in areas along left and right edges of image.
	while(row < clippedTransitionRow1) {
	  // Fill in left edge.
	  int column = startColumn;
	  int outputColumn = 0;
	  while(column < clippedTransitionColumn0) {
	    result(outputRow, outputColumn) = static_cast<OutputType>(
	      innerProduct2D<
	      AccumulatorType, AccumulatorType, KernelType, SignalType>(
		kernel, signal, 0, kernel.rows(),
		transitionColumn0 - column, kernel.columns(),
		row - transitionRow0, 0)
	      + integrateArea<AccumulatorType>(
		accumulatedKernel, 0, kernel.rows(),
		0, transitionColumn0 - column));
	    ++column;
	    ++outputColumn;
	  }
	    
	  // Fill in right edge.
	  column = clippedTransitionColumn1;
	  outputColumn +=
	    clippedTransitionColumn1 - clippedTransitionColumn0;
	  while(column < stopColumn) {
	    result(outputRow, outputColumn) = static_cast<OutputType>(
	      innerProduct2D<
	      AccumulatorType, AccumulatorType, KernelType, SignalType>(
		kernel, signal, 0, kernel.rows(),
		0, signal.columns() - column + kColOverTwo,
		row - transitionRow0, column - transitionColumn0)
	      + integrateArea<AccumulatorType>(
		accumulatedKernel, 0, kernel.rows(),
		signal.columns() - column + kColOverTwo, kernel.columns()));
	    ++column;
	    ++outputColumn;
	  }
	  ++row;
	  ++outputRow;
	}
	
	// Fill in areas along bottom of image.
        while(row < stopRow) {
	  // Fill in bottom left corner.
	  int column = startColumn;
	  int outputColumn = 0;
	  while(column < clippedTransitionColumn0) {
	    result(outputRow, outputColumn) = static_cast<OutputType>(
	      innerProduct2D<
	      AccumulatorType, AccumulatorType, KernelType, SignalType>(
		kernel, signal, 0, signal.rows() - row + kRowOverTwo,
		transitionColumn0 - column, kernel.columns(),
		row - transitionRow0, 0)
	      + integrateArea<AccumulatorType>(
		accumulatedKernel,
		signal.rows() - row + kRowOverTwo, kernel.rows(),
		0, kernel.columns())
	      + integrateArea<AccumulatorType>(
		accumulatedKernel, 0, signal.rows() - row + kRowOverTwo,
		0, transitionColumn0 - column));
	    ++column;
	    ++outputColumn;
	  }

	  // Fill in bottom edge.
	  while(column < clippedTransitionColumn1) {
	    result(outputRow, outputColumn) = static_cast<OutputType>(
	      innerProduct2D<
	      AccumulatorType, AccumulatorType, KernelType, SignalType>(
		kernel, signal, 0, signal.rows() - row + kRowOverTwo,
		0, kernel.columns(),
		row - transitionRow0, column - transitionColumn0)
	      + integrateArea<AccumulatorType>(
		accumulatedKernel,
		signal.rows() - row + kRowOverTwo, kernel.rows(),
		0, kernel.columns()));
	    ++column;
	    ++outputColumn;
	  }

	  // Fill in bottom right corner.
	  while(column < stopColumn) {
	    result(outputRow, outputColumn) = static_cast<OutputType>(
	      innerProduct2D<
	      AccumulatorType, AccumulatorType, KernelType, SignalType>(
		kernel, signal, 0, signal.rows() - row + kRowOverTwo,
		0, signal.columns() - column + kColOverTwo,
		row - transitionRow0, column - transitionColumn0)
	      + integrateArea<AccumulatorType>(
		accumulatedKernel,
		signal.rows() - row + kRowOverTwo, kernel.rows(),
		0, kernel.columns())
	      + integrateArea<AccumulatorType>(
		accumulatedKernel, 0, signal.rows() - row + kRowOverTwo,
		signal.columns() - column + kColOverTwo, kernel.columns()));
	    ++column;
	    ++outputColumn;
	  }
	  ++row;
	  ++outputRow;
	}

	// Fill in the middle of the image.
	Index2D newCorner0(kRowOverTwo, kColOverTwo);
	Index2D newCorner1(static_cast<int>(signal.rows()) - kRowOverTwo,
			   static_cast<int>(signal.columns()) - kColOverTwo);
	Index2D resultCorner0(static_cast<int>(kernel.rows()) - 1,
			      static_cast<int>(kernel.columns()) - 1);
        correlate2DCommon<OutputType, AccumulatorType, KernelType, SignalType>(
	  kernel, signal, result, newCorner0, newCorner1, resultCorner0);
        return result;
      }


      template <class OutputType, class AccumulatorType,
		class KernelType, class SignalType>
      Array2D<OutputType>
      correlate2DReflectSignal(const Array2D<KernelType>& kernel,
			       const Array2D<SignalType>& signal,
			       const Index2D& corner0,
			       const Index2D& corner1)
      {
	// Brace yourself...

	size_t outputRows =
	  corner1.getRow() - corner0.getRow();
	size_t outputColumns =
	  corner1.getColumn() - corner0.getColumn();
        Array2D<OutputType> result(outputRows, outputColumns);

	// Make a flipped (upside down) kernel for coding convenience.
	Array2D<KernelType> kernelUD(kernel.rows(), kernel.columns());
	for(size_t rowIndex = 0; rowIndex < kernel.rows(); ++rowIndex) {
	  std::copy(
	    kernel.rowBegin(rowIndex), kernel.rowEnd(rowIndex),
	    kernelUD.rowBegin(kernel.rows() - rowIndex - 1));
	}

	// Constants specified without reference to signal or result.
	const int kRowOverTwo = static_cast<int>(kernel.rows()) / 2;
	const int kColOverTwo = static_cast<int>(kernel.columns()) / 2;
	
	// Constants specified with respect to rows & columns in
	// argument signal.
	const int startRow = corner0.getRow();
	const int stopRow = corner1.getRow();
	const int transitionRow0 = kRowOverTwo;
	const int transitionRow1 =
	  static_cast<int>(signal.rows()) - transitionRow0;
	const int startColumn = corner0.getColumn();
	const int stopColumn = corner1.getColumn();
	const int transitionColumn0 = kColOverTwo;
	const int transitionColumn1 =
	  static_cast<int>(signal.columns()) - transitionColumn0;
	const int clippedTransitionRow0 =
	  brick::common::clip(transitionRow0, startRow, stopRow);
	const int clippedTransitionRow1 =
	  brick::common::clip(transitionRow1, startRow, stopRow);
	const int clippedTransitionColumn0 =
	  brick::common::clip(transitionColumn0, startColumn, stopColumn);
	const int clippedTransitionColumn1 =
	  brick::common::clip(transitionColumn1, startColumn, stopColumn);

	// Fill in areas along top of image.
	int row = startRow;
	int outputRow = 0;
        while(row < clippedTransitionRow0) {
	  // Intermediate variable keeping track of the first kernel
	  // row that overlaps valid image data.
	  size_t kernelStartRow = transitionRow0 - row;

	  // Fill in top left corner.
	  int column = startColumn;
	  int outputColumn = 0;
	  while(column < clippedTransitionColumn0) {
	    AccumulatorType elementValue = correlate2DLeftReflection<
	      AccumulatorType, AccumulatorType, KernelType, SignalType>(
		kernel, signal, kernelStartRow, kernel.rows(), 0,
		column + kColOverTwo + 1);
	    elementValue += correlate2DLeftReflection<
	      AccumulatorType, AccumulatorType, KernelType, SignalType>(
		kernelUD, signal,
		kernel.rows() - kernelStartRow, kernel.rows(), 0,
		column + kColOverTwo + 1);
	    result(outputRow, outputColumn) =
	      static_cast<OutputType>(elementValue);
	    ++column;
	    ++outputColumn;
	  }

	  // Fill in top edge.
	  while(column < clippedTransitionColumn1) {
	    AccumulatorType elementValue = innerProduct2D<
	      AccumulatorType, AccumulatorType, KernelType, SignalType>(
		kernel, signal, kernelStartRow, kernel.rows(),
		0, kernel.columns(), 0, column - transitionColumn0);
	    elementValue += innerProduct2D<
	      AccumulatorType, AccumulatorType, KernelType, SignalType>(
		kernelUD, signal, kernel.rows() - kernelStartRow,
		kernel.rows(), 0, kernel.columns(), 0,
		column - transitionColumn0);
	    result(outputRow, outputColumn) =
	      static_cast<OutputType>(elementValue);
	    ++column;
	    ++outputColumn;
	  }

	  // Fill in top right corner.
	  while(column < stopColumn) {
	    AccumulatorType elementValue = correlate2DRightReflection<
	      AccumulatorType, AccumulatorType, KernelType, SignalType>(
		kernel, signal, kernelStartRow, kernel.rows(), 0,
		signal.columns() - column + kColOverTwo);
	    elementValue += correlate2DRightReflection<
	      AccumulatorType, AccumulatorType, KernelType, SignalType>(
	      kernelUD, signal,
	      kernel.rows() - kernelStartRow, kernel.rows(), 0,
	      signal.columns() - column + kColOverTwo);
	    result(outputRow, outputColumn) =
	      static_cast<OutputType>(elementValue);
	    ++column;
	    ++outputColumn;
	  }
	  ++row;
	  ++outputRow;
	}

	// Fill in areas along left and right edges of image.
	while(row < clippedTransitionRow1) {

	  // Fill in left edge.
	  int column = startColumn;
	  int outputColumn = 0;
	  while(column < clippedTransitionColumn0) {
	    result(outputRow, outputColumn) = correlate2DLeftReflection<
	      OutputType, AccumulatorType, KernelType, SignalType>(
		kernel, signal, 0, kernel.rows(), row - transitionRow0,
		column + kColOverTwo + 1);
	    ++column;
	    ++outputColumn;
	  }

	  // Fill in right edge.
	  column = clippedTransitionColumn1;
	  outputColumn +=
	    clippedTransitionColumn1 - clippedTransitionColumn0;
	  while(column < stopColumn) {
	    result(outputRow, outputColumn) = correlate2DRightReflection<
	      OutputType, AccumulatorType, KernelType, SignalType>(
		kernel, signal, 0, kernel.rows(), row - transitionRow0,
		signal.columns() - column + kColOverTwo);
	    ++column;
	    ++outputColumn;
	  }
	  ++row;
	  ++outputRow;
	}

	// Fill in areas along bottom of image.
        while(row < stopRow) {
	  // Intermediate variable keeping track of the last kernel
	  // row that overlaps valid image data.
	  size_t kernelStopRow = signal.rows() - row + kRowOverTwo;
	  size_t kernelStopRowUD = kernel.rows() - kernelStopRow;

	  // Fill in bottom left corner.
	  int column = startColumn;
	  int outputColumn = 0;
	  while(column < clippedTransitionColumn0) {
	    AccumulatorType elementValue = correlate2DLeftReflection<
	      AccumulatorType, AccumulatorType, KernelType, SignalType>(
		kernel, signal, 0, kernelStopRow,
		signal.rows() - kernelStopRow,
		column + kColOverTwo + 1);
	    elementValue += correlate2DLeftReflection<
	      AccumulatorType, AccumulatorType, KernelType, SignalType>(
		kernelUD, signal,
		0, kernelStopRowUD, signal.rows() - kernelStopRowUD,
		column + kColOverTwo + 1);
	    result(outputRow, outputColumn) =
	      static_cast<OutputType>(elementValue);
	    ++column;
	    ++outputColumn;
	  }

	  // Fill in bottom edge.
	  while(column < clippedTransitionColumn1) {
	    AccumulatorType elementValue = innerProduct2D<
	      AccumulatorType, AccumulatorType, KernelType, SignalType>(
		kernel, signal, 0, kernelStopRow,
		0, kernel.columns(), signal.rows() - kernelStopRow,
		column - transitionColumn0);
	    elementValue += innerProduct2D<
	      AccumulatorType, AccumulatorType, KernelType, SignalType>(
		kernelUD, signal, 0, kernelStopRowUD,
		0, kernel.columns(), signal.rows() - kernelStopRowUD,
		column - transitionColumn0);
	    result(outputRow, outputColumn) =
	      static_cast<OutputType>(elementValue);
	    ++column;
	    ++outputColumn;
	  }

	  // Fill in bottom right corner.
	  while(column < stopColumn) {
	    AccumulatorType elementValue = correlate2DRightReflection<
	      AccumulatorType, AccumulatorType, KernelType, SignalType>(
		kernel, signal, 0, kernelStopRow,
		signal.rows() - kernelStopRow,
		signal.columns() + kColOverTwo - column);
	    elementValue += correlate2DRightReflection<
	      AccumulatorType, AccumulatorType, KernelType, SignalType>(
		kernelUD, signal, 0, kernelStopRowUD,
		signal.rows() - kernelStopRowUD,
		signal.columns() + kColOverTwo - column);
	    result(outputRow, outputColumn) =
	      static_cast<OutputType>(elementValue);
	    ++column;
	    ++outputColumn;
	  }
	  ++row;
	  ++outputRow;
	}

	// Fill in the middle of the image.
	Index2D newCorner0(kRowOverTwo, kColOverTwo);
	Index2D newCorner1(static_cast<int>(signal.rows()) - kRowOverTwo,
			   static_cast<int>(signal.columns()) - kColOverTwo);
	Index2D resultCorner0(static_cast<int>(kernel.rows()) - 1,
			      static_cast<int>(kernel.columns()) - 1);
        correlate2DCommon<OutputType, AccumulatorType, KernelType, SignalType>(
	  kernel, signal, result, newCorner0, newCorner1, resultCorner0);
        return result;
      }


      template <class OutputType, class AccumulatorType,
		class KernelType, class SignalType>
      Array2D<OutputType>
      correlate2DWrapSignal(const Array2D<KernelType>& kernel,
			    const Array2D<SignalType>& signal,
			    const Index2D& corner0,
			    const Index2D& corner1)
      {
	// Brace yourself...
	
	size_t outputRows =
	  corner1.getRow() - corner0.getRow();
	size_t outputColumns =
	  corner1.getColumn() - corner0.getColumn();
        Array2D<OutputType> result(outputRows, outputColumns);

	// Constants specified without reference to signal or result.
	const int kRowOverTwo = static_cast<int>(kernel.rows()) / 2;
	const int kColOverTwo = static_cast<int>(kernel.columns()) / 2;
	
	// Constants specified with respect to rows & columns in
	// argument signal.
	const int startRow = corner0.getRow();
	const int stopRow = corner1.getRow();
	const int transitionRow0 = kRowOverTwo;
	const int transitionRow1 =
	  static_cast<int>(signal.rows()) - transitionRow0;
	const int startColumn = corner0.getColumn();
	const int stopColumn = corner1.getColumn();
	const int transitionColumn0 = kColOverTwo;
	const int transitionColumn1 =
	  static_cast<int>(signal.columns()) - transitionColumn0;
	const int clippedTransitionRow0 =
	  brick::common::clip(transitionRow0, startRow, stopRow);
	const int clippedTransitionRow1 =
	  brick::common::clip(transitionRow1, startRow, stopRow);
	const int clippedTransitionColumn0 =
	  brick::common::clip(transitionColumn0, startColumn, stopColumn);
	const int clippedTransitionColumn1 =
	  brick::common::clip(transitionColumn1, startColumn, stopColumn);

	// Fill in areas along top of image.
	int row = startRow;
	int outputRow = 0;
	Array1D<AccumulatorType> accumulator(result.columns());
        while(row < clippedTransitionRow0) {

	  // Intermediate variable to simplify things later.
	  size_t kernelStartRow = transitionRow0 - row;

	  // Fill in top left corner.  Do this in two passes to reduce
	  // cache misses.
	  int outputColumn = 0;
	  for(int column = startColumn; column < clippedTransitionColumn0;
	      ++column) {
	    // Add upper left and upper right contributions.
	    accumulator[outputColumn] = correlate2DLeftRightWrap<
	      AccumulatorType, AccumulatorType, KernelType, SignalType>(
		kernel, signal, kernelStartRow, kernel.rows(), 0,
		column + kColOverTwo + 1);
	    ++outputColumn;
	  }
	  outputColumn = 0;
	  for(int column = startColumn; column < clippedTransitionColumn0;
	      ++column) {
	    // Add lower left and lower right contributions.
	    accumulator[outputColumn] += correlate2DLeftRightWrap<
	      AccumulatorType, AccumulatorType, KernelType, SignalType>(
		kernel, signal, 0, kernelStartRow,
		signal.rows() - kernelStartRow, column + kColOverTwo + 1);
	    result(outputRow, outputColumn) =
	      static_cast<OutputType>(accumulator[outputColumn]);
	    ++outputColumn;
	  }

	  // Fill in top edge.  Do this in two passes to reduce
	  // cache misses.
	  outputColumn = clippedTransitionColumn0 - startColumn;
	  for(int column = clippedTransitionColumn0;
	      column < clippedTransitionColumn1;
	      ++column) {
	    // Add the upper contribution.
	    accumulator[outputColumn] = innerProduct2D<
	      AccumulatorType, AccumulatorType, KernelType, SignalType>(
		kernel, signal, kernelStartRow, kernel.rows(),
		0, kernel.columns(), 0, column - transitionColumn0);
	    ++outputColumn;
	  }
	  outputColumn = clippedTransitionColumn0 - startColumn;
	  for(int column = clippedTransitionColumn0;
	      column < clippedTransitionColumn1;
	      ++column) {
	    // Add the lower contribution.
	    accumulator[outputColumn] += innerProduct2D<
	      AccumulatorType, AccumulatorType, KernelType, SignalType>(
		kernel, signal, 0, kernelStartRow,
		0, kernel.columns(),
		signal.rows() - kernelStartRow, column - transitionColumn0);
	    result(outputRow, outputColumn) =
	      static_cast<OutputType>(accumulator[outputColumn]);
	    ++outputColumn;
	  }

	  // Fill in top right corner.  Because we're wrapping, this
	  // corner has the same element values as a shifted version
	  // of the top left corner.  If we calculated the top left
	  // corner, we simply copy it.
	  outputColumn = clippedTransitionColumn1 - startColumn;
	  for(int column = clippedTransitionColumn1; column < stopColumn;
	      ++column) {
	    int donorColumn =
	      outputColumn - static_cast<int>(signal.columns());
	    if(donorColumn >= 0) {
	      result(outputRow, outputColumn) = result(outputRow, donorColumn);
	    } else {
	      // Add upper left and upper right contributions.
	      accumulator[outputColumn] = correlate2DLeftRightWrap<
		AccumulatorType, AccumulatorType, KernelType, SignalType>(
		  kernel, signal, kernelStartRow, kernel.rows(), 0,
		  column - signal.columns() + kColOverTwo + 1);
	    }
	    ++outputColumn;
	  }
	  outputColumn = clippedTransitionColumn1 - startColumn;
	  for(int column = clippedTransitionColumn1; column < stopColumn;
	      ++column) {
	    int donorColumn = outputColumn - static_cast<int>(signal.columns());
	    if(donorColumn < 0) {
	      // Add lower left and lower right contributions.
	      accumulator[outputColumn] += correlate2DLeftRightWrap<
		AccumulatorType, AccumulatorType, KernelType, SignalType>(
		  kernel, signal, 0, kernelStartRow,
		  signal.rows() - kernelStartRow,
		  column - signal.columns() + kColOverTwo + 1);
	      result(outputRow, outputColumn) =
		static_cast<OutputType>(accumulator[outputColumn]);
	    }
	    ++outputColumn;
	  }
	  ++row;
	  ++outputRow;
	}

	// Fill in areas along left and right edges of image.
        while(row < clippedTransitionRow1) {

	  // Fill in left edge.
	  size_t outputColumn = 0;
	  for(int column = startColumn; column < clippedTransitionColumn0;
	      ++column) {
	    result(outputRow, outputColumn) = correlate2DLeftRightWrap<
	      OutputType, AccumulatorType, KernelType, SignalType>(
		kernel, signal, 0, kernel.rows(), row - transitionRow0,
		column + kColOverTwo + 1);
	    ++outputColumn;
	  }

	  // Fill in right edge.
	  outputColumn = clippedTransitionColumn1 - startColumn;
	  for(int column = clippedTransitionColumn1; column < stopColumn;
	      ++column) {
	    int donorColumn = (static_cast<int>(outputColumn)
			       - static_cast<int>(signal.columns()));
	    if(donorColumn >= 0) {
	      result(outputRow, outputColumn) = result(outputRow, donorColumn);
	    } else {
	      result(outputRow, outputColumn) = correlate2DLeftRightWrap<
		OutputType, AccumulatorType, KernelType, SignalType>(
		  kernel, signal, 0, kernel.rows(), row - transitionRow0,
		  column - signal.columns() + kColOverTwo + 1);
	    }
	    ++outputColumn;
	  }
	  ++row;
	  ++outputRow;
	}

	// Fill in areas along bottom of image.
        while(row < stopRow) {
	  int donorRow = row - outputRow - static_cast<int>(signal.rows());

	  if(donorRow >= 0) {
	    std::copy(result.rowBegin(donorRow), result.rowEnd(donorRow),
		      result.rowBegin(outputRow));
	  } else {

	    // Intermediate variable to simplify things later.
	    size_t kernelStartRow = transitionRow0 - row + signal.rows();

	    // Fill in bottom left corner.  Do this in two passes to reduce
	    // cache misses.
	    int outputColumn = 0;
	    for(int column = startColumn; column < clippedTransitionColumn0;
		++column) {
	      // Add upper left and upper right contributions.
	      accumulator[outputColumn] = correlate2DLeftRightWrap<
		AccumulatorType, AccumulatorType, KernelType, SignalType>(
		  kernel, signal, kernelStartRow, kernel.rows(), 0,
		  column + kColOverTwo + 1);
	      ++outputColumn;
	    }
	    outputColumn = 0;
	    for(int column = startColumn; column < clippedTransitionColumn0;
		++column) {
	      // Add lower left and lower right contributions.
	      accumulator[outputColumn] += correlate2DLeftRightWrap<
		AccumulatorType, AccumulatorType, KernelType, SignalType>(
		  kernel, signal, 0, kernelStartRow,
		  signal.rows() - kernelStartRow, column + kColOverTwo + 1);
	      result(outputRow, outputColumn) =
		static_cast<OutputType>(accumulator[outputColumn]);
	      ++outputColumn;
	    }

	    // Fill in bottom edge.  Do this in two passes to reduce
	    // cache misses.
	    outputColumn = clippedTransitionColumn0 - startColumn;
	    for(int column = clippedTransitionColumn0;
		column < clippedTransitionColumn1;
		++column) {
	      // Add the upper contribution.
	      accumulator[outputColumn] = innerProduct2D<
		AccumulatorType, AccumulatorType, KernelType, SignalType>(
		  kernel, signal, kernelStartRow, kernel.rows(),
		  0, kernel.columns(), 0, column - transitionColumn0);
	      ++outputColumn;
	    }
	    outputColumn = clippedTransitionColumn0 - startColumn;
	    for(int column = clippedTransitionColumn0;
		column < clippedTransitionColumn1;
		++column) {
	      // Add the lower contribution.
	      accumulator[outputColumn] += innerProduct2D<
		AccumulatorType, AccumulatorType, KernelType, SignalType>(
		  kernel, signal, 0, kernelStartRow,
		  0, kernel.columns(),
		  signal.rows() - kernelStartRow, column - transitionColumn0);
	      result(outputRow, outputColumn) =
		static_cast<OutputType>(accumulator[outputColumn]);
	      ++outputColumn;
	    }

	    // Fill in bottom right corner.  Because we're wrapping, this
	    // corner has the same element values as a shifted version
	    // of the bottom left corner.  If we calculated the bottom left
	    // corner, we simply copy it.
	    outputColumn = clippedTransitionColumn1 - startColumn;
	    for(int column = clippedTransitionColumn1; column < stopColumn;
		++column) {
	      int donorColumn = outputColumn - static_cast<int>(signal.columns());
	      if(donorColumn >= 0) {
		result(outputRow, outputColumn) =
		  result(outputRow, donorColumn);
	      } else {
		// Add upper left and upper right contributions.
		accumulator[outputColumn] = correlate2DLeftRightWrap<
		  AccumulatorType, AccumulatorType, KernelType, SignalType>(
		    kernel, signal, kernelStartRow, kernel.rows(), 0,
		    column - signal.columns() + kColOverTwo + 1);
	      }
	      ++outputColumn;
	    }
	    outputColumn = clippedTransitionColumn1 - startColumn;
	    for(int column = clippedTransitionColumn1; column < stopColumn;
		++column) {
	      int donorColumn = outputColumn - static_cast<int>(signal.columns());
	      if(donorColumn < 0) {
		// Add lower left and lower right contributions.
		accumulator[outputColumn] += correlate2DLeftRightWrap<
		  AccumulatorType, AccumulatorType, KernelType, SignalType>(
		    kernel, signal, 0, kernelStartRow,
		    signal.rows() - kernelStartRow,
		    column - signal.columns() + kColOverTwo + 1);
		result(outputRow, outputColumn) =
		  static_cast<OutputType>(accumulator[outputColumn]);
	      }
	      ++outputColumn;
	    }
	  }
	  ++row;
	  ++outputRow;
	}

	// Fill in the middle of the image.
	Index2D newCorner0(kRowOverTwo, kColOverTwo);
	Index2D newCorner1(static_cast<int>(signal.rows()) - kRowOverTwo,
			   static_cast<int>(signal.columns()) - kColOverTwo);
	Index2D resultCorner0(static_cast<int>(kernel.rows()) - 1,
			      static_cast<int>(kernel.columns()) - 1);
        correlate2DCommon<OutputType, AccumulatorType, KernelType, SignalType>(
	  kernel, signal, result, newCorner0, newCorner1, resultCorner0);
        return result;
      }

      
      template <class KernelType>
      Array2D<KernelType>
      reverseKernel(const Array2D<KernelType>& kernel) {
	Array2D<KernelType> reversedKernel(kernel.rows(), kernel.columns());
	for(size_t rowIndex = 0; rowIndex < kernel.rows(); ++rowIndex) {
	  std::reverse_copy(
	    kernel.rowBegin(rowIndex), kernel.rowEnd(rowIndex),
	    reversedKernel.rowBegin(kernel.rows() - rowIndex - 1));
	}
	return reversedKernel;
      }
      
                   
    } // namespace privateCode
    /// @endcond


    template <class OutputType, class AccumulatorType,
	      class KernelType, class SignalType>
    inline Array2D<OutputType>
    convolve2D(const Array2D<KernelType>& kernel,
	       const Array2D<SignalType>& signal,
	       ConvolutionStrategy strategy,
	       ConvolutionROI roi)
    {
      Array2D<KernelType> reversedKernel = privateCode::reverseKernel(kernel);
      return correlate2D<OutputType, AccumulatorType, KernelType, SignalType>(
	reversedKernel, signal, strategy, roi);
    }
    

    template <class OutputType, class AccumulatorType,
	      class KernelType, class SignalType, class FillType>
    inline Array2D<OutputType>
    convolve2D(const Array2D<KernelType>& kernel,
	       const Array2D<SignalType>& signal,
	       ConvolutionStrategy strategy,
	       ConvolutionROI roi,
	       const FillType& fillValue)
    {
      Array2D<KernelType> reversedKernel = privateCode::reverseKernel(kernel);
      return correlate2D<OutputType, AccumulatorType, KernelType, SignalType>(
	reversedKernel, signal, strategy, roi, fillValue);
    }
    

    template <class OutputType, class AccumulatorType,
	      class KernelType, class SignalType>
    inline Array2D<OutputType>
    convolve2D(const Array2D<KernelType>& kernel,
	       const Array2D<SignalType>& signal,
	       ConvolutionStrategy strategy,
	       const Index2D& corner0,
	       const Index2D& corner1)
    {
      Array2D<KernelType> reversedKernel = privateCode::reverseKernel(kernel);
      return correlate2D<OutputType, AccumulatorType, KernelType, SignalType>(
	reversedKernel, signal, strategy, corner0, corner1);
    }
    

    template <class OutputType, class AccumulatorType,
	      class KernelType, class SignalType, class FillType>
    inline Array2D<OutputType>
    convolve2D(const Array2D<KernelType>& kernel,
	       const Array2D<SignalType>& signal,
	       ConvolutionStrategy strategy,
	       const Index2D& corner0,
	       const Index2D& corner1,
	       const FillType& fillValue)
    {
      Array2D<KernelType> reversedKernel = privateCode::reverseKernel(kernel);
      return correlate2D<OutputType, AccumulatorType, KernelType, SignalType>(
	reversedKernel, signal, strategy, corner0, corner1, fillValue);
    }
    

    template <class OutputType, class AccumulatorType,
	      class KernelType, class SignalType>
    Array2D<OutputType>
    correlate2D(const Array2D<KernelType>& kernel,
		const Array2D<SignalType>& signal,
		ConvolutionStrategy strategy,
		ConvolutionROI roi)
    {
      switch(roi) {
      case BRICK_CONVOLVE_ROI_SAME:
      {
	Index2D corner0(0, 0);
	Index2D corner1(static_cast<int>(signal.rows()), static_cast<int>(signal.columns()));
	return correlate2D<OutputType, AccumulatorType, KernelType, SignalType>(
	  kernel, signal, strategy, corner0, corner1);
        break;
      }
      case BRICK_CONVOLVE_ROI_VALID:
      {
	Index2D corner0(static_cast<int>(kernel.rows()) / 2, static_cast<int>(kernel.columns()) / 2);
	Index2D corner1(static_cast<int>(signal.rows()) - static_cast<int>(corner0.getRow()),
			static_cast<int>(signal.columns()) - static_cast<int>(corner0.getColumn()));
	return correlate2D<OutputType, AccumulatorType, KernelType, SignalType>(
	  kernel, signal, strategy, corner0, corner1);
        break;
      }
      case BRICK_CONVOLVE_ROI_FULL:
      {
	Index2D corner0(-(static_cast<int>(kernel.rows()) / 2),
			-(static_cast<int>(kernel.columns()) / 2));
	Index2D corner1(static_cast<int>(signal.rows()) - static_cast<int>(corner0.getRow()),
			static_cast<int>(signal.columns()) - static_cast<int>(corner0.getColumn()));
	return correlate2D<OutputType, AccumulatorType, KernelType, SignalType>(
	  kernel, signal, strategy, corner0, corner1);
        break;
      }
      default:
        BRICK_THROW(brick::common::LogicException, "correlate2D()",
                  "Illegal value for roi argument.");
        break;
      }
      return Array2D<OutputType>();
    }
    

    template <class OutputType, class AccumulatorType,
	      class KernelType, class SignalType, class FillType>
    Array2D<OutputType>
    correlate2D(const Array2D<KernelType>& kernel,
		const Array2D<SignalType>& signal,
		ConvolutionStrategy strategy,
		ConvolutionROI roi,
		const FillType& fillValue)
    {
      switch(roi) {
      case BRICK_CONVOLVE_ROI_SAME:
      {
	Index2D corner0(0, 0);
	Index2D corner1(static_cast<int>(signal.rows()), static_cast<int>(signal.columns()));
	return correlate2D<OutputType, AccumulatorType, KernelType, SignalType, FillType>(
	  kernel, signal, strategy, corner0, corner1, fillValue);
        break;
      }
      case BRICK_CONVOLVE_ROI_VALID:
      {
	Index2D corner0(static_cast<int>(kernel.rows()) / 2, static_cast<int>(kernel.columns()) / 2);
	Index2D corner1(static_cast<int>(signal.rows()) - corner0.getRow(),
			static_cast<int>(signal.columns()) - corner0.getColumn());
	return correlate2D<OutputType, AccumulatorType, KernelType, SignalType, FillType>(
	  kernel, signal, strategy, corner0, corner1, fillValue);
        break;
      }
      case BRICK_CONVOLVE_ROI_FULL:
      {
	Index2D corner0(-(static_cast<int>(kernel.rows()) / 2),
			-(static_cast<int>(kernel.columns()) / 2));
	Index2D corner1(static_cast<int>(signal.rows()) - corner0.getRow(),
			static_cast<int>(signal.columns()) - corner0.getColumn());
	return correlate2D<OutputType, AccumulatorType, KernelType, SignalType, FillType>(
	  kernel, signal, strategy, corner0, corner1, fillValue);
        break;
      }
      default:
        BRICK_THROW(brick::common::LogicException, "correlate2D()",
                  "Illegal value for roi argument.");
        break;
      }
      return Array2D<OutputType>();
    }

    
    template <class OutputType, class AccumulatorType,
	      class KernelType, class SignalType>
    Array2D<OutputType>
    correlate2D(const Array2D<KernelType>& kernel,
		const Array2D<SignalType>& signal,
		ConvolutionStrategy strategy,
		const Index2D& corner0,
		const Index2D& corner1)
    {
      if(kernel.rows() % 2 != 1) {
        BRICK_THROW(brick::common::ValueException, "correlate2D()",
                  "Argument kernel must have an odd number of rows.");
      }
      if(kernel.columns() % 2 != 1) {
        BRICK_THROW(brick::common::ValueException, "correlate2D()",
                  "Argument kernel must have an odd number of columns.");
      }
      // Note(xxx): is the following check necessary?
      if(kernel.rows() > signal.rows()) {
        BRICK_THROW(brick::common::ValueException, "correlate2D()",
                  "Argument kernel must not have more rows than "
                  "argument signal.");
      }
      if(kernel.columns() > signal.columns()) {
        BRICK_THROW(brick::common::ValueException, "correlate2D()",
                  "Argument kernel must not have more columns than "
                  "argument signal.");
      }

      switch(strategy) {
      case BRICK_CONVOLVE_TRUNCATE_RESULT:
        return privateCode::correlate2DTruncateResult<
	  OutputType, AccumulatorType, KernelType, SignalType>(kernel, signal);
        break;
      case BRICK_CONVOLVE_PAD_RESULT:
      case BRICK_CONVOLVE_PAD_SIGNAL:
        BRICK_THROW(brick::common::ValueException, "correlate2D()",
                  "The specified convolution strategy requires that a "
		  "fill value be specified.");
      case BRICK_CONVOLVE_ZERO_PAD_SIGNAL:
        return privateCode::correlate2DZeroPadSignal<
	  OutputType, AccumulatorType, KernelType, SignalType>(
	    kernel, signal, corner0, corner1);
        break;
      case BRICK_CONVOLVE_REFLECT_SIGNAL:
        return privateCode::correlate2DReflectSignal<
	  OutputType, AccumulatorType, KernelType, SignalType>(
	    kernel, signal, corner0, corner1);
        break;
      case BRICK_CONVOLVE_WRAP_SIGNAL:
        return privateCode::correlate2DWrapSignal<
	  OutputType, AccumulatorType, KernelType, SignalType>(
	    kernel, signal, corner0, corner1);
        break;
      default:
        BRICK_THROW(brick::common::LogicException, "correlate2D()",
                  "Illegal value for strategy argument.");
        break;
      }
      return Array2D<OutputType>();
    }
    

    template <class OutputType, class AccumulatorType,
	      class KernelType, class SignalType, class FillType>
    Array2D<OutputType>
    correlate2D(const Array2D<KernelType>& kernel,
		const Array2D<SignalType>& signal,
		ConvolutionStrategy strategy,
		const Index2D& corner0,
		const Index2D& corner1,
		const FillType& fillValue)
    {
      if(kernel.rows() % 2 != 1) {
        BRICK_THROW(brick::common::ValueException, "correlate2D()",
                  "Argument kernel must have an odd number of rows.");
      }
      if(kernel.columns() % 2 != 1) {
        BRICK_THROW(brick::common::ValueException, "correlate2D()",
                  "Argument kernel must have an odd number of columns.");
      }
      // Note(xxx): is the following check necessary?
      if(kernel.rows() > signal.rows()) {
        BRICK_THROW(brick::common::ValueException, "correlate2D()",
                  "Argument kernel must not have more rows than "
                  "argument signal.");
      }
      if(kernel.columns() > signal.columns()) {
        BRICK_THROW(brick::common::ValueException, "correlate2D()",
                  "Argument kernel must not have more columns than "
                  "argument signal.");
      }
    
      switch(strategy) {
      case BRICK_CONVOLVE_TRUNCATE_RESULT:
        return privateCode::correlate2DTruncateResult<
	  OutputType, AccumulatorType, KernelType, SignalType>(kernel, signal);
        break;
      case BRICK_CONVOLVE_PAD_RESULT:
        return privateCode::correlate2DPadResult<
	  OutputType, AccumulatorType, KernelType, SignalType>(
	    kernel, signal, corner0, corner1,
	    static_cast<OutputType>(fillValue));
        break;
      case BRICK_CONVOLVE_PAD_SIGNAL:
        return privateCode::correlate2DPadSignal<
	  OutputType, AccumulatorType, KernelType, SignalType>(
	    kernel, signal, corner0, corner1,
	    static_cast<SignalType>(fillValue));
        break;
      case BRICK_CONVOLVE_ZERO_PAD_SIGNAL:
        return privateCode::correlate2DZeroPadSignal<
	  OutputType, AccumulatorType, KernelType, SignalType>(
	    kernel, signal, corner0, corner1);
        break;
      case BRICK_CONVOLVE_REFLECT_SIGNAL:
        return privateCode::correlate2DReflectSignal<
	  OutputType, AccumulatorType, KernelType, SignalType>(
	    kernel, signal, corner0, corner1);
        break;
      case BRICK_CONVOLVE_WRAP_SIGNAL:
        return privateCode::correlate2DWrapSignal<
	  OutputType, AccumulatorType, KernelType, SignalType>(
	    kernel, signal, corner0, corner1);
        break;
      default:
        BRICK_THROW(brick::common::LogicException, "correlate2D()",
                  "Illegal value for strategy argument.");
        break;
      }
      return Array2D<OutputType>();
    }
    
  } // namespace numeric

} // namespace brick

#endif /* #ifndef BRICK_NUMERIC_CONVOLVE2D_IMPL_HH */
