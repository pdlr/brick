/**
***************************************************************************
* @file brick/numeric/convolve1D_impl.hh
*
* Header file defining 1D correlation and convolution functions.
*
* Copyright (C) 2006-2007, 2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_CONVOLVE1D_IMPL_HH
#define BRICK_NUMERIC_CONVOLVE1D_IMPL_HH

// This file is included by convolve1D.hh, and should not be directly included
// by user code, so no need to include convolve1D.hh here.
// 
// #include <brick/numeric/convolve1D.hh>

#include <algorithm> // For std::reverse_copy()
#include <numeric> // For std::partial_sum()

#include <brick/common/functional.hh> // for clip()
#include <brick/common/functional.hh> // for clip()
#include <brick/numeric/numericTraits.hh>

namespace brick {

  namespace numeric {

    /// @cond privateCode    
    namespace privateCode {

      template <class OutputType, class KernelType, class SignalType>
      Array1D<OutputType>
      accumulateKernel(const Array1D<KernelType>& kernel,
		       SignalType fillValue)
      {
	typedef typename brick::numeric::ArithmeticTraits<
          KernelType, KernelType>::SumType KernelSumType;
	Array1D<KernelSumType> accumulatedKernel(kernel.size() + 1);
	accumulatedKernel[0] = static_cast<KernelSumType>(0);
	std::partial_sum(kernel.begin(), kernel.end(),
			 accumulatedKernel.begin() + 1);

	Array1D<OutputType> result(accumulatedKernel.size());
	for(size_t index0 = 0; index0 < result.size(); ++index0) {
	  result[index0] = static_cast<OutputType>(
	    accumulatedKernel[index0]
	    * static_cast<OutputType>(fillValue));
	}
	return result;
      }


      template <class OutputType, class KernelType, class SignalType>
      void
      correlate1DCommon(
	const Array1D<KernelType>& kernel,
	typename Array1D<SignalType>::const_iterator signalBegin,
	typename Array1D<SignalType>::const_iterator signalEnd,
	typename Array1D<OutputType>::iterator resultIterator)
      {
        typedef typename Array1D<SignalType>::const_iterator SignalIterator;
        typedef typename Array1D<KernelType>::const_iterator KernelIterator;

        SignalIterator signalIterator = signalBegin;
        while(signalIterator != signalEnd) {
          // Calculate one dot-product between the kernel and the
          // current section of the signal.
          OutputType dotProduct = static_cast<OutputType>(0);
          SignalIterator signalIteratorCopy = signalIterator;
          KernelIterator kernelIterator = kernel.begin();
          while(kernelIterator != kernel.end()) {
            dotProduct += static_cast<OutputType>(
	      *(kernelIterator++)
	      * static_cast<OutputType>(*(signalIteratorCopy++)));;
          }
          *resultIterator = dotProduct;
          ++resultIterator;
          ++signalIterator;
        }
      }


      template <class OutputType, class KernelType, class SignalType>
      Array1D<OutputType>
      correlate1DTruncateResult(const Array1D<KernelType>& kernel,
				const Array1D<SignalType>& signal)
      {
        Array1D<OutputType> result(signal.size() - kernel.size() + 1);
        correlate1DCommon<OutputType, KernelType, SignalType>(
	  kernel, signal.begin(), signal.end() - kernel.size() + 1,
	  result.begin());
        return result;
      }


      template <class OutputType, class KernelType, class SignalType>
      Array1D<OutputType>
      correlate1DPadResult(const Array1D<KernelType>& kernel,
			   const Array1D<SignalType>& signal,
			   int boundary0,
			   int boundary1,
			   OutputType fillValue)
      {
        Array1D<OutputType> result(boundary1 - boundary0);

	// Check for degenerate case.
	if(kernel.size() > signal.size()) {
	  result = fillValue;
	  return result;
	}

	// Establish significant coordinates in the output image.

	// Constants specified without reference to signal or result.
	const int kSizeOverTwo = static_cast<int>(kernel.size()) / 2;
	
	// Constants specified with respect to argument signal.
	const int transitionIndex0 = kSizeOverTwo;
	const int transitionIndex1 = static_cast<int>(signal.size()) - kSizeOverTwo;
	const int clippedTransitionIndex0 =
	  brick::common::clip(transitionIndex0, boundary0, boundary1);
	const int clippedTransitionIndex1 =
	  brick::common::clip(transitionIndex1, boundary0, boundary1);

	// Constants specified with respect to result.
	const int resultTransitionIndex0 =
	  clippedTransitionIndex0 - boundary0;

	int inputIndex = boundary0;
	int outputIndex = 0;
	while(inputIndex < clippedTransitionIndex0) {
	  result[outputIndex] = fillValue;
	  ++inputIndex;
	  ++outputIndex;
	}
	inputIndex = clippedTransitionIndex1;
	outputIndex += clippedTransitionIndex1 - clippedTransitionIndex0;
	while(inputIndex < boundary1) {
	  result[outputIndex] = fillValue;
	  ++inputIndex;
	  ++outputIndex;
	}

        correlate1DCommon<OutputType, KernelType, SignalType>(
	  kernel,
	  signal.begin() + clippedTransitionIndex0 - kSizeOverTwo,
	  signal.begin() + clippedTransitionIndex1 - kSizeOverTwo,
	  result.begin() + resultTransitionIndex0);
        return result;
      }
    
                   
      template <class OutputType, class KernelType, class SignalType>
      Array1D<OutputType>
      correlate1DZeroPadSignal(const Array1D<KernelType>& kernel,
			       const Array1D<SignalType>& signal,
			       int boundary0,
			       int boundary1)
      {
        Array1D<OutputType> result(boundary1 - boundary0);

	// Constants specified without reference to signal or result.
	const int kSizeOverTwo = static_cast<int>(kernel.size()) / 2;
	
	// Constants specified with respect to argument signal.
	const int transitionIndex0 = kSizeOverTwo;
	const int transitionIndex1 = static_cast<int>(signal.size()) - kSizeOverTwo;
	const int clippedTransitionIndex0 =
	  brick::common::clip(transitionIndex0, boundary0, boundary1);
	const int clippedTransitionIndex1 =
	  brick::common::clip(transitionIndex1, boundary0, boundary1);

	// Constants specified with respect to result.
	const int resultTransitionIndex0 =
	  clippedTransitionIndex0 - boundary0;

	int inputIndex = boundary0;
	int outputIndex = 0;
	while(inputIndex < clippedTransitionIndex0) {
          OutputType dotProduct0 = static_cast<OutputType>(0);
          int kernelIndex = clippedTransitionIndex0 - inputIndex;
          size_t signalIndex = boundary0 + kSizeOverTwo;
          while(kernelIndex < static_cast<int>(kernel.size())) {
            dotProduct0 += static_cast<OutputType>(
	      kernel[kernelIndex]
	      * static_cast<OutputType>(signal[signalIndex]));
            ++kernelIndex;
            ++signalIndex;
          }
          result[outputIndex] = dotProduct0;
	  ++inputIndex;
	  ++outputIndex;
        }

	inputIndex = clippedTransitionIndex1;
	outputIndex += clippedTransitionIndex1 - clippedTransitionIndex0;
	while(inputIndex < boundary1) {
          OutputType dotProduct0 = static_cast<OutputType>(0);
          int kernelStopIndex = boundary1 - inputIndex;
          int kernelIndex = 0;
          size_t signalIndex = inputIndex - kSizeOverTwo;
          while(kernelIndex < kernelStopIndex) {
            dotProduct0 += static_cast<OutputType>(
	      kernel[kernelIndex]
	      * static_cast<OutputType>(signal[signalIndex]));
            ++kernelIndex;
            ++signalIndex;
          }
          result[outputIndex] = dotProduct0;
	  ++inputIndex;
	  ++outputIndex;
        }
	
        correlate1DCommon<OutputType, KernelType, SignalType>(
	  kernel,
	  signal.begin() + clippedTransitionIndex0 - kSizeOverTwo,
	  signal.begin() + clippedTransitionIndex1 - kSizeOverTwo,
	  result.begin() + resultTransitionIndex0);
        return result;
      }


      template <class OutputType, class KernelType, class SignalType>
      Array1D<OutputType>
      correlate1DPadSignal(const Array1D<KernelType>& kernel,
			   const Array1D<SignalType>& signal,
			   int boundary0,
			   int boundary1,
			   SignalType fillValue)
      {
        Array1D<OutputType> result(boundary1 - boundary0);

	// Precompute numbers which will make it easy to total up the
	// portions of the kernel which overlap padded areas of the
	// signal.
	Array1D<OutputType> accumulatedKernel = accumulateKernel<
	  OutputType, KernelType, SignalType>(kernel, fillValue);
	
	// Constants specified without reference to signal or result.
	const int kSizeOverTwo = static_cast<int>(kernel.size()) / 2;
	
	// Constants specified with respect to argument signal.
	const int transitionIndex0 = kSizeOverTwo;
	const int transitionIndex1 = static_cast<int>(signal.size()) - kSizeOverTwo;
	const int clippedTransitionIndex0 =
	  brick::common::clip(transitionIndex0, boundary0, boundary1);
	const int clippedTransitionIndex1 =
	  brick::common::clip(transitionIndex1, boundary0, boundary1);

	// Constants specified with respect to result.
	const int resultTransitionIndex0 =
	  clippedTransitionIndex0 - boundary0;

	int inputIndex = boundary0;
	int outputIndex = 0;
	while(inputIndex < clippedTransitionIndex0) {
          int kernelIndex = clippedTransitionIndex0 - inputIndex;
          size_t signalIndex = boundary0 + kSizeOverTwo;
          OutputType dotProduct0 = accumulatedKernel[kernelIndex];
          while(kernelIndex < static_cast<int>(kernel.size())) {
            dotProduct0 += static_cast<OutputType>(
	      kernel[kernelIndex]
	      * static_cast<OutputType>(signal[signalIndex]));
            ++kernelIndex;
            ++signalIndex;
          }
          result[outputIndex] = dotProduct0;
	  ++inputIndex;
	  ++outputIndex;
        }

	inputIndex = clippedTransitionIndex1;
	outputIndex += clippedTransitionIndex1 - clippedTransitionIndex0;
	while(inputIndex < boundary1) {
          int kernelStopIndex = boundary1 - inputIndex;
          int kernelIndex = 0;
          size_t signalIndex = inputIndex - kSizeOverTwo;
          OutputType dotProduct0 = (accumulatedKernel[kernel.size()]
				    - accumulatedKernel[kernelStopIndex]);
          while(kernelIndex < kernelStopIndex) {
            dotProduct0 += static_cast<OutputType>(
	      kernel[kernelIndex]
	      * static_cast<OutputType>(signal[signalIndex]));
            ++kernelIndex;
            ++signalIndex;
          }
          result[outputIndex] = dotProduct0;
	  ++inputIndex;
	  ++outputIndex;
        }
	
        correlate1DCommon<OutputType, KernelType, SignalType>(
	  kernel,
	  signal.begin() + clippedTransitionIndex0 - kSizeOverTwo,
	  signal.begin() + clippedTransitionIndex1 - kSizeOverTwo,
	  result.begin() + resultTransitionIndex0);
        return result;
      }


      template <class OutputType, class KernelType, class SignalType>
      Array1D<OutputType>
      correlate1DReflectSignal(const Array1D<KernelType>& kernel,
			       const Array1D<SignalType>& signal,
			       int boundary0,
			       int boundary1)
      {
        Array1D<OutputType> result(boundary1 - boundary0);

	// Constants specified without reference to signal or result.
	const int kSizeOverTwo = static_cast<int>(kernel.size()) / 2;
	
	// Constants specified with respect to argument signal.
	const int transitionIndex0 = kSizeOverTwo;
	const int transitionIndex1 = static_cast<int>(signal.size()) - kSizeOverTwo;
	const int clippedTransitionIndex0 =
	  brick::common::clip(transitionIndex0, boundary0, boundary1);
	const int clippedTransitionIndex1 =
	  brick::common::clip(transitionIndex1, boundary0, boundary1);

	// Constants specified with respect to result.
	const int resultTransitionIndex0 =
	  clippedTransitionIndex0 - boundary0;

	int inputIndex = boundary0;
	int outputIndex = 0;
	while(inputIndex < clippedTransitionIndex0) {
          OutputType dotProduct0 = static_cast<OutputType>(0);

          int kernelIndex = clippedTransitionIndex0 - inputIndex;
          size_t signalIndex = boundary0 + kSizeOverTwo;
          while(kernelIndex < static_cast<int>(kernel.size())) {
            dotProduct0 += static_cast<OutputType>(
	      kernel[kernelIndex]
	      * static_cast<OutputType>(signal[signalIndex]));
            ++kernelIndex;
            ++signalIndex;
	  }
          kernelIndex = clippedTransitionIndex0 - inputIndex - 1;
          signalIndex = boundary0 + kSizeOverTwo;
          while(kernelIndex >= 0) {
            dotProduct0 += static_cast<OutputType>(
	      kernel[kernelIndex]
	      * static_cast<OutputType>(signal[signalIndex]));
            --kernelIndex;
            ++signalIndex;
	  }
          result[outputIndex] = dotProduct0;
	  ++inputIndex;
	  ++outputIndex;
	}


	inputIndex = clippedTransitionIndex1;
	outputIndex += clippedTransitionIndex1 - clippedTransitionIndex0;
	while(inputIndex < boundary1) {
          OutputType dotProduct0 = static_cast<OutputType>(0);

          int kernelStopIndex = boundary1 - inputIndex;
          int kernelIndex = 0;
          size_t signalIndex = inputIndex - kSizeOverTwo;
          while(kernelIndex < kernelStopIndex) {
            dotProduct0 += static_cast<OutputType>(
	      kernel[kernelIndex]
	      * static_cast<OutputType>(signal[signalIndex]));
            ++kernelIndex;
            ++signalIndex;
          }
	  --signalIndex;
          while(kernelIndex < static_cast<int>(kernel.size())) {
            dotProduct0 += static_cast<OutputType>(
	      kernel[kernelIndex]
	      * static_cast<OutputType>(signal[signalIndex]));
            ++kernelIndex;
            --signalIndex;
          }
          result[outputIndex] = dotProduct0;
	  ++inputIndex;
	  ++outputIndex;
        }
	
        correlate1DCommon<OutputType, KernelType, SignalType>(
	  kernel,
	  signal.begin() + clippedTransitionIndex0 - kSizeOverTwo,
	  signal.begin() + clippedTransitionIndex1 - kSizeOverTwo,
	  result.begin() + resultTransitionIndex0);
	return result;
      }


      template <class OutputType, class KernelType, class SignalType>
      Array1D<OutputType>
      correlate1DWrapSignal(const Array1D<KernelType>& kernel,
			    const Array1D<SignalType>& signal,
			    int boundary0,
			    int boundary1)
      {
	Array1D<OutputType> result(boundary1 - boundary0);

	// Constants specified without reference to signal or result.
	const int kSizeOverTwo = static_cast<int>(kernel.size()) / 2;
	
	// Constants specified with respect to argument signal.
	const int transitionIndex0 = kSizeOverTwo;
	const int transitionIndex1 = static_cast<int>(signal.size()) - kSizeOverTwo;
	const int clippedTransitionIndex0 =
	  brick::common::clip(transitionIndex0, boundary0, boundary1);
	const int clippedTransitionIndex1 =
	  brick::common::clip(transitionIndex1, boundary0, boundary1);

	// Constants specified with respect to result.
	const int resultTransitionIndex0 =
	  clippedTransitionIndex0 - boundary0;

	int inputIndex = boundary0;
	int outputIndex = 0;
	while(inputIndex < clippedTransitionIndex0) {
          OutputType dotProduct0 = static_cast<OutputType>(0);

          int kernelIndex = clippedTransitionIndex0 - inputIndex;
          size_t signalIndex = boundary0 + kSizeOverTwo;
          while(kernelIndex < static_cast<int>(kernel.size())) {
            dotProduct0 += static_cast<OutputType>(
	      kernel[kernelIndex]
	      * static_cast<OutputType>(signal[signalIndex]));
            ++kernelIndex;
            ++signalIndex;
	  }
	  kernelIndex = clippedTransitionIndex0 - inputIndex - 1;
	  signalIndex = signal.size() - 1;
          while(kernelIndex >= 0) {
            dotProduct0 += static_cast<OutputType>(
	      kernel[kernelIndex]
	      * static_cast<OutputType>(signal[signalIndex]));
            --kernelIndex;
            --signalIndex;
	  }
          result[outputIndex] = dotProduct0;
	  ++inputIndex;
	  ++outputIndex;
	}


	inputIndex = clippedTransitionIndex1;
	outputIndex += clippedTransitionIndex1 - clippedTransitionIndex0;
	while(inputIndex < boundary1) {
          OutputType dotProduct0 = static_cast<OutputType>(0);

          int kernelStopIndex = boundary1 - inputIndex;
          int kernelIndex = 0;
          size_t signalIndex = inputIndex - kSizeOverTwo;
          while(kernelIndex < kernelStopIndex) {
            dotProduct0 += static_cast<OutputType>(
	      kernel[kernelIndex]
	      * static_cast<OutputType>(signal[signalIndex]));
            ++kernelIndex;
            ++signalIndex;
          }
	  signalIndex = 0;
          while(kernelIndex < static_cast<int>(kernel.size())) {
            dotProduct0 += static_cast<OutputType>(
	      kernel[kernelIndex]
	      * static_cast<OutputType>(signal[signalIndex]));
            ++kernelIndex;
            ++signalIndex;
          }
          result[outputIndex] = dotProduct0;
	  ++inputIndex;
	  ++outputIndex;
        }

        correlate1DCommon<OutputType, KernelType, SignalType>(
	  kernel,
	  signal.begin() + clippedTransitionIndex0 - kSizeOverTwo,
	  signal.begin() + clippedTransitionIndex1 - kSizeOverTwo,
	  result.begin() + resultTransitionIndex0);
	return result;
      }

    } // namespace privateCode
    /// @endcond


    template <class OutputType, class KernelType, class SignalType>
    inline Array1D<OutputType>
    convolve1D(const Array1D<KernelType>& kernel,
	       const Array1D<SignalType>& signal,
	       ConvolutionStrategy strategy,
	       ConvolutionROI roi)
    {
      Array1D<KernelType> reversedKernel(kernel.size());
      std::reverse_copy(kernel.begin(), kernel.end(), reversedKernel.begin());
      return correlate1D<OutputType, KernelType, SignalType>(
	reversedKernel, signal, strategy, roi);
    }
    

    template <class OutputType, class KernelType, class SignalType,
	      class FillType>
    inline Array1D<OutputType>
    convolve1D(const Array1D<KernelType>& kernel,
	       const Array1D<SignalType>& signal,
	       ConvolutionStrategy strategy,
	       ConvolutionROI roi,
	       const FillType& fillValue)
    {
      Array1D<KernelType> reversedKernel(kernel.size());
      std::reverse_copy(kernel.begin(), kernel.end(), reversedKernel.begin());
      return correlate1D<OutputType, KernelType, SignalType>(
	reversedKernel, signal, strategy, roi, fillValue);
    }
    

    template <class OutputType, class KernelType, class SignalType>
    inline Array1D<OutputType>
    convolve1D(const Array1D<KernelType>& kernel,
	       const Array1D<SignalType>& signal,
	       ConvolutionStrategy strategy,
	       int boundary0,
	       int boundary1)
    {
      Array1D<KernelType> reversedKernel(kernel.size());
      std::reverse_copy(kernel.begin(), kernel.end(), reversedKernel.begin());
      return correlate1D<OutputType, KernelType, SignalType>(
	reversedKernel, signal, strategy, boundary0, boundary1);
    }
    

    template <class OutputType, class KernelType, class SignalType,
	      class FillType>
    inline Array1D<OutputType>
    convolve1D(const Array1D<KernelType>& kernel,
	       const Array1D<SignalType>& signal,
	       ConvolutionStrategy strategy,
	       int boundary0,
	       int boundary1,
	       const FillType& fillValue)
    {
      Array1D<KernelType> reversedKernel(kernel.size());
      std::reverse_copy(kernel.begin(), kernel.end(), reversedKernel.begin());
      return correlate1D<OutputType, KernelType, SignalType>(
	reversedKernel, signal, strategy, boundary0, boundary1, fillValue);
    }
    

    template <class OutputType, class KernelType, class SignalType>
    Array1D<OutputType>
    correlate1D(const Array1D<KernelType>& kernel,
		const Array1D<SignalType>& signal,
		ConvolutionStrategy strategy,
		ConvolutionROI roi)
    {
      switch(roi) {
      case BRICK_CONVOLVE_ROI_SAME:
      {
	int boundary0 = 0;
	int boundary1 = static_cast<int>(signal.size());
	return correlate1D<OutputType, KernelType, SignalType>(
	  kernel, signal, strategy, boundary0, boundary1);
        break;
      }
      case BRICK_CONVOLVE_ROI_VALID:
      {
	int boundary0 = static_cast<int>(kernel.size()) / 2;
	int boundary1 = static_cast<int>(signal.size()) - boundary0;
	return correlate1D<OutputType, KernelType, SignalType>(
	  kernel, signal, strategy, boundary0, boundary1);
        break;
      }
      case BRICK_CONVOLVE_ROI_FULL:
      {
	int boundary0 = -static_cast<int>(kernel.size()) / 2;
	int boundary1 = static_cast<int>(signal.size()) - boundary0;
	return correlate1D<OutputType, KernelType, SignalType>(
	  kernel, signal, strategy, boundary0, boundary1);
        break;
      }
      default:
        BRICK_THROW(brick::common::LogicException, "correlate1D()",
                    "Illegal value for roi argument.");
        break;
      }
      return Array1D<OutputType>();
    }
    

    template <class OutputType, class KernelType, class SignalType,
	      class FillType>
    Array1D<OutputType>
    correlate1D(const Array1D<KernelType>& kernel,
		const Array1D<SignalType>& signal,
		ConvolutionStrategy strategy,
		ConvolutionROI roi,
		const FillType& fillValue)
    {
      switch(roi) {
      case BRICK_CONVOLVE_ROI_SAME:
      {
	int boundary0 = 0;
	int boundary1 = static_cast<int>(signal.size());
	return correlate1D<OutputType, KernelType, SignalType, FillType>(
	  kernel, signal, strategy, boundary0, boundary1, fillValue);
        break;
      }
      case BRICK_CONVOLVE_ROI_VALID:
      {
	int boundary0 = static_cast<int>(kernel.size()) / 2;
	int boundary1 = static_cast<int>(signal.size()) - boundary0;
	return correlate1D<OutputType, KernelType, SignalType, FillType>(
	  kernel, signal, strategy, boundary0, boundary1, fillValue);
        break;
      }
      case BRICK_CONVOLVE_ROI_FULL:
      {
	int boundary0 = -static_cast<int>(kernel.size()) / 2;
	int boundary1 = static_cast<int>(signal.size()) - boundary0;
	return correlate1D<OutputType, KernelType, SignalType, FillType>(
	  kernel, signal, strategy, boundary0, boundary1, fillValue);
        break;
      }
      default:
        BRICK_THROW(brick::common::LogicException, "correlate1D()",
                    "Illegal value for roi argument.");
        break;
      }
      return Array1D<OutputType>();
    }

    
    template <class OutputType, class KernelType, class SignalType>
    Array1D<OutputType>
    correlate1D(const Array1D<KernelType>& kernel,
		const Array1D<SignalType>& signal,
		ConvolutionStrategy strategy,
		int boundary0,
		int boundary1)
    {
      if(kernel.size() % 2 != 1) {
        BRICK_THROW(brick::common::ValueException, "correlate1D()",
                    "Argument kernel must have an odd number of elements.");
      }
      // Note(xxx): is the following check necessary?
      if(kernel.size() > signal.size()) {
        BRICK_THROW(brick::common::ValueException, "correlate1D()",
                    "Argument kernel must not have more elements than "
                    "argument signal.");
      }

      switch(strategy) {
      case BRICK_CONVOLVE_TRUNCATE_RESULT:
        return privateCode::correlate1DTruncateResult<
	  OutputType, KernelType, SignalType>(kernel, signal);
        break;
      case BRICK_CONVOLVE_PAD_RESULT:
      case BRICK_CONVOLVE_PAD_SIGNAL:
        BRICK_THROW(brick::common::ValueException, "correlate1D()",
                    "The specified convolution strategy requires that a "
                    "fill value be specified.");
      case BRICK_CONVOLVE_ZERO_PAD_SIGNAL:
        return privateCode::correlate1DZeroPadSignal<
	  OutputType, KernelType, SignalType>(
	    kernel, signal, boundary0, boundary1);
        break;
      case BRICK_CONVOLVE_REFLECT_SIGNAL:
        return privateCode::correlate1DReflectSignal<
	  OutputType, KernelType, SignalType>(
	    kernel, signal, boundary0, boundary1);
        break;
      case BRICK_CONVOLVE_WRAP_SIGNAL:
        return privateCode::correlate1DWrapSignal<
	  OutputType, KernelType, SignalType>(
	    kernel, signal, boundary0, boundary1);
        break;
      default:
        BRICK_THROW(brick::common::LogicException, "correlate1D()",
                    "Illegal value for strategy argument.");
        break;
      }
      return Array1D<OutputType>();
    }
    

    template <class OutputType, class KernelType, class SignalType,
	      class FillType>
    Array1D<OutputType>
    correlate1D(const Array1D<KernelType>& kernel,
		const Array1D<SignalType>& signal,
		ConvolutionStrategy strategy,
		int boundary0,
		int boundary1,
		const FillType& fillValue)
    {
      if(kernel.size() % 2 != 1) {
        BRICK_THROW(brick::common::ValueException, "correlate1D()",
                    "Argument kernel must have an odd number of elements.");
      }
      // Note(xxx): is the following check necessary?
      if(kernel.size() > signal.size()) {
        BRICK_THROW(brick::common::ValueException, "correlate1D()",
                    "Argument kernel must not have more elements than "
                    "argument signal.");
      }

      switch(strategy) {
      case BRICK_CONVOLVE_TRUNCATE_RESULT:
        return privateCode::correlate1DTruncateResult<
	  OutputType, KernelType, SignalType>(kernel, signal);
        break;
      case BRICK_CONVOLVE_PAD_RESULT:
        return privateCode::correlate1DPadResult<
	  OutputType, KernelType, SignalType>(
	    kernel, signal, boundary0, boundary1,
	    static_cast<OutputType>(fillValue));
        break;
      case BRICK_CONVOLVE_PAD_SIGNAL:
        return privateCode::correlate1DPadSignal<
	  OutputType, KernelType, SignalType>(
	    kernel, signal, boundary0, boundary1,
	    static_cast<SignalType>(fillValue));
        break;
      case BRICK_CONVOLVE_ZERO_PAD_SIGNAL:
        return privateCode::correlate1DZeroPadSignal<
	  OutputType, KernelType, SignalType>(
	    kernel, signal, boundary0, boundary1);
        break;
      case BRICK_CONVOLVE_REFLECT_SIGNAL:
        return privateCode::correlate1DReflectSignal<
	  OutputType, KernelType, SignalType>(
	    kernel, signal, boundary0, boundary1);
        break;
      case BRICK_CONVOLVE_WRAP_SIGNAL:
        return privateCode::correlate1DWrapSignal<
	  OutputType, KernelType, SignalType>(
	    kernel, signal, boundary0, boundary1);
        break;
      default:
        BRICK_THROW(brick::common::LogicException, "correlate1D()",
                    "Illegal value for strategy argument.");
        break;
      }
      return Array1D<OutputType>();
    }

  } // namespace numeric

} // namespace brick

#endif /* #ifndef BRICK_NUMERIC_CONVOLVE1D_IMPL_HH */
