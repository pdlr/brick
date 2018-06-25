/**
***************************************************************************
* @file brick/numeric/convolve2D.hh
*
* Header file declaring 1D correlation and convolution functions.
*
* Copyright (C) 2006-2007, 2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_CONVOLVE2D_HH
#define BRICK_NUMERIC_CONVOLVE2D_HH

#include <brick/numeric/array2D.hh>
#include <brick/numeric/convolutionStrategy.hh>
#include <brick/numeric/index2D.hh>

namespace brick {

  namespace numeric {

    /** Unstable: interface subject to change. **/
    template <class OutputType, class AccumulatorType,
	      class KernelType, class SignalType>
    inline Array2D<OutputType>
    convolve2D(const Array2D<KernelType>& kernel,
	       const Array2D<SignalType>& signal,
	       ConvolutionStrategy strategy = BRICK_CONVOLVE_PAD_RESULT,
	       ConvolutionROI roi = BRICK_CONVOLVE_ROI_SAME);


    /** Unstable: interface subject to change. **/
    template <class OutputType, class AccumulatorType,
	      class KernelType, class SignalType, class FillType>
    inline Array2D<OutputType>
    convolve2D(const Array2D<KernelType>& kernel,
	       const Array2D<SignalType>& signal,
	       ConvolutionStrategy strategy,
	       ConvolutionROI roi,
	       const FillType& fillValue);


    /** Unstable: interface subject to change. **/
    template <class OutputType, class AccumulatorType,
	      class KernelType, class SignalType>
    Array2D<OutputType>
    convolve2D(const Array2D<KernelType>& kernel,
	       const Array2D<SignalType>& signal,
	       ConvolutionStrategy strategy,
	       const Index2D& corner0,
	       const Index2D& corner1);


    /** Unstable: interface subject to change. **/
    template <class OutputType, class AccumulatorType,
	      class KernelType, class SignalType, class FillType>
    Array2D<OutputType>
    convolve2D(const Array2D<KernelType>& kernel,
	       const Array2D<SignalType>& signal,
	       ConvolutionStrategy strategy,
	       const Index2D& corner0,
	       const Index2D& corner1,
	       const FillType& fillValue);


    /** Unstable: interface subject to change. **/
    template <class OutputType, class AccumulatorType,
	      class KernelType, class SignalType>
    inline Array2D<OutputType>
    correlate2D(const Array2D<KernelType>& kernel,
		const Array2D<SignalType>& signal,
		ConvolutionStrategy strategy = BRICK_CONVOLVE_PAD_RESULT,
		ConvolutionROI roi = BRICK_CONVOLVE_ROI_SAME);


    /** Unstable: interface subject to change. **/
    template <class OutputType, class AccumulatorType,
	      class KernelType, class SignalType, class FillType>
    Array2D<OutputType>
    correlate2D(const Array2D<KernelType>& kernel,
		const Array2D<SignalType>& signal,
		ConvolutionStrategy strategy,
		ConvolutionROI roi,
		const FillType& fillValue);


    /** Unstable: interface subject to change. **/
    template <class OutputType, class AccumulatorType,
	      class KernelType, class SignalType>
    Array2D<OutputType>
    correlate2D(const Array2D<KernelType>& kernel,
		const Array2D<SignalType>& signal,
		ConvolutionStrategy strategy,
		const Index2D& corner0,
		const Index2D& corner1);


    /** Unstable: interface subject to change. **/
    template <class OutputType, class AccumulatorType,
	      class KernelType, class SignalType, class FillType>
    Array2D<OutputType>
    correlate2D(const Array2D<KernelType>& kernel,
		const Array2D<SignalType>& signal,
		ConvolutionStrategy strategy,
		const Index2D& corner0,
		const Index2D& corner1,
		const FillType& fillValue);


  } // namespace numeric

} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/numeric/convolve2D_impl.hh>

#endif /* #ifndef BRICK_NUMERIC_CONVOLVE2D_HH */
