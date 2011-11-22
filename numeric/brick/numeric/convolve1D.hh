/**
***************************************************************************
* @file brick/numeric/convolve1D.hh
*
* Header file declaring 1D correlation and convolution functions.
*
* Copyright (C) 2006-2007, 2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_CONVOLVE1D_HH
#define BRICK_NUMERIC_CONVOLVE1D_HH

#include <brick/numeric/array1D.hh>
#include <brick/numeric/convolutionStrategy.hh>

namespace brick {

  namespace numeric {
    
    /** Unstable: interface subject to change. **/
    template <class OutputType, class KernelType, class SignalType>
    inline Array1D<OutputType>
    convolve1D(const Array1D<KernelType>& kernel,
	       const Array1D<SignalType>& signal,
	       ConvolutionStrategy strategy = BRICK_CONVOLVE_PAD_RESULT,
	       ConvolutionROI roi = BRICK_CONVOLVE_ROI_SAME);


    /** Unstable: interface subject to change. **/
    template <class OutputType, class KernelType, class SignalType,
	      class FillType>
    inline Array1D<OutputType>
    convolve1D(const Array1D<KernelType>& kernel,
	       const Array1D<SignalType>& signal,
	       ConvolutionStrategy strategy,
	       ConvolutionROI roi,
	       const FillType& fillValue);
  

    /** Unstable: interface subject to change. **/
    template <class OutputType, class KernelType, class SignalType>
    Array1D<OutputType>
    convolve1D(const Array1D<KernelType>& kernel,
	       const Array1D<SignalType>& signal,
	       ConvolutionStrategy strategy,
	       int boundary0,
	       int boundary1);

    
    /** Unstable: interface subject to change. **/
    template <class OutputType, class KernelType, class SignalType,
	      class FillType>
    Array1D<OutputType>
    convolve1D(const Array1D<KernelType>& kernel,
	       const Array1D<SignalType>& signal,
	       ConvolutionStrategy strategy,
	       int boundary0,
	       int boundary1,
	       const FillType& fillValue);

    
    /** Unstable: interface subject to change. **/
    template <class OutputType, class KernelType, class SignalType>
    inline Array1D<OutputType>
    correlate1D(const Array1D<KernelType>& kernel,
		const Array1D<SignalType>& signal,
		ConvolutionStrategy strategy = BRICK_CONVOLVE_PAD_RESULT,
		ConvolutionROI roi = BRICK_CONVOLVE_ROI_SAME);


    /** Unstable: interface subject to change. **/
    template <class OutputType, class KernelType, class SignalType,
	      class FillType>
    Array1D<OutputType>
    correlate1D(const Array1D<KernelType>& kernel,
		const Array1D<SignalType>& signal,
		ConvolutionStrategy strategy,
		ConvolutionROI roi,
		const FillType& fillValue);

    
    /** Unstable: interface subject to change. **/
    template <class OutputType, class KernelType, class SignalType>
    Array1D<OutputType>
    correlate1D(const Array1D<KernelType>& kernel,
		const Array1D<SignalType>& signal,
		ConvolutionStrategy strategy,
		int boundary0,
		int boundary1);

    
    /** Unstable: interface subject to change. **/
    template <class OutputType, class KernelType, class SignalType,
	      class FillType>
    Array1D<OutputType>
    correlate1D(const Array1D<KernelType>& kernel,
		const Array1D<SignalType>& signal,
		ConvolutionStrategy strategy,
		int boundary0,
		int boundary1,
		const FillType& fillValue);

  } // namespace numeric

} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/numeric/convolve1D_impl.hh>

#endif /* #ifndef BRICK_NUMERIC_CONVOLVE1D_HH */
