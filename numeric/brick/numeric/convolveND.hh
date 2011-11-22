/**
***************************************************************************
* @file brick/numeric/convolveND.hh
*
* Header file declaring founctions for doing N-dimensional convolution
* and correlation.
*
* Copyright (C) 2008, 2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
***************************************************************************
**/

#ifndef BRICK_NUMERIC_CONVOLVEND_HH
#define BRICK_NUMERIC_CONVOLVEND_HH

#include <iostream>
#include <string>
#include <brick/common/exception.hh>
#include <brick/numeric/arrayND.hh>
#include <brick/numeric/convolutionStrategy.hh>


namespace brick {

  namespace numeric {
    

    /** Unstable: interface subject to change. **/
    template <class OutputType, class AccumulatorType,
	      class KernelType, class SignalType, size_t Dimension>
    ArrayND<Dimension, OutputType>
    convolve(const Array1D<KernelType>& kernel,
             const ArrayND<Dimension, SignalType>& signal,
             size_t axis,
             ConvolutionStrategy strategy = BRICK_CONVOLVE_PAD_RESULT,
             ConvolutionROI roi = BRICK_CONVOLVE_ROI_SAME);

    
  } // namespace numeric

} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/numeric/convolveND_impl.hh>

#endif /* #ifdef BRICK_NUMERIC_CONVOLVEND_HH */
