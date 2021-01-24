/**
***************************************************************************
* @file brick/numeric/convolveND_impl.hh
*
* Header file defining functions for doing N-dimensional convolution
* and correlation.
*
* Copyright (C) 2008, 2011 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
***************************************************************************
**/

#ifndef BRICK_NUMERIC_CONVOLVEND_IMPL_HH
#define BRICK_NUMERIC_CONVOLVEND_IMPL_HH

// This file is included by convolveND.hh, and should not be directly included
// by user code, so no need to include convolveND.hh here.
//
// #include <brick/numeric/convolveND.hh>

#include <iostream>
#include <string>
#include <brick/common/exception.hh>
#include <brick/numeric/arrayND.hh>
#include <brick/numeric/convolutionStrategy.hh>


namespace brick {

  namespace numeric {

    namespace privateCode {

      /// @cond privateCode
      inline bool
      convolveNDUpdatePosition(Array1D<size_t>& position,
                               Array1D<size_t> const& lowerBound,
                               Array1D<size_t> const& upperBound)
      {
        for(size_t ii = position.size() - 1; ii < position.size(); --ii) {
          if(++(position[ii]) < upperBound[ii]) {
            return true;
          }
          position[ii] = lowerBound[ii];
        }
        return false;
      }

    } // namespace privateCode
    /// @endcond


    // Unstable: interface subject to change.
    template <class OutputType, class AccumulatorType,
	      class KernelType, class SignalType, size_t Dimension>
    ArrayND<Dimension, OutputType>
    convolve(const Array1D<KernelType>& kernel,
             const ArrayND<Dimension, SignalType>& signal,
             size_t axis,
             ConvolutionStrategy strategy,
             ConvolutionROI roi)
    {
      Array1D<AccumulatorType> reversedKernel(kernel.size());
      for(size_t ii = 0; ii < kernel.size(); ++ii) {
        reversedKernel[kernel.size() - ii - 1] =
          static_cast<AccumulatorType>(kernel[ii]);
      }

      ArrayND<Dimension, OutputType> result(signal.getShape());
      result = OutputType(0);

      if(strategy != BRICK_CONVOLVE_PAD_RESULT) {
        BRICK_THROW(brick::common::NotImplementedException, "convolve()",
                  "Only ConvolutionStrategy BRICK_CONVOLVE_PAD_RESULT is "
                  "currently implemented.");
      }
      if(roi != BRICK_CONVOLVE_ROI_SAME) {
        BRICK_THROW(brick::common::NotImplementedException, "convolve()",
                  "Only ConvolutionROI BRICK_CONVOLVE_ROI_SAME is "
                  "currently implemented.");
      }

      // We'll be using 1D indexing into the ArrayND instance.  This
      // stride tells us how much we need to increment the index in
      // order to move one element along the requested axis, and
      // resultOffset tells us the offset between the first element of
      // the convolution kernel and the center of the convolution
      // kernel.
      size_t stride = signal.getStride(axis);
      size_t resultOffset = (reversedKernel.size() / 2) * stride;

      // This array keeps track of where we are in the ND array.  It
      // will be stepped in raster order through each position at
      // which the convolution kernel is applied to signal.
      Array1D<size_t> elementPosition(Dimension);
      elementPosition = 0;

      // Establish upper and lower bounds on elementPosition for valid
      // indexing.
      Array1D<size_t> positionLowerBound = elementPosition.copy();
      Array1D<size_t> positionUpperBound = signal.getShape().copy();
      positionUpperBound[axis] -= (reversedKernel.size() - 1);

      do {
        size_t signalIndex = signal.flattenIndex(elementPosition);
        size_t resultIndex = signalIndex + resultOffset;

        AccumulatorType accumulator = AccumulatorType(0);
        for(size_t ii = 0; ii < reversedKernel.size(); ++ii) {
          accumulator += (reversedKernel[ii]
                          * static_cast<AccumulatorType>(signal(signalIndex)));
          signalIndex += stride;
        }
        result(resultIndex) = static_cast<OutputType>(accumulator);
      } while(privateCode::convolveNDUpdatePosition(
                elementPosition, positionLowerBound, positionUpperBound));
      return result;
    }


  } // namespace numeric

} // namespace brick

#endif /* #ifdef BRICK_NUMERIC_CONVOLVEND_IMPL_HH */
