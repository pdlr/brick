/**
***************************************************************************
* @file brick/numeric/fft.hh
*
* Header file declaring functions that implement fast fourier
* transforms.  These functions are not intended to replace the
* excellent and highly-optimized alternative implementations available
* elsewhere.  They are only to provide a convenient substitute when
* other libraries are not available, and fast execution is not
* important.
*
* Copyright (C) 2017 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_NUMERIC_FFT_HH
#define BRICK_NUMERIC_FFT_HH

#include <brick/numeric/array1D.hh>

namespace brick {

  namespace numeric {

    /// Warning: this interface is not yet stable.
    template <class ComplexType>
    Array1D<ComplexType>
    computeFFT(Array1D<ComplexType> const& inputSignal);
    
  } // namespace numeric

} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/numeric/fft_impl.hh>

#endif /* #ifndef BRICK_NUMERIC_FFT_HH */
