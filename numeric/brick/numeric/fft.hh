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

    /**
     * Warning: this interface is not yet stable.
     *
     * This function currently implements an un-optimized version of
     * the Cooley-Tukey radix-2 recursive decimation-in-time FFT
     * algorithm.  Over time it may grow to dispatch to other FFT
     * algorithms as well, based on the input signal or a
     * configuration argument.
     *
     * The goal here isn't to make an FFT implementation that competes
     * with with the more optimized versions available, just to have a
     * quick and easy FFT for use when other libraries aren't handy.
     *
     * @param inputSignal This argument is the complex-valued signal
     * from which to compute the Fourier transform.  For now, the
     * number of elements in this sequence must be a power of two.
     *
     * @return The return value is the Discrete Fourier transform of
     * argument inputSignal.  Assuming there are N elements in
     * inputSignal, the return value will also have N elements.  The
     * nth element of the return value is the Fourier coefficient for
     * frequency 2*pi*n/N radians per sample, where n starts at zero
     * and goes to (N - 1).  Note that elements with n > (N/2) are
     * above Nyquist, and in real signals will alias back to
     * frequencies symmetrically distributed around n = N/2.
     */
    template <class ComplexType>
    Array1D<ComplexType>
    computeFFT(Array1D<ComplexType> const& inputSignal);

  } // namespace numeric

} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/numeric/fft_impl.hh>

#endif /* #ifndef BRICK_NUMERIC_FFT_HH */
