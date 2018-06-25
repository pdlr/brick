/**
***************************************************************************
* @file brick/numeric/fft_impl.hh
*
* Header file defining inline functions and function templates that
* are declared in brick/numeric/fft.hh.
*
* Copyright (C) 2017 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_NUMERIC_FFT_IMPL_HH
#define BRICK_NUMERIC_FFT_IMPL_HH

// This file is included by fft.hh, and should not be directly
// included by user code, so no need to include fft.hh here.
//
// #include <brick/numeric/fft.hh>

#include <brick/common/constants.hh>
#include <brick/common/exception.hh>
#include <brick/common/mathFunctions.hh>

namespace brick {

  namespace numeric {

    namespace privateCode {

      template <class ComplexType>
      Array1D<ComplexType>
      computeRadix2Twiddles(std::size_t const count)
      {
        typedef typename ComplexType::value_type FloatType;
        Array1D<ComplexType> twiddleFactors(count);
        for(std::size_t ii = 0; ii < count; ++ii) {
          FloatType exponent(ii * (-brick::common::constants::twoPi / count));
          twiddleFactors[ii] = ComplexType(
            brick::common::cosine(exponent),
            brick::common::sine(exponent));
        }
        return twiddleFactors;
      }


      bool
      isPowerOfTwo(std::size_t signalLength)
      {
        double exponent = std::log(double(signalLength)) / std::log(2.0);
        int integerExponent = static_cast<int>(exponent + 0.5);
        std::size_t checkValue = std::size_t(1) << integerExponent;
        return checkValue == signalLength;
      }


      template <class ComplexType>
      void
      recursiveRadix2FFT(Array1D<ComplexType>& result,
                         std::size_t const outputIndex,
                         Array1D<ComplexType> const& inputSignal,
                         std::size_t const inputIndex,
                         std::size_t const count,
                         std::size_t const stride,
                         Array1D<ComplexType> const& twiddleFactors)
      {
        // Handle the trivial case.
        if(count == 1) {
          result[outputIndex] = inputSignal[inputIndex];
          return;
        }

        // Divide and conquer.
        std::size_t countOver2 = count / 2;
        std::size_t strideTimes2 = stride * 2;

        // Compute the FFT of the even elements.
        recursiveRadix2FFT(result, outputIndex,
                           inputSignal, inputIndex,
                           countOver2, strideTimes2, twiddleFactors);

        // Compute the FFT of the odd elements.
        recursiveRadix2FFT(result, outputIndex + countOver2,
                           inputSignal, inputIndex + stride,
                           countOver2, strideTimes2, twiddleFactors);

        // Apply a frequency-dependent phase shift ("twiddle factor")
        // to the odd FFT, and combine with the Even FFT to get the
        // whole.
        for(std::size_t ii = 0; ii < countOver2; ++ii) {
          std::size_t evenIndex = outputIndex + ii;
          std::size_t oddIndex = evenIndex + countOver2;

          ComplexType& evenPart = result[evenIndex];
          ComplexType& oddPart = result[oddIndex];

          ComplexType temp0 = evenPart;
          ComplexType temp1 = twiddleFactors[ii * stride] * oddPart;
          evenPart = temp0 + temp1;
          oddPart = temp0 - temp1;
        }
      }

    } // namespace privateCode


    template <class ComplexType>
    Array1D<ComplexType>
    computeFFT(Array1D<ComplexType> const& inputSignal)
    {
      if(not privateCode::isPowerOfTwo(inputSignal.size())) {
        BRICK_THROW(brick::common::NotImplementedException,
                    "computeFFT()",
                    "Length of index signal must currently be a power of two.");
      }

      // This code implements the recursive radix-2 decimation-in-time
      // algorithm of Cooley and Tukey.  This algorithm recognizes
      // that the naive FFT calculation (which is O(N^2), where N
      // is the number of elements of the input signal), can be
      // written as the sum of two sub-FFTs.  The first of these is
      // the FFT of the even elements of the input sequence, and the
      // second is the FFT of the odd elements of the input sequence,
      // with a frequency-dependent phase offset.  These phase offsets
      // are traditionally called "twiddle factors," and are computed
      // up-front.
      //
      // These two sub-FFTs have half as many elements (each) as the
      // whole FFT.  FFTs are periodic, though, so only N/2 FFT values
      // need be computed for each sub-FFT.  Repeating this trick for
      // each sub-FFT gives a recursive algorithm with O(N*log(N))
      // complexity.
      Array1D<ComplexType> result(inputSignal.size());
      std::size_t inputIndex = 0;
      std::size_t outputIndex = 0;
      std::size_t inputStride = 1;
      std::size_t count = inputSignal.size();
      Array1D<ComplexType> twiddleFactors =
        privateCode::computeRadix2Twiddles<ComplexType>(count);
      privateCode::recursiveRadix2FFT(
        result, outputIndex, inputSignal, inputIndex, count, inputStride,
        twiddleFactors);

      return result;
    }

  } // namespace numeric

} // namespace brick


#endif /* #ifndef BRICK_NUMERIC_FFT_IMPL_HH */
