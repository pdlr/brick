/**
***************************************************************************
* @file brick/numeric/sampledFunctions_impl.hh
*
* Header file defining the inline and template functions declared in
* sampledFunctions.hh.
*
* Copyright (C) 2006-2012 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_NUMERIC_SAMPLEDFUNCTIONS_IMPL_HH
#define BRICK_NUMERIC_SAMPLEDFUNCTIONS_IMPL_HH

// This file is included by sampledFunctions.hh, and should not be
// directly included by user code, so no need to include
// sampledFunctions.hh here.
//
// #include <brick/numeric/sampledFunctions.hh>

#include <cmath>
#include <brick/common/constants.hh>
#include <brick/common/mathFunctions.hh>
#include <brick/numeric/functional.hh>
#include <brick/numeric/utilities.hh>

namespace brick {

  namespace numeric {

    // This function returns a Blackman-Harris window of the specified
    // size.
    template <class Type>
    Array1D<Type>
    getBlackmanHarrisWindow1D(size_t windowSize)
    {
      Array1D<Type> result(windowSize);
      getBlackmanHarrisWindow1D<Type>(result.begin(), result.end());
      return result;
    }


    // This function generates a Blackman-Harris window of the specified
    // size.
    template <class Iter, class Type>
    void
    getBlackmanHarrisWindow1D(Iter beginIter, Iter endIter)
    {
      const unsigned int windowSize = static_cast<unsigned int>(
        endIter - beginIter);
      const Type piOverNMinus1 = static_cast<Type>(
        brick::common::constants::pi / (windowSize - 1));
      const Type k0 = 2 * piOverNMinus1;
      const Type k1 = 4 * piOverNMinus1;
      const Type k2 = 6 * piOverNMinus1;
      const Type c0 = static_cast<Type>(-0.48829);
      const Type c1 = static_cast<Type>(0.14128);
      const Type c2 = static_cast<Type>(-0.01168);
      const Type c3 = static_cast<Type>(0.35875);
      unsigned int ii = 0;
      while(beginIter != endIter) {
        *beginIter = (c0 * brick::common::cosine(k0 * ii)
                      + c1 * brick::common::cosine(k1 * ii)
                      + c2 * brick::common::cosine(k2 * ii)
                      + c3);
        ++ii;
        ++beginIter;
      }
    }


    // This function returns an array in which the elements are
    // sampled from a 1D Gaussian.
    template <class Type>
    Array1D<Type>
    getGaussian1D(double sigma, size_t size, bool normalize)
    {
      if(size == 0) {
        size = static_cast<size_t>(6.0 * sigma + 0.5);
        if(size % 2 == 0) {
          ++size;
        }
      }
      Array1D<Type> result(size);
      getGaussian1D(result.begin(), result.end(), sigma, normalize);
      return result;
    }


    // This function generates a sampled a 1D Gaussian.
    template <class Iter, class Type>
    void
    getGaussian1D(Iter beginIter, Iter endIter, Type sigma, bool normalize)
    {
      const unsigned int windowSize = static_cast<unsigned int>(
        endIter - beginIter);
      const Type increment = static_cast<Type>(1.0);
      Type xx = static_cast<Type>((1.0 - windowSize)/2.0);

      if(sigma <= 0.0) {
        sigma = windowSize / static_cast<Type>(6.0);
      }
      Gaussian1DFunctor<Type> functor(sigma);

      if(normalize) {
        Type runningTotal = static_cast<Type>(0.0);
        Iter myIter = beginIter;
        while(myIter != endIter) {
          *myIter = functor(xx);
          runningTotal += *myIter;
          xx += increment;
          ++myIter;
        }
        while(beginIter != endIter) {
          *beginIter /= runningTotal;
          ++beginIter;
        }
      } else {
        while(beginIter != endIter) {
          *beginIter = functor(xx);
          xx += increment;
          ++beginIter;
        }
      }
    }


    // This function returns a Hamming window of the specified size.
    template <class Type>
    Array1D<Type>
    getHammingWindow1D(size_t windowSize)
    {
      Array1D<Type> result(windowSize);
      getHammingWindow1D<Type>(result.begin(), result.end());
      return result;
    }


    // This function generates a Hamming window of the specified size.
    template <class Type, class Iter>
    void
    getHammingWindow1D(Iter beginIter, Iter endIter)
    {
      const unsigned int windowSize = static_cast<unsigned int>(
        endIter - beginIter);
      const Type piOverNMinus1 = static_cast<Type>(
        brick::common::constants::pi / (windowSize - 1));
      const Type k0 = 2 * piOverNMinus1;
      const Type c0 = static_cast<Type>(-0.46);
      const Type c1 = static_cast<Type>(0.54);

      unsigned int ii = 0;
      while(beginIter != endIter) {
        *beginIter = (c0 * brick::common::cosine(k0 * ii) + c1);
        ++ii;
        ++beginIter;
      }
    }


    // This function returns a Hann window of the specified size.
    template <class Type>
    Array1D<Type>
    getHannWindow1D(size_t windowSize)
    {
      Array1D<Type> result(windowSize);
      getHannWindow1D<Type>(result.begin(),result.end());
      return result;
    }


    // This function generates a Hann window of the specified size.
    template <class Iter, class Type>
    void
    getHannWindow1D(Iter beginIter, Iter endIter)
    {
      const unsigned int windowSize = static_cast<unsigned int>(
        endIter - beginIter);
      const Type piOverNMinus1 = static_cast<Type>(
        brick::common::constants::pi / (windowSize - 1));
      const Type k0 = 2 * piOverNMinus1;
      const Type c0 = static_cast<Type>(-0.5);
      const Type c1 = static_cast<Type>(0.5);

      unsigned int ii = 0;
      while(beginIter != endIter) {
        *beginIter = (c0 * brick::common::cosine(k0 * ii) + c1);
        ++ii;
        ++beginIter;
      }
    }

  } // namespace numeric

} // namespace brick

#endif /* #ifndef BRICK_NUMERIC_SAMPLEDFUNCTIONS_IMPL_HH */
