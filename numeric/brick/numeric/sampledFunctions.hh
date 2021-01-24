/**
***************************************************************************
* @file brick/numeric/sampledFunctions.hh
*
* Header file declaring functions that return sampled versions of
* common functions, such as gaussians, etc.
*
* Copyright (C) 2006-2012 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_NUMERIC_SAMPLEDFUNCTIONS_HH
#define BRICK_NUMERIC_SAMPLEDFUNCTIONS_HH

#include <brick/numeric/array1D.hh>

namespace brick {

  namespace numeric {

    /**
     * This function returns a Blackman-Harris window of the specified
     * size.  The Blackman-Harris window is a low-resolution
     * (high-dynamic-range) window.  Use this if you want low
     * side-lobes.
     *
     * @param windowSize This argument specifies the size of the
     * required window in samples.
     *
     * @return The return value is an array containing the sampled
     * window function.
     */
    template <class Type>
    Array1D<Type>
    getBlackmanHarrisWindow1D(size_t windowSize);


    /**
     * This function generates a Blackman-Harris window of the specified
     * size.  The Blackman-Harris window is a low-resolution
     * (high-dynamic-range) window.  Use this if you want low
     * side-lobes.
     *
     * @param beginIter This argument and the next specify a sequence
     * into which the sampled window values will be copied.
     *
     * @param endIter This argument and the previous specify a
     * sequence into which the sampled window values will be copied.
     */
    template <class Iter, class Type>
    void
    getBlackmanHarrisWindow1D(Iter beginIter, Iter endIter);


    /**
     * This function returns an array in which the elements are
     * sampled from a 1D Gaussian.  The maximum of the Gaussian is at
     * the center of the returned array, and the distance between
     * adjacent samples (elements of the array) is assumed to be
     * exactly 1.0.
     *
     * @param sigma This argument specifies the desired standard
     * deviation of the Gaussian.
     *
     * @param size This argument specifies how many elements should be
     * in the returned array.  If size is set to 0, the returned value
     * will have the smallest odd number of elements that is greater
     * than or equal to 6.0 * sigma.
     *
     * @param normalize This argument specifies whether, after the
     * elements of the array have been computed, the array should be
     * rescaled so that its elements sum to 1.0.
     *
     * @return The return value is an array containing the sampled
     * values.
     */
    template <class Type>
    Array1D<Type>
    getGaussian1D(double sigma,
                  size_t size = 0,
                  bool normalize = false);


    /**
     * This function generates a sampled a 1D Gaussian.  The maximum
     * of the Gaussian is at the center of the returned array, and the
     * distance between adjacent samples (elements of the array) is
     * assumed to be exactly 1.0.
     *
     * @param beginIter This argument and the next specify a sequence
     * into which the sampled window values will be copied.
     *
     * @param endIter This argument and the previous specify a
     * sequence into which the sampled window values will be copied.
     *
     * @param sigma This argument specifies the desired standard
     * deviation of the Gaussian.  If this argument is less than or
     * equal to 0.0, then it will be automatically set so that the the
     * entire sequence spans approximately 6*sigma.
     *
     * @param normalize This argument specifies whether, after the
     * Gaussian has been computed, its elements should be rescaled so
     * that they sum to 1.0.
     */
    template <class Iter, class Type>
    void
    getGaussian1D(Iter beginIter, Iter endIter,
                  Type sigma = -1.0, bool normalize = false);


    /**
     * This function returns a Hamming window of the specified size.
     * The Hamming window is a moderate-resolution window.  Use this
     * if you can tolerate not-too-high side-lobes and need better
     * resolution than you get with a Blackman-Harris window.
     *
     * @param windowSize This argument specifies the size of the
     * required window in samples.
     *
     * @return The return value is an array containing the sampled
     * window function.
     */
    template <class Type>
    Array1D<Type>
    getHammingWindow1D(size_t windowSize);


    /**
     * This function generates a Hamming window of the specified size.
     * The Hamming window is a moderate-resolution window.  Use this
     * if you can tolerate not-too-high side-lobes and need better
     * resolution than you get with a Blackman-Harris window.
     *
     * @param beginIter This argument and the next specify a sequence
     * into which the sampled window values will be copied.
     *
     * @param endIter This argument and the previous specify a
     * sequence into which the sampled window values will be copied.
     */
    template <class Type, class Iter>
    void
    getHammingWindow1D(Iter beginIter, Iter endIter);


    /**
     * This function returns a Hann window of the specified size.
     * The Hann window is a moderate-resolution window with higher
     * dynamic range (and lower resolution) than the Hamming window.
     *
     * @param windowSize This argument specifies the size of the
     * required window in samples.
     *
     * @return The return value is an array containing the sampled
     * window function.
     */
    template <class Type>
    Array1D<Type>
    getHannWindow1D(size_t windowSize);


    /**
     * This function generates a Hann window of the specified size.
     * The Hann window is a moderate-resolution window with higher
     * dynamic range (and lower resolution) than the Hamming window.
     *
     * @param beginIter This argument and the next specify a sequence
     * into which the sampled window values will be copied.
     *
     * @param endIter This argument and the previous specify a
     * sequence into which the sampled window values will be copied.
     */
    template <class Iter, class Type>
    void
    getHannWindow1D(Iter beginIter, Iter endIter);

  } // namespace numeric

} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/numeric/sampledFunctions_impl.hh>

#endif /* #ifndef BRICK_NUMERIC_SAMPLEDFUNCTIONS_HH */
