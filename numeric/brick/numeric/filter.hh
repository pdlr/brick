/**
***************************************************************************
* @file brick/numeric/filter.hh
*
* Header file declaring functions for performing filtering of arrays.
*
* Copyright (C) 2012 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_NUMERIC_FILTER_HH
#define BRICK_NUMERIC_FILTER_HH

#include <brick/numeric/array1D.hh>
#include <brick/numeric/array2D.hh>

namespace brick {

  namespace numeric {

    /**
     * This function filters each column of an array by convolving
     * (correlating, actually) with the specified filter kernel.  The
     * result is returned through the pre-allocated outputArray
     * argument.  It is an error if the size of outputArray differs
     * from inputArray.  The returned array is the same size as the
     * input array, and regions of the returned array for which the
     * convolution result is not valid (because the convolution kernel
     * extends off the edge of the array) will be left untouched.  The
     * convolution will be performed using the math with precision
     * specified by KernelType, and the result will be cast to
     * OutputType.
     *
     * @param outputArray The result of the filtering will be returned
     * through this argument, which must have the same number of rows
     * and columns as inputArray.
     *
     * @param inputArray This argument is the Array to be filtered.
     *
     * @param kernel This argument specifies the correlation kernel to
     * be used for filtering.
     */
    template<class InputType, class OutputType, class KernelType>
    void
    filterColumns(
      Array2D<OutputType>& outputArray,
      Array2D<InputType> const& inputArray,
      Array1D<KernelType> const& kernel);


    /**
     * This function filters each column of an array by convolving
     * (correlating, actually) with the specified filter kernel.  The
     * result is returned through the pre-allocated outputArray
     * argument.  It is an error if the size of outputArray differs
     * from inputArray.  The returned array is the same size as the
     * input array, and regions of the returned array for which the
     * convolution result is not valid (because the convolution kernel
     * extends off the edge of the array) will be left untouched.  The
     * convolution will be performed using the math with precision
     * specified by KernelType, and the result will be cast to
     * OutputType.
     *
     * @param outputArray The result of the filtering will be returned
     * through this argument, which must have the same number of rows
     * and columns as inputArray.
     *
     * @param inputArray This argument is the Array to be filtered.
     *
     * @param kernel This argument specifies the correlation kernel to
     * be used for filtering.
     */
    template<class InputType, class OutputType, class KernelType>
    void
    filterRows(
      Array2D<OutputType>& outputArray,
      Array2D<InputType> const& inputArray,
      Array1D<KernelType> const& kernel);

  } // namespace numeric

} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/numeric/filter_impl.hh>

#endif /* #ifndef BRICK_NUMERIC_FILTER_HH */
