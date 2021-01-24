/**
***************************************************************************
* @file brick/numeric/boxIntegrator2D_impl.hh
*
* Header file defining the inline and template functions declared in
* boxIntegrator2D.hh.
*
* Copyright (C) 2006,2012 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_NUMERIC_BOXINTEGRATOR2D_IMPL_HH
#define BRICK_NUMERIC_BOXINTEGRATOR2D_IMPL_HH

// This file is included by boxIntegrator2D.hh, and should not be
// directly included by user code, so no need to include
// boxIntegrator2D.hh here.
//
// #include <brick/numeric/boxIntegrator2D.hh>

#include <brick/common/functional.hh>

namespace brick {

  namespace numeric {


    // This constructor performs almost no work, and simply
    // initializes the class instance to a "zero" state.
    template <class Type0, class Type1>
    BoxIntegrator2D<Type0, Type1>::
    BoxIntegrator2D()
      : m_cache(),
        m_corner0(0, 0)
    {
      // Empty.
    }


    // This constructor initializes the class instance, and then
    // passes its arguments to member function setArray() in order
    // to precompute integral information.
    template <class Type0, class Type1>
    BoxIntegrator2D<Type0, Type1>::
    BoxIntegrator2D(const Array2D<Type0>& inputArray)
      : m_cache(),
        m_corner0(0, 0)
    {
      this->setArray(inputArray);
    }


    // This constructor initializes the class instance, and then
    // passes its arguments to member function setArray() in order
    // to precompute integral information.
    template <class Type0, class Type1>
    template <class Functor>
    BoxIntegrator2D<Type0, Type1>::
    BoxIntegrator2D(const Array2D<Type0>& inputArray, Functor functor)
      : m_cache(),
        m_corner0(0, 0)
    {
      this->setArray(inputArray, functor);
    }


    // This constructor works just like the single argument version
    // of setArray, with the exception that the pre-integration is
    // performed over only a rectangular sub-array.
    template <class Type0, class Type1>
    BoxIntegrator2D<Type0, Type1>::
    BoxIntegrator2D(const Array2D<Type0>& inputArray,
                    const Index2D& corner0,
                    const Index2D& corner1)
      : m_cache(),
        m_corner0(0, 0)
    {
      this->setArray(inputArray, corner0, corner1);
    }


    // The copy constructor does a shallow copy of its argument,
    // however behavior of the class should be indistinguishable
    // from a deep copy (operations on each of the original and copy
    // will not affect the other).
    template <class Type0, class Type1>
    BoxIntegrator2D<Type0, Type1>::
    BoxIntegrator2D(const BoxIntegrator2D& other)
      : m_cache(other.m_cache),
        m_corner0(other.m_corner0)
    {
      // Empty.
    }


    // This member function returns the integral over the
    // rectangular region with corner0 and corner1 at its diagonally
    // opposite corners.
    template <class Type0, class Type1>
    Type1
    BoxIntegrator2D<Type0, Type1>::
    getIntegral(const Index2D& corner0, const Index2D& corner1)
    {
      return (m_cache(corner1.getRow(), corner1.getColumn())
              - m_cache(corner1.getRow(), corner0.getColumn())
              - m_cache(corner0.getRow(), corner1.getColumn())
              + m_cache(corner0.getRow(), corner0.getColumn()));
    }


    // This member function works just like member function
    // getIntegral(const Index2D&, const Index2D&), with two
    // exceptions: it is very slightly slower; and corner indexing
    // is slightly more sophisticated.
    template <class Type0, class Type1>
    Type1
    BoxIntegrator2D<Type0, Type1>::
    getIntegral(const Index2D& corner0, const Index2D& corner1,
                bool /* dummy */)
    {
      int row0 = corner0.getRow() - m_corner0.getRow();
      int row1 = corner1.getRow() - m_corner0.getRow();
      int column0 = corner0.getColumn() - m_corner0.getColumn();
      int column1 = corner1.getColumn() - m_corner0.getColumn();
      return (m_cache(row1, column1)
              - m_cache(row1, column0)
              - m_cache(row0, column1)
              + m_cache(row0, column0));
    }


    // This member function returns the raw 2D integral from which box
    // integration is performed.
    template <class Type0, class Type1>
    Type1
    BoxIntegrator2D<Type0, Type1>::
    getRawIntegral(size_t row, size_t column)
    {
      return m_cache(row, column);
    }


    // This member function works just like member function
    // getRawIntegral(size_t, size_t), but uses raster-indexing.
    template <class Type0, class Type1>
    Type1
    BoxIntegrator2D<Type0, Type1>::
    getRawIntegral(size_t index0)
    {
      return m_cache(index0);
    }


    // This member function discards and previously cached integral
    // information, performs a double integral over the input array,
    // and caches the result for future use in computing
    // sub-integrals.
    template <class Type0, class Type1>
    void
    BoxIntegrator2D<Type0, Type1>::
    setArray(const Array2D<Type0>& inputArray)
    {
      this->setArray(inputArray, common::StaticCastFunctor<Type0, Type1>());
    }


    // This member function discards and previously cached integral
    // information, applies the specified functor to each element of
    // the input array, performs a double integral over the input
    // array, and caches the result for future use in computing
    // sub-integrals.
    template <class Type0, class Type1>
    template <class Functor>
    void
    BoxIntegrator2D<Type0, Type1>::
    setArray(const Array2D<Type0>& inputArray, Functor functor)
    {
      m_corner0.setValue(0, 0);
      this->fillCache(
        inputArray.begin(), inputArray.rows(), inputArray.columns(),
        inputArray.columns(), functor);
    }


    // This constructor works just like the single argument version
    // of setArray, with the exception that the pre-integration is
    // performed over only a rectangular sub-array.
    template <class Type0, class Type1>
    void
    BoxIntegrator2D<Type0, Type1>::
    setArray(const Array2D<Type0>& inputArray,
             const Index2D& corner0,
             const Index2D& corner1)
    {
      this->setArray(inputArray, corner0, corner1,
                     common::StaticCastFunctor<Type0, Type1>());
    }


    // This member function works just like the three-argument
    // version of setArray, with the exception that the the
    // specified functor is applied to each element of the array
    // before integration.
    template <class Type0, class Type1>
    template <class Functor>
    void
    BoxIntegrator2D<Type0, Type1>::
    setArray(const Array2D<Type0>& inputArray,
             const Index2D& corner0,
             const Index2D& corner1,
             Functor functor)
    {
      int row0 = corner0.getRow();
      int row1 = corner1.getRow();
      int column0 = corner0.getColumn();
      int column1 = corner1.getColumn();
      if(row1 < row0) {
        std::swap(row0, row1);
      }
      if(column1 < column0) {
        std::swap(column0, column1);
      }

      int roiRows = row1 - row0;
      int roiColumns = column1 - column0;

      m_corner0.setValue(row0, column0);
      this->fillCache(
        inputArray.begin() + (row0 * inputArray.columns() + column0),
        roiRows, roiColumns, inputArray.columns(), functor);
    }


    // This constructor works just like the single argument version
    // of setArray, with the exception that the pre-integration is
    // performed over only a rectangular sub-array.
    template <class Type0, class Type1>
    template <class Functor>
    void
    BoxIntegrator2D<Type0, Type1>::
    fillCache(typename Array2D<Type0>::const_iterator inIter,
              int roiRows,
              int roiColumns,
              int inputArrayColumns,
              Functor functor)
    {
      int rowIncrement = inputArrayColumns - roiColumns;
      m_cache.reinit(roiRows + 1, roiColumns + 1);

      // First row of cache represents boxes with zero height, so
      // values are identically zero.
      std::fill(m_cache.rowBegin(0), m_cache.rowEnd(0),
                static_cast<Type1>(functor(0)));

      // Pre-integrate first row.
      typename Array2D<Type0>::const_iterator endIter = inIter + roiColumns;
      typename Array2D<Type1>::iterator outIter =
        m_cache.begin() + m_cache.columns();
      *(outIter++) = static_cast<Type1>(functor(0));
      while(inIter != endIter) {
        *outIter = *(outIter - 1) + static_cast<Type1>(functor(*inIter));
        ++outIter;
        ++inIter;
      }

      // Pre-integrate all remaining rows.
      typename Array2D<Type1>::iterator previousRowIter =
        m_cache.begin() + m_cache.columns();
      for(int row = 1; row < roiRows; ++row) {
        inIter += rowIncrement;
        Type1 rowSum = static_cast<Type1>(0);
        *(outIter++) = rowSum;
        ++previousRowIter;
        for(int column = 0; column < roiColumns; ++column) {
          rowSum += static_cast<Type1>(functor(*inIter));
          *outIter = rowSum + (*previousRowIter);
          ++inIter;
          ++outIter;
          ++previousRowIter;
        }
      }
    }

  } // namespace numeric

} // namespace brick

#endif /* #ifndef BRICK_NUMERIC_BOXINTEGRATOR2D_IMPL_HH */
