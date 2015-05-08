/**
***************************************************************************
* @file brick/numeric/bilinearInterpolator_impl.hh
*
* Header file defining inline and template functions for the
* BilinearInterpolator class.
*
* Copyright (C) 1999-2007,2011,2015 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_BILINEARINTERPOLATOR_IMPL_HH
#define BRICK_NUMERIC_BILINEARINTERPOLATOR_IMPL_HH

// This file is included by bilinearInterpolator.hh, and should not be
// directly included by user code, so no need to include
// bilinearInterpolator.hh here.
// 
// #include <brick/numeric/bilinearInterpolator.hh>

#include <cmath>

#include <brick/common/mathFunctions.hh>

namespace brick {

  namespace numeric {
    
    template <class TypeIn, class TypeOut, class FloatType>
    inline TypeOut
    BilinearInterpolator<TypeIn, TypeOut, FloatType>::
    operator()(FloatType row, FloatType column) const
    {
      this->checkBounds(row, column);

      FloatType temporaryFloat;
      FloatType x1;
      brick::common::splitFraction(column, temporaryFloat, x1);
      size_t i0 = static_cast<size_t>(temporaryFloat);
      
      FloatType y1;
      brick::common::splitFraction(row, temporaryFloat, y1);
      size_t j0 = static_cast<size_t>(temporaryFloat);
      
      FloatType x0 = 1.0 - x1;
      FloatType y0 = 1.0 - y1;
      size_t index0 = m_array.columns() * j0 + i0;
      return TypeOut(
        x0 * (y0 * static_cast<FloatType>(this->m_array(index0))
              + y1 * static_cast<FloatType>(
                this->m_array(index0 + m_array.columns())))
        + x1 * (y0 * static_cast<FloatType>(this->m_array(index0 + 1))
                + y1 * static_cast<FloatType>(
                  this->m_array(index0 + m_array.columns() + 1))));
    }


    template <class TypeIn, class TypeOut, class FloatType>
    inline void
    BilinearInterpolator<TypeIn, TypeOut, FloatType>::
    checkBounds(FloatType
#ifdef BRICK_NUMERIC_CHECKBOUNDS
                row
#endif /* #ifdef BRICK_NUMERIC_CHECKBOUNDS */
                , FloatType
#ifdef BRICK_NUMERIC_CHECKBOUNDS
                column
#endif /* #ifdef BRICK_NUMERIC_CHECKBOUNDS */
      ) const
    {
#ifdef BRICK_NUMERIC_CHECKBOUNDS
      int row0 = static_cast<int>(row);
      if(row0 < 0
         || (row0 == std::numeric_limits<int>::max())
         || (static_cast<FloatType>(row0 + 1)
             >= static_cast<FloatType>(m_array.rows()))) {
        std::ostringstream message;
        message << "Row index " << row << " is invalid for a(n) "
                << m_array.rows() << " x " << m_array.columns() << " array.";
        BRICK_THROW(brick::common::IndexException,
                    "BilinearInterpolator::checkBounds()",
                    message.str().c_str());
      }
      int column0 = static_cast<int>(column);
      if(column0 < 0
         || (column0 == std::numeric_limits<int>::max())
         || (static_cast<FloatType>(column0 + 1)
             >= static_cast<FloatType>(m_array.columns()))) {
        std::ostringstream message;
        message << "Column index " << column << " is too high "
                << "to interpolate within  a(n) "
                << m_array.rows() << " x " << m_array.columns() << " array.";
        BRICK_THROW(brick::common::IndexException,
                    "BilinearInterpolator::checkBounds()",
                    message.str().c_str());
      }
#endif /* #ifdef BRICK_NUMERIC_CHECKBOUNDS */
    }

  } // namespace numeric

} // namespace brick

#endif /* #ifndef BRICK_NUMERIC_BILINEARINTERPOLATOR_IMPL_HH */
