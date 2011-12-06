/**
***************************************************************************
* @file brick/numeric/bilinearInterpolator_impl.hh
*
* Header file defining inline and template functions for the
* BilinearInterpolator class.
*
* Copyright (C) 1999-2007,2011 David LaRose, dlr@cs.cmu.edu
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

namespace brick {

  namespace numeric {
    
    template <class TypeIn, class TypeOut>
    inline TypeOut
    BilinearInterpolator<TypeIn, TypeOut>::
    operator()(double row, double column)
    {
      this->checkBounds(row, column);
    
      double tmpDbl;
      double x1 = modf(column, &tmpDbl);
      size_t i0 = static_cast<size_t>(tmpDbl);
      double y1 = modf(row, &tmpDbl);
      size_t j0 = static_cast<size_t>(tmpDbl);
      double x0 = 1.0 - x1;
      double y0 = 1.0 - y1;
      size_t index0 = m_array.columns() * j0 + i0;
      return TypeOut(
        x0 * (y0 * this->m_array(index0)
              + y1 * this->m_array(index0 + m_array.columns()))
        + x1 * (y0 * this->m_array(index0 + 1)
                + y1 * this->m_array(index0 + m_array.columns() + 1)));
    }


    template <class TypeIn, class TypeOut>
    inline void
    BilinearInterpolator<TypeIn, TypeOut>::
    checkBounds(double row, double column) const
    {
#ifdef BRICK_NUMERIC_CHECKBOUNDS
      int row0 = static_cast<int>(row);
      if(row0 < 0 || (row0 + 1) >= static_cast<int>(m_array.rows())) {
        std::ostringstream message;
        message << "Row index " << row << " is invalid for a(n) "
                << m_array.rows() << " x " << m_array.columns() << " array.";
        BRICK_THROW(IndexException, "BilinearInterpolator::checkBounds()",
                  message.str().c_str());
      }
      int column0 = static_cast<int>(column);
      if(column0 < 0 || (column0 + 1) >= static_cast<int>(m_array.columns())) {
        std::ostringstream message;
        message << "Column index " << column << " is invalid for a(n) "
                << m_array.rows() << " x " << m_array.columns() << " array.";
        BRICK_THROW(IndexException, "BilinearInterpolator::checkBounds()",
                  message.str().c_str());
      }
#endif
    }

  } // namespace numeric

} // namespace brick

#endif /* #ifndef BRICK_NUMERIC_BILINEARINTERPOLATOR_IMPL_HH */
