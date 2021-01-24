/**
***************************************************************************
* @file brick/numeric/stencil2D_impl.hh
*
* Header file defining the Stencil2D class template.
*
* Copyright (C) 2006, 2011 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_NUMERIC_STENCIL2D_IMPL_HH
#define BRICK_NUMERIC_STENCIL2D_IMPL_HH

// This file is included by stencil2D.hh, and should not be directly included
// by user code, so no need to include stencil2D.hh here.
//
// #include <brick/numeric/stencil2D.hh>

#include <cmath>

namespace brick {

  namespace numeric {

    template <class Type, int Size>
    Stencil2D<Type, Size>::
    Stencil2D()
      : m_basePtr(0),
        m_rows(0),
        m_columns(0),
        m_targetSize(0),
        m_numberOfElements(0),
        m_ptr(0),
        m_offsetArray(),
        m_incrementArray(),
        m_patternArray()
    {
      for(size_t elementIndex = 0; elementIndex < Size; ++elementIndex) {
        m_offsetArray[elementIndex] = 0;
        m_incrementArray[elementIndex] = 0;
        m_patternArray[elementIndex] = Index2D(0, 0);
      }
    }


    template <class Type, int Size>
    Stencil2D<Type, Size>::
    Stencil2D(size_t rows, size_t columns)
      : m_basePtr(0),
        m_rows(0),
        m_columns(0),
        m_targetSize(0),
        m_numberOfElements(rows * columns),
        m_ptr(0),
        m_offsetArray(),
        m_incrementArray(),
        m_patternArray()
    {
      if(rows * columns > Size) {
        std::ostringstream message;
        message << "Requested size, (" << rows << ", " << columns << ") "
                << "is too big for a " << Size << " element stencil.";
        BRICK_THROW(brick::common::ValueException, "Stencil2D::Stencil2D()",
                  message.str().c_str());
      }

      // Record the rows x columns ROI.
      size_t index0 = 0;
      for(int row = 0; row < static_cast<int>(rows); ++row) {
        for(int column = 0; column < static_cast<int>(columns); ++column) {
          m_offsetArray[index0] = 0;
          m_incrementArray[index0] = 0;
          m_patternArray[index0] = Index2D(row, column);
          ++index0;
        }
      }

      // Fill in any unused elements with zeros.
      while(index0 < Size) {
        m_offsetArray[index0] = 0;
        m_incrementArray[index0] = 0;
        m_patternArray[index0] = Index2D(0, 0);
        ++index0;
      }
    }


    template <class Type, int Size>
    Stencil2D<Type, Size>::
    Stencil2D(const Array2D<bool>& pattern)
      : m_basePtr(0),
        m_rows(0),
        m_columns(0),
        m_targetSize(0),
        m_numberOfElements(0),
        m_ptr(0),
        m_offsetArray(),
        m_incrementArray(),
        m_patternArray()
    {
      this->setPattern(pattern);
    }


    template <class Type, int Size>
    Stencil2D<Type, Size>::
    Stencil2D(const std::vector<Index2D>& pattern)
      : m_basePtr(0),
        m_rows(0),
        m_columns(0),
        m_targetSize(0),
        m_numberOfElements(pattern.size()),
        m_ptr(0),
        m_offsetArray(),
        m_incrementArray(),
        m_patternArray()
    {
      if(pattern.size() > Size) {
        std::ostringstream message;
        message << "Requested pattern has too many elements for a "
                << Size << " element stencil.";
        BRICK_THROW(brick::common::ValueException, "Stencil2D::Stencil2D()",
                  message.str().c_str());
      }
      std::fill(m_offsetArray.begin(), m_offsetArray.end(), 0);
      std::fill(m_incrementArray.begin(), m_incrementArray.end(), 0);
      std::copy(pattern.begin(), pattern.end(), m_patternArray.begin());
      for(size_t index0 = pattern.size(); index0 < Size; ++index0) {
        m_patternArray[index0] = Index2D(0, 0);
      }
    }


    template <class Type, int Size>
    Type&
    Stencil2D<Type, Size>::
    getReference(size_t elementNumber)
    {
      this->checkBounds(elementNumber);
      return *(m_ptr + m_offsetArray[elementNumber]);
    }


    template <class Type, int Size>
    Type
    Stencil2D<Type, Size>::
    getValue(size_t elementNumber) const
    {
      this->checkBounds(elementNumber);
      return *(m_ptr + m_offsetArray[elementNumber]);
    }


    template <class Type, int Size>
    void
    Stencil2D<Type, Size>::
    goTo(size_t index)
    {
      m_ptr = m_basePtr + index;
    }


    template <class Type, int Size>
    void
    Stencil2D<Type, Size>::
    goTo(size_t row, size_t column)
    {
      m_ptr = m_basePtr + row * m_columns + column;
    }


    template <class Type, int Size>
    void
    Stencil2D<Type, Size>::
    setPattern(const Array2D<bool>& pattern)
    {
      // Record which of the elements in pattern are true.
      size_t index0 = 0;
      for(int row = 0; row < static_cast<int>(pattern.rows()); ++row) {
        for(int column = 0; column < static_cast<int>(pattern.columns());
            ++column) {
          if(pattern(row, column)) {
            if(index0 >= Size) {
              std::ostringstream message;
              message << "Requested pattern has too many elements for a "
                      << Size << " element stencil.";
              BRICK_THROW(brick::common::ValueException,
                          "Stencil2D::Stencil2D()",
                          message.str().c_str());
            }
            m_offsetArray[index0] = 0;
            m_incrementArray[index0] = 0;
            m_patternArray[index0] = Index2D(row, column);
            ++index0;
          }
        }
      }
      m_numberOfElements = index0;

      // Fill in any unused elements with zeros.
      while(index0 < Size) {
        m_offsetArray[index0] = 0;
        m_incrementArray[index0] = 0;
        m_patternArray[index0] = Index2D(0, 0);
        ++index0;
      }
    }


    template <class Type, int Size> template <class Type2>
    void
    Stencil2D<Type, Size>::
    setTarget(const Array2D<Type2>& target)
    {
      m_basePtr = target.data();
      m_rows = target.rows();
      m_columns = target.columns();
      m_targetSize = static_cast<int>(target.size());
      m_ptr = m_basePtr;

      for(size_t elementIndex = 0; elementIndex < m_numberOfElements;
          ++elementIndex) {
        Index2D targetIndex = m_patternArray[elementIndex];
        m_offsetArray[elementIndex] =
          static_cast<int>(targetIndex.getColumn() + targetIndex.getRow() * m_columns);
      }
      for(size_t elementIndex = 0; elementIndex < m_numberOfElements - 1;
          ++elementIndex) {
        m_incrementArray[elementIndex] =
          m_offsetArray[elementIndex + 1] - m_offsetArray[elementIndex];
      }
      m_incrementArray[m_numberOfElements - 1] =
        m_targetSize - m_offsetArray[m_numberOfElements - 1];
    }


    template <class Type, int Size> template <class Type2>
    void
    Stencil2D<Type, Size>::
    setTarget(Array2D<Type2>& target)
    {
      m_basePtr = target.data();
      m_rows = target.rows();
      m_columns = target.columns();
      m_targetSize = target.size();
      m_ptr = m_basePtr;

      for(size_t elementIndex = 0; elementIndex < m_numberOfElements;
          ++elementIndex) {
        Index2D targetIndex = m_patternArray[elementIndex];
        m_offsetArray[elementIndex] =
          targetIndex.getColumn() + targetIndex.getRow() * m_columns;
      }
      for(size_t elementIndex = 0; elementIndex < m_numberOfElements - 1;
          ++elementIndex) {
        m_incrementArray[elementIndex] =
          m_offsetArray[elementIndex + 1] - m_offsetArray[elementIndex];
      }
      m_incrementArray[m_numberOfElements - 1] =
        m_targetSize - m_offsetArray[m_numberOfElements - 1];
    }


    template <class Type, int Size>
    inline void
    Stencil2D<Type, Size>::
    checkBounds(size_t
#ifdef BRICK_NUMERIC_CHECKBOUNDS
                elementNumber
#endif /* #ifdef BRICK_NUMERIC_CHECKBOUNDS */
      ) const
    {
#ifdef BRICK_NUMERIC_CHECKBOUNDS
      if(elementNumber >= m_numberOfElements) {
        std::ostringstream message;
        message << "Element number " << elementNumber
                << " is invalid for a(n) " << m_numberOfElements
                << " element stencil. ";
        BRICK_THROW(brick::common::IndexException,
                    "Stencil2D::checkBounds(size_t)",
                    message.str().c_str());
      }

      Type* resultPtr = m_ptr + this->m_offsetArray[elementNumber];
      if(resultPtr < m_basePtr) {
        BRICK_THROW(brick::common::IndexException,
                    "Stencil2D::checkBounds(size_t)",
                    "Pointer has retreated too far.  "
                    "Stencil no longer addresses valid memory.");
      }

      if(resultPtr >= m_basePtr + m_targetSize) {
        BRICK_THROW(brick::common::IndexException,
                    "Stencil2D::checkBounds(size_t)",
                    "Pointer has advanced too far.  "
                    "Stencil no longer addresses valid memory.");
      }
#endif /* #ifdef BRICK_NUMERIC_CHECKBOUNDS */
    }

  } // namespace numeric

} // namespace brick

#endif /* #ifndef BRICK_NUMERIC_STENCIL2D_IMPL_HH */
