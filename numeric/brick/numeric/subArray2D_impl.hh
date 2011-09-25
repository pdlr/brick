/**
***************************************************************************
* @file brick/numeric/subArray2D.hh 
* Header file declaring SubArray2D class template.
*
* Copyright (C) 2001-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/


#ifndef BRICK_NUMERIC_SUBARRAY2D_IMPL_HH
#define BRICK_NUMERIC_SUBARRAY2D_IMPL_HH

// This file is included by subArray2D.hh, and should not be directly included
// by user code, so no need to include subArray2D.hh here.
// 
// #include <brick/numeric/subArray2D.hh>

#include <cmath>
#include <brick/common/mathFunctions.hh>   // For absoluteValue().
#include <brick/common/stridedPointer.hh>

namespace brick {

  namespace numeric {
    
    template <class Type>
    SubArray2D<Type>::
    SubArray2D(const Array2D<Type>& source)
      : m_source(source),
        m_startRow(0),
        m_stopRow(source.rows()),
        m_rowStride(1),
        m_rows(source.rows()),
        m_startColumn(0),
        m_stopColumn(source.columns()),
        m_columnStride(1),
        m_columns(source.columns())
    {
      // Empty
    }


    template <class Type>
    SubArray2D<Type>::
    SubArray2D(const Array2D<Type>& source, const Slice& rowSlice,
               const Slice& columnSlice)
      : m_source(source),
        m_startRow(rowSlice.start()),
        m_stopRow(rowSlice.stop()),
        m_rowStride(rowSlice.stride()),
        m_rows(0),
        m_startColumn(columnSlice.start()),
        m_stopColumn(columnSlice.stop()),
        m_columnStride(columnSlice.stride()),
        m_columns(0)
    {
      // It's convenient to be able to specify "last row/column" somehow.
      // Setting row1/column1 to zero will wrap to source.rows()/source.columns()
      // if appropriate.
      if((this->m_stopRow == 0) && (this->m_rowStride > 0)) {
        this->m_stopRow = static_cast<int>(source.rows());
      }
      if((this->m_stopColumn == 0) && (this->m_columnStride > 0)) {
        this->m_stopColumn = static_cast<int>(source.columns());
      }
    
      // Negative indexing is also super-convenient
      while(this->m_startRow < 0) {
        this->m_startRow += static_cast<int>(source.rows());
      }
      while(this->m_stopRow < 0) {
        this->m_stopRow += static_cast<int>(source.rows());
      }
      while(this->m_startColumn < 0) {
        this->m_startColumn += static_cast<int>(source.columns());
      }
      while(this->m_stopColumn < 0) {
        this->m_stopColumn += static_cast<int>(source.columns());
      }

      // Now compute sizes (yuck).
      this->m_rows = ((this->m_stopRow - this->m_startRow)
                      / this->m_rowStride); // integer division
      if(this->m_rows < 0) {
        this->m_rows = 0; 
      } else {
        // Can't think of a better way to do this.
        if(brick::common::absoluteValue(
             static_cast<int>(this->m_rows) * this->m_rowStride)
           < brick::common::absoluteValue(
             this->m_stopRow - this->m_startRow)) {
          ++this->m_rows;
        }
      }
      this->m_columns = ((this->m_stopColumn - this->m_startColumn)
                         / this->m_columnStride); // integer division
      if(this->m_columns < 0) {
        this->m_columns = 0; 
      } else {
        // Can't think of a better way to do this.
        if(brick::common::absoluteValue(
             static_cast<int>(this->m_columns) * this->m_columnStride)
           < brick::common::absoluteValue(
             this->m_stopColumn - this->m_startColumn)) {
          ++this->m_columns;
        }
      }
      // Make sure we won't read/write outside of the array.
      this->checkArray2DSize(source);
    }

    
    template <class Type>
    SubArray2D<Type>::
    SubArray2D(const SubArray2D<Type> &other)
      : m_source(other.m_source),
        m_startRow(other.m_startRow),
        m_stopRow(other.m_stopRow),
        m_rowStride(other.m_rowStride),
        m_rows(other.m_rows),
        m_startColumn(other.m_startColumn),
        m_stopColumn(other.m_stopColumn),
        m_columnStride(other.m_columnStride),
        m_columns(other.m_columns)
    {
      // Empty
    }


    template <class Type>
    SubArray2D<Type>::operator Array2D<Type>() const
    {
      Array2D<Type> returnVal(this->rows(), this->columns());
      subArray(returnVal) = *this;
      return returnVal;
    }


    template <class Type>
    SubArray2D<Type>& SubArray2D<Type>::
    operator=(const SubArray2D<Type>& other)
    {
      this->checkSubArray2DSize(other);
      // if(this->chooseMajorAxis() == 0) {
      if(1) {
        return this->copyRowMajor(other);
      } else {
        return this->copyColumnMajor(other);
      }
    }


    template <class Type>
    SubArray2D<Type>& SubArray2D<Type>::
    operator=(Type value)
    {
      int outRow = m_startRow;
      while(outRow < m_stopRow) {
        brick::common::StridedPointer<Type> outPtr0(
          this->m_source.data(outRow, this->m_startColumn),
                  this->m_columnStride);
        brick::common::StridedPointer<Type> outPtr1 =
          outPtr0 + this->columns();
        std::fill(outPtr0, outPtr1, value);
        outRow += this->m_rowStride;
      }
      return *this;
    }

    
    template <class Type>
    SubArray2D<Type>& SubArray2D<Type>::
    operator+=(const SubArray2D<Type>& other)
    {
      this->checkSubArray2DSize(other);
      int inRow = other.m_startRow;
      int outRow = this->m_startRow;
      while(inRow < other.m_stopRow) {
        brick::common::StridedPointer<Type> thisDataPtr0(
          this->m_source.data(outRow, this->m_startColumn),
                       this->m_columnStride);
        brick::common::StridedPointer<Type> thisDataPtr1 =
          thisDataPtr0 + this->columns();
        brick::common::StridedPointer<const Type> otherDataPtr0(
          other.m_source.data(inRow, other.m_startColumn),
          other.m_columnStride);
        std::transform(thisDataPtr0, thisDataPtr1, otherDataPtr0, thisDataPtr0,
                       std::plus<Type>());
        inRow += other.m_rowStride;
        outRow += this->m_rowStride;
      }
      return *this;
    }

    
    template <class Type>
    SubArray2D<Type>& SubArray2D<Type>::
    operator-=(const SubArray2D<Type>& other)
    {
      this->checkSubArray2DSize(other);
      int inRow = other.m_startRow;
      int outRow = this->m_startRow;
      while(inRow < other.m_stopRow) {
        brick::common::StridedPointer<Type> thisDataPtr0(
          this->m_source.data(outRow, this->m_startColumn),
          this->m_columnStride);
        brick::common::StridedPointer<Type> thisDataPtr1 =
          thisDataPtr0 + this->columns();
        brick::common::StridedPointer<const Type> otherDataPtr0(
          other.m_source.data(inRow, other.m_startColumn),
          other.m_columnStride);
        std::transform(thisDataPtr0, thisDataPtr1, otherDataPtr0, thisDataPtr0,
                       std::minus<Type>());
        inRow += other.m_rowStride;
        outRow += this->m_rowStride;
      }
      return *this;
    }

    
    template <class Type>
    inline void SubArray2D<Type>::
    checkArray2DSize(const Array2D<Type>& other) const
    {
#ifdef BRICK_NUMERIC_CHECKBOUNDS
      if((m_startRow < 0) || (m_startRow >= static_cast<int>(other.rows()))) {
        std::ostringstream message;
        message << "Invalid start row: " << m_startRow << std::endl;
        BRICK_THROW(IndexException, "SubArray2D::checkArray2DSize()",
                  message.str().c_str());
      }
      if(m_rowStride > 0) {
        if(m_stopRow > static_cast<int>(other.rows())) {
          std::ostringstream message;
          message << "Invalid stop row: " << m_stopRow << std::endl;
          BRICK_THROW(IndexException, "SubArray2D::checkArray2DSize()",
                    message.str().c_str());
        }
      } else if(m_rowStride < 0) {
        if(m_stopRow < -1) {
          std::ostringstream message;
          message << "Invalid stop row: " << m_stopRow << std::endl;
          BRICK_THROW(IndexException, "SubArray2D::checkArray2DSize()",
                    message.str().c_str());
        }
      } else {
        // m_rowStride == 0
        BRICK_THROW(IndexException, "SubArray2D::checkArray2DSize()",
                  "Invalid row stride: 0\n");
      }
      if((m_startColumn < 0)
         || (m_startColumn >= static_cast<int>(other.columns()))) {
        std::ostringstream message;
        message << "Invalid start column: " << m_startColumn << std::endl;
        BRICK_THROW(IndexException, "SubArray2D::checkArray2DSize()",
                  message.str().c_str());
      }
      if(m_columnStride > 0) {
        if(m_stopColumn > static_cast<int>(other.columns())) {
          std::ostringstream message;
          message << "Invalid stop column: " << m_stopColumn << std::endl;
          BRICK_THROW(IndexException, "SubArray2D::checkArray2DSize()",
                    message.str().c_str());
        }
      } else if(m_columnStride < 0) {
        if(m_stopColumn < -1) {
          std::ostringstream message;
          message << "Invalid stop column: " << m_stopColumn << std::endl;
          BRICK_THROW(IndexException, "SubArray2D::checkArray2DSize()",
                    message.str().c_str());
        }
      } else {
        // m_columnStride == 0
        BRICK_THROW(IndexException, "SubArray2D::checkArray2DSize()",
                  "Invalid column stride: 0\n");
      }
#endif
    }

    
    template <class Type>
    inline void SubArray2D<Type>::
    checkSubArray2DSize(const SubArray2D<Type>& other) const
    {
#ifdef BRICK_NUMERIC_CHECKBOUNDS
      if(other.rows() != this->rows()) {
        std::ostringstream message;
        message << "Row mismatch: " << other.rows() << " vs. "
                << this->rows() << std::endl;
        BRICK_THROW(IndexException, "SubArray2D::checkSubArray2DSize()",
                  message.str().c_str());
      }
      if(other.columns() != this->columns()) {
        std::ostringstream message;
        message << "Column mismatch: " << other.columns() << " vs. "
                << this->columns() << std::endl;
        BRICK_THROW(IndexException, "SubArray2D::checkSubArray2DSize()",
                  message.str().c_str());
      }
#endif
    }

    
    template <class Type>
    SubArray2D<Type>& SubArray2D<Type>::
    copyColumnMajor(const SubArray2D<Type>& other)
    {
      int inColumn = other.m_startColumn;
      int outColumn = this->m_startColumn;
      while(inColumn < other.m_stopColumn) {
        brick::common::StridedPointer<const Type> inPtr0(
          other.m_source.data(other.m_startRow, inColumn),
          other.m_rowStride * other.m_source.getRowStep());
        brick::common::StridedPointer<const Type> inPtr1 =
          inPtr0 + other.rows();
        brick::common::StridedPointer<Type> outPtr0(
          this->m_source.data(this->m_startRow, outColumn),
          this->m_rowStride * this->m_source.getRowStep());
        std::copy(inPtr0, inPtr1, outPtr0);
        inColumn += other.m_columnStride;
        outColumn += this->m_columnStride;
      }
      return *this;
    }
  
    template <class Type>
    SubArray2D<Type>& SubArray2D<Type>::
    copyRowMajor(const SubArray2D<Type>& other)
    {
      int inRow = other.m_startRow;
      int outRow = this->m_startRow;
      while(inRow < other.m_stopRow) {
        brick::common::StridedPointer<const Type> inPtr0(
          other.m_source.data(inRow, other.m_startColumn),
          other.m_columnStride);
        brick::common::StridedPointer<const Type> inPtr1 =
          inPtr0 + other.columns();
        brick::common::StridedPointer<Type> outPtr0(
          this->m_source.data(outRow, this->m_startColumn),
          this->m_columnStride);
        std::copy(inPtr0, inPtr1, outPtr0);
        inRow += other.m_rowStride;
        outRow += this->m_rowStride;
      }
      return *this;
    }

    /* Non-member functions */

    template <class Type>
    std::ostream& operator<<(std::ostream& stream,
                             const SubArray2D<Type>& subArray)
    {
      Array2D<Type> array(subArray.rows(), subArray.columns());
      subArray2D(array) = subArray;
      stream << "SubArray2D([[";
      for(int row = 0; row < array.rows(); ++row) {
        for(int column = 0; column < array.columns() - 1; ++column) {
          stream << array(row, column) << ", ";
        }
        stream << array(row, array.columns() - 1) << "],\n";
        if(row != array.rows() - 1) {
          stream << "          ";
        }
      }
      stream.flush();
      return stream;
    }

  } // namespace numeric

} // namespace brick

#endif /* #ifdef BRICK_NUMERIC_SUBARRAY2D_IMPL_HH */
