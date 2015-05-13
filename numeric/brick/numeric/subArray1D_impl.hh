/**
***************************************************************************
* @file brick/numeric/subArray1D.hh
*
* Header file declaring SubArray1D class template.
*
* Copyright (C) 2001-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/


#ifndef BRICK_NUMERIC_SUBARRAY1D_IMPL_HH
#define BRICK_NUMERIC_SUBARRAY1D_IMPL_HH

// This file is included by subArray1D.hh, and should not be directly included
// by user code, so no need to include subArray1D.hh here.
// 
// #include <brick/numeric/subArray1D.hh>

#include <cmath>
#include <brick/common/exception.hh>
#include <brick/common/mathFunctions.hh>   // For absoluteValue().
#include <brick/common/stridedPointer.hh>

namespace brick {

  namespace numeric {
    
    template <class Type>
    SubArray1D<Type>::
    SubArray1D(const Array1D<Type>& source)
      : m_source(source),
        m_start(0),
        m_stop(source.size()),
        m_stride(1),
        m_size(source.size())
    {
      // Empty
    }

    
    template <class Type>
    SubArray1D<Type>::
    SubArray1D(const Array1D<Type>& source, const Slice& slice0)
      : m_source(source),
        m_start(slice0.start()),
        m_stop(slice0.stop()),
        m_stride(slice0.stride()),
        m_size(0)
    {
      // It's convenient to be able to specify "last element" somehow.
      // Setting the stop value to zero will wrap to source.size()
      // if appropriate.
      if((this->m_stop == 0) && (this->m_stride > 0)) {
        this->m_stop = source.size();
      }
    
      // Negative indexing is also super-convenient
      while(this->m_start < 0) {
        this->m_start += source.size();
      }
      while(this->m_stop < 0) {
        this->m_stop += source.size();
      }

      // Now compute sizes (yuck).
      int intSize = ((this->m_stop - this->m_start)
                     / this->m_stride); // integer division
      if(intSize < 0) {
        this->m_size = 0; 
      } else {
        this->m_size = intSize;
        // Can't think of a better way to do this.
        if(common::absoluteValue(static_cast<int>(this->m_size)
                                 * this->m_stride)
           < common::absoluteValue(this->m_stop - this->m_start)) {
          ++this->m_size;
        }
      }

      // Make sure we won't read/write outside of the array.
      this->checkArray1DSize(source);
    }

  
    template <class Type>
    SubArray1D<Type>::
    SubArray1D(const SubArray1D<Type> &other)
      : m_source(other.m_source),
        m_start(other.m_start),
        m_stop(other.m_stop),
        m_stride(other.m_stride),
        m_size(other.m_size)
    {
      // Empty
    }

  
    template <class Type>
    SubArray1D<Type>::operator Array1D<Type>() const
    {
      Array1D<Type> returnVal(this->size());
      subArray(returnVal) = *this;
      return returnVal;
    }

  
    template <class Type>
    SubArray1D<Type>& SubArray1D<Type>::
    operator=(const SubArray1D<Type>& other)
    {
      this->checkSubArray1DSize(other);
      return this->copy(other);
    }

  
    template <class Type>
    SubArray1D<Type>& SubArray1D<Type>::
    operator=(Type value)
    {
      brick::common::StridedPointer<Type> outPtr0(
        this->m_source.data(m_start), this->m_stride);
      brick::common::StridedPointer<Type> outPtr1 = outPtr0 + this->size();
      std::fill(outPtr0, outPtr1, value);
      return *this;
    }

  
    template <class Type>
    SubArray1D<Type>& SubArray1D<Type>::
    operator+=(const SubArray1D<Type>& other)
    {
      this->checkSubArray1DSize(other);
      brick::common::StridedPointer<Type> thisDataPtr0(
        this->m_source.data(this->m_start), this->m_stride);
      brick::common::StridedPointer<Type> thisDataPtr1 =
        thisDataPtr0 + this->size();
      brick::common::StridedPointer<const Type>
        otherDataPtr0(other.m_source.data(other.m_start), other.m_stride);
      std::transform(thisDataPtr0, thisDataPtr1, otherDataPtr0, thisDataPtr0,
                     std::plus<Type>());
      return *this;
    }

  
    template <class Type>
    SubArray1D<Type>& SubArray1D<Type>::
    operator-=(const SubArray1D<Type>& other)
    {
      this->checkSubArray1DSize(other);
      brick::common::StridedPointer<Type> thisDataPtr0(
        this->m_source.data(this->m_start), this->m_stride);
      brick::common::StridedPointer<Type> thisDataPtr1 =
        thisDataPtr0 + this->size();
      brick::common::StridedPointer<const Type>
        otherDataPtr0(other.m_source.data(other.m_start), other.m_stride);
      std::transform(thisDataPtr0, thisDataPtr1, otherDataPtr0, thisDataPtr0,
                     std::minus<Type>());
      return *this;
    }

    
    template <class Type>
    inline void SubArray1D<Type>::
    checkArray1DSize(const Array1D<Type>&
#ifdef BRICK_NUMERIC_CHECKBOUNDS
                     other
#endif /* #ifdef BRICK_NUMERIC_CHECKBOUNDS */
      ) const
    {
#ifdef BRICK_NUMERIC_CHECKBOUNDS
      if((m_start < 0) || (m_start >= static_cast<int>(other.size()))) {
        std::ostringstream message;
        message << "Invalid start index: " << m_start << std::endl;
        BRICK_THROW(brick::common::IndexException,
                    "SubArray1D::checkArray1DSize()",
                    message.str().c_str());
      }
      if(m_stride > 0) {
        if(m_stop > static_cast<int>(other.size())) {
          std::ostringstream message;
          message << "Invalid stop index: " << m_stop << std::endl;
          BRICK_THROW(brick::common::IndexException,
                      "SubArray1D::checkArray1DSize()",
                      message.str().c_str());
        }
      } else if(m_stride < 0) {
        if(m_stop < -1) {
          std::ostringstream message;
          message << "Invalid stop index: " << m_stop << std::endl;
          BRICK_THROW(brick::common::IndexException,
                      "SubArray1D::checkArray1DSize()",
                      message.str().c_str());
        }
      } else {
        // m_stride == 0
        BRICK_THROW(brick::common::IndexException,
                    "SubArray1D::checkArray1DSize()",
                    "Invalid stride: 0");
      }
#endif /* #ifdef BRICK_NUMERIC_CHECKBOUNDS */
    }

  
    template <class Type>
    inline void SubArray1D<Type>::
    checkSubArray1DSize(const SubArray1D<Type>&
#ifdef BRICK_NUMERIC_CHECKBOUNDS
                        other
#endif /* #ifdef BRICK_NUMERIC_CHECKBOUNDS */
      ) const
    {
#ifdef BRICK_NUMERIC_CHECKBOUNDS
      if(other.size() != this->size()) {
        std::ostringstream message;
        message << "Size mismatch: " << other.size() << " vs. "
                << this->size() << std::endl;
        BRICK_THROW(brick::common::IndexException,
                    "SubArray1D::checkSubArray1DSize()",
                    message.str().c_str());
      }
#endif /* #ifdef BRICK_NUMERIC_CHECKBOUNDS */
    }

  
    template <class Type>
    SubArray1D<Type>& SubArray1D<Type>::
    copy(const SubArray1D<Type>& other)
    {
      brick::common::StridedPointer<const Type>
        inPtr0(other.m_source.data(other.m_start), other.stride());
      brick::common::StridedPointer<const Type> inPtr1 = inPtr0 + other.size();
      brick::common::StridedPointer<Type>
        outPtr0(this->m_source.data(this->m_start), this->m_stride);
      std::copy(inPtr0, inPtr1, outPtr0);
      return *this;
    }


    /* Non-member functions */

    template <class Type>
    std::ostream& operator<<(std::ostream& stream,
                             const SubArray1D<Type>& subArray)
    {
      Array1D<Type> array(subArray.size());
      subArray1D(array) = subArray;
      stream << "SubArray1D([";
      for(int index0 = 0; index < array.size() - 1; ++index0) {
        stream << array(index0) << ", ";
      }
      stream << array(array.size() - 1) << "])";
      stream.flush();
      return stream;
    }

  } // namespace numeric

} // namespace brick

#endif /* #ifdef BRICK_NUMERIC_SUBARRAY1D_IMPL_HH */
