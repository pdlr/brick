/**
***************************************************************************
* @file brick/numeric/array2D.hh
*
* Header file defining Array2D class template.
*
* Copyright (C) 2001-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_ARRAY2D_IMPL_HH
#define BRICK_NUMERIC_ARRAY2D_IMPL_HH

// This file is included by array2D.hh, and should not be directly included
// by user code, so no need to include array2D.hh here.
// 
// #include <brick/numeric/array2D.hh>

#include <algorithm>
#include <functional>
#include <sstream>
#include <vector>
#include <brick/common/expect.hh>
#include <brick/common/functional.hh>
#include <brick/numeric/functional.hh>
#include <brick/numeric/numericTraits.hh>

namespace brick {

  namespace numeric {
    
    // Static constant describing how the string representation of an
    // Array2D should start.
    template <class Type>
    const std::string&
    Array2D<Type>::
    ioIntro()
    {
      static const std::string intro = "Array2D(";
      return intro;
    }

    // Static constant describing how the string representation of an
    // Array2D should end.
    template <class Type>
    const char&
    Array2D<Type>::
    ioOutro()
    {
      static const char outro = ')';
      return outro;
    }


    // Static constant describing how the the data portion of the
    // string representation of an Array1D should start.
    template <class Type>
    const char&
    Array2D<Type>::
    ioOpening()
    {
      static const char opening = '[';
      return opening;
    }


    // Static constant describing how the the data portion of the
    // string representation of an Array2D should end.
    template <class Type>
    const char&
    Array2D<Type>::
    ioClosing()
    {
      static const char closing = ']';
      return closing;
    }


    // Static constant describing how individual elements should be
    // separated in the string representation of Array2D.
    template <class Type>
    const char&
    Array2D<Type>::
    ioSeparator()
    {
      static const char separator = ',';
      return separator;
    }

    // Non-static member functions below.

    template <class Type>
    Array2D<Type>::
    Array2D()
      : m_rows(0),
        m_columns(0),
        m_rowStep(0),
        m_size(0),
        m_storageSize(0),
        m_dataPtr(0),
        m_referenceCount(0)
    {
      // Empty
    }


    template <class Type>
    Array2D<Type>::
    Array2D(size_t arrayRows, size_t arrayColumns, size_t rowStep)
      : m_rows(arrayRows),
        m_columns(arrayColumns),
        m_rowStep(rowStep ? rowStep : arrayColumns),
        m_size(0),           // This will be set in the call to allocate().
        m_storageSize(0),    // This will be set in the call to allocate().
        m_dataPtr(0),        // This will be set in the call to allocate().
        m_referenceCount(0)  // This will be set in the call to allocate().
    {
      this->allocate(arrayRows, arrayColumns, rowStep);
    }

  
    // Construct from an initialization string.
    template <class Type>
    Array2D<Type>::
    Array2D(const std::string& inputString, size_t rowStep)
      : m_rows(0),
        m_columns(0),
        m_rowStep(0),
        m_size(0),
        m_storageSize(0),
        m_dataPtr(0),
        m_referenceCount(0)
    {
      // We'll use the stream input operator to parse the string.
      std::istringstream inputStream(inputString);

      // Now read the string into an array.
      Array2D<Type> inputArray;
      inputStream >> inputArray;
      if(!inputStream) {
        std::ostringstream message;
        message << "Couldn't parse input string: \"" << inputString << "\".";
        BRICK_THROW(common::ValueException,
                    "Array2D::Array2D(const std::string&)",
                    message.str().c_str());                 
      }

      // If all went well, copy into *this.
      if(m_rowStep == 0) {
        *this = inputArray;
      } else {
        this->reinit(inputArray.rows(), inputArray.columns(), rowStep);
        this->copy(inputArray);
      }
    }


    /* When copying from a Array2D do a shallow copy */
    /* Update reference count if the array we're copying has */
    /* valid data. */
    template <class Type>
    Array2D<Type>::
    Array2D(const Array2D<Type>& source)
      : m_rows(source.m_rows),
        m_columns(source.m_columns),
        m_rowStep(source.m_rowStep),
        m_size(source.m_size),
        m_storageSize(source.m_storageSize),
        m_dataPtr(source.m_dataPtr),
        m_referenceCount(source.m_referenceCount)
    {
      // Empty.
    }

  
    /* Here's a constructor for getting image data into the array */
    /* cheaply. */
    template <class Type>
    Array2D<Type>::
    Array2D(size_t arrayRows, size_t arrayColumns, Type* const dataPtr,
            size_t rowStep)
      : m_rows(arrayRows),
        m_columns(arrayColumns),
        m_rowStep(rowStep ? rowStep : arrayColumns),
        m_size(arrayRows * arrayColumns),
        m_storageSize(arrayRows * m_rowStep),
        m_dataPtr(dataPtr),
        m_referenceCount(0)
    {
      // Empty
    }

  
    // Construct an array around external data that was allocated by
    // an Array?D instance.
    template <class Type>
    Array2D<Type>::
    Array2D(size_t arrayRows, size_t arrayColumns, Type* const dataPtr,
            common::ReferenceCount const& referenceCount)
      : m_rows(arrayRows),
        m_columns(arrayColumns),
        m_rowStep(arrayColumns),
        m_size(arrayRows * arrayColumns),
        m_storageSize(m_size),
        m_dataPtr(dataPtr),
        m_referenceCount(referenceCount)
    {
      // Empty.
    }


    // Construct an array around external data that was allocated by
    // an Array?D instance.
    template <class Type>
    Array2D<Type>::
    Array2D(size_t arrayRows, size_t arrayColumns, Type* const dataPtr,
            size_t rowStep, common::ReferenceCount const& referenceCount)
      : m_rows(arrayRows),
        m_columns(arrayColumns),
        m_rowStep(rowStep),
        m_size(arrayRows*arrayColumns),
        m_storageSize(arrayRows * m_rowStep),
        m_dataPtr(dataPtr),
        m_referenceCount(referenceCount)
    {
      // Empty.
    }
    
  
    template <class Type>
    Array2D<Type>::
    ~Array2D()
    {
      deAllocate();
    }

  
    template <class Type>
    inline void Array2D<Type>::
    checkDimension(size_t
#ifdef BRICK_NUMERIC_CHECKBOUNDS
                   arrayRows
#endif /* #ifdef BRICK_NUMERIC_CHECKBOUNDS */
                   , size_t
#ifdef BRICK_NUMERIC_CHECKBOUNDS
                   arrayColumns
#endif /* #ifdef BRICK_NUMERIC_CHECKBOUNDS */
      ) const
    {
#ifdef BRICK_NUMERIC_CHECKBOUNDS
      if(arrayRows != this->rows()
         || arrayColumns != this->columns()) {
        std::ostringstream message;
        message << "Size mismatch: required dimension is ("
                << arrayRows << ", " << arrayColumns << ") "
                << " while *this has dimension "
                << this->rows() << ", " << this->columns() << ").";
        BRICK_THROW(common::IndexException, "Array2D::checkDimension()",
                    message.str().c_str());
      }
#endif /* #ifdef BRICK_NUMERIC_CHECKBOUNDS */
    }


    template <class Type>
    Array2D<Type> Array2D<Type>::
    copy() const
    {
      Array2D<Type> newArray(m_rows, m_columns, m_rowStep);
      newArray.copy(*this);
      return newArray;
    }

    
    template <class Type> template <class Type2>
    void Array2D<Type>::
    copy(const Array2D<Type2>& source)
    {
      if(source.size() != m_size) {
        std::ostringstream message;
        message << "Mismatched array sizes. Source array has "
                << source.size() << " elements, while destination array has "
                << m_size << " elements.";
        BRICK_THROW(common::ValueException, "Array2D::copy(const Array2D&)",
                    message.str().c_str());
      }
      if(m_size != 0) {
        if(this->isContiguous() && source.isContiguous()) {
          this->copy(source.data());          
        } else {
          for(unsigned int rr = 0; rr < this->rows(); ++rr) {
            this->getRow(rr).copy(source.getRow(rr));
          }
        }
      }
    }


    template <class Type> template <class Type2>
    void Array2D<Type>::
    copy(const Type2* dataPtr)
    {
      if(dataPtr == 0) {
        BRICK_THROW(common::ValueException, "Array2D::copy(const Type2*)",
                    "Argument is a NULL pointer.");
      }
      if(this->isContiguous()) {
        std::transform(dataPtr, dataPtr + m_size, m_dataPtr,
                       StaticCastFunctor<Type2, Type>());        
      } else {
        for(unsigned int rr = 0; rr < this->rows(); ++rr) {
          this->getRow(rr).copy(dataPtr);
          dataPtr += m_columns;
        }
      }
    }


    // Return an Array2D instance referencing a subset of *this.
    template <class Type>
    Array2D<Type>
    Array2D<Type>::
    getRegion(Index2D const& corner0, Index2D const& corner1)
    {
      if(corner0.getRow()    >= corner1.getRow() ||
         corner0.getColumn() >= corner1.getColumn()) {
        BRICK_THROW(common::ValueException, "Array2D::getRegion()",
                    "Argument corner0 must be above and to the left of "
                    "corner1.  In other words, corner0.getRow() must return "
                    "a value smaller than is returned by corner1.getRow(), "
                    "and corner0.getRow() must return a value smaller than "
                    "is returned by corner1.getRow().");
      }

      if(corner0.getRow()    < 0                                 ||
         corner0.getColumn() < 0                                 ||
         corner0.getRow()    > static_cast<int>(this->m_rows)    ||
         corner0.getColumn() > static_cast<int>(this->m_columns) ||
         corner1.getRow()    < 0                                 ||
         corner1.getColumn() < 0                                 ||
         corner1.getRow()    > static_cast<int>(this->m_rows)    ||
         corner1.getColumn() > static_cast<int>(this->m_columns)) {
        BRICK_THROW(common::IndexException, "Array2D::getRegion()",
                    "Arguments corner0 and/or corner1 are out of bounds.");
      }

      unsigned int regionRows = corner1.getRow() - corner0.getRow();
      unsigned int regionColumns = corner1.getColumn() - corner0.getColumn();
      return Array2D<Type>(regionRows, regionColumns, this->getData(corner0),
                           this->getRowStep());
    }

    
    template <class Type>
    Array1D<Type> Array2D<Type>::
    getRow(size_t index)
    {
      this->checkBounds(index, 0);
      return Array1D<Type>(this->columns(),
                           m_dataPtr + (index * this->m_rowStep));
    }

  
    template <class Type>
    const Array1D<Type> Array2D<Type>::
    getRow(size_t index) const
    {
      this->checkBounds(index, 0);
      return Array1D<Type>(this->columns(),
                           m_dataPtr + (index * this->m_rowStep));
    }


    template <class Type>
    const Array1D<Type> Array2D<Type>::
    ravel() const
    {
      if(!(this->isContiguous())) {
        BRICK_THROW(common::StateException, "Array2D::ravel()",
                    "Can't flatten a non-contiguous array.");
      }
      if(this->isReferenceCounted()) {
        return Array1D<Type>(m_size, m_dataPtr, m_referenceCount);
      }
      return Array1D<Type>(m_size, m_dataPtr);
    }

    
    template <class Type>
    Array1D<Type> Array2D<Type>::
    ravel()
    {
      if(!(this->isContiguous())) {
        BRICK_THROW(common::StateException, "Array2D::ravel()",
                    "Can't flatten a non-contiguous array.");
      }
      if(this->isReferenceCounted()) {
        return Array1D<Type>(m_size, m_dataPtr, m_referenceCount);
      }
      return Array1D<Type>(m_size, m_dataPtr);
    }


    template <class Type>
    std::istream&
    Array2D<Type>::
    readFromStream(std::istream& inputStream)
    {
      // Most of the time, InputType will be the same as Type.
      // TextOutputType is the type one would use to write out a value
      // of Type.  Not surprisingly, this is the type you need to use
      // when reading those values back in.
      typedef typename NumericTraits<Type>::TextOutputType InputType;

      // If stream is in a bad state, we can't read from it.
      if (!inputStream){
        return inputStream;
      }
    
      // It's a lot easier to use a try block than to be constantly
      // testing whether the IO has succeeded, so we tell inputStream to
      // complain if anything goes wrong.
      std::ios_base::iostate oldExceptionState = inputStream.exceptions();
      inputStream.exceptions(
        std::ios_base::badbit | std::ios_base::failbit | std::ios_base::eofbit);

      // Now on with the show.
      try{
        common::Expect::FormatFlag flags = common::Expect::SkipWhitespace;

        // Skip any preceding whitespace.
        inputStream >> common::Expect("", flags);
              
        // We won't require the input format to start with "Array2D(", but
        // if it does we read it here.
        bool foundIntro = false;
        if(inputStream.peek() == ioIntro()[0]) {
          foundIntro = true;
          inputStream >> common::Expect(ioIntro(), flags);
        }

        // OK.  We've dispensed with the intro.  What's left should be of
        // the format "[row, row, row, ...]".  We require the square
        // brackets to be there.
        inputStream >> common::Expect(&(ioOpening()), 1, flags);

        // Read the data.  We'll use the Array1D<Type> stream operator to
        // read each row.
        Array1D<Type> inputValue;
        std::vector< Array1D<Type> > inputBuffer;
        while(1) {
          // Read the next row.
          inputStream >> inputValue;
          inputBuffer.push_back(inputValue);

          // Read the separator, or else the closing character.
          char inChar = 0;
          inputStream >> inChar;
          if(inChar == ioClosing()) {
            // Found a closing.  Stop here.
            break;
          }
          if(inChar != ioSeparator()) {
            // Missing separator?  Fail here.
            inputStream.clear(std::ios_base::failbit);
          }
        }
    
        // If we found an intro, we expect the corresponding outro.
        if(foundIntro) {
          inputStream >> common::Expect(&(ioOutro()), 1, flags);
        }

        // Now we're done with all of the parsing, verify that all rows
        // have the same length.
        size_t arrayRows = inputBuffer.size();
        size_t arrayColumns = ((inputBuffer.size() != 0)
                               ? inputBuffer[0].size() : 0);
        for(size_t index = 1; index < arrayRows; ++index) {
          if(inputBuffer[index].size() != arrayColumns) {
            // Inconsistent row lengths!  Fail here.
            inputStream.clear(std::ios_base::failbit);
          }
        }

        // And finally, copy the data.
        this->reinit(arrayRows, arrayColumns);
        for(size_t index = 0; index < arrayRows; ++index) {
          std::copy(inputBuffer[index].begin(), inputBuffer[index].end(),
                    this->begin() + (index * arrayColumns));
        }
      } catch(std::ios_base::failure) {
        // Empty
      }
      inputStream.exceptions(oldExceptionState);
      return inputStream;
    }
  

    template <class Type>
    void Array2D<Type>::
    reinit(size_t arrayRows, size_t arrayColumns, size_t rowStep)
    {
      this->allocate(arrayRows, arrayColumns, rowStep);
    }


    template <class Type>
    void Array2D<Type>::
    reinitIfNecessary(size_t arrayRows, size_t arrayColumns, size_t rowStep)
    {
      if(this->getStorageSize()
         != this->getStorageSize(arrayRows, arrayColumns, rowStep)) {
        this->reinit(arrayRows, arrayColumns, rowStep);
      } else {
        if((this->rows() != arrayRows)
           || (this->columns() != arrayColumns)) {
          this->reshape(arrayRows, arrayColumns, rowStep);
        }
      }
    }


    /* After reshaping, matrix is still row major order */
    template <class Type>
    void Array2D<Type>::
    reshape(int arrayRows, int arrayColumns, int rowStep)
    {
      // If one axis is specified as -1, it will be automatically 
      // chosen to match the number of elements in the array.
      if((arrayRows == -1) && (arrayColumns != 0)) {
        arrayRows = static_cast<int>(this->size()) / arrayColumns;
      } else if((arrayColumns == -1) && (arrayRows != 0)) {
        arrayColumns = static_cast<int>(this->size()) / arrayRows;
      }
      if(rowStep <= 0) {
        rowStep = arrayColumns;
      }
      if((arrayRows * arrayColumns) != static_cast<int>(this->size())) {
        std::ostringstream message;
        message << "Can't reshape a(n) " << this->size()
                << " element array to have " << arrayRows << " rows and "
                << arrayColumns << " columns.";
        BRICK_THROW(common::ValueException, "Array2D::reshape()",
                    message.str().c_str());
      }
      if((arrayRows * rowStep) != static_cast<int>(this->getStorageSize())) {
        std::ostringstream message;
        message << "Although the values for arrayRows (" << arrayRows
                << ") and arrayColumns (" << arrayColumns << ") are ok, "
                << "the specified value for rowStep (" << rowStep
                << ") doesn't match the the size of the available "
                << "storage (" << this->getStorageSize() << ").";
        BRICK_THROW(common::ValueException, "Array2D::reshape()",
                    message.str().c_str());
      }
      m_rows = arrayRows;
      m_columns = arrayColumns;
      m_rowStep = rowStep;
    }

  
    template <class Type>
    Array1D<size_t> Array2D<Type>::
    shape() const
    {
      Array1D<size_t> rc(2);
      rc(0) = this->rows();
      rc(1) = this->columns();
      return rc;
    }


    template <class Type>
    size_t Array2D<Type>::
    shape(size_t axis) const
    {
      size_t result;
      switch(axis) {
      case 0:
        result = this->rows();
        break;
      case 1:
        result = this->columns();
        break;
      default:
        std::ostringstream message;
        message << "Invalid Axis: "<< axis << ".";
        BRICK_THROW(common::ValueException, "Array2D::shape(size_t)",
                    message.str().c_str());
        result = 0;
        break;
      }
      return result;
    }


    template <class Type>
    Array2D<Type>& Array2D<Type>::
    operator=(Type value)
    {
      if(this->isContiguous()) {
        std::fill(m_dataPtr, m_dataPtr + m_size, value);
      } else {
        for(unsigned int rr = 0; rr < this->rows(); ++rr) {
          this->getRow(rr) = value;
        }
      }
      return *this;
    }

  
    template <class Type>
    Array2D<Type>& Array2D<Type>::
    operator=(const Array2D<Type>& source)
    {
      // Check for self-assignment
      if(&source != this) {
        this->deAllocate();
        m_rows = source.m_rows;
        m_columns = source.m_columns;
        m_rowStep = source.m_rowStep;
        m_size = source.m_size;
        m_storageSize = source.m_storageSize;
        m_dataPtr = source.m_dataPtr;
        m_referenceCount = source.m_referenceCount;
      }
      return *this;
    }


    template <class Type> template <class Type2>
    Array2D<Type>&
    Array2D<Type>::
    operator+=(const Array2D<Type2>& arg)
    {
      if(m_size != arg.size()) {
        std::ostringstream message;
        message << "Mismatched array sizes. Argument array is "
                << arg.rows() << " x " << arg.columns()
                << ", while destination array is "
                << m_rows << " x " << m_columns << ".";
        BRICK_THROW(common::ValueException, "Array2D::operator+=()",
                    message.str().c_str());
      }
      if(this->isContiguous()) {
        std::transform(m_dataPtr, m_dataPtr + m_size, arg.data(), m_dataPtr,
                       std::plus<Type>());
      } else {
        for(unsigned int rr = 0; rr < this->rows(); ++rr) {
          this->getRow(rr) += arg.getRow(rr);
        }
      }      
      return *this;
    }


    template <class Type> template <class Type2>
    Array2D<Type>&
    Array2D<Type>::
    operator-=(const Array2D<Type2>& arg)
    {
      if(m_size != arg.size()) {
        std::ostringstream message;
        message << "Mismatched array sizes. Argument array is "
                << arg.rows() << " x " << arg.columns()
                << ", while destination array is "
                << m_rows << " x " << m_columns << ".";
        BRICK_THROW(common::ValueException, "Array2D::operator-=()",
                    message.str().c_str());
      }
      if(this->isContiguous()) {
        std::transform(m_dataPtr, m_dataPtr + m_size, arg.data(), m_dataPtr,
                       std::minus<Type>());
      } else {
        for(unsigned int rr = 0; rr < this->rows(); ++rr) {
          this->getRow(rr) -= arg.getRow(rr);
        }
      }      
      return *this;
    }


    template <class Type> template <class Type2>
    Array2D<Type>&
    Array2D<Type>::
    operator*=(const Array2D<Type2>& arg)
    {
      if(m_size != arg.size()) {
        std::ostringstream message;
        message << "Mismatched array sizes. Argument array is "
                << arg.rows() << " x " << arg.columns()
                << ", while destination array is "
                << m_rows << " x " << m_columns << ".";
        BRICK_THROW(common::ValueException, "Array2D::operator*=()",
                    message.str().c_str());
      }
      if(this->isContiguous()) {
        std::transform(m_dataPtr, m_dataPtr + m_size, arg.data(), m_dataPtr,
                       std::multiplies<Type>());
      } else {
        for(unsigned int rr = 0; rr < this->rows(); ++rr) {
          this->getRow(rr) *= arg.getRow(rr);
        }
      }      
      return *this;
    }


    template <class Type> template <class Type2>
    Array2D<Type>&
    Array2D<Type>::
    operator/=(const Array2D<Type2>& arg)
    {
      if(m_size != arg.size()) {
        std::ostringstream message;
        message << "Mismatched array sizes. Argument array is "
                << arg.rows() << " x " << arg.columns()
                << ", while destination array is "
                << m_rows << " x " << m_columns << ".";
        BRICK_THROW(common::ValueException, "Array2D::operator/=()",
                    message.str().c_str());
      }
      if(this->isContiguous()) {
        std::transform(m_dataPtr, m_dataPtr + m_size, arg.data(), m_dataPtr,
                       std::divides<Type>());
      } else {
        for(unsigned int rr = 0; rr < this->rows(); ++rr) {
          this->getRow(rr) /= arg.getRow(rr);
        }
      }      
      return *this;
    }


    template <class Type>
    Array2D<Type>&
    Array2D<Type>::
    operator*=(Type arg)
    {
      if(this->isContiguous()) {
        std::transform(m_dataPtr, m_dataPtr + m_size, m_dataPtr,
                       std::bind2nd(std::multiplies<Type>(), arg));
      } else {
        for(unsigned int rr = 0; rr < this->rows(); ++rr) {
          this->getRow(rr) *= arg;
        }
      }      
      return *this;
    }


    template <class Type>
    Array2D<Type>&
    Array2D<Type>::
    operator/=(Type arg)
    {
      if(this->isContiguous()) {
        std::transform(m_dataPtr, m_dataPtr + m_size, m_dataPtr,
                       std::bind2nd(std::divides<Type>(), arg));
      } else {
        for(unsigned int rr = 0; rr < this->rows(); ++rr) {
          this->getRow(rr) /= arg;
        }
      }      
      return *this;
    }


    template <class Type>
    Array2D<Type>&
    Array2D<Type>::
    operator+=(Type arg)
    {
      if(this->isContiguous()) {
        std::transform(m_dataPtr, m_dataPtr + m_size, m_dataPtr,
                       std::bind2nd(std::plus<Type>(), arg));
      } else {
        for(unsigned int rr = 0; rr < this->rows(); ++rr) {
          this->getRow(rr) += arg;
        }
      }      
      return *this;
    }


    template <class Type>
    Array2D<Type>&
    Array2D<Type>::
    operator-=(Type arg)
    {
      if(this->isContiguous()) {
        std::transform(m_dataPtr, m_dataPtr + m_size, m_dataPtr,
                       std::bind2nd(std::minus<Type>(), arg));
      } else {
        for(unsigned int rr = 0; rr < this->rows(); ++rr) {
          this->getRow(rr) -= arg;
        }
      }      
      return *this;
    }


    template <class Type>
    Array2D<Type> Array2D<Type>::
    transpose() const
    {
      Array2D<Type> newMx(m_columns, m_rows);

      // Waiting for row & column iterators
      Type *tPtr0 = newMx.m_dataPtr;
      for(size_t j = 0; j < m_columns; ++j) {
        // const Type *tPtr1 = this->data(0, j);
        const Type *tPtr1 = this->data(j);
        for(size_t ii = 0; ii < m_rows; ++ii) {
          *tPtr0 = *tPtr1;
          ++tPtr0;
          tPtr1 += m_rowStep;
        }
      }
      return newMx;
    }


    template <class Type>
    void Array2D<Type>::
    allocate(size_t arrayRows, size_t arrayColumns, size_t rowStep)
    {
      // Make sure to release any shared resources.
      this->deAllocate();
      
      // Check array size.  It doesn't make sense to allocate memory
      // for a zero size array.
      m_rows = arrayRows;
      m_columns = arrayColumns;
      m_rowStep = rowStep != 0 ? rowStep : arrayColumns;
      m_size = m_rows * m_columns;
      m_storageSize = m_rows * m_rowStep;
      if(m_storageSize > 0) {
        // Allocate data storage.  new() should throw an exception if
        // we run out of memory.
        m_dataPtr = new(Type[m_storageSize]);

        // Set reference count to show that exactly one Array is pointing
        // to this data.
        m_referenceCount.reset(1);
      }
      return;
    }


    template <class Type>
    inline void Array2D<Type>::
    checkBounds(size_t
#ifdef BRICK_NUMERIC_CHECKBOUNDS
                index
#endif /* #ifdef BRICK_NUMERIC_CHECKBOUNDS */
      ) const
    {
#ifdef BRICK_NUMERIC_CHECKBOUNDS
      if(index >= m_size) {
        std::ostringstream message;
        message << "Index " << index << " is invalid for a(n) " << m_rows
                << " x " << m_columns << " array.";
        BRICK_THROW(common::IndexException, "Array2D::checkBounds(size_t)",
                    message.str().c_str());
      }
#endif /* #ifdef BRICK_NUMERIC_CHECKBOUNDS */
    }


    template <class Type>
    inline void Array2D<Type>::
    checkBounds(size_t
#ifdef BRICK_NUMERIC_CHECKBOUNDS
                rowIndex
#endif /* #ifdef BRICK_NUMERIC_CHECKBOUNDS */
                , size_t
#ifdef BRICK_NUMERIC_CHECKBOUNDS
                columnIndex
#endif /* #ifdef BRICK_NUMERIC_CHECKBOUNDS */
      ) const
    {
#ifdef BRICK_NUMERIC_CHECKBOUNDS
      if(rowIndex >= m_rows) {
        std::ostringstream message;
        message << "Row index " << rowIndex << " is invalid for a(n) "
                << m_rows << " x " << m_columns << " array.";
        BRICK_THROW(common::IndexException,
                    "Array2D::checkBounds(size_t, size_t)",
                    message.str().c_str());
      }
      if(columnIndex >= m_columns) {
        std::ostringstream message;
        message << "Column index " << columnIndex << " is invalid for a(n) "
                << m_rows << " x " << m_columns << " array.";
        BRICK_THROW(common::IndexException,
                    "Array2D::checkBounds(size_t, size_t)",
                    message.str().c_str());
      }
#endif
    }


    template <class Type>
    void Array2D<Type>::
    deAllocate()
    {
      // Are we responsible for deallocating the contents of this array?
      if(m_referenceCount.isCounted()) {
        // If yes, are we currently the only array pointing to this data?
        if(!m_referenceCount.isShared()) {
          // If yes, then delete the data.
          delete[] m_dataPtr;
        }
      }
      // Abandon our pointers to data.  Reference would take care of
      // itself, but it's cleaner conceptually to wipe it here, rather
      // than waiting for a subsequent call to ~ReferenceCount() or
      // ReferenceCount.allocate().
      m_dataPtr = 0;
      m_size = 0;
      m_storageSize = 0;
      m_rows = 0;
      m_columns = 0;
      m_referenceCount.reset(0);
    }

    /* ======== Non-member functions ======== */

    template <class Type>
    inline Array2D<Type>
    sqrt(const Array2D<Type>& array0)
    {
      return squareRoot(array0);
    }

    // This function returns an Array2D instance of the same shape and
    // element type as its input, in which each element contains the
    // square root of the corresponding element of the input array.
    template <class Type>
    Array2D<Type>
    squareRoot(const Array2D<Type>& array0)
    {
      Array2D<Type> result(array0.rows(), array0.columns(),
                           array0.getRowStep());
      if(array0.isContiguous()) {
        std::transform(array0.begin(), array0.end(),
                       result.begin(), SquareRootFunctor<Type>());
      } else {
        for(unsigned int row = 0; row < array0.rows(); ++row) {
          Array1D<Type> inRow = array0.getRow(row);
          Array1D<Type> outRow = result.getRow(row);
          std::transform(inRow.begin(), inRow.end(), outRow.begin(),
                         SquareRootFunctor<Type>());
        }
      }      
        
      return result;
    }

     
    template <class Type>
    Array2D<Type> operator+(const Array2D<Type>& array0,
                            const Array2D<Type>& array1)
    {
      if((array0.rows() != array1.rows())
         || (array0.columns() != array1.columns())) {
        std::ostringstream message;
        message << "Array sizes do not match.  Array0 is " << array0.rows()
                << " x " << array0.columns() << ", while array1 is "
                << array1.rows() << " x " << array1.columns() << ".";
        BRICK_THROW(common::ValueException, "Array2D::operator+()",
                    message.str().c_str());
      }

      Array2D<Type> result(array0.rows(), array0.columns(),
                           array0.getRowStep());
      if(array0.isContiguous() && array1.isContiguous()) {
        std::transform(array0.begin(), array0.end(), array1.begin(),
                       result.begin(), std::plus<Type>());
      } else {
        for(unsigned int row = 0; row < array0.rows(); ++row) {
          Array1D<Type> in0Row = array0.getRow(row);
          Array1D<Type> in1Row = array1.getRow(row);
          Array1D<Type> outRow = result.getRow(row);
          std::transform(in0Row.begin(), in0Row.end(), in1Row.begin(),
                         outRow.begin(), std::plus<Type>());
        }
      }
      return result;
    }
    
      
    template <class Type>
    Array2D<Type> operator-(const Array2D<Type>& array0,
                            const Array2D<Type>& array1)
    {
      if((array0.rows() != array1.rows())
         || (array0.columns() != array1.columns())) {
        std::ostringstream message;
        message << "Array sizes do not match.  Array0 is " << array0.rows()
                << " x " << array0.columns() << ", while array1 is "
                << array1.rows() << " x " << array1.columns() << ".";
        BRICK_THROW(common::ValueException, "Array2D::operator-()",
                    message.str().c_str());
      }
      Array2D<Type> result(array0.rows(), array0.columns(),
                           array0.getRowStep());
      if(array0.isContiguous() && array1.isContiguous()) {
        std::transform(array0.begin(), array0.end(), array1.begin(),
                       result.begin(), std::minus<Type>());
      } else {
        for(unsigned int row = 0; row < array0.rows(); ++row) {
          Array1D<Type> in0Row = array0.getRow(row);
          Array1D<Type> in1Row = array1.getRow(row);
          Array1D<Type> outRow = result.getRow(row);
          std::transform(in0Row.begin(), in0Row.end(), in1Row.begin(),
                         outRow.begin(), std::minus<Type>());
        }
      }
      return result;
    }


    template <class Type>
    Array2D<Type> operator*(const Array2D<Type>& array0,
                            const Array2D<Type>& array1)
    {
      if((array0.rows() != array1.rows())
         || (array0.columns() != array1.columns())) {
        std::ostringstream message;
        message << "Array sizes do not match.  Array0 is " << array0.rows()
                << " x " << array0.columns() << ", while array1 is "
                << array1.rows() << " x " << array1.columns() << ".";
        BRICK_THROW(common::ValueException, "Array2D::operator*()", message.str().c_str());
      }
      Array2D<Type> result(array0.rows(), array0.columns(),
                           array0.getRowStep());
      if(array0.isContiguous() && array1.isContiguous()) {
        std::transform(array0.begin(), array0.end(), array1.begin(),
                       result.begin(), std::multiplies<Type>());
      } else {
        for(unsigned int row = 0; row < array0.rows(); ++row) {
          Array1D<Type> in0Row = array0.getRow(row);
          Array1D<Type> in1Row = array1.getRow(row);
          Array1D<Type> outRow = result.getRow(row);
          std::transform(in0Row.begin(), in0Row.end(), in1Row.begin(),
                         outRow.begin(), std::multiplies<Type>());
        }
      }
      return result;
    }

    template <class Type>
    Array2D<Type> operator/(const Array2D<Type>& array0,
                            const Array2D<Type>& array1)
    {
      if((array0.rows() != array1.rows())
         || (array0.columns() != array1.columns())) {
        std::ostringstream message;
        message << "Array sizes do not match.  Array0 is " << array0.rows()
                << " x " << array0.columns() << ", while array1 is "
                << array1.rows() << " x " << array1.columns() << ".";
        BRICK_THROW(common::ValueException, "Array2D::operator/()",
                    message.str().c_str());
      }
      Array2D<Type> result(array0.rows(), array0.columns(),
                           array0.getRowStep());
      if(array0.isContiguous() && array1.isContiguous()) {
        std::transform(array0.begin(), array0.end(), array1.begin(),
                       result.begin(), std::divides<Type>());
      } else {
        for(unsigned int row = 0; row < array0.rows(); ++row) {
          Array1D<Type> in0Row = array0.getRow(row);
          Array1D<Type> in1Row = array1.getRow(row);
          Array1D<Type> outRow = result.getRow(row);
          std::transform(in0Row.begin(), in0Row.end(), in1Row.begin(),
                         outRow.begin(), std::divides<Type>());
        }
      }
      return result;
    }

    template <class Type>
    Array2D<Type> operator+(const Array2D<Type>& array0, Type scalar)
    {
      Array2D<Type> result(array0.rows(), array0.columns(),
                           array0.getRowStep());
      if(array0.isContiguous()) {
        std::transform(array0.begin(), array0.end(), result.begin(),
                       std::bind2nd(std::plus<Type>(), scalar));
      } else {
        for(unsigned int row = 0; row < array0.rows(); ++row) {
          Array1D<Type> inRow = array0.getRow(row);
          Array1D<Type> outRow = result.getRow(row);
          std::transform(inRow.begin(), inRow.end(), outRow.begin(), 
                         std::bind2nd(std::plus<Type>(), scalar));
        }
      }
      return result;
    }

    template <class Type>
    Array2D<Type> operator-(const Array2D<Type>& array0, Type scalar)
    {
      Array2D<Type> result(array0.rows(), array0.columns(),
                           array0.getRowStep());
      if(array0.isContiguous()) {
        std::transform(array0.begin(), array0.end(), result.begin(),
                       std::bind2nd(std::minus<Type>(), scalar));
      } else {
        for(unsigned int row = 0; row < array0.rows(); ++row) {
          Array1D<Type> inRow = array0.getRow(row);
          Array1D<Type> outRow = result.getRow(row);
          std::transform(inRow.begin(), inRow.end(), outRow.begin(), 
                         std::bind2nd(std::minus<Type>(), scalar));
        }
      }
      return result;
    }

    template <class Type>
    Array2D<Type> operator*(const Array2D<Type>& array0, Type scalar)
    {
      Array2D<Type> result(array0.rows(), array0.columns(),
                           array0.getRowStep());
      if(array0.isContiguous()) {
        std::transform(array0.begin(), array0.end(), result.begin(),
                       std::bind2nd(std::multiplies<Type>(), scalar));
      } else {
        for(unsigned int row = 0; row < array0.rows(); ++row) {
          Array1D<Type> inRow = array0.getRow(row);
          Array1D<Type> outRow = result.getRow(row);
          std::transform(inRow.begin(), inRow.end(), outRow.begin(), 
                         std::bind2nd(std::multiplies<Type>(), scalar));
        }
      }
      return result;
    }

    template <class Type>
    Array2D<Type> operator/(const Array2D<Type>& array0, Type scalar)
    {
      Array2D<Type> result(array0.rows(), array0.columns(),
                           array0.getRowStep());
      if(array0.isContiguous()) {
        std::transform(array0.begin(), array0.end(), result.begin(),
                       std::bind2nd(std::divides<Type>(), scalar));
      } else {
        for(unsigned int row = 0; row < array0.rows(); ++row) {
          Array1D<Type> inRow = array0.getRow(row);
          Array1D<Type> outRow = result.getRow(row);
          std::transform(inRow.begin(), inRow.end(), outRow.begin(), 
                         std::bind2nd(std::divides<Type>(), scalar));
        }
      }
      return result;
    }

    template <class Type>
    inline Array2D<Type> operator+(Type scalar, const Array2D<Type>& array0)
    {
      return array0 + scalar;
    }

    template <class Type>
    inline Array2D<Type> operator*(Type scalar, const Array2D<Type>& array0)
    {
      return array0 * scalar;
    }


    // Elementwise comparison of an Array2D with a constant.
    template <class Type>
    Array2D<bool>
    operator==(const Array2D<Type>& array0, const Type arg)
    {
      Array2D<bool> result(array0.rows(), array0.columns(),
                           array0.getRowStep());
      if(array0.isContiguous()) {
        std::transform(array0.begin(), array0.end(), result.data(),
                       std::bind2nd(std::equal_to<Type>(), arg));
      } else {
        for(unsigned int row = 0; row < array0.rows(); ++row) {
          Array1D<Type> inRow = array0.getRow(row);
          Array1D<bool> outRow = result.getRow(row);
          std::transform(inRow.begin(), inRow.end(), outRow.begin(), 
                         std::bind2nd(std::equal_to<Type>(), arg));
        }
      }
      return result;
    }

    
    // Elementwise comparison of an Array2D with another array.
    template <class Type>
    Array2D<bool>
    operator==(const Array2D<Type>& array0, const Array2D<Type>& array1)
    {
      array0.checkDimension(array1.rows(), array1.columns());
      Array2D<bool> result(array0.rows(), array0.columns(),
                           array0.getRowStep());
      if(array0.isContiguous() && array1.isContiguous()) {
        std::transform(array0.begin(), array0.end(), array1.begin(),
                       result.begin(), std::equal_to<Type>());
      } else {
        for(unsigned int row = 0; row < array0.rows(); ++row) {
          Array1D<Type> in0Row = array0.getRow(row);
          Array1D<Type> in1Row = array1.getRow(row);
          Array1D<bool> outRow = result.getRow(row);
          std::transform(in0Row.begin(), in0Row.end(), in1Row.begin(),
                         outRow.begin(), std::equal_to<Type>());
        }
      }
      return result;
    }

  
    template <class Type>
    Array2D<bool> operator>(const Array2D<Type>& array0, Type arg)
    {
      Array2D<bool> result(array0.rows(), array0.columns(),
                           array0.getRowStep());
      if(array0.isContiguous()) {
        std::transform(array0.begin(), array0.end(), result.begin(),
                       std::bind2nd(std::greater<Type>(), arg));
      } else {
        for(unsigned int row = 0; row < array0.rows(); ++row) {
          Array1D<Type> inRow = array0.getRow(row);
          Array1D<bool> outRow = result.getRow(row);
          std::transform(inRow.begin(), inRow.end(), outRow.begin(), 
                         std::bind2nd(std::greater<Type>(), arg));
        }
      }
      return result;
    }

    template <class Type>
    Array2D<bool> operator<(const Array2D<Type>& array0, Type arg)
    {
      Array2D<bool> result(array0.rows(), array0.columns(),
                           array0.getRowStep());
      if(array0.isContiguous()) {
        std::transform(array0.begin(), array0.end(), result.begin(),
                       std::bind2nd(std::less<Type>(), arg));
      } else {
        for(unsigned int row = 0; row < array0.rows(); ++row) {
          Array1D<Type> inRow = array0.getRow(row);
          Array1D<bool> outRow = result.getRow(row);
          std::transform(inRow.begin(), inRow.end(), outRow.begin(), 
                         std::bind2nd(std::less<Type>(), arg));
        }
      }
      return result;
    }

    template <class Type>
    Array2D<bool> operator>=(const Array2D<Type>& array0, Type arg)
    {
      Array2D<bool> result(array0.rows(), array0.columns(),
                           array0.getRowStep());
      if(array0.isContiguous()) {
        std::transform(array0.begin(), array0.end(), result.begin(),
                       std::bind2nd(std::greater_equal<Type>(), arg));
      } else {
        for(unsigned int row = 0; row < array0.rows(); ++row) {
          Array1D<Type> inRow = array0.getRow(row);
          Array1D<bool> outRow = result.getRow(row);
          std::transform(inRow.begin(), inRow.end(), outRow.begin(), 
                         std::bind2nd(std::greater_equal<Type>(), arg));
        }
      }
      return result;
    }

    template <class Type>
    Array2D<bool> operator<=(const Array2D<Type>& array0, Type arg)
    {
      Array2D<bool> result(array0.rows(), array0.columns(),
                           array0.getRowStep());
      if(array0.isContiguous()) {
        std::transform(array0.begin(), array0.end(), result.begin(),
                       std::bind2nd(std::less_equal<Type>(), arg));
      } else {
        for(unsigned int row = 0; row < array0.rows(); ++row) {
          Array1D<Type> inRow = array0.getRow(row);
          Array1D<bool> outRow = result.getRow(row);
          std::transform(inRow.begin(), inRow.end(), outRow.begin(), 
                         std::bind2nd(std::less_equal<Type>(), arg));
        }
      }
      return result;
    }

    template <class Type>
    std::ostream& operator<<(std::ostream& stream, const Array2D<Type>& array0)
    {
      // Most of the time, OutputType will be the same as Type.
      typedef typename NumericTraits<Type>::TextOutputType OutputType;

      stream << "Array2D([[";
      for(size_t row = 0; row < array0.rows(); ++row) {
        if (array0.columns() > 0) {
          for(size_t column = 0; column < array0.columns() - 1; ++column) {
            stream << static_cast<OutputType>(array0(row, column)) << ", ";
          }
          stream << static_cast<OutputType>(array0(row, array0.columns() - 1));
          if(row != array0.rows() - 1) {
            stream << "],\n";
            stream << "         [";
          }
        }
      }	
      stream << "]])";
      stream.flush();
      return stream;
    }

  
    template <class Type>
    std::istream& operator>>(std::istream& inputStream, Array2D<Type>& array0)
    {
      return array0.readFromStream(inputStream);
    }
  
  } // namespace numeric

} // namespace brick

#endif  /* #ifndef BRICK_NUMERIC_ARRAY2D_IMPL_HH */

