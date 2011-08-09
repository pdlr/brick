/**
***************************************************************************
* @file brick/numeric/array3D.hh
*
* Header file defining Array3D class template.
*
* Copyright (C) 2001-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_ARRAY3D_IMPL_HH
#define BRICK_NUMERIC_ARRAY3D_IMPL_HH

// This file is included by array3D.hh, and should not be directly included
// by user code, so no need to include array3D.hh here.
// 
// #include <brick/numeric/array3D.hh>

#include <algorithm>
#include <sstream>
#include <numeric>
#include <functional>
#include <brick/numeric/numericTraits.hh>

namespace brick {

  namespace numeric {
    
    // Static constant describing how the string representation of an
    // Array3D should start.
    template <class Type>
    const std::string&
    Array3D<Type>::
    ioIntro()
    {
      static const std::string intro = "Array3D(";
      return intro;
    }

    // Static constant describing how the string representation of an
    // Array3D should end.
    template <class Type>
    const char&
    Array3D<Type>::
    ioOutro()
    {
      static const char outro = ')';
      return outro;
    }


    // Static constant describing how the the data portion of the
    // string representation of an Array1D should start.
    template <class Type>
    const char&
    Array3D<Type>::
    ioOpening()
    {
      static const char opening = '[';
      return opening;
    }


    // Static constant describing how the the data portion of the
    // string representation of an Array3D should end.
    template <class Type>
    const char&
    Array3D<Type>::
    ioClosing()
    {
      static const char closing = ']';
      return closing;
    }


    // Static constant describing how individual elements should be
    // separated in the string representation of Array3D.
    template <class Type>
    const char&
    Array3D<Type>::
    ioSeparator()
    {
      static const char separator = ',';
      return separator;
    }


    // Non-static member functions below.

    template <class Type>
    Array3D<Type>::
    Array3D()
      : m_shape0(0),
        m_shape1(0),
        m_shape2(0),
        m_shape1Times2(0),
        m_size(0),
        m_dataPtr(0),
        m_refCountPtr(0),
        m_isAllocated(false)
    {
      // Empty.
    }

  
    template <class Type>
    Array3D<Type>::
    Array3D(size_t arrayShape0, size_t arrayShape1, size_t arrayShape2)
      : m_shape0(arrayShape0),
        m_shape1(arrayShape1),
        m_shape2(arrayShape2),
        m_shape1Times2(0), // This will be set in the call to allocate().
        m_size(0),         // This will be set in the call to allocate().
        m_dataPtr(0),      // This will be set in the call to allocate().
        m_refCountPtr(0),  // This will be set in the call to allocate().
        m_isAllocated(false)
    {
      this->allocate();
    }

  
    // Construct from an initialization string.
    template <class Type>
    Array3D<Type>::
    Array3D(const std::string& inputString)
      : m_shape0(0),
        m_shape1(0),
        m_shape2(0),
        m_shape1Times2(0),
        m_size(0),
        m_dataPtr(0),
        m_refCountPtr(0),
        m_isAllocated(false)
    {
      // We'll use the stream input operator to parse the string.
      std::istringstream inputStream(inputString);

      // Now read the string into an array.
      Array3D<Type> inputArray;
      inputStream >> inputArray;
      if(!inputStream) {
        std::ostringstream message;
        message << "Couldn't parse input string: \"" << inputString << "\".";
        BRICK_THROW(common::ValueException, "Array3D::Array3D(const std::string&)",
                   message.str().c_str());                 
      }

      // If all went well, copy into *this.
      *this = inputArray;
    }


  
    /* When copying from a Array3D do a shallow copy */
    /* Update reference count if the array we're copying has */
    /* valid data. */
    template <class Type>
    Array3D<Type>::
    Array3D(const Array3D<Type>& source)
      : m_shape0(source.m_shape0),
        m_shape1(source.m_shape1),
        m_shape2(source.m_shape2),
        m_shape1Times2(source.m_shape1 * source.m_shape2),
        m_size(source.m_size),
        m_dataPtr(source.m_dataPtr),
        m_refCountPtr(source.m_refCountPtr),
        m_isAllocated(source.m_isAllocated)
    {
      if(m_isAllocated) {
        ++(*m_refCountPtr);
      }
    }


    /* Here's a constructor for getting image data into the array */
    /* cheaply. */
    template <class Type>
    Array3D<Type>::
    Array3D(size_t arrayShape0, size_t arrayShape1, size_t arrayShape2, Type* const dataPtr)
      : m_shape0(arrayShape0),
        m_shape1(arrayShape1),
        m_shape2(arrayShape2),
        m_shape1Times2(arrayShape1 * arrayShape2),
        m_size(arrayShape0 * arrayShape1 * arrayShape2),
        m_dataPtr(dataPtr),
        m_refCountPtr(0),
        m_isAllocated(false)
    {
      // empty
    }

  
    template <class Type>
    Array3D<Type>::
    ~Array3D()
    {
      deAllocate();
    }


    template <class Type>
    inline void Array3D<Type>::
    checkDimension(
#ifdef BRICK_NUMERIC_CHECKBOUNDS
      size_t arrayShape0, size_t arrayShape1, size_t arrayShape2
#else /* #ifdef BRICK_NUMERIC_CHECKBOUNDS */
      size_t, size_t, size_t
#endif /* #ifdef BRICK_NUMERIC_CHECKBOUNDS [...] #else */
      ) const
    {
#ifdef BRICK_NUMERIC_CHECKBOUNDS
      if(arrayShape0 != this->shape0()
         || arrayShape1 != this->shape1()
         || arrayShape2 != this->shape2()){
        std::ostringstream message;
        message << "Size mismatch: required dimension is ("
                << arrayShape0 << ", " << arrayShape1 << ", " << arrayShape2 << ") "
                << " while *this has dimension "
                << this->shape0() << ", " << this->shape1() << ", "
                << this->shape2() << ") ";
        BRICK_THROW(common::IndexException, "Array3D::checkDimension()",
                  message.str().c_str());
      }
#endif
    }


    template <class Type>
    Array3D<Type> Array3D<Type>::
    copy() const
    {
      Array3D<Type> newArray(m_shape0, m_shape1, m_shape2);
      newArray.copy(*this);
      return newArray;
    }

  
    template <class Type> template <class Type2>
    void Array3D<Type>::
    copy(const Array3D<Type2>& source)
    {
      if(source.size() != m_size) {
        std::ostringstream message;
        message << "Mismatched array sizes. Source array has "
                << source.size() << " elements, while destination array has "
                << m_size << " elements.";
        BRICK_THROW(common::ValueException, "Array3D::copy(const Array3D&)",
                   message.str().c_str());
      }
      if(m_size != 0) {
        this->copy(source.data());
      }
    }

  
    template <class Type> template <class Type2>
    void Array3D<Type>::
    copy(const Type2* dataPtr)
    {
      if (dataPtr == 0) {
        BRICK_THROW(common::ValueException, "Array3D::copy(const Type2*)",
                  "Argument is a NULL pointer.");
      }
      std::copy(dataPtr, dataPtr + m_size, m_dataPtr);
    }

  
    // This member function sets the value of the array from an input
    // stream.
    template <class Type>
    std::istream&
    Array3D<Type>::
    readFromStream(std::istream& inputStream)
    {
      // Most of the time, InputType will be the same as Type.
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

        // We won't require the input format to start with "Array3D(", but
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
        Array2D<Type> inputValue;
        std::vector< Array2D<Type> > inputBuffer;
        while(1) {
          // Read the next row.
          inputStream >> inputValue;
          inputBuffer.push_back(inputValue);

          // Read the separator, or else the closing character.
          char inChar = 0;
          inputStream >> inChar;
          if(inChar == this->ioClosing()) {
            // Found a closing.  Stop here.
            break;
          }
          if(inChar != this->ioSeparator()) {
            // Missing separator?  Fail here.
            inputStream.clear(std::ios_base::failbit);
          }
        }
    
        // If we found an intro, we expect the corresponding outro.
        if(foundIntro) {
          inputStream >> common::Expect(&(ioOutro()), 1, flags);
        }

        // Now we're done with all of the parsing, verify that all slices
        // have the same number of rows and columns.
        size_t arrayShape0 = inputBuffer.size();
        size_t arrayShape1 = ((arrayShape0 != 0) ? inputBuffer[0].rows() : 0);
        size_t arrayShape2 = ((arrayShape0 != 0) ? inputBuffer[0].columns() : 0);
        for(size_t index = 1; index < arrayShape0; ++index) {
          if((inputBuffer[index].rows() != arrayShape1)
             || (inputBuffer[index].columns() != arrayShape2)) {
            // Inconsistent slice sizes!  Fail here.
            inputStream.clear(std::ios_base::failbit);
          }
        }

        // And finally, copy the data.
        size_t sliceSize = arrayShape1 * arrayShape2;
        this->reinit(arrayShape0, arrayShape1, arrayShape2);
        for(size_t index = 0; index < arrayShape0; ++index) {
          std::copy(inputBuffer[index].begin(), inputBuffer[index].end(),
                    this->begin() + (index * sliceSize));
        }
      } catch(std::ios_base::failure) {
        // Empty
      }
      inputStream.exceptions(oldExceptionState);
      return inputStream;
    }

  
    template <class Type>
    void Array3D<Type>::
    reinit(size_t arrayShape0, size_t arrayShape1, size_t arrayShape2)
    {
      this->deAllocate();
      this->m_shape0 = arrayShape0;
      this->m_shape1 = arrayShape1;
      this->m_shape2 = arrayShape2;
      this->allocate();
    }


    template <class Type>
    void Array3D<Type>::
    reinitIfNecessary(size_t arrayShape0, size_t arrayShape1,
                      size_t arrayShape2)
    {
      if(this->size() != arrayShape0 * arrayShape1 * arrayShape2) {
        this->reinit(arrayShape0, arrayShape1, arrayShape2);
      } else {
        if(this->shape0() != arrayShape0
           || this->shape1() != arrayShape1) {
          this->reshape(arrayShape0, arrayShape1, arrayShape2);
        }
      }
    }


    template <class Type>
    void Array3D<Type>::
    reshape(int arrayShape0, int arrayShape1, int arrayShape2)
    {
      if ((arrayShape0 * arrayShape1 * arrayShape2 != 0)){
        // If one axis is specified as -1, it will be automatically 
        // chosen to match the number of elements in the array.
        if((arrayShape0 == -1) && (arrayShape1 != -1) && (arrayShape2 != -1)) {
          arrayShape0 = static_cast<int>(this->size()) / (arrayShape1 * arrayShape2);
        } else if((arrayShape1 == -1) && (arrayShape0 != -1) && (arrayShape2 != -1)) {
          arrayShape1 = static_cast<int>(this->size()) / (arrayShape0 * arrayShape2);
        } else if((arrayShape2 == -1) && (arrayShape1 != -1) && (arrayShape0 != -1)) {
          arrayShape2 = static_cast<int>(this->size()) / (arrayShape1 * arrayShape0);
        }
      }
      if((arrayShape0 * arrayShape1 * arrayShape2) != static_cast<int>(this->size())) {
        std::ostringstream message;
        message << "Can't reshape a(n) " << this->size()
                << " element array to be " << arrayShape0 << " x " << arrayShape1
                << " x " << arrayShape2 << ".";
        BRICK_THROW(common::ValueException, "Array3D::reshape()", message.str().c_str());
      }

      m_shape0 = static_cast<size_t>(arrayShape0);
      m_shape1 = static_cast<size_t>(arrayShape1);
      m_shape2 = static_cast<size_t>(arrayShape2);
      m_shape1Times2 = arrayShape1 * arrayShape2;
    }

  
    template <class Type>
    Array1D<size_t> Array3D<Type>::
    shape() const
    {
      Array1D<size_t> rc(3);
      rc(0) = this->shape0();
      rc(1) = this->shape1();
      rc(2) = this->shape2();
      return rc;
    }


    template <class Type>
    size_t Array3D<Type>::
    shape(size_t axis) const
    {
      switch(axis) {
      case 0:
        return this->shape0();
        break;
      case 1:
        return this->shape1();
        break;
      case 2:
        return this->shape2();
        break;
      default:
        std::ostringstream message;
        message << "Invalid Axis: "<< axis << ".";
        BRICK_THROW(common::ValueException, "Array3D::shape(size_t)",
                  message.str().c_str());
        break;
      }
      return 0;                     // To keep the darn compiler happy.
    }


    template <class Type>
    Array2D<Type>
    Array3D<Type>::
    slice(size_t index0)
    {
      this->checkBounds(index0, 0, 0);
      return Array2D<Type>(
        m_shape1, m_shape2, m_dataPtr + (index0 * m_shape1Times2));
    }

  
    template <class Type>
    const Array2D<Type>
    Array3D<Type>::
    slice(size_t index0) const
    {
      this->checkBounds(index0, 0, 0);
      return Array2D<Type>(
        m_shape1, m_shape2, m_dataPtr + (index0 * m_shape1Times2));
    }

  
    template <class Type>
    Array3D<Type>& Array3D<Type>::
    operator=(const Array3D<Type>& source)
    {
      // Check for self-assignment
      if(&source != this) {
        this->deAllocate();
        m_shape0 = source.m_shape0;
        m_shape1 = source.m_shape1;
        m_shape2 = source.m_shape2;
        m_shape1Times2 = source.m_shape1Times2;
        m_size = source.m_size;
        m_dataPtr = source.m_dataPtr;
        m_refCountPtr = source.m_refCountPtr;
        m_isAllocated = source.m_isAllocated;
        if(m_isAllocated) {
          ++(*m_refCountPtr);
        }
      }
      return *this;
    }


    template <class Type>
    Array3D<Type>& Array3D<Type>::
    operator=(Type value)
    {
      std::fill(m_dataPtr, m_dataPtr + m_size, value);
      return *this;
    }


    template <class Type> template <class Type2>
    Array3D<Type>& Array3D<Type>::
    operator*=(const Array3D<Type2>& arg)
    {
      if(m_size != arg.size()) {
        std::ostringstream message;
        message << "Mismatched array sizes. Argument array has "
                << arg.size() << " elements, while destination array has "
                << m_size << " elements.";
        BRICK_THROW(common::ValueException, "Array3D::operator*=()",
                  message.str().c_str());
      }
      std::transform(m_dataPtr, m_dataPtr + m_size, arg.data(), m_dataPtr,
                     std::multiplies<Type>());
      return *this;
    }


    template <class Type> template <class Type2>
    Array3D<Type>& Array3D<Type>::
    operator/=(const Array3D<Type2>& arg)
    {
      if(m_size != arg.size()) {
        std::ostringstream message;
        message << "Mismatched array sizes. Argument array has "
                << arg.size() << " elements, while destination array has "
                << m_size << " elements.";
        BRICK_THROW(common::ValueException, "Array3D::operator/=()",
                  message.str().c_str());
      }
      std::transform(m_dataPtr, m_dataPtr + m_size, arg.data(), m_dataPtr,
                     std::divides<Type>());
      return *this;
    }


    template <class Type> template <class Type2>
    Array3D<Type>& Array3D<Type>::
    operator+=(const Array3D<Type2>& arg)
    {
      if(m_size != arg.size()) {
        std::ostringstream message;
        message << "Mismatched array sizes. Argument array has "
                << arg.size() << " elements, while destination array has "
                << m_size << " elements.";
        BRICK_THROW(common::ValueException, "Array3D::operator+=()",
                  message.str().c_str());
      }
      std::transform(m_dataPtr, m_dataPtr + m_size, arg.data(), m_dataPtr,
                     std::plus<Type>());
      return *this;
    }


    template <class Type> template <class Type2>
    Array3D<Type>& Array3D<Type>::
    operator-=(const Array3D<Type2>& arg)
    {
      if(m_size != arg.size()) {
        std::ostringstream message;
        message << "Mismatched array sizes. Argument array has "
                << arg.size() << " elements, while destination array has "
                << m_size << " elements.";
        BRICK_THROW(common::ValueException, "Array3D::operator-=()",
                  message.str().c_str());
      }
      std::transform(m_dataPtr, m_dataPtr + m_size, arg.data(), m_dataPtr,
                     std::minus<Type>());
      return *this;
    }


    template <class Type>
    Array3D<Type>& Array3D<Type>::
    operator+=(const Type arg)
    {
      std::transform(m_dataPtr, m_dataPtr + m_size, m_dataPtr,
                     std::bind2nd(std::plus<Type>(), arg));
      return *this;
    }


    template <class Type>
    Array3D<Type>& Array3D<Type>::
    operator-=(const Type arg)
    {
      std::transform(m_dataPtr, m_dataPtr + m_size, m_dataPtr,
                     std::bind2nd(std::minus<Type>(), arg));
      return *this;
    }


    template <class Type>
    Array3D<Type>& Array3D<Type>::
    operator*=(const Type arg)
    {
      std::transform(m_dataPtr, m_dataPtr + m_size, m_dataPtr,
                     std::bind2nd(std::multiplies<Type>(), arg));
      return *this;
    }


    template <class Type>
    Array3D<Type>& Array3D<Type>::
    operator/=(const Type arg)
    {
      std::transform(m_dataPtr, m_dataPtr + m_size, m_dataPtr,
                     std::bind2nd(std::divides<Type>(), arg));
      return *this;
    }


    template <class Type>
    void Array3D<Type>::
    allocate()
    {
      m_shape1Times2  = m_shape1 * m_shape2;
      m_size = m_shape0 * m_shape1 * m_shape2;
      if(m_shape0 > 0 && m_shape1 > 0 && m_shape2 > 0) {
        m_dataPtr = new(Type[m_size]); // should throw an exeption
        m_refCountPtr = new size_t;	 // if we're out of memory.
        *m_refCountPtr = 1;
        m_isAllocated = true;
        return;
      }
      m_dataPtr = 0;
      m_refCountPtr = 0;
      m_isAllocated = false;
      return;
    }


    template <class Type>
    inline void Array3D<Type>::
    checkBounds(
#ifdef BRICK_NUMERIC_CHECKBOUNDS
      size_t index
#else /* #ifdef BRICK_NUMERIC_CHECKBOUNDS */
      size_t
#endif /* #ifdef BRICK_NUMERIC_CHECKBOUNDS [...] #else */
      ) const
    {
#ifdef BRICK_NUMERIC_CHECKBOUNDS
      if(index >= m_size) {
        std::ostringstream message;
        message << "Index " << index << " is invalid for a(n) " << m_size
                << " element array.";
        BRICK_THROW(common::IndexException, "Array3D::checkBounds(size_t)",
                  message.str().c_str());
      }
#endif /* #ifdef BRICK_NUMERIC_CHECKBOUNDS */
    }

  
    template <class Type>
    inline void Array3D<Type>::
    checkBounds(
#ifdef BRICK_NUMERIC_CHECKBOUNDS
      size_t index0, size_t index1, size_t index2
#else /* #ifdef BRICK_NUMERIC_CHECKBOUNDS */
      size_t, size_t, size_t
#endif /* #ifdef BRICK_NUMERIC_CHECKBOUNDS [...] #else */
      ) const
    {
#ifdef BRICK_NUMERIC_CHECKBOUNDS
      if(index0 >= m_shape0) {
        std::ostringstream message;
        message << "index0 should be less than " << m_shape0
                << ", but is actually " << index0 << ".";
        BRICK_THROW(common::IndexException, "Array3D::checkBounds()",
                   message.str().c_str());
      }
      if(index1 >= m_shape1) {
        std::ostringstream message;
        message << "index1 should be less than " << m_shape1
                << ", but is actually " << index1 << ".";
        BRICK_THROW(common::IndexException, "Array3D::checkBounds()",
                   message.str().c_str());
      }
      if(index2 >= m_shape2) {
        std::ostringstream message;
        message << "index2 should be less than " << m_shape2
                << ", but is actually " << index2 << ".";
        BRICK_THROW(common::IndexException, "Array3D::checkBounds()",
                   message.str().c_str());
      }
#endif
    }

  
    template <class Type>
    void Array3D<Type>::
    deAllocate()
    {
      if(m_isAllocated == true) {
        if(--(*m_refCountPtr) == 0) {
          delete[] m_dataPtr;
          delete m_refCountPtr;
          m_isAllocated = false;
          m_dataPtr = 0;
          m_refCountPtr = 0;
        }
      } else {
        m_dataPtr = 0;
        m_refCountPtr = 0;
      }
    }

  
    /* Non-member functions which will ultimately wind up in a different file */
//   template <class Type> 
//   Type maximum(const Array3D<Type>& array)
//   {
//     return *std::max_element(array.data(), array.data() + array.size());
//   }

//   template <class Type> 
//   Type minimum(const Array3D<Type>& array)
//   {
//     return *std::min_element(array.data(), array.data() + array.size());
//   }

//   template <class Type> 
//   Type sum(const Array3D<Type>& array)
//   {
//     return std::accumulate(array.data(), array.data() + array.size(),
//                            static_cast<Type>(0));
//   }
  
    template <class Type>
    Array3D<Type> operator+(const Array3D<Type>& array0,
                            const Array3D<Type>& array1)
    {
      if((array0.shape0() != array1.shape0())
         || (array0.shape1() != array1.shape1())	
         || (array0.shape2() != array1.shape2())) {
        std::ostringstream message;
        message << "Array sizes do not match.  Array0 is "
                << array0.shape0() << " x " << array0.shape1()
                << " x " << array0.shape2()
                << ", while array1 is "
                << array1.shape0() << " x " << array1.shape1()
                << " x " << array1.shape2() << ".";
        BRICK_THROW(common::ValueException, "Array3D::operator+()", message.str().c_str());
      }
      Array3D<Type> result(array0.shape0(), array0.shape1(), array0.shape2());
      std::transform(array0.begin(), array0.end(), array1.begin(),
                     result.begin(), std::plus<Type>());
      return result;
    }


    template <class Type>
    Array3D<Type> operator-(const Array3D<Type>& array0,
                            const Array3D<Type>& array1)
    {
      if((array0.shape0() != array1.shape0())
         || (array0.shape1() != array1.shape1())	
         || (array0.shape2() != array1.shape2())) {
        std::ostringstream message;
        message << "Array sizes do not match.  Array0 is "
                << array0.shape0() << " x " << array0.shape1()
                << " x " << array0.shape2()
                << ", while array1 is "
                << array1.shape0() << " x " << array1.shape1()
                << " x " << array1.shape2() << ".";
        BRICK_THROW(common::ValueException, "Array3D::operator-()", message.str().c_str());
      }
      Array3D<Type> result(array0.shape0(), array0.shape1(), array0.shape2());
      std::transform(array0.begin(), array0.end(), array1.begin(),
                     result.begin(), std::minus<Type>());
      return result;
    }

  
    template <class Type>
    Array3D<Type> operator*(const Array3D<Type>& array0,
                            const Array3D<Type>& array1)
    {
      if((array0.shape0() != array1.shape0())
         || (array0.shape1() != array1.shape1())	
         || (array0.shape2() != array1.shape2())) {
        std::ostringstream message;
        message << "Array sizes do not match.  Array0 is "
                << array0.shape0() << " x " << array0.shape1()
                << " x " << array0.shape2()
                << ", while array1 is "
                << array1.shape0() << " x " << array1.shape1()
                << " x " << array1.shape2() << ".";
        BRICK_THROW(common::ValueException, "Array3D::operator*()", message.str().c_str());
      }
      Array3D<Type> result(array0.shape0(), array0.shape1(), array0.shape2());
      std::transform(array0.begin(), array0.end(), array1.begin(),
                     result.begin(), std::multiplies<Type>());
      return result;
    }
  

    template <class Type>
    Array3D<Type> operator/(const Array3D<Type>& array0,
                            const Array3D<Type>& array1)
    {
      if((array0.shape0() != array1.shape0())
         || (array0.shape1() != array1.shape1())	
         || (array0.shape2() != array1.shape2())) {
        std::ostringstream message;
        message << "Array sizes do not match.  Array0 is "
                << array0.shape0() << " x " << array0.shape1()
                << " x " << array0.shape2()
                << ", while array1 is "
                << array1.shape0() << " x " << array1.shape1()
                << " x " << array1.shape2() << ".";
        BRICK_THROW(common::ValueException, "Array3D::operator/()", message.str().c_str());
      }
      Array3D<Type> result(array0.shape0(), array0.shape1(), array0.shape2());
      std::transform(array0.begin(), array0.end(), array1.begin(), result.begin(),
                     std::divides<Type>());
      return result;
    }


    template <class Type>
    Array3D<Type> operator+(const Array3D<Type>& array0, Type scalar)
    {
      Array3D<Type> result(array0.shape0(), array0.shape1(), array0.shape2());
      std::transform(array0.begin(), array0.end(), result.begin(),
                     std::bind2nd(std::plus<Type>(), scalar));
      return result;
    }


    template <class Type>
    Array3D<Type> operator-(const Array3D<Type>& array0, Type scalar)
    {
      Array3D<Type> result(array0.shape0(), array0.shape1(), array0.shape2()); 
      std::transform(array0.begin(), array0.end(), result.begin(),
                     std::bind2nd(std::minus<Type>(), scalar));
      return result;
    }


    template <class Type>
    Array3D<Type> operator*(const Array3D<Type>& array0, Type scalar)
    {
      Array3D<Type> result(array0.shape0(), array0.shape1(), array0.shape2()); 
      std::transform(array0.begin(), array0.end(), result.begin(),
                     std::bind2nd(std::multiplies<Type>(), scalar));
      return result;
    }


    template <class Type>
    Array3D<Type> operator/(const Array3D<Type>& array0, Type scalar)
    {
      Array3D<Type> result(array0.shape0(), array0.shape1(), array0.shape2()); 
      std::transform(array0.begin(), array0.end(), result.begin(),
                     std::bind2nd(std::divides<Type>(), scalar));
      return result;
    }


    template <class Type>
    inline Array3D<Type> operator+(Type scalar, const Array3D<Type>& array0)
    {
      return array0 + scalar;
    }


    template <class Type>
    inline Array3D<Type> operator*(Type scalar, const Array3D<Type>& array0)
    {
      return array0 * scalar;
    }


    // Elementwise comparison of an Array3D with a constant.
    template <class Type>
    Array3D<bool>
    operator==(const Array3D<Type>& array0, const Type arg)
    {
      Array3D<bool> result(array0.shape0(), array0.shape1(), array0.shape2());
      std::transform(array0.begin(), array0.end(), result.data(),
                     std::bind2nd(std::equal_to<Type>(), arg));
      return result;
    }

    
    // Elementwise comparison of an Array3D with another array.
    template <class Type>
    Array3D<bool>
    operator==(const Array3D<Type>& array0, const Array3D<Type>& array1)
    {
      array0.checkDimension(array1.shape0(), array1.shape1(), array1.shape2());
      Array3D<bool> result(array0.shape0(), array0.shape1(), array0.shape2());
      std::transform(array0.begin(), array0.end(), array1.begin(),
                     result.begin(), std::equal_to<Type>());
      return result;
    }

    
    template <class Type>
    Array3D<bool> operator>(const Array3D<Type>& array0, Type arg)
    {
      Array3D<bool> result(array0.shape0(), array0.shape1(), array0.shape2()); 
      std::transform(array0.begin(), array0.end(), result.begin(),
                     std::bind2nd(std::greater<Type>(), arg));
      return result;
    }


    template <class Type>
    Array3D<bool> operator<(const Array3D<Type>& array0, Type arg)
    {
      Array3D<bool> result(array0.shape0(), array0.shape1(), array0.shape2()); 
      std::transform(array0.begin(), array0.end(), result.begin(),
                     std::bind2nd(std::less<Type>(), arg));
      return result;
    }

  
    template <class Type>
    Array3D<bool> operator>=(const Array3D<Type>& array0, Type arg)
    {
      Array3D<bool> result(array0.shape0(), array0.shape1(), array0.shape2()); 
      std::transform(array0.begin(), array0.end(), result.begin(),
                     std::bind2nd(std::greater_equal<Type>(), arg));
      return result;
    }

  
    template <class Type>
    Array3D<bool> operator<=(const Array3D<Type>& array0, Type arg)
    {
      Array3D<bool> result(array0.shape0(), array0.shape1(), array0.shape2()); 
      std::transform(array0.begin(), array0.end(), result.begin(),
                     std::bind2nd(std::less_equal<Type>(), arg));
      return result;
    }


    // This operator outputs a text representation of an Array3D
    // instance to a std::ostream.
    template <class Type>  
    std::ostream& operator<<(std::ostream& stream, const Array3D<Type>& array0)
    {
      // Most of the time, OutputType will be the same as Type.
      typedef typename NumericTraits<Type>::TextOutputType OutputType;

      size_t index0, index1, index2;
      stream << "Array3D([";
      for(index0 = 0; index0 < array0.shape0(); ++index0) {
        if(index0 != 0) {
          stream << "         ";
        }
        stream << "[";
        for(index1 = 0; index1 < array0.shape1(); ++index1) {
          if(index1 != 0) {
            stream << "          ";	
          }
          stream << "[";
          for(index2 = 0; index2 < array0.shape2(); ++index2) {
            stream << static_cast<OutputType>(array0(index0, index1, index2));
            if(index2 != array0.shape2() - 1) {
              stream << ", ";
            }
          }
          stream << "]";
          if(index1 != array0.shape1() - 1) {
            stream << ",\n";
          }	
        }
        stream << "]";
        if(index0 != array0.shape0() - 1) {
          stream << ",\n";
        }
      }
      stream << "])\n";
      stream.flush();
      return stream;
    }


    template <class Type>
    std::istream&
    operator>>(std::istream& inputStream, Array3D<Type>& array0)
    {
      return array0.readFromStream(inputStream);
    }
  
  } // namespace numeric

} // namespace brick

#endif /* #ifdef BRICK_NUMERIC_ARRAY3D_IMPL_HH */
