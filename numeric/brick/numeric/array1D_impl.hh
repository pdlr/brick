/**
***************************************************************************
* @file brick/numeric/array1D_impl.hh
*
* Header file defining Array1D inline and template functions.
*
* Copyright (C) 2001-2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_ARRAY1D_IMPL_HH
#define BRICK_NUMERIC_ARRAY1D_IMPL_HH

// This file is included by array1D.hh, and should not be directly included
// by user code, so no need to include array1D.hh here.
// 
// #include <brick/numeric/array1D.hh>

#include <algorithm>
#include <sstream>
#include <vector>
#include <brick/common/expect.hh>
#include <brick/numeric/numericTraits.hh>

namespace brick {

  namespace numeric {
    
    // Static constant describing how the string representation of an
    // Array1D should start.
    template <class Type>
    const std::string&
    Array1D<Type>::
    ioIntro()
    {
      static const std::string intro = "Array1D(";
      return intro;
    }


    // Static constant describing how the string representation of an
    // Array1D should end.
    template <class Type>
    const char&
    Array1D<Type>::
    ioOutro()
    {
      static const char outro = ')';
      return outro;
    }


    // Static constant describing how the the data portion of the
    // string representation of an Array1D should start.
    template <class Type>
    const char&
    Array1D<Type>::
    ioOpening()
    {
      static const char opening = '[';
      return opening;
    }


    // Static constant describing how the the data portion of the
    // string representation of an Array1D should end.
    template <class Type>
    const char&
    Array1D<Type>::
    ioClosing()
    {
      static const char closing = ']';
      return closing;
    }


    // Static constant describing how individual elements should be
    // separated in the string representation of Array1D.
    template <class Type>
    const char&
    Array1D<Type>::
    ioSeparator()
    {
      static const char separator = ',';
      return separator;
    }

    
    // Non-static member functions below.

    template <class Type>
    Array1D<Type>::
    Array1D()
      : m_size(0),
        m_dataPtr(0),
        m_referenceCount(0)
    {
      // Empty.
    }

  
    template <class Type>
    Array1D<Type>::
    Array1D(size_t arraySize)
      : m_size(0),          // This will be set in the call to allocate().
        m_dataPtr(0),       // This will be set in the call to allocate().
        m_referenceCount(0) // This will be set in the call to allocate().
    {
      this->allocate(arraySize);
    }


    // Construct from an initialization string.
    template <class Type>
    Array1D<Type>::
    Array1D(const std::string& inputString)
      : m_size(0),
        m_dataPtr(0),
        m_referenceCount(0)
    {
      // We'll use the stream input operator to parse the string.
      std::istringstream inputStream(inputString);

      // Now read the string into an array.
      Array1D<Type> inputArray;
      inputStream >> inputArray;
      if(!inputStream) {
        std::ostringstream message;
        message << "Couldn't parse input string: \"" << inputString << "\".";
        BRICK_THROW(common::ValueException,
                    "Array1D::Array1D(const std::string&)",
                    message.str().c_str());                 
      }

      // If all went well, copy into *this.
      *this = inputArray;
    }

  
    /* When copying from an Array1D do a shallow copy */
    /* Update reference count if the array we're copying has */
    /* valid data. */
    template <class Type>
    Array1D<Type>::
    Array1D(const Array1D<Type>& source)
      : m_size(source.m_size),
        m_dataPtr(source.m_dataPtr),
        m_referenceCount(source.m_referenceCount)
    {
      // Empty.
    }


    template <class Type>
    Array1D<Type>::
    Array1D(size_t arraySize, Type* const dataPtr)
      : m_size(arraySize),
        m_dataPtr(dataPtr),
        m_referenceCount(0)
    {
      // Empty.
    }


    template <class Type>
    Array1D<Type>::
    Array1D(size_t arraySize, Type* const dataPtr,
            common::ReferenceCount const& referenceCount)
      : m_size(arraySize),
        m_dataPtr(dataPtr),
        m_referenceCount(referenceCount)
    {
      // Empty.
    }


    template <class Type>
    Array1D<Type>::
    Array1D(std::initializer_list<Type> initializer)
      : m_size(0),          // This will be set in the call to allocate().
        m_dataPtr(0),       // This will be set in the call to allocate().
        m_referenceCount(0) // This will be set in the call to allocate().
    {
      this->allocate(initializer.size());
      std::copy(initializer.begin(), initializer.end(), this->begin());
    }

    
    template <class Type>
    Array1D<Type>::~Array1D()
    {
      deAllocate();
    }


    template <class Type>
    inline void Array1D<Type>::
    checkDimension(size_t
#ifdef BRICK_NUMERIC_CHECKBOUNDS
                   arraySize
#endif /* #ifdef BRICK_NUMERIC_CHECKBOUNDS */
      ) const
    {
#ifdef BRICK_NUMERIC_CHECKBOUNDS
      if(arraySize != this->size()) {
        std::ostringstream message;
        message << "Size mismatch: required size is " << arraySize
                << " while *this has dimension " << this->size() << ".";
        BRICK_THROW(common::IndexException, "Array1D::checkDimension()",
                  message.str().c_str());
      }
#endif /* #ifdef BRICK_NUMERIC_CHECKBOUNDS */
    }


    template <class Type>
    Array1D<Type> Array1D<Type>::
    copy() const
    {
      Array1D<Type> newArray(m_size);
      newArray.copy(*this);
      return newArray;
    }

  
    template <class Type> template <class Type2>
    void Array1D<Type>::
    copy(const Array1D<Type2>& source)
    {
      if(source.size() != m_size) {
        std::ostringstream message;
        message << "Mismatched array sizes. Source array has "
                << source.size() << " elements, while destination array has "
                << m_size << " elements.";
        BRICK_THROW(common::ValueException, "Array1D::copy(const Array1D&)",
                    message.str().c_str());
      }
      if(m_size != 0) {
        this->copy(source.getData());
      }
    }

  
    template <class Type> template <class Type2>
    void
    Array1D<Type>::
    copy(const Type2* dataPtr)
    {
      if(dataPtr == 0) {
        BRICK_THROW(common::ValueException, "Array1D::copy(const Type2*)",
                    "Argument is a NULL pointer.");
      }
      std::copy(dataPtr, dataPtr + m_size, m_dataPtr);
    }


    template <class Type>
    Array1D<Type>& Array1D<Type>::
    operator=(Type val)
    {
      std::fill(m_dataPtr, m_dataPtr + m_size, val);
      return *this;
    }


    // This member function sets the value of the array from an input
    // stream.
    template <class Type>
    std::istream&
    Array1D<Type>::
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
        common::Expect::FormatFlag flags = common::Expect::SkipWhitespace();
      
        // Skip any preceding whitespace.
        inputStream >> common::Expect("", flags);

        // We won't require the input format to start with "Array1D(", but
        // if it does we read it here.
        bool foundIntro = false;
        if(inputStream.peek() == ioIntro()[0]) {
          foundIntro = true;
          inputStream >> common::Expect(ioIntro(), flags);
        }

        // OK.  We've dispensed with the intro.  What's left should be of
        // the format "[#, #, #, ...]".  We require the square brackets to
        // be there.
        inputStream >> common::Expect(&(ioOpening()), 1, flags);

        // Read the data.
        InputType inputValue;
        std::vector<Type> inputBuffer;
        while(1) {
          // Read the next value.
          inputStream >> inputValue;
          inputBuffer.push_back(static_cast<Type>(inputValue));

          // Read the separator, or else the closing character.
          char inChar = 0;
          inputStream >> inChar;
          if(inChar == ioClosing()) {
            // Found a closing.  Stop here.
            break;
          }
          if(inChar != ioSeparator()) {
            // Missing separator.  Fail here.
            inputStream.clear(std::ios_base::failbit);
          }
        }

        // If we found an intro, we expect the corresponding outro.
        if(foundIntro) {
          inputStream >> common::Expect(&(ioOutro()), 1, flags);
        }

        // Now we're done with all of the parsing.  Copy the data to *this.
        this->reinit(inputBuffer.size());
        std::copy(inputBuffer.begin(), inputBuffer.end(), this->begin());

      } catch(std::ios_base::failure) {
        // Empty
      }
      inputStream.exceptions(oldExceptionState);
      return inputStream;
    }
  

    template <class Type>
    void Array1D<Type>::
    reinit(size_t arraySize)
    {
      this->allocate(arraySize);
    }
  

    template <class Type>
    void Array1D<Type>::
    reinitIfNecessary(size_t arraySize)
    {
      if(this->size() != arraySize) {
        this->reinit(arraySize);
      }
    }
  

    template <class Type>
    Array1D<Type>& Array1D<Type>::
    operator=(const Array1D<Type>& source)
    {
      // Check for self-assignment
      if(&source != this) {
        this->deAllocate();
        m_size = source.m_size;
        m_dataPtr = source.m_dataPtr;
        m_referenceCount = source.m_referenceCount;
      }
      return *this;
    }


    template <class Type> template <class Type2>
    Array1D<Type>&
    Array1D<Type>::
    operator+=(const Array1D<Type2>& arg)
    {
      if(m_size != arg.size()) {
        std::ostringstream message;
        message << "Mismatched array sizes. Argument array has "
                << arg.size() << " elements, while destination array has "
                << m_size << " elements.";
        BRICK_THROW(common::ValueException, "Array1D::operator+=(const Array1D&)",
                     message.str().c_str());
      }
      std::transform(m_dataPtr, m_dataPtr + m_size, arg.data(), m_dataPtr,
                     std::plus<Type>());
      return *this;
    }


    template <class Type>
    Array1D<Type>&
    Array1D<Type>::
    operator+=(const Type arg)
    {
      std::transform(m_dataPtr, m_dataPtr + m_size, m_dataPtr,
                     std::bind2nd(std::plus<Type>(), arg));
      return *this;
    }


    template <class Type> template <class Type2>
    Array1D<Type>&
    Array1D<Type>::
    operator-=(const Array1D<Type2>& arg)
    {
      if(m_size != arg.size()) {
        std::ostringstream message;
        message << "Mismatched array sizes. Argument array has "
                << arg.size() << " elements, while destination array has "
                << m_size << " elements.";
        BRICK_THROW(common::ValueException,
                    "Array1D::operator-=(const Array1D&)",
                    message.str().c_str());
      }
      std::transform(m_dataPtr, m_dataPtr + m_size, arg.data(), m_dataPtr,
                     std::minus<Type>());
      return *this;
    }


    template <class Type>
    Array1D<Type>&
    Array1D<Type>::
    operator-=(const Type arg)
    {
      std::transform(m_dataPtr, m_dataPtr + m_size, m_dataPtr,
                     std::bind2nd(std::minus<Type>(), arg));
      return *this;
    }


    template <class Type> template <class Type2>
    Array1D<Type>&
    Array1D<Type>::
    operator*=(const Array1D<Type2>& arg)
    {
      if(m_size != arg.size()) {
        std::ostringstream message;
        message << "Mismatched array sizes. Argument array has "
                << arg.size() << " elements, while destination array has "
                << m_size << " elements.";
        BRICK_THROW(common::ValueException, "Array1D::operator*=(const Array1D&)",
                     message.str().c_str());
      }
      std::transform(m_dataPtr, m_dataPtr + m_size, arg.data(), m_dataPtr,
                     std::multiplies<Type>());
      return *this;
    }


    template <class Type>
    Array1D<Type>&
    Array1D<Type>::
    operator*=(const Type arg)
    {
      std::transform(m_dataPtr, m_dataPtr + m_size, m_dataPtr,
                     std::bind2nd(std::multiplies<Type>(), arg));
      return *this;
    }


    template <class Type> template <class Type2>
    Array1D<Type>&
    Array1D<Type>::
    operator/=(const Array1D<Type2>& arg)
    {
      if(m_size != arg.size()) {
        std::ostringstream message;
        message << "Mismatched array sizes. Argument array has "
                << arg.size() << " elements, while destination array has "
                << m_size << " elements.";
        BRICK_THROW(common::ValueException, "Array1D::operator/=(const Array1D&)",
                    message.str().c_str());
      }
      std::transform(m_dataPtr, m_dataPtr + m_size, arg.data(), m_dataPtr,
                     std::divides<Type>());
      return *this;
    }


    template <class Type>
    Array1D<Type>&
    Array1D<Type>::
    operator/=(const Type arg)
    {
      std::transform(m_dataPtr, m_dataPtr + m_size, m_dataPtr,
                     std::bind2nd(std::divides<Type>(), arg));
      return *this;
    }


    template <class Type>
    Array1D<Type>&
    Array1D<Type>::
    operator<<=(int numberOfBits)
    {
      for(int ii = this->size() - 1; ii >= 0; --ii) {
        this->m_dataPtr[ii] <<= numberOfBits;
      }      
      return *this;
    }

    
    template <class Type>
    Array1D<Type>&
    Array1D<Type>::
    operator>>=(int numberOfBits)
    {
      for(int ii = this->size() - 1; ii >= 0; --ii) {
        this->m_dataPtr[ii] >>= numberOfBits;
      }      
      return *this;
    }

    
    template <class Type>
    void Array1D<Type>::
    allocate(size_t arraySize)
    {
      // Make sure to release any shared resources.
      this->deAllocate();
      
      // Check array size.  It doesn't make sense to allocate memory
      // for a zero size array.
      m_size = arraySize;
      if(m_size > 0) {
        // Allocate data storage.  new() should throw an exception if
        // we run out of memory.
        m_dataPtr = new Type[m_size];
 
        // Set reference count to show that exactly one Array is pointing
        // to this data.
        m_referenceCount.reset(1);
      }
      return;
    }


    template <class Type>
    inline void Array1D<Type>::
    checkBounds(size_t
#ifdef BRICK_NUMERIC_CHECKBOUNDS
                index
#endif /* #ifdef BRICK_NUMERIC_CHECKBOUNDS */
      ) const
    {
#ifdef BRICK_NUMERIC_CHECKBOUNDS
      if(index >= m_size) {
        std::ostringstream message;
        message << "Index " << index << " is invalid for a(n) " << m_size
                << " element array.";
        BRICK_THROW(common::IndexException, "Array1D::checkBounds()",
                  message.str().c_str());
      }
#endif /* #ifdef BRICK_NUMERIC_CHECKBOUNDS */
    }


    template <class Type>
    void Array1D<Type>::
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
      m_referenceCount.reset(0);
    }

    /* Non-member functions which should maybe wind up in a different file. */
  
    template <class Type>
    Array1D<Type> operator+(const Array1D<Type>& array0,
                            const Array1D<Type>& array1)
    {
      if(array0.size() != array1.size()) {
        std::ostringstream message;
        message << "Array sizes do not match.  Array0 has " << array0.size()
                << " elements, while array1 has " << array1.size()
                << " elements.";
        BRICK_THROW(common::ValueException, "Array1D::operator+()", message.str().c_str());
      }
      Array1D<Type> result(array0.size());
      std::transform(array0.begin(), array0.end(), array1.begin(),
                     result.begin(), std::plus<Type>());
      return result;
    }


    template <class Type>
    Array1D<Type> operator-(const Array1D<Type>& array0,
                            const Array1D<Type>& array1)
    {
      if(array0.size() != array1.size()) {
        std::ostringstream message;
        message << "Array sizes do not match.  Array0 has " << array0.size()
                << " elements, while array1 has " << array1.size()
                << " elements.";
        BRICK_THROW(common::ValueException, "Array1D::operator-()", message.str().c_str());
      }
      Array1D<Type> result(array0.size());
      std::transform(array0.begin(), array0.end(), array1.begin(),
                     result.begin(), std::minus<Type>());
      return result;
    }


    template <class Type>
    Array1D<Type> operator*(const Array1D<Type>& array0,
                            const Array1D<Type>& array1)
    {
      if(array0.size() != array1.size()) {
        std::ostringstream message;
        message << "Array sizes do not match.  Array0 has " << array0.size()
                << " elements, while array1 has " << array1.size()
                << " elements.";
        BRICK_THROW(common::ValueException, "Array1D::operator*()", message.str().c_str());
      }
      Array1D<Type> result(array0.size());
      std::transform(array0.begin(), array0.end(), array1.begin(),
                     result.begin(), std::multiplies<Type>());
      return result;
    }


    template <class Type>
    Array1D<Type> operator/(const Array1D<Type>& array0,
                            const Array1D<Type>& array1)
    {
      if(array0.size() != array1.size()) {
        std::ostringstream message;
        message << "Array sizes do not match.  Array0 has " << array0.size()
                << " elements, while array1 has " << array1.size()
                << " elements.";
        BRICK_THROW(common::ValueException, "Array1D::operator/()", message.str().c_str());
      }
      Array1D<Type> result(array0.size());
      std::transform(array0.begin(), array0.end(), array1.begin(),
                     result.begin(), std::divides<Type>());
      return result;
    }


    template <class Type>
    Array1D<Type> operator+(const Array1D<Type>& array0, Type scalar)
    {
      Array1D<Type> result(array0.size());
      std::transform(array0.begin(), array0.end(), result.begin(),
                     std::bind2nd(std::plus<Type>(), scalar));
      return result;
    }


    template <class Type>
    Array1D<Type> operator-(const Array1D<Type>& array0, Type scalar)
    {
      Array1D<Type> result(array0.size());
      std::transform(array0.begin(), array0.end(), result.begin(),
                     std::bind2nd(std::minus<Type>(), scalar));
      return result;
    }


    template <class Type>
    Array1D<Type> operator*(const Array1D<Type>& array0, Type scalar)
    {
      Array1D<Type> result(array0.size());
      std::transform(array0.begin(), array0.end(), result.begin(),
                     std::bind2nd(std::multiplies<Type>(), scalar));
      return result;
    }


    template <class Type>
    Array1D<Type> operator/(const Array1D<Type>& array0, Type scalar)
    {
      Array1D<Type> result(array0.size());
      std::transform(array0.begin(), array0.end(), result.begin(),
                     std::bind2nd(std::divides<Type>(), scalar));
      return result;
    }


    template <class Type>
    inline Array1D<Type> operator+(Type scalar, const Array1D<Type>& array0)
    {
      return array0 + scalar;
    }


    template <class Type>
    Array1D<Type> operator-(Type scalar, const Array1D<Type>& array0)
    {
      Array1D<Type> result(array0.size());
      std::transform(array0.begin(), array0.end(), result.begin(),
                     std::bind1st(std::minus<Type>(), scalar));
      return result;
    }

  
    template <class Type>
    inline Array1D<Type> operator*(Type scalar, const Array1D<Type>& array0)
    {
      return array0 * scalar;
    }

  
    template <class Type>
    Array1D<Type> operator/(Type scalar, const Array1D<Type>& array0)
    {
      Array1D<Type> result(array0.size());
      std::transform(array0.begin(), array0.end(), result.begin(),
                     std::bind1st(std::divides<Type>(), scalar));
      return result;
    }
  

    // Elementwise comparison of an Array1D with a constant.
    template <class Type>
    Array1D<bool>
    operator==(const Array1D<Type>& array0, const Type arg)
    {
      Array1D<bool> result(array0.size());
      std::transform(array0.begin(), array0.end(), result.data(),
                     std::bind2nd(std::equal_to<Type>(), arg));
      return result;
    }

    
    // Elementwise comparison of an Array1D with another array.
    template <class Type>
    Array1D<bool>
    operator==(const Array1D<Type>& array0, const Array1D<Type>& array1)
    {
      array0.checkDimension(array1.size());
      Array1D<bool> result(array0.size());
      std::transform(array0.begin(), array0.end(), array1.begin(),
                     result.begin(), std::equal_to<Type>());
      return result;
    }

  
    template <class Type>
    Array1D<bool>
    operator>(const Array1D<Type>& array0, const Type arg)
    {
      Array1D<bool> result(array0.size());
      std::transform(array0.begin(), array0.end(), result.begin(),
                     std::bind2nd(std::greater<Type>(), arg));
      return result;
    }

  
    template <class Type>
    Array1D<bool>
    operator>=(const Array1D<Type>& array0, const Type arg)
    {
      Array1D<bool> result(array0.size());
      std::transform(array0.begin(), array0.end(), result.begin(),
                     std::bind2nd(std::greater_equal<Type>(), arg));
      return result;
    }


    template <class Type>
    Array1D<bool>
    operator<(const Array1D<Type>& array0, const Type arg)
    {
      Array1D<bool> result(array0.size());
      std::transform(array0.begin(), array0.end(), result.begin(),
                     std::bind2nd(std::less<Type>(), arg));
      return result;
    }


    template <class Type>
    Array1D<bool>
    operator<=(const Array1D<Type>& array0, const Type arg)
    {
      Array1D<bool> result(array0.size());
      std::transform(array0.begin(), array0.end(), result.begin(),
                     std::bind2nd(std::less_equal<Type>(), arg));
      return result;
    }


    template <class Type>
    std::ostream& operator<<(std::ostream& stream, const Array1D<Type>& array0)
    {
      // Most of the time, OutputType will be the same as Type.
      typedef typename NumericTraits<Type>::TextOutputType OutputType;
    
      if (!stream){
        BRICK_THROW(common::IOException,
                    "operator<<(std::ostream&, const Array1D&)",
                    "Invalid stream\n");
      }

      size_t index;
      stream << "Array1D([";
      if (array0.size() > 0){
        for(index = 0; index < array0.size() - 1; ++index) {
          stream << static_cast<OutputType>(array0(index)) << ", ";
        }
        stream << static_cast<OutputType>(array0(index));
      }
      stream << "])";
      return stream;
    }

    // Sets the value of an Array1D instance from a std::istream.
    template <class Type>
    std::istream& operator>>(std::istream& inputStream, Array1D<Type>& array0)
    {
      return array0.readFromStream(inputStream);
    }

  } // namespace numeric

} // namespace brick

#endif /* #ifndef BRICK_NUMERIC_ARRAY1D_IMPL_HH */
