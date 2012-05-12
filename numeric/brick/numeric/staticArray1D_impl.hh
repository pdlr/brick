/**
***************************************************************************
* @file brick/numeric/staticArray1D_impl.hh
*
* Header file defining inline and template functions declared in
* StaticArray1D.hh.
*
* Copyright (C) 2001-2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_STATICARRAY1D_IMPL_HH
#define BRICK_NUMERIC_STATICARRAY1D_IMPL_HH

// This file is included by staticArray1D.hh, and should not be
// directly included by user code, so no need to include
// staticArray1D.hh here.
// 
// #include <brick/numeric/staticArray1D.hh>

#include <algorithm>
#include <sstream>
#include <vector>
#include <brick/common/expect.hh>
#include <brick/numeric/numericTraits.hh>

namespace brick {

  namespace numeric {
    
    // Static constant describing how the string representation of an
    // StaticArray1D should start.
    template <class Type>
    const std::string&
    StaticArray1D<Type, Size>::
    ioIntro()
    {
      static const std::string intro = "StaticArray1D(";
      return intro;
    }

  
    // Static constant describing how the string representation of an
    // StaticArray1D should end.
    template <class Type>
    const std::string&
    StaticArray1D<Type, Size>::
    ioOutro()
    {
      static const std::string outro = ")";
      return outro;
    }

  
    // Non-static member functions below.

    template <class Type>
    StaticArray1D<Type, Size>::
    StaticArray1D()
      : m_dataArray()
    {
      // Empty.
    }

  
    // Construct from an initialization string.
    template <class Type>
    StaticArray1D<Type, Size>::
    StaticArray1D(const std::string& inputString)
      : m_dataArray()
    {
      // We'll use the stream input operator to parse the string.
      std::istringstream inputStream(inputString);

      // Now read the string into a staticArray.
      StaticArray1D<Type, Size> inputStaticArray;
      inputStream >> inputStaticArray;
      if(!inputStream) {
        std::ostringstream message;
        message << "Couldn't parse input string: \"" << inputString << "\".";
        BRICK_THROW3(ValueException,
                   "StaticArray1D::StaticArray1D(const std::string&)",
                   message.str().c_str());                 
      }

      // If all went well, copy into *this.
      *this = inputStaticArray;
    }

  
    // When copying from a StaticArray1D do a deep copy.
    template <class Type>
    StaticArray1D<Type, Size>::
    StaticArray1D(const StaticArray1D<Type, Size>& source)
      : m_dataArray()
    {
      std::copy(source.begin(), source.end(), m_dataArray);
    }

  
    template <class Type>
    StaticArray1D<Type, Size>::
    ~StaticArray1D()
    {
      // Empty.
    }

  
    template <class Type>
    StaticArray1D<Type, Size>& StaticArray1D<Type, Size>::
    operator=(Type val)
    {
      std::fill(this->begin(), this->end(), val);
      return *this;
    }


    template <class Type>
    StaticArray1D<Type, Size>& StaticArray1D<Type, Size>::
    operator=(const StaticArray1D<Type, Size>& source)
    {
      std::copy(source.begin(), source.end(), this->begin());
      return *this;
    }

  
    template <class Type> template <class Type2>
    StaticArray1D<Type, Size>&
    StaticArray1D<Type, Size>::
    operator+=(const StaticArray1D<Type2, Size>& arg)
    {
      std::transform(this->begin, this->end(), arg.begin(), this->begin(),
                     std::plus<Type>());
      return *this;
    }

  
    template <class Type>
    StaticArray1D<Type, Size>&
    StaticArray1D<Type, Size>::
    operator+=(const Type arg)
    {
      std::transform(this->begin(), this->end(), this->begin(),
                     std::bind2nd(std::plus<Type>(), arg));
      return *this;
    }

  
    template <class Type> template <class Type2>
    StaticArray1D<Type, Size>&
    StaticArray1D<Type, Size>::
    operator-=(const StaticArray1D<Type2, Size>& arg)
    {
      std::transform(this->begin(), this->end(), arg.begin(), this->begin(),
                     std::minus<Type>());
      return *this;
    }

  
    template <class Type>
    StaticArray1D<Type, Size>&
    StaticArray1D<Type, Size>::
    operator-=(const Type arg)
    {
      std::transform(this->begin(), this->end(), this->begin(),
                     std::bind2nd(std::minus<Type>(), arg));
      return *this;
    }

  
    template <class Type> template <class Type2>
    StaticArray1D<Type, Size>&
    StaticArray1D<Type, Size>::
    operator*=(const StaticArray1D<Type2, Size>& arg)
    {
      std::transform(this->begin(), this->end(), arg.begin(), this->begin(),
                     std::multiplies<Type>());
      return *this;
    }

  
    template <class Type>
    StaticArray1D<Type, Size>&
    StaticArray1D<Type, Size>::
    operator*=(const Type arg)
    {
      std::transform(this->begin(), this->end(), this->begin(),
                     std::bind2nd(std::multiplies<Type>(), arg));
      return *this;
    }

  
    template <class Type> template <class Type2>
    StaticArray1D<Type, Size>&
    StaticArray1D<Type, Size>::
    operator/=(const StaticArray1D<Type2, Size>& arg)
    {
      std::transform(this->begin(), this->end(), arg.begin(), this->begin(),
                     std::divides<Type>());
      return *this;
    }

  
    template <class Type>
    StaticArray1D<Type, Size>&
    StaticArray1D<Type, Size>::
    operator/=(const Type arg)
    {
      std::transform(this->begin(), this->end(), this->begin(),
                     std::bind2nd(std::divides<Type>(), arg));
      return *this;
    }

  
    template <class Type>
    StaticArray1D<bool, Size> StaticArray1D<Type, Size>::
    operator>(const Type arg)
      const
    {
      StaticArray1D<bool, Size> result(m_size);
      std::transform(this->begin(), this->end(), result.begin(),
                     std::bind2nd(std::greater<Type>(), arg));
      return result;
    }

  
    template <class Type>
    StaticArray1D<bool, Size> StaticArray1D<Type, Size>::
    operator>=(const Type arg)
      const
    {
      StaticArray1D<bool, Size> result(m_size);
      std::transform(this->begin(), this->end(), result.begin(),
                     std::bind2nd(std::greater_equal<Type>(), arg));
      return result;
    }

  
    template <class Type>
    StaticArray1D<bool, Size> StaticArray1D<Type, Size>::
    operator<(const Type arg)
      const
    {
      StaticArray1D<bool, Size> result(m_size);
      std::transform(this->begin(), this->end(), result.begin(),
                     std::bind2nd(std::less<Type>(), arg));
      return result;
    }

  
    template <class Type>
    StaticArray1D<bool, Size> StaticArray1D<Type, Size>::
    operator<=(const Type arg)
      const
    {
      StaticArray1D<bool, Size> result(m_size);
      std::transform(this->begin(), this->end(), result.begin(),
                     std::bind2nd(std::less_equal<Type>(), arg));
      return result;
    }

  
    template <class Type>
    inline void StaticArray1D<Type, Size>::
    checkBounds(size_t index) const
    {
#ifdef BRICK_NUMERIC_CHECKBOUNDS
      if(index < 0 || index >= Size) {
        std::ostringstream message;
        message << "Index " << index << " is invalid for a(n) " << m_size
                << " element staticArray.";
        BRICK_THROW(IndexException, "StaticArray1D::checkBounds()",
                  message.str().c_str());
      }
#endif
    }


    /* ========== Non-member functions =========== */
  
    template <class Type>
    StaticArray1D<Type, Size>
    operator+(const StaticArray1D<Type, Size>& staticArray0,
              const StaticArray1D<Type, Size>& staticArray1)
    {
      StaticArray1D<Type, Size> result(staticArray0.size());
      std::transform(staticArray0.begin(), staticArray0.end(),
                     staticArray1.begin(),
                     result.begin(), std::plus<Type>());
      return result;
    }

  
    template <class Type>
    StaticArray1D<Type, Size>
    operator-(const StaticArray1D<Type, Size>& staticArray0,
              const StaticArray1D<Type, Size>& staticArray1)
    {
      StaticArray1D<Type, Size> result(staticArray0.size());
      std::transform(staticArray0.begin(), staticArray0.end(),
                     staticArray1.begin(),
                     result.begin(), std::minus<Type>());
      return result;
    }

  
    template <class Type>
    StaticArray1D<Type, Size>
    operator*(const StaticArray1D<Type, Size>& staticArray0,
              const StaticArray1D<Type, Size>& staticArray1)
    {
      StaticArray1D<Type, Size> result(staticArray0.size());
      std::transform(staticArray0.begin(), staticArray0.end(),
                     staticArray1.begin(),
                     result.begin(), std::multiplies<Type>());
      return result;
    }

  
    template <class Type>
    StaticArray1D<Type, Size>
    operator/(const StaticArray1D<Type, Size>& staticArray0,
              const StaticArray1D<Type, Size>& staticArray1)
    {
      StaticArray1D<Type, Size> result(staticArray0.size());
      std::transform(staticArray0.begin(), staticArray0.end(),
                     staticArray1.begin(),
                     result.begin(), std::divides<Type>());
      return result;
    }

  
    template <class Type>
    StaticArray1D<Type, Size>
    operator+(const StaticArray1D<Type, Size>& staticArray0, Type scalar)
    {
      StaticArray1D<Type, Size> result(staticArray0.size());
      std::transform(staticArray0.begin(), staticArray0.end(), result.begin(),
                     std::bind2nd(std::plus<Type>(), scalar));
      return result;
    }

  
    template <class Type>
    StaticArray1D<Type, Size>
    operator-(const StaticArray1D<Type, Size>& staticArray0, Type scalar)
    {
      StaticArray1D<Type, Size> result(staticArray0.size());
      std::transform(staticArray0.begin(), staticArray0.end(), result.begin(),
                     std::bind2nd(std::minus<Type>(), scalar));
      return result;
    }


    template <class Type>
    StaticArray1D<Type, Size>
    operator*(const StaticArray1D<Type, Size>& staticArray0, Type scalar)
    {
      StaticArray1D<Type, Size> result(staticArray0.size());
      std::transform(staticArray0.begin(), staticArray0.end(), result.begin(),
                     std::bind2nd(std::multiplies<Type>(), scalar));
      return result;
    }


    template <class Type>
    StaticArray1D<Type, Size>
    operator/(const StaticArray1D<Type, Size>& staticArray0, Type scalar)
    {
      StaticArray1D<Type, Size> result(staticArray0.size());
      std::transform(staticArray0.begin(), staticArray0.end(), result.begin(),
                     std::bind2nd(std::divides<Type>(), scalar));
      return result;
    }


    template <class Type>
    inline StaticArray1D<Type, Size>
    operator+(Type scalar, const StaticArray1D<Type, Size>& staticArray0)
    {
      return staticArray0 + scalar;
    }


    template <class Type>
    StaticArray1D<Type, Size>
    operator-(Type scalar, const StaticArray1D<Type, Size>& staticArray0)
    {
      StaticArray1D<Type, Size> result(staticArray0.size());
      std::transform(staticArray0.begin(), staticArray0.end(), result.begin(),
                     std::bind1st(std::minus<Type>(), scalar));
      return result;
    }


    template <class Type>
    inline StaticArray1D<Type, Size>
    operator*(Type scalar, const StaticArray1D<Type, Size>& staticArray0)
    {
      return staticArray0 * scalar;
    }


    template <class Type>
    StaticArray1D<Type, Size>
    operator/(Type scalar, const StaticArray1D<Type, Size>& staticArray0)
    {
      StaticArray1D<Type, Size> result(staticArray0.size());
      std::transform(staticArray0.begin(), staticArray0.end(), result.begin(),
                     std::bind1st(std::divides<Type>(), scalar));
      return result;
    }
  

    template <class Type>
    std::ostream&
    operator<<(std::ostream& stream,
               const StaticArray1D<Type, Size>& staticArray0)
    {
      // Most of the time, OutputType will be the same as Type.
      typedef typename NumericTraits<Type>::TextOutputType OutputType;
    
      if (!stream){
        BRICK_THROW(IOException, "operator<<(std::ostream&, const StaticArray1D&)",
                  "Invalid stream\n");
      }

      size_t index;
      stream << "StaticArray1D([";
      if (staticArray0.size() > 0){
        for(index = 0; index < staticArray0.size() - 1; ++index) {
          stream << static_cast<OutputType>(staticArray0(index)) << ", ";
        }
        stream << static_cast<OutputType>(staticArray0(index));
      }
      stream << "])";
      return stream;
    }


    // Sets the value of an StaticArray1D instance from a std::istream.
    template <class Type>
    std::istream&
    operator>>(std::istream& inputStream,
               StaticArray1D<Type, Size>& staticArray0)
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
      
        // We won't require the input format to start with
        // "StaticArray1D(", but if it does we read it here.
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
        std::vector<Type> inputBuffer(Size);
        for(size_t index0 = 0; index0 < Size - 1; ++index0) {
          // Read the next value.
          inputStream >> inputValue;
          inputBuffer[index0] = static_cast<Type>(inputValue);

          // Read the separator.
          char inChar = 0;
          inputStream >> inChar;
          if(inChar != ioSeparator) {
            // Missing separator.  Fail here.
            inputStream.clear(std::ios_base::failbit);
          }
        }

        // Read the final value.
        inputStream >> inputValue;
        inputBuffer[Size - 1] = static_cast<Type>(inputValue);

        // Read the closing character.
        char inChar = 0;
        inputStream >> inChar;
        if(inChar != ioClosing) {
          // Missing Closing.  Fail here.
          inputStream.clear(std::ios_base::failbit);
        }

        // If we found an intro, we expect the corresponding outro.
        if(foundIntro) {
          inputStream >> common::Expect(&(ioOutro()), 1, flags);
        }

        // Now we're done with all of the parsing.  Copy the data to *this.
        std::copy(inputBuffer.begin(), inputBuffer.end(), this->begin());

      } catch(std::ios_base::failure) {
        // Empty
      }
      inputStream.exceptions(oldExceptionState);
      return inputStream;
    }

  } // namespace numeric

} // namespace brick

#endif /* #ifdef BRICK_NUMERIC_STATICARRAY1D_IMPL_HH */

