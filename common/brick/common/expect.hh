/**
***************************************************************************
* @file brick/common/inputStream.hh
*
* Header file declaring convenience functions and classes for working
* with input streams.
*
* Copyright (C) 2004-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_COMMON_INPUTSTREAM_HH
#define BRICK_COMMON_INPUTSTREAM_HH

#include <cstring>
#include <iostream>

namespace brick {

  namespace common {

    template <class Type = unsigned int>
    class Enum {
    public:
      Enum(Type value = 0) : m_value(value) {}
      ~Enum() {}

      operator bool() const {return static_cast<bool>(m_value);}
      
      Enum operator&(Enum<Type> const& other) const {
        return Enum<Type>(m_value & other.m_value);
      }
      Enum operator|(Enum<Type> const& other) const {
        return Enum<Type>(m_value | other.m_value);
      }
      

      Type getValue();

    private:
      Type m_value;
    };
    
    
    /**
     ** This class reads from an inputstream and checks to make sure
     ** the input matches expectations.  Use it like this:
     **
     ** @code
     **   // Strict: only pass if the exact string is read.
     **   myInputStream >> Expect("Number attending:") >> attendeeCount;
     **
     **   // Less strict: accept leading whitespace without issues.
     **   myInputStream >> Expect("Number attending:", Expect::SkipWhitespace)
     **                 >> attendeeCount;
     **
     **   // Don't actually look at the input.  Just read an
     **   // appropriate number of characters.
     **   myInputStream >> Expect("Number attending:", Expect::Sloppy)
     **                 >> attendeeCount;
     **
     **   // First skip any leading whitespace, then read the required
     **   // number of input characters without actually looking at
     **   // them.
     **   myInputStream >> Expect("Number attending:",
     **                           Expect::Sloppy | Expect::SkipWhitespace)
     **                 >> attendeeCount;
     **   
     ** @endCode
     **
     ** If the input doesn't match expectation, failbit will be set on
     ** the input stream.
     **/
    class Expect {
    public:

# if 0
      /**
       ** This local type represents available modifiers to Expect's
       ** behavior.  See the usage example in the class documentation
       ** for more information.
       **/
      typedef enum {
        /// Default behavior.  Expect target exactly.
        NoFlag         = 0x00,

        /// Remove any preceding whitespace from the stream input
        /// before comparing against the expectation.
        SkipWhitespace = 0x01,

        /// Simply read the correct number of bytes from the input
        /// stream, but don't actually verify the that the bytes have the
        /// right value.
        Sloppy         = 0x02,
      } FormatEnum;

      /**
       ** This type gives us a way to represent the result of applying
       ** bitwise operations to FormatEnum values.
       **/
      struct FormatFlag {
        unsigned int skipWhitespace : 1;
        unsigned int sloppy : 1;

        FormatFlag(unsigned int formatEnum)
          : skipWhitespace(formatEnum & SkipWhitespace),
            sloppy(formatEnum & Sloppy)
          {}
      };

#endif

      struct FormatFlag : public Enum<unsigned int> {
        FormatFlag(unsigned int value) : Enum<unsigned int>(value){}
        FormatFlag(Enum<unsigned int> const& value)
          : Enum<unsigned int>(value){}
      };
      static const FormatFlag NoFlag;
      static const FormatFlag SkipWhitespace;
      static const FormatFlag Sloppy;
      
      
      /**
       * Constructor.
       */
      Expect(std::string const& expectation, FormatFlag formatFlag = NoFlag);

      
      /**
       * Constructor.
       */
      Expect(char const* expectationPtr, FormatFlag formatFlag = NoFlag);


      /**
       * Constructor.
       */
      Expect(char const* expectationPtr, size_t length,
             FormatFlag formatFlag = NoFlag);
      

      /**
       * Destructor.
       */
      virtual
      ~Expect();


      std::istream&
      expect(std::istream& stream) const;
      
    private:

      // Skips to the next non-whitespace character.
      std::istream&
      skipWhiteSpace(std::istream& stream) const;
      

      static const size_t s_chunkSize = 1024;

      FormatFlag m_formatFlag;
      char const* m_target;
      unsigned int m_targetLength;
    };


    std::istream&
    operator>>(std::istream& stream, Expect const& expect);

  } // namespace common
  
} // namespace brick


#endif /* #ifndef BRICK_COMMON_INPUTSTREAM_HH */
