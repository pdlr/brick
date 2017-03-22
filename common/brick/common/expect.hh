/**
***************************************************************************
* @file brick/common/expect.hh
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
    
    /**
     ** This class reads from an inputstream and checks to make sure
     ** the input matches expectations.  Use it like this:
     **
     ** @code
     **   // Strict: only pass if the exact string is read.
     **   myInputStream >> Expect("Number attending:") >> attendeeCount;
     **
     **   // Less strict: accept leading whitespace without issues.
     **   Expect::FormatFlag flags = Expect::SkipWhitespace
     **   myInputStream >> Expect("Number attending:", flags) >> attendeeCount;
     **
     **   // Don't actually look at the input.  Just read an
     **   // appropriate number of characters.
     **   myInputStream >> Expect("Number attending:", Expect::Sloppy)
     **                 >> attendeeCount;
     **
     **   // First skip any leading whitespace, then read the required
     **   // number of input characters without actually looking at
     **   // them.
     **   Expect::FormatFlag flags = Expect::SkipWhitespace | Expect::Sloppy;
     **   myInputStream >> Expect("Number attending:", flags) >> attendeeCount;
     **   
     ** @endCode
     **
     ** If the input doesn't match expectation, failbit will be set on
     ** the input stream.
     **/
    class Expect {
    public:

      /**
       ** This local type represents available modifiers to Expect's
       ** behavior.  See the usage example in the class documentation
       ** for more information.  We'd very much like to use an enum
       ** here, but we want to allow bitwise operations on these
       ** modifiers, and that gets tricky with enums (the result of
       ** the bitwise operation is an integral type, which you can't
       ** cast back to an enum).  We don't want to use an integral
       ** type directly, because that would introduce an ambiguity in
       ** the Expect constructors (Expect(char const*, size_t) vs
       ** Expect(char const*, FormatFlag)).  Instead, we create this
       ** surrogate.
       **/
      struct FormatFlag {
        unsigned int value;
        explicit FormatFlag(unsigned int arg = 0) : value(arg) {}
        operator unsigned int() const {return value;}
      };

      // Ugh.  We go through a contortion below (using static functions
      // instead of static members) to avoid the static initialization
      // order fiasco.
      
      /// Default behavior.  Expect target exactly.
      static FormatFlag NoFlag() {
        static FormatFlag returnValue(0x0);
        return returnValue;
      }        

      /// Remove any preceding whitespace from the stream input
      /// before comparing against the expectation.
      static FormatFlag SkipWhitespace() {
        static FormatFlag returnValue(0x1);
        return returnValue;
      }

      /// Simply read the correct number of bytes from the input
      /// stream, but don't actually verify the that the bytes have
      /// the right value.  This is useful if you need to read from
      /// disk as fast as possible.
      static FormatFlag Sloppy() {
        static FormatFlag returnValue(0x2);
        return returnValue;
      }
      
      
      /**
       * Constructor.
       */
      Expect(std::string const& expectation,
             FormatFlag formatFlag = Expect::NoFlag());

      
      /**
       * Constructor.
       */
      Expect(char const* expectationPtr,
             FormatFlag formatFlag = Expect::NoFlag());


      /**
       * Constructor.
       */
      Expect(char const* expectationPtr, size_t length,
             FormatFlag formatFlag = Expect::NoFlag());
      

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
