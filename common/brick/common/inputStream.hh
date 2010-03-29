/**
***************************************************************************
* @file brick/common/inputStream.hh
*
* Header file declaring InputStream class.
*
* Copyright (C) 2004-2010 David LaRose, dlr@cs.cmu.edu
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
     ** The InputStream class is used to add convenience functions to
     ** existing istreams.  For example, it provides the expect() member
     ** function, which helps with reading strictly formatted input.
     ** Ideally, this would inherit most of its member functions from
     ** StreamType, but alas that requires access to stuff that isn't
     ** available in the the public std::istream interface.
     **/
    template <class StreamType=std::istream>
    class InputStream
    {
    public:
      /* ============ Typedefs "inherited" from StreamType ============ */
      /// This typedef simply mirrors StreamType::char_type.
      typedef typename StreamType::char_type char_type;


      /* ================== New typedefs and constants ================== */
    
      /**
       ** These bit flag constants are used to control InputStream
       ** behavior.  You pass them as constructor arguments, or else as
       ** arguments to the (not yet implemented) setFlags member
       ** function.
       **/
      typedef unsigned int FormatFlag;

    
      /**
       * NO_FLAG is the default state.
       */
      static const FormatFlag NO_FLAG         = 0x00;

    
      /**
       * SKIP_WHITESPACE indicates that expect() member function calls
       * should remove any preceding whitespace from the stream input
       * before comparing against the expected input.
       */
      static const FormatFlag SKIP_WHITESPACE = 0x01;

    
      /**
       * SLOPPY_EXPECT indicates that expect() member function calls
       * should simply read the correct number of bytes from the input
       * stream, but not actually verify the that the bytes have the
       * right value.
       */
      static const FormatFlag SLOPPY_EXPECT   = 0x02;

    
      /** 
       * This constructor pretends it's copying an input stream.
       * Unfortunately, the istream copy constructor is private, so
       * instead we just take a reference.  You might use this
       * constructor like this:
       *
       *   InputStream myStream(std::cin,
       *                        (InputStream::SKIP_WHITESPACE
       *                         | InputStream::SLOPPY_EXPECT));
       *
       * @param inputStream This is the stream to be copied.
       * 
       * @param flags This argument specifies flags that influence the
       * behavior of the InputStream instance.
       */
      InputStream(StreamType& inputStream,
                  FormatFlag flags=NO_FLAG)
        : m_istream(inputStream),
          m_skipWhiteSpace((flags & SKIP_WHITESPACE) != 0),
          m_sloppyExpect((flags & SLOPPY_EXPECT) != 0) {}


      /** 
       * Destructor.
       */
      ~InputStream() {}

      
      /* =================== Conversion operators =================== */

      /** 
       * This operator converts to void* by dispatching to the
       * StreamType conversion operator.
       * 
       * @return The return value is a void*, or 0 if this->fail().
       */
      operator void*() const {
        return m_istream;
      }


      /** 
       * This operator converts to StreamType reference.
       * 
       * @return The return value is a reference to the StreamType
       * we're pretending we're derived from.
       */
      operator StreamType&() const {
        return m_istream;
      }


      /* ============ Operators "inherited" from StreamType ============ */

      /** 
       * This operator simply dispatches to the corresponding
       * StreamType operator.
       * 
       * @param typeInstance This argument is passed directly to the
       * corresponding istream operator.
       *
       * @return The return value is a reference to *this.
       */
      template <class Type>
      InputStream&
      operator>>(Type& typeInstance) {
        m_istream >> typeInstance;
        return *this;
      }


      /* ======== Member functions "inherited" from StreamType ======== */

      /** 
       * This member function simply dispatches to the corresponding
       * StreamType member function.
       * 
       * @return The return value is simply the return value of the
       * corresponding StreamType member function.
       */
      bool
      bad() {
        return m_istream.bad();
      }

    
      /** 
       * This member function simply dispatches to the corresponding
       * StreamType member function.
       * 
       * @param state This argument is passed directly to the
       * corresponding istream member function.
       */
      void
      clear(typename StreamType::iostate state=StreamType::goodbit) {
        m_istream.clear(state);
      }


      /** 
       * This member function simply dispatches to the corresponding
       * StreamType member function.
       * 
       * @return The return value is simply the return value of the
       * corresponding StreamType member function.
       */
      bool
      eof() {
        return m_istream.eof();
      }

    
      /** 
       * This member function simply dispatches to the corresponding
       * StreamType member function.
       * 
       * @return The return value is simply the return value of the
       * corresponding StreamType member function.
       */
      bool
      fail() {
        return m_istream.fail();
      }

    
      /** 
       * This member function simply dispatches to the corresponding
       * StreamType member function.
       * 
       * @return The return value is simply the return value of the
       * corresponding StreamType member function.
       */
      typename StreamType::int_type
      get() {
        return m_istream.get();
      }


      /** 
       * This member function reads no more than the (numberOfCharacters
       * - 1) characters into an inputBuffer, stopping at the first
       * newline character, and terminating with '\0'.
       * 
       * @param inputBuffer This argument specifies buffer into which
       * characters will be read.
       * 
       * @param numberOfCharacters This argument specifies the size of
       * the buffer.
       * 
       * @return 
       */
      InputStream
      getline(char_type* inputBuffer, std::streamsize numberOfCharacters) {
        return m_istream.getline(inputBuffer, numberOfCharacters);
      }
    
    
      /** 
       * This member function simply dispatches to the corresponding
       * StreamType member function.
       * 
       * @return The return value is simply the return value of the
       * corresponding StreamType member function.
       */
      typename StreamType::int_type
      peek() {
        return m_istream.peek();
      }

    
      /* =================== New member functions =================== */

      /** 
       * This member function reads one character from the input and sets
       * failbit if doesn't match the specified value.
       * 
       * @param inputChar This argument defines what the upcoming stream
       * input character should be.
       *
       * @return The return value is a reference to *this.
       */
      InputStream&
      expect(typename StreamType::char_type inputChar) {
        typename StreamType::char_type readChar = 0;
        *this >> readChar;
        if(readChar != inputChar) {
          this->clear(std::ios_base::failbit);
        }
        return *this;
      }


      /** 
       * This member function reads characters from the input and sets
       * failbit if they don't match the specified string.
       * 
       * @param inputString This argument defines what the upcoming
       * stream input characters should be.
       *
       * @return The return value is a reference to *this.
       */
      InputStream&
      expect(const std::string& inputString) {
        return this->expect(inputString.c_str(), inputString.size());
      }


      /** 
       * This member function reads characters from the input and sets
       * failbit if they don't match the specified string.  The second
       * argument is optional.  If it is not provided, the number of
       * characters to read and compare will be determined by calling
       * strlen() on the first argument.  If the second argument is
       * provided, then it specifies the number of characters to read
       * and compare.  Note that this second mode does not respect '\0'
       * termination of the input string, so this call:
       *
       *   inputStream.expect("hello world", 240); // Catastrophe!
       *
       * will likely cause problems by reading past the end of the input
       * string.
       * 
       * @param inputCString This argument defines what the upcoming
       * stream input characters should be.
       *
       * @param stringSize This argument, if specified, indicates how
       * many characters to read and compare.
       * 
       * @return The return value is a reference to *this.
       */
      InputStream&
      expect(const char* inputCString, size_t stringSize=0) {
        // We need an area into which to read.
        typename StreamType::char_type inputBuffer[s_chunkSize + 1];

        // Remove leading white space, if required.
        if(m_skipWhiteSpace) {
          this->skipWhiteSpace();
        }
      
        // inputString might be very long, so we keep count of how much
        // we've read.
        size_t readCount = 0;

        // If the caller hasn't specified how long the string is, we
        // find out using strlen().
        if(stringSize == 0) {
          stringSize = strlen(inputCString);        
        }
      
        while(readCount < stringSize) {
          // Figure out how much to read.
          size_t thisChunkCount = stringSize - readCount;

          // But be sure not to overflow the buffer.
          if(thisChunkCount > s_chunkSize) {
            thisChunkCount = s_chunkSize;
          }

          // Read the appropriate number of characters.
          m_istream.read(inputBuffer,
                         static_cast<std::streamsize>(thisChunkCount));
          
          // Quit if we weren't able to read enough bytes.
          if(static_cast<size_t>(m_istream.gcount()) != thisChunkCount) {
            this->clear(std::ios_base::failbit);
            break;
          }
          
          // Now evaluate what we've read.
          bool matchFlag = true;
          if(!m_sloppyExpect) {
            // We use memcmp() instead of strncmp() because the caller
            // might conceivably expect a string with a '\0' in the middle.
            matchFlag = (
              std::memcmp(reinterpret_cast<const void*>(inputCString
                                                        + readCount),
                          reinterpret_cast<const void*>(inputBuffer),
                          thisChunkCount * sizeof(char)) == 0);
          }

          // Quit if it wasn't what we expected.
          if(!matchFlag) {
            this->clear(std::ios_base::failbit);
            break;
          }
        
          // Increment our count of how many characters read.
          readCount += thisChunkCount;
        }
        return *this;
      }
    

      /** 
       * This member function skips characters up to the next non-whitespace.
       */
      void
      skipWhiteSpace() {
        typename StreamType::char_type inputChar;
        m_istream >> inputChar;
        m_istream.putback(inputChar);
      }

    
    private:
      static const size_t s_chunkSize = 1024;
      StreamType& m_istream;
      bool m_skipWhiteSpace;
      bool m_sloppyExpect;
    };

  } // namespace common
  
} // namespace brick


#endif /* #ifndef BRICK_COMMON_INPUTSTREAM_HH */
