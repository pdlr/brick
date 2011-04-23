/**
***************************************************************************
* @file brick/common/inputStream.cc
*
* Header file defining convenience functions and classes for working
* with input streams.
*
* Copyright (C) 2004-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <iostream>
#include <brick/common/inputStream.hh>

namespace brick {

  namespace common {

    const Expect::FormatFlag Expect::NoFlag(0x0);
    const Expect::FormatFlag Expect::SkipWhitespace(0x1);
    const Expect::FormatFlag Expect::Sloppy(0x2);
    
    
    // Constructor.
    Expect::
    Expect(std::string const& expectation, FormatFlag formatFlag)
      : m_formatFlag(formatFlag),
        m_target(expectation.c_str()),
        m_targetLength(expectation.size())
    {
      // Empty.
    }

      
    // Constructor.
    Expect::
    Expect(char const* expectationPtr, FormatFlag formatFlag)
      : m_formatFlag(formatFlag),
        m_target(expectationPtr),
        m_targetLength(strlen(expectationPtr))
    {
      // Empty.
    }

    
    // Constructor.
    Expect::
    Expect(char const* expectationPtr, size_t length, FormatFlag formatFlag)
      : m_formatFlag(formatFlag),
        m_target(expectationPtr),
        m_targetLength(length)
    {
      // Empty.
    }

    
    // Destructor.
    Expect::
    ~Expect()
    {
      // Empty.
    }
      

    std::istream&
    Expect::
    expect(std::istream& stream) const
    {
      // We need an area into which to read.
      std::istream::char_type inputBuffer[s_chunkSize + 1];

      // Remove leading white space, if required.
      if((m_formatFlag & SkipWhitespace) != 0) {
        this->skipWhiteSpace(stream);
      }
      
      // inputString might be very long, so we keep count of how much
      // we've read.
      size_t readCount = 0;
      while(readCount < m_targetLength) {
        // Figure out how much to read.
        size_t thisChunkCount = m_targetLength - readCount;

        // But be sure not to overflow the buffer.
        if(thisChunkCount > s_chunkSize) {
          thisChunkCount = s_chunkSize;
        }

        // Read the appropriate number of characters.
        stream.read(inputBuffer, static_cast<std::streamsize>(thisChunkCount));
          
        // Quit if we weren't able to read enough bytes.
        if(static_cast<size_t>(stream.gcount()) != thisChunkCount) {
          stream.clear(std::ios_base::failbit);
          break;
        }
          
        // Now evaluate what we've read.
        bool matchFlag = true;
        if((m_formatFlag & Sloppy) == 0) {
          // We use memcmp() instead of strncmp() because the caller
          // might conceivably expect a string with a '\0' in the middle.
          matchFlag = (std::memcmp(reinterpret_cast<const void*>(
                                     m_target + readCount),
                                   reinterpret_cast<const void*>(inputBuffer),
                                   thisChunkCount * sizeof(char))
                       == 0);
        }

        // Quit if it wasn't what we expected.
        if(!matchFlag) {
          stream.clear(std::ios_base::failbit);
          break;
        }
        
        // Increment our count of how many characters read.
        readCount += thisChunkCount;
      }
      return stream;
    }

    
    // Skips to the next non-whitespace character.
    std::istream&
    Expect::
    skipWhiteSpace(std::istream& stream) const {
      std::istream::char_type inputChar;
      stream >> inputChar;
      stream.putback(inputChar);
      return stream;
    }
    
    
    std::istream&
    operator>>(std::istream& stream, Expect const& expect)
    {
      return expect.expect(stream);
    }

  } // namespace common
  
} // namespace brick
