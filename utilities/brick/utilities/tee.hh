/**
***************************************************************************
* @file brick/utilities/tee.hh
*
* Header file declaring a class for sending out put to multiple streams
* at once.
*
* Copyright (C) 2017 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_UTILITIES_TEE_HH
#define BRICK_UTILITIES_TEE_HH

#include <ostream>
#include <streambuf>
#include <vector>

namespace brick {

  namespace utilities {

    /**
     ** This output stream buffer class template is used internally by
     ** class template Tee<CharType, Traits>.  It maintains a vector
     ** of client streambuffers, which it calls "sinks," and directs
     ** all of its output to each sink.
     **/
    template <class CharType, class Traits = std::char_traits<CharType> >
    class TeeBuffer : public std::basic_streambuf<CharType>
    {
    public:
      typedef typename std::basic_streambuf<CharType>::int_type int_type;

      /// After default construction, output goes nowhere.
      TeeBuffer() :
        std::basic_streambuf<CharType>(),
        m_sinkVector() {}

      /// Add a "sink" streambuffer.  Output will get copied to *sinkPtr.
      void add(std::basic_streambuf<CharType>* sinkPtr) {
        m_sinkVector.push_back(sinkPtr);
      }

    private:
      // Disallow copying.
      TeeBuffer(TeeBuffer const&) = delete;
      TeeBuffer(TeeBuffer const&&) = delete;
      TeeBuffer& operator=(TeeBuffer const&) = delete;
      TeeBuffer& operator=(TeeBuffer const&&) = delete;
      
      // Private member functions.

      /// Override std::basic_streambuf<CharType>::overflow().
      int_type overflow(int_type nextChar) {
        int_type returnValue = 0;
        for(auto sinkPtr : m_sinkVector) {
          int_type putcResult = sinkPtr->sputc(nextChar);
          if(Traits::eq_int_type(putcResult, Traits::eof())) {
            returnValue = Traits::eof();
          }
        }
        return returnValue;
      }

      // Private data members.
      std::vector<std::basic_streambuf<CharType>*> m_sinkVector;
      
    }; // class TeeBuffer


    /**
     ** This class allows writing the same thing to multiple output
     ** streams.  It is named after the unix "tee" command, which
     ** sends its standard input to both standard output and a log
     ** file.  Use it like this:
     **
     ** @code
     **   std::ofstream outputFileStream("outputFile");
     **   Tee<char> logStream;
     **   logStream.add(outputFileStream);
     **   logStream.add(std::cout);
     **   logStream << "This is a log message" << std::endl;
     * @endcode
     *
     * The code snippet above will print "This is a log message" on
     * std::cout, and also write it to the file "outputFile."
     **/
    template <class CharType = char, class Traits = std::char_traits<CharType> >
    class Tee : public std::basic_ostream<CharType>
    {
    public:

      /** 
       * After default construction, output goes nowhere.  Output directed
       * to a Tee instance in it's default state is simply ignored.
       */
      Tee() :
        std::basic_ostream<CharType>(),
        m_buffer()
        {
          std::ostream::rdbuf(&m_buffer);
        }

      /** 
       * Arrange for output sent to *this to be implicitly sent to
       * another stream.  If you call add() with multiple ostreams,
       * then any stream output sent to *this will be directed to each
       * of those ostreams.
       * 
       * @param sink This argument is the ostream that should receive
       * output from *this.
       */
      void add(std::basic_ostream<CharType>& sink) {
        sink.flush();
        m_buffer.add(sink.rdbuf());
      }

    private:
      // Disallow copying.
      Tee(Tee const&) = delete;
      Tee(Tee const&&) = delete;
      Tee& operator=(Tee const&) = delete;
      Tee& operator=(Tee const&&) = delete;
      
      // Private member functions.

      // Private data members.
      TeeBuffer<CharType, Traits> m_buffer;
      
    }; // class Tee
    
  } // namespace utilities
    
} // namespace brick

#endif /* #ifndef BRICK_UTILITIES_TEE_HH */
