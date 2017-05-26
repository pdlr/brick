/**
***************************************************************************
* @file brick/utilities/tee.hh
*
* Header file declaring Tee class.
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
    
    template <class CharType, class Traits = std::char_traits<CharType> >
    class TeeBuffer : public std::basic_streambuf<CharType>
    {
    public:
      typedef typename std::basic_streambuf<CharType>::int_type int_type;
      
      TeeBuffer() :
        std::basic_streambuf<CharType>(),
        m_sinkVector() {}

      void add(std::basic_streambuf<CharType>* sinkPtr) {
        m_sinkVector.push_back(sinkPtr);
      }

    private:
      // Disallow copying.
      TeeBuffer(TeeBuffer const&) = delete;
      // xxx move constructor.
      TeeBuffer& operator=(TeeBuffer const&) = delete;
      // xxx move operator.
      
      // Private member functions.

      /// Override std::basic_streambuf<CharType>::overflow().
      int_type overflow(int_type nextChar) {
        for(auto sinkPtr : m_sinkVector) {
          sinkPtr->sputc(nextChar);
        }
        return nextChar;
      }

      // Private data members.
      std::vector<std::basic_streambuf<CharType>*> m_sinkVector;
      
    }; // class TeeBuffer


    template <class CharType = char, class Traits = std::char_traits<CharType> >
    class Tee : public std::basic_ostream<CharType>
    {
    public:
      
      Tee() :
        std::basic_ostream<CharType>(),
        m_buffer()
        {
          std::ostream::rdbuf(&m_buffer);
        }

      void add(std::basic_ostream<CharType>& sink) {
        sink.flush();
        m_buffer.add(sink.rdbuf());
      }

    private:
      // Disallow copying.
      Tee(Tee const&) = delete;
      // xxx move constructor.
      Tee& operator=(Tee const&) = delete;
      // xxx move operator.
      
      // Private member functions.

      // Private data members.
      TeeBuffer<CharType> m_buffer;
      
    }; // class Tee
    
  } // namespace utilities
    
} // namespace brick

#endif /* #ifndef BRICK_UTILITIES_TEE_HH */
