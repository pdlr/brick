/**
***************************************************************************
* @file brick/numeric/extremumRecorder.hh
*
* Header file declaring ExtremumRecorder class Template.
*
* Copyright (C) 2008, 2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_EXTREMUMRECORDER_HH
#define BRICK_NUMERIC_EXTREMUMRECORDER_HH

#include <limits>

namespace brick {

  namespace numeric {
    
    /**
     ** A simple class template to help you avoid writing "if(myVar >
     ** maxVal) {maxVal = myVar; bestIndex = ii;} all over your code.
     ** Please see also the specializations MaxRecorder and
     ** MinRecorder.
     **/
    template<class Type, class Payload, class Functor>
    class ExtremumRecorder {
    public:

      /** 
       * This constructor explicitly initializes the stored extremum.
       * For example, if using operator>() as functor, you would want
       * to pass a very small number as the constructor argument, so
       * that the first call to member function test() would set a new
       * extremum (you might get this very small number using
       * std::numeric_limits<Type>::min(), or
       * -std::numeric_limits<Type>::max(), depending on the
       * signedness of Type).  This constructor also explicitly sets
       * the payload.
       */
      ExtremumRecorder(Type const& maxValue, Payload const& payload,
                       Functor functor)
        : m_extremeValue(maxValue), m_payload(payload), m_functor(functor) {}


      
      /** 
       * This member function is an alias for member function test().
       * 
       * @param value This argument will be compared to the saved
       * extreme value using the saved functor, and then copied into
       * the saved maximum if the functor returns true.
       * 
       * @param payload This argument is copied into the saved payload
       * if and only if the saved extreme value was updated.
       * 
       * @return The return value is true if the saved extreme value
       * was updated, false otherwise.
       */
      bool
      evaluate(const Type& value, const Payload& payload) {
        return this->test(value, payload);
      }


      /** 
       * This member function compares its first argument with the
       * saved extreme value, and updates the saved extremeValue (and
       * payload) if the specified functor returns true when called
       * with the arguments(value, previousExtremeValue).  For
       * example, if the functor is simply operator>(), this has the
       * effect of updating the saved value only if the new value is
       * larger.
       * 
       * @param value This argument will be compared to the saved
       * extreme value using the saved functor, and then copied into
       * the saved maximum if the functor returns true.
       * 
       * @param payload This argument is copied into the saved payload
       * if and only if the saved extreme value was updated.
       * 
       * @return The return value is true if the saved extreme value
       * was updated, false otherwise.
       */
      bool
      test(const Type& value, const Payload& payload) {
        if(m_functor(value, m_extremeValue)) {
          m_extremeValue = value;
          m_payload = payload;
          return true;
        }
        return false;
      }


      /** 
       * This member function returns the saved maximum.
       * 
       * @return The return value is the largest value passed to
       * member function test() since construction, or since the last
       * call to member function reset().
       */
      const Type&
      getExtremeValue() {return m_extremeValue;}


      /** 
       * This member function returns the saved "payload".
       * 
       * @return The return value is the payload associated with the
       * "most extreme" value passed to member function test() since
       * construction, or since the last call to member function
       * reset().
       */
      const Payload&
      getPayload() {return m_payload;}


      /** 
       * This member function resets the saved extreme value just as if
       * *this were freshly constructed.
       */
      void
      reset(Type const& value, Payload const& payload) {
        m_extremeValue = value;
        m_payload = payload;
      }

    private:

      Type m_extremeValue;
      Payload m_payload;
      Functor m_functor;

    };
    
  } // namespace numeric

} // namespace brick

#endif // #ifdef BRICK_NUMERIC_EXTREMUMRECORDER_HH
