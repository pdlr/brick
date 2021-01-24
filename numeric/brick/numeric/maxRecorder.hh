/**
***************************************************************************
* @file brick/numeric/extremumRecorder.hh
*
* Header file declaring MaxRecorder class Template.
*
* Copyright (C) 2008, 2011 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_MAXRECORDER_HH
#define BRICK_NUMERIC_MAXRECORDER_HH

#include <functional>
#include <limits>
#include <brick/numeric/extremumRecorder.hh>

namespace brick {

  namespace numeric {

    /**
     ** A simple class template to help you avoid writing "if(myVar >
     ** maxVal) {maxVal = myVar; bestIndex = ii;} all over your code.
     **/
    template<class Type, class Payload>
    class MaxRecorder
      : public ExtremumRecorder< Type, Payload, std::greater<Type> >
    {
    public:

      /**
       * Default constructor initializes max to a very small number
       * (either std::numeric_limits<Type>::min(),
       * -std::numeric_limits<Type>::max(), depending on the
       * signedness of Type), and sets the payload to it's default
       * value.
       */
      MaxRecorder()
        : ExtremumRecorder< Type, Payload, std::greater<Type> >(
            Type(), Payload(), std::greater<Type>()) {
        this->reset();
      }


      /**
       * This constructor allows the user to explicitly set the starting
       * maxValue, and the default payload.
       */
      MaxRecorder(Type const& maxValue, Payload const& payload)
        : ExtremumRecorder< Type, Payload, std::greater<Type> >(
            maxValue, payload, std::greater<Type>()) {}


      // Member function test(), and its alias evaluate(), are
      // inherited from ExtremumRecorder.
      //
      // bool evaluate(const Type& value, const Payload& payload);
      // bool test(const Type& value, const Payload& payload);


      /**
       * This member function returns the saved maximum.
       *
       * @return The return value is the largest value passed to
       * member function test() since construction, or since the last
       * call to member function reset().
       */
      const Type&
      getMax() {return this->getMaximum();}


      /**
       * This member function returns the saved maximum.
       *
       * @return The return value is the largest value passed to
       * member function test() since construction, or since the last
       * call to member function reset().
       */
      const Type&
      getMaximum() {return this->getExtremeValue();}


      // Member function getPayload() is inherited from
      // ExtremumRecorder.
      //
      // const Payload& getPayload();


      /**
       * This member function resets the saved maximum just as if
       * *this were freshly constructed with the default constructor.
       */
      void
      reset() {
        if(std::numeric_limits<Type>::is_signed) {
          this->ExtremumRecorder< Type, Payload, std::greater<Type> >::reset(
            -std::numeric_limits<Type>::max(), Payload());
        } else {
          this->ExtremumRecorder< Type, Payload, std::greater<Type> >::reset(
            std::numeric_limits<Type>::min(), Payload());
        }
      }
    };

  } // namespace numeric

} // namespace brick

#endif // #ifdef BRICK_NUMERIC_MAXRECORDER_HH
