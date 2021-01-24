/**
***************************************************************************
* @file brick/numeric/extremumRecorder.hh
*
* Header file declaring MinRecorder class Template.
*
* Copyright (C) 2008, 2011 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_MINRECORDER_HH
#define BRICK_NUMERIC_MINRECORDER_HH

#include <functional>
#include <limits>
#include <brick/numeric/extremumRecorder.hh>

namespace brick {

  namespace numeric {

    /**
     ** A simple class template to help you avoid writing "if(myVar <
     ** minVal) {minVal = myVar; bestIndex = ii;} all over your code.
     **/
    template<class Type, class Payload>
    class MinRecorder
      : public ExtremumRecorder< Type, Payload, std::less<Type> >
    {
    public:

      /**
       * Default constructor initializes min to a very large number
       * (std::numeric_limits<Type>::max()), and sets the payload to
       * it's default value.
       */
      MinRecorder()
        : ExtremumRecorder< Type, Payload, std::less<Type> >(
            Type(), Payload(), std::less<Type>()) {
        this->reset();
      }


      /**
       * This constructor allows the user to explicitly set the starting
       * minValue, and the default payload.
       */
      MinRecorder(Type const& minValue, Payload const& payload)
        : ExtremumRecorder< Type, Payload, std::less<Type> >(
            minValue, payload, std::less<Type>()) {}


      // Member function test(), and its alias evaluate(), are
      // inherited from ExtremumRecorder.
      //
      // bool evaluate(const Type& value, const Payload& payload);
      // bool test(const Type& value, const Payload& payload);


      /**
       * This member function returns the saved minimum.
       *
       * @return The return value is the largest value passed to
       * member function test() since construction, or since the last
       * call to member function reset().
       */
      const Type&
      getMin() {return this->getMinimum();}


      /**
       * This member function returns the saved minimum.
       *
       * @return The return value is the largest value passed to
       * member function test() since construction, or since the last
       * call to member function reset().
       */
      const Type&
      getMinimum() {return this->getExtremeValue();}


      // Member function getPayload() is inherited from
      // ExtremumRecorder.
      //
      // const Payload& getPayload();


      /**
       * This member function resets the saved minimum just as if
       * *this were freshly constructed with the default constructor.
       */
      void
      reset() {
        this->ExtremumRecorder< Type, Payload, std::less<Type> >::reset(
          std::numeric_limits<Type>::max(), Payload());
      }
    };

  } // namespace numeric

} // namespace brick

#endif // #ifdef BRICK_NUMERIC_MAXRECORDER_HH
