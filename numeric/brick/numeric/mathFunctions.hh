/**
***************************************************************************
* @file brick/numeric/mathFunctions.hh
*
* Header file declaring general math function templates to replace
* abs(), sin(), cos(), etc in portable code.  This builds on
* brick::common::mathFunctions.hh, but uses explicit function
* overloading rather than genereric templates.  The advantage of
* explicit overloading is that the generic templated versions of math
* functions are often inappropriate for unanticipated types.  If you
* use brick::numeric functions, you can be confident that the generic
* math functions won't be arbitrarily applied to types for which they
* aren't appropriate.
*
* Copyright (C) 2015, David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_MATHFUNCTIONS_HH
#define BRICK_NUMERIC_MATHFUNCTIONS_HH

#include <brick/common/mathFunctions.hh>

namespace brick {

  namespace numeric {


/**
 ** This macro makes it easy to dispatch from a brick::numeric
 ** overloaded function to the function template of the same name in
 ** brick::common.  It is used below to generate families of math
 ** functions.
 **/
#define BRICK_NUMERIC_DECLARE_MATHFUNCTION_1(FunctionName, Type)               \
    inline Type FunctionName(Type arg) {                                       \
      return brick::common::FunctionName<Type>(arg);                           \
    }

#define BRICK_NUMERIC_DECLARE_MATHFUNCTION_2(FunctionName, Type)               \
    inline Type FunctionName(Type arg0, Type arg1) {                           \
      return brick::common::FunctionName<Type>(arg0, arg1);                    \
    }
    
#define BRICK_NUMERIC_DECLARE_MATHFUNCTION_FAMILY_1(FunctionName)              \
    BRICK_NUMERIC_DECLARE_MATHFUNCTION_1(FunctionName, brick::common::Int8)    \
    BRICK_NUMERIC_DECLARE_MATHFUNCTION_1(FunctionName, brick::common::Int16)   \
    BRICK_NUMERIC_DECLARE_MATHFUNCTION_1(FunctionName, brick::common::Int32)   \
    BRICK_NUMERIC_DECLARE_MATHFUNCTION_1(FunctionName, brick::common::Int64)   \
    BRICK_NUMERIC_DECLARE_MATHFUNCTION_1(FunctionName, brick::common::UInt8)   \
    BRICK_NUMERIC_DECLARE_MATHFUNCTION_1(FunctionName, brick::common::UInt16)  \
    BRICK_NUMERIC_DECLARE_MATHFUNCTION_1(FunctionName, brick::common::UInt32)  \
    BRICK_NUMERIC_DECLARE_MATHFUNCTION_1(FunctionName, brick::common::UInt64)  \
    BRICK_NUMERIC_DECLARE_MATHFUNCTION_1(FunctionName, brick::common::Float32) \
    BRICK_NUMERIC_DECLARE_MATHFUNCTION_1(FunctionName, brick::common::Float64)

#define BRICK_NUMERIC_DECLARE_MATHFUNCTION_FAMILY_2(FunctionName)              \
    BRICK_NUMERIC_DECLARE_MATHFUNCTION_2(FunctionName, brick::common::Int8)    \
    BRICK_NUMERIC_DECLARE_MATHFUNCTION_2(FunctionName, brick::common::Int16)   \
    BRICK_NUMERIC_DECLARE_MATHFUNCTION_2(FunctionName, brick::common::Int32)   \
    BRICK_NUMERIC_DECLARE_MATHFUNCTION_2(FunctionName, brick::common::Int64)   \
    BRICK_NUMERIC_DECLARE_MATHFUNCTION_2(FunctionName, brick::common::UInt8)   \
    BRICK_NUMERIC_DECLARE_MATHFUNCTION_2(FunctionName, brick::common::UInt16)  \
    BRICK_NUMERIC_DECLARE_MATHFUNCTION_2(FunctionName, brick::common::UInt32)  \
    BRICK_NUMERIC_DECLARE_MATHFUNCTION_2(FunctionName, brick::common::UInt64)  \
    BRICK_NUMERIC_DECLARE_MATHFUNCTION_2(FunctionName, brick::common::Float32) \
    BRICK_NUMERIC_DECLARE_MATHFUNCTION_2(FunctionName, brick::common::Float64)
    
    // Please see brick/common/mathFunctions.hh for documentation of
    // math functions.

    // Declarations of single-argument overloads.
    BRICK_NUMERIC_DECLARE_MATHFUNCTION_FAMILY_1(absoluteValue);
    BRICK_NUMERIC_DECLARE_MATHFUNCTION_FAMILY_1(arccosine);
    BRICK_NUMERIC_DECLARE_MATHFUNCTION_FAMILY_1(arcsine);
    BRICK_NUMERIC_DECLARE_MATHFUNCTION_FAMILY_1(cosine);
    BRICK_NUMERIC_DECLARE_MATHFUNCTION_FAMILY_1(sine);
    BRICK_NUMERIC_DECLARE_MATHFUNCTION_FAMILY_1(squareRoot);
    BRICK_NUMERIC_DECLARE_MATHFUNCTION_FAMILY_1(tangent);

    // Declarations of two-argument overloads.
    BRICK_NUMERIC_DECLARE_MATHFUNCTION_FAMILY_2(arctangent2);

    // Declarations for functions that don't fit a common signature.
    inline void splitFraction(brick::common::Float32 arg,
                              brick::common::Float32& integerPart,
                              brick::common::Float32& fractionalPart) {
      brick::common::splitFraction(arg, integerPart, fractionalPart);
    }

    inline void splitFraction(brick::common::Float64 arg,
                              brick::common::Float64& integerPart,
                              brick::common::Float64& fractionalPart) {
      brick::common::splitFraction(arg, integerPart, fractionalPart);
    }
    
  } // namespace numeric

} // namespace brick

#endif /* #ifndef BRICK_NUMERIC_MATHFUNCTIONS_HH */
