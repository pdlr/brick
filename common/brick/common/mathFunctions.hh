/**
***************************************************************************
* @file brick/common/mathFunctions.hh
*
* Header file declaring general math function templates to replace
* abs(), sin(), cos(), etc in portable code.
*
* Copyright (C) 2010, David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_COMMON_MATHFUNCTIONS_HH
#define BRICK_COMMON_MATHFUNCTIONS_HH

#include <cmath>
#include <cstdlib>

namespace brick {

  namespace common {

    /** 
     * This function template takes the place of std::abs(),
     * std::fabs(), std::fabsl(), etc., hopefully enabling generic
     * code.
     * 
     * @param arg This argument will have its absolute value computed.
     * 
     * @return The return value is the absolute value of arg.
     */
    template <class Type>
    inline Type absoluteValue(Type arg);

    
    /** 
     * This function template takes the place of std::atan2(),
     * std::atan2f(), std::atan2l(), etc., hopefully enabling generic
     * code.
     * 
     * @param yy This argument proportional to sin(theta).
     * 
     * @param xx This argument proportional to cos(theta).
     * 
     * @return The return value is the principal of the arc tangent of
     * (yy/xx).
     */
    template <class Type>
    inline Type arcTangent2(Type yy, Type xx);


    /** 
     * This function template takes the place of std::cos(),
     * std::cosf(), std::cosl(), etc., hopefully enabling generic
     * code.
     * 
     * @param arg This argument, expressed in radians, will have its
     * cosine computed.
     * 
     * @return The return value is the cosine of arg.
     */
    template <class Type>
    inline Type cosine(Type arg);


    /** 
     * This function template takes the place of std::sin(),
     * std::sinf(), std::sinl(), etc., hopefully enabling generic
     * code.
     * 
     * @param arg This argument, expressed in radians, will have its
     * sine computed.
     * 
     * @return The return value is the sine of arg.
     */
    template <class Type>
    inline Type sine(Type arg);


    /** 
     * This function template takes the place of std::sqrt(),
     * std::sqrtf(), std::sqrtl(), etc., hopefully enabling generic
     * code.
     * 
     * @param arg This argument will have its non-negative square root
     * computed.
     * 
     * @return The return value is the non-negative square root of arg.
     */
    template <class Type>
    inline Type squareRoot(Type arg);

    
    /** 
     * This function template takes the place of std::tan(),
     * std::tanf(), std::tanl(), etc., hopefully enabling generic
     * code.
     * 
     * @param arg This argument, expressed in radians, will have its
     * tangent computed.
     * 
     * @return The return value is the tangent of arg.
     */
    template <class Type>
    inline Type tangent(Type arg);

  } // namespace common

} // namespace brick


/* ==================== Implementation follows ==================== */

namespace brick {

  namespace common {

    template <class Type>
    inline Type absoluteValue(Type arg) {return (arg >= Type(0)) ? arg : -arg;}

    template<>
    inline int absoluteValue(int arg) {return std::abs(arg);}

    template<>
    inline long int absoluteValue(long int arg) {return std::labs(arg);}

    /* ======== Need to figure out how to test for llabs() ======== */
    /* ======== and imaxabs() availability                 ======== */
    // template<>
    // inline long long int absoluteValue(long long int arg) {
    //   return std::llabs(arg);
    // }
    // 
    // template<>
    // inline intmax_t absoluteValue(intmax_t arg) {return std::imaxabs(arg);}

    /* ======== Need to figure out how to test for fabsf() ======== */
    /* ======== and fabsl() availability                   ======== */
    // template<>
    // inline float absoluteValue(float arg) {return std::fabsf(arg);}

    template<>
    inline double absoluteValue(double arg) {return std::fabs(arg);}
    
    // template<>
    // inline long double absoluteValue(long double arg) {
    //   return std::fabsl(arg);
    // }

    
    template <class Type>
    inline Type arcTangent2(Type yy, Type xx) {
      return static_cast<Type>(
        std::atan2(static_cast<double>(yy), static_cast<double>(xx)));
    }


    /* ======== Need to figure out how to test for atan2f()  ======== */
    /* ======== and atan2l() availability                    ======== */
    // template <>
    // inline float arcTangent2(float arg) {return std::atan2f(yy, xx);}
    //
    // template <>
    // inline long double arcTangent2(long double) {
    //   return std::atan2l(yy, xx);
    // }


    template <class Type>
    inline Type cosine(Type arg) {
      return static_cast<Type>(std::cos(static_cast<double>(arg)));
    }

    /* ======== Need to figure out how to test for cosf()  ======== */
    /* ======== and cosl() availability                    ======== */
    // template<>
    // inline float cosine(float arg) {return std::cosf(arg);}
    //
    // template<>
    // inline long double cosine(long double arg) {return std::cosl(arg);}


    template <class Type>
    inline Type sine(Type arg) {
      return static_cast<Type>(std::sin(static_cast<double>(arg)));
    }

    /* ======== Need to figure out how to test for sinf()  ======== */
    /* ======== and sinl() availability                    ======== */
    // template<>
    // inline float sine(float arg) {return std::sinf(arg);}
    // 
    // template<>
    // inline long double sine(long double arg) {return std::sinl(arg);}


    template <class Type>
    inline Type squareRoot(Type arg) {
      return static_cast<Type>(std::sqrt(static_cast<double>(arg)));
    }

    /* ======== Need to figure out how to test for sqrtf() ======== */
    /* ======== and sqrtl availability                     ======== */
    // template<>
    // inline float squareRoot(float arg) {return std::sqrtf(arg);}
    // 
    // template<>
    // inline long double squareRoot(long double arg) {return std::sqrtl(arg);}


    template <class Type>
    inline Type tangent(Type arg) {
      return static_cast<Type>(std::tan(static_cast<double>(arg)));
    }

    /* ======== Need to figure out how to test for tanf()  ======== */
    /* ======== and tanl() availability                    ======== */
    // template<>
    // inline float tangent(float arg) {return std::tanf(arg);}
    // 
    // template<>
    // inline long double tangent(long double arg) {return std::tanl(arg);}

    
  } // namespace common

} // namespace brick


#endif /* #ifndef BRICK_COMMON_MATHFUNCTIONS_HH */
