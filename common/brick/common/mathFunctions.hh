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
     * This function template takes the place of std::acos(),
     * std::acosf(), std::acosl(), etc., hopefully enabling generic
     * code.
     *
     * @param arg This argument will have its arccosine computed.
     *
     * @return The return value is the arccosine of arg, expressed in
     * radians.
     */
    template <class Type>
    inline Type arccosine(Type arg);


    /**
     * This function template takes the place of std::asin(),
     * std::asinf(), std::asinl(), etc., hopefully enabling generic
     * code.
     *
     * @param arg This argument will have its arcsine computed.
     *
     * @return The return value is the arcsine of arg, expressed in
     * radians.
     */
    template <class Type>
    inline Type arcsine(Type arg);


    /**
     * This function template takes the place of std::atan2(),
     * std::atan2f(), std::atan2l(), etc., hopefully enabling generic
     * code.
     *
     * @param yy This argument proportional to sin(theta).
     *
     * @param xx This argument proportional to cos(theta).
     *
     * @return The return value is the principal value of the arc
     * tangent of (yy/xx), expressed in radians.
     */
    template <class Type>
    inline Type arctangent2(Type yy, Type xx);


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
     * This function template takes the place of std::log(),
     * std::logf(), std::logl(), etc., hopefully enabling generic
     * code.
     *
     * @param arg This argument will have its natural logarithm
     * computed.
     *
     * @return The return value is the natural log of arg.
     */
    template <class Type>
    inline Type logarithm(Type arg);


    /**
     * This function template takes the place of std::ceil(),
     * std::ceilf(), std::ceill(), etc., hopefully enabling generic
     * code.
     *
     * @param arg This argument will be rounded down to the nearest
     * integer.
     *
     * @return The return value is the ceil of arg.
     */
    template <class Type>
    inline Type roundToCeiling(Type arg);


    /**
     * This function template takes the place of std::floor(),
     * std::floorf(), std::floorl(), etc., hopefully enabling generic
     * code.
     *
     * @param arg This argument will be rounded down to the nearest
     * integer.
     *
     * @return The return value is the floor of arg.
     */
    template <class Type>
    inline Type roundToFloor(Type arg);


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
     * This function template takes the place of std::modf(),
     * std::modff(), std::modfl(), etc., hopefully enabling generic
     * code.
     *
     * @param arg This argument is the value to be split into integer
     * and fractional parts.
     *
     * @param integerPart This argument returns the integer part of
     * the argument.
     *
     * @param fractionalPart This argument returns the fractional part
     * of the argument.
     */
    template <class Type>
    inline void splitFraction(Type arg, Type& integerPart,
                              Type& fractionalPart);


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
    inline Type arccosine(Type arg) {
      return static_cast<Type>(std::acos(static_cast<double>(arg)));
    }

    /* ======== Need to figure out how to test for cosf()  ======== */
    /* ======== and cosl() availability                    ======== */
    // template<>
    // inline float arccosine(float arg) {return std::acosf(arg);}
    //
    // template<>
    // inline long double arccosine(long double arg) {return std::acosl(arg);}


    template <class Type>
    inline Type arcsine(Type arg) {
      return static_cast<Type>(std::asin(static_cast<double>(arg)));
    }

    /* ======== Need to figure out how to test for sinf()  ======== */
    /* ======== and sinl() availability                    ======== */
    // template<>
    // inline float arcsine(float arg) {return std::asinf(arg);}
    //
    // template<>
    // inline long double arcsine(long double arg) {return std::asinl(arg);}


    template <class Type>
    inline Type arctangent2(Type yy, Type xx) {
      return static_cast<Type>(
        std::atan2(static_cast<double>(yy), static_cast<double>(xx)));
    }


    /* ======== Need to figure out how to test for atan2f()  ======== */
    /* ======== and atan2l() availability                    ======== */
    // template <>
    // inline float arctangent2(float arg) {return std::atan2f(yy, xx);}
    //
    // template <>
    // inline long double arctangent2(long double) {
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
    inline Type logarithm(Type arg) {
      return static_cast<Type>(std::log(static_cast<double>(arg)));
    }

    /* ======== Need to figure out how to test for logf()  ======== */
    /* ======== and logl() availability                    ======== */
    // template<>
    // inline float logarithm(float arg) {return std::cosf(arg);}
    //
    // template<>
    // inline long double logarithm(long double arg) {return std::cosl(arg);}


    template <class Type>
    inline Type roundToCeiling(Type arg) {
      return static_cast<Type>(std::ceil(static_cast<double>(arg)));
    }

    /* ======== Need to figure out how to test for ceilf()  ======== */
    /* ======== and ceill() availability                    ======== */
    // template<>
    // inline float roundToCeiling(float arg) {return std::ceilf(arg);}
    //
    // template<>
    // inline long double roundToCeiling(long double arg) {
    //   return std::ceill(arg);
    // }


    template <class Type>
    inline Type roundToFloor(Type arg) {
      return static_cast<Type>(std::floor(static_cast<double>(arg)));
    }

    /* ======== Need to figure out how to test for floorf()  ======== */
    /* ======== and floorl() availability                    ======== */
    // template<>
    // inline float roundToFloor(float arg) {return std::floorf(arg);}
    //
    // template<>
    // inline long double roundToFloor(long double arg) {
    //   return std::floorl(arg);
    // }


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
    inline void splitFraction(Type arg, Type& integerPart,
                              Type& fractionalPart)
    {
      // std::round() won't be available until C++11.
      integerPart = static_cast<Type>(
        (arg >= static_cast<Type>(0.0))
        ? std::floor(static_cast<double>(arg))
        : std::ceil(static_cast<double>(arg)));
      fractionalPart = arg - integerPart;
    }


    template <>
    inline void splitFraction(double arg, double& integerPart,
                              double& fractionalPart)
    {
      fractionalPart = std::modf(arg, &integerPart);
    }


    /* ======== Need to figure out how to test for modff() ======== */
    /* ======== and modfl availability                     ======== */
    // template <>
    // inline void splitFraction(float arg, float& integerPart,
    //                           float& fractionalPart)
    // {
    //   fractionalPart = std::modff(arg, &fractionalPart);
    // }
    //
    // template <>
    // inline void splitFraction(long double arg, long double& integerPart,
    //                           long double& fractionalPart)
    // {
    //   fractionalPart = std::modfl(arg, &fractionalPart);
    // }


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
