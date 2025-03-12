/**
***************************************************************************
* @file brick/numeric/functional.hh
*
* Header file declaring numeric functors.
*
* Copyright (C) 2003-2011, David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <cmath>
#include <functional>
#include <brick/numeric/numericTraits.hh>

#ifndef BRICK_NUMERIC_FUNCTIONAL_HH
#define BRICK_NUMERIC_FUNCTIONAL_HH

namespace brick {

  namespace numeric {

    /**
     ** Functor template which computes the value of a Gaussian
     ** evaluated at its argument.
     **/
    template <class Type>
    struct Gaussian1DFunctor
    {
      /**
       * The constructor initializes the Gaussian to zero mean and the
       * specified standard deviation.
       *
       * @param sigma This argument specifies the standard deviation of
       * the Gaussian.
       */
      Gaussian1DFunctor(Type sigma = 1.0)
        : m_k0(1.0 / (std::sqrt(2.0 * M_PI) * sigma)),
          m_k1(-1.0 / (2.0 * sigma * sigma)) {}


      /**
       * Compute the value of the Gaussian at the specified input value,
       * and return the result.
       *
       * @param input The point at which to evaluate the Guassian.
       */
      inline Type
      operator()(const Type& input) {
        return m_k0 * std::exp(input * input * m_k1);
      }

    private:
      Type m_k0;
      Type m_k1;
    };


    /**
     ** Functor template which computes the natural logarithm of its
     ** argument (using std::log(), if appropriate).
     **/
    template <class Type>
    struct LogFunctor
    {
      /**
       * Compute the log of the input argument, and return the result.
       *
       * @param input The natural log of this argument will be computed.
       * @return The natural log of the argument.
       */
      inline Type
      operator()(const Type& input) {
        return std::log(input);
      }
    };


    /**
     ** Functor template which uses static_cast to convert instances of
     ** one type into instances of another, but does the right thing
     ** with regard to rounding, so that the difference between the
     ** input and the returned value is minimized.
     **/
    template <class TypeIn, class TypeOut>
    struct NumericTypeConversionFunctor
    {
      /**
       * This operator returns the output value which minimizes the
       * difference between input and output.  In general, this takes
       * some smarts, so this may need specializing later.
       *
       * @param input Will be converted to type TypeOut.
       * @return The result of the conversion.
       */
      inline TypeOut
      operator()(const TypeIn& input) {
        // This "if" should optimize away, since traits are known at compile
        // time.
        if(!NumericTraits<TypeIn>().isIntegral()
           && NumericTraits<TypeOut>().isIntegral()) {
          return static_cast<TypeOut>(input + 0.5);
        }
        return static_cast<TypeOut>(input);
      }
    };


    /**
     ** Functor template which computes the square root of its argument
     ** (using std::sqrt(), if appropriate).
     **/
    template <class Type>
    struct SquareRootFunctor
    {
      /**
       * Take the square root of the input argument, and return the result.
       *
       * @param input This argument will be passed to a square root function.
       * @return The result of the square root call.
       */
      inline Type
      operator()(const Type& input) {
        return std::sqrt(input);
      }
    };

  } // namespace numeric

}  // namespace brick


/* ================= Specializations =================== */

namespace brick {

  namespace numeric {

  // None

  } // namespace numeric

}  // namespace brick

#endif // #ifndef BRICK_NUMERIC_FUNCTIONAL_HH
