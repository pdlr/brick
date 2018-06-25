/**
**********************************************************************
* @file brick/optimization/gradientFunction.hh
*
* Header file declaring GradientFunction class template.
*
* Copyright (C) 2003-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
**********************************************************************
**/

#ifndef BRICK_OPTIMIZATION_GRADIENTFUNCTION_HH
#define BRICK_OPTIMIZATION_GRADIENTFUNCTION_HH

#include <functional>
#include <brick/numeric/derivativeRidders.hh>

namespace brick {

  namespace optimization {

    /**
     ** The GradientFunction class template is derived from
     ** std::unary_function, and adds one additional member function
     ** for computing the gradient of the function output with respect
     ** to the argument.  If accuracy is important, consider using
     ** GradientFunctionRidders instead of GradientFunction.
     **
     ** Template argument Functor is assumed to be a subclass of
     ** std::unary_function.  Functor::result_type is assumed to be a
     ** continuous scalar type.  Functor::argument_type is assumed to be
     ** a vector or 1D array type which supports the following
     ** interface:
     **
     **   argument_type(size_t N): construct an N-element vector.
     **   size_t size(): return the number of elements in the vector.
     **   element_type& operator[](size_t i): return a reference to the
     **                        (i)th element of the array.
     **   const element_type& operator[](size_t i) const: return a
     **                        const reference to the (i)th element of
     **                        the array.
     **
     ** It is further assumed that element type of argument_type is a
     ** continuous scalar, and can be implicitly cast to and from
     ** the type specified by template argument Scalar.
     **
     ** Template argument Scalar specifies the precision with which
     ** internal calculations will be conducted.
     **
     ** Here's a usage example:
     **
     ** @code
     **   typedef GradientFunction<MyObjectiveFunction> GradFunctor;
     **   MyObjectiveFunction functor();
     **   GradFunctor gradFunctor(functor);
     **   OptimizerBFGS<GradFunctor> optimizer(gradFunctor);
     **   optimizer.setStartPoint(myStartPoint);
     **   myResult = optimizer.optimum();
     ** @endcode
     **/
    template <class Functor, class Scalar = double>
    class GradientFunction
      : public std::unary_function<typename Functor::argument_type,
                                   typename Functor::result_type>
    {
    public:
      /**
       * Constructor.
       *
       * @param functor This argument is the function object to be
       * adapted.
       *
       * @param epsilon If the gradient() method is not overridden in a
       * subclass, the gradient will be computed by using symmetric
       * divided differences with a total step size of 2 * epsilon.
       */
      GradientFunction(const Functor& functor, Scalar epsilon=1.0e-6) :
        m_functor(functor), m_epsilon(epsilon) {
        if(epsilon == 0.0) {
          BRICK_THROW(brick::common::ValueException,
		      "GradientFunction::GradientFunction()",
		      "Invalid value (0.0) for argument epsilon.");
        }
      }


      /**
       * Destructor.
       */
      virtual ~GradientFunction() {}


      /**
       * This method numerically approximates the gradient of
       * this->operator() by divided differences.  This method should
       * often be overridden by a subclass.
       *
       * This function will throw ValueException if you set
       * constructor argument epsilon small enough that, when added to
       * the elements of theta, it gets completely rounded away.
       *
       * @param theta The point around which to compute the gradient.
       * @return The computed gradient.
       */
      typename Functor::argument_type
      gradient(const typename Functor::argument_type& theta);


      /**
       * This method evaluates the function value a the specified point.
       *
       * @param theta The point at which to evaluate the function.
       * @return The function value at theta.
       */
      typename Functor::result_type
      operator()(const typename Functor::argument_type& theta) {
        return m_functor(theta);
      }

    private:
      Functor m_functor;
      Scalar m_epsilon;

    }; // class GradientFunction


    /**
     ** The GradientFunctionRidders class template is just like
     ** GradientFunction, from which it is derived, but uses Ridders's
     ** method to estimate derivatives.  This takes about 10x as long
     ** as the naive method implemented by GradientFunction, but is
     ** much more accurate.  See GradientFunction for more documentation.
     **/
    template <class Functor, class Scalar = double>
    class GradientFunctionRidders
      : public GradientFunction<Functor, Scalar>
    {
    public:

      /**
       * Constructor.
       *
       * @param functor This argument is the function object to be
       * adapted.
       *
       * @param stepBound This argument specifies the largest finite
       * difference to used in gradient() computation.  See
       * numeric::DerivativeRidders for more information.
       */
      GradientFunctionRidders(const Functor& functor,
                              Scalar stepBound = 0.01,
                              Scalar errorTolerance = 1.0e-5);


      /**
       * Destructor.
       */
      virtual ~GradientFunctionRidders() {}


      /**
       * This method numerically approximates the gradient of
       * this->operator().
       *
       * This function will throw ValueException if the estimated
       * error in the returned derivative is larger than the tolerance
       * specified in the constructor.  Usually this happens because
       * constructor argument stepBound was too large.
       */
      typename Functor::argument_type
      gradient(const typename Functor::argument_type& theta);

    private:

      typedef brick::numeric::NDimensionalFunctorAdapter<Functor, Scalar>
        FunctorAdaptor;

      Scalar m_errorTolerance;
      brick::numeric::DerivativeRidders<FunctorAdaptor> m_ridders;

      // Experimental members to allow auto-tuning of initial Ridders
      // step.
      Scalar m_maxStepBound;
      Scalar m_minStepBound;
      typename Functor::argument_type m_stepBounds;

    }; // class GradientFunctionRidders

  } // namespace optimization

} // namespace brick


/* ======= Declarations to maintain compatibility with legacy code. ======= */

namespace brick {

  using optimization::GradientFunction;

} // namespace brick


/*******************************************************************
 * Member function definitions follow.  This would be a .C file
 * if it weren't templated.
 *******************************************************************/

namespace brick {

  namespace optimization {

    // Numerically approximate the gradient of this->operator() by
    // divided differences.
    template <class Functor, class Scalar>
    typename Functor::argument_type
    GradientFunction<Functor, Scalar>::
    gradient(const typename Functor::argument_type& theta)
    {
      // Create some vectors to use as input to operator()().
      typename Functor::argument_type thetaMinus(theta.size());
      typename Functor::argument_type thetaPlus(theta.size());
      // Return value must be a vector, too, so use argument_type.
      typename Functor::argument_type result(theta.size());
      // Initialize arguments.
      for(size_t index = 0; index < theta.size(); ++index) {
        thetaMinus[index] = theta[index];
        thetaPlus[index] = theta[index];
      }
      // Now compute each partial derivative.
      for(size_t index = 0; index < theta.size(); ++index) {
        // Set up the difference.
        thetaMinus[index] = theta[index] - m_epsilon;
        thetaPlus[index] = theta[index] + m_epsilon;
        // Compute divided difference.
        typename Functor::result_type valueMinus =
          m_functor.operator()(thetaMinus);
        typename Functor::result_type valuePlus =
          m_functor.operator()(thetaPlus);

        // Dodge some roundoff error by finding out precisely what the
        // difference in arguments turned out to be.
        Scalar delta = thetaPlus[index] - thetaMinus[index];
        if(delta == 0.0) {
          BRICK_THROW(brick::common::ValueException,
		      "GradientFunction::gradient()",
		      "Difference over which gradent is computed rounded "
		      "to zero.");
        }
        result[index] = (valuePlus - valueMinus) / delta;

        // Reset argument values.
        thetaMinus[index] = theta[index];
        thetaPlus[index] = theta[index];
      }
      return result;
    }


    /* ===== Definitions for GradientFunctionRidders ==== */

    template <class Functor, class Scalar>
    GradientFunctionRidders<Functor, Scalar>::
    GradientFunctionRidders(const Functor& functor,
                            Scalar stepBound,
                            Scalar errorTolerance)
      : GradientFunction<Functor, Scalar>(functor),
        m_errorTolerance(errorTolerance),
        m_ridders(FunctorAdaptor(functor), stepBound),
        m_maxStepBound(stepBound),
        m_minStepBound(1.0e-6), // Should allow this to be set by user.
        m_stepBounds()
    {
      // Empty.
    }


    template <class Functor, class Scalar>
    typename Functor::argument_type
    GradientFunctionRidders<Functor, Scalar>::
    gradient(const typename Functor::argument_type& theta)
    {
      if(m_stepBounds.size() != theta.size()) {
        m_stepBounds = typename Functor::argument_type(theta.size());
        for(unsigned int ii = 0; ii < theta.size(); ++ii) {
          m_stepBounds[ii] = m_maxStepBound;
        }
      }

      // Return value must be a vector, so use argument_type.
      typename Functor::argument_type result(theta.size());
      typename Functor::argument_type errorEstimateVector(theta.size());
      for(unsigned int ii = 0; ii < theta.size(); ++ii) {
        result[ii] = Scalar(0.0);
        errorEstimateVector[ii] = Scalar(0.0);
      }

      // Make a local copy of theta that we can change with impunity.
      typename Functor::argument_type zeroPoint(theta.size());
      for(size_t index = 0; index < theta.size(); ++index) {
        zeroPoint[index] = theta[index];
      }

      // Compute each partial derivative.
      for(size_t index = 0; index < theta.size(); ++index) {

        // Adjust bound upward, if allowable.
        if(m_stepBounds[index] < m_maxStepBound) {
          m_stepBounds[index] *= 10.0;
          if(m_stepBounds[index] > m_maxStepBound) {
            m_stepBounds[index] = m_maxStepBound;
          }
        }

        // Loop until good result.
        while(1) {
          // We want m_ridders to evaluate gradient around theta, but
          // we can't just use theta as the zero point because that
          // would mean evaluating the adapted (1D) functor around
          // 0.0, which would fake out m_ridders when it's trying to
          // assess the scale of theta, leading to a numerically less
          // acceptable result.  For this reason, we zero out the
          // index'th element of zeroPoint, and then tell m_ridders to
          // evaluate around theta[index].
          zeroPoint[index] = 0.0;

          Scalar errorValue;
          m_ridders.getFunctor().setTargetDimension(index);
          m_ridders.getFunctor().setZeroPoint(zeroPoint);
          m_ridders.setStepBound(m_stepBounds[index]);
          result[index] = m_ridders.estimateDerivative(
            theta[index], errorValue);
          errorEstimateVector[index] = errorValue;

          // Fix our change to zeroPoint.
          zeroPoint[index] = theta[index];

          if(errorValue > m_errorTolerance) {
            m_stepBounds[index] /= 10.0;
            if(m_stepBounds[index] < m_minStepBound) {
              std::ostringstream message;
              message << "Out-of-bounds error reported by DerivativeRidders.\n"
                      << "Theta is " << theta << "\n"
                      << "Result (so far) is " << result << "\n"
                      << "StepBounds are: " << m_stepBounds << "\n"
                      << "Reported error is " << errorEstimateVector;
              BRICK_THROW(brick::common::ValueException,
			  "GradientFunctionRidders::gradient()",
			  message.str().c_str());
            }
            continue;
          }

          break;
        }
      }
      return result;
    }

  } // namespace optimization

} // namespace brick

#endif /* #ifndef BRICK_OPTIMIZATION_GRADIENTFUNCTION_HH */
