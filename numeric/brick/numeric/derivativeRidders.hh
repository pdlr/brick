/**
***************************************************************************
* @file brick/numeric/derivativeRidders.hh
*
* Header file declaring code to approximate derivatives using Ridders's
* method.
*
* (C) Copyright 2010-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying LICENSE file for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_DERIVATIVERIDDERS_HH
#define BRICK_NUMERIC_DERIVATIVERIDDERS_HH

#include <functional>
#include <brick/numeric/array2D.hh>

namespace brick {

  namespace numeric {

    /**
     ** This class implements the algorithm described in [1] for
     ** numerically approximating first derivatives, which applies the
     ** Romberg method to numerical differentiation.  The paper also
     ** describes how to approximate higher order derivatives, but we
     ** haven't implemented that yet.  The implementation here
     ** includes termination criterion, test of estimated error, and
     ** table-based implementation from [2].
     **
     ** Template argument Functor is the the type of the functor that
     ** will be differentiated.  Template argument Scalar specifies
     ** what type will be used for internal calculations.
     **
     ** Note: we can make this run faster by following the
     ** implementation in [1].  We should be able to get better
     ** accuracy by making the finite difference step be exactly
     ** representable as described in [2], 230.  Perhaps we should
     ** build in an automatic check of estimated error, with restart,
     ** so as to minimize the number of failed estimates (see
     ** documentation for constructor parameter "stepBound").
     **
     ** C. J. F. Ridders, Accurate Computation of F'(x) and
     ** F'(x)F''(x), Advances in Engineering Software, Vol 4, No 2, p
     ** 75, 1982.
     **
     ** [2] W. H. Press et al., Numerical Recipes, The Art of
     ** Scientific Computing, Third Edition, Cambridge University
     ** Press, 2007.
     **/
    template <class Functor, class Scalar = typename Functor::argument_type>
    class DerivativeRidders {
    public:
      
      /** 
       * Constructor.
       * 
       * @param functor The "()" operator of this argument should
       * compute the function to be numerically differentiated.
       * 
       * @param stepBound This argument specifies the largest step to
       * use when computing approximate derivatives.  Setting it too
       * large will cause DerivativeRidders to terminate early with
       * large error estimate.  Setting it too small hurts accuracy.
       * Choose a value an order of magnitude or so smaller than the
       * "characteristic scale" (nominally f / f'' , corresponding to
       * the amount of parameter variation required to generate a
       * significant change in the function value) of functor.  Or
       * leave it at the default, don't panic, and check the error
       * estimate using member function getErrorEstimate().
       * 
       * @param stepDecreaseFactor This argument specifies by what
       * factor to divide the finite difference step at each
       * iteration.  Ritters' paper sets it to 2.0, Press and Flannery
       * use 1.4.  Possibly 1.4 gives better estimates, but requires
       * more iterations to terminate.
       * 
       * @param tableauSize This argument bounds how many finite
       * difference approximations will be used as input to the
       * Romberg method.  A each derivative calculation will take at
       * most twice this many function evaluations.
       * 
       * @param backtrackFactor This argument corresponds to Press's
       * and Flannery's "SAFE" constant.  Iteration will terminate
       * when estimated error is greater than or equal to this number
       * times the best error estimate seen so far.
       */
      explicit
      DerivativeRidders(Functor const& functor,
                        Scalar stepBound = 1.0e-2,
                        Scalar stepDecreaseFactor = 2.0,
                        unsigned int tableauSize = 10,
                        Scalar backtrackFactor = 2.0);

      
      /**
       * Copy constructor does a deep copy.
       */
      DerivativeRidders(DerivativeRidders<Functor, Scalar> const& other);
      
      
      /**
       * Destructor.
       */
      ~DerivativeRidders() {};


      /** 
       * Assignment operator does a deep copy.
       * 
       * @param other This argument is the instance to be copied.
       * 
       * @return The return value is a reference to *this.
       */
      DerivativeRidders<Functor, Scalar>&
      operator=(DerivativeRidders<Functor, Scalar> const& other);
      

      /** 
       * This member function estimates the first derivative of
       * functor at the specified argument value, using at most 2 *
       * tableauSize function evaluations (tableauSize is a
       * constructor argument).
       * 
       * @param argument This argument specifies at what parameter
       * value to estimate the first derivative.  This is the "x" in
       * f'(x).
       * 
       * @param errorEstimate This argument returns a rough indication
       * of the error of the returned derivative value.  It is
       * computed by comparing adjacent estimates in the Romberg
       * table.
       * 
       * @return The return value is the estimated first derivative.
       */
      Scalar
      estimateDerivative(typename Functor::argument_type const& argument,
                         Scalar& errorEstimate);


      /** 
       * This member function provides access to the saved copy of
       * functor.  It is useful in case you want to tweak some of its
       * parameters between calls to estimateDerivative().
       * 
       * @return The return value is a reference to the stored functor.
       */
      Functor&
      getFunctor() {return m_functor;}


      /** 
       * Provides retrospective access to the constructed tableau.
       * This function is provided for debugging only, and may be
       * removed once the class is stable.
       * 
       * @return The return value is a deep copy of the tableau
       * constructed in the most recent call to estimateDerivative().
       */
      numeric::Array2D<double>
      getTableau() {return m_tableau.copy();}


      void
      setStepBound(Scalar stepBound) {m_stepBound = stepBound;}
      
    private:

      Scalar m_backtrackFactor;
      Functor m_functor;
      Scalar m_rombergFactor;
      Scalar m_stepBound;
      Scalar m_stepDecreaseFactor;
      Array2D<Scalar> m_tableau;

    };
    

    /**
     ** This class takes a scalar-valued N-Dimensional function and
     ** wraps it so that you can interact with it as if it took a
     ** scalar argument.  This is useful if you want to compute its
     ** gradient using the DerivativeRidders class template.
     **
     ** Template argument Functor is the the type of the functor that
     ** will be adapted.  Template argument Scalar specifies what
     ** should normally be available as typename
     ** Functor::argument_type::element_type.  It is included in case
     ** your argument type doesn't define and element_type typedef.
     **/
    template <class Functor, class Scalar>
    class NDimensionalFunctorAdapter
      : public std::unary_function<Scalar, typename Functor::result_type>
    {
    public:
      
      /** 
       * Constructor.
       * 
       * @param functor This argument is the functor to be adapted.
       * Functor::operator()() must take exactly one argument, and
       * must not modify its argument, or this class will not work.
       * 
       * @param zeroPoint This argument specifies the point in
       * parameter space around which functor will be evaluated.
       * Calling this->operator()(0) will return the value of
       * functor(zeroPoint).
       * 
       * @param dimension This argument specifies along which
       * dimension to evaluate functor.  See the documentation for
       * NDimensionalFunctorAdapter::operator()(Scalar) for more
       * information.
       */
      NDimensionalFunctorAdapter(
        Functor const& functor,
        typename Functor::argument_type const& zeroPoint,
        unsigned int dimension);


      /** 
       * Single argument constructor.  If you use this, you must call
       * member functions setZeroPoint() and setTargetDimension()
       * before calling this->operator()().
       * 
       * @param functor This argument is the functor to be adapted.
       * Functor::operator()() must take exactly one argument, and
       * must not modify its argument, or this class will not work.
       */
      explicit
      NDimensionalFunctorAdapter(Functor const& functor);


      /** 
       * Copy constructor does a deep copy.
       * 
       * @param other This argument is the instance to be copied.
       */
      NDimensionalFunctorAdapter(
        NDimensionalFunctorAdapter<Functor, Scalar> const& other);
      
      
      /**
       * Destructor.
       */
      virtual
      ~NDimensionalFunctorAdapter() {};


      /** 
       * Assignment operator does a deep copy.
       * 
       * @param other This argument is the instance to be copied.
       * 
       * @return The return value is a reference to *this.
       */
      NDimensionalFunctorAdapter<Functor, Scalar>&
      operator=(NDimensionalFunctorAdapter<Functor, Scalar> const& other);

      
      /** 
       * This member function evaluates functor (see constructor
       * arguments) and returns the result.  Calling
       * this->operator()(x) will return the value of
       * functor(newPoint), where
       * 
       *   newPoint[i] == zeroPoint[i], for i != dimension,
       *
       * and
       *
       *   newPoint[dimension] == zeroPoint[dimension] + x.
       *
       * See constructor argument documentation for descriptions of
       * dimension and zeroPoint.
       *
       * @param argument specifies at which value to evaluate functor.
       * 
       * @return The return value is the calculated function value.
       */
      inline typename Functor::result_type
      operator()(Scalar argument);


      /** 
       * This member function specifies along which axis the
       * N-Dimensional function will be evaluated.
       * 
       * @param dimension This argument specifies along which
       * dimension to evaluate functor.  See the documentation for
       * NDimensionalFunctorAdapter::operator()(Scalar) for more
       * information.
       */
      void
      setTargetDimension(unsigned int dimension) {m_dimension = dimension;}

      
      /** 
       * This member function sets the point around which the function
       * will be evaluated.
       * 
       * @param zeroPoint This argument specifies the point in
       * parameter space around which functor will be evaluated.
       * Calling this->operator()(0) will return the value of
       * functor(zeroPoint).
       */
      void
      setZeroPoint(typename Functor::argument_type const& zeroPoint);
      
    private:

      typename Functor::argument_type m_argument;
      unsigned int m_dimension;
      Functor m_functor;
      typename Functor::argument_type m_zeroPoint;
      
    };
    
  } // namespace numeric

} // namespace brick


/*******************************************************************
 * Definitions of templates and inline functions follow.
 *******************************************************************/

#include <limits>
#include <brick/common/mathFunctions.hh>

namespace brick {

  namespace numeric {
  
    template <class Functor, class Scalar>
    DerivativeRidders<Functor, Scalar>::
    DerivativeRidders(Functor const& functor,
                      Scalar stepBound,
                      Scalar stepDecreaseFactor,
                      unsigned int tableauSize,
                      Scalar backtrackFactor)
      : m_backtrackFactor(backtrackFactor),
        m_functor(functor),
        m_rombergFactor(stepDecreaseFactor * stepDecreaseFactor),
        m_stepBound(stepBound),
        m_stepDecreaseFactor(stepDecreaseFactor),
        m_tableau(tableauSize, tableauSize)
    {
      if(stepBound == 0) {
        BRICK_THROW(common::ValueException,
                  "DerivativeRidders::DerivativeRidders()",
                  "Argument stepBound must be nonzero.");
      }
    }


    template <class Functor, class Scalar>
    DerivativeRidders<Functor, Scalar>::
    DerivativeRidders(DerivativeRidders<Functor, Scalar> const& other)
      : m_backtrackFactor(other.m_backtrackFactor),
        m_functor(other.m_functor),
        m_rombergFactor(other.m_rombergFactor),
        m_stepBound(other.m_stepBound),
        m_stepDecreaseFactor(other.m_stepDecreaseFactor),
        m_tableau(other.m_tableau.copy())
    {
      // Empty.
    }


    template <class Functor, class Scalar>
    DerivativeRidders<Functor, Scalar>&
    DerivativeRidders<Functor, Scalar>::
    operator=(DerivativeRidders<Functor, Scalar> const& other)
    {
      if(&other != this) {
        m_backtrackFactor = other.m_backtrackFactor;
        m_functor = other.m_functor;
        m_rombergFactor = other.m_rombergFactor;
        m_stepBound = other.m_stepBound;
        m_stepDecreaseFactor = other.m_stepDecreaseFactor;
        m_tableau = other.m_tableau.copy();
      }
      return *this;
    }

    
    template <class Functor, class Scalar>
    Scalar
    DerivativeRidders<Functor, Scalar>::
    estimateDerivative(typename Functor::argument_type const& argument,
                       Scalar& errorEstimate)
    {
      Scalar result;
      errorEstimate = std::numeric_limits<Scalar>::max();

      // xxx
      m_tableau = Scalar(0);
      
      // Start with the largest step.
      Scalar currentStep = m_stepBound;
      for(unsigned int ii = 0; ii < m_tableau.columns(); ++ii) {

        // Compute derivative estimate at the current step size.
        typename Functor::argument_type x0 = argument - currentStep;
        typename Functor::argument_type x1 = argument + currentStep;

        // xxx Diverging from the methods of Ridders and Press & Flannery here.
        // m_tableau(0, ii) = ((m_functor(x1) - m_functor(x0))
        //                     / (2 * currentStep));
        m_tableau(0, ii) = ((m_functor(x1) - m_functor(x0))
                            / (x1 - x0));
        
        // Now fill in table of extrapolated estimates.
        Scalar weight = m_rombergFactor;
        for(unsigned int jj = 1; jj <= ii; ++jj) {
          Scalar prevEstimateThisScale = m_tableau(jj - 1, ii);
          Scalar prevEstimatePrevScale = m_tableau(jj - 1, ii - 1);
          Scalar thisEstimate =
            ((prevEstimateThisScale * weight - prevEstimatePrevScale)
             / (weight - 1));
          m_tableau(jj, ii) = thisEstimate;

          // Extrapolations that appear to give a better estimate than
          // the current front-runner become the current front-runner.
          Scalar thisErrorEstimate = std::max(
            common::absoluteValue(thisEstimate - prevEstimateThisScale),
            common::absoluteValue(thisEstimate - prevEstimatePrevScale));
          if(thisErrorEstimate <= errorEstimate) {
            errorEstimate = thisErrorEstimate;
            result = thisEstimate;
          }

          // Next time we'll use a bigger weight to reflect that we've
          // moved one level down in the Romberg method.
          weight *= m_rombergFactor;
        }

        // Check to make sure error isn't way out of line (probably
        // due to numerical precision.  If so, quit now.
        if(ii != 0) {
          Scalar scaleErrorEstimate =
            common::absoluteValue(
              m_tableau(ii, ii) - m_tableau(ii - 1, ii - 1));
          if(scaleErrorEstimate >= m_backtrackFactor * errorEstimate) {
            break;
          }
        }
        currentStep /= m_stepDecreaseFactor;
      }
      return result;
    }


    /* == Definitions for NDimensionalFunctorAdapter == */

    template <class Functor, class Scalar>
    NDimensionalFunctorAdapter<Functor, Scalar>::
    NDimensionalFunctorAdapter(
      Functor const& functor,
      typename Functor::argument_type const& zeroPoint,
      unsigned int dimension)
      : m_argument(),
        m_dimension(0),
        m_functor(functor),
        m_zeroPoint()
    {
      this->setTargetDimension(dimension);
      this->setZeroPoint(zeroPoint);
    }

      
    template <class Functor, class Scalar>
    NDimensionalFunctorAdapter<Functor, Scalar>::
    NDimensionalFunctorAdapter(Functor const& functor)
      : m_argument(),
        m_dimension(0),
        m_functor(functor),
        m_zeroPoint()
    {
      // Empty.
    }

    
    // Copy constructor does a deep copy.
    template <class Functor, class Scalar>
    NDimensionalFunctorAdapter<Functor, Scalar>::
    NDimensionalFunctorAdapter(
      NDimensionalFunctorAdapter<Functor, Scalar> const& other)
      : m_argument(),
        m_dimension(0),
        m_functor(other.m_functor),
        m_zeroPoint()
    {
      this->setTargetDimension(other.m_dimension);
      this->setZeroPoint(other.m_zeroPoint);
    }

    
    // Assignment operator does a deep copy.
    template <class Functor, class Scalar>
    NDimensionalFunctorAdapter<Functor, Scalar>&
    NDimensionalFunctorAdapter<Functor, Scalar>::
    operator=(NDimensionalFunctorAdapter<Functor, Scalar> const& other)
    {
      if(&other != this) {
        m_functor = other.m_functor;
        this->setTargetDimension(other.m_dimension);
        this->setZeroPoint(other.m_zeroPoint);
      }
      return *this;
    }

    
    template <class Functor, class Scalar>
    typename Functor::result_type
    NDimensionalFunctorAdapter<Functor, Scalar>::
    operator()(Scalar argument) {
      m_argument[m_dimension] = m_zeroPoint[m_dimension] + argument;
      return m_functor(m_argument);
    }


    // This member function sets the point around which the function
    // will be evaluated.
    template <class Functor, class Scalar>
    void
    NDimensionalFunctorAdapter<Functor, Scalar>::
    setZeroPoint(typename Functor::argument_type const& zeroPoint)
    {
      // This looks ugly and wasteful, but it has two advantages: it
      // works for types that have shallow-copy semantics; and it
      // works for types that don't define an element_type typedef.
      m_argument = typename Functor::argument_type(zeroPoint.size());
      m_zeroPoint = typename Functor::argument_type(zeroPoint.size());
      for(unsigned int ii = 0; ii < m_zeroPoint.size(); ++ii) {
        m_argument[ii] = zeroPoint[ii];
        m_zeroPoint[ii] = zeroPoint[ii];
      }
    }
    
  } // namespace numeric

} // namespace brick

#endif /* #ifndef BRICK_NUMERIC_DERIVATIVERIDDERS_HH */
