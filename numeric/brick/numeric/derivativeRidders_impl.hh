/**
***************************************************************************
* @file brick/numeric/derivativeRidders.hh
*
* Implementation of inline and template functions for the
* DerivativeRidders class template.
*
* (C) Copyright 2010-2012 David LaRose, dlr@cs.cmu.edu
* See accompanying LICENSE file for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_DERIVATIVERIDDERS_IMPL_HH
#define BRICK_NUMERIC_DERIVATIVERIDDERS_IMPL_HH

// This file is included by derivativeRidders.hh, and should not be
// directly included by user code, so no need to include
// derivativeRidders.hh here.
//
// #include <brick/numeric/derivativeRidders.hh>

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
      Scalar result(0);
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

#endif /* #ifndef BRICK_NUMERIC_DERIVATIVERIDDERS_IMPL_HH */
