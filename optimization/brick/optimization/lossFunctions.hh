/**
**********************************************************************
* @file brick/optimization/lossFunctions.hh
*
* Header file declaring loss function classes to enable robust
* optimization.
*
* These classes are typically used in robust estimation, where they
* implement "M-estimators."  You might set this up as a maximum likelihood
* problem framework with the model
*
* @code
*   P = Product_over_i[exp(-L(x_i / sigma_i))]
* @endcode
*
* Where L(.) is the loss function in question, the x_i are independent
* error terms, and the sigma_i are their covariances.  In this
* framework, we typically solve for the maximum likelihood estimate by
* taking the logarithm, then setting the first derivative of the
* result to zero.  In solving this equation, the first derivative of
* the loss function is important (and the actual value of L(.) isn't).
* We call the first derivative of the loss function the "weight," and
* the classes in this file expose it through member function
* getWeight().
*
* In these classes, the value of L(.) is available through member
* function getValue(), although it is often inefficient to compute.
*
* Copyright (C) 2018 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
**********************************************************************
**/

#ifndef BRICK_OPTIMIZATION_LOSSFUNCTIONS_HH
#define BRICK_OPTIMIZATION_LOSSFUNCTIONS_HH

#include <functional>

namespace brick {

  namespace optimization {

    /**
     ** Implements a loss function based on the Cauchy distribution
     ** (AKA Lorentzian distribution), which applies decreasing weight
     ** to measurements as they become more deviant.  For a loss
     ** function that completely disregards outliers, see
     ** LossFunctionTukeyBiweight.  For a loss function that applies
     ** the same weight to all outliers (in contrast to normal L2
     ** loss, which gives increasing weight to more deviant
     ** measurements), see LossFunctionHuber.
     **/
    template <class FloatType>
    class LossFunctionCauchy
      : public std::unary_function<FloatType, FloatType>
    {
    public:
      /** 
       * Constructs a loss function based on the Cauchy distribution:
       *
       * @code
       *   P(x) ~= 1 / (1 + 0.5 * (x/sigma)^2)
       * @endcode
       *
       * where x is the random variable, and sigma is the variance.
       *
       * Since we're using this in a maximum likelihood estimation
       * framework, we care about the log of this distribution.  Here
       * we'll switch variables to z, which has presumably already
       * been normalized by sigma.
       *
       * @code
       *   L(z) = ln(1 + z^2 / 2)
       * @endcode
       *
       * This function has first derivative as follows:
       *
       * @code
       *   psi(z) = z / (1 + z^2 / 2)
       * @endcode
       */
      LossFunctionCauchy();

      
      /** 
       * Applies the loss function, the function L(z) described in the
       * constructor comment.  This implements a function that is
       * convex (downward), but with less and less change as you get
       * farther from zero.  Note that it is unusual to call this
       * function.  If you are using this class in maximum likelihood
       * estimation, you probably want member function getWeight()
       * instead.
       *
       * @param argument This is the input value at which to calculate
       * the value of loss function.
       * 
       * @return The return value is the resulting value of the loss
       * function.
       */
      FloatType getValue(FloatType const& argument);


      /** 
       * Returns the first derivative of this->getValue(), the
       * function psi(z) described in the constructor comment.  This
       * is useful in Gauss-Newton iteration when setting the first
       * derivative (with respect to optimization parameters) of the
       * total error to zero.  This is sometimes called the "weight"
       * of the loss function.
       * 
       * @param argument This argument is the input value at which to
       * calculate the first derivative of the loss function.
       * 
       * @return The return value the first derivative with respect to
       * the input argument.
       */
      FloatType getWeight(FloatType const& argument);


      /**
       * Returns the square root of the loss function value.  This is
       * useful if working with code that stuffs residuals directly
       * into a sum-of-squares error function.  It lets you return
       * residuals that secretly (and inefficiently) implement an
       * M-estimator, invisible to the calling context.
       */
      FloatType getL2Equivalent(FloatType const& argument);
      
    private:
    };

    
    /**
     ** Implements the Huber loss function [1], which is quadratic for
     ** arguments with magnitude less than or equal to a user
     ** specified constant, and linear for arguments with magnitude
     ** greater than that constant.  This loss function essentially
     ** says that all deviant measurements get the same relative
     ** weight (in contrast to normal L2 loss, which gives increasing
     ** weight to more deviant measurements).  For a loss function
     ** that gives decreasing weight to more deviant measurements, see
     ** LossFunctionLorentzian.  For a loss function that completely
     ** disregards deviant measurements, see
     ** LossFunctionTukeyBiweight.
     **
     ** [1] Huber, Peter J. (1964). "Robust Estimation of a Location
     ** Parameter". Annals of Statistics. 53 (1): 73–101.
     **/
    template <class FloatType>
    class LossFunctionHuber
      : public std::unary_function<FloatType, FloatType>
    {
    public:
      /** 
       * Constructs a Huber loss function, specifying the desired
       * transition point between quadratic and linear shape.
       *
       * @code
       *   L(z) = (1/2) * z^2,             for |z| <= delta
       *          delta * (|z| - delta/2), otherwise.
       * @endcode
       *
       * @param delta This argument is the argument magnitude at which
       * the shape of the loss function switches from quadratic to
       * linear.
       */
      LossFunctionHuber(FloatType const& delta = FloatType(1.0));

      
      /** 
       * Applies the loss function.  This implements a convex
       * (downward) function that de-weights outliers when compared to
       * L2 loss.  Note that if you're using this in a maximum
       * likelihood estimation problem, you probably want member
       * function getWeight() instead.
       *
       * @param argument This is the input value at which to calculate
       * the value of loss function.
       * 
       * @return The return value is the resulting value of the loss
       * function.
       */
      FloatType getValue(FloatType const& argument);


      /** 
       * Returns the first derivative of this->getValue() with
       * respect to its parameter.  This is useful in Gauss-Newton
       * iteration when setting the first derivative (with respect to
       * optimization parameters) of the total error to zero.  This is
       * sometimes called the "weight" of the loss function.
       * 
       * @param argument This argument is the input value at which to
       * calculate the first derivative of the loss function.
       * 
       * @return The return value the first derivative with respect to
       * the input argument.
       */
      FloatType getWeight(FloatType const& argument);


      /**
       * Returns the square root of the loss function value.  This is
       * useful if working with code that stuffs residuals directly
       * into a sum-of-squares error function.  It lets you return
       * residuals that secretly (and inefficiently) implement an
       * M-estimator, invisible to the calling context.
       *
       * In the rest of this comment, f(z) is the return value of
       * this->getL2Equivalent().
       *
       * @code
       *   # The square of f(z) is the loss function itself.
       *   f(z)^2 = L(z)
       *
       *   # For |z| <= delta, the derivation is simple.
       *   f(z)^2 = (1/2) * z^2,                    |z| <= delta
       *   f(z)   = (1/sqrt(2)) * z
       *   f(z)   = (sqrt(2)/2) * z
       *
       *   # For |z| > delta, we can just punt can compute square root.
       *   f(z)^2 = delta * (|z| - delta/2),        |z| < delta
       *   f(z)   = sqrt(delta * |z| - delta^2 / 2)
       * @endcode
       */
      FloatType getL2Equivalent(FloatType const& argument);
      
    private:
      FloatType const m_delta;
      FloatType const m_deltaSquaredOverTwo;
    };
    

    /**
     ** Implements the Pseudo-Huber loss function [2], which
     ** approximates the Huber loss function, but has continuous
     ** higher-order derivatives.
     **
     ** Charbonnier, P.; Blanc-Feraud, L.; Aubert, G.; Barlaud,
     ** M. (1997). "Deterministic edge-preserving regularization in
     ** computed imaging". IEEE Trans. Image Processing. 6 (2):
     ** 298–311
     **/
    template <class FloatType>
    class LossFunctionPseudoHuber
      : public std::unary_function<FloatType, FloatType>
    {
    public:
      /** 
       * Constructs a pseudo-Huber loss function, specifying the delta
       * parameter that controlls the asymptotic shape of the
       * function.
       *
       * @code
       *   L(z) = delta^2 * (sqrt(1 + (z/delta)^2) - 1)
       * @endcode
       * 
       * @param delta This argument specifies the value of delta in
       * the equation described above.
       */
      LossFunctionPseudoHuber(FloatType const& delta = FloatType(1.0));

      
      /** 
       * Applies the loss function.  This implements a convex
       * (downward) function that de-weights outliers when compared to
       * L2 loss.  Note that you probably want member function
       * getWeight() instead.
       *
       * @param argument This is the input value at which to calculate
       * the value of loss function.
       * 
       * @return The return value is the resulting value of the loss
       * function.
       */
      FloatType getValue(FloatType const& argument);


      /** 
       * Returns the first derivative of this->getValue() with
       * respect to its parameter.  This is useful in Gauss-Newton
       * iteration when setting the first derivative (with respect to
       * optimization parameters) of the total error to zero.  This is
       * sometimes called the "weight" of the loss function.
       * 
       * @code
       *   L(z) = delta^2 * (sqrt(1 + (z/delta)^2) - 1)
       *   d/dz(L(z)) = delta^2 * d/dz(sqrt(1 + (z/delta)^2) - 1)
       *   d/dz(L(z)) = delta^2 * d/dz(sqrt(1 + (z^2 / delta^2)))
       *   d/dz(L(z)) = delta^2 * (1/2)(1 / sqrt(1 + (z^2 / delta^2))
       *                * d/dz(1 + (z^2 / delta^2))
       *   d/dz(L(z)) = delta^2 * (1/2)(1 / sqrt(1 + (z^2 / delta^2))
       *                * 2/(delta^2) * z
       *   d/dz(L(z)) = (z / sqrt(1 + (z^2 / delta^2)))
       * @endcode
       * 
       * @param argument This argument is the input value at which to
       * calculate the first derivative of the loss function.
       * 
       * @return The return value the first derivative with respect to
       * the input argument.
       */
      FloatType getWeight(FloatType const& argument);


      /**
       * Returns the square root of the loss function value.  This is
       * useful if working with code that stuffs residuals directly
       * into a sum-of-squares error function.  It lets you return
       * residuals that secretly (and inefficiently) implement an
       * M-estimator, invisible to the calling context.
       */
      FloatType getL2Equivalent(FloatType const& argument);

    private:
      FloatType const m_deltaSquared;
    };
    

    /**
     ** Implements the Tukey biweight loss function, which applies
     ** zero weight to outlier measurements.  For a loss function that
     ** gives decreasing weight to more deviant measurements -- but
     ** doesn't completely disregard them -- see
     ** LossFunctionLorentzian.  For a loss function that applies the
     ** same weight to all outliers (in contrast to normal L2 loss,
     ** which gives increasing weight to more deviant measurements),
     ** see LossFunctionHuber.
     **/
    template <class FloatType>
    class LossFunctionTukeyBiweight
      : public std::unary_function<FloatType, FloatType>
    {
    public:
      /** 
       * Constructs a Tukey biweight loss function, which has first
       * derivative as follows:
       *
       * @code
       *   psi(z) = z * (1 - z^2 / c^2)^2,     |z| <= c
       *            0                          otherwise
       * @endcode
       *
       * where c is a user supplied constant.
       * 
       * @param cc This argument specifies the scale of the loss
       * function.  Note that the optimal value of cc for normally
       * distributed errors is 6.0.

       */
      LossFunctionTukeyBiweight(FloatType const& cc = FloatType(6.0));

      
      /** 
       * Applies the loss function, This implements a function that is
       * convex (downward) between -cc and cc (remember that cc was
       * the constructor argument), and constant elsewhere.  Note that
       * it is unusual to call this function.  You probably want
       * member function getWeight() instead.
       *
       * This value was found by integrating psi(z).
       *
       * @code
       *   G(z) = z^6/(6*c^4) - z^4/(2*c^2) + z^2/2
       * @endcode
       *
       * where G(z) is the indefinite integral over z of psi(z).
       *
       * Making this definite over the range of -c to x, where x is
       * arbitrarily chosen, we have:
       *
       * @code
       *   L(x) = (x^2 - c^2)^3 / (6 * c^4) + constant
       * @endcode
       *
       * With the constraint, L(0) = 0, we easily find the constant.
       *
       * @code
       *   L(x) = (x^2 - c^2)^3 / (6 * c^4) + c^2 / 6
       * @endcode
       *
       * @param argument This is the input value at which to calculate
       * the value of loss function.
       * 
       * @return The return value is the resulting value of the loss
       * function.
       */
      FloatType getValue(FloatType const& argument);


      /** 
       * Returns the first derivative of this->getValue() with
       * respect to its parameter.  This is useful in Gauss-Newton
       * iteration when setting the first derivative (with respect to
       * optimization parameters) of the total error to zero.  This is
       * sometimes called the "weight" of the loss function.
       * 
       * @param argument This argument is the input value at which to
       * calculate the first derivative of the loss function.
       * 
       * @return The return value the first derivative with respect to
       * the input argument.
       */
      FloatType getWeight(FloatType const& argument);


      /**
       * Returns the square root of the loss function value.  This is
       * useful if working with code that stuffs residuals directly
       * into a sum-of-squares error function.  It lets you return
       * residuals that secretly (and inefficiently) implement an
       * M-estimator, invisible to the calling context.
       */
      FloatType getL2Equivalent(FloatType const& argument);
      
    private:
      FloatType const m_c;
      FloatType const m_cSquared;
    };

    
  } // namespace optimization

} // namespace brick


/*******************************************************************
 * Member function definitions follow.  This would be a .C file
 * if it weren't templated.
 *******************************************************************/

#include <brick/common/constants.hh>
#include <brick/numeric/differentiableScalar.hh>
#include <brick/numeric/mathFunctions.hh>

namespace brick {

  namespace optimization {


    // Constructs a loss function based on the Cauchy distribution.
    template <class FloatType>
    LossFunctionCauchy<FloatType>::
    LossFunctionCauchy()
      : std::unary_function<FloatType, FloatType>()
    {
      // Empty.
    }

      
    // Applies the loss function, the function L(z) described in the
    // constructor comment.
    template <class FloatType>
    FloatType
    LossFunctionCauchy<FloatType>::
    getValue(FloatType const& argument)
    {
      return brick::numeric::logarithm(
        FloatType(1.0) + argument * argument / FloatType(2.0));
    }


    // Returns the first derivative of this->getValue(), the
    // function psi(z) described in the constructor comment.
    template <class FloatType>
    FloatType
    LossFunctionCauchy<FloatType>::
    getWeight(FloatType const& argument)
    {
      return argument / (FloatType(1.0) + argument * argument / FloatType(2.0));
    }      


    // Returns the square root of the loss function value.
    template <class FloatType>
    FloatType
    LossFunctionCauchy<FloatType>::
    getL2Equivalent(FloatType const& argument)
    {
      return brick::numeric::squareRoot(this->getValue(argument));
    }

    
    // Constructs a Huber loss function, specifying the desired
    // transition point between quadratic and linear shape.
    template <class FloatType>
    LossFunctionHuber<FloatType>::
    LossFunctionHuber(FloatType const& delta)
      : std::unary_function<FloatType, FloatType>(),
        m_delta(delta),
        m_deltaSquaredOverTwo(delta * delta * FloatType(0.5))
    {
      // Empty.
    }

    
    // Applies the loss function.
    template <class FloatType>
    FloatType
    LossFunctionHuber<FloatType>::
    getValue(FloatType const& argument)
    {
      FloatType argMagnitude = brick::numeric::absoluteValue(argument);
      if(argMagnitude <= m_delta) {
        return argument * argument * FloatType(0.5);
      }
      return (m_delta * argMagnitude - m_deltaSquaredOverTwo);
    }


    // Returns the first derivative of this->getValue() with
    // respect to its parameter.
    template <class FloatType>
    FloatType
    LossFunctionHuber<FloatType>::
    getWeight(FloatType const& argument)
    {
      FloatType argMagnitude = brick::numeric::absoluteValue(argument);
      if(argMagnitude <= m_delta) {
        return argument;
      } else if(argument < FloatType(0.0)) {
        return -m_delta;
      }
      return m_delta;
    }
    

    // Returns the square root of the loss function value.
    template <class FloatType>
    FloatType
    LossFunctionHuber<FloatType>::
    getL2Equivalent(FloatType const& argument)
    {
      FloatType argMagnitude = brick::numeric::absoluteValue(argument);
      if(argMagnitude <= m_delta) {
        return brick::common::constants::rootTwoOverTwo * argument;
      }
      return brick::numeric::squareRoot(
        m_delta * argMagnitude - m_deltaSquaredOverTwo);
    }
    
    
    // Constructs a pseudo-Huber loss function, specifying the internal
    // parameter delta.
    template <class FloatType>
    LossFunctionPseudoHuber<FloatType>::
    LossFunctionPseudoHuber(FloatType const& delta)
      : std::unary_function<FloatType, FloatType>(),
        m_deltaSquared(delta * delta)
    {
      // Empty.
    }


    // Applies the loss function.
    template <class FloatType>
    FloatType
    LossFunctionPseudoHuber<FloatType>::
    getValue(FloatType const& argument)
    {
      return m_deltaSquared *
        (brick::numeric::squareRoot(
          FloatType(1.0) + argument * argument / m_deltaSquared)
         - FloatType(1.0));
    }


    // Returns the first derivative of this->getValue() with
    // respect to its parameter.
    template <class FloatType>
    FloatType
    LossFunctionPseudoHuber<FloatType>::
    getWeight(FloatType const& argument)
    {
      return argument / brick::numeric::squareRoot(
        FloatType(1.0) + argument * argument / m_deltaSquared);
    }

    
    // Returns the square root of the loss function value.
    template <class FloatType>
    FloatType
    LossFunctionPseudoHuber<FloatType>::
    getL2Equivalent(FloatType const& argument)
    {
      return brick::numeric::squareRoot(
        this->getValue(argument));
    }


    // Constructs at Tukey biweight loss function.
    template <class FloatType>
    LossFunctionTukeyBiweight<FloatType>::
    LossFunctionTukeyBiweight(FloatType const& cc)
      : std::unary_function<FloatType, FloatType>(),
        m_c(cc),
        m_cSquared(cc * cc)
    {
      // Empty.
    }
      

    // Applies the loss function, This implements a function that is
    // convex (downward) between -cc and cc (remember that cc was
    // the constructor argument), and constant elsewhere.
    template <class FloatType>
    FloatType
    LossFunctionTukeyBiweight<FloatType>::
    getValue(FloatType const& argument)
    {
      FloatType integrationConstant = m_cSquared / FloatType(6.0);
      FloatType argMagnitude = brick::numeric::absoluteValue(argument);
      if(argMagnitude <= m_c) {
        FloatType x2MinusC2 = argument * argument - m_cSquared;
        FloatType numerator = x2MinusC2 * x2MinusC2 * x2MinusC2;
        FloatType denominator = FloatType(6) * m_cSquared * m_cSquared;
        return numerator / denominator + integrationConstant;
      }
      return integrationConstant;
    }


    // Returns the first derivative of this->getValue() with
    // respect to its parameter.
    template <class FloatType>
    FloatType
    LossFunctionTukeyBiweight<FloatType>::
    getWeight(FloatType const& argument)
    {
      FloatType argMagnitude = brick::numeric::absoluteValue(argument);
      if(argMagnitude <= m_c) {
        FloatType temp0 = FloatType(1.0) - argument * argument / m_cSquared;
        return argument * temp0 * temp0;
      }
      return FloatType(0.0);
    }


    template <class FloatType>
    FloatType
    LossFunctionTukeyBiweight<FloatType>::
    getL2Equivalent(FloatType const& argument)
    {
      return brick::numeric::squareRoot(this->getValue(argument));
    }
    
  } // namespace optimization

} // namespace brick

#endif /* #ifndef BRICK_OPTIMIZATION_LOSSFUNCTIONS_HH */
