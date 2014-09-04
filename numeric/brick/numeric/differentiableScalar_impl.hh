/**
***************************************************************************
* @file brick/numeric/differentiableScalar_impl.hh
*
* Header file defining inline and template functions for the
* DifferentiableScalar class.
*
* Copyright (C) 2014 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_NUMERIC_DIFFERENTIABLESCALAR_IMPL_HH
#define BRICK_NUMERIC_DIFFERENTIABLESCALAR_IMPL_HH

// This file is included by differentiableScalar.hh, and should not be
// directly included by user code, so no need to include
// differentiableScalar.hh here.
// 
// #include <brick/numeric/differentiableScalar.hh>

#include <brick/common/exception.hh>
#include <brick/common/mathFunctions.hh>

namespace brick {

  namespace numeric {

    // The default constructor makes a differentiableScalar with value
    // 0, first partial derivative equal to 1, and all other partials
    // equal to 0.
    template <class Type, uint32_t Dimension>
    DifferentiableScalar<Type, Dimension>::
    DifferentiableScalar()
      : m_value(0.0),
        m_partials()
    {
      if(Dimension == 0) {
        BRICK_THROW(brick::common::ValueException,
                    "DifferentiableScalar::DifferentiableScalar()",
                    "Template argument Dimension must be greater than "
                    "or equal to 1.");
      }
      m_partials[0] = 1.0;
      std::fill(this->m_partials + 1, this->m_partials + Dimension, Type(0.0));
    }


    // This default constructor makes a differentiableScalar with a
    // user specified value, first partial derivative equal to 1, and
    // all other partials equal to 0.
    template <class Type, uint32_t Dimension>
    DifferentiableScalar<Type, Dimension>::
    DifferentiableScalar(Type value)
      : m_value(value),
        m_partials()
    {
      if(Dimension == 0) {
        BRICK_THROW(brick::common::ValueException,
                    "DifferentiableScalar::DifferentiableScalar()",
                    "Template argument Dimension must be greater than "
                    "or equal to 1.");
      }
      m_partials[0] = 1.0;
      std::fill(this->m_partials + 1, this->m_partials + Dimension, Type(0.0));
    }


    // This default constructor makes a differentiableScalar with a
    // user specified value and partial derivatives explicitly set by
    // the second constructor argument.
    template <class Type, uint32_t Dimension>
    template <class Iter>
    DifferentiableScalar<Type, Dimension>::
    DifferentiableScalar(Type value, Iter partialsBegin)
      : m_value(value),
        m_partials()
    {
      if(Dimension == 0) {
        BRICK_THROW(brick::common::ValueException,
                    "DifferentiableScalar::DifferentiableScalar()",
                    "Template argument Dimension must be greater than "
                    "or equal to 1.");
      }
      std::copy(partialsBegin, partialsBegin + Dimension,
                &(this->m_partials[0]));
    }

      
    // The copy constructor does a deep copy.
    template <class Type, uint32_t Dimension>
    DifferentiableScalar<Type, Dimension>::
    DifferentiableScalar(DifferentiableScalar<Type, Dimension> const& other)
      : m_value(other.m_value),
        m_partials()
    {
      std::copy(&(other.m_partials[0]), &(other.m_partials[0]) + Dimension,
                &(this->m_partials[0]));
    }
      

    // The assigment operator does a deep copy.
    template <class Type, uint32_t Dimension>
    DifferentiableScalar<Type, Dimension>&
    DifferentiableScalar<Type, Dimension>::
    operator=(DifferentiableScalar<Type, Dimension> const& other)
    {
      if(&other != this) {
        this->m_value = other.m_value;
        std::copy(&(other.m_partials[0]), &(other.m_partials[0]) + Dimension,
                  &(this->m_partials[0]));
      }
      return *this;
    }
      

    // This operator multiplies the differentiableScalar by another
    // differentiableScalar.
    template <class Type, uint32_t Dimension>
    DifferentiableScalar<Type, Dimension>&
    DifferentiableScalar<Type, Dimension>::
    operator*=(DifferentiableScalar<Type, Dimension> const& other)
    {
      // Apply the chain rule.
      for(uint32_t ii = 0; ii < Dimension; ++ii) {
        this->m_partials[ii] = (this->m_value * other.m_partials[ii]
                                + other.m_value * this->m_partials[ii]);
      }
      this->m_value *= other.m_value;
      return *this;
    }
    

    // This operator divides the differentiableScalar by another
    // differentiableScalar.
    template <class Type, uint32_t Dimension>
    DifferentiableScalar<Type, Dimension>&
    DifferentiableScalar<Type, Dimension>::
    operator/=(DifferentiableScalar<Type, Dimension> const& other)
    {
      // Apply the quotient rule.  Note: the quotient rule follows
      // from the chain rule.  If you start with g(x) = f(x)h(x) and
      // apply the chain rule, you can then rearrange to find the
      // derivative of f(x) = g(x)/h(x).  This gives you:
      //
      //   f'(x) = (g'(x)h(x) - g(x)h'(x)) / (h(x))^2
      Type otherSquared = other.m_value * other.m_value;
      for(uint32_t ii = 0; ii < Dimension; ++ii) {
        this->m_partials[ii] =
          (this->m_partials[ii] * other.m_value
           - this->m_value * other.m_partials[ii]) / otherSquared;
      }
      this->m_value /= other.m_value;
      return *this;
    }
    

    // This operator adds another differentiableScalar to *this.
    template <class Type, uint32_t Dimension>
    DifferentiableScalar<Type, Dimension>&
    DifferentiableScalar<Type, Dimension>::
    operator+=(DifferentiableScalar<Type, Dimension> const& other)
    {
      for(uint32_t ii = 0; ii < Dimension; ++ii) {
        this->m_partials[ii] += other.m_partials[ii];
      }
      this->m_value += other.m_value;
      return *this;
    }

    
    // This operator subtracts another differentiableScalar from *this.
    template <class Type, uint32_t Dimension>
    DifferentiableScalar<Type, Dimension>&
    DifferentiableScalar<Type, Dimension>::
    operator-=(DifferentiableScalar<Type, Dimension> const& other)
    {
      for(uint32_t ii = 0; ii < Dimension; ++ii) {
        this->m_partials[ii] -= other.m_partials[ii];
      }
      this->m_value -= other.m_value;
      return *this;
    }


    // This member function is equivalent to
    // this->setPartialDerivative(0, value).
    template <class Type, uint32_t Dimension>
    Type const&
    DifferentiableScalar<Type, Dimension>::
    setDerivative(Type const& derivative)
    {
      return m_partials[0] = derivative;
    }

      
    // This member function sets the first derivative of *this
    // with respect to the "ii"-th parameter.
    template <class Type, uint32_t Dimension>
    Type const&
    DifferentiableScalar<Type, Dimension>::
    setPartialDerivative(uint32_t ii, Type const& derivative)
    {
      return m_partials[ii] = derivative;
    }

    
    /* ============ Non-member function definitions ============ */

  
    // This operator multiplies two DifferentiableScalar instances,
    // applying the chain rule to compute the derivatives of the
    // result.
    template<class Type, uint32_t Dimension>
    DifferentiableScalar<Type, Dimension>
    operator*(DifferentiableScalar<Type, Dimension> const& arg0,
              DifferentiableScalar<Type, Dimension> const& arg1)
    {
      DifferentiableScalar<Type, Dimension> result(arg0);
      result *= arg1;
      return result;
    }


    // This operator divides two DifferentiableScalar instances,
    // applying the chain rule to compute the derivatives of the
    // result.
    template<class Type, uint32_t Dimension>
    DifferentiableScalar<Type, Dimension>
    operator/(DifferentiableScalar<Type, Dimension> const& arg0,
              DifferentiableScalar<Type, Dimension> const& arg1)
    {
      DifferentiableScalar<Type, Dimension> result(arg0);
      result /= arg1;
      return result;
    }      


    // This operator adds two DifferentiableScalar instances.
    template<class Type, uint32_t Dimension>
    DifferentiableScalar<Type, Dimension>
    operator+(DifferentiableScalar<Type, Dimension> const& arg0,
              DifferentiableScalar<Type, Dimension> const& arg1)
    {
      DifferentiableScalar<Type, Dimension> result(arg0);
      result += arg1;
      return result;
    }
  

    // This operator subtracts two DifferentiableScalar instances.
    template<class Type, uint32_t Dimension>
    DifferentiableScalar<Type, Dimension>
    operator-(DifferentiableScalar<Type, Dimension> const& arg0,
              DifferentiableScalar<Type, Dimension> const& arg1)
    {
      DifferentiableScalar<Type, Dimension> result(arg0);
      result -= arg1;
      return result;
    }


    // This operator computes the cosine of a DifferentiableScalar
    // instance, with partial derivatives.
    template<class Type, uint32_t Dimension>
    DifferentiableScalar<Type, Dimension>
    cosine(DifferentiableScalar<Type, Dimension> const& arg0)
    {
      DifferentiableScalar<Type, Dimension> result(
        brick::common::cosine(arg0.getValue()));

      // Derivative of cos(f(x)) is -sin(f(x)) * f'(x).
      Type sineValue = brick::common::sine(arg0.getValue());
      for(uint32_t ii = 0; ii < Dimension; ++ii) {
        result.setPartialDerivative(
          ii, -sineValue * arg0.getPartialDerivative(ii));
      }

      return result;
    }

    
    // This operator computes the sine of a DifferentiableScalar
    // instance, with partial derivatives.
    template<class Type, uint32_t Dimension>
    DifferentiableScalar<Type, Dimension>
    sine(DifferentiableScalar<Type, Dimension> const& arg0)
    {
      DifferentiableScalar<Type, Dimension> result(
        brick::common::sine(arg0.getValue()));

      // Derivative of sin(f(x)) is cos(f(x)) * f'(x).
      Type cosineValue = brick::common::cosine(arg0.getValue());
      for(uint32_t ii = 0; ii < Dimension; ++ii) {
        result.setPartialDerivative(
          ii, cosineValue * arg0.getPartialDerivative(ii));
      }

      return result;
    }
    
  } // namespace numeric

} // namespace brick

#endif /* #ifndef BRICK_NUMERIC_DIFFERENTIABLESCALAR_IMPL_HH */
