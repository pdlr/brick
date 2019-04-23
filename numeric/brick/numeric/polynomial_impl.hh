/**
***************************************************************************
* @file brick/numeric/polynomial_impl.hh
*
* Header file defining inline and template functions for the
* Polynomial class.
*
* Copyright (C) 2006,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_NUMERIC_POLYNOMIAL_IMPL_HH
#define BRICK_NUMERIC_POLYNOMIAL_IMPL_HH

// This file is included by polynomial.hh, and should not be directly included
// by user code, so no need to include polynomial.hh here.
//
// #include <brick/numeric/polynomial.hh>

#include <algorithm>

namespace brick {

  namespace numeric {


    // The default constructor makes a constant polynomial: p(x) = 1.
    template <class Type>
    Polynomial<Type>::
    Polynomial()
      : m_coefficientArray(1)
    {
      m_coefficientArray[0] = Type(1.0);
    }


    // This constructor makes a constant polynomial:
    template <class Type>
    Polynomial<Type>::
    Polynomial(Type coefficient0)
      : m_coefficientArray(1)
    {
      m_coefficientArray[0] = coefficient0;
    }


    // This constructor makes a first order polynomial:
    template <class Type>
    Polynomial<Type>::
    Polynomial(Type coefficient1, Type coefficient0)
      : m_coefficientArray(2)
    {
      m_coefficientArray[0] = coefficient0;
      m_coefficientArray[1] = coefficient1;
    }


    // This constructor makes a first order polynomial:
    template <class Type>
    Polynomial<Type>::
    Polynomial(Type coefficient2, Type coefficient1, Type coefficient0)
      : m_coefficientArray(3)
    {
      m_coefficientArray[0] = coefficient0;
      m_coefficientArray[1] = coefficient1;
      m_coefficientArray[2] = coefficient2;
    }


    // This constructor makes a polynomial of arbitrary order.
    template <class Type>
    Polynomial<Type>::
    Polynomial(Array1D<Type> const& coefficients)
      : m_coefficientArray(coefficients.copy())
    {
      // Empty.
    }


#if BRICK_NUMERIC_POLYNOMIAL_TEMPLATED_CONSTRUCTORS
    //  This constructor makes a polynomial of arbitrary order.
    template <class Type>
    template <class Iter>
    Polynomial<Type>::
    Polynomial(Iter beginIter, Iter endIter, bool isSequence)
      : m_coefficientArray(endIter - beginIter)
    {
      std::copy(beginIter, endIter, m_coefficientArray.begin());
    }
#endif /* #if BRICK_NUMERIC_POLYNOMIAL_TEMPLATED_CONSTRUCTORS */

    // The copy constructor does a deep copy.
    template <class Type>
    Polynomial<Type>::
    Polynomial(Polynomial<Type> const& other)
      : m_coefficientArray(other.m_coefficientArray.copy())
    {
      // Empty.
    }


    // This member function returns the coefficients of the
    // polynomial.
    template <class Type>
    Array1D<Type>
    Polynomial<Type>::
    getCoefficientArray() const
    {
      return m_coefficientArray.copy();
    }


    // The assigment operator does a deep copy.
    template <class Type>
    Polynomial<Type>&
    Polynomial<Type>::
    operator=(Polynomial<Type> const& other)
    {
      if(&other != this) {
        m_coefficientArray = other.m_coefficientArray.copy();
      }
      return *this;
    }


    // This operator evaluates the polynomial.
    template <class Type>
    Type
    Polynomial<Type>::
    operator()(Type xValue) const
    {
      Type result = m_coefficientArray[0];
      if(m_coefficientArray.size() > 1) {
        result += xValue * m_coefficientArray[1];
      }
      Type xToTheN = xValue;
      for(size_t index0 = 2; index0 < m_coefficientArray.size(); ++index0) {
        xToTheN *= xValue;
        result += xToTheN * m_coefficientArray[index0];
      }
      return result;
    }


    // This operator multiplies the polynomial by another polynomial.
    template <class Type>
    Polynomial<Type>&
    Polynomial<Type>::
    operator*=(Polynomial<Type> const& other)
    {
      size_t oldOrder = this->getOrder();
      size_t otherOrder = other.getOrder();
      size_t newOrder = oldOrder + otherOrder;
      Array1D<Type> newCoefficientArray(newOrder + 1);
      newCoefficientArray = Type(0);

      // We're going to use oldOrder and otherOrder as bounds for the
      // loops which calculate the new coefficients.  Increment them
      // now to avoid writing "oldOrder + 1" in the loop termination
      // condition.
      ++oldOrder;
      ++otherOrder;
      for(size_t index0 = 0; index0 < oldOrder; ++index0) {
        for(size_t index1 = 0; index1 < otherOrder; ++index1) {
          newCoefficientArray[index0 + index1] +=
            m_coefficientArray[index0] * other.m_coefficientArray[index1];
        }
      }
      m_coefficientArray = newCoefficientArray;
      return *this;
    }


    // This operator adds another polynomial to *this.
    template <class Type>
    Polynomial<Type>&
    Polynomial<Type>::
    operator+=(Polynomial<Type> const& other)
    {
      size_t largerSize = m_coefficientArray.size();
      size_t smallerSize = other.m_coefficientArray.size();
      if(smallerSize > largerSize) {
        std::swap(smallerSize, largerSize);
      }
      Array1D<Type> newCoefficientArray(largerSize);

      size_t index0 = 0;
      while(index0 < smallerSize) {
        newCoefficientArray[index0] =
          m_coefficientArray[index0] + (other.m_coefficientArray)[index0];
        ++index0;
      }
      while(index0 < m_coefficientArray.size()) {
        newCoefficientArray[index0] = m_coefficientArray[index0];
        ++index0;
      }
      while(index0 < other.m_coefficientArray.size()) {
        newCoefficientArray[index0] = other.m_coefficientArray[index0];
        ++index0;
      }
      m_coefficientArray = newCoefficientArray;
      return *this;
    }


    // This operator subtracts another polynomial from *this.
    template <class Type>
    Polynomial<Type>&
    Polynomial<Type>::
    operator-=(Polynomial<Type> const& other)
    {
      size_t largerSize = m_coefficientArray.size();
      size_t smallerSize = other.m_coefficientArray.size();
      if(smallerSize > largerSize) {
        std::swap(smallerSize, largerSize);
      }
      Array1D<Type> newCoefficientArray(largerSize);

      size_t index0 = 0;
      while(index0 < smallerSize) {
        newCoefficientArray[index0] =
          m_coefficientArray[index0] - (other.m_coefficientArray)[index0];
        ++index0;
      }
      while(index0 < m_coefficientArray.size()) {
        newCoefficientArray[index0] = m_coefficientArray[index0];
        ++index0;
      }
      while(index0 < other.m_coefficientArray.size()) {
        newCoefficientArray[index0] = -(other.m_coefficientArray[index0]);
        ++index0;
      }
      m_coefficientArray = newCoefficientArray;
      return *this;
    }



    /* ============ Non-member function definitions ============ */

    // This operator multiplies two Polynomial instances.
    template<class Type>
    Polynomial<Type>
    operator*(Polynomial<Type> const& arg0, Polynomial<Type> const& arg1)
    {
      Polynomial<Type> result = arg0;
      result *= arg1;
      return result;
    }


    // This operator adds two Polynomial instances.
    template<class Type>
    Polynomial<Type>
    operator+(Polynomial<Type> const& arg0, Polynomial<Type> const& arg1)
    {
      Polynomial<Type> result = arg0;
      result += arg1;
      return result;
    }


    // This operator subtracts two Polynomial instances.
    template<class Type>
    Polynomial<Type>
    operator-(Polynomial<Type> const& arg0, Polynomial<Type> const& arg1)
    {
      Polynomial<Type> result = arg0;
      result -= arg1;
      return result;
    }


  } // namespace numeric

} // namespace brick

#endif /* #ifndef BRICK_NUMERIC_POLYNOMIAL_IMPL_HH */
