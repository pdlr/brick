/**
***************************************************************************
* @file brick/numeric/solveQuadratic_impl.hh
*
* Implementation of function templates for solving quadratic
* polynomial equations of a single variable.
*
* Copyright (C) 2001-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

// This file is included by solveQuadratic.hh, and should not be
// directly included by user code, so no need to include
// solveQuadratic.hh here.
// 
// #include <brick/numeric/solveQuadratic.hh>

#ifndef BRICK_NUMERIC_SOLVEQUADRATIC_IMPL_HH
#define BRICK_NUMERIC_SOLVEQUADRATIC_IMPL_HH

#include <cmath>

namespace brick {

  namespace numeric {

    // This function computes the real roots of the quadratic polynomial
    // c0*x^2 + c1*x + c0 = 0.
    template <class Type>
    bool
    solveQuadratic(Type c0, Type c1, Type c2,
                   Type& root0, Type& root1)
    {
      Type ss = (c1 * c1) - (static_cast<Type>(4.0) * c0 * c2);
      if(ss < static_cast<Type>(0.0)) {
        return false;
      }
      Type rt = std::sqrt(ss);
      if(c1 < static_cast<Type>(0.0)) {
        rt = -rt;
      }
      Type qq = static_cast<Type>(-0.5) * (rt + c1);
      root0 = qq / c0;
      root1 = c2 / qq;
      return true;
    }



    // This function computes the (possibly complex) roots of the
    // quadratic polynomial c0*x^2 + c1*x + c2 = 0.
    template <class Type>
    void
    solveQuadratic(Type c0, Type c1, Type c2,
                   std::complex<Type>& root0, std::complex<Type>& root1)
    {
      Type ss = (c1 * c1) - (static_cast<Type>(4.0) * c0 * c2);
      std::complex<Type> rt;
      if(ss >= static_cast<Type>(0.0)) {
        rt.real() = std::sqrt(ss);
        rt.imag() = static_cast<Type>(0.0);
      } else {
        rt.real() = static_cast<Type>(0.0);
        rt.imag() = std::sqrt(-ss);
      }
      if(c1 < static_cast<Type>(0.0)) {
        rt = -rt;
      }
      std::complex<Type> qq = static_cast<Type>(-0.5) * (rt + c1);
      root0 = qq / c0;
      root1 = c2 / qq;
    }


    // This function computes the roots of the quadratic polynomial
    // x^2 + c0*x + c1 = 0, where c0 and c1 are complex.
    template <class Type>
    void
    solveQuadratic(std::complex<Type> c0, std::complex<Type> c1,
                   std::complex<Type>& root0, std::complex<Type>& root1)
    {
      std::complex<Type> ss = (c0 * c0) - (static_cast<Type>(4.0) * c1);
      std::complex<Type> rt = std::sqrt(ss);
      if((c0.real() * rt.real() + c0.imag() * rt.imag())
         < static_cast<Type>(0.0)) {
        rt = -rt;
      }
      std::complex<Type> qq = static_cast<Type>(-0.5) * (rt + c0);
      root0 = qq;

      // If qq is zero, then c0 and c1 are zero, and we know the
      // second root.
      if(qq.real() || qq.imag()) {
        root1 = c1 / qq;
      } else {
        root1.real() = Type(0.0);
        root1.imag() = Type(0.0);
      }
    } 
    
  } // namespace numeric

} // namespace brick

#endif /* #ifndef BRICK_NUMERIC_UTILITIES_HH */
