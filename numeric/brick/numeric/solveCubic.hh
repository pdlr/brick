/**
***************************************************************************
* @file brick/numeric/solveCubic.hh
*
* Header file declaring a function template for solving cubic
* polynomial equations of a single variable.
*
* Copyright (C) 2001-2009,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_SOLVECUBIC_HH
#define BRICK_NUMERIC_SOLVECUBIC_HH

#include <brick/common/complexNumber.hh>

namespace brick {

  namespace numeric {

    /**
     * This function computes the real roots of the cubic polynomial
     * x^3 + c0*x^2 + c1*x + c2 = 0.
     *
     * @param c0 This argument is the quadratic coefficient of the
     * polynomial.
     *
     * @param c1 This argument is the linear coefficient of the
     * polynomial.
     *
     * @param c2 This argument is the constant coefficient of the
     * polynomial.
     *
     * @param root0 This reference argument is used to return the
     * first real root of the polynomial.
     *
     * @param root1 If the polynomial has three real roots, this
     * reference argument is used to return the second root.
     *
     * @param root2 If the polynomial has three real roots, this
     * argument is used to return the third root.
     *
     * @return If the polynomial has three real roots, the return
     * value is true.  If the polynomial has only one real root, the
     * return value is false, and arguments root1 and root2 are not
     * changed.
     */
    template <class Type>
    bool
    solveCubic(Type c0, Type c1, Type c2,
               Type& root0, Type& root1, Type& root2);


    /**
     * This function computes the (possibly complex) roots of the
     * cubic polynomial x^3 + c0*x^2 + c1*x + c2 = 0.  As of mid-2015,
     * the std::complex API is not stable, so we use
     * brick::common::ComplexNumber instead.
     *
     * @param c0 This argument is the quadratic coefficient of the
     * polynomial.
     *
     * @param c1 This argument is the linear coefficient of the
     * polynomial.
     *
     * @param c2 This argument is the constant coefficient of the
     * polynomial.
     *
     * @param root0 This reference argument is used to return the
     * first root of the polynomial.
     *
     * @param root1 This reference argument is used to return the
     * second root of the polynomial.
     *
     * @param root2 This reference argument is used to return the
     * third root of the polynomial.
     */
    template <class Type>
    void
    solveCubic(Type c0, Type c1, Type c2,
               brick::common::ComplexNumber<Type>& root0,
               brick::common::ComplexNumber<Type>& root1,
               brick::common::ComplexNumber<Type>& root2);


  } // namespace numeric

} // namespace brick

// Include file containing definitions of inline and template
// functions.
#include <brick/numeric/solveCubic_impl.hh>

#endif /* #ifndef BRICK_NUMERIC_SOLVECUBIC_HH */
