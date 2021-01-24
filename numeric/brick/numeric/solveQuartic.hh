/**
***************************************************************************
* @file brick/numeric/solveQuartic.hh
*
* Header file declaring a function for solving quartic polynomial
* equations of a single variable.
*
* Copyright (C) 2001-2009 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_SOLVEQUARTIC_HH
#define BRICK_NUMERIC_SOLVEQUARTIC_HH

#include <brick/common/complexNumber.hh>

namespace brick {

  namespace numeric {

    /**
     * This function computes the (possibly complex) roots of the
     * quartic polynomial x^4 + c0*x^3 + c1*x^2 + c2*x + c3 = 0.
     *
     * @param c0 This argument is the cubic  coefficient of the
     * polynomial.
     *
     * @param c1 This argument is the quadratic coefficient of the
     * polynomial.
     *
     * @param c2 This argument is the linear coefficient of the
     * polynomial.
     *
     * @param c3 This argument is the constant coefficient of the
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
     *
     * @param root3 This reference argument is used to return the
     * fourth root of the polynomial.
     */
    template <class Type>
    void
    solveQuartic(Type c0, Type c1, Type c2, Type c3,
                 brick::common::ComplexNumber<Type>& root0,
                 brick::common::ComplexNumber<Type>& root1,
                 brick::common::ComplexNumber<Type>& root2,
                 brick::common::ComplexNumber<Type>& root3);


  } // namespace numeric

} // namespace brick

// Include file containing definitions of inline and template
// functions.
#include <brick/numeric/solveQuartic_impl.hh>

#endif /* #ifndef BRICK_NUMERIC_SOLVEQUARTIC_HH */
