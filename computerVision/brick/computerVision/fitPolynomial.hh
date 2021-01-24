/**
***************************************************************************
* @file brick/computerVision/fitPolynomial.hh
*
* Header file declaring functions for estimating polynomial coefficients.
*
* Copyright (C) 2019 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
***************************************************************************
**/

#ifndef BRICK_COMPUTERVISION_FITPOLYNOMIAL_HH
#define BRICK_COMPUTERVISION_FITPOLYNOMIAL_HH

#include <brick/numeric/polynomial.hh>

namespace brick {

  namespace computerVision {

    /**
     * This function computes the best polynomial fit between the two input
     * sequences.  That is, it solves for scalars c_0, c_1, ... that minimize
     * the quantity
     *
     *   e = sum_over_j((sum_over_i(c_i * x_j^i) - y_j)^2),
     *
     * where x_j and y_j are the elements of the two input sequences.  The
     * resulting coeffients are returned as a numeric::Polynomial.  Put
     * another way, this function solves for the polynomial P() that
     * most nearly (in the least squares sense) satisfies the system
     * of equations:
     *
     *   P(x_j) = y_j
     *
     * @param xBegin This iterator, along with argument xEnd, specify the
     * sequence of "X" input values.
     *
     * @param xEnd This iterator, along with argument xBegin, specify the
     * sequence of "X" input values.
     *
     * @param yBegin This iterator specifies the beginning of a sequence of
     * "Y" input values.  If there are fewer elements in this sequence than
     * there are in the X input sequence, then the result is undefined.
     *
     * @param order This argument specifies the order of the returned
     * polynomial.  Set this to 1 for linear, 2 for quadratic, etc.
     *
     * @return The return value is the best fit polynomial.
     */
    template <class FloatType, class Iter0Type, class Iter1Type>
    brick::numeric::Polynomial<FloatType>
    fitPolynomial(Iter0Type xBegin, Iter0Type xEnd, Iter1Type yBegin,
                  std::size_t order);

  } // namespace computerVision

} // namespace brick

// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/fitPolynomial_impl.hh>

#endif // #ifndef BRICK_COMPUTERVISION_FITPOLYNOMIAL_HH
