/**
***************************************************************************
* @file brick/numeric/solveQuadratic.hh
*
* Header file declaring a function template for solving quadratic
* polynomial equations of a single variable.
*
* Copyright (C) 2001-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_SOLVEQUADRATIC_HH
#define BRICK_NUMERIC_SOLVEQUADRATIC_HH

#include <complex>

namespace brick {

  namespace numeric {
    
    /** 
     * This function computes the real roots of the quadratic polynomial
     * c0*x^2 + c1*x + c2 = 0.
     *
     * Note that the two well known versions of the quadratic formula:
     * 
     *   x = (-c1 +|- sqrt(c1**2 - 4c0*c2)) / (2*c0)
     * 
     * and
     * 
     *   x = 2*c2 / (-c1 +|- sqrt(c1**2 - 4*c0*c2))
     * 
     * both tend to be inaccurate when c0 and/or c2 are small, since
     * then the quantity (-c1 +|- sqrt(c1**2 - 4*c0*c2)) gets very
     * small and loses significance.  Instead we use the form
     * advocated by Press and Flannery (Numerical Recipes):
     *
     *   q = -(1/2)(c1 + sgn(c1)*sqrt(c1**2 - 4*c1*c2))
     *   x1 = q/c0, x2 = c2/q
     *
     * This is just the same as using both of the first two forms, each
     * one for the root at which it is most numerically stable.
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
     * @param root0 If the polynomial has real roots, this reference
     * argument is used to return the first root.
     * 
     * @param root1 If the polynomial has real roots, this reference
     * argument is used to return the second root.
     * 
     * @return If the polynomial has real roots, the return value is
     * true.  If the polynomial does not have real roots, the return
     * value is false, and arguments root0 and root1 are not changed.
     */
    template <class Type>
    bool
    solveQuadratic(Type c0, Type c1, Type c2,
                   Type& root0, Type& root1);

  
    /** 
     * This function computes the (possibly complex) roots of the
     * quadratic polynomial c0*x^2 + c1*x + c2 = 0.
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
     */
    template <class Type>
    void
    solveQuadratic(Type c0, Type c1, Type c2,
                   std::complex<Type>& root0, std::complex<Type>& root1);


    /** 
     * This function computes the roots of the quadratic polynomial
     * x^2 + c0*x + c1 = 0, where c0 and c1 are complex.
     * 
     * @param c0 This argument is the linear coefficient of the
     * polynomial.
     * 
     * @param c1 This argument is the constant coefficient of the
     * polynomial.
     * 
     * @param root0 This reference argument is used to return the
     * first root of the polynomial.
     * 
     * @param root1 This reference argument is used to return the
     * second root of the polynomial.
     */
    template <class Type>
    void
    solveQuadratic(std::complex<Type> c0, std::complex<Type> c1,
                   std::complex<Type>& root0, std::complex<Type>& root1);
    
  
  } // namespace numeric

} // namespace brick

// Include file containing definitions of inline and template
// functions.
#include <brick/numeric/solveQuadratic_impl.hh>

#endif /* #ifndef BRICK_NUMERIC_SOLVEQUADRATIC_HH */
