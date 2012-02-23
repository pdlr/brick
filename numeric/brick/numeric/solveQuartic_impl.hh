/**
***************************************************************************
* @file brick/numeric/solveQuartic_impl.hh
*
* Implementation of function templates for solving quartic polynomial
* equations of a single variable.
*
* Copyright (C) 2001-2009,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_SOLVEQUARTIC_IMPL_HH
#define BRICK_NUMERIC_SOLVEQUARTIC_IMPL_HH

// This file is included by solveQuartic.hh, and should not be
// directly included by user code, so no need to include
// solveQuartic.hh here.
// 
// #include <brick/numeric/solveQuartic.hh>

#include <cmath>
#include <brick/common/constants.hh>
#include <brick/numeric/solveCubic.hh>
#include <brick/numeric/solveQuadratic.hh>

namespace brick {

  namespace numeric {

    // This function computes the (possibly complex) roots of the
    // quartic polynomial x^4 + c0*x^3 + c1*x^2 + c2*x + c3 = 0.
    template <class Type>
    void
    solveQuartic(Type c0, Type c1, Type c2, Type c3,
                 std::complex<Type>& root0, std::complex<Type>& root1,
                 std::complex<Type>& root2, std::complex<Type>& root3)
    {
      // This solution method follows the "Quick and memorable
      // solution from first principles" algorithm contributed to
      // http://en.wikipedia.org/wiki/Quartic_equation by user
      // Goatchurch on Nov. 6, 2006.
      //
      // Start with the change of variables x = u - c0/4, to get a
      // depressed quartic:
      //
      //   u^4 + alpha*u^2 + beta*u + gamma = 0,
      //
      //   alpha = (-3*c0^2 / 8) + c1,
      //   beta = (c0^3 / 8) - ((c0 * c1) / 2) + c2,
      //   gamma = (-3*c0^4 / 256) + ((c0^2 * c1) / 16) - ((c0 * c2) / 4) + c3.

      Type c0Squared = c0 * c0;
      Type c0Cubed = c0Squared * c0;
      Type minus3C0Squared = Type(-3.0) * c0Squared;

      Type alpha = (minus3C0Squared / Type(8.0)) + c1;
      Type beta = (c0Cubed / Type(8.0)) - ((c0 * c1) / Type(2.0)) + c2;
      Type gamma = (((minus3C0Squared * c0Squared) / Type(256.0))
                    + ((c0Squared * c1) / Type(16.0)) - ((c0 * c2) / Type(4.0))
                    + c3);

      // We'll factor the depressed quartic into the product of two
      // quadratics:
      //
      //   (u^2 + p*u + q) * (u^2 + r*u + s) = 0,
      //
      // where r = -p (so that the cubed term drops out).
      //
      // Setting the two polynomials in u equal, we have:
      //
      //   alpha = s - p^2 + q
      //   beta = p * (s - q)
      //   gamma = sq
      // 
      // From which we get:
      //
      //   alpha + p^2 = s + q
      //   (alpha + p^2)^2 = (s + q)^2
      // 
      //   beta / p = s - q
      //   (beta / p)^2 = (s - q)^2
      //
      //   (alpha + p^2)^2 - (beta / p)^2 = (s + q)^2 - (s - q)^2
      //   = 4sq = 4*gamma
      //
      // Substituting P = p^2 gives:
      //
      //   P^2 + 2*alpha*P + alpha^2 - (beta^2 / P) = 4*gamma
      // 
      // or equivalently
      // 
      //   P^3 + 2*alpha*P^2 +(alpha^2 - 4*gamma)*P - beta^2 = 0
      //
      // which we solve for P.

      Type k0 = Type(2.0) * alpha;
      Type k1 = alpha * alpha - Type(4.0) * gamma;
      Type k2 = -(beta * beta);
      std::complex<Type> r0;
      std::complex<Type> r1;
      std::complex<Type> r2;
      solveCubic(k0, k1, k2, r0, r1, r2);

      // Pick the largest root of the polynomial.
      Type mag2R0 = r0.real() * r0.real(); // R1 is always real.
      Type mag2R1 = r1.real() * r1.real() + r1.imag() * r1.imag();
      Type mag2R2 = r2.real() * r2.real() + r2.imag() * r2.imag();
      if(mag2R1 > mag2R0) {
        std::swap(r1, r0);
        std::swap(mag2R1, mag2R0);
      }
      if(mag2R2 > mag2R0) {
        std::swap(r2, r0);
        std::swap(mag2R2, mag2R0);
      }

      // Recover p.
      std::complex<Type> pp = std::sqrt(r0);

      // Recover s and q from the equations
      //
      //   alpha + p^2 = s + q
      //   beta / p = s - q
      std::complex<Type> alphaPlusPpSquared = alpha + pp * pp;
      std::complex<Type> betaOverPp =
        mag2R0 ? (beta / pp) : std::complex<Type>(Type(0.0), Type(0.0));
      std::complex<Type> ss = (alphaPlusPpSquared + betaOverPp) / Type(2.0);
      std::complex<Type> qq = (alphaPlusPpSquared - betaOverPp) / Type(2.0);

      // Now find the roots of the polynomial in u, and change
      // variables back to x:
      // 
      //   (u^2 + p*u + q) * (u^2 + r*u + s) = 0,
      //   x = u - c0/4
      
      solveQuadratic(pp, qq, r0, r1);
      root0 = r0 - c0 / Type(4.0);
      root1 = r1 - c0 / Type(4.0);

      solveQuadratic(-pp, ss, r0, r1);
      root2 = r0 - c0 / Type(4.0);
      root3 = r1 - c0 / Type(4.0);
    }

  } // namespace numeric

} // namespace brick

#endif /* #ifndef BRICK_NUMERIC_SOLVEQUARTIC_IMPL_HH */
