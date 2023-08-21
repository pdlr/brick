/**
***************************************************************************
* @file brick/numeric/solveCubic_impl.hh
*
* Implementation of functions for solving cubic polynomial equations
* of a single variable.
*
* Copyright (C) 2001-2009,2012 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_SOLVECUBIC_IMPL_HH
#define BRICK_NUMERIC_SOLVECUBIC_IMPL_HH

// This file is included by solveCubic.hh, and should not be
// directly included by user code, so no need to include
// solveCubic.hh here.
//
// #include <brick/numeric/solveCubic.hh>

#include <cmath>
#include <brick/common/constants.hh>
#include <brick/numeric/mathFunctions.hh>
#include <brick/numeric/solveQuadratic.hh>

#include <iostream>
#include <vector>

namespace brick {

  namespace numeric {

    // This function computes the real roots of the cubic polynomial
    // x^3 + c0*x^2 + c1*x + c2 = 0.
    template <class Type>
    bool
    solveCubic(Type c0, Type c1, Type c2,
               Type& root0, Type& root1, Type& root2)
    {
      // We follow the formulation in Press et al, "Numerical Recipes,
      // The Art of Scientific Computing," third edition, Cambridge
      // University Press, 2007.

      bool returnValue = true;

      Type c0Squared = c0 * c0;
      Type qq = ((c0Squared - (Type(3.0) * c1)) / Type(9.0));
      Type rr = ((Type(2.0) * c0Squared * c0
                  - Type(9.0) * c0 * c1
                  + Type(27.0) * c2)
                 / Type(54.0));

      Type rrSquared = rr * rr;
      Type qqCubed = qq * qq * qq;
      Type c0OverThree = c0 / Type(3.0);
      if(rrSquared < qqCubed) {
        // Looks like we have three real roots.
        Type theta = std::acos(rr / squareRoot(qqCubed));
        Type minusTwoRootQq = Type(-2.0) * Type(squareRoot(qq));
        Type twoPi = Type(2.0 * brick::common::constants::pi);

        root0 = (minusTwoRootQq * std::cos(theta / Type(3.0))
                 - c0OverThree);
        root1 = (minusTwoRootQq * std::cos((theta + twoPi) / Type(3.0))
                 - c0OverThree);
        root2 = (minusTwoRootQq * std::cos((theta - twoPi) / Type(3.0))
                 - c0OverThree);
      } else {
        // Looks like we have some complex roots.
        bool signRr = rr > Type(0.0);
        Type absRr = signRr ? rr : -rr;
        Type aa = std::pow(absRr + squareRoot(rrSquared - qqCubed), 1.0 / 3.0);
        if(signRr) {
          aa = -aa;
        }

        Type bb = (aa == Type(0.0)) ? Type(0.0) : (qq / aa);

        root0 = (aa + bb) - c0OverThree;
        returnValue = false;
      }
      return returnValue;
    }



    // This function computes the (possibly complex) roots of the
    // cubic polynomial x^3 + c0*x^2 + c1*x + c2 = 0.
    template <class Type>
    void
    solveCubic(Type c0, Type c1, Type c2,
               brick::common::ComplexNumber<Type>& root0,
               brick::common::ComplexNumber<Type>& root1,
               brick::common::ComplexNumber<Type>& root2)
    {
      // We follow the formulation in Press et al, "Numerical Recipes,
      // The Art of Scientific Computing," third edition, Cambridge
      // University Press, 2007.
      Type c0Squared = c0 * c0;
      Type qq = ((c0Squared - (Type(3.0) * c1)) / Type(9.0));
      Type rr = ((Type(2.0) * c0Squared * c0
                  - Type(9.0) * c0 * c1
                  + Type(27.0) * c2)
                 / Type(54.0));

      Type rrSquared = rr * rr;
      Type qqCubed = qq * qq * qq;
      Type c0OverThree = c0 / Type(3.0);
      if(rrSquared < qqCubed) {
        // Looks like we have three real roots.
        Type theta = std::acos(rr / squareRoot(qqCubed));
        Type minusTwoRootQq = Type(-2.0) * Type(squareRoot(qq));
        Type twoPi = Type(2.0 * brick::common::constants::pi);

        root0.setValue((minusTwoRootQq * std::cos(theta / Type(3.0))
                        - c0OverThree),
                       0.0);
        root1.setValue((minusTwoRootQq * std::cos((theta + twoPi) / Type(3.0))
                        - c0OverThree),
                       0.0);
        root2.setValue((minusTwoRootQq * std::cos((theta - twoPi) / Type(3.0))
                        - c0OverThree),
                       0.0);
      } else {
        // Looks like we have some complex roots.
        bool signRr = rr > Type(0.0);
        Type absRr = signRr ? rr : -rr;
        Type aa = std::pow(absRr + squareRoot(rrSquared - qqCubed), 1.0 / 3.0);
        if(signRr) {
          aa = -aa;
        }

        Type bb = (aa == Type(0.0)) ? Type(0.0) : (qq / aa);

        root0.setValue(((aa + bb) - c0OverThree), 0.0);
        root1.setValue((Type(-0.5) * (aa + bb) - c0OverThree),
                       (Type(squareRoot(3.0)) / Type(2.0)) * (aa - bb));
        root2.setValue(root1.getRealPart(), -(root1.getImaginaryPart()));
      }
    }

    // This function is not part of the public interface, and is not tested.
    template<class Type>
    Type
    cubify(Type const& c0, Type const& c1, Type const& c2,
           Type const& c3, Type const& xx)
    {
      return (c0 * xx * xx * xx + c1 * xx * xx + c2 * xx + c3);
    }

    // This function is not part of the public interface, and is not tested.
    // It is part of an effort to build a solver that accurately finds the
    // roots of cubic polynomials with very small initial coefficients.  But
    // that effort hasn't run to completion yet.
    //
    // This function iteratively computes the (possibly complex) roots
    // of the cubic polynomial c0*x^3 + c1*x^2 + c2*x + c3 = 0.
    template <class Type>
    void
    solveCubicIterative(Type const& c0, Type const& c1, Type const& c2,
                        Type const& c3, Type const& /* epsilon */,
                        brick::common::ComplexNumber<Type>& root0,
                        brick::common::ComplexNumber<Type>& root1,
                        brick::common::ComplexNumber<Type>& root2)
    {
      // First let's figure some things out about the cubic.
      // If f(x) = c0*x^3 + c1*x^2 + c2*x + c3
      // then f'(x) = 3*c0*x^2 + 2*c1*x + c2
      //      f''(x) = 6*c0*x + 2*c1
      //      f'''(x) = 6*c0
      //
      // If c0 > 0, then f''(x) slopes from bottom left to top right, and
      // f'(x) is convex (i.e., concave up).  Similarly, if c0 < 0, then
      // f'(x) is concave down.
      //
      // The central extremum of f'(x) happens when f'' is 0.  This is at
      //   x_c = -c1 / (3*c0).
      //
      // At this point, the value of f'(x) is
      //   f'(x_c) = -c1^2/(3*c0) + c2
      //
      // This means if c0 and (-c1^2/(3*c0) + c2) have the same sign, then
      // f'(x) never crosses the X axis, and has no real roots.  This means
      // that f(x) is monotonically increasing or decreasing, and has one
      // real root.
      
      // Assuming small value for c0, find candidate roots by neglecting
      // the cubic term.  Warning(xxx): This is dangerous.
      brick::common::ComplexNumber<Type> candidate0;
      brick::common::ComplexNumber<Type> candidate1;
      solveQuadratic(c1, c2, c3, candidate0, candidate1);

      candidate0.setValue(-2020275687.0, 0.0);

      std::cout << "Candidate0: " << candidate0 << " -- "
                << cubify(c0, c1, c2, c3, candidate0.getRealPart())
                << std::endl;
      std::cout << "Candidate1: " << candidate1 << " -- "
                << cubify(c0, c1, c2, c3, candidate1.getRealPart())
                << std::endl;

      std::vector<brick::common::ComplexNumber<Type> > candidateVector;
      candidateVector.push_back(candidate0);
      candidateVector.push_back(candidate1);

      for(auto x0 : candidateVector) {
        // Change variables to expand around this root.  That is, we
        // write a polynomial in a variable, epsilon, such that:
        //
        //   k0 * epsilon^3 + k1 * epsilon^2 + k2 * epsilon + k3
        //   = (c0 * (x0 + epsilon)^3 + c1 * (x0 + epsilon)^2
        //      + c2 * (x0 + epsilon) + c3.
        //
        // We could work out all the algebra to solve for k0, k1,
        // etc., or we could just use Taylor series expansion.
        //
        //  f(x0) = c0 * x0^3 /* Because c1 * x0^2 + c2 * x0 + c3 = 0. */
        //  f'(x0) = 3*c0*x0^2 + 2*c1*x0 + c2
        //  f''(x0) = 6*c0*x0 + 2*c1
        //  f'''(x0) = 6*c0
        //
        //  f(x0 + epsilon) = (f(x0) + epsilon * f'(x0)
        //                     + (epsilon^2/2) * f''(x0)
        //                     + (epsilon^3/6) * f'''(x0))
        //                  = (c0*x0^3 + epsilon * (3*c0*x0^2 + 2*c1*x0 + c2)
        //                     + (epsilon^2/2) * (6*c0*x0 + 2*c1)
        //                     + (epsilon^3/6) * (6*c0)
        //                  = (c0*epsilon^3 + (3*c0*x0 + c1)*epsilon^2
        //                     + (3*c0*x0^2 + 2*c1*x0 + c2)*epsilon
        //                     + c0*x0^3

        // The roots of the cubic equation in epsilon correspond to
        // roots of the input cubic.  But the equation in epsilon has
        // the same leading coefficient, so it's subject to the same
        // numerical concerns.  However, for the root(s) with a small
        // value of We have a very strong argument for neglecting the
        // cubic term:
        //   - We deliberately set the problem up so that epsilon is small.
        //   - We multiply this tiny value by c0, which is is itself assumed
        //     to be small.
        // The resulting estimates of the "close" root(s) should be better
        // that the ones we started with.

        auto k0 = 3.0 * c0 * x0 + c1;
        auto k1 = 3.0 * c0 * x0 * x0 + 2 * c1 * x0 + c2;
        auto k2 = c0 * x0 * x0 * x0;
        brick::common::ComplexNumber<Type> epsilon0;
        brick::common::ComplexNumber<Type> epsilon1;
        solveQuadratic(k1/k0, k2/k0, epsilon0, epsilon1);
        auto candidate00 = epsilon0 + x0;
        auto candidate01 = epsilon1 + x0;
        std::cout << "Candidate00: (" << epsilon0 << ") "
                  << candidate00 << " -- "
                  << cubify(c0, c1, c2, c3, candidate00.getRealPart())
                  << std::endl;
        std::cout << "Candidate01: (" << epsilon1 << ") "
                  << candidate01 << " -- "
                  << cubify(c0, c1, c2, c3, candidate01.getRealPart())
                  << std::endl;
        
      }

      root0 = candidate0;
      root1 = candidate1;
      root2 = candidate1;
    }

  } // namespace numeric

} // namespace brick

#endif /* #ifndef BRICK_NUMERIC_SOLVECUBIC_IMPL_HH */
