/**
***************************************************************************
* @file brick/numeric/solveCubic_impl.hh
*
* Implementation of functions for solving cubic polynomial equations
* of a single variable.
*
* Copyright (C) 2001-2009,2012 David LaRose, dlr@cs.cmu.edu
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

  } // namespace numeric

} // namespace brick

#endif /* #ifndef BRICK_NUMERIC_SOLVECUBIC_IMPL_HH */
