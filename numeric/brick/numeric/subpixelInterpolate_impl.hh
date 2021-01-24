/**
***************************************************************************
* @file brick/computerVision/subpixelInterpolate_impl.hh
*
* Header file defining inline and template routines for fitting
* quadratics to local patches of images or arrays, and using those
* quadratics to find local maxima or minima with subpixel precision.
*
* Copyright (C) 2009,2012 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_NUMERIC_SUBPIXELINTERPOLATE_IMPL_HH
#define BRICK_NUMERIC_SUBPIXELINTERPOLATE_IMPL_HH

// This file is included by subPixelInterpolate.hh, and should not be
// directly included by user code, so no need to include
// subPixelInterpolate.hh here.
//
// #include <brick/numeric/subPixelInterpolate.hh>

#include <cmath>

namespace brick {

  namespace numeric {

    // This function solves for the coefficients of a 2D quadratic,
    // given the function values on a 3x3 grid around the origin.
    template <class Type0, class Type1>
    void
    getQuadraticCoefficients3x3(
      Type0 value00, Type0 value01, Type0 value02,
      Type0 value10, Type0 value11, Type0 value12,
      Type0 value20, Type0 value21, Type0 value22,
      Type1& k0, Type1& k1, Type1& k2, Type1& k3, Type1& k4, Type1& k5)
    {
      // Given f(x0, y0) = value00, we can rearrange the quadratic
      // equation like this:
      //
      // @code
      //   [x0*x0, 2*x0*y0, y0*y0, x0, y0, 1] [k0, k1, k2, k3, k5, k5]^T = value00
      // @endcode
      //
      // We're always going to look at a 3x3 neighborhood around (0,
      // 0), so this gives us 9 simultaneous equations that are
      // quadratic in x, y, but linear in the parameters for which
      // we're solving.
      //
      // @code
      //   [     ][k0]   [value00]
      //   [  A  ][k1] = [value01]
      //   [     ][k2]   [value02]
      //          [k3]   [value10]
      //          [k4]   [value11]
      //          [k5]   [value12]
      //                 [value20]
      //                 [value21]
      //                 [value22]
      // @endcode
      //
      // Where A is a 9x6 matrix:
      //
      // @code
      //       [1,  2,  1, -1, -1,  1]
      //       [0,  0,  1,  0, -1,  1]
      //       [1, -2,  1,  1, -1,  1]
      //       [1,  0,  0, -1,  0,  1]
      //   A = [0,  0,  0,  0,  0,  1]
      //       [1,  0,  0,  1,  0,  1]
      //       [1, -2,  1, -1,  1,  1]
      //       [0,  0,  1,  0,  1,  1]
      //       [1,  2,  1,  1,  1,  1]
      // @endcode
      //
      // We simply hardcode the pseudoInverse here.
      k0 = (0.16666666666666663 * value00
            + -0.33333333333333331 * value01
            + 0.16666666666666663 * value02
            + 0.16666666666666663 * value10
            + -0.33333333333333337 * value11
            + 0.16666666666666663 * value12
            + 0.16666666666666663 * value20
            + -0.33333333333333331 * value21
            + 0.16666666666666663 * value22);
      k1 = (0.125 * value00
            + -0.125 * value02
            + -0.125 * value20
            + 0.125 * value22);
      k2 = (0.16666666666666663 * value00
            + 0.16666666666666663 * value01
            + 0.16666666666666663 * value02
            + -0.33333333333333331 * value10
            + -0.33333333333333331 * value11
            + -0.33333333333333331 * value12
            + 0.16666666666666663 * value20
            + 0.16666666666666663 * value21
            + 0.16666666666666663 * value22);
      k3 = (-0.16666666666666666 * value00
            + 0.16666666666666666 * value02
            + -0.16666666666666666 * value10
            + 0.16666666666666666 * value12
            + -0.16666666666666666 * value20
            + 0.16666666666666666 * value22);
      k4 = (-0.16666666666666666 * value00
            + -0.16666666666666666 * value01
            + -0.16666666666666666 * value02
            + 0.16666666666666666 * value20
            + 0.16666666666666666 * value21
            + 0.16666666666666666 * value22);
      k5 = (-0.11111111111111116 * value00
            + 0.22222222222222227 * value01
            + -0.11111111111111116 * value02
            + 0.22222222222222221 * value10
            + 0.55555555555555558 * value11
            + 0.22222222222222221 * value12
            + -0.11111111111111116 * value20
            + 0.22222222222222227 * value21
            + -0.11111111111111116 * value22);
    }


    template <class Type0, class Type1>
    bool
    subpixelInterpolate(Type1 centerRowCoord, Type1 centerColumnCoord,
                        Type0 value00, Type0 value01, Type0 value02,
                        Type0 value10, Type0 value11, Type0 value12,
                        Type0 value20, Type0 value21, Type0 value22,
                        Type1& extremumRowCoord, Type1& extremumColumnCoord,
                        Type1& extremeValue)
    {
      // First, we solve for the coefficients of a quadratic in the
      // neighborhood of (centerRowCoord, centerColumnCoord), with the
      // parameters of the quadratic being location (rowCoord,
      // columnCoord), and the scalar quadratic values being provided
      // in the argument list (value00, value01, etc.) above.
      //
      // @code
      //   valueRC = [R, C] [k0, k1] [R] + [R, C] [k3] + k5
      //                    [k1, k2] [C]          [k4]
      // @endcode
      //
      // This function call solves for k0 through k5 in the above
      // equation.
      Type1 k0, k1, k2, k3, k4, k5;
      getQuadraticCoefficients3x3(
        value00, value01, value02,
        value10, value11, value12,
        value20, value21, value22,
        k0, k1, k2, k3, k4, k5);

      // Now we have the parameters of the quadratic, we need to find
      // its extremum.  That is, the place where it's derivitive goes
      // to zero.  We need to solve this equation for x and y:
      //
      // @code
      //   2[k0, k1][x] + [k3] = [0]
      //    [k1, k2][y]   [k4]   [0]
      // @endcode
      //
      // We simply hardcode a cofactor method inverse below, starting
      // with finding the determinant.
      Type1 determinant = k0 * k2 - k1 * k1;

      // Before we continue with the matrix inversion, consider
      // whether we're going to find a min, a max, or a saddle point.
      // Essentially, we need to know if matrix [k0,k1; k1, k2] has
      // one positive eigenvalue and one negative eigenvalue (in which
      // case we'll find a saddle point).  If we were to solve for
      // those eigenvalues, we'd be trying to find values of lambda
      // such that
      //
      // @code
      //   determinant([k0 - lambda, k1         ]) == 0
      //               [         k1, k2 - lambda]
      // @endcode
      //
      // If you write this out, you get a quadratic equation:
      //
      //   lambda^2 + (k0 + k2) * lambda + (k0 * k2 - k1^2) = 0
      //
      // the signs of the two solutions to this equation are different
      // if and only if (k0 * k2 - k1^2) is negative (check out the
      // documentation of solveQuadratic() in utilities.h to see why
      // this is true), but this is just the determinant we just computed.
      //
      // So, if determinant is zero, we can't do the matrix inverse.
      // If determinant < 0, we'll only find a saddle point.  In both
      // cases, we simply return false.
      if(determinant <= 0) {
        return false;
      }

      // Now continue with the matrix inversion, and compute the
      // position of the extremum.
      Type1 detTimes2 = determinant * 2.0;
      Type1 inverseMx00 = k2 / detTimes2;
      Type1 inverseMx01 = -k1 / detTimes2;
      Type1 inverseMx11 = k0 / detTimes2;

      Type1 maxX = -k3 * inverseMx00 - k4 * inverseMx01;
      Type1 maxY = -k3 * inverseMx01 - k4 * inverseMx11;

      extremumRowCoord = centerRowCoord + maxY;
      extremumColumnCoord = centerColumnCoord + maxX;

      // Now that we have the location of the extremum, use the
      // quadratic equation to estimate its "real" value.
      extremeValue =
        (k0 * maxX * maxX + 2 * k1 * maxX * maxY + k2 * maxY * maxY
         + k3 * maxX + k4 * maxY + k5);

      return true;
    }



    // Given pixel values in a 3x1 arrangement centerPosition, this
    // function fits a quadratic to the values and returns the
    // location and interpolated value of the extremum (max or min
    // value) of the quadratic.
    template <class Type0, class Type1>
    bool
    subpixelInterpolate(Type1 centerPosition,
                        Type0 value0, Type0 value1, Type0 value2,
                        Type1& extremumPosition,
                        Type1& extremeValue)
    {
      // First, we solve for the coefficients of a quadratic in the
      // neighborhood of centerPosition, with the parameters of the
      // quadratic being position, and the scalar quadratic values
      // being provided in the argument list (value0, value1, value2)
      // above.
      //
      // @code
      //   valueP = k0*p*p + k1*p + k2
      // @endcode
      //
      // To make this easier, we solve as if centerPosition were 0.
      // We'll correct this later.  Assuming centerPosition is zero,
      // we can write the above equation three times, one for each of
      // the input values, assuming coodinates of -1, 0, 1.  Writing
      // this in matrix form gives:
      //
      // @code
      //   [1, -1, 1][k0]   [value0]
      //   [0,  0, 1][k1] = [value1]
      //   [1,  1, 1][k2]   [value2]
      // @endcode
      //
      // Solving this using the Moore-Penrose pseudoinverse gives:
      //
      // @code
      //   [k0]   [ 0.5, -1.0, 0.5][value0]
      //   [k1] = [-0.5,  0.0, 0.5][value1]
      //   [k2]   [ 0.0,  1.0, 0.0][value2]
      // @endcode
      Type1 k0 =  0.5 * value0 - value1 + 0.5 * value2;
      Type1 k1 = -0.5 * value0 + 0.5 * value2;
      Type1 k2 = value1;

      // Now we have the parameters of the quadratic, we need to find
      // its extremum.  That is, the place where it's derivitive goes
      // to zero.  We need to solve this equation for p:
      //
      // @code
      //   2*k0*p + k1 = 0
      // @code
      //
      // If k0 is zero, we can't solve this.  Instead we return false.
      if(k0 == 0) {
        return false;
      }
      Type1 maxP = -k1 / (2 * k0);

      // Now correct the "assume centerPosition is zero" dodge we
      // introduced above.
      extremumPosition = centerPosition + maxP;

      // Now that we have the location of the extremum, use the
      // quadratic equation to estimate its "real" value.
      extremeValue = k0 * maxP * maxP + k1 * maxP + k2;

      return true;
    }

  } // namespace numeric

} // namespace brick

#endif /* #ifndef BRICK_NUMERIC_SUBPIXELINTERPOLATE_IMPL_HH */
