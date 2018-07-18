/**
***************************************************************************
* @file brick/computerVision/subpixelInterpolate.hh
*
* Header file declaring routines for fitting quadratics to local
* patches of images or arrays, and using those quadratics to find
* local maxima or minima with subpixel precision.
*
* Copyright (C) 2009,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_NUMERIC_SUBPIXELINTERPOLATE_HH
#define BRICK_NUMERIC_SUBPIXELINTERPOLATE_HH


namespace brick {

  namespace numeric {

    /**
     * This function solves for the coefficients of a 2D quadratic,
     * given the function values on a 3x3 grid around the origin.  It
     * is used internally by subpixelInterpolate(), but is exposed
     * here in case its useful in other contexts.  This function is
     * completely hardcoded, so it should be reasonably fast.
     *
     * Note that this function talks in (x, y) coordinates, while
     * subpixelInterpolate() talks in (row, column) coordinates.  Our
     * normal convention is that row number goes up with increasing y
     * coordinate and column number goes up with increasing x
     * coordinate, so this makes the valueYX arguments match nicely
     * with the valueRowColumn arguments of subpixelInterpolate(), but
     * it's something to be aware of if you don't follow this
     * convention.
     *
     * More detail about what this routine does: assuming you have
     * the function
     *
     * @code
     *   f(x, y) = [x, y] [k0, k1] [x] + [x, y] [k3] + k5
     *                    [k1, k2] [y]          [k4]
     * @endcode
     *
     * and are given the values of f(x, y) for the 3x3 grid x = {-1,
     * 0, 1}, y = {-1, 0, 1}, this function returns the values of k0
     * ... k5 that most nearly satisfy the equation (in the least
     * squares sense).
     *
     * @param value00 This argument is the value of f(-1, -1).
     *
     * @param value01 This argument is the value of f(0, -1).
     *
     * @param value02 This argument is the value of f(1, -1).
     *
     * @param value10 This argument is the value of f(-1, 0).
     *
     * @param value11 This argument is the value of f(0, 0).
     *
     * @param value12 This argument is the value of f(1, 0).
     *
     * @param value20 This argument is the value of f(-1, 1).
     *
     * @param value21 This argument is the value of f(0, 1).
     *
     * @param value22 This argument is the value of f(1, 1).
     *
     * @param k0 This argument returns a recovered parameter of the
     * quadratic function.
     *
     * @param k1 This argument returns a recovered parameter of the
     * quadratic function.
     *
     * @param k2 This argument returns a recovered parameter of the
     * quadratic function.
     *
     * @param k3 This argument returns a recovered parameter of the
     * quadratic function.
     *
     * @param k4 This argument returns a recovered parameter of the
     * quadratic function.
     *
     * @param k5 This argument returns a recovered parameter of the
     * quadratic function.
     */
    template <class Type0, class Type1>
    inline void
    getQuadraticCoefficients3x3(
      Type0 value00, Type0 value01, Type0 value02,
      Type0 value10, Type0 value11, Type0 value12,
      Type0 value20, Type0 value21, Type0 value22,
      Type1& k0, Type1& k1, Type1& k2, Type1& k3, Type1& k4, Type1& k5);


    /**
     * Given pixel values in a 3x3 array around (centerRow,
     * centerColumn), this function fits a quadratic to the array
     * values and returns the location and interpolated value of the
     * extremum (max or min value) of the quadratic.  Note that when
     * choosing the values of parameters centerRowCoord and
     * centerColumnCoord, you have to think carefully about your pixel
     * coordinate system.  It is necessary to decide whether a row
     * coordinate of 6.0 means "the center of the sixth row," or "On
     * the border between the fifth and sixth rows".  We encourage the
     * convention that the pixel at integer coordinates [i, j] has
     * corners at (i, j), (i, j+1), (i+1, j), and (i+1, j+1).  Under
     * this convention, the center of pixel [i, j] is at row
     * coordinate (i + 0.5), column coordinate (j + 0.5).
     *
     * @param centerRowCoord This argument is the row coordinate of the
     * center of the center pixel of the 3x3 neighborhood.
     *
     * @param centerColumnCoord This argument is the column coordinate of the
     * center of the center pixels of the 3x3 neighborhood.
     *
     * @param value00 This argument is the pixel value at (centerRowCoord -
     * 1, centerColumnCoord - 1).
     *
     * @param value01 This argument is the pixel value at (centerRowCoord -
     * 1, centerColumnCoord).
     *
     * @param value02 This argument is the pixel value at (centerRowCoord -
     * 1, centerColumnCoord + 1).
     *
     * @param value10 This argument is the pixel value at (centerRowCoord,
     * centerColumnCoord - 1).
     *
     * @param value11 This argument is the pixel value at (centerRowCoord,
     * centerColumnCoord).
     *
     * @param value12 This argument is the pixel value at (centerRowCoord,
     * centerColumnCoord + 1).
     *
     * @param value20 This argument is the pixel value at (centerRowCoord +
     * 1, centerColumnCoord - 1).
     *
     * @param value21 This argument is the pixel value at (centerRowCoord +
     * 1, centerColumnCoord).
     *
     * @param value22 This argument is the pixel value at (centerRowCoord +
     * 1, centerColumnCoord + 1).
     *
     * @param extremumRowCoord This argument returns the (subpixel) row
     * coordinate of the recovered extremum.
     *
     * @param extremumColumnCoord This argument returns the (subpixel)
     * column coordinate of the recovered extremum.
     *
     * @param extremeValue This argument returns the interpolated
     * pixel value at the recovered extremum.
     *
     * @return The return value is true if the extremum is well
     * defined, false otherwise.  Note that this function does not
     * check to see whether the location is extremum is reasonable.
     * You'll have to check for yourself to make sure that it doesn't
     * think the extremum is miles from (centerRowCoord,
     * centerColumnCoord).  A reasonable check is to make sure it's
     * within the original 3x3 neighborhood.
     */
    template <class Type0, class Type1>
    bool
    subpixelInterpolate(Type1 centerRowCoord, Type1 centerColumnCoord,
                        Type0 value00, Type0 value01, Type0 value02,
                        Type0 value10, Type0 value11, Type0 value12,
                        Type0 value20, Type0 value21, Type0 value22,
                        Type1& extremumRowCoord, Type1& extremumColumnCoord,
                        Type1& extremeValue);


    /**
     * Given pixel values in a 3x1 arrangement around centerPosition,
     * this function fits a quadratic to the values and returns the
     * location and interpolated value of the extremum (max or min
     * value) of the quadratic.  Note that when choosing the values of
     * parameter centerPosition, you have to think carefully about
     * your pixel coordinate system.  It is necessary to decide
     * whether a coordinate of 6.0 means "the center of the sixth
     * pixel," or "On the border between the fifth and sixth pixels".
     * We encourage the convention that the pixel at integer
     * coordinate i has boundaries at i and (i+1).  Under this
     * convention, the center of pixel i is at coordinate (i + 0.5).
     *
     * @param centerPosition This argument is the coordinate of the
     * center of the center pixel of the 3x1 neighborhood.
     *
     * @param value0 This argument is the pixel value at
     * (centerPosition - 1).
     *
     * @param value1 This argument is the pixel value at
     * centerPosition.
     *
     * @param value2 This argument is the pixel value at
     * (centerPosition + 1).
     *
     * @param extremumPosition This argument returns the (subpixel)
     * position of the recovered extremum.
     *
     * @param extremeValue This argument returns the interpolated
     * pixel value at the recovered extremum.
     *
     * @return The return value is true if the extremum is well
     * defined, false otherwise.  Note that this function does not
     * check to see whether the location is extremum is reasonable.
     * You'll have to check for yourself to make sure that it doesn't
     * think the extremum is miles from centerPosition.
     */
    template <class Type0, class Type1>
    bool
    subpixelInterpolate(Type1 centerPosition,
                        Type0 value0, Type0 value1, Type0 value2,
                        Type1& extremumPosition,
                        Type1& extremeValue);

  } // namespace numeric

} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/numeric/subpixelInterpolate_impl.hh>


#endif /* #ifndef BRICK_NUMERIC_SUBPIXELINTERPOLATE_HH */
