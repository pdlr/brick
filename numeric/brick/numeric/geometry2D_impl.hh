/**
***************************************************************************
* @file brick/numeric/geometry2D.cc
*
* Source file file declaring useful functions for dlrNumeric.
*
* Copyright (C) 2001-2007 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/


#ifndef BRICK_NUMERIC_GEOMETRY2D_IMPL_HH
#define BRICK_NUMERIC_GEOMETRY2D_IMPL_HH

// This file is included by geometry2D.hh, and should not be directly included
// by user code, so no need to include geometry2D.hh here.
// 
// #include <brick/numeric/geometry2D.hh>

#include <brick/numeric/utilities.hh>

namespace brick {

  namespace numeric {

    template <class Type>
    bool
    bilaterate(Vector2D<Type> const& point0, Vector2D<Type> const& point1,
               Type range0, Type range1,
               Vector2D<Type>& intersection0, Vector2D<Type>& intersection1)
    {
      // We use the traditional approach of transforming the points so
      // that point0 is at the origin and point1 is on the X axis,
      // solving this new bilateration, and then transforming back to
      // the original coordinate system.

      Vector2D<Type> baselineVector(point1 - point0);
    
      // Coordinates for the translated input points will be (0, 0), and
      // (x0, 0);
      Type x0 = magnitude<Type>(baselineVector);
    
      // Compute direction vectors for the new X and Y axis.
      Vector2D<Type> e0(baselineVector / x0);
      Vector2D<Type> e1(-e0.y(), e0.x());

      // In our new coordinate system, we can write x^2 + y^2 = r0^2,
      // and (x - x0)^2 + y^2 = r1^2.
      //
      // Subtracting these two equations gives 2*x0*x - x0^2 = r0^2 -
      // r1^2, which gives us x = (r0^2 - r1^2 + x0^2) / (2*x0).
      Type r0Squared = range0 * range0;
      Type xValue = (r0Squared - range1 * range1 + x0 * x0) / (2 * x0);

      // Substituting this into x^2 + y^2 = r0^2, we have y = sqrt(r0^2 - x^2).
      Type r0SquaredMinusXSquared = r0Squared - xValue * xValue;
      if(r0SquaredMinusXSquared < 0.0) {
        return false;
      }
      Type yValue = std::sqrt(r0SquaredMinusXSquared);

      // Now we tranform back into the original coordinate system.
      Vector2D<Type> xComponent = xValue * e0;
      Vector2D<Type> yComponent = yValue * e1;
      intersection0 = point0 + xComponent + yComponent;
      intersection1 = point0 + xComponent - yComponent;
      return true;
    }
  
  } // namespace numeric

} // namespace brick

#endif /* #ifndef BRICK_NUMERIC_GEOMETRY2D_IMPL_HH */
