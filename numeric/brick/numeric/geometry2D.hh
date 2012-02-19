/**
***************************************************************************
* @file brick/numeric/geometry2D.hh
*
* Header file declaring some 2D geometric routines for dlrNumeric.
*
* Copyright (C) 2007,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_GEOMETRY2D_HH
#define BRICK_NUMERIC_GEOMETRY2D_HH

#include <brick/numeric/vector2D.hh>

namespace brick {

  namespace numeric {
    
    /** 
     * This function finds a 2D point, given two other points and
     * distance between each of those points and the point to be
     * recovered.  That is, it answers the question "what point is at
     * range0 from point0 and also at range1 from point1?"  Note that
     * there may be two points which satisfy this constraint, or there
     * may be none if the sum of range0 and range1 is less than the
     * distance between point0 and point1.
     * 
     * @param point0 This argument specifies the first known point.
     * 
     * @param point1 This argument specifies the second known point.
     * 
     * @param range0 This argument specifies the distance from point0 to
     * the target point.
     * 
     * @param range1 This argument specifies the distance from point1 to
     * the target point.
     * 
     * @param intersection0 This argument returns the first of the two
     * points which is at range0 from point0 and at range1 from point1.
     * If no solution is found, then this argument will not be touched.
     * 
     * @param intersection1 This argument returns the second of the two
     * points which is at range0 from point0 and at range1 from point1.
     * In degenerate cases, intersection1 will be identical to
     * intersection0.  If no solution is found, then this argument will
     * not be touched.
     * 
     * @return The return value is true if a solution was found, false
     * otherwise.
     */
    template <class Type>
    bool
    bilaterate(Vector2D<Type> const& point0, Vector2D<Type> const& point1,
               Type range0, Type range1,
               Vector2D<Type>& intersection0, Vector2D<Type>& intersection1);

  } // namespace numeric

} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/numeric/geometry2D_impl.hh>

#endif /* #ifndef BRICK_NUMERIC_GEOMETRY2D_HH */
