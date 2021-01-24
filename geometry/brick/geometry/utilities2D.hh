/**
***************************************************************************
* @file brick/geometry/utilities2D.hh
*
* Header file declaring some 2D geometric utilities for finding
* intersects, etc.
*
* Copyright (C) 2009 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_GEOMETRY_UTILITIES2D_HH
#define BRICK_GEOMETRY_UTILITIES2D_HH

#include <brick/numeric/transform2D.hh>
#include <brick/numeric/vector2D.hh>
#include <brick/geometry/lineSegment2D.hh>
#include <brick/geometry/ray2D.hh>
#include <brick/geometry/triangle2D.hh>


namespace brick {

  namespace geometry {

    template <class Type>
    bool
    checkIntersect(LineSegment2D<Type> const& lineSegment0,
                   LineSegment2D<Type> const& lineSegment1);


    template <class Type>
    bool
    checkIntersect(LineSegment2D<Type> const& lineSegment0,
                   LineSegment2D<Type> const& lineSegment1,
                   brick::numeric::Vector2D<Type>& intersect);


    template <class Type>
    bool
    checkIntersect(Ray2D<Type> const& ray,
                   LineSegment2D<Type> const& lineSegment);


    template <class Type>
    bool
    checkIntersect(Ray2D<Type> const& ray,
                   LineSegment2D<Type> const& lineSegment,
                   brick::numeric::Vector2D<Type>& intersect);


    template <class Type>
    bool
    checkIntersect(Ray2D<Type> const& ray,
                   LineSegment2D<Type> const& lineSegment,
                   double& lambda);


    template <class Type>
    bool
    checkIntersect(Ray2D<Type> const& ray,
                   LineSegment2D<Type> const& lineSegment,
                   brick::numeric::Vector2D<Type>& intersect,
                   double& lambda);


    /**
     * Return the centroid of triangle, which is coincident with the
     * intersection of its three medians.
     *
     * @param triangle The triangle of which the centroid will be
     * computed.
     *
     * @return The return value is the centroid.
     */
    template <class Type>
    brick::numeric::Vector2D<Type>
    getCentroid(Triangle2D<Type> const& triangle);


    template <class Type>
    brick::numeric::Vector2D<Type>
    findClosestPoint(brick::numeric::Vector2D<Type> const& point,
                     Ray2D<Type> const& ray);


    template <class Type>
    LineSegment2D<Type>
    operator*(brick::numeric::Transform2D<Type> const& transform,
              LineSegment2D<Type> const& inputSegment);


    template <class Type>
    Ray2D<Type>
    operator*(brick::numeric::Transform2D<Type> const& transform,
              Ray2D<Type> const& inputRay);


  } // namespace utilities

} // namespace brick

// Include definitions of inline and template functions.
#include <brick/geometry/utilities2D_impl.hh>

#endif /* #ifndef BRICK_GEOMETRY_UTILITIES2D_HH */
