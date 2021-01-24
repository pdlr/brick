/**
***************************************************************************
* @file brick/geometry/utilities3D.hh
*
* Header file declaring some 3D geometric utilities for finding
* intersects, etc.
*
* WARNING: Some of the functions in this file are defined only for
* Float64 type, and will remain so until the brickLinearAlgebra
* library is extended to handle types besides Float64.
*
* Copyright (C) 2007 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_GEOMETRY_UTILITIES3D_HH
#define BRICK_GEOMETRY_UTILITIES3D_HH

#include <brick/numeric/vector3D.hh>
#include <brick/numeric/transform3D.hh>
#include <brick/geometry/circle3D.hh>
#include <brick/geometry/plane3D.hh>
#include <brick/geometry/ray3D.hh>
#include <brick/geometry/triangle3D.hh>

namespace brick {

  namespace geometry {

    template <class Type>
    bool
    checkIntersect(Ray3D<Type> const& ray, Triangle3D<Type> const& triangle);


    template <class Type>
    bool
    checkIntersect(Ray3D<Type> const& ray, Triangle3D<Type> const& triangle,
                   brick::numeric::Vector3D<Type>& intersect);


    template <class Type>
    bool
    checkIntersect(Ray3D<Type> const& ray, Triangle3D<Type> const& triangle,
                   Type& lambda);


    template <class Type>
    bool
    checkIntersect(Ray3D<Type> const& ray, Triangle3D<Type> const& triangle,
                   brick::numeric::Vector3D<Type>& intersect, Type& lambda);


    template <class Type>
    brick::numeric::Vector3D<Type>
    findIntersect(Ray3D<Type> const& ray, Plane3D<Type> const& plane);


    template <class Type>
    brick::numeric::Vector3D<Type>
    findIntersect(Ray3D<Type> const& ray, Plane3D<Type> const& plane,
                  Type& distance);


    template <class Type>
    brick::numeric::Vector3D<Type>
    findIntersect(Ray3D<Type> const& ray0, Ray3D<Type> const& ray1,
                  Type& distance0, Type& distance1, Type& residual);


    template <class Type>
    Circle3D<Type>
    operator*(brick::numeric::Transform3D<Type> const& transform,
              Circle3D<Type> const& inputCircle);


    template <class Type>
    Plane3D<Type>
    operator*(brick::numeric::Transform3D<Type> const& transform,
              Plane3D<Type> const& inputPlane);


    template <class Type>
    Ray3D<Type>
    operator*(brick::numeric::Transform3D<Type> const& transform,
              Ray3D<Type> const& inputRay);


    template <class Type>
    Triangle3D<Type>
    operator*(brick::numeric::Transform3D<Type> const& transform,
              Triangle3D<Type> const& inputTriangle);

  } // namespace utilities

} // namespace brick


// Include definitions of inline and template functions.
#include <brick/geometry/utilities3D_impl.hh>

#endif /* #ifndef BRICK_GEOMETRY_UTILITIES3D_HH */
