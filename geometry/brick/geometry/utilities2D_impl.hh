/**
***************************************************************************
* @file brick/utilities/utilities2D.cc
*
* Source file defining some 2D geometric utilities for finding
* intersects, etc.
*
* Copyright (C) 2009-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_GEOMETRY_UTILITIES2D_IMPL_HH
#define BRICK_GEOMETRY_UTILITIES2D_IMPL_HH

// This file is included by circle2D.hh, and should not be directly included
// by user code, so no need to include circle2D.hh here.
// 
// #include <brick/geometry/utilities2D.hh>
#include <brick/numeric/utilities.hh>

// Anonymous namespace for local symbols.
namespace {

  template <class Type>
  bool
  solve2By2LinearSystem(Type const& a00, Type const& a01,
                        Type const& a10, Type const& a11,
                        Type const& b0, Type const& b1,
                        Type& x0, Type& x1) {
    // Use cofactor inversion.
    Type const& determinant = a00 * a11 - a01 * a10;
    if(determinant == 0.0) {
      return false;
    }
    x0 = (a11 * b0 - a01 * b1) / determinant;
    x1 = (a00 * b1 - a10 * b0) / determinant;
    return true;
  }

} // namespace


namespace brick {

  namespace geometry {

    template <class Type>
    bool
    checkIntersect(LineSegment2D<Type> const& lineSegment0,
                   LineSegment2D<Type> const& lineSegment1)
    {
      brick::numeric::Vector2D<Type> dummy;
      return checkIntersect(lineSegment0, lineSegment1, dummy);
    }


    template <class Type>
    bool
    checkIntersect(LineSegment2D<Type> const& lineSegment0,
                   LineSegment2D<Type> const& lineSegment1,
                   brick::numeric::Vector2D<Type>& intersect)
    {
      // The point at which lineSegment0 intersects the line of
      // lineSegment1 satisfies this equation:
      //
      //    v_0 + alpha * (v_1 - v_0) = w_0 + beta * (w_1 - w_0)
      //
      // where:
      //
      //    v_0 and v_1 are the vertices of the first line segment.
      //
      //    w_0 and w_1 are the vertices of the second line segment.
      //
      //    o and d are the origin and direction of the ray, respectively.
      //
      //    alpha and beta are scalar free parameters.
      //
      // We rearrange this equation:
      //
      //   A * [alpha] = b
      //       [ beta]
      //
      // where A is the 2x2 matrix [(v_0 - v_1), (w_1 - w_0)], and b is the 2
      // element vector [v_0 - w_0].
      
      // First find matrix A.
      brick::numeric::Vector2D<Type> lineDirection0 =
        lineSegment0.getVertex0() - lineSegment0.getVertex1();
      brick::numeric::Vector2D<Type> lineDirection1 =
        lineSegment1.getVertex0() - lineSegment1.getVertex1();
      double const& A_00 = lineDirection0.x();
      double const& A_01 = -lineDirection1.x();
      double const& A_10 = lineDirection0.y();
      double const& A_11 = -lineDirection1.y();

      // Now find b.  We'll call it bb.
      brick::numeric::Vector2D<Type> bb =
        lineSegment0.getVertex0() - lineSegment1.getVertex0();

      // Solve for alpha and beta.
      double alpha;
      double beta;
      if(!solve2By2LinearSystem(A_00, A_01, A_10, A_11, bb.x(), bb.y(),
                                alpha, beta)) {
        return false;
      }
      
      // The line segment runs from alpha = 0 to alpha = 1 and from
      // beta = 0 to beta = 1.  Check this here.
      if((alpha < 0.0) || (alpha >= 1.0) || (beta < 0.0) || (beta >= 1.0)) {
        return false;
      }

      // All done.  Now fill in the return parameters.
      intersect = lineSegment0.getVertex0() - alpha * lineDirection0;
      return true;
    }

    
    template <class Type>
    bool
    checkIntersect(Ray2D<Type> const& ray,
                   LineSegment2D<Type> const& lineSegment)
    {
      double dummy0;
      brick::numeric::Vector2D<Type> dummy1;
      return checkIntersect(ray, lineSegment, dummy1, dummy0);
    }


    template <class Type>
    bool
    checkIntersect(Ray2D<Type> const& ray,
                   LineSegment2D<Type> const& lineSegment,
                   brick::numeric::Vector2D<Type>& intersect)
    {
      double dummy;
      return checkIntersect(ray, lineSegment, intersect, dummy);
    }
    

    template <class Type>
    bool
    checkIntersect(Ray2D<Type> const& ray,
                   LineSegment2D<Type> const& lineSegment,
                   double& lambda)
    {
      brick::numeric::Vector2D<Type> dummy;
      return checkIntersect(ray, lineSegment, dummy, lambda);
    }


    template <class Type>
    bool
    checkIntersect(Ray2D<Type> const& ray,
                   LineSegment2D<Type> const& lineSegment,
                   brick::numeric::Vector2D<Type>& intersect, double& lambda)
    {
      // The point at which the ray intersects the line of the
      // line segment satisfies this equation:
      //
      //    v_0 + alpha * (v_1 - v_0) = o + beta * d
      //
      // where:
      //
      //    v_0 and v_1 are the vertices of the line segment.
      //
      //    o and d are the origin and direction of the ray, respectively.
      //
      //    alpha and beta are scalar free parameters.
      //
      // We rearrange this equation:
      //
      //   A * [alpha] = b
      //       [ beta]
      //
      // where A is the 2x2 matrix [(v_0 - v_1), d], and b is the 2
      // element vector [v_0 - o].
      
      // First find matrix A.
      brick::numeric::Vector2D<Type> lineDirection =
        lineSegment.getVertex0() - lineSegment.getVertex1();
      brick::numeric::Vector2D<Type> const& rayDirection =
        ray.getDirectionVector();
      double const& A_00 = lineDirection.x();
      double const& A_01 = rayDirection.x();
      double const& A_10 = lineDirection.y();
      double const& A_11 = rayDirection.y();

      // Now find b.  We'll call it bb.
      brick::numeric::Vector2D<Type> bb =
        lineSegment.getVertex0() - ray.getOrigin();

      // Solve for alpha and beta.
      double alpha;
      double beta;
      if(!solve2By2LinearSystem(A_00, A_01, A_10, A_11, bb.x(), bb.y(),
                                alpha, beta)) {
        return false;
      }

      // The line segment runs from alpha = 0 to alpha = 1.  Check
      // this here.
      if((alpha < 0.0) || (alpha >= 1.0)) {
        return false;
      }

      // Check that the intersection is with the ray, not with the
      // fictional part of the ray that extends back past the ray
      // origin.
      if(beta < 0.0) {
        return false;
      }
      
      // All done.  Now fill in the return parameters.
      lambda = beta;
      intersect = ray.getOrigin() + beta * ray.getDirectionVector();
      return true;
    }

    
    // Return the centroid of triangle, which is coincident with the
    // intersection of its three medians.
    template <class Type>
    brick::numeric::Vector2D<Type>
    getCentroid(Triangle2D<Type> const& triangle)
    {
      // First find the three medians.
      LineSegment2D<Type> median0(
        triangle.getVertex0(),
        (triangle.getVertex1() + triangle.getVertex2()) / Type(2.0));
      LineSegment2D<Type> median1(
        triangle.getVertex1(),
        (triangle.getVertex2() + triangle.getVertex1()) / Type(2.0));
      LineSegment2D<Type> median2(
        triangle.getVertex2(),
        (triangle.getVertex0() + triangle.getVertex0()) / Type(2.0));

      // Find three intersections.
      numeric::Vector2D<Type> intersect0;
      numeric::Vector2D<Type> intersect1;
      numeric::Vector2D<Type> intersect2;
      checkIntersect(median0, median1, intersect0);
      checkIntersect(median1, median2, intersect1);
      checkIntersect(median2, median0, intersect2);

      // Return the average.
      return (intersect0 + intersect1 + intersect2) / Type(3.0);
    }

    
    template <class Type>
    brick::numeric::Vector2D<Type>
    findClosestPoint(brick::numeric::Vector2D<Type> const& point,
                     Ray2D<Type> const& ray)
    {
      brick::numeric::Vector2D<Type> v1 = point - ray.getOrigin();
      double kk = brick::numeric::dot<Type>(v1, ray.getDirectionVector());
      return ray.getOrigin() + kk * ray.getDirectionVector();
    }


    template <class Type>
    LineSegment2D<Type>
    operator*(brick::numeric::Transform2D<Type> const& transform,
              LineSegment2D<Type> const& inputSegment)
    {
      brick::numeric::Vector2D<Type> newVertex0 =
        transform * inputSegment.getVertex0();
      brick::numeric::Vector2D<Type> newVertex1 =
        transform * inputSegment.getVertex1();
      return LineSegment2D<Type>(newVertex0, newVertex1);
    }
    
    
    template <class Type>
    Ray2D<Type>
    operator*(brick::numeric::Transform2D<Type> const& transform,
              Ray2D<Type> const& inputRay)
    {
      brick::numeric::Vector2D<Type> newOrigin =
        transform * inputRay.getOrigin();
      brick::numeric::Vector2D<Type> newEndpoint =
        transform * (inputRay.getOrigin() + inputRay.getDirectionVector());
      return Ray2D<Type>(newOrigin, newEndpoint - newOrigin, false);
    }
    
    
  } // namespace utilities
    
} // namespace brick

#endif /* #ifndef BRICK_GEOMETRY_UTILITIES2D_IMPL_HH */
