/**
***************************************************************************
* @file brick/geometry/utilities3D_impl.hh
*
* Source file defining some 3D geometric utilities for finding
* intersects, etc.
*
* Copyright (C) 2007 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_GEOMETRY_UTILITIES3D_IMPL_HH
#define BRICK_GEOMETRY_UTILITIES3D_IMPL_HH

// This file is included by utilities3D.hh, and should not be directly
// included by user code, so no need to include utilities3D.hh here.
//
// #include <brick/geometry/utilities3D.hh>

#include <sstream>
#include <brick/common/exception.hh>
#include <brick/common/types.hh>
#include <brick/linearAlgebra/linearAlgebra.hh>
#include <brick/numeric/array2D.hh>
#include <brick/numeric/numericTraits.hh>
#include <brick/numeric/mathFunctions.hh>
#include <brick/numeric/utilities.hh>

namespace brick {

  namespace geometry {

    namespace privateCode {

      template <class FloatType>
      brick::numeric::Array2D<FloatType>
      pseudoinverse(brick::numeric::Array2D<FloatType> const& A)
      {
        brick::numeric::Array2D<FloatType> ATranspose = A.transpose();
        brick::numeric::Array2D<FloatType> ATA =
          brick::numeric::matrixMultiply<FloatType>(ATranspose, A);
        if(ATA.rows() != 2 || ATA.columns() != 2) {
          BRICK_THROW(brick::common::NotImplementedException,
                      "geometry::privateCode::psuedoInverse()",
                      "Currently only matrices of rank 2 are supported.");
        }

        // Use cofactor inverse, rather than lapack, because it
        // generalizes to more floating point types.
        FloatType determinant = ATA(0, 0) * ATA(1, 1) - ATA(0, 1) * ATA(1, 0);
        if(brick::numeric::absoluteValue(determinant)
           < brick::numeric::NumericTraits<FloatType>::epsilon()) {

          BRICK_THROW(brick::common::ValueException,
                      "geometry::privateCode::psuedoInverse()",
                      "Matrix is not full rank");
        }

        brick::numeric::Array2D<FloatType> ATAInv(2, 2);
        ATAInv(0, 0) = ATA(1, 1) / determinant;
        ATAInv(0, 1) = -(ATA(0, 1) / determinant);
        ATAInv(1, 0) = -(ATA(1, 0) / determinant);
        ATAInv(1, 1) = ATA(0, 0) / determinant;
        return brick::numeric::matrixMultiply<FloatType>(ATAInv, ATranspose);
      }

    } // namespace privateCode


    template <class Type>
    bool
    checkIntersect(Ray3D<Type> const& ray, Triangle3D<Type> const& triangle)
    {
      Type dummy0;
      brick::numeric::Vector3D<Type> dummy1;
      return checkIntersect(ray, triangle, dummy1, dummy0);
    }


    template <class Type>
    bool
    checkIntersect(Ray3D<Type> const& ray, Triangle3D<Type> const& triangle,
                   brick::numeric::Vector3D<Type>& intersect)
    {
      Type dummy;
      return checkIntersect(ray, triangle, intersect, dummy);
    }


    template <class Type>
    bool
    checkIntersect(Ray3D<Type> const& ray, Triangle3D<Type> const& triangle,
                   Type& lambda)
    {
      brick::numeric::Vector3D<Type> dummy;
      return checkIntersect(ray, triangle, dummy, lambda);
    }


    template <class Type>
    bool
    checkIntersect(Ray3D<Type> const& ray, Triangle3D<Type> const& triangle,
                   brick::numeric::Vector3D<Type>& intersect, Type& lambda)
    {
      // The point at which the ray intersects the plane of the
      // triangle satisfies this equation:
      //
      //   v_0 + alpha_0 * e_0 + alpha_1 * e_1 = o + beta * d,
      //
      // where:
      //
      //    v_0, v_1, and v_2 are the three vertices of the triangle.
      //
      //    e_0 = v_1 - v_0 is the first leg of the triangle.
      //
      //    e_1 = v_2 - v_0 is the second leg of the triangle.
      //
      //    o and d are the origin and direction of the ray, respectively.
      //
      //    alpha_0, alpha_1, and beta are scalars.
      //
      // We rearrange this equation and solve for alpha_0, alpha_1, and beta.
      brick::common::Float64 bufferA[9];
      brick::common::Float64 bufferB[3];

      brick::numeric::Array2D<brick::common::Float64> AMatrix(3, 3, bufferA);
      brick::numeric::Array1D<brick::common::Float64> bVector(3, bufferB);

      brick::numeric::Vector3D<Type> e_0 =
        triangle.getVertex1() - triangle.getVertex0();
      brick::numeric::Vector3D<Type> e_1 =
        triangle.getVertex2() - triangle.getVertex0();

      AMatrix[0] = e_0.x();
      AMatrix[3] = e_0.y();
      AMatrix[6] = e_0.z();
      AMatrix[1] = e_1.x();
      AMatrix[4] = e_1.y();
      AMatrix[7] = e_1.z();
      AMatrix[2] = -ray.getDirectionVector().x();
      AMatrix[5] = -ray.getDirectionVector().y();
      AMatrix[8] = -ray.getDirectionVector().z();

      bVector[0] = ray.getOrigin().x() - triangle.getVertex0().x();
      bVector[1] = ray.getOrigin().y() - triangle.getVertex0().y();
      bVector[2] = ray.getOrigin().z() - triangle.getVertex0().z();

      try {
        brick::linearAlgebra::linearSolveInPlace(AMatrix, bVector);
      } catch(brick::common::ValueException const&) {
        // Singular matrix indicates degenerate triangle, or ray
        // parallel to the plane of the triangle.
        return false;
      }

      if(bVector[0] < 0.0 || bVector[1] < 0.0) {
        // Negative values for alpha_0 or alpha_1 indicate that the
        // intersection point is outside the triangle.
        return false;
      }

      // At this point we know that the intersection lies between
      // sides e_0 and e_1.  We need to make sure that this point also
      // lies within the bounds of the third side.  We do this by
      // requiring that the cross products with two adjacent sides be
      // antiparallel.
      brick::numeric::Vector3D<Type> intersectionPoint =
        ray.getOrigin() + bVector[2] * ray.getDirectionVector();
      brick::numeric::Vector3D<Type> offsetFromVertex1 =
        triangle.getVertex1() - intersectionPoint;
      brick::numeric::Vector3D<Type> e_2 =
        triangle.getVertex2() - triangle.getVertex1();

      brick::numeric::Vector3D<Type> crossProduct0 =
        brick::numeric::cross(offsetFromVertex1, -e_0);
      brick::numeric::Vector3D<Type> crossProduct1 =
        brick::numeric::cross(offsetFromVertex1, e_2);
      Type dotProduct =
        brick::numeric::dot<Type>(crossProduct0, crossProduct1);

      if(dotProduct >= 0.0) {
        return false;
      }

      lambda = bVector[2];
      intersect = intersectionPoint;
      return true;
    }


    template <class Type>
    brick::numeric::Vector3D<Type>
    findIntersect(Ray3D<Type> const& ray, Plane3D<Type> const& plane)
    {
      Type dummy;
      return findIntersect(ray, plane, dummy);
    }


    template <class Type>
    brick::numeric::Vector3D<Type>
    findIntersect(Ray3D<Type> const& ray, Plane3D<Type> const& plane, Type& distance)
    {
      brick::common::Float64 bufferA[9];
      brick::common::Float64 bufferB[3];

      brick::numeric::Array2D<Type> AMatrix(3, 3, bufferA);
      brick::numeric::Array1D<Type> bVector(3, bufferB);

      AMatrix[0] = ray.getDirectionVector().x();
      AMatrix[3] = ray.getDirectionVector().y();
      AMatrix[6] = ray.getDirectionVector().z();
      AMatrix[1] = plane.getDirectionVector0().x();
      AMatrix[4] = plane.getDirectionVector0().y();
      AMatrix[7] = plane.getDirectionVector0().z();
      AMatrix[2] = plane.getDirectionVector1().x();
      AMatrix[5] = plane.getDirectionVector1().y();
      AMatrix[8] = plane.getDirectionVector1().z();

      bVector[0] = plane.getOrigin().x() - ray.getOrigin().x();
      bVector[1] = plane.getOrigin().y() - ray.getOrigin().y();
      bVector[2] = plane.getOrigin().z() - ray.getOrigin().z();

      try {
        brick::linearAlgebra::linearSolveInPlace(AMatrix, bVector);
      } catch(brick::common::ValueException const&) {
        std::ostringstream message;
        message << "Unable to find intersection of " << ray << " with "
                << plane << ".  Perhaps the ray is parallel to the plane.";
        BRICK_THROW(brick::common::ValueException, "findIntersect()", message.str().c_str());
      }

      distance = bVector[0];
      return ray.getOrigin() + distance * ray.getDirectionVector();
    }


    template <class Type>
    brick::numeric::Vector3D<Type>
    findIntersect(Ray3D<Type> const& ray0, Ray3D<Type> const& ray1,
                  Type& distance0, Type& distance1, Type& residual)
    {
      // Let z_0 and v_0 be the origin and direction vector of ray0,
      // respectively.  Let z_1 and v_1 be the origin and direction
      // vector of ray1.  Every point, p0, on ray0 can be expressed as:
      //
      // @verbatim
      //   p0 = z_0 + a_0 * v_0
      // @endverbatim
      //
      // Where a_0 is a scalar.  We can write points on ray1
      // similarly.  We're looking for the scalars a_0 and a_1 that
      // minimize (in the least squares sense) `the length of the
      // residual vector (p0 - p1).
      //
      // @verbatim
      //   p0 - p1 = (z_0 + a_0 * v_0) - (z_1 + a_1 * v_1) ~= 0
      // @endverbatim
      //
      // This is easily rearranged into a homogeneous linear system of
      // equations:
      //
      // @verbatim
      //   [v_0; -v_1] * [a_0; a_1]^T = [z_1 - z_0]
      // @endverbatim

      brick::numeric::Array2D<Type> AMatrix(3, 2);
      AMatrix(0, 0) = ray0.getDirectionVector().x();
      AMatrix(1, 0) = ray0.getDirectionVector().y();
      AMatrix(2, 0) = ray0.getDirectionVector().z();
      AMatrix(0, 1) = -ray1.getDirectionVector().x();
      AMatrix(1, 1) = -ray1.getDirectionVector().y();
      AMatrix(2, 1) = -ray1.getDirectionVector().z();

      brick::numeric::Array1D<Type> bVector(3);
      bVector[0] = ray1.getOrigin().x() - ray0.getOrigin().x();
      bVector[1] = ray1.getOrigin().y() - ray0.getOrigin().y();
      bVector[2] = ray1.getOrigin().z() - ray0.getOrigin().z();

      brick::numeric::Array2D<Type> APinv;
      try {
        APinv = privateCode::pseudoinverse(AMatrix);
      } catch(brick::common::ValueException ) {
        BRICK_THROW(brick::common::ValueException, "findIntersect()",
                  "Trouble inverting matrix.  "
                  "Perhaps input rays are parallel.");
      }
      brick::numeric::Array1D<Type> parameters =
        brick::numeric::matrixMultiply<Type>(APinv, bVector);

      distance0 = parameters[0];
      distance1 = parameters[1];
      brick::numeric::Vector3D<Type> point0 =
        ray0.getOrigin() + distance0 * ray0.getDirectionVector();
      brick::numeric::Vector3D<Type> point1 =
        ray1.getOrigin() + distance1 * ray1.getDirectionVector();
      residual = brick::numeric::magnitude<Type>(point1 - point0) / 2.0;
      return Type(0.5) * (point0 + point1);
    }


    template <class Type>
    Circle3D<Type>
    operator*(brick::numeric::Transform3D<Type> const& transform,
              Circle3D<Type> const& inputCircle)
    {
      return Circle3D<Type>(
        transform * inputCircle.getOrigin(),
        transform.rotate(inputCircle.getBasisVector(0)),
        transform.rotate(inputCircle.getBasisVector(1)),
        true);
    }


    template <class Type>
    Plane3D<Type>
    operator*(brick::numeric::Transform3D<Type> const& transform,
              Plane3D<Type> const& inputPlane)
    {
      brick::numeric::Vector3D<Type> newOrigin =
        transform * inputPlane.getOrigin();
      brick::numeric::Vector3D<Type> newEndPoint0 =
        transform * (inputPlane.getOrigin() + inputPlane.getDirectionVector0());
      brick::numeric::Vector3D<Type> newEndPoint1 =
        transform * (inputPlane.getOrigin() + inputPlane.getDirectionVector1());
      return Plane3D<Type>(newOrigin, newEndPoint0, newEndPoint1);
    }


    template <class Type>
    Ray3D<Type>
    operator*(brick::numeric::Transform3D<Type> const& transform,
              Ray3D<Type> const& inputRay)
    {
      brick::numeric::Vector3D<Type> newOrigin =
        transform * inputRay.getOrigin();
      brick::numeric::Vector3D<Type> newEndpoint =
        transform * (inputRay.getOrigin() + inputRay.getDirectionVector());
      return Ray3D<Type>(newOrigin, newEndpoint - newOrigin, false);
    }


    template <class Type>
    Triangle3D<Type>
    operator*(brick::numeric::Transform3D<Type> const& transform,
              Triangle3D<Type> const& inputTriangle)
    {
      brick::numeric::Vector3D<Type> newVertex0 =
        transform * inputTriangle.getVertex0();
      brick::numeric::Vector3D<Type> newVertex1 =
        transform * inputTriangle.getVertex1();
      brick::numeric::Vector3D<Type> newVertex2 =
        transform * inputTriangle.getVertex2();
      return Triangle3D<Type>(newVertex0, newVertex1, newVertex2);
    }


  } // namespace utilities

} // namespace brick

#endif /* #ifndef BRICK_GEOMETRY_UTILITIES3D_IMPL_HH */
