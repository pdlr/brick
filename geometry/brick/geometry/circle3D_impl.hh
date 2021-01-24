/**
***************************************************************************
* @file brick/geometry/circle3D_impl.hh
*
* Source file defining the Circle3D class template.
*
* Copyright (C) 2014 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_GEOMETRY_CIRCLE3D_IMPL_HH
#define BRICK_GEOMETRY_CIRCLE3D_IMPL_HH

// This file is included by circle3D.hh, and should not be directly included
// by user code, so no need to include circle3D.hh here.
//
// #include <brick/geometry/circle3D.hh>

#include <brick/common/mathFunctions.hh>
#include <brick/numeric/utilities.hh>

namespace brick {

  namespace geometry {

    // The default constructor initializes to the unit circle in the
    // X-Y plane.
    template <class Type>
    Circle3D<Type>::
    Circle3D()
      : m_origin(0.0, 0.0, 0.0),
        m_radius(1.0),
        m_basisVector0(1.0, 0.0, 0.0),
        m_basisVector1(0.0, 1.0, 0.0)
    {
      // Empty.
    }


    // This constructor initializes the circle using explicitly
    // specified values.
    template <class Type>
    Circle3D<Type>::
    Circle3D(brick::numeric::Vector3D<Type> const& origin,
             brick::numeric::Vector3D<Type> const& basisVector0,
             brick::numeric::Vector3D<Type> const& basisVector1,
             bool isOrthonormalized)
      : m_origin(origin),
        m_radius(numeric::magnitude<Type>(basisVector0)),
        m_basisVector0(basisVector0),
        m_basisVector1(basisVector1)
    {
      if(!isOrthonormalized) {
        // Make sure m_basisVector1 is perpendicular to basisVector0.
        Type dotProduct = numeric::dot<Type>(m_basisVector1, m_basisVector0);
        dotProduct /= (m_radius * m_radius);
        m_basisVector1 -= dotProduct * m_basisVector0;

        // Make sure m_basisVector1 has length equal to m_radius.
        Type length = numeric::magnitude<Type>(m_basisVector1);
        m_basisVector1 *= m_radius / length;
      }
    }


    // The copy constructor deep copies its argument.
    template <class Type>
    Circle3D<Type>::
    Circle3D(Circle3D<Type> const& source)
      : m_origin(source.m_origin),
        m_radius(source.m_radius),
        m_basisVector0(source.m_basisVector0),
        m_basisVector1(source.m_basisVector1)
    {
      // Empty.
    }


    // The assignment operator deep copies its argument.
    template <class Type>
    Circle3D<Type>&
    Circle3D<Type>::
    operator=(Circle3D<Type> const& source)
    {
      if(&source != this) {
        m_origin = source.m_origin;
        m_radius = source.m_radius;
        m_basisVector0 = source.m_basisVector0;
        m_basisVector1 = source.m_basisVector1;
      }
      return *this;
    }


    // This member function returns the geometric center of the circle.
    template <class Type>
    brick::numeric::Vector3D<Type> const&
    Circle3D<Type>::
    getBasisVector(brick::common::UInt32 basisIndex) const
    {
      if(0 == basisIndex) {
        return this->m_basisVector0;
      }
      return this->m_basisVector1;
    }


    // This member function returns a point on the perimiter of the
    // circle, rotated a user-specified angle around the circle,
    // from an arbitrary-but-unchanging start point.
    template <class Type>
    brick::numeric::Vector3D<Type>
    Circle3D<Type>::
    getPerimeterPoint(Type angle) const
    {
      return (m_origin
              + brick::common::cosine(angle) * m_basisVector0
              + brick::common::sine(angle) * m_basisVector1);
    }


    /* ======= Non-member functions. ======= */

    // template <class Type>
    // std::ostream&
    // operator<<(std::ostream& stream, const Circle3D<Type>& circle)
    // {
    //   stream << "Circle3D{ "
    //          << circle.getOrigin() << ", "
    //          << circle.getRadius() << ", "
    //          << circle.getPerimiterPoint(0.0) << ", "
    //          << circle.getPerimiterPoint(common::constants::piOverTwo)
    //          << " }";
    //   return stream;
    // }

  } // namespace geometry

} // namespace brick

#endif /* #ifndef BRICK_GEOMETRY_CIRCLE3D_IMPL_HH */
