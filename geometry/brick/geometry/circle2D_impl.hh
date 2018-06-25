/**
***************************************************************************
* @file brick/geometry/circle2D_impl.hh
*
* Source file defining the Circle2D class template.
*
* Copyright (C) 2008-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_GEOMETRY_CIRCLE2D_IMPL_HH
#define BRICK_GEOMETRY_CIRCLE2D_IMPL_HH

// This file is included by circle2D.hh, and should not be directly included
// by user code, so no need to include circle2D.hh here.
//
// #include <brick/geometry/circle2D.hh>

namespace brick {

  namespace geometry {

    // The default constructor initializes to the unit circle.
    template <class Type>
    Circle2D<Type>::
    Circle2D()
      : m_origin(0.0, 0.0),
        m_radius(1.0)
    {
      // Empty.
    }


    // This constructor initializes the circle using explicitly
    // specified values.
    template <class Type>
    Circle2D<Type>::
    Circle2D(brick::numeric::Vector2D<Type> const& origin,
             Type const& radius)
      : m_origin(origin),
        m_radius(radius)
    {
      if(radius < 0.0) {
        m_radius = -m_radius;
      }
    }


    // The copy constructor deep copies its argument.
    template <class Type>
    Circle2D<Type>::
    Circle2D(Circle2D<Type> const& source)
      : m_origin(source.m_origin),
        m_radius(source.m_radius)
    {
      // Empty.
    }


    // The assignment operator deep copies its argument.
    template <class Type>
    Circle2D<Type>&
    Circle2D<Type>::
    operator=(Circle2D<Type> const& source)
    {
      if(&source != this) {
        m_origin = source.m_origin;
        m_radius = source.m_radius;
      }
      return *this;
    }


    /* ======= Non-member functions. ======= */

    template <class Type>
    std::ostream&
    operator<<(std::ostream& stream, const Circle2D<Type>& circle)
    {
      stream << "Circle2D{ "
             << circle.getOrigin() << ", "
             << circle.getRadius() << " }";
      return stream;
    }

  } // namespace geometry

} // namespace brick

#endif /* #ifndef BRICK_GEOMETRY_CIRCLE2D_IMPL_HH */
