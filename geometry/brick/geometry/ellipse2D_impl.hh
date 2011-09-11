/**
***************************************************************************
* @file brick/geometry/ellipse2D_impl.hh
*
* Source file defining the Ellipse2D class.
*
* Copyright (C) 2008 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_GEOMETRY_ELLIPSE2D_IMPL_HH
#define BRICK_GEOMETRY_ELLIPSE2D_IMPL_HH

// This file is included by ellipse2D.hh, and should not be directly included
// by user code, so no need to include ellipse2D.hh here.
// 
// #include <brick/geometry/ellipse2D.hh>

namespace brick {

  namespace geometry {
    
    // The default constructor initializes to the unit circle.
    template <class Type>
    Ellipse2D<Type>::
    Ellipse2D()
      : m_origin(0.0, 0.0),
        m_semimajorAxis(1.0, 0.0),
        m_semiminorAxis(0.0, 1.0)
    {
      // Empty.
    }

    
    // This constructor initializes the ellipse using explicitly
    // specified values.
    template <class Type>
    Ellipse2D<Type>::
    Ellipse2D(brick::numeric::Vector2D<Type> const& origin,
              brick::numeric::Vector2D<Type> const& semimajorAxis,
              Type ratio)
      : m_origin(origin),
        m_semimajorAxis(semimajorAxis),
        m_semiminorAxis(semimajorAxis.y() * ratio, semimajorAxis.x() * ratio)
    {
      if(ratio > 1.0) {
        std::swap(m_semimajorAxis, m_semiminorAxis);
      }
    }

    
    // The copy constructor deep copies its argument.
    template <class Type>
    Ellipse2D<Type>::
    Ellipse2D(Ellipse2D<Type> const& source)
      : m_origin(source.m_origin),
        m_semimajorAxis(source.m_semimajorAxis),
        m_semiminorAxis(source.m_semiminorAxis)
    {
      // Empty.
    }


    // The assignment operator deep copies its argument.
    template <class Type>
    Ellipse2D<Type>&
    Ellipse2D<Type>::
    operator=(Ellipse2D<Type> const& source)
    {
      if(&source != this) {
        m_origin = source.m_origin;
        m_semimajorAxis = source.m_semimajorAxis;
        m_semiminorAxis = source.m_semiminorAxis;
      }
      return *this;
    }


    /* ======= Non-member functions. ======= */

    std::ostream&
    operator<<(std::ostream& stream, const Ellipse2D& ellipse)
    {
      stream << "Ellipse2D{ "
             << ellipse.getOrigin() << ", "
             << ellipse.getSemimajorAxis() << ", "
             << ellipse.getSemiminorAxis() << " }";
      return stream;
    }
    
  } // namespace geometry
    
} // namespace brick

#endif /* #ifndef BRICK_GEOMETRY_CIRCLE2D_IMPL_HH */
