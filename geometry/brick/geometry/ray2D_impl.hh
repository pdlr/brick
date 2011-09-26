/**
***************************************************************************
* @file brick/geometry/ray2D_impl.hh
*
* Source file defining the Ray2D class template.
*
* Copyright (C) 2007-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_GEOMETRY_RAY2D_IMPL_HH
#define BRICK_GEOMETRY_RAY2D_IMPL_HH

// This file is included by ray2D.hh, and should not be directly included
// by user code, so no need to include ray2D.hh here.
// 
// #include <brick/geometry/ray2D.hh>
#include <brick/numeric/utilities.hh>

namespace brick {

  namespace geometry {
    
    // The default constructor initializes to the ray that starts
    // at the origin and points along the X axis.
    template <class Type>
    Ray2D<Type>::
    Ray2D()
      : m_origin(0.0, 0.0), m_direction(1.0, 0.0)
    {
      // Empty.
    }

    
    // This constructor initializes the ray using a point and a
    // direction vector.
    template <class Type>
    Ray2D<Type>::
    Ray2D(brick::numeric::Vector2D<Type> const& point,
          brick::numeric::Vector2D<Type> const& direction,
          bool normalize)
      : m_origin(point),
        m_direction(direction)
    {
      if(normalize) {
        m_direction /= brick::numeric::magnitude<Type>(m_direction);
      }
    }

    
    // This constructor initializes the ray according to the
    // equation Ax + By + C = 0.
    template <class Type>
    Ray2D<Type>::
    Ray2D(Type const& aa, Type const& bb, Type const& cc)
      : m_origin(-(aa * cc), -(bb * cc)),
        m_direction(-bb, aa)
    {
      Type kk = aa * aa + bb * bb;
      m_origin /= kk;
      m_direction /= brick::numeric::magnitude<Type>(m_direction);
    }

    
    // The copy constructor deep copies its argument.
    template <class Type>
    Ray2D<Type>::
    Ray2D(Ray2D<Type> const& source)
      : m_origin(source.m_origin),
        m_direction(source.m_direction)
    {
      // Empty.
    }


    // Destructor.
    template <class Type>
    Ray2D<Type>::
    ~Ray2D()
    {
      // Empty.
    }

      
    // The assignment operator deep copies its argument.
    template <class Type>
    Ray2D<Type>&
    Ray2D<Type>::
    operator=(Ray2D<Type> const& source)
    {
      if(&source != this) {
        m_origin = source.m_origin;
        m_direction = source.m_direction;
      }
      return *this;
    }


    // This member function returns the direction of the ray.
    template <class Type>
    brick::numeric::Vector2D<Type> const&
    Ray2D<Type>::
    getDirection() const
    {
      return m_direction;
    }
      

    // This member function returns the direction of the ray.
    template <class Type>
    brick::numeric::Vector2D<Type> const&
    Ray2D<Type>::
    getDirectionVector() const
    {
      return this->getDirection();
    }
      

    // This member function returns the start point of the ray.
    template <class Type>
    brick::numeric::Vector2D<Type> const&
    Ray2D<Type>::
    getOrigin() const
    {
      return m_origin;
    }


    /* ======= Non-member functions. ======= */

    template <class Type>
    std::ostream&
    operator<<(std::ostream& stream, Ray2D<Type> const& ray)
    {
      stream << "Ray2D { "
             << ray.getOrigin() << ", "
             << ray.getDirectionVector() << " }";
      return stream;
    }
    
  } // namespace geometry
    
} // namespace brick

#endif /* #ifndef BRICK_GEOMETRY_RAY2D_IMPL_HH */
