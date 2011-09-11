/**
***************************************************************************
* @file brick/geometry/ray3D.cc
*
* Source file defining the Ray3D class.
*
* Copyright (C) 2007 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/geometry/ray3D.hh>
#include <brick/numeric/utilities.hh>

namespace dnum = brick::numeric;

namespace brick {

  namespace geometry {
    
      
    // The default constructor initializes to the ray that starts
    // at the origin and points along the X axis.
    template <class Type>
    Ray3D<Type>::
    Ray3D()
      : m_origin(0.0, 0.0, 0.0), m_direction(1.0, 0.0, 0.0)
    {
      // Empty.
    }

    
    // This constructor initializes the ray using a point and a
    // direction vector.
    template <class Type>
    Ray3D<Type>::
    Ray3D(brick::numeric::Vector3D<Type> const& point,
          brick::numeric::Vector3D<Type> const& direction,
          bool normalize = true)
      : m_origin(point),
        m_direction(direction)
    {
      if(normalize) {
        m_direction /= dnum::magnitude(m_direction);
      }
    }

    
    // The copy constructor deep copies its argument.
    template <class Type>
    Ray3D<Type>::
    Ray3D(Ray3D<Type> const& source)
      : m_origin(source.m_origin),
        m_direction(source.m_direction)
    {
      // Empty.
    }


    // Destructor.
    template <class Type>
    Ray3D<Type>::
    ~Ray3D()
    {
      // Empty.
    }

      
    // The assignment operator deep copies its argument.
    template <class Type>
    Ray3D<Type>&
    Ray3D<Type>::
    operator=(Ray3D<Type> const& source)
    {
      if(&source != this) {
        m_origin = source.m_origin;
        m_direction = source.m_direction;
      }
      return *this;
    }
      

    // This member function returns the direction of the ray.
    template <class Type>
    brick::numeric::Vector3D<Type> const&
    Ray3D<Type>::
    getDirection() const
    {
      return m_direction;
    }
      

    // This member function returns the direction of the ray.
    template <class Type>
    brick::numeric::Vector3D<Type> const&
    Ray3D<Type>::
    getDirectionVector() const
    {
      return this->getDirection();
    }
      

    // This member function returns the start point of the ray.
    template <class Type>
    brick::numeric::Vector3D<Type> const&
    Ray3D<Type>::
    getOrigin() const
    {
      return m_origin;
    }


    // The copy constructor deep copies its argument.
    template <class Type>
    Ray3D<Type>::
    Ray3D(const Ray3D& source)
      : m_origin(source.m_origin),
        m_direction(source.m_direction)
    {
      // Empty.
    }


    /* ======= Non-member functions. ======= */

    template <class Type>
    std::ostream&
    operator<<(std::ostream& stream, const Ray3D<Type>& ray)
    {
      stream << "Ray3D{ "
             << ray.getOrigin() << ", "
             << ray.getDirectionVector() << " }";
      return stream;
    }
    
  } // namespace geometry
    
} // namespace brick
