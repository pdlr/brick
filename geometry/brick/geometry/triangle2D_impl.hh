/**
***************************************************************************
* @file brick/geometry/triangle2D_impl.hh
*
* Source file defining the Triangle2D class template.
*
* Copyright (C) 2014 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_GEOMETRY_TRIANGLE2D_IMPL_HH
#define BRICK_GEOMETRY_TRIANGLE2D_IMPL_HH

// This file is included by triangle2D.hh, and should not be directly included
// by user code, so no need to include triangle2D.hh here.
// 
// #include <brick/geometry/triangle2D.hh>

#include <brick/common/expect.hh>
#include <brick/common/mathFunctions.hh>
#include <brick/numeric/utilities.hh>

namespace brick {

  namespace geometry {
    
    // The default constructor initializes to a triangle in the X-Y
    // plane.
    template <class Type>
    Triangle2D<Type>::
    Triangle2D()
      : m_vertex0(0.0, 0.0),
        m_vertex1(1.0, 0.0),
        m_vertex2(0.0, 1.0)
    {
      // Empty.
    }

    
    // This constructor initializes the triangle using three points.
    template <class Type>
    Triangle2D<Type>::
    Triangle2D(brick::numeric::Vector2D<Type> const& vertex0,
               brick::numeric::Vector2D<Type> const& vertex1,
               brick::numeric::Vector2D<Type> const& vertex2)
      : m_vertex0(vertex0),
        m_vertex1(vertex1),
        m_vertex2(vertex2)
    {
      // Empty.
    }

    
    // The copy constructor deep copies its argument.
    template <class Type>
    Triangle2D<Type>::
    Triangle2D(Triangle2D<Type> const& source)
      : m_vertex0(source.m_vertex0),
        m_vertex1(source.m_vertex1),
        m_vertex2(source.m_vertex2)
    {
      // Empty.
    }


    // Destructor.
    template <class Type>
    Triangle2D<Type>::
    ~Triangle2D()
    {
      // Empty.
    }


    // The assignment operator deep copies its argument.
    template <class Type>
    Triangle2D<Type>&
    Triangle2D<Type>::
    operator=(Triangle2D<Type> const& source)
    {
      if(&source != this) {
        m_vertex0 = source.m_vertex0;
        m_vertex1 = source.m_vertex1;
        m_vertex2 = source.m_vertex2;
      }
      return *this;
    }


    // Return the total area of the triangle
    template <class Type>
    Type
    Triangle2D<Type>::
    getArea() const
    {
      // Calculate the area one way.
      Type area0 =
        (this->m_vertex0.getX() * (this->m_vertex1.getY()
                                   - this->m_vertex2.getY())
         + this->m_vertex1.getX() * (this->m_vertex2.getY()
                                     - this->m_vertex0.getY())
         + this->m_vertex2.getX() * (this->m_vertex0.getY()
                                     - this->m_vertex1.getY())) * 0.5;

      // Calculate the area another way.
      Type area1 =
        (this->m_vertex0.getY() * (this->m_vertex1.getX()
                                   - this->m_vertex2.getX())
         + this->m_vertex1.getY() * (this->m_vertex2.getX()
                                     - this->m_vertex0.getX())
         + this->m_vertex2.getY() * (this->m_vertex0.getX()
                                     - this->m_vertex1.getX())) * 0.5;

      // Take the average.
      return (brick::common::absoluteValue(area0)
              + brick::common::absoluteValue(area1)) * 0.5;
    }


    // This member function returns the one of the three vertices
    // that define the triangle.
    template <class Type>
    brick::numeric::Vector2D<Type> const&
    Triangle2D<Type>::
    getVertex0() const
    {
      return m_vertex0;
    }
      

    // This member function returns the one of the three vertices
    // that define the triangle.
    template <class Type>
    brick::numeric::Vector2D<Type> const&
    Triangle2D<Type>::
    getVertex1() const
    {
      return m_vertex1;
    }
      

    // This member function returns the one of the three vertices
    // that define the triangle.
    template <class Type>
    brick::numeric::Vector2D<Type> const&
    Triangle2D<Type>::
    getVertex2() const
    {
      return m_vertex2;
    }
      

    /* ======= Non-member functions. ======= */

    template <class Type>
    std::ostream&
    operator<<(std::ostream& stream, Triangle2D<Type> const& triangle)
    {
      stream << "Triangle2D { "
             << triangle.getVertex0() << ", "
             << triangle.getVertex1() << ", "
             << triangle.getVertex2() << " }";
      return stream;
    }


    template <class Type>
    std::istream&
    operator>>(std::istream& stream, Triangle2D<Type>& triangle)
    {
      // If stream is in a bad state, we can't read from it.
      if (!stream){
        return stream;
      }
    
      // It's a lot easier to use a try block than to be constantly
      // testing whether the IO has succeeded, so we tell inputStream to
      // complain if anything goes wrong.
      std::ios_base::iostate oldExceptionState = stream.exceptions();
      stream.exceptions(
        std::ios_base::badbit | std::ios_base::failbit | std::ios_base::eofbit);

      // Now on with the show.
      try{
        brick::numeric::Vector2D<Type> vertex0;
        brick::numeric::Vector2D<Type> vertex1;
        brick::numeric::Vector2D<Type> vertex2;
        
        brick::common::Expect::FormatFlag flags =
          brick::common::Expect::SkipWhitespace;
        stream >> brick::common::Expect("Triangle2D", flags)
               >> brick::common::Expect("{", flags)
               >> vertex0
               >> brick::common::Expect(",", flags)
               >> vertex1
               >> brick::common::Expect(",", flags)
               >> vertex2
               >> brick::common::Expect("}");

        triangle.setValue(vertex0, vertex1, vertex2);
        
      } catch(std::ios_base::failure) {
        // Empty
      }
      stream.exceptions(oldExceptionState);
      return stream;
    }

    
  } // namespace utilities
    
} // namespace brick

#endif /* #ifndef BRICK_GEOMETRY_TRIANGLE2D_IMPL_HH */
