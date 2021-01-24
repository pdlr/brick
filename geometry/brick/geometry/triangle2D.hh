/**
***************************************************************************
* @file brick/geometry/triangle2D.hh
*
* Header file declaring the Triangle2D class template.
*
* Copyright (C) 2014 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_GEOMETRY_TRIANGLE2D_HH
#define BRICK_GEOMETRY_TRIANGLE2D_HH

#include <iostream>
#include <brick/numeric/vector2D.hh>

namespace brick {

  namespace geometry {

    /**
     ** The Triangle2D class represents a triangle in 2D space.
     **/
    template <class Type>
    class Triangle2D {
    public:

      /**
       * The default constructor initializes to a right triangle with
       * two of its legs having length 1.0 and running from the origin
       * along the X and Y axes.
       */
      Triangle2D();


      /**
       * This constructor initializes the triangle using three points.
       *
       * @param vertex0 This argument is one of the three points that
       * define the triangle.
       *
       * @param vertex1 This argument is one of the three points that
       * define the triangle.
       *
       * @param vertex2 This argument is one of the three points that
       * define the triangle.
       */
      Triangle2D(brick::numeric::Vector2D<Type> const& vertex0,
                 brick::numeric::Vector2D<Type> const& vertex1,
                 brick::numeric::Vector2D<Type> const& vertex2);


      /**
       * The copy constructor deep copies its argument.
       *
       * @param source This argument is the class instance to be copied.
       */
      Triangle2D(Triangle2D<Type> const& source);


      /**
       * Destructor.
       */
      ~Triangle2D();


      /**
       * The assignment operator deep copies its argument.
       *
       * @param source This argument is the class instance to be copied.
       *
       * @return The return value is a reference to *this.
       */
      Triangle2D<Type>&
      operator=(Triangle2D<Type> const& source);


      /**
       * Return the total area of the triangle
       *
       * @return The return value is the area of the triangle.
       */
      Type
      getArea() const;


      /**
       * This member function returns the one of the three vertices
       * that define the triangle.
       *
       * @return The return value is a Vector2D representing the
       * first vertex of the triangle.
       */
      brick::numeric::Vector2D<Type> const&
      getVertex0() const;


      /**
       * This member function returns the one of the three vertices
       * that define the triangle.
       *
       * @return The return value is a Vector2D representing the
       * second vertex of the triangle.
       */
      brick::numeric::Vector2D<Type> const&
      getVertex1() const;


      /**
       * This member function returns the one of the three vertices
       * that define the triangle.
       *
       * @return The return value is a Vector2D representing the
       * first vertex of the triangle.
       */
      brick::numeric::Vector2D<Type> const&
      getVertex2() const;


      Triangle2D<Type>&
      setValue(brick::numeric::Vector2D<Type> const& vertex0,
               brick::numeric::Vector2D<Type> const& vertex1,
               brick::numeric::Vector2D<Type> const& vertex2) {
        m_vertex0 = vertex0;
        m_vertex1 = vertex1;
        m_vertex2 = vertex2;
        return *this;
      }

    private:

      // Private data members.
      brick::numeric::Vector2D<Type> m_vertex0;
      brick::numeric::Vector2D<Type> m_vertex1;
      brick::numeric::Vector2D<Type> m_vertex2;

    }; // class Triangle2D



    /* ======= Non-member functions. ======= */

    template <class Type>
    std::ostream&
    operator<<(std::ostream& stream, Triangle2D<Type> const& triangle);


    template <class Type>
    std::istream&
    operator>>(std::istream& stream, Triangle2D<Type>& triangle);


  } // namespace utilities

} // namespace brick


// Include definitions of inline and template functions.
#include <brick/geometry/triangle2D_impl.hh>

#endif /* #ifndef BRICK_GEOMETRY_TRIANGLE2D_HH */
