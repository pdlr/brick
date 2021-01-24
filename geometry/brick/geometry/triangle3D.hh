/**
***************************************************************************
* @file brick/geometry/triangle3D.hh
*
* Header file declaring the Triangle3D class template.
*
* Copyright (C) 2008-2011 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_GEOMETRY_TRIANGLE3D_HH
#define BRICK_GEOMETRY_TRIANGLE3D_HH

#include <iostream>
#include <brick/numeric/vector3D.hh>

namespace brick {

  namespace geometry {

    /**
     ** The Triangle3D class represents a triangle in 3D space.
     **/
    template <class Type>
    class Triangle3D {
    public:

      /**
       * The default constructor initializes to a triangle in the X-Y
       * plane.
       */
      Triangle3D();


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
      Triangle3D(brick::numeric::Vector3D<Type> const& vertex0,
                 brick::numeric::Vector3D<Type> const& vertex1,
                 brick::numeric::Vector3D<Type> const& vertex2);


      /**
       * The copy constructor deep copies its argument.
       *
       * @param source This argument is the class instance to be copied.
       */
      Triangle3D(Triangle3D<Type> const& source);


      /**
       * Destructor.
       */
      ~Triangle3D();


      /**
       * The assignment operator deep copies its argument.
       *
       * @param source This argument is the class instance to be copied.
       *
       * @return The return value is a reference to *this.
       */
      Triangle3D<Type>&
      operator=(Triangle3D<Type> const& source);


      /**
       * This member function returns the one of the three vertices
       * that define the triangle.
       *
       * @return The return value is a Vector3D representing the
       * first vertex of the triangle.
       */
      brick::numeric::Vector3D<Type> const&
      getVertex0() const;


      /**
       * This member function returns the one of the three vertices
       * that define the triangle.
       *
       * @return The return value is a Vector3D representing the
       * second vertex of the triangle.
       */
      brick::numeric::Vector3D<Type> const&
      getVertex1() const;


      /**
       * This member function returns the one of the three vertices
       * that define the triangle.
       *
       * @return The return value is a Vector3D representing the
       * first vertex of the triangle.
       */
      brick::numeric::Vector3D<Type> const&
      getVertex2() const;


      Triangle3D<Type>&
      setValue(brick::numeric::Vector3D<Type> const& vertex0,
               brick::numeric::Vector3D<Type> const& vertex1,
               brick::numeric::Vector3D<Type> const& vertex2) {
        m_vertex0 = vertex0;
        m_vertex1 = vertex1;
        m_vertex2 = vertex2;
        return *this;
      }

    private:

      // Private data members.
      brick::numeric::Vector3D<Type> m_vertex0;
      brick::numeric::Vector3D<Type> m_vertex1;
      brick::numeric::Vector3D<Type> m_vertex2;

    }; // class Triangle3D



    /* ======= Non-member functions. ======= */

    template <class Type>
    std::ostream&
    operator<<(std::ostream& stream, Triangle3D<Type> const& triangle);


    template <class Type>
    std::istream&
    operator>>(std::istream& stream, Triangle3D<Type>& triangle);


  } // namespace utilities

} // namespace brick


// Include definitions of inline and template functions.
#include <brick/geometry/triangle3D_impl.hh>

#endif /* #ifndef BRICK_GEOMETRY_TRIANGLE3D_HH */
