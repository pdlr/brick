/**
***************************************************************************
* @file brick/geometry/circle2D.hh
*
* Header file declaring the Circle2D class template.
*
* Copyright (C) 2008-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_GEOMETRY_CIRCLE2D_HH
#define BRICK_GEOMETRY_CIRCLE2D_HH

#include <iostream>
#include <brick/numeric/vector2D.hh>

namespace brick {

  namespace geometry {

    /**
     ** The Circle2D class represents a circle in 2D space.
     **/
    template <class Type>
    class Circle2D {
    public:

      /**
       * The default constructor initializes to the unit circle.
       */
      inline
      Circle2D();


      /**
       * This constructor initializes the circle using explicitly
       * specified values.
       *
       * @param origin This argument specifies the position of the
       * geometric center of the circle.
       *
       * @param radius This argument specifies the radius of the
       * circle.  If the value of radius is less than 0.0, it will be
       * multiplied by -1.
       */
      inline
      Circle2D(brick::numeric::Vector2D<Type> const& origin,
               Type const& radius);


      /**
       * The copy constructor deep copies its argument.
       *
       * @param source This argument is the class instance to be copied.
       */
      inline
      Circle2D(Circle2D<Type> const& source);


      /**
       * Destructor.
       */
      ~Circle2D() {}


      /**
       * The assignment operator deep copies its argument.
       *
       * @param source This argument is the class instance to be copied.
       *
       * @return The return value is a reference to *this.
       */
      inline Circle2D<Type>&
      operator=(Circle2D<Type> const& source);


      /**
       * This member function returns the geometric center of the circle.
       *
       * @return The return value is the point at the center of the
       * circle.
       */
      brick::numeric::Vector2D<Type> const&
      getOrigin() const {return m_origin;}


      /**
       * This member function returns the radius of the circle.
       *
       * @return The return value is the radius of the circle.
       */
      Type const&
      getRadius() const {return m_radius;}


    private:
      // Private member functions.

      // Private data members.
      brick::numeric::Vector2D<Type> m_origin;
      Type m_radius;

    }; // class Circle2D



    /* ======= Non-member functions. ======= */

    std::ostream&
    operator<<(std::ostream& stream, Circle2D<Type> const& circle);

  } // namespace geometry

} // namespace brick


// Include definitions of inline and template functions.
#include <brick/geometry/circle2D_impl.hh>

#endif /* #ifndef BRICK_GEOMETRY_CIRCLE2D_HH */
