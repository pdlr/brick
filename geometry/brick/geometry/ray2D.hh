/**
***************************************************************************
* @file brick/geometry/ray2D.hh
*
* Header file declaring Ray2D class template.
*
* Copyright (C) 2007-2011 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_GEOMETRY_RAY2D_HH
#define BRICK_GEOMETRY_RAY2D_HH

#include <iostream>
#include <brick/numeric/vector2D.hh>

namespace brick {

  namespace geometry {

    /**
     ** The Ray2D class represents a ray in 2D space.
     **/
    template <class Type>
    class Ray2D {
    public:

      /**
       * The default constructor initializes to the ray that starts
       * at the origin and points along the X axis.
       */
      Ray2D();


      /**
       * This constructor initializes the ray using a point and a direction.
       *
       * @param point This argument specifies the start point of the ray.
       *
       * @param direction This argument specifies the direction of the ray.
       *
       * @param normalized If the direction vector is already
       * normalized to unit length, then you can save some computation
       * by setting this argument to false.
       */
      Ray2D(brick::numeric::Vector2D<Type> const& point,
            brick::numeric::Vector2D<Type> const& direction,
            bool normalize = true);


      /**
       * This constructor initializes the ray according to the
       * equation Ax + By + C = 0.  After construction, the ray origin
       * is set to [-aa*cc / k, -bb*cc / k], where k = aa*aa + bb*bb,
       * and the ray direction is set parallel to [-bb, aa].
       *
       * @param aa This argument is the "A" coefficient of the line equation.
       *
       * @param bb This argument is the "B" coefficient of the line equation.
       *
       * @param cc This argument is the "C" coefficient of the line equation.
       */
      Ray2D(Type const& aa, Type const& bb, Type const& cc);


      /**
       * The copy constructor deep copies its argument.
       *
       * @param source This argument is the class instance to be
       * copied.
       */
      Ray2D(Ray2D<Type> const& source);


      /**
       * Destructor.
       */
      ~Ray2D();


      /**
       * The assignment operator deep copies its argument.
       *
       * @param source This argument is the class instance to be
       * copied.
       *
       * @return The return value is a reference to *this.
       */
      Ray2D<Type>&
      operator=(Ray2D<Type> const& source);


      /**
       * This member function returns the direction of the ray.
       *
       * @return The return value is a const reference to Vector2D
       * representing the direction of the ray.
       */
      brick::numeric::Vector2D<Type> const&
      getDirection() const;


      /**
       * This member function returns the direction of the ray.
       *
       * @return The return value is a const reference to Vector2D
       * representing the direction of the ray.
       */
      brick::numeric::Vector2D<Type> const&
      getDirectionVector() const;


      /**
       * This member function returns the start point of the ray.
       *
       * @return The return value is a const reference to Vector2D
       * representing the start point of the ray.
       */
      brick::numeric::Vector2D<Type> const&
      getOrigin() const;


   private:
      // Private member functions.

      // Private data members.
      brick::numeric::Vector2D<Type> m_origin;
      brick::numeric::Vector2D<Type> m_direction;

    }; // class Ray2D



    /* ======= Non-member functions. ======= */

    template <class Type>
    std::ostream&
    operator<<(std::ostream& stream, Ray2D<Type> const& ray);


  } // namespace geometry

} // namespace brick


// Include definitions of inline and template functions.
#include <brick/geometry/ray2D_impl.hh>

#endif /* #ifndef BRICK_GEOMETRY_RAY2D_HH */
