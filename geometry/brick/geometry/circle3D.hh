/**
***************************************************************************
* @file brick/geometry/circle3D.hh
*
* Header file declaring the Circle3D class template.
*
* Copyright (C) 2014 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_GEOMETRY_CIRCLE3D_HH
#define BRICK_GEOMETRY_CIRCLE3D_HH

#include <iostream>
#include <brick/numeric/vector3D.hh>

namespace brick {

  namespace geometry {

    /**
     ** The Circle3D class represents a 2D circle in 3D space.
     **/
    template <class Type>
    class Circle3D {
    public:

      /**
       * The default constructor initializes to the unit circle in the
       * X-Y plane.
       */
      inline
      Circle3D();


      /**
       * This constructor initializes the circle using explicitly
       * specified values.
       *
       * @param origin This argument specifies the position of the
       * geometric center of the circle.
       *
       * @param basisVector0 This argument specifies a vector that
       * points from the origin to the perimeter of the circle.  The
       * length of this vector defines the radius of the circle, and
       * its orientation specifies a "zero point" along the perimeter.
       * Calls to this->getPerimiterPoint(0.0) will return (origin +
       * basisVector0).
       *
       * @param basisVector1 This argument is used only to constrain
       * the plane in which the circle lies.  It must not be zero
       * length, and must not be parallel to basisVector1.
       *
       * @param isOrthonormalized This argument allows the user to
       * save some computation if basisVector1 is already
       * perpendicular to basisVector0, and has length equal to the
       * radius of the circle.  Be careful with this argument: if you
       * set it to true, and basisVector1 is not already perpendicular
       * to basisVector0, or basisVector1 is not the same length as
       * basisVector0, then the resulting Circle3D instance won't
       * actually describe a circle.
       */
      inline
      Circle3D(brick::numeric::Vector3D<Type> const& origin,
               brick::numeric::Vector3D<Type> const& basisVector0,
               brick::numeric::Vector3D<Type> const& basisVector1,
               bool isOrthonormalized = false);


      /**
       * The copy constructor deep copies its argument.
       *
       * @param source This argument is the class instance to be copied.
       */
      inline
      Circle3D(Circle3D<Type> const& source);


      /**
       * Destructor.
       */
      ~Circle3D() {}


      /**
       * The assignment operator deep copies its argument.
       *
       * @param source This argument is the class instance to be copied.
       *
       * @return The return value is a reference to *this.
       */
      inline Circle3D<Type>&
      operator=(Circle3D<Type> const& source);


      /**
       * This member function returns the geometric center of the circle.
       *
       * @return The return value is the point at the center of the
       * circle.
       */
      brick::numeric::Vector3D<Type> const&
      getBasisVector(brick::common::UInt32 basisIndex) const;


      /**
       * This member function returns the geometric center of the circle.
       *
       * @return The return value is the point at the center of the
       * circle.
       */
      brick::numeric::Vector3D<Type> const&
      getOrigin() const {return m_origin;}


      /**
       * This member function returns a point on the perimiter of the
       * circle, rotated a user-specified angle around the circle,
       * from an arbitrary-but-unchanging start point.
       *
       * @param angle This argument specifies which point on the
       * perimeter of the circle will be returned.
       *
       * @return The return value is the point at the center of the
       * circle.
       */
      inline brick::numeric::Vector3D<Type>
      getPerimeterPoint(Type angle) const;


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
      brick::numeric::Vector3D<Type> m_origin;
      Type m_radius;
      brick::numeric::Vector3D<Type> m_basisVector0;
      brick::numeric::Vector3D<Type> m_basisVector1;

    }; // class Circle3D



    /* ======= Non-member functions. ======= */

    // std::ostream&
    // operator<<(std::ostream& stream, Circle3D<Type> const& circle);

  } // namespace geometry

} // namespace brick


// Include definitions of inline and template functions.
#include <brick/geometry/circle3D_impl.hh>

#endif /* #ifndef BRICK_GEOMETRY_CIRCLE3D_HH */
