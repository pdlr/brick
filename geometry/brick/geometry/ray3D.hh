/**
***************************************************************************
* @file brick/geometry/ray3D.hh
*
* Header file declaring Ray3D class template.
*
* Copyright (C) 2007 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_GEOMETRY_RAY3D_HH
#define BRICK_GEOMETRY_RAY3D_HH

#include <iostream>
#include <brick/numeric/vector3D.hh>

namespace brick {

  namespace geometry {
    
    /**
     ** The Ray3D class represents a ray in 3D space.
     **/
    template <class Type>
    class Ray3D {
    public:
      
      /** 
       * The default constructor initializes to the ray that starts
       * at the origin and points along the X axis.
       */
      Ray3D();

    
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
      Ray3D(brick::numeric::Vector3D<Type> const& point,
            brick::numeric::Vector3D<Type> const& direction,
            bool normalize = true);

    
      /** 
       * The copy constructor deep copies its argument.
       * 
       * @param source This argument is the class instance to be
       * copied.
       */
      Ray3D(Ray3D<Type> const& source);


      /** 
       * Destructor.
       */
      ~Ray3D();

      
      /** 
       * The assignment operator deep copies its argument.
       * 
       * @param source This argument is the class instance to be
       * copied.
       * 
       * @return The return value is a reference to *this.
       */
      Ray3D<Type>&
      operator=(Ray3D<Type> const& source);


      /** 
       * This member function returns the direction of the ray.
       * 
       * @return The return value is a const reference to Vector3D
       * representing the direction of the ray.
       */
      brick::numeric::Vector3D<Type> const&
      getDirection() const;
      

      /** 
       * This member function returns the direction of the ray.
       * 
       * @return The return value is a const reference to Vector3D
       * representing the direction of the ray.
       */
      brick::numeric::Vector3D<Type> const&
      getDirectionVector() const;
      

      /** 
       * This member function returns the start point of the ray.
       * 
       * @return The return value is a const reference to Vector3D
       * representing the start point of the ray.
       */
      brick::numeric::Vector3D<Type> const&
      getOrigin() const;


   private:
      // Private member functions.

      // Private data members.
      brick::numeric::Vector3D<Type> m_origin;
      brick::numeric::Vector3D<Type> m_direction;

    }; // class Ray3D


    /* ======= Non-member functions. ======= */

    template <class Type>
    std::ostream&
    operator<<(std::ostream& stream, Ray3D<Type> const& ray);
    
    
  } // namespace geometry
    
} // namespace brick


// Include definitions of inline and template functions.
#include <brick/geometry/ray3D_impl.hh>
    
#endif /* #ifndef BRICK_GEOMETRY_RAY3D_HH */
