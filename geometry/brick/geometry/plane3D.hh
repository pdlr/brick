/**
***************************************************************************
* @file brick/geometry/plane3D.hh
*
* Header file declaring the Plane3D class template.
*
* Copyright (C) 2007-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_GEOMETRY_PLANE3D_HH
#define BRICK_GEOMETRY_PLANE3D_HH

#include <iostream>
#include <brick/numeric/vector3D.hh>

// TBD - do we really need this include?
// #include <brick/numeric/utilities.hh>

namespace brick {

  namespace geometry {
    
    /**
     ** The Plane3D class represents a plane in 3D space.
     **/
    template <class Type>
    class Plane3D {
    public:
      
      /** 
       * The default constructor initializes to the X-Y plane.
       */
      Plane3D();

    
      /** 
       * This constructor initializes the plane using three points.
       * 
       * @param point0 This argument is one of the three points that
       * define the plane.
       * 
       * @param point1 This argument is one of the three points that
       * define the plane.
       * 
       * @param point2 This argument is one of the three points that
       * define the plane.
       * 
       * @param orthonormalize If the two vectors (point1 - point0)
       * and (point2 - point0) are orthonormal, then you can save some
       * computation by setting this argument to false.
       */
      Plane3D(Vector3D<Type> const& point0,
              Vector3D<Type> const& point1,
              Vector3D<Type> const& point2,
              bool orthonormalize = true);

    
      /** 
       * This constructor initializes the plane using a collection of
       * points.  The resulting plane minimizes the sum-of-squares
       * residual.
       * 
       * @param beginIterator This argument points to the first
       * Vector3D instance in the sequence.
       * 
       * @param endIterator This argument points to one past the end
       * of the input sequence.
       * 
       * @param inlierProportion This parameter specifies what
       * proportion of the input is expected to actually lie on the
       * plane.  The number of points to ignore is proportional to
       * (1.0 - inlierPercentage).
       */
      template <class Iterator> 
      Plane3D(Iterator beginIterator,
              Iterator endIterator,
              double inlierPercentage = 1.0);

    
      /** 
       * The copy constructor deep copies its argument.
       * 
       * @param source This argument is the class instance to be copied.
       */
      Plane3D(Plane3D<Type> const& source);


      /** 
       * Destructor.
       */
      ~Plane3D() {}


      /** 
       * The assignment operator deep copies its argument.
       * 
       * @param source This argument is the class instance to be copied.
       * 
       * @return The return value is a reference to *this.
       */
      Plane3D<Type>&
      operator=(Plane3D<Type> const& source);


      /** 
       * This member function returns one of a pair of orthonormal
       * direction vectors that span the plane.  Note that there are
       * infinitely many such pairs, so the vector returned by
       * getDirectionVector0() is in some sense arbitrary.  It will,
       * however, be consistent through the lifetime of the Plane3D
       * instance.
       * 
       * @return The return value is a Vector3D representing the
       * first direction of the pair.
       */
      Vector3D<Type> const&
      getDirectionVector0() const;
      

      /** 
       * This member function returns one of a pair of orthonormal
       * direction vectors that span the plane.  Note that there are
       * infinitely many such pairs, so the vector returned by
       * getDirectionVector1() is in some sense arbitrary.  It will,
       * however, be consistent through the lifetime of the Plane3D
       * instance.
       * 
       * @return The return value is a Vector3D representing the
       * second direction of the pair.
       */
      Vector3D<Type> const&
      getDirectionVector1() const;


      // xxx
      Vector3D<Type>
      getNormal() const;
      

      /** 
       * This member function returns one of the infinitely many
       * points on the plane that could serve as the origin of a 2D
       * coordinate system.
       * 
       * @return The return value is a Vector3D representing the
       * returned point.
       */
      Vector3D<Type> const&
      getOrigin() const;


      // xxx Should this be moved to utilities3D.hh?
      double
      findDistance(Vector3D<Type> const& point);
        
    private:
      // Private member functions.
      template <class Iterator> 
      void
      estimateFromSequence(Iterator beginIterator,
                           Iterator endIterator);


      // Private data members.
      Vector3D<Type> m_origin;
      Vector3D<Type> m_directionVector0;
      Vector3D<Type> m_directionVector1;

    }; // class Plane3D



    /* ======= Non-member functions. ======= */

    std::ostream&
    operator<<(std::ostream& stream, Plane3D<Type> const& plane);
    
    
  } // namespace utilities
    
} // namespace brick


// Include definitions of inline and template functions.
#include <brick/geometry/plane3D_impl.hh>

#endif /* #ifndef BRICK_GEOMETRY_PLANE3D_HH */
