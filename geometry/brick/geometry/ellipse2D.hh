/**
***************************************************************************
* @file brick/geometry/ellipse2D.hh
*
* Header file declaring the Ellipse2D class template.
*
* Copyright (C) 2008-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_GEOMETRY_ELLIPSE2D_HH
#define BRICK_GEOMETRY_ELLIPSE2D_HH

#include <iostream>
#include <brick/numeric/array1D.hh>
#include <brick/numeric/vector2D.hh>

namespace brick {

  namespace geometry {
    
    /**
     ** The Ellipse2D class represents a planar ellipse in 2D space.
     **/
    template <class Type>
    class Ellipse2D {
    public:
      
      /** 
       * The default constructor initializes to the unit circle.
       */
      inline
      Ellipse2D();

      
      /** 
       * This constructor initializes the ellipse using explicitly
       * specified values.
       * 
       * @param origin This argument specifies the position of the
       * geometric center of the ellipse.
       * 
       * @param majorAxis This argument represents a vector pointing
       * from the center of the ellipse to one of the two points on
       * the boundary of the ellipse that is farthest from the center.
       * 
       * @param ratio This argument specifies the length of the minor
       * axis as a proportion of the lenth of the major axis.  It must
       * be less than or equal to 1.0.
       */
      inline
      Ellipse2D(brick::numeric::Vector2D<Type> const& origin,
                brick::numeric::Vector2D<Type> const& semimajorAxis,
                Type ratio);

    
      /** 
       * The copy constructor deep copies its argument.
       * 
       * @param source This argument is the class instance to be copied.
       */
      inline
      Ellipse2D(Ellipse2D<Type> const& source);


      /** 
       * Destructor.
       */
      ~Ellipse2D() {}


      /** 
       * The assignment operator deep copies its argument.
       * 
       * @param source This argument is the class instance to be copied.
       * 
       * @return The return value is a reference to *this.
       */
      inline Ellipse2D<Type>&
      operator=(Ellipse2D<Type> const& source);


      /** 
       * Estimate ellipse parameters from a series of points on the
       * ellipse.  After calling this member function, *this will
       * match the input points as closely as possible.  This member
       * function currently implements the closed form algorithm of
       * [1], which is based on Fitzgibbon's formulation [2].  Note
       * that this algorithm minimizes RMS error in an algebraic
       * distance of the ellipse parameterization, rather than
       * minimizaing in Euclidean space, leading to a bias in the
       * resulting estimate.  If this bothers us, we can implement an
       * iterative method.
       *
       * [1] Halir, R., and Flusser, J., "Numerically Stable Direct
       * Least Squares Fitting Of Ellipses." 1998.
       *
       * [2] Fitzgibbon, A. W., Pilu, M, and Fischer, R. B., "Direct
       * Least Squares Fitting of Ellipses." Proc. of the 13th
       * International Conference on Pattern Recognition, pp 253–257,
       * 1996.
       *
       * @param beginIter This argument is an iterator pointing to the
       * beginning of a sequence of Vector2D<Type> instances.
       * 
       * @param endIter This argument is an iterator pointing to the
       * end of a sequence of Vector2D<Type> instances.
       */
      template<class IterType>
      void
      estimate(IterType beginIter, IterType endIter);

      
      /** 
       * This member function returns the geometric center of the ellipse.
       * 
       * @return The return value is the point at the centroid of the
       * ellipse.
       */
      brick::numeric::Vector2D<Type> const&
      getOrigin() const {return m_origin;}
      

      /** 
       * This member function returns a vector pointing from the
       * center of the ellipse to the point on the edge of the ellipse
       * that is farthest from the center.  Note that there are two
       * such farthest points on opposite sides of the ellipse.  The
       * vector returned by this member function will remain
       * consistent for the life of the ellipse, and will reflect the
       * semimajor axis specified as a constructor argument (if the
       * three-argument constructor was used).
       * 
       * @return The return value is a vector pointing along the
       * semimajor axis of the ellipse.
       */
      brick::numeric::Vector2D<Type> const&
      getSemimajorAxis() const {return m_semimajorAxis;}
      

      /** 
       * This member function returns a vector pointing from the
       * center of the ellipse to the point on the edge of the ellipse
       * that is closest to the center.  Note that there are two such
       * closest points on opposite sides of the ellipse.  The vector
       * returned by this member function will remain consistent for
       * the life of the ellipse, and will be normally be rotated 90
       * degrees counterclockwise from the vector returned by
       * getSemiMajorAxis().
       * 
       * @return The return value is a vector pointing along the
       * semiminor axis of the ellipse.
       */
      brick::numeric::Vector2D<Type> const&
      getSemiminorAxis() const {return m_semimajorAxis;}
      

    private:
      // Private member functions.

      // Convert from implicit ellipse representation to trigonometric
      void
      convertAlgebraicToTrigonometric(
        brick::numeric::Array1D<Type> const& algebraicParameters,
        brick::numeric::Vector2D<Type>& origin,
        brick::numeric::Vector2D<Type>& semimajorAxis,
        brick::numeric::Vector2D<Type>& semiminorAxis);
      
      
      // Private data members.
      brick::numeric::Vector2D<Type> m_origin;
      brick::numeric::Vector2D<Type> m_semimajorAxis;
      brick::numeric::Vector2D<Type> m_semiminorAxis;

    }; // class Ellipse2D



    /* ======= Non-member functions. ======= */

    template <class Type>
    std::ostream&
    operator<<(std::ostream& stream, Ellipse2D<Type> const& ellipse);
    
  } // namespace geometry
    
} // namespace brick


// Include definitions of inline and template functions.
#include <brick/geometry/ellipse2D_impl.hh>

#endif /* #ifndef BRICK_GEOMETRY_ELLIPSE2D_HH */
