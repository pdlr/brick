/**
***************************************************************************
* @file brick/numeric/amanatidesWoo2D.hh
*
* Header file declaring AmanatidesWoo2D class.
*
* Copyright (C) 2004-2007,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_AMANATIDESWOO2D_HH
#define BRICK_NUMERIC_AMANATIDESWOO2D_HH

#include <iostream>
#include <limits>
#include <brick/common/exception.hh>
#include <brick/numeric/vector2D.hh>
#include <brick/numeric/transform2D.hh>
#include <brick/numeric/utilities.hh>
#include <brick/numeric/amanatidesWoo2DIterator.hh>

namespace brick {

  namespace numeric {
    
    /**
     ** This class implements the Fast Voxel Traversal Algorithm of
     ** Amanatides and Woo [1] for 2D arrays.  The algorithm is for
     ** efficient ray tracing through regular grids.  For convenience,
     ** there are two coordinate systems associated with this class:
     **
     ** - The pixel coordinate system is an integer-based coordinate
     ** system associated with the array.
     **
     ** - The world coordinate system, which can differ from the pixel
     ** coordinate system only by translation, rotation, and
     ** non-isotropic scaling, but is in other respects arbitrary.
     **
     ** The pixel coordinate system has its origin at the corner of the
     ** 2D array, has its first coordinate (which we'll call U) parallel
     ** to the rows of the array, and has its second coordinate (which
     ** we'll call V) parallel to the columns of the array.  As you move
     ** along any row of the 2D array, the transition from the first
     ** column to the second occurs at U == 1.0, the transition from the
     ** second column to the third occurs at U == 2.0, and so on.
     ** Similarly as you move along any column of the array, the
     ** transition from the first row to the second occurs at V == 1.0,
     ** and the transition from the second row to the third occurs at V
     ** == 2.0.
     **
     ** Once constructed, you can use the AmanatidesWoo2D class to
     ** access pixels along a straight line through the 2D array using
     ** STL-style iteration.  That is, the member functions
     ** AmanatidesWoo2D::begin() returns an iterator which references
     ** subsequent pixels along the line, and will become equal (==) to
     ** the iterator returned by AmanatidesWoo2D::end() when the line
     ** passes out of the 2D array.  Here's example usage:
     **
     ** @code
     **   typedef AmanatidesWoo2D< Array2D<int> >::iterator awIterator;
     **   Transform2D<double> pixelFromWorld;
     **   AmanatidesWoo2D< Array2D<int> > rayTracer(
     **     myArray2D, pixelFromWorld, rayOrigin, rayDirection, false);
     **   awIterator iterator0 = rayTracer.begin(); 
     **   awIterator iterator1 = rayTracer.end();
     **   for(; iterator0 != iterator1; ++iterator0) {
     **     *iterator0 = 0;
     **   }
     ** @endcode
     **
     ** [1] John Amanatides and Andrew Woo, "A Fast Voxel Traversal
     ** Algorithm for Ray Tracing", Eurographics ’87, 1987, pp 3-10.
     **/
    template <class ARRAY2D, class FLOAT_TYPE = double, class INT_TYPE = int>
    class AmanatidesWoo2D {
    public:
      /* ================== Public typedefs ================== */

      /**
       ** This typedef specifies the type which will be returned by
       ** member functions begin() and end().  AmanatidesWoo2DIterator
       ** is a forward iterator, and has unusual semantics for
       ** operator==().  Be sure to understand them before using this
       ** class.
       **/
      typedef AmanatidesWoo2DIterator<ARRAY2D> iterator;

      // Iteration over const arrays will be done using a different
      // class.
      // typedef AmanatidesWoo2DConstIterator<ARRAY2D> const_iterator;
    
      /* ================== Public methods ================== */

      /** 
       * This constructor specifies all of the internal state of the
       * AmanatidesWoo2D class.  After construction, the AmanatidesWoo2D
       * instance is ready to be used.
       * 
       * @param data This argument specifies the 2D array over which to
       * iterate.
       *
       * @param pixelTworld This parameter specifies a coordinate
       * transformation which takes world coordinates and converts them
       * into pixel coordinates.
       *     
       * @param rayOrigin This parameter specifies the starting point
       * of the ray to be traced, expressed in world coordinates.  If
       * you'd rather express rayOrigin and rayDirection in pixel
       * coordinates, simply set argument pixelTworld to the identity
       * transform.  Note that rayOrigin does not have to lie inside the
       * pixel array.
       *     
       * @param rayDirection This parameter specifies the direction of
       * the ray to be traced, expressed in world coordinates.  In other
       * words, any point on the line of interest can be written as
       * (rayOrigin + t * rayDirection) for some value of t.  If
       * you'd rather express rayOrigin and rayDirection in pixel
       * coordinates, simply set argument pixelTworld to the identity
       * transform.
       *     
       * @param downstreamOnly This boolean argument specifies whether
       * pixels "upstream" of rayOrigin are to be included in the pixel
       * traversal.  If downstreamOnly is set to false, pixels which
       * intersect the line (rayOrigin + t * rayDirection) will be
       * included in the iteration even if they intersect the line at
       * negative values of t.
       */
      AmanatidesWoo2D(ARRAY2D& data,
                      const Transform2D<FLOAT_TYPE>& pixelTworld,
                      const Vector2D<FLOAT_TYPE>& rayOrigin,
                      const Vector2D<FLOAT_TYPE>& rayDirection,
                      bool downstreamOnly=true);

    
      /** 
       * This is the copy constructor.  After copying, the new
       * AmanatidesWoo2D instance and the copied instance both reference
       * the same pixel array.
       * 
       * @param source The AmanatidesWoo2D instance to be copied.
       */
      AmanatidesWoo2D(const AmanatidesWoo2D& source);

      /** 
       * This is the destructor.  It destroys the AmanatidesWoo2D
       * instance and cleans up any allocated memory.
       */
      ~AmanatidesWoo2D();

      /** 
       * This member function returns an iterator which references the
       * first pixel along the traced ray.  If constructor argument
       * downstreamOnly was set to false, this iterator will either
       * point to a pixel on the very edge of the array, or (in the case
       * that the ray does not intersect the pixel array) be equal to
       * this->end().  If constructor argument downstreamOnly was set to
       * true, then this iterator will point to an edge pixel, or point
       * to an interior pixel (in the case that rayOrigin lies within
       * the boundaries of the pixel array), or be equal to this->end()
       * (in the case that the ray does not intersect the pixel array,
       * and in the case that rayOrigin lies outside the pixel array and
       * rayDirection points away from the pixel array).
       * 
       * @return The return value is an iterator pointing to the first
       * pixel in the array, or else an iterator for which (iter ==
       * this->end()) is true.
       */
      iterator
      begin();

      /** 
       * This member function returns an iterator which references an
       * invalid pixel, and which will be equal (==) to the iterator
       * returned by member function begin() when that iterator has
       * fully traversed the line of pixels.
       * 
       * @return The return value is an iterator pointing to an invalid
       * pixel, for use as the end iterator in the comparison clause of
       * a loop: "while(rayIterator != endIterator) {++rayIterator; ...;}".
       */
      iterator
      end();

      /** 
       * This member function returns a reference to the array object
       * over which iteration is performed.
       * 
       * @return The return value is a reference to the data array.
       */
      ARRAY2D&
      getData() {return m_data;}
    
      /** 
       * This member function returns true if the iterator returned by
       * member function begin() will point to a valid pixel.  In other
       * words, if constructor argument downstreamOnly was set to
       * false, the return value of validIntersection() indicates
       * whether the ray specified in the constructor actually
       * intersects the pixel array.  If constructor argument
       * downstreamOnly was set to true, the return value of
       * validIntersection whether the "downstream" half of the ray
       * actually intersects the pixel array.
       * 
       * @return A boolean indicating whether or not (this->begin() ==
       * this->end()) will be true.
       */
      inline bool
      validIntersection();

      /** 
       * The assignment operator copies its argument.  After copying,
       * *this and the copied instance both reference the same pixel
       * array.
       * 
       * @param source The AmanatidesWoo2D instance to be copied.
       * @return The return value is a reference to *this.
       */
      AmanatidesWoo2D&
      operator=(const AmanatidesWoo2D& source);
  
    private:

      /** 
       * This private member function computes the parameters tEntry and
       * tExit such that (rayOrigin + tEntry * rayDirection) is very
       * first point of intersection between the ray and the pixel
       * array, and (rayOrigin + tExit * rayDirection) is very last of
       * intersection between the ray and the pixel array.
       * 
       * @param rayOriginPixel This parameter specifies the starting
       * point of the ray to be traced, expressed in pixel coordinates.
       * @param rayDirectionPixel This parameters specifies the
       * direction of the ray to be traced, expressed in pixel
       * coordinates.
       * @param m_data This parameter is a reference the pixel array.
       * @return The return value is a std::pair<FLOAT_TYPE, FLOAT_TYPE> in
       * which the first value is tEntry and the second value is tExit.
       */
      std::pair<FLOAT_TYPE, FLOAT_TYPE>
      findEntryAndExitPoints(const Vector2D<FLOAT_TYPE>& rayOriginPixel,
                             const Vector2D<FLOAT_TYPE>& rayDirectionPixel,
                             const ARRAY2D& data);

      /** 
       * This function is not documented because it will change very
       * soon.
       * 
       * @param rayOrigin 
       * @param rayDirection 
       * @param bVector 
       * @param cConstant 
       * @param defaultValue 
       * @return 
       */
      FLOAT_TYPE
      findIntersection(const Vector2D<FLOAT_TYPE>& rayOrigin,
                       const Vector2D<FLOAT_TYPE>& rayDirection,
                       const Vector2D<FLOAT_TYPE>& bVector,
                       FLOAT_TYPE cConstant,
                       FLOAT_TYPE defaultValue);

      /// This data member is a reference the array over which to
      /// iterate.
      ARRAY2D& m_data;

      /// This data member indicates the U coordinate of the first pixel
      /// along the path.
      INT_TYPE m_initialU;

      /// This data member indicates the V coordinate of the first pixel
      /// along the path.
      INT_TYPE m_initialV;

      /// This data member indicates whether the U coordinate will be
      /// increasing (m_stepU == 1) or decreasing (m_stepU == -1) as we
      /// travel along the ray path.
      INT_TYPE m_stepU;

      /// This data member indicates whether the V coordinate will be
      /// increasing (m_stepV == 1) or decreasing (m_stepV == -1) as we
      /// travel along the ray path.
      INT_TYPE m_stepV;
    
      /// This data member indicates what increment of ray parameter t
      /// moves us exactly 1 pixel in the U direction.
      FLOAT_TYPE m_tDeltaU;

      /// This data member indicates what increment of ray parameter t
      /// moves us exactly 1 pixel in the V direction.
      FLOAT_TYPE m_tDeltaV;
    
      /// This data member indicates the value of ray parameter t at
      /// which the ray first enters a pixel with U coordinate not equal
      /// to m_initialU.
      FLOAT_TYPE m_tMaxU;

      /// This data member indicates the value of ray parameter t at
      /// which the ray first enters a pixel with V coordinate not equal
      /// to m_initialV.
      FLOAT_TYPE m_tMaxV;

      /// This data member indicates the value of ray parameter t at
      /// which we start following the ray.  It frequently corresponds
      /// to the point at which the ray enters the pixelated space, but
      /// may also be set to zero if the ray origin is within the
      /// pixelated space and constructor argument downstreamOnly was
      /// set to true.
      FLOAT_TYPE m_tStart;

      /// This data member indicates whether the portion of the ray we
      /// care about actually intersects the specified pixel array.
      bool m_validIntersection;
    };

  } // namespace numeric

} // namespace brick

// Include file containing definitions of inline and template
// functions.
#include <brick/numeric/amanatidesWoo2D_impl.hh>

#endif /* #ifndef BRICK_NUMERIC_AMANATIDESWOO2D_HH */
