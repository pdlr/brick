/**
***************************************************************************
* @file brick/numeric/amanatidesWoo3D.hh
*
* Header file declaring AmanatidesWoo3D class.
*
* Copyright (C) 2004-2007,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_AMANATIDESWOO3D_HH
#define BRICK_NUMERIC_AMANATIDESWOO3D_HH

#include <iostream>
#include <limits>
#include <brick/common/exception.hh>
#include <brick/numeric/vector3D.hh>
#include <brick/numeric/transform3D.hh>
#include <brick/numeric/amanatidesWoo3DIterator.hh>

namespace brick {

  namespace numeric {
    
    /**
     ** This class implements the Fast Voxel Traversal Algorithm of
     ** Amanatides and Woo [Ref] for 3D arrays.  The algorithm is for
     ** efficient ray tracing through regular grids.  For convenience,
     ** there are two coordinate systems associated with this class:
     **
     ** - The voxel coordinate system is an integer-based coordinate
     ** system associated with the array.
     **
     ** - The world coordinate system, which can differ from the voxel
     ** coordinate system only by translation, rotation, and
     ** non-isotropic scaling, but is in other respects arbitrary.
     **
     ** The voxel coordinate system has its origin at the corner of the
     ** 3D array, has its first coordinate (which we'll call U) parallel
     ** to the rows of the array, has its second coordinate (which we'll
     ** call V) parallel to the columns of the array, and has its third
     ** coordinate (which we'll call W) parallel to the direction of
     ** increasing slice number.  As you move along any row of any slice
     ** in the 3D array, the transition from the first column to the
     ** second occurs at U == 1.0, the transition from the second column
     ** to the third occurs at U == 2.0, and so on.  Similarly as you
     ** move along any column of any slice of the array, the transition
     ** from the first row to the second occurs at V == 1.0, and the
     ** transition from the second row to the third occurs at V ==
     ** 2.0. At any row, column coordinate, as you move through
     ** subsequent slices of the array, the transition from the first
     ** slice to the second occurs at W == 1.0, and the transition from
     ** the second slice to the third occurs at W == 2.0, and so on.
     **
     ** Once constructed, you can use the AmanatidesWoo3D class to
     ** access voxels along a straight line through the 3D array using
     ** STL-style iteration.  That is, the member functions
     ** AmanatidesWoo3D::begin() returns an iterator which references
     ** subsequent voxels along the line, and will become equal (==) to
     ** the iterator returned by AmanatidesWoo3D::end() when the line
     ** passes out of the 3D array.  Here's an example usage, which will
     ** probably format terribly in the Doxygen-generated doc:
     **
     **   typedef AmanatidesWoo3D< Array3D<int> >::iterator awIterator;
     **
     **   AmanatidesWoo3D< Array3D<int> > rayTracer(
     **     myArray3D, rayOrigin, rayDirection, false);
     **
     **   awIterator iterator0 = rayTracer.begin(); 
     **
     **   awIterator iterator1 = rayTracer.end();
     **
     **   for(; iterator0 != iterator1; ++iterator0) {
     **
     **     *iterator0 = 0;
     **
     **   }
     **/
    template <class ARRAY3D>
    class AmanatidesWoo3D {
    public:
      /* ================== Public typedefs ================== */

      /**
       ** This typedef specifies the type which will be returned by
       ** member functions begin() and end().  AmanatidesWoo3DIterator
       ** is a forward iterator, and has unusual semantics for
       ** operator==().  Be sure to understand them before using this
       ** class.
       **/
      typedef AmanatidesWoo3DIterator<ARRAY3D> iterator;

      // Iteration over const arrays will be done using a different
      // class.
      // typedef AmanatidesWoo3DConstIterator<ARRAY3D> const_iterator;
    
      /* ================== Public methods ================== */

      /** 
       * This constructor specifies all of the internal state of the
       * AmanatidesWoo3D class.  After construction, the AmanatidesWoo3D
       * instance is ready to be used.
       * 
       * @param data This argument specifies the 3D array over which to
       * iterate.
       *
       * @param voxelTworld This parameter specifies a coordinate
       * transformation which takes world coordinates and converts them
       * into voxel coordinates.
       *     
       * @param rayOrigin This parameter specifies the starting point
       * of the ray to be traced, expressed in world coordinates.  If
       * you'd rather express rayOrigin and rayDirection in voxel
       * coordinates, simply set argument voxelTworld to the identity
       * transform.  Note that rayOrigin does not have to lie inside the
       * voxel array.
       *     
       * @param rayDirection This parameter specifies the direction of
       * the ray to be traced, expressed in world coordinates.  In other
       * words, any point on the line of interest can be written as
       * (rayOrigin + t * rayDirection) for some value of t.  If
       * you'd rather express rayOrigin and rayDirection in voxel
       * coordinates, simply set argument voxelTworld to the identity
       * transform.
       *     
       * @param downstreamOnly This boolean argument specifies whether
       * voxels "upstream" of rayOrigin are to be included in the voxel
       * traversal.  If downstreamOnly is set to false, voxels which
       * intersect the line (rayOrigin + t * rayDirection) will be
       * included in the iteration even if they intersect the line at
       * negative values of t.
       */
      AmanatidesWoo3D(ARRAY3D& data,
                      const Transform3D<double>& voxelTworld,
                      const Vector3D<double>& rayOrigin,
                      const Vector3D<double>& rayDirection,
                      bool downstreamOnly=true);

    
      /** 
       * This is the copy constructor.  After copying, the new
       * AmanatidesWoo3D instance and the copied instance both reference
       * the same voxel array.
       * 
       * @param source The AmanatidesWoo3D instance to be copied.
       */
      AmanatidesWoo3D(const AmanatidesWoo3D& source);

      /** 
       * This is the destructor.  It destroys the AmanatidesWoo3D
       * instance and cleans up any allocated memory.
       */
      ~AmanatidesWoo3D();

      /** 
       * This member function returns an iterator which references the
       * first voxel along the traced ray.  If constructor argument
       * downstreamOnly was set to false, this iterator will either
       * point to a voxel on the very edge of the array, or (in the case
       * that the ray does not intersect the voxel array) be equal to
       * this->end().  If constructor argument downstreamOnly was set to
       * true, then this iterator will point to an edge voxel, or point
       * to an interior voxel (in the case that rayOrigin lies within
       * the boundaries of the voxel array), or be equal to this->end()
       * (in the case that the ray does not intersect the voxel array,
       * and in the case that rayOrigin lies outside the voxel array and
       * rayDirection points away from the voxel array).
       * 
       * @return The return value is an iterator pointing to the first
       * voxel in the array, or else an iterator for which (iter ==
       * this->end()) is true.
       */
      iterator
      begin();

      /** 
       * This member function returns an iterator which references an
       * invalid voxel, and which will be equal (==) to the iterator
       * returned by member function begin() when that iterator has
       * fully traversed the line of voxels.
       * 
       * @return The return value is an iterator pointing to an invalid
       * voxel, for use as the end iterator in the comparison clause of
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
      ARRAY3D&
      getData() {return m_data;}

      /** 
       * This member function returns true if the iterator returned by
       * member function begin() will point to a valid voxel.  In other
       * words, if constructor argument downstreamOnly was set to
       * false, the return value of validIntersection() indicates
       * whether the ray specified in the constructor actually
       * intersects the voxel array.  If constructor argument
       * downstreamOnly was set to true, the return value of
       * validIntersection whether the "downstream" half of the ray
       * actually intersects the voxel array.
       * 
       * @return A boolean indicating whether or not (this->begin() ==
       * this->end()) will be true.
       */
      inline bool
      validIntersection();
    
      /** 
       * The assignment operator copies its argument.  After copying,
       * *this and the copied instance both reference the same voxel
       * array.
       * 
       * @param source The AmanatidesWoo3D instance to be copied.
       * @return The return value is a reference to *this.
       */
      AmanatidesWoo3D&
      operator=(const AmanatidesWoo3D& source);
  
    private:

      /** 
       * This private member function computes the parameters tEntry and
       * tExit such that (rayOrigin + tEntry * rayDirection) is very
       * first point of intersection between the ray and the voxel
       * array, and (rayOrigin + tExit * rayDirection) is very last of
       * intersection between the ray and the voxel array.
       * 
       * @param rayOriginVoxel This parameter specifies the starting
       * point of the ray to be traced, expressed in voxel coordinates.
       * @param rayDirectionVoxel This parameters specifies the
       * direction of the ray to be traced, expressed in voxel
       * coordinates.
       * @param m_data This parameter is a reference the voxel array.
       * @return The return value is a std::pair<double, double> in
       * which the first value is tEntry and the second value is tExit.
       */
      std::pair<double, double>
      findEntryAndExitPoints(const Vector3D<double>& rayOriginVoxel,
                             const Vector3D<double>& rayDirectionVoxel,
                             const ARRAY3D& data);

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
      double
      findIntersection(const Vector3D<double>& rayOrigin,
                       const Vector3D<double>& rayDirection,
                       const Vector3D<double>& bVector,
                       double cConstant,
                       double defaultValue);

      /// This data member is a reference the array over which to
      /// iterate.
      ARRAY3D& m_data;

      /// This data member indicates the U coordinate of the first voxel
      /// along the path.
      int m_initialU;

      /// This data member indicates the V coordinate of the first voxel
      /// along the path.
      int m_initialV;

      /// This data member indicates the W coordinate of the first voxel
      /// along the path.
      int m_initialW;
    
      /// This data member indicates whether the U coordinate will be
      /// increasing (m_stepU == 1) or decreasing (m_stepU == -1) as we
      /// travel along the ray path.
      int m_stepU;

      /// This data member indicates whether the V coordinate will be
      /// increasing (m_stepV == 1) or decreasing (m_stepV == -1) as we
      /// travel along the ray path.
      int m_stepV;
    
      /// This data member indicates whether the W coordinate will be
      /// increasing (m_stepW == 1) or decreasing (m_stepW == -1) as we
      /// travel along the ray path.
      int m_stepW;
    
      /// This data member indicates what increment of ray parameter t
      /// moves us exactly 1 voxel in the U direction.
      double m_tDeltaU;

      /// This data member indicates what increment of ray parameter t
      /// moves us exactly 1 voxel in the V direction.
      double m_tDeltaV;
    
      /// This data member indicates what increment of ray parameter t
      /// moves us exactly 1 voxel in the W direction.
      double m_tDeltaW;
    
      /// This data member indicates the value of ray parameter t at
      /// which the ray first enters a voxel with U coordinate not equal
      /// to m_initialU.
      double m_tMaxU;

      /// This data member indicates the value of ray parameter t at
      /// which the ray first enters a voxel with V coordinate not equal
      /// to m_initialV.
      double m_tMaxV;

      /// This data member indicates the value of ray parameter t at
      /// which the ray first enters a voxel with W coordinate not equal
      /// to m_initialW.
      double m_tMaxW;

      /// This data member indicates the value of ray parameter t at
      /// which we start following the ray.  It frequently corresponds
      /// to the point at which the ray enters the voxelated space, but
      /// may also be set to zero if the ray origin is within the
      /// voxelated space and constructor argument downstreamOnly was
      /// set to true.
      double m_tStart;

      /// This data member indicates whether the portion of the ray we
      /// care about actually intersects the specified voxel array.
      bool m_validIntersection;
    };

  } // namespace numeric

} // namespace brick

// Include file containing definitions of inline and template
// functions.
#include <brick/numeric/amanatidesWoo3D_impl.hh>

#endif /* #ifndef BRICK_NUMERIC_AMANATIDESWOO3D_HH */
