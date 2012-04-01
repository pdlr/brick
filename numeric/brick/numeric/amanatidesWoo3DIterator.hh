/**
***************************************************************************
* @file brick/numeric/amanatidesWoo3DIterator.hh
*
* Header file declaring AmanatidesWoo3DIterator class.
*
* Copyright (C) 2004-2007,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_AMANATIDESWOO3DITERATOR_HH
#define BRICK_NUMERIC_AMANATIDESWOO3DITERATOR_HH

#include <iostream>
#include <brick/common/exception.hh>

namespace brick {

  namespace numeric {
    
    /**
     ** This class provides access to the elements of a data array along
     ** a straight path, and does the actual work of Amanatides and
     ** Woo's fast voxel traversal algorithm.  Typically, an
     ** AmanatidesWoo3DIterator instance will be created by an
     ** AmanatidesWoo3D object in order to access a line of voxels
     ** specified through the AmanatidesWoo3D class interface.  The user
     ** will probably never need to directly construct an
     ** AmanatidesWoo3DIterator.  For more information on the fast voxel
     ** traversal algorithm of Amanatides and Woo, please refer to
     ** [ref].
     **/
    template <class ARRAY3D, class FLOAT_TYPE = double, class INT_TYPE = int>
    class AmanatidesWoo3DIterator
      : public std::iterator<std::forward_iterator_tag,
                             typename ARRAY3D::value_type>
    {
    public:
      /** 
       * The class constructor is initialized with all of the internal
       * variables of the voxel traversal algorithm.
       * 
       * @param data This parameter is a reference to the 3D data over
       * which to iterate.
       * @param startU This parameter specifies the starting U coordinate
       * (column) in the voxel data.  Its value must lie in the range
       * [0..N), where N is the number of columns in parameter 'data'.
       * @param startV This parameter specifies the starting V coordinate
       * (row) in the voxel data.  Its value must lie in the range
       * [0..M), where M is the number of rows in parameter 'data'.
       * @param startW This parameter specifies the starting W coordinate
       * (row) in the voxel data.  Its value must lie in the range
       * [0..M), where M is the number of rows in parameter 'data'.
       * @param stepU This parameter specifies the increment by which
       * the U coordinate changes as we move along the direction of the
       * ray.  It must be either 1 or -1.
       * @param stepV This parameter specifies the increment by which
       * the V coordinate changes as we move along the direction of the
       * ray.  It must be either 1 or -1.
       * @param stepW This parameter specifies the increment by which
       * the W coordinate changes as we move along the direction of the
       * ray.  It must be either 1 or -1.
       * @param tMaxU This parameter specifies the value of ray
       * parameter 't' at which the ray passes from the current column
       * into the next column.  Parameter 't' is described in the
       * documentation for class AmanatidesWoo3D.
       * @param tMaxV This parameter specifies the value of ray
       * parameter 't' at which the ray passes from the current row into
       * the next row.  Parameter 't' is described in the documentation
       * for class AmanatidesWoo3D.
       * @param tMaxW This parameter specifies the value of ray
       * parameter 't' at which the ray passes from the current slice into
       * the next slice.  Parameter 't' is described in the documentation
       * for class AmanatidesWoo3D.
       * @param tDeltaU This parameter specifies the increment to ray
       * parameter 't' which moves one exactly one column width to the
       * left or right, where left and right describe the directions of
       * the negative and positive U axis, respectively.  Parameter 't'
       * is described in the documentation for class AmanatidesWoo3D.
       * @param tDeltaV This parameter specifies the increment to ray
       * parameter 't' which moves one exactly one row width up or down,
       * where up and down describe the directions of the negative and
       * positive V axis, respectively.  Parameter 't' is described in
       * the documentation for class AmanatidesWoo3D.
       * @param tDeltaW This parameter specifies the increment to ray
       * parameter 't' which moves one exactly one slice width up or down,
       * where up and down describe the directions of the negative and
       * positive Z axis, respectively.  Parameter 't' is described in
       * the documentation for class AmanatidesWoo3D.
       * @param tStart This parameter specifies the value of ray
       * parameter 't' at the very beginning point of the iteration.
       */
      AmanatidesWoo3DIterator(ARRAY3D& data,
                              INT_TYPE startU, INT_TYPE startV, INT_TYPE startW,
                              INT_TYPE stepU, INT_TYPE stepV, INT_TYPE stepW,
                              FLOAT_TYPE tMaxU, FLOAT_TYPE tMaxV, FLOAT_TYPE tMaxW,
                              FLOAT_TYPE tDeltaU, FLOAT_TYPE tDeltaV, FLOAT_TYPE tDeltaW,
                              FLOAT_TYPE tStart);

      /** 
       * Copy constructor.
       * 
       * @param source This argument specifies the AmanatidesWoo3D
       * instance to be copied.
       */
      AmanatidesWoo3DIterator(const AmanatidesWoo3DIterator& source);

      /** 
       * Destructor.
       */
      ~AmanatidesWoo3DIterator() {};

      /** 
       * This method returns the ray parameter t at which the ray being
       * followed passes into the current voxel.  In other words, the
       * value t such that (rayOrigin + t * rayDirection) is the point
       * of entry into the current voxel.
       * 
       * @return The return value is the value of t at which the ray
       * passes into the current voxel.
       */
      FLOAT_TYPE
      tEntry() {return m_tEntry;}

      /** 
       * This method returns the ray parameter t at which the ray being
       * followed passes out of the current voxel.  In other words, the
       * value t such that (rayOrigin + t * rayDirection) is the point
       * of exit from the current voxel.  Invoking this method carries a
       * computational cost of 1 FLOAT_TYPE comparison.
       * 
       * @return The return value is the value of t at which the ray
       * passes out of the current voxel.
       */
      FLOAT_TYPE
      tExit() {return std::min(m_tMaxU, std::min(m_tMaxV, m_tMaxW));}

      /** 
       * This method returns the U coordinate of the current voxel.
       * 
       * @return The return value is the U coordinate of the current
       * voxel.
       */
      INT_TYPE
      U() {return m_U;}

      /** 
       * This method returns the V coordinate of the current voxel.
       * 
       * @return The return value is the V coordinate of the current
       * voxel.
       */
      INT_TYPE
      V() {return m_V;}

      /** 
       * This method returns the W coordinate of the current voxel.
       * 
       * @return The return value is the W coordinate of the current
       * voxel.
       */
      INT_TYPE
      W() {return m_W;}

      /** 
       * This operator returns a reference to the Array3D element at the
       * current voxel.  With each increment of the
       * AmanatidesWoo3DIterator instance, this operator will return a
       * reference to the next voxel along the ray.
       * 
       * @return The return value is a reference the the relevant
       * Array3D element.
       */
      inline typename ARRAY3D::value_type& // element_type?
      operator*();

      /** 
       * This operator returns a pointer to the Array3D element at the
       * current voxel.  With each increment of the
       * AmanatidesWoo3DIterator instance, this operator will return a
       * pointer to the next voxel along the ray.
       * 
       * @return The return value is a pointer the the relevant
       * Array3D element.
       */
      inline typename ARRAY3D::value_type* // element_type?
      operator->();

      /** 
       * The pre-increment operator increments the iterator so that it
       * points to the next voxel along the path.
       * 
       * @return The return value is a reference to *this.
       */
      inline AmanatidesWoo3DIterator&
      operator++();	             // prefix

      /** 
       * The post-increment operator increments the iterator so that it
       * points to the next voxel along the path.  It differs from the
       * pre-increment operator in its return value.  Traditionally,
       * post-increment is a little slower than pre-increment.
       * 
       * @param dummy This parameter is a dummy which indicates to the
       * compiler that this operation is post-increment (rather than
       * pre-increment).
       *
       * @return The return value is a copy of *this which was generated
       * before the increment.
       */
      inline AmanatidesWoo3DIterator
      operator++(int dummy);                 // postfix

      /** 
       * This is the assignment operator.  It copies the value of its
       * argument into *this.
       * 
       * @param source This argument specifies the
       * AmanatidesWoo3DIterator instance to be copied.
       * @return The return value is a reference to *this.
       */
      AmanatidesWoo3DIterator&
      operator=(const AmanatidesWoo3DIterator& source);

      /** 
       * The equality operator returns true if both the argument and
       * *this currently reference a valid voxel, or if both the
       * argument and *this currently reference an invalid voxel.  In all
       * other cases the return is false.
       *
       * NOTE: This behavior is not exactly what you'd expect for an
       * equality operator.  references the same voxel as the argument.
       * 
       * @param other This argument is a second AmanatidesWoo3DIterator
       * instance that is to be compared with *this.
       *
       * @return The return value is true if *this and other both
       * reference in-bounds voxels, or if *this and other both
       * reference out-of-bounds voxels, false otherwise.
       */
      inline bool
      operator==(const AmanatidesWoo3DIterator& other);

      /** 
       * The equality operator returns false if both the argument and
       * *this currently reference a valid voxel, or if both the
       * argument and *this currently reference an invalid voxel.  In all
       * other cases the return is true.
       *
       * NOTE: This behavior is not exactly what you'd expect for an
       * equality operator.  references the same voxel as the argument.
       * 
       * @param other This argument is a second AmanatidesWoo3DIterator
       * instance that is to be compared with *this.
       *
       * @return The return value is false if *this and other both
       * reference in-bounds voxels, or if *this and other both
       * reference out-of-bounds voxels, true otherwise.
       */
      inline bool
      operator!=(const AmanatidesWoo3DIterator& other);

    private:
      ARRAY3D& m_data;
      bool m_inBounds;
      INT_TYPE m_stepU;
      INT_TYPE m_stepV;
      INT_TYPE m_stepW;
      FLOAT_TYPE m_tDeltaU;
      FLOAT_TYPE m_tDeltaV;
      FLOAT_TYPE m_tDeltaW;
      FLOAT_TYPE m_tEntry;
      FLOAT_TYPE m_tMaxU;
      FLOAT_TYPE m_tMaxV;
      FLOAT_TYPE m_tMaxW;
      INT_TYPE m_U;
      INT_TYPE m_uLimit;
      INT_TYPE m_V;
      INT_TYPE m_vLimit;
      INT_TYPE m_W;
      INT_TYPE m_wLimit;
    };

  } // namespace numeric

} // namespace brick

// Include file containing definitions of inline and template
// functions.
#include <brick/numeric/amanatidesWoo3DIterator_impl.hh>

#endif /* #ifndef BRICK_NUMERIC_AMANATIDESWOO3DITERATOR_HH */
