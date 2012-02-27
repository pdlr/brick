/**
***************************************************************************
* @file brick/numeric/amanatidesWoo2DIterator.hh
*
* Header file declaring AmanatidesWoo2DIterator class.
*
* Copyright (C) 2004-2007,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_AMANATIDESWOO2DITERATOR_HH
#define BRICK_NUMERIC_AMANATIDESWOO2DITERATOR_HH

#include <iostream>
#include <brick/common/exception.hh>

namespace brick {

  namespace numeric {
    
    /**
     ** This class provides access to the elements of a data array along
     ** a straight path, and does the actual work of Amanatides and
     ** Woo's fast voxel traversal algorithm.  Typically, an
     ** AmanatidesWoo2DIterator instance will be created by an
     ** AmanatidesWoo2D object in order to access a line of pixels
     ** specified through the AmanatidesWoo2D class interface.  The user
     ** will probably never need to directly construct an
     ** AmanatidesWoo2DIterator.  For more information on the fast voxel
     ** traversal algorithm of Amanatides and Woo, please refer to
     ** [ref].
     **/
    template <class ARRAY2D, class FLOAT_TYPE = double, class INT_TYPE = int>
    class AmanatidesWoo2DIterator
      : public std::iterator<std::forward_iterator_tag,
                             typename ARRAY2D::value_type>
    {
    public:
      /** 
       * The class constructor is initialized with all of the internal
       * variables of the voxel traversal algorithm.
       * 
       * @param data This parameter is a reference to the 2D data over
       * which to iterate.
       * @param startU This parameter specifies the starting U coordinate
       * (column) in the pixel data.  Its value must lie in the range
       * [0..N), where N is the number of columns in parameter 'data'.
       * @param startV This parameter specifies the starting V coordinate
       * (row) in the pixel data.  Its value must lie in the range
       * [0..M), where M is the number of rows in parameter 'data'.
       * @param stepU This parameter specifies the increment by which
       * the U coordinate changes as we move along the direction of the
       * ray.  It must be either 1 or -1.
       * @param stepV This parameter specifies the increment by which
       * the V coordinate changes as we move along the direction of the
       * ray.  It must be either 1 or -1.
       * @param tMaxU This parameter specifies the value of ray
       * parameter 't' at which the ray passes from the current column
       * into the next column.  Parameter 't' is described in the
       * documentation for class AmanatidesWoo2D.
       * @param tMaxV This parameter specifies the value of ray
       * parameter 't' at which the ray passes from the current row into
       * the next row.  Parameter 't' is described in the documentation
       * for class AmanatidesWoo2D.
       * @param tDeltaU This parameter specifies the increment to ray
       * parameter 't' which moves one exactly one column width to the
       * left or right, where left and right describe the directions of
       * the negative and positive U axis, respectively.  Parameter 't'
       * is described in the documentation for class AmanatidesWoo2D.
       * @param tDeltaV This parameter specifies the increment to ray
       * parameter 't' which moves one exactly one row width up or down,
       * where up and down describe the directions of the negative and
       * positive V axis, respectively.  Parameter 't' is described in
       * the documentation for class AmanatidesWoo2D.
       * @param tStart This parameter specifies the value of ray
       * parameter 't' at the very beginning point of the iteration.
       */
      AmanatidesWoo2DIterator(ARRAY2D& data,
                              INT_TYPE startU, INT_TYPE startV,
                              INT_TYPE stepU, INT_TYPE stepV,
                              FLOAT_TYPE tMaxU, FLOAT_TYPE tMaxV,
                              FLOAT_TYPE tDeltaU, FLOAT_TYPE tDeltaV,
                              FLOAT_TYPE tStart);

      /** 
       * Copy constructor.
       * 
       * @param source This argument specifies the AmanatidesWoo2D
       * instance to be copied.
       */
      AmanatidesWoo2DIterator(const AmanatidesWoo2DIterator& source);

      /** 
       * Destructor.
       */
      ~AmanatidesWoo2DIterator() {};

      /** 
       * This method returns the ray parameter t at which the ray being
       * followed passes into the current pixel.  In other words, the
       * value t such that (rayOrigin + t * rayDirection) is the point
       * of entry into the current pixel.
       * 
       * @return The return value is the value of t at which the ray
       * passes into the current pixel.
       */
      FLOAT_TYPE
      tEntry() {return m_tEntry;}

      /** 
       * This method returns the ray parameter t at which the ray being
       * followed passes out of the current pixel.  In other words, the
       * value t such that (rayOrigin + t * rayDirection) is the point
       * of exit from the current pixel.  Invoking this method carries a
       * computational cost of 1 FLOAT_TYPE comparison.
       * 
       * @return The return value is the value of t at which the ray
       * passes out of the current pixel.
       */
      FLOAT_TYPE
      tExit() {return std::min(m_tMaxU, m_tMaxV);}

      /** 
       * This method returns the U coordinate of the current pixel.
       * 
       * @return The return value is the U coordinate of the current
       * pixel.
       */
      INT_TYPE
      U() {return m_U;}

      /** 
       * This method returns the V coordinate of the current pixel.
       * 
       * @return The return value is the V coordinate of the current
       * pixel.
       */
      INT_TYPE
      V() {return m_V;}

      /** 
       * This operator returns a reference to the Array2D element at the
       * current pixel.  With each increment of the
       * AmanatidesWoo2DIterator instance, this operator will return a
       * reference to the next pixel along the ray.
       * 
       * @return The return value is a reference the the relevant
       * Array2D element.
       */
      inline typename ARRAY2D::value_type& // element_type?
      operator*();

      /** 
       * This operator returns a pointer to the Array2D element at the
       * current pixel.  With each increment of the
       * AmanatidesWoo2DIterator instance, this operator will return a
       * pointer to the next pixel along the ray.
       * 
       * @return The return value is a pointer the the relevant
       * Array2D element.
       */
      inline typename ARRAY2D::value_type* // element_type?
      operator->();

      /** 
       * The pre-increment operator increments the iterator so that it
       * points to the next pixel along the path.
       * 
       * @return The return value is a reference to *this.
       */
      inline AmanatidesWoo2DIterator&
      operator++();	             // prefix

      /** 
       * The post-increment operator increments the iterator so that it
       * points to the next pixel along the path.  It differs from the
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
      inline AmanatidesWoo2DIterator
      operator++(int dummy);                 // postfix

      /** 
       * This is the assignment operator.  It copies the value of its
       * argument into *this.
       * 
       * @param source This argument specifies the
       * AmanatidesWoo2DIterator instance to be copied.
       * @return The return value is a reference to *this.
       */
      AmanatidesWoo2DIterator&
      operator=(const AmanatidesWoo2DIterator& source);

      /** 
       * The equality operator returns true if both the argument and
       * *this currently reference a valid pixel, or if both the
       * argument and *this currently reference an invalid pixel.  In all
       * other cases the return is false.
       *
       * NOTE: This behavior is not exactly what you'd expect for an
       * equality operator.  references the same pixel as the argument.
       * 
       * @param other This argument is a second AmanatidesWoo2DIterator
       * instance that is to be compared with *this.
       *
       * @return The return value is true if *this and other both
       * reference in-bounds pixels, or if *this and other both
       * reference out-of-bounds pixels, false otherwise.
       */
      inline bool
      operator==(const AmanatidesWoo2DIterator& other);

      /** 
       * The equality operator returns false if both the argument and
       * *this currently reference a valid pixel, or if both the
       * argument and *this currently reference an invalid pixel.  In all
       * other cases the return is true.
       *
       * NOTE: This behavior is not exactly what you'd expect for an
       * equality operator.  references the same pixel as the argument.
       * 
       * @param other This argument is a second AmanatidesWoo2DIterator
       * instance that is to be compared with *this.
       *
       * @return The return value is false if *this and other both
       * reference in-bounds pixels, or if *this and other both
       * reference out-of-bounds pixels, true otherwise.
       */
      inline bool
      operator!=(const AmanatidesWoo2DIterator& other);

    private:
      ARRAY2D& m_data;
      bool m_inBounds;
      INT_TYPE m_stepU;
      INT_TYPE m_stepV;
      FLOAT_TYPE m_tDeltaU;
      FLOAT_TYPE m_tDeltaV;
      FLOAT_TYPE m_tEntry;
      FLOAT_TYPE m_tMaxU;
      FLOAT_TYPE m_tMaxV;
      INT_TYPE m_U;
      INT_TYPE m_uLimit;
      INT_TYPE m_V;
      INT_TYPE m_vLimit;
    };

  } // namespace numeric

} // namespace brick

// Include file containing definitions of inline and template
// functions.
#include <brick/numeric/amanatidesWoo2DIterator_impl.hh>

#endif /* #ifndef BRICK_NUMERIC_AMANATIDESWOO2DITERATOR_HH */
