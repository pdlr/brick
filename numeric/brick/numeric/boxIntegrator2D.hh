/**
***************************************************************************
* @file brick/numeric/boxIntegrator2D.hh
*
* Header file declaring the BoxIntegrator2D class.
*
* Copyright (C) 2006,2012 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_NUMERIC_BOXINTEGRATOR2D_HH
#define BRICK_NUMERIC_BOXINTEGRATOR2D_HH

#include <brick/numeric/array2D.hh>
#include <brick/numeric/index2D.hh>


namespace brick {

  namespace numeric {

    /**
     ** This class provides an efficient way integrate over
     ** rectangular regions of an Array2D instance.  Computational
     ** cost is incurred almost entirely in the constructor (or
     ** alternatively, the setArray() method), which has complexity
     ** O(N), where N is the number of elements in the 2D array.
     ** Member function getIntegral() returns the sum over a
     ** rectangular region of the array, and has complexity O(1).
     ** Template argument Type0 specifies the element type of the
     ** input array, while template argument Type1 specifies the
     ** output type, as well as the type used for internal sums.
     **/
    template <class Type0, class Type1>
    class BoxIntegrator2D {
    public:

      /**
       * This constructor performs almost no work, and simply
       * initializes the class instance to a "zero" state.
       */
      inline
      BoxIntegrator2D();


      /**
       * This constructor initializes the class instance, and then
       * passes its arguments to member function setArray() in order
       * to precompute integral information.
       *
       * @param inputArray This argument specifies the array over
       * which to integrate.  The data in the array is not retained,
       * so it is safe to de-allocate the array after the constructor
       * call has completed.
       */
      explicit
      BoxIntegrator2D(const Array2D<Type0>& inputArray);


      /**
       * This constructor initializes the class instance, and then
       * passes its arguments to member function setArray() in order
       * to precompute integral information.
       *
       * @param inputArray This argument specifies the array over
       * which to integrate.  The data in the array is not retained,
       * so it is safe to de-allocate the array after the constructor
       * call has completed.
       *
       * @param functor This single-argument functor will be applied
       * to each element of the array before integration.
       */
      template <class Functor>
      explicit
      BoxIntegrator2D(const Array2D<Type0>& inputArray, Functor functor);


      /**
       * This constructor works just like the single argument version
       * of setArray, with the exception that the pre-integration is
       * performed over only a rectangular sub-array.  This is useful
       * if you know you will only need to compute integrals over a
       * portion of the input array.
       *
       * @param inputArray This argument specifies the array over
       * which to integrate.  The data in the array is not retained,
       * so it is safe to de-allocate the array after the constructor
       * call has completed.
       *
       * @param corner0 This argument, along with corner1, defines the
       * sub-region over which to integrate.
       *
       * @param corner1 This argument, along with corner0, defines the
       * sub-region over which to integrate.
       */
      BoxIntegrator2D(const Array2D<Type0>& inputArray,
                      const Index2D& corner0,
                      const Index2D& corner1);


      /**
       * The copy constructor does a shallow copy of its argument,
       * however behavior of the class should be indistinguishable
       * from a deep copy (operations on each of the original and copy
       * will not affect the other).
       *
       * @param other This argument is the BoxIntegrator2D instance to
       * be copied.
       */
      BoxIntegrator2D(const BoxIntegrator2D& other);


      /**
       * This member function returns the integral over the
       * rectangular region with corner0 and corner1 at its diagonally
       * opposite corners.  It is very efficient, requiring only 4 2D
       * indexing operations and 3 add operations to compute the
       * result.
       *
       * For most intuitive operation, the elements of corner0 should
       * be either uniformly smaller or uniformly larger than the
       * elements of corner1.  If this is the case, then the returned
       * integral will include the elements located on the smaller of
       * corner0 and corner1, as well as on the two adjacent edges,
       * but will not include the element at the larger of corner0 and
       * corner1, and its two adjacent edges.
       *
       * If one element of corner0 is larger than the corresponding
       * element of corner1, while the other element is smaller than
       * the corresponding element of corner1, then the result will be
       * the same as returned by the following code snipped.
       *
       * @code
       * Index2D corner2(min(corner0.getU(), corner1.getU()),
       *                 min(corner0.getV(), corner1.getV()));
       * Index2D corner3(max(corner0.getU(), corner1.getU()),
       *                 max(corner0.getV(), corner1.getV()));
       * Type2 result = -1 * integrator.getIntegral(corner2, corner3);
       * @endcode
       *
       * @param corner0 This argument specifies one of the region
       * corners.  If the array was set using the single-argument
       * constructor or single-argument version of setArray, then
       * corner0 specifies the corner as an index into the entire
       * array.  If the array was set using the three-argument
       * constructor or three-argument version of setArray, then
       * corner0 specifies the corner relative to region over which
       * pre-integration was performed.
       *
       * @param corner1 This argument specifies the region corner
       * diagonally opposite to corner0.  If the array was set using
       * the single-argument constructor or single-argument version of
       * setArray, then corner1 specifies the corner as an index into
       * the entire array.  If the array was set using the
       * three-argument constructor or three-argument version of
       * setArray, then corner1 specifies the corner relative to
       * region over which pre-integration was performed.
       *
       * @return The return value is the integral over the specified
       * region.
       */
      Type1
      getIntegral(const Index2D& corner0, const Index2D& corner1);


      /**
       * This member function works just like member function
       * getIntegral(const Index2D&, const Index2D&), with two
       * exceptions: it is very slightly slower; and corner indexing
       * is slightly more sophisticated.  Please see the descriptions
       * of argments corner0 and corner1 for more information on this
       * revised indexing.  Please see the description of member
       * function getIntegral(const Index2D&, const Index2D&) for more
       * information on the integration.
       *
       * @param corner0 This argument specifies one of the region
       * corners.  The corner is specified as an index into the entire
       * array, regardless of whether the single-argument or
       * three-argument version of the constructor (or setArray()
       * method) was used.
       *
       * @param corner1 This argument specifies the region corner
       * diagonally opposit to corner0.  The corner is specified as an
       * index into the entire array, regardless of whether the
       * single-argument or three-argument version of the constructor
       * (or setArray() method) was used.
       *
       * @param dummy This argument is not used.  Its only role is to
       * differentiate this member function from the two argument
       * version of getIntegral().
       *
       * @return The return value is the integral over the specified
       * region.
       */
      Type1
      getIntegral(const Index2D& corner0, const Index2D& corner1, bool dummy);


      /**
       * This member function returns the raw 2D integral from which
       * box integration is performed.  It is normally used only by
       * other functions in the dlr_libs suite of libraries, but is a
       * stable part of the BoxIntegrator interface.
       *
       * @param row This argument specifies the row (relative to the
       * region over which pre-integration was performed) at which to
       * sample the raw integral.
       *
       * @param column This argument specifies the column (relative to
       * the region over which pre-integration was performed) at which
       * to sample the raw integral.
       *
       * @return The return value is the value of the raw pre-integral
       * at the specified location.
       */
      Type1
      getRawIntegral(size_t row, size_t column);


      /**
       * This member function is just like getRawIntegral(size_t,
       * size_t), but uses raster-order single indexing.
       *
       * @param index0 This argument specifies the location at which
       * to sample the raw integral.
       *
       * @return The return value is the value of the raw pre-integral
       * at the specified location.
       */
      Type1
      getRawIntegral(size_t index0);


      /**
       * This member function discards and previously cached integral
       * information, performs a double integral over the input array,
       * and caches the result for future use in computing
       * sub-integrals.
       *
       * @param inputArray This argument specifies the array over
       * which to integrate.  The data in the array is not retained,
       * so it is safe to de-allocate the array after the constructor
       * call has completed.
       */
      void
      setArray(const Array2D<Type0>& inputArray);


      /**
       * This member function discards and previously cached integral
       * information, applies the specified functor to each element of
       * the input array, performs a double integral over the input
       * array, and caches the result for future use in computing
       * sub-integrals.
       *
       * @param inputArray This argument specifies the array over
       * which to integrate.  The data in the array is not retained,
       * so it is safe to de-allocate the array after the constructor
       * call has completed.
       *
       * @param functor This single-argument functor will be applied
       * to each element of the array before integration.
       */
      template <class Functor>
      void
      setArray(const Array2D<Type0>& inputArray, Functor functor);


      /**
       * This member functionworks just like the single argument
       * version of setArray, with the exception that the
       * pre-integration is performed over only a rectangular
       * sub-array.  This is useful if you know you will only need to
       * compute integrals over a portion of the input array.
       *
       * @param inputArray This argument specifies the array over
       * which to integrate.  The data in the array is not retained,
       * so it is safe to de-allocate the array after the constructor
       * call has completed.
       *
       * @param corner0 This argument, along with corner1, defines the
       * sub-region over which to integrate.
       *
       * @param corner1 This argument, along with corner0, defines the
       * sub-region over which to integrate.
       */
      void
      setArray(const Array2D<Type0>& inputArray,
               const Index2D& corner0,
               const Index2D& corner1);


      /**
       * This member function works just like the three-argument
       * version of setArray, with the exception that the the
       * specified functor is applied to each element of the array
       * before integration.
       *
       * @param inputArray This argument specifies the array over
       * which to integrate.  The data in the array is not retained,
       * so it is safe to de-allocate the array after the constructor
       * call has completed.
       *
       * @param corner0 This argument, along with corner1, defines the
       * sub-region over which to integrate.
       *
       * @param corner1 This argument, along with corner0, defines the
       * sub-region over which to integrate.
       *
       * @param functor This single-argument functor will be applied
       * to each element of the array before integration.
       */
      template <class Functor>
      void
      setArray(const Array2D<Type0>& inputArray,
               const Index2D& corner0,
               const Index2D& corner1,
               Functor functor);


    protected:


      /**
       * This protected member function does the actual work of
       * pre-integrating the input array.
       *
       * @param inIter This argument is an iterator pointing to the
       * upper-left corner of the region to be integrated.
       *
       * @param roiRows This argument specifies how many rows there
       * are in the region to be integrated.
       *
       * @param roiColumns This argument specifies how many columns
       * there are in the region to be integrated.
       *
       * @param inputArrayColumns This argument how many columns there
       * are in the input array.
       */
      template <class Functor>
      void
      fillCache(typename Array2D<Type0>::const_iterator inIter,
                int roiRows,
                int roiColumns,
                int inputArrayColumns,
                Functor functor);


      Array2D<Type1> m_cache;
      Index2D m_corner0;
    };


  } // namespace numeric

} // namespace brick

// Include file containing definitions of inline and template
// functions.
#include <brick/numeric/boxIntegrator2D_impl.hh>

#endif /* #ifndef BRICK_NUMERIC_BOXINTEGRATOR2D_HH */
