/**
***************************************************************************
* @file brick/numeric/subArray1D.hh
*
* Header file declaring SubArray1D class template.
*
* Copyright (C) 2001-2011 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/


#ifndef BRICK_NUMERIC_SUBARRAY1D_HH
#define BRICK_NUMERIC_SUBARRAY1D_HH

#include <brick/numeric/array1D.hh>
#include <brick/numeric/slice.hh>

namespace brick {

  namespace numeric {

    /**
     ** Header file defining a simple SubArray class to work with Array1D.h
     ** The goal here is simplicity.  This is not intended to be a full
     ** interface.  Just enough to let you select and copy subArrays.
     ** Use this class via the function subArray(), like this:
     **
     **  subArray(array0) = subArray(array1, Slice(0, 6, 2));
     **
     ** Note that this class has deep copy semantics.
     **/
    template <class Type>
    class SubArray1D {
    public:

      /**
       * The single argument constructs a subarray referencing every
       * element of the source array.
       *
       * @param source This argument specifies the Array1D into which to index.
       */
      SubArray1D(const Array1D<Type>& source);


      /**
       * This constructor permits slicing of the source array.  The
       * resulting subarray references only those elements of the source
       * array indicated by argument slice0.
       *
       * @param source This argument specifies the Array1D into which to
       * index.
       *
       * @param slice0 This argument specifies which elements to include
       * in the subArray.
       */
      SubArray1D(const Array1D<Type>& source, const Slice& slice0);


      /**
       * This is the copy constructor.  After construction, *this will
       * reference the same elements (of the same source array) as the
       * copied SubArray1D instance.
       *
       * @param other This argument specifies the SubArray1D instance to
       * be copied
       */
      SubArray1D(const SubArray1D<Type> &other);


      /**
       * Destroys the SubArray1D instance.
       */
      virtual
      ~SubArray1D() {}


      /**
       * This conversion operator generates an Array1D instance from a
       * SubArray1D.  The element values from the *this are copied into
       * the newly created Array1D instance.
       *
       * @return An Array1D instance containing copies of the elements
       * referenced by *this.
       */
      operator Array1D<Type>() const;


      /**
       * This member function returns the number of elements referenced by
       * *this.  For example, if a SubArray1D references every third element
       * of an 18 element Array1D, then its size() method will return 6.
       *
       * @return The number of elements referenced by *this.
       */
      inline size_t
      size() const {return this->m_size;}


      /**
       * This member function returns the index of the first element
       * referenced by *this.  For example, if a SubArray1D references
       * every third element of an 18 element Array1D, starting from
       * element 4, then its start() method will return 4.
       *
       * @return The index of the first element referenced by *this.
       */
      inline size_t
      start() const {return this->m_start;}


      /**
       * This member function returns the spacing of the elements
       * referenced by *this.  For example, if a SubArray1D references
       * every third element of an Array1D, then its stride() method
       * will return 3.
       *
       * @return The spacing of the elements referenced by *this.
       */
      inline size_t
      stride() const {return this->m_stride;}


      /**
       * This method returns an Array1D instance referencing the same
       * memory as the Array1D instance from which *this was
       * constructed.
       *
       * @return An Array1D instance which references the same memory as
       * the Array1D instance from which *this was constructed.
       */
      inline Array1D<Type>
      getArray() const {return this->m_source;}


      /**
       * This assignment operator performs a deep copy, copying each
       * element from other into *this.  Note that this modifies the
       * corresponding elements of the array from which *this was
       * created.
       *
       * @param other This argument specifies the SubArray1D instance to
       * be copied
       *
       * @return A reference to *this.
       */
      SubArray1D<Type>&
      operator=(const SubArray1D<Type>& other);


      /**
       * This assignment operator sets each element reference by *this
       * to the specified value.  Note that this modifies the
       * corresponding elements of the array from which *this was
       * created.
       *
       * @param value This argument specifies the value to which the
       * array elements should be set.
       *
       * @return A reference to *this.
       */
      SubArray1D<Type>&
      operator=(Type value);

      /**
       * This operator increments each element of *this by the
       * corresponding element of its argument.  Note that this modifies
       * the corresponding elements of the array from which *this was
       * created.
       *
       * @param other This argument specifies a SubArray1D instance, the
       * elements of which will be added to the elements of *this.
       *
       * @return A reference to *this.
       */
      SubArray1D<Type>&
      operator+=(const SubArray1D<Type>& other);

      /**
       * This operator decrements each element of *this by the
       * corresponding element of its argument.  Note that this modifies
       * the corresponding elements of the array from which *this was
       * created.
       *
       * @param other This argument specifies a SubArray1D instance, the
       * elements of which will be subtracted from the elements of *this.
       *
       * @return A reference to *this.
       */
      SubArray1D<Type>&
      operator-=(const SubArray1D<Type>& other);

    private:
      inline void
      checkArray1DSize(const Array1D<Type>& other) const;

      inline void
      checkSubArray1DSize(const SubArray1D<Type>& other) const;

      SubArray1D<Type>&
      copy(const SubArray1D<Type>& other);

      Array1D<Type> m_source;
      int m_start;
      int m_stop;
      int m_stride;
      size_t m_size;
    };


    /* Non-member functions */

    /**
     * This is a convenience function for constructing SubArray1D
     * instances which reference every element of the source array.  Use
     * this function when you want to modify every element of an Array1D
     * instance.
     *
     * @param source This argument specifies the Array1D into which to index.
     *
     * @return A SubArray1D instance referencing every element of the
     * source array.
     */
    template <class Type>
    inline SubArray1D<Type>
    subArray(const Array1D<Type>& source) {return SubArray1D<Type>(source);}


    /**
     * This is a convenience function for constructing SubArray1D
     * instances which reference only selected elements of the source
     * array.  Use this function when you want to modify only some of
     * the elements of an Array1D instance.  For example, you might use
     * the following to copy the first, sixth, eleventh, and sixteenth
     * elements of array1 into the second, sixth, tenth, and fourteenth
     * elements of array0:
     *
     *  subArray(array0, Slice(1, 15, 4)) = subArray(array1, Slice(0, 16, 5));
     *
     * @param source This argument specifies the Array1D into which to
     * index.
     *
     * @param rowSlice A slice instance indicating the target elements.
     *
     * @return A SubArray1D instance referencing the selected elements
     * of the source array.
     */
    template <class Type>
    inline SubArray1D<Type>
    subArray(const Array1D<Type>& source,
             const Slice& rowSlice) {
      return SubArray1D<Type>(source, rowSlice);
    }


    /**
     * This stream output operator sends a text representation of the
     * SubArray1D instance to the supplied stream instance.
     *
     * @param stream The stream to which data should be written.
     *
     * @param subArray0 The SubArray1D instance to be written.
     *
     * @return A reference to argument stream.
     */
    template <class Type>
    std::ostream&
    operator<<(std::ostream& stream, const SubArray1D<Type>& subArray0);

  } // namespace numeric

} // namespace brick

#include <brick/numeric/subArray1D_impl.hh>

#endif /* #ifdef BRICK_NUMERIC_SUBARRAY1D_HH */
