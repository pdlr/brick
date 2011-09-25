/**
***************************************************************************
* @file brick/numeric/subArray2D.hh 
* Header file declaring SubArray2D class template.
*
* Copyright (C) 2001-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/


#ifndef BRICK_NUMERIC_SUBARRAY2D_HH
#define BRICK_NUMERIC_SUBARRAY2D_HH

#include <brick/numeric/array2D.hh>
#include <brick/numeric/slice.hh>

namespace brick {

  namespace numeric {
    
    /**
     ** Header file defining a simple SubArray class to work with Array2D.h
     ** The goal here is simplicity.  This is not intended to be a full
     ** interface.  Just enough to let you select and copy subArrays.
     ** Use this class via the function subArray(), like this:
     **
     **  subArray(array0) = subArray(array1, Slice(0, 6, 2), Slice(1, 4));
     **
     ** or
     **
     **  Array2D<double> array0 = subArray(array1, Slice(0, 6, 2), Slice(1, 4));
     **
     ** Note that this class has deep copy semantics.
     **/
    template <class Type>
    class SubArray2D {
    public:
      /** 
       * The single argument constructs a subarray referencing every
       * element of the source array.
       * 
       * @param source This argument specifies the Array2D into which to index.
       */
      SubArray2D(const Array2D<Type>& source);

      /** 
       * This constructor permits slicing of the source array.  The
       * resulting subarray references only those elements of the source
       * array which lie in rows indicated by rowSlice and also in
       * columns indicated by columnSlice.
       * 
       * @param source This argument specifies the Array2D into which to
       * index.
       *
       * @param rowSlice This argument specifies which rows to include
       * in the subArray.
       *
       * @param columnSlice This argument specifies which columns to
       * include in the subArray.
       */
      SubArray2D(const Array2D<Type>& source, const Slice& rowSlice,
                 const Slice& columnSlice);

      /** 
       * This is the copy constructor.  After construction, *this will
       * reference the same elements (of the same source array) as the
       * copied SubArray2D instance.
       * 
       * @param other This argument specifies the SubArray2D instance to
       * be copied
       */
      SubArray2D(const SubArray2D<Type> &other);

      /** 
       * Destroys the SubArray2D instance.
       */
      virtual
      ~SubArray2D() {}

      /** 
       * This conversion operator generates an Array2D instance from a
       * SubArray2D.  The element values from the *this are copied into
       * the newly created Array2D instance.
       * 
       * @return An Array2D instance containing copies of the elements
       * referenced by *this.
       */
      operator Array2D<Type>() const;

      /** 
       * This member function returns the number of rows referenced by
       * *this.  For example, if a SubArray2D references every third row
       * of an 18 row Array2D, then its rows() method will return 6.
       * 
       * @return The number of rows referenced by *this.
       */
      inline size_t
      rows() const {return this->m_rows;}

      /** 
       * This member function returns the number of columns referenced
       * by *this.  For example, if a SubArray2D references every third
       * column of an 18 column Array2D, then its columns() method will
       * return 6.
       * 
       * @return The number of columns referenced by *this.
       */
      inline size_t
      columns() const {return this->m_columns;}

      /** 
       * This member function returns the index of the first row
       * referenced by *this.  For example, if a SubArray2D references
       * every third row of an 18 row Array2D, starting from row 4, then
       * its startRow() method will return 4.
       * 
       * @return The index of the first row referenced by *this.
       */
      inline size_t
      startRow() const {return this->m_startRow;}

      /** 
       * This member function returns the index of the first column
       * referenced by *this.  For example, if a SubArray2D references
       * every third column of an 18 column Array2D, starting from
       * column 4, then its startColumn() method will return 4.
       * 
       * @return The index of the first column referenced by *this.
       */
      inline size_t
      startColumn() const {return this->m_startColumn;}
    
      /** 
       * This member function returns the spacing of the rows referenced
       * by *this.  For example, if a SubArray2D references every third
       * row of an Array2D, then its rowStride() method will return 3.
       * 
       * @return The spacing of the rows referenced by *this.
       */
      inline size_t
      rowStride() const {return this->m_rowStride;}

      /** 
       * This member function returns the spacing of the columns
       * referenced by *this.  For example, if a SubArray2D references
       * every third column of an Array2D, then its columnStride()
       * method will return 3.
       * 
       * @return The spacing of the columns referenced by *this.
       */
      inline size_t
      columnStride() const {return this->m_columnStride;}

      /** 
       * This method returns an Array2D instance referencing the same
       * memory as the Array2D instance from which *this was
       * constructed.
       * 
       * @return An Array2D instance which references the same memory as
       * the Array2D instance from which *this was constructed.
       */
      inline Array2D<Type>
      getArray() const {return this->m_source;}
  
      /** 
       * This assignment operator performs a deep copy, copying each
       * element from other into *this.  Note that this modifies the
       * corresponding elements of the array from which *this was
       * created.
       * 
       * @param other This argument specifies the SubArray2D instance to
       * be copied
       * @return A reference to *this.
       */
      SubArray2D<Type>&
      operator=(const SubArray2D<Type>& other);

      /** 
       * This assignment operator sets each element reference by *this
       * to the specified value.  Note that this modifies the
       * corresponding elements of the array from which *this was
       * created.
       * 
       * @param value This argument specifies the value to which the
       * array elements should be set.
       * @return A reference to *this.
       */
      SubArray2D<Type>& operator=(Type value);

      /** 
       * This operator increments each element of *this by the
       * corresponding element of its argument.  Note that this modifies
       * the corresponding elements of the array from which *this was
       * created.
       * 
       * @param other This argument specifies a SubArray2D instance, the
       * elements of which will be added to the elements of *this.
       * @return A reference to *this.
       */
      SubArray2D<Type>&
      operator+=(const SubArray2D<Type>& other);

      /** 
       * This operator decrements each element of *this by the
       * corresponding element of its argument.  Note that this modifies
       * the corresponding elements of the array from which *this was
       * created.
       * 
       * @param other This argument specifies a SubArray2D instance, the
       * elements of which will be subtracted from the elements of *this.
       * @return A reference to *this.
       */
      SubArray2D<Type>&
      operator-=(const SubArray2D<Type>& other);

    private:

      inline void
      checkArray2DSize(const Array2D<Type>& other) const;

      inline void
      checkSubArray2DSize(const SubArray2D<Type>& other) const;

      SubArray2D<Type>&
      copyColumnMajor(const SubArray2D<Type>& other);

      SubArray2D<Type>&
      copyRowMajor(const SubArray2D<Type>& other);

      Array2D<Type> m_source;
      int m_startRow;
      int m_stopRow;
      int m_rowStride;
      /// number of rows
      size_t m_rows;
      int m_startColumn;
      int m_stopColumn;
      int m_columnStride;
      /// number of columns
      size_t m_columns;
    };


    /* Non-member functions */

    /** 
     * This is a convenience function for constructing SubArray2D
     * instances which reference every element of the source array.  Use
     * this function when you want to modify every element of an Array2D
     * instance.  For example, you might use the following to copy the 9
     * elements at the intersections of rows [0, 2, 4] and columns [1,
     * 2, 3] of array1 into an appropriately sized array0:
     *
     *  subArray(array0) = subArray(array1, Slice(0, 6, 2), Slice(1, 4));
     * 
     * @param source This argument specifies the Array2D into which to index.
     * @return A SubArray2D instance referencing every element of the
     * source array.
     */
    template <class Type>
    inline SubArray2D<Type>
    subArray(const Array2D<Type>& source) {return SubArray2D<Type>(source);}

  
    /** 
     * This is a convenience function for constructing SubArray2D
     * instances which reference only selected elements of the source
     * array.  Use this function when you want to modify only some of
     * the elements of an Array2D instance.  For example, you might use
     * the following to copy the 9 elements at the intersections of rows
     * [0, 2, 4] and columns [1, 2, 3] of array1 into the elements at
     * the intersections of rows [2, 4, 6] and columns [3, 6, 9] of
     * array0:
     *
     *  subArray(array0, Slice(2, 8, 2), Slice(3, 12, 3)
     *    = subArray(array1, Slice(0, 6, 2), Slice(1, 4));
     * 
     * @param source This argument specifies the Array2D into which to index.
     * @param rowSlice A slice instance indicating the target rows.
     * @param columnSlice A slice instance indicating the target columns.
     * @return A SubArray2D instance referencing the elements of the
     * source array which lie at the intersections of the selected rows
     * and columns.
     */
    template <class Type>
    inline SubArray2D<Type>
    subArray(const Array2D<Type>& source,
             const Slice& rowSlice,
             const Slice& columnSlice) {
      return SubArray2D<Type>(source, rowSlice, columnSlice);
    }

  
    /** 
     * This function is just like subArray(const Array2D<Type>&, const
     * Slice&, const Slice&), above except that it constructs a
     * SubArray2D instance which references only selected elements _of a
     * particular row_ in the source array.  For example, you might use
     * the following to copy the 3 elements at the intersections of row
     * 2 and columns [1, 2, 3] of array1 into the elements at the
     * intersections of rows 4 and columns [3, 6, 9] of array0:
     *
     *  subArray(array0, 6, Slice(3, 12, 3)
     *    = subArray(array1, 2, Slice(1, 4));
     * 
     * @param source This argument specifies the Array2D into which to index.
     * @param row An index indicating the target row.
     * @param columnSlice A slice instance indicating the target columns.
     * @return A SubArray2D instance referencing the elements of the
     * source array which lie at the intersections of the target row
     * and the selected columns.
     */
    template <class Type>
    inline SubArray2D<Type>
    subArray(const Array2D<Type>& source,
             const int row,
             const Slice& columnSlice) {
      return SubArray2D<Type>(source, Slice(row, row+1), columnSlice);
    }

  
    /** 
     * This function is just like subArray(const Array2D<Type>&, const
     * Slice&, const Slice&), above except that it constructs a
     * SubArray2D instance which references only selected elements _of a
     * particular column_ in the source array.  For example, you might
     * use the following to copy the 3 elements at the intersections of
     * rows [0, 2, 4] and column 2 of array1 into the elements at the
     * intersections of rows [2, 4, 6] and column 6 array0:
     *
     *  subArray(array0, Slice(2, 8, 2), 6)
     *    = subArray(array1, Slice(0, 6, 2), 2);
     * 
     * @param source This argument specifies the Array2D into which to index.
     * @param rowSlice A slice instance indicating the target rows.
     * @param column An index indicating the target column.
     * @return A SubArray2D instance referencing the elements of the
     * source array which lie at the intersections of the selected rows
     * and the target column.
     */
    template <class Type>
    inline SubArray2D<Type>
    subArray(const Array2D<Type>& source,
             const Slice& rowSlice,
             int column) {
      return SubArray2D<Type>(source, rowSlice, Slice(column, column + 1));
    }

  
    /** 
     * This stream output operator sends a text representation of the
     * SubArray2D instance to the supplied stream instance.
     * 
     * @param stream The stream to which data should be written.
     * @param subArray0 The SubArray2D instance to be written.
     * @return A reference to argument stream.
     */
    template <class Type>
    std::ostream&
    operator<<(std::ostream& stream, const SubArray2D<Type>& subArray0);

  } // namespace numeric

} // namespace brick

#include <brick/numeric/subArray1D_impl.hh>

#endif /* #ifdef BRICK_NUMERIC_SUBARRAY2D_HH */
