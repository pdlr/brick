/**
***************************************************************************
* @file brick/sparse/array2D_impl.hh
*
* Header file defining a sparse matrix class template.
*
* Copyright (C) 2014 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_SPARSE_ARRAY2D_IMPL_HH
#define BRICK_SPARSE_ARRAY2D_IMPL_HH

// This file is included by array2D.hh, and should not be directly included
// by user code, so no need to include array2D.hh here.
//
// #include <brick/sparse/array2D.hh>

namespace brick {

  namespace sparse {

    // Default constructor initializes to zero size.
    template <class Type>
    Array2D<Type>::
    Array2D()
    {
      //xxx
    }


      /**
       * Constructs a "arrayRows x arrayColumns" element sparse array.
       *
       * @param arrayRows Number of rows in the array after successful
       * construction.
       *
       * @param arrayColumns Number of columns in the array after successful
       * construction.
       */
      Array2D(size_t arrayRows, size_t arrayColumns);


      /**
       * The copy constructor does a shallow copy.  The newly created
       * array points to the same data as copied array.
       *
       * @param source The Array2D instance to be copied.
       */
      Array2D(const Array2D<Type> &source);


      /**
       * Destroys the Array2D instance and deletes the internal data
       * store if no remaining arrays point to it.
       */
      virtual
      ~Array2D();


      /**
       * Reset the array to zero size, abandoning all contents.  This is
       * equivalent to this->reinit(0, 0);
       */
      void
      clear() {this->reinit(0, 0);}


      /**
       * Returns the array dimension along the "second" axis (the number
       * of columns in the array.)  This is synonymous with shape(1), but
       * may execute faster.
       *
       * @return Number of columns.
       */
      size_t
      columns() const {return this->getColumns();}


      /**
       * Optionally (if and only if BRICK_NUMERIC_CHECKBOUNDS is
       * defined) throw an exception if the shape of *this is
       * different than specified.
       *
       * @param arrayRows The required number of rows.
       *
       * @param arrayColumns The required number of columns.
       */
      inline void
      checkDimension(size_t arrayRows, size_t arrayColumns) const;


      /**
       * Allocates a new array and deep copies the contents of *this.
       *
       * @return A new array that is a (deep) copy of *this.
       */
      Array2D<Type>
      copy() const;


      /**
       * Deep copies the contents of source.  It is an error if source
       * does not have the same shape as *this.
       *
       * @param source The array to be copied.
       *
       * @exception ValueException thrown when array sizes do not match.
       */
      template <class Type2> void
      copy(const Array2D<Type2>& source);


      /**
       * This is an alias for member function isEmpty(), It returns
       * true if the array instance contains no elements.  It has
       * complexity O(1).
       *
       * @return The return value indicates whether or not the array is
       * empty.
       */
      bool
      empty() const {return this->isEmpty();}


      /**
       * Returns the array dimension along the "second" axis (the number
       * of columns in the array.)  This is synonymous with shape(1), but
       * may execute faster.
       *
       * @return Number of columns.
       */
      size_t
      getColumns() const {return m_columns;}


      /**
       * This member function returns a specific element of the array
       * by value.
       *
       * @param rowIndex This argument and the next specify which element
       * value should be returned.
       *
       * @param columnIndex This argument and the previous specify which
       * element value should be returned.
       *
       * @return The return value is a copy of the requested element.
       */
      Type
      getElement(size_t rowIndex, size_t columnIndex) const;


      /**
       * Returns the the internal ReferenceCount instance by
       * reference.  This member function is included to support
       * certain tests, and should generally not be used in client
       * code.
       *
       * @return A reference to the internal referenceCount instance.
       */
      common::ReferenceCount const&
      getReferenceCount() const {return m_referenceCount;}


      /**
       * Returns the array dimension along the "first" axis (the number
       * of rows in the array.)  This is synonymous with shape(0), but
       * may execute faster.
       *
       * @return Number of rows.
       */
      size_t
      getRows() const {return m_rows;}


      /**
       * Returns the number of elements in the array.  This is the
       * product of rows() and columns().
       *
       * @return Number of elements.
       */
      size_t
      getSize() const {return m_size;}


      /**
       * Returns true if the array instance contains no elements.  It
       * has complexity O(1).
       *
       * @return The return value indicates whether or not the array is
       * empty.
       */
      bool
      isEmpty() const {return this->getSize() == 0;}


      /**
       * This member function sets the value of the array from an input
       * stream.  The array is modified only if the read was successful,
       * otherwise the array is not modified, and failbit is set in the
       * stream state.  Because of this nice behavior, readFromStream is
       * quite slow.
       *
       * @param inputStream This is the stream from which to read the
       * array.
       *
       * @return The return value is a reference to inputStream.
       */
      std::istream&
      readFromStream(std::istream& inputStream);


      /**
       * Changes the shape of the array and reallocates storage.  The
       * current array contents are lost.  After a successful call to
       * reinit(), the array will be arrayRows x arrayColumns.
       *
       * @param arrayRows Requested row dimension.
       *
       * @param arrayColumns Requested column dimension.
       */
      void
      reinit(size_t arrayRows, size_t arrayColumns);


      /**
       * Behaves just like member function reinit() unless *this
       * already has the requested shape, in which case this member
       * function does nothing, or *this has the requested size (but a
       * different shape), in which case this member function behaves
       * just like member function resize().
       *
       * @param arrayRows Requested row dimension.
       *
       * @param arrayColumns Requested column dimension.
       */
      inline void
      reinitIfNecessary(size_t arrayRows, size_t arrayColumns);


      /**
       * This is an alias for member function getRows().  It returns
       * the array dimension along the "first" axis (the number of
       * rows in the array.)  This is synonymous with shape(0), but
       * may execute faster.
       *
       * @return Number of rows.
       */
      size_t
      rows() const {return this->getRows();}


      /**
       * This member function sets the value of a specific element of
       * the array.
       *
       * @param row This argument and the next specify which element
       * value should be set.
       *
       * @param column This argument and the previous specify which
       * element value should be set.
       *
       * @param value This argument will be copied into the selected
       * array element.
       *
       * @return The return value is a reference to the newly created
       * array element.
       */
      Type&
      setElement(size_t rowIndex, size_t columnIndex, const Type& value) {
        return this->operator()(rowIndex, columnIndex) = value;
      }


      /**
       * Returns a 2 element Array1D containing the dimensions of *this
       * along each axis.  That is, the first element of the returned
       * array will be this->rows(), and the second element will be
       * this->columns().
       *
       * @return Array describing the shape of *this.
       */
      Array1D<size_t>
      shape() const;


      /**
       * Returns the array dimension along the the axis indicated by
       * index.  This is synonymous with shape()[axis], but may execute
       * faster.
       *
       * @param axis Selected axis.
       *
       * @exception ValueException thrown on invalid axis parameter.
       *
       * @return Dimension along the selected axis.
       */
      size_t
      shape(size_t axis) const;


      /**
       * This is a synonym for member function getSize().  It returns
       * the number of elements in the array.  This is the product of
       * rows() and columns().
       *
       * @return Number of elements.
       */
      size_t
      size() const {return this->getSize();}


      /**
       * Assignment operator shallow copies the contents of source.
       * After the copy, both arrays reference the same data.
       *
       * @param source The Array2D instance to be copied.
       *
       * @return Reference to *this.
       */
      Array2D<Type>&
      operator=(const Array2D<Type>& source);


      /**
       * Returns a specific element of the array by reference.
       *
       * @param rowIndex Row of the selected element.
       *
       * @param columnIndex Column of the selected element.
       *
       * @return Reference to the array element.
       */
      Type&
      operator()(size_t rowIndex, size_t columnIndex) {
        this->checkBounds(rowIndex, columnIndex);
        return m_dataPtr[columnIndex + rowIndex * m_rowStep];
      }


      /**
       * Returns a specific element of the array by value.
       *
       * @param rowIndex Row of the selected element.
       *
       * @param column Column of the selected element.
       *
       * @return Value of the array element.
       */
      Type
      operator()(size_t rowIndex, size_t columnIndex) const {
        this->checkBounds(rowIndex, columnIndex);
        return m_dataPtr[columnIndex + rowIndex * m_rowStep];
      }


      /**
       * Increments each element of *this by the value of the
       * corresponding element of arg.
       *
       * @param arg Array2D of values to be added to the elements of
       * *this.
       *
       * @exception ValueException thrown when array sizes differ.
       *
       * @return Reference to *this.
       */
      template <class Type2> Array2D<Type>&
      operator+=(const Array2D<Type2>& arg);


      /**
       * Decrements each element of *this by the value of the
       * corresponding element of arg.
       *
       * @param arg Array2D of values to be subtracted from the elements
       * of *this.
       *
       * @exception ValueException thrown when array sizes differ.
       *
       * @return Reference to *this.
       */
      template <class Type2> Array2D<Type>&
      operator-=(const Array2D<Type2>& arg);


      /**
       * Multiplies each element of *this by the value of the
       * corresponding element of arg.
       *
       * @param arg Array2D of values by which the elements of *this are
       * to be multiplied.
       *
       * @exception ValueException thrown when array sizes differ.
       *
       * @return Reference to *this.
       */
      template <class Type2> Array2D<Type>&
      operator*=(const Array2D<Type2>& arg);


      /**
       * Divides each element of *this by the value of the
       * corresponding element of arg.
       *
       * @param arg Array2D of values by which the elements of *this are
       * to be divided.
       *
       * @exception ValueException thrown when array sizes differ.
       *
       * @return Reference to *this.
       */
      template <class Type2> Array2D<Type>&
      operator/=(const Array2D<Type2>& arg);


      /**
       * Increments each element of *this by a constant.
       *
       * @param arg Value by which array elements will be incremented.
       *
       * @return Reference to *this.
       */
      Array2D<Type>&
      operator+=(Type arg);


      /**
       * Decrements each element of *this by a constant.
       *
       * @param arg Value by which array elements will be decremented.
       *
       * @return Reference to *this.
       */
      Array2D<Type>&
      operator-=(Type arg);


      /**
       * Multiplies each element of *this by a constant.
       *
       * @param arg Value by which array elements will be multiplied.
       *
       * @return Reference to *this.
       */
      Array2D<Type>&
      operator*=(Type arg);


      /**
       * Divides each element of *this by a constant.
       *
       * @param arg Value by which array elements will be divided.
       *
       * @return Reference to *this.
       */
      Array2D<Type>&
      operator/=(Type arg);


      /**
       * For integral Types, left-shifts each element of *this by the
       * specified number of bits.  This is equivalent to multiplying
       * by 2**numberOfBits.
       *
       * @param numberOfBits This argument specifies how many bits to
       * shift.
       *
       * @return Reference to *this, after the operation has been
       * completed..
       */
      Array2D<Type>&
      operator<<=(int numberOfBits);


      /**
       * For integral Types, right-shifts each element of *this by the
       * specified number of bits.  This is equivalent to dividing by
       * 2**numberOfBits.
       *
       * @param numberOfBits This argument specifies how many bits to
       * shift.
       *
       * @return Reference to *this, after the operation has been
       * completed..
       */
      Array2D<Type>&
      operator>>=(int numberOfBits);

    private:

      typedef brick::numeric::Index2D Key;
      struct Index2DComparator
        : public std::binary_function<Key, Key, bool>
      {
        bool operator()(Key const& key0, Key const& key1) {
          if(key0.getRow() > key1.getRow()) {return true;}
          if(key0.getColumn() > key1.getColumn()) {return true;}
          return false;
        }
      };


      /* ******** Private member functions ******** */
      /**
       * Allocate memory for array data and initialize reference count.
       */
      void
      allocate(size_t arrayRows, size_t arrayColumns);


      /**
       * Optionally (if and only if BRICK_NUMERIC_CHECKBOUNDS is
       * defined) throw an exception if either row or column is out of
       * range for this array.
       *
       * @param row This argument is the row index to check.
       *
       * @param column This argument  is the column index to check.
       */inline void
      checkBounds(size_t row, size_t column) const;


       /**
       * Release the memory for array data by decrementing the reference
       * count, or by returning the memory to the heap if it is not
       * referenced by any remaining Array instances.
       */
      void
      deAllocate();


      // Constants to help with formatting.  We use the initialization
      // on first use paradigm for the string constants to avoid
      // headaches.

      /**
       * Static constant describing how the string representation of an
       * Array2D should start.
       */
      inline static const std::string& ioIntro();

      /**
       * Static constant describing how the string representation of an
       * Array2D should end.
       */
      inline static const char& ioOutro();

      /**
       * Static constant describing how the the data portion of the
       * string representation of an Array2D should start.
       */
      inline static const char& ioOpening();

      /**
       * Static constant describing how the the data portion of the
       * string representation of an Array2D should end.
       */
      inline static const char& ioClosing();

      /**
       * Static constant describing how individual elements should be
       * separated in the string representation of Array2D.
       */
      inline static const char& ioSeparator();


      /* ********Private data members******** */
      size_t m_rows;
      size_t m_columns;
      size_t m_size;
      std::map<Key, Type>* m_dictionaryPtr;
      common::ReferenceCount m_referenceCount;
    };


    /* Non-member functions that will ultimately wind up in a different file */

    /**
     * This function returns an array that is the same size as its
     * argument, and in which the value of each element is the square
     * root of the corresponding element of the argument.
     *
     * @param array0 The square root calculation will be performed for
     * each element of this array.
     *
     * @return An array of square root values.
     */
    template <class Type>
    Array2D<Type>
    squareRoot(const Array2D<Type>& array0);


    /**
     * This function returns an array that is the same size as its
     * argument, and in which the value of each element is the square
     * root of the corresponding element of the argument.  It is really
     * just an alias for the function squareRoot(), below.
     *
     * @param array0 The square root calculation will be performed for
     * each element of this array.
     *
     * @return An array of square root values.
     */
    template <class Type>
    inline Array2D<Type>
    sqrt(const Array2D<Type>& array0);


    /**
     * Elementwise addition of Array2D instances.
     *
     * @param array0 First argument for addition.
     *
     * @param array1 Second argument for addition.
     *
     * @exception ValueException thrown when array sizes differ.
     *
     * @return Array2D instance in which the value of each element is
     * the sum of the values of the corresponding elements of the two
     * Array2D arguments.
     */
    template <class Type>
    Array2D<Type>
    operator+(const Array2D<Type>& array0,
              const Array2D<Type>& array1);


    /**
     * Elementwise subtraction of Array2D instances.
     *
     * @param array0 First argument for subtraction.
     *
     * @param array1 Second argument for subtraction.
     *
     * @exception ValueException thrown when array sizes differ.
     *
     * @return Array2D instance in which the value of each element is
     * the difference of the values of the corresponding elements of the two
     * Array2D arguments.
     */
    template <class Type>
    Array2D<Type>
    operator-(const Array2D<Type>& array0,
              const Array2D<Type>& array1);


    /**
     * Elementwise multiplication of Array2D instances.
     *
     * @param array0 First argument for multiplication.
     *
     * @param array1 Second argument for multiplication.
     *
     * @exception ValueException thrown when array sizes differ.
     *
     * @return Array2D instance in which the value of each element is
     * the product of the values of the corresponding elements of the two
     * Array2D arguments.
     */
    template <class Type>
    Array2D<Type>
    operator*(const Array2D<Type>& array0,
              const Array2D<Type>& array1);


    /**
     * Elementwise division of Array2D instances.
     *
     * @param array0 First argument for division.
     *
     * @param array1 Second argument for division.
     *
     * @exception ValueException thrown when array sizes differ.
     *
     * @return Array2D instance in which the value of each element is
     * the dividend of the values of the corresponding elements of the two
     * Array2D arguments.
     */
    template <class Type>
    Array2D<Type>
    operator/(const Array2D<Type>& array0,
              const Array2D<Type>& array1);


    /**
     * Addition of Array2D and scalar.
     *
     * @param array0 Array2D argument of the addition.
     *
     * @param scalar Scalar argument of the addition.
     *
     * @return Array2D instance in which the value of each element is
     * the sum of the corresponding element of the Array2D argument and
     * the scalar argument.
     */
    template <class Type>
    Array2D<Type>
    operator+(const Array2D<Type>& array0, Type scalar);


    /**
     * Subtraction of Array2D and scalar.
     *
     * @param array0 Array2D argument of the subtraction.
     *
     * @param scalar Scalar argument of the subtraction.
     *
     * @return Array2D instance in which the value of each element is
     * the difference of the corresponding element of the Array2D
     * argument and the scalar argument.
     */
    template <class Type>
    Array2D<Type>
    operator-(const Array2D<Type>& array0, Type scalar);


    /**
     * Multiplication of Array2D and scalar.
     *
     * @param array0 Array2D argument of the multiplication.
     *
     * @param scalar Scalar argument of the multiplication.
     *
     * @return Array2D instance in which the value of each element is
     * the product of the corresponding element of the Array2D argument
     * and the scalar argument.
     */
    template <class Type>
    Array2D<Type>
    operator*(const Array2D<Type>& array0, Type scalar);


    /**
     * Division of Array2D and scalar.
     *
     * @param array0 Array2D argument of the division.
     *
     * @param scalar Scalar argument of the division.
     *
     * @return Array2D instance in which the value of each element is
     * the difference of the corresponding element of the Array2D
     * argument and the scalar argument.
     */
    template <class Type>
    Array2D<Type>
    operator/(const Array2D<Type>& array0, Type scalar);


    /**
     * Addition of scalar and Array2D.
     *
     * @param scalar Scalar argument of the addition.
     *
     * @param array0 Array2D argument of the addition.
     *
     * @return Array2D instance in which the value of each element is
     * the sum of the scalar argument and the corresponding element of
     * the Array2D argument.
     */
    template <class Type>
    inline Array2D<Type>
    operator+(Type scalar, const Array2D<Type>& array0);


    /**
     * Multiplication of scalar and Array2D.
     *
     * @param scalar Scalar argument of the multiplication.
     *
     * @param array0 Array2D argument of the multiplication.
     *
     * @return Array2D instance in which the value of each element is
     * the product of the scalar argument and the corresponding element
     * of the Array2D argument.
     */
    template <class Type>
    inline Array2D<Type>
    operator*(Type scalar, const Array2D<Type>& array0);


    /**
     * Elementwise comparison of an Array2D with a constant.
     *
     * @param array0 An Array2D instance.
     *
     * @param arg Value to which the elements of array0 should be
     * compared.
     *
     * @return An Array2D<bool> in which each element has value "true"
     * if the corresponding element of array0 is equal to arg.
     */
    template <class Type>
    Array2D<bool>
    operator==(const Array2D<Type>& array0, const Type arg);


    /**
     * Elementwise comparison of an Array2D with another array.
     *
     * @param array0 An Array2D instance.
     *
     * @param array1 A second Array2D instance with the same size as
     * array0.
     *
     * @return An Array2D<bool> in which each element has value "true"
     * if the corresponding element of array0 is equal to the
     * corresponding element of array1.
     */
    template <class Type>
    Array2D<bool>
    operator==(const Array2D<Type>& array0, const Array2D<Type>& array1);


    /**
     * Elementwise comparison of Array2D with a constant.
     *
     * @param array0 An Array2D instance.
     *
     * @param arg Value to which the elements of array0 should be
     * compared.
     *
     * @return An Array2D<bool> in which each element has value "true"
     * if the corresponding element of array0 is greater than arg.
     */
    template <class Type>
    Array2D<bool>
    operator>(const Array2D<Type>& array0, Type arg);


    /**
     * Elementwise comparison of Array2D with a constant.
     *
     * @param array0 An Array2D instance.
     *
     * @param arg Value to which the elements of array0 should be
     * compared.
     *
     * @return An Array2D<bool> in which each element has value "true"
     * if the corresponding element of array0 is less than arg.
     */
    template <class Type>
    Array2D<bool>
    operator<(const Array2D<Type>& array0, Type arg);


    /**
     * Elementwise comparison of Array2D with a constant.
     *
     * @param array0 An Array2D instance.
     *
     * @param arg Value to which the elements of array0 should be
     * compared.
     *
     * @return An Array2D<bool> in which each element has value "true"
     * if the corresponding element of array0 is greater than or equal
     * to arg.
     */
    template <class Type>
    Array2D<bool>
    operator>=(const Array2D<Type>& array0, Type arg);


    /**
     * Elementwise comparison of Array2D with a constant.
     *
     * @param array0 An Array2D instance.
     *
     * @param arg Value to which the elements of array0 should be
     * compared.
     *
     * @return An Array2D<bool> in which each element has value "true"
     * if the corresponding element of array0 is less than or equal to
     * arg.
     */
    template <class Type>
    Array2D<bool>
    operator<=(const Array2D<Type>& array0, Type arg);


    /**
     * Outputs a text representation of an Array2D instance to a
     * std::ostream.  The output format looks like this:
     *
     * Array2D([[1, 2, 4, 8, 16],
     *          [5, 1, 6, 7, 2]])
     *
     * Where the array elements are output using
     * operator<<(std::ostream&, const Type&)
     * and each element is separated from its neighbors by a comma and
     * whitespace.
     *
     * @param stream Reference to the the output stream.
     *
     * @param array0 Const reference to the Array2D to be output.
     *
     * @return Reference to output stream.
     */
    template <class Type>
    std::ostream&
    operator<<(std::ostream& stream, const Array2D<Type>& array0);


    /**
     * Sets the value of an Array2D instance from a std::istream.
     * The input format is as described for
     * operator<<(std::ostream&, const Array2D<Type>&) above.
     *
     * @param stream Reference to the the input stream.
     *
     * @param array0 Reference to the Array2D that will take the input.
     *
     * @return Reference to input stream.
     */
    template <class Type>
    std::istream&
    operator>>(std::istream& stream, Array2D<Type>& array0);

  } // namespace sparse

} // namespace brick

#include <brick/sparse/array2D_impl.hh>

#endif /* #ifdef BRICK_SPARSE_ARRAY2D_IMPL_HH */
