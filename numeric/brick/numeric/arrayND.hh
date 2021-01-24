/**
***************************************************************************
* @file brick/numeric/arrayND.hh
*
* Header file declaring ArrayND class template.
*
* Copyright (C) 2008, 2011 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
***************************************************************************
**/

#ifndef BRICK_NUMERIC_ARRAYND_HH
#define BRICK_NUMERIC_ARRAYND_HH

#include <iostream>
#include <string>
#include <brick/common/exception.hh>
#include <brick/numeric/array1D.hh>

namespace brick {


  namespace numeric {

    /**
     ** @warning WARNING: this class is very immature, is not well
     ** tested, and has an interface that is subject to change.
     **
     ** The ArrayND class template represents a multi-dimensional
     ** array of arbitrary type.  This class has internal reference
     ** counting.
     **
     ** Element storage in the ArrayND class is contiguous along the
     ** last axis, then the next-to-last axis, etc.  To make this
     ** clearer, you can think of the 2-dimensional array in the
     ** example below as having 3 rows and 4 columns, with the
     ** elements stored in raster order, starting from the top-left
     ** corner of the array.
     **
     ** @code
     **   size_t const shape[] = {3, 4};
     **   ArrayND<2, int> exampleArray3x4(shape);
     ** @endcode
     **
     ** In the current implementation, there is no provision for
     ** padding between rows of the array (for example, to align on
     ** 16-byte boundaries).  We may add support for this later,
     ** however the default behavior will always be to have no padding
     ** unless it is explicitly enabled.  If your code counts on
     ** contiguous storage in an ArrayND instance passed in from code
     ** you do not control, you can make sure there's no padding by
     ** confirming that the isContiguous() method returns true.
     **
     ** IMPORTANT: This class does _shallow_ copies by default.  If you
     ** type:
     **
     **       array1 = array2;
     **
     ** then array1 and array2 point to the same data.  To do a deep copy,
     ** type
     **
     **       array1.copy(array2);
     **
     ** or
     **
     **       array1 = array2.copy();
     **
     ** The advantage of the first form is that it doesn't involve allocating
     ** memory.  The advantage of the second form is that there's no error if
     ** array1 and array2 have different shapes.
     **/
    template <size_t Dimension, class Type>
    class ArrayND {
    public:
      /* ======== Public typedefs ======== */

      /**
       ** Typedef for value_type describes the contents of the array.
       **/
      typedef Type value_type;

      /**
       ** Typedef for iterator type helps with standard library interface.
       **/
      typedef Type* iterator;

      /**
       ** Typedef for const_iterator type helps with standard library
       ** interface.
       **/
      typedef const Type* const_iterator;


      /* ======== Public member functions ======== */

      /**
       * Default constructor initializes to zero size.
       */
      ArrayND();


      /**
       * Constructs an array with a user-specified shape.  For
       * example, the following code constructs a 2-dimensional
       * ArrayND instance with 5 rows and 6 columns:
       *
       * @code
       *   size_t shape[] = {5, 6};
       *   ArrayND<2, int> myArray(2, shape);
       * @endcode
       *
       * Here's a second example, which creates a 3-dimensional
       * ArrayND instance with 2 levels of 5 rows and 6 columns:
       *
       * @code
       *   size_t const shape[] = {2, 5, 6};
       *   ArrayND<3, int> myArray(3, shape);
       * @endcode
       *
       * @param dimension This argument specifies the dimension
       * (number of axes) in the array.  Set this to 1 if you want a
       * vector, 2 if you want a matrix, 3 if you want a 3D array,
       * etc.
       *
       * @param shape This array specifies the size (number of
       * elements along each axis) of the new ArrayND instance.
       */
      explicit
      ArrayND(size_t dimension, size_t const shape[]);


      /**
       * Constructs an array with a user-specified shape.  For
       * example, the following code constructs a 2-dimensional
       * ArrayND instance with 5 rows and 6 columns:
       *
       * @code
       *   ArrayND<2, int> myArray(Array1D<size_t>("[5, 6]"));
       * @endcode
       *
       * Here's a second example, which creates a 3-dimensional
       * ArrayND instance with 2 levels of 5 rows and 6 columns:
       *
       * @code
       *   ArrayND<3, int> myArray(Array1D<size_t>("[2, 5, 6]"));
       * @endcode
       *
       * @param shape This array specifies the size (number of
       * elements along each axis) of the new ArrayND instance.
       */
      explicit
      ArrayND(const Array1D<size_t>& shape);


      /**
       * The copy constructor does a shallow copy.  The newly created
       * array points to the same data as copied array.
       *
       * @param source The ArrayND<> instance to be copied.
       */
      ArrayND(const ArrayND<Dimension, Type> &source);


      /**
       * Construct from an initialization string.  Ultimately, this
       * constructor will be very flexible about interpreting string
       * syntax, but right now it only accepts a secret format known
       * only to the developers of this library.
       *
       * @param inputString The argument specifies the string from which
       * the array will be constructed.
       */
      explicit
      ArrayND(const std::string& inputString);


      /**
       * Destroys the ArrayND instance and deletes the internal data
       * store if no remaining arrays point to it.
       */
      ~ArrayND();


      /**
       * Return begin() iterator for Standard Library algorithms.
       *
       * @return Iterator pointing to the first element of the array.
       */
      iterator
      begin() {return m_storage.begin();}


      /**
       * Return begin() const_iterator for Standard Library algorithms.
       *
       * @return Const iterator pointing to the first element of the
       * array.
       */
      const_iterator
      begin() const {return m_storage.begin();}


      /**
       * Optionally throw an exception if the shape of *this is
       * different than specified.
       *
       * @param size The required array size.
       */
      inline void
      checkDimension(Array1D<size_t> const& shape) const;


      /**
       * Reset the array to zero size, abandoning all contents.  This is
       * equivalent to this->reinit(0);
       */
      void
      clear() {m_shape.reinit(0); m_storage.reinit(0);}


      /**
       * Allocate a new array and deep copy the contents of *this.
       *
       * @return A new array which is a (deep) copy of *this.
       */
      ArrayND<Dimension, Type>
      copy() const;


      /**
       * Deep copies the contents of source.  It is an error if source
       * does not have the same shape as *this.
       *
       * @param source The array to be copied.
       *
       * @exception ValueException on incompatible array sizes
       *
       * @see void ArrayND<size_t, Type>::copy(const Type2* dataPtr)
       */
      template <class Type2> void
      copy(const ArrayND<Dimension, Type2>& source);


      /**
       * Copies elements from dataPtr.  There must be valid data at all
       * addresses from dataPtr to (dataPtr + this->size());
       *
       * @param dataPtr Pointer to the data to be copied.
       */
      template <class Type2> void
      copy(const Type2* dataPtr);


      /**
       * Returns a pointer to the internal data store.  This is ugly
       * but often necessary for interfacing with external libraries.
       * In the current implementation, data is stored contiguously.
       * At some point, we may add support for padding to byte-align
       * subsequent rows.  In this event, any ArrayND instance for
       * which the returned data is not contiguous will also return
       * false from member function isContiguous().
       *
       * @return Pointer to the internal data store.
       */
      Type*
      data() {return m_storage.data();}


      /**
       * This version of data(void) is appropriate for const ArrayND,
       * and returns a pointer-to-const.
       *
       * @return Const pointer to the internal data store.
       *
       * @see Type* ArrayND<size_t, Type>::data()
       */
      const Type*
      data() const {return m_storage.data();}


      /**
       * This member function returns true if the array instance
       * contains no elements.  It has complexity O(1).
       *
       * @return The return value indicates whether or not the array is
       * empty.
       */
      bool
      empty() const {return m_storage.empty();}


      /**
       * Return end() iterator for Standard Library algorithms.
       *
       * @return Iterator pointing just past the last element of
       * the array.
       */
      iterator
      end() {return m_storage.end();}


      /**
       * Return end() const_iterator for Standard Library algorithms.
       *
       * @return Const iterator pointing just past the last element of
       * the array.
       */
      const_iterator
      end() const {return m_storage.end();}


      /**
       * This member function converts a multi-dimensional index (such
       * as would be passed to getElement(Array1D<size_t> const&))
       * into a scalar index (such as would be passed to
       * getElement(size_t).
       *
       * @param indexArray This argument is the array to be converted.
       *
       * @return The return value is the equivalent 1D index.
       */
      size_t
      flattenIndex(const Array1D<size_t>& indexArray) const;


      /**
       * This member function returns a specific element of the array
       * by value.
       *
       * @param index0 This argument specifies which element value
       * should be returned.
       *
       * @return The return value is a copy of the requested element.
       */
      inline Type
      getElement(size_t index0) const {return this->operator()(index0);}


      /**
       * This member function returns a specific element of the array
       * by value.
       *
       * @param index0 This argument specifies which element value
       * should be returned.
       *
       * @return The return value is a copy of the requested element.
       */
      inline Type
      getElement(Array1D<size_t> const& index0) const {
        return this->operator()(index0);
      }


      /**
       * This member function returns a (reference to an) array
       * describing the shape of *this.  Note that Array1D copy
       * semantics are _shallow_, so you'll actually get an array that
       * references memory inside *this.  If you copy (shallow copy)
       * this array, and then modify its contents, you will break the
       * internal state of *this.
       *
       * @return The return value is an array describing the number of
       * elements along each axis of *this.
       */
      inline const Array1D<size_t>&
      getShape() const {return m_shape;}


      /**
       * This member function establishes the relationship between
       * single-indexing (for example using getElement(size_t)) and
       * array indexing (for example using getElement(Array1D<size_t>
       * const&)) by returning the single-index offset between
       * adjacent elements along the specified axis.
       *
       * @param axis This argument specifies the axis along which we
       * wish to move.
       *
       * @return The return value is the number by which the argument
       * to getElement(size_t) must be incremented in order to move
       * one element along the specified axis.
       */
      size_t
      getStride(size_t axis) const {return m_strideArray[axis];}


      /**
       * This member function returns an array in which each element
       * contains the distance between adjacent elements along the
       * corresponding axis.  See member function getStride() for more
       * information.  Note that Array1D copy semantics are _shallow_,
       * so you'll actually get an array that references memory inside
       * *this.  If you copy (shallow copy) this array, and then
       * modify its contents, you will break the internal state of
       * *this.
       *
       * @return The return value is an array in which the ii-th
       * element is equal to the result of this->getStride(ii).
       */
      Array1D<size_t> const&
      getStrideArray() const {return m_strideArray;}


      /**
       * Indicates whether the internal data array is being managed (and
       * reference counted) by *this.  This member function is only
       * needed in very unusual circumstances.
       *
       * @return The return value is a bool indicating whether the internal
       * data is being managed by *this.
       */
      bool
      isAllocated() const {return true;}


      /**
       * Returns an Array1D, with size equal to this->size(), which
       * references the same data as *this.  In other words, ravel()
       * returns a flattened version of *this.
       *
       * @return Array1D referencing the same data as *this.
       */
      Array1D<Type>
      ravel() {return m_storage;}


      /**
       * Returns an Array1D, with size equal to this->size(), which
       * references the same data as *this.  In other words, ravel()
       * returns a flattened version of *this.
       *
       * @return Array1D referencing the same data as *this.
       */
      Array1D<Type> const
      ravel() const {return m_storage;}


      /**
       * Changes the shape of the array and reallocates storage.
       * The current array contents are lost.
       *
       * @param size Number of elements in the array after reallocation.
       */
      void
      reinit(Array1D<size_t> const& shape);


      /**
       * This member function sets the value of a specific element of
       * the array.
       *
       * @param index0 This argument specifies which element value
       * should be set.
       *
       * @param value This argument will be copied into the selected
       * array element.
       *
       * @return The return value is a reference to the newly set
       * array element.
       */
      Type&
      setElement(size_t index0, const Type& value) {
        return this->operator()(index0) = value;
      }


      /**
       * This member function sets the value of a specific element of
       * the array.
       *
       * @param index0 This argument specifies which element value
       * should be set.
       *
       * @param value This argument will be copied into the selected
       * array element.
       *
       * @return The return value is a reference to the newly set
       * array element.
       */
      Type&
      setElement(Array1D<size_t> const& index0, const Type& value) {
        return this->operator()(this->flattenIndex(index0)) = value;
      }


      /**
       * Returns the number of elements in the array.
       *
       * @return The number of elements in the array.
       */
      size_t
      size() const {return m_storage.size();}


      /**
       * Assignment operator shallow copies the contents of source.
       * After the copy, both arrays reference the same data.
       *
       * @param source The ArrayND instance to be copied.
       *
       * @return Reference to *this.
       */
      ArrayND<Dimension, Type>&
      operator=(const ArrayND<Dimension, Type>& source);


      /**
       * Assign value to every element in the array.
       *
       * @param value The value to be copied.
       *
       * @return Reference to *this.
       */
      inline ArrayND<Dimension, Type>&
      operator=(Type value);


      /**
       * Returns the (index)th element of the array by reference.
       *
       * @param index Indicates the selected element.
       *
       * @return Reference to the (index)th element of the array.
       */
      Type&
      operator()(size_t index0) {return m_storage(index0);}


      /**
       * Returns the (index)th element of the array by value.
       *
       * @param index Indicates the selected element.
       *
       * @return Value of the (index)th element of the array.
       */
      Type const& operator()(size_t index0) const {return m_storage(index0);}


      /**
       * Returns the (index)th element of the array by reference.
       *
       * @param index Indicates the selected element.
       *
       * @return Reference to the (index)th element of the array.
       */
      Type&
      operator()(Array1D<size_t> const& index0) {
        return m_storage[this->flattenIndex(index0)];
      }


      /**
       * Returns the (index)th element of the array by value.
       *
       * @param index Indicates the selected element.
       *
       * @return Value of the (index)th element of the array.
       */
      Type const& operator()(Array1D<size_t> const& index0) const {
        return m_storage[this->flattenIndex(index0)];
      }


      /**
       * Returns the (index)th element of the array by reference.
       * Synonymous with operator()().
       *
       * @param index Indicates the selected element.
       *
       * @return Reference to the (index)th element of the array.
       */
      Type& operator[](size_t index) {return this->operator()(index);}


      /**
       * Returns the (index)th element of the array by value.
       * Synonymous with operator()() const.
       *
       * @param index Indicates the selected element.
       *
       * @return Value of the (index)th element of the array.
       */
      Type const&
      operator[](size_t index) const {return this->operator()(index);}


      /**
       * Increments each element of *this by the value of the
       * corresponding element of arg.
       *
       * @param arg ArrayND of values to be added to the elements of
       * *this.
       *
       * @exception ValueException on incompatible array sizes
       *
       * @return Reference to *this.
       */
      template <class Type2>
      ArrayND<Dimension, Type>&
      operator+=(const ArrayND<Dimension, Type2>& arg);


      /**
       * Increments each element of *this by a constant.
       *
       * @param arg Value by which array elements will be incremented.
       *
       * @return Reference to *this.
       */
      ArrayND<Dimension, Type>&
      operator+=(const Type arg);


      /**
       * Decrements each element of *this by the value of the
       * corresponding element of arg.
       *
       * @param arg ArrayND of values to be subtracted from the elements
       * of *this.
       *
       * @exception ValueException on incompatible array sizes
       *
       * @return Reference to *this.
       */
      template <class Type2>
      ArrayND<Dimension, Type>&
      operator-=(const ArrayND<Dimension, Type2>& arg);


      /**
       * Decrements each element of *this by a constant.
       *
       * @param arg Value by which array elements will be decremented.
       *
       * @return Reference to *this.
       */
      ArrayND<Dimension, Type>&
      operator-=(const Type arg);


      /**
       * Multiplies each element of *this by the value of the
       * corresponding element of arg.
       *
       * @param arg ArrayND of values by which the elements of *this
       * should be multiplied.
       *
       * @exception ValueException on incompatible array sizes
       *
       * @return Reference to *this.
       */
      template <class Type2>
      ArrayND<Dimension, Type>&
      operator*=(const ArrayND<Dimension, Type2>& arg);


      /**
       * Multiplies each element of *this by a constant.
       *
       * @param arg Value by which array elements will be multiplied.
       *
       * @return Reference to *this.
       */
      ArrayND<Dimension, Type>&
      operator*=(const Type arg);


      /**
       * Divides each element of *this by the value of the corresponding
       * element of arg.
       *
       * @param arg ArrayND of values by which the elements of *this
       * should be divided.
       *
       * @exception ValueException on incompatible array sizes
       *
       * @return Reference to *this.
       */
      template <class Type2>
      ArrayND<Dimension, Type>&
      operator/=(const ArrayND<Dimension, Type2>& arg);


      /**
       * Divides each element of *this by a constant.
       *
       * @param arg Value by which array elements will be divided.
       *
       * @return Reference to *this.
       */
      ArrayND<Dimension, Type>&
      operator/=(const Type arg);

    private:
      /* ======== Private member functions ======== */

      size_t computeSize(Array1D<size_t> const& shape) const;
      Array1D<size_t> computeStride(Array1D<size_t> const& shape) const;


      /* ======== Private data members ======== */
      Array1D<size_t> m_shape;
      Array1D<Type> m_storage;
      Array1D<size_t> m_strideArray;
    };


    /* Non-member functions which should maybe wind up in a different file */

    /**
     * Elementwise addition of ArrayND instances.
     *
     * @param array0 First argument for addition.
     *
     * @param arrayN Second argument for addition.
     *
     * @exception ValueException on incompatible array sizes
     *
     * @return ArrayND instance in which the value of each element is
     * the sum of the values of the corresponding elements of the two
     * ArrayND arguments.
     */
    template <size_t Dimension, class Type>
    ArrayND<Dimension, Type>
    operator+(const ArrayND<Dimension, Type>& array0,
              const ArrayND<Dimension, Type>& arrayN);


    /**
     * Elementwise subtraction of ArrayND instances.
     *
     * @param array0 First argument for subtraction.
     *
     * @param arrayN Second argument for subtraction.
     *
     * @exception ValueException on incompatible array sizes
     *
     * @return ArrayND instance in which the value of each element is
     * the difference of the values of the corresponding elements of the two
     * ArrayND arguments.
     */
    template <size_t Dimension, class Type>
    ArrayND<Dimension, Type>
    operator-(const ArrayND<Dimension, Type>& array0,
              const ArrayND<Dimension, Type>& arrayN);


    /**
     * Elementwise multiplication of ArrayND instances.
     *
     * @param array0 First argument for multiplication.
     *
     * @param arrayN Second argument for multiplication.
     *
     * @exception ValueException on incompatible array sizes
     *
     * @return ArrayND instance in which the value of each element is
     * the product of the values of the corresponding elements of the two
     * ArrayND arguments.
     */
    template <size_t Dimension, class Type>
    ArrayND<Dimension, Type>
    operator*(const ArrayND<Dimension, Type>& array0,
              const ArrayND<Dimension, Type>& arrayN);


    /**
     * Elementwise division of ArrayND instances.
     *
     * @param array0 First argument for division.
     *
     * @param arrayN Second argument for division.
     *
     * @exception ValueException on incompatible array sizes
     *
     * @return ArrayND instance in which the value of each element is
     * the dividend of the values of the corresponding elements of the two
     * ArrayND arguments.
     */
    template <size_t Dimension, class Type>
    ArrayND<Dimension, Type>
    operator/(const ArrayND<Dimension, Type>& array0,
              const ArrayND<Dimension, Type>& arrayN);


    /**
     * Addition of ArrayND and scalar.
     *
     * @param array ArrayND argument of the addition.
     *
     * @param scalar Scalar argument of the addition.
     *
     * @return ArrayND instance in which the value of each element is
     * the sum of the corresponding element of the ArrayND argument and
     * the scalar argument.
     */
    template <size_t Dimension, class Type>
    ArrayND<Dimension, Type>
    operator+(const ArrayND<Dimension, Type>& array, Type scalar);


    /**
     * Subtraction of ArrayND and scalar.
     *
     * @param array0 ArrayND argument of the subtraction.
     *
     * @param scalar Scalar argument of the subtraction.
     *
     * @return ArrayND instance in which the value of each element is
     * the difference of the corresponding element of the ArrayND
     * argument and the scalar argument.
     */
    template <size_t Dimension, class Type>
    ArrayND<Dimension, Type>
    operator-(const ArrayND<Dimension, Type>& array0, Type scalar);


    /**
     * Multiplication of ArrayND and scalar.
     *
     * @param array0 ArrayND argument of the multiplication.
     *
     * @param scalar Scalar argument of the multiplication.
     *
     * @return ArrayND instance in which the value of each element is
     * the product of the corresponding element of the ArrayND argument
     * and the scalar argument.
     */
    template <size_t Dimension, class Type>
    ArrayND<Dimension, Type>
    operator*(const ArrayND<Dimension, Type>& array0, Type scalar);


    /**
     * Division of ArrayND and scalar.
     *
     * @param array0 ArrayND argument of the division.
     *
     * @param scalar Scalar argument of the division.
     *
     * @return ArrayND instance in which the value of each element is
     * the difference of the corresponding element of the ArrayND
     * argument and the scalar argument.
     */
    template <size_t Dimension, class Type>
    ArrayND<Dimension, Type>
    operator/(const ArrayND<Dimension, Type>& array0, Type scalar);


    /**
     * Addition of scalar and ArrayND.
     *
     * @param scalar Scalar argument of the addition.
     *
     * @param array0 ArrayND argument of the addition.
     *
     * @return ArrayND instance in which the value of each element is
     * the sum of the scalar argument and the corresponding element of
     * the ArrayND argument.
     */
    template <size_t Dimension, class Type>
    inline ArrayND<Dimension, Type>
    operator+(Type scalar, const ArrayND<Dimension, Type>& array0);


    /**
     * Subtraction of scalar and ArrayND.
     *
     * @param scalar Scalar argument of the subtraction.
     *
     * @param array0 ArrayND argument of the subtraction.
     *
     * @return ArrayND instance in which the value of each element is
     * the difference of the scalar argument and the corresponding
     * element of the ArrayND argument.
     */
    template <size_t Dimension, class Type>
    inline ArrayND<Dimension, Type>
    operator-(Type scalar, const ArrayND<Dimension, Type>& array0);


    /**
     * Multiplication of scalar and ArrayND.
     *
     * @param scalar Scalar argument of the multiplication.
     *
     * @param array0 ArrayND argument of the multiplication.
     *
     * @return ArrayND instance in which the value of each element is
     * the product of the scalar argument and the corresponding element
     * of the ArrayND argument.
     */
    template <size_t Dimension, class Type>
    inline ArrayND<Dimension, Type>
    operator*(Type scalar, const ArrayND<Dimension, Type>& array0);


    /**
     * Division of scalar and ArrayND.
     *
     * @param scalar Scalar argument of the division.
     *
     * @param array0 ArrayND argument of the division.
     *
     * @return ArrayND instance in which the value of each element is
     * the dividend of the scalar argument and the corresponding element
     * of the ArrayND argument.
     */
    template <size_t Dimension, class Type>
    inline ArrayND<Dimension, Type>
    operator/(Type scalar, const ArrayND<Dimension, Type>& array0);


    /**
     * Elementwise comparison of an ArrayND with a constant.
     *
     * @param array0 An ArrayND instance.
     *
     * @param arg Value to which the elements of array0 should be
     * compared.
     *
     * @return An ArrayND<bool> in which each element has value "true"
     * if the corresponding element of array0 is equal to arg.
     */
    template <size_t Dimension, class Type>
    ArrayND<Dimension, bool>
    operator==(const ArrayND<Dimension, Type>& array0, const Type arg);


    /**
     * Elementwise comparison of an ArrayND with another array.
     *
     * @param array0 An ArrayND instance.
     *
     * @param arrayN A second ArrayND instance with the same size as
     * array0.
     *
     * @return An ArrayND<bool> in which each element has value "true"
     * if the corresponding element of array0 is equal to the
     * corresponding element of arrayN.
     */
    template <size_t Dimension, class Type>
    ArrayND<Dimension, bool>
    operator==(const ArrayND<Dimension, Type>& array0,
               const ArrayND<Dimension, Type>& arrayN);


    /**
     * Elementwise comparison of an ArrayND with a constant.
     *
     * @param array0 An ArrayND instance.
     *
     * @param arg Value to which array elements should be compared.
     *
     * @return An ArrayND<bool> in which each element has value "true"
     * if the corresponding element of array0 is greater than arg.
     */
    template <size_t Dimension, class Type>
    ArrayND<Dimension, bool>
    operator>(const ArrayND<Dimension, Type>& array0, const Type arg);


    /**
     * Elementwise comparison of an ArrayND with a constant.
     *
     * @param array0 An ArrayND instance.
     *
     * @param arg Value to which array elements should be compared.
     *
     * @return An ArrayND<bool> in which each element has value "true"
     * if the corresponding element of array0 is greater than or equal
     * to arg.
     */
    template <size_t Dimension, class Type>
    ArrayND<Dimension, bool>
    operator>=(const ArrayND<Dimension, Type>& array0, const Type arg);


    /**
     * Elementwise comparison of an ArrayND with a constant.
     *
     * @param array0 An ArrayND instance.
     *
     * @param arg Value to which array elements should be compared.
     *
     * @return An ArrayND<bool> in which each element has value "true"
     * if the corresponding element of array0 is less than arg.
     */
    template <size_t Dimension, class Type>
    ArrayND<Dimension, bool>
    operator<(const ArrayND<Dimension, Type>& array0, const Type arg);


    /**
     * Elementwise comparison of an ArrayND with a constant.
     *
     * @param array0 An ArrayND instance.
     *
     * @param arg Value to which array elements should be compared.
     *
     * @return An ArrayND<bool> in which each element has value "true"
     * if the corresponding element of array0 is less than or equal to
     * arg.
     */
    template <size_t Dimension, class Type>
    ArrayND<Dimension, bool>
    operator<=(const ArrayND<Dimension, Type>& array0, const Type arg);


    template <size_t Dimension, class Type>
    std::ostream& operator<<(std::ostream& stream,
                             const ArrayND<Dimension, Type>& array0);


    template <size_t Dimension, class Type>
    std::istream& operator>>(std::istream& stream,
                             ArrayND<Dimension, Type> & array0);

  } // namespace numeric

} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/numeric/arrayND_impl.hh>

#endif /* #ifdef BRICK_NUMERIC_ARRAYND_HH */
