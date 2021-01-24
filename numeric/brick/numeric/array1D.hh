/**
***************************************************************************
* @file brick/numeric/array1D.hh
*
* Header file declaring Array1D class.
*
* Copyright (C) 2001-2011 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_ARRAY1D_HH
#define BRICK_NUMERIC_ARRAY1D_HH

#include <iostream>
#include <string>
#include <brick/common/exception.hh>
#include <brick/common/referenceCount.hh>

namespace brick {


  /**
   ** This namespace contains code for 1D, 2D, and 3D arrays,
   ** matrices, coordinate transformations, rotation conversions, and
   ** much more.
   **/
  namespace numeric {

    /**
     ** The Array1D class template represents a 1D array of arbitrary type.
     ** This class has internal reference counting.
     **
     ** IMPORTANT: This class does _shallow_ copies by default.  If you
     ** type:
     **
     ** @code
     **       array1 = array2;
     ** @endcode
     **
     ** then array1 and array2 point to the same data.  To do a deep copy,
     ** type
     **
     ** @code
     **       array1.copy(array2);
     ** @endcode
     **
     ** or
     **
     ** @code
     **       array1 = array2.copy();
     ** @endcode
     **
     ** The advantage of the first form is that it doesn't involve allocating
     ** memory.  The advantage of the second form is that there's no error if
     ** array1 and array2 have different shapes/sizes.
     **/
    template <class Type>
    class Array1D {
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
      Array1D();


      /**
       * Constructs an "arraySize" element array.
       *
       * @param arraySize Number of elements in the array after successful
       * construction.
       */
      explicit
      Array1D(size_t arraySize);


      /**
       * Constructs an Array1D by parsing an input string.  For example,
       * the code line
       *
       *   Array1D<double> testArray("[1.0, 2.0, 3.0]");
       *
       * would construct a 3 element array.
       *
       * @param inputString A formatted string which could be reasonably
       * expected to parse into the desired array values.
       */
      explicit
      Array1D(const std::string& inputString);


      /**
       * The copy constructor does a shallow copy.  The newly created
       * array points to the same data as copied array.
       *
       * @param source The Array1D<> instance to be copied.
       */
      Array1D(const Array1D<Type> &source);


      /**
       * Construct an array around external data.  Arrays constructed
       * in this way will not implement reference counting, and will
       * not delete dataPtr when done.
       *
       * @param arraySize Number of elements in the array after construction.
       *
       * @param dataPtr A C-style array of Type into which the newly
       * constructed Array1D should index.
       */
      Array1D(size_t arraySize, Type* const dataPtr);


      /**
       * Construct an array around external data which was allocated by
       * an Array?D instance.  Arrays constructed in this way _do_
       * implement reference counting, and will delete dataPtr when
       * done.  This constructor is provided primarily so that higher
       * dimensionality array classes can return Array1D instances which
       * reference their data without being friend classes.  Caveat
       * emptor.
       *
       * @param arraySize Number of elements in the array after construction.
       *
       * @param dataPtr A C-style array of Type into which the newly
       * constructed Array1D should index.
       *
       * @param referenceCountPtr ReferenceCount instance indicating
       * the number of Array classes currently using dataPtr.
       */
      Array1D(size_t arraySize, Type* const dataPtr,
              common::ReferenceCount const& referenceCount);


      /**
       * Construct an array using an initializer list.  This lets users
       * build arrays using the following syntax:
       *
       * @code
       *    Array1D<double> {0.5, 2.0, 11.0, -0.1};
       * @endCode
       *
       * @param initializer This argument is generated automatically
       * by the compiler.
       */
      Array1D(std::initializer_list<Type> initializer);


      /**
       * Destroys the Array1D instance and deletes the internal data
       * store if no remaining arrays point to it.
       */
      ~Array1D();


      /**
       * Return begin() iterator for Standard Library algorithms.
       *
       * @return Iterator pointing to the first element of the array.
       */
      iterator
      begin() {return m_dataPtr;}


      /**
       * Return begin() const_iterator for Standard Library algorithms.
       *
       * @return Const iterator pointing to the first element of the
       * array.
       */
      const_iterator
      begin() const {return m_dataPtr;}


      /**
       * Optionally throw an exception if the shape of *this is
       * different than specified.
       *
       * @param arraySize The required array size.
       */
      inline void
      checkDimension(size_t arraySize) const;


      /**
       * Reset the array to zero size, abandoning all contents.  This is
       * equivalent to this->reinit(0);
       */
      void
      clear() {this->reinit(0);}


      /**
       * Allocate a new array and deep copy the contents of *this.
       *
       * @return A new array which is a (deep) copy of *this.
       */
      Array1D<Type>
      copy() const;


      /**
       * Deep copies the contents of source.  It is an error if source
       * does not have the same size as *this.
       *
       * @param source The array to be copied.
       *
       * @exception ValueException on incompatible array sizes
       *
       * @see void Array1D<Type>::copy(const Type2* dataPtr)
       */
      template <class Type2> void
      copy(const Array1D<Type2>& source);


      /**
       * Copies elements from dataPtr.  There must be valid data at all
       * addresses from dataPtr to (dataPtr + this->size());
       *
       * @param dataPtr Pointer to the data to be copied.
       */
      template <class Type2> void
      copy(const Type2* dataPtr);


      /**
       * This is an alias for member function getData().  It returns a
       * pointer to the internal data store.  This is ugly but often
       * necessary for interfacing with external libraries.  Data is
       * stored contiguously.
       *
       * @return Pointer to the internal data store.
       */
      Type*
      data() {return m_dataPtr;}


      /**
       * This is an alias for member function getData() const.  This
       * version of data(void) is appropriate for const Array1D, and
       * returns a pointer-to-const.
       *
       * @return Const pointer to the internal data store.
       */
      const Type*
      data() const {return m_dataPtr;}


      /**
       * This is an alias for member function getData(size_t).  Just
       * like data(void), which is documented above, but returns a
       * pointer to the (index)th element instead of the first
       * element.  Note that "data(index)" is synonymous with "data()
       * + index" .
       *
       * @param index Indicates to which element of the array the return
       * value should point.
       *
       * @return Pointer to the (index)th element of the array.
       */
      Type*
      data(size_t index) {
        this->checkBounds(index);
        return m_dataPtr + index;
      }


      /**
       * This is an alias for member function getData(size_t) const.
       * This version of data(size_t) is appropriate for const
       * Array1D, and returns a pointer-to-const.
       *
       * @param index Indicates to which element of the array the return
       * value should point.
       *
       * @return Const pointer to the (index)th element of the array.
       */
      const Type*
      data(size_t index) const {
        this->checkBounds(index);
        return m_dataPtr + index;
      }


      /**
       * This is an alias for member function isEmpty(), and returns
       * true if the array instance contains no elements.  It has
       * complexity O(1).
       *
       * @return The return value indicates whether or not the array is
       * empty.
       */
      bool
      empty() const {return this->size() == 0;}


      /**
       * Return end() iterator for Standard Library algorithms.
       *
       * @return Iterator pointing just past the last element of
       * the array.
       */
      iterator
      end() {return m_dataPtr + m_size;}


      /**
       * Return end() const_iterator for Standard Library algorithms.
       *
       * @return Const iterator pointing just past the last element of
       * the array.
       */
      const_iterator
      end() const {return m_dataPtr + m_size;}


      /**
       * Returns a pointer to the internal data store.  This is ugly
       * but often necessary for interfacing with external libraries.
       * Data is stored contiguously.
       *
       * @return Pointer to the internal data store.
       */
      Type*
      getData() {return m_dataPtr;}


      /**
       * This version of getData(void) is appropriate for const
       * Array1D, and returns a pointer-to-const.
       *
       * @return Const pointer to the internal data store.
       */
      const Type*
      getData() const {return m_dataPtr;}


      /**
       * Just like getData(void), which is documented above, but
       * returns a pointer to the (index)th element instead of the
       * first element.  Note that "getData(index)" is synonymous with
       * "getData() + index" .
       *
       * @param index Indicates to which element of the array the return
       * value should point.
       *
       * @return Pointer to the (index)th element of the array.
       */
      Type*
      getData(size_t index) {
        this->checkBounds(index);
        return m_dataPtr + index;
      }


      /**
       * This version of getData(size_t) is appropriate for const
       * Array1D, and returns a pointer-to-const.
       *
       * @param index Indicates to which element of the array the return
       * value should point.
       *
       * @return Const pointer to the (index)th element of the array.
       */
      const Type*
      getData(size_t index) const {
        this->checkBounds(index);
        return m_dataPtr + index;
      }


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
       * Returns the number of elements in the array.  This is a synonym
       * for getSize().
       *
       * @return The number of elements in the array.
       */
      size_t
      getLength() const {return this->getSize();}


      /**
       * Returns the number of elements in the array.
       *
       * @return The number of elements in the array.
       */
      size_t
      getSize() const {return m_size;}


      /**
       * Returns the the internal ReferenceCount instance by
       * reference.  This member function is included to support
       * certain tests, and should genrally not be used in client
       * code.
       *
       * @return A reference to the internal referenceCount instance.
       */
      common::ReferenceCount const&
      getReferenceCount() const {return m_referenceCount;}


      /**
       * This member function returns true if the array instance
       * contains no elements.  It has complexity O(1).
       *
       * @return The return value indicates whether or not the array is
       * empty.
       */
      bool
      isEmpty() const {return this->size() == 0;}


      /**
       * Indicates whether the internal data array is being managed (and
       * reference counted) by *this.  This member function is only
       * needed in very unusual circumstances.
       *
       * @return The return value is a bool indicating whether the internal
       * data is being managed by *this.
       */
      bool
      isReferenceCounted() const {return m_referenceCount.isCounted();}


      /**
       * Returns the number of elements in the array.  This is a synonym
       * for getSize().
       *
       * @return The number of elements in the array.
       */
      size_t
      length() const {return this->getSize();}


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
       * Changes the shape of the array and reallocates storage.
       * The current array contents are lost.
       *
       * @param size Number of elements in the array after reallocation.
       */
      void
      reinit(size_t arraySize);


      /**
       * Behaves just like member function reinit() unless
       * this->size() matches the requested size, in which case this
       * member function does nothing.
       *
       * @param size Required number of elements in the array.
       */
      inline void
      reinitIfNecessary(size_t arraySize);


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
       * Returns the number of elements in the array.  This is a synonym
       * for getSize().
       *
       * @return The number of elements in the array.
       */
      size_t
      size() const {return this->getSize();}


      /**
       * Assignment operator shallow copies the contents of source.
       * After the copy, both arrays reference the same data.
       *
       * @param source The Array1D instance to be copied.
       *
       * @return Reference to *this.
       */
      Array1D<Type>&
      operator=(const Array1D<Type>& source);


      /**
       * Assign value to every element in the array.
       *
       * @param value The value to be copied.
       *
       * @return Reference to *this.
       */
      Array1D<Type>&
      operator=(Type value);


      /**
       * Returns the (index)th element of the array by reference.
       *
       * @param index Indicates the selected element.
       *
       * @return Reference to the (index)th element of the array.
       */
      Type& operator()(size_t index) {
        this->checkBounds(index);
        return m_dataPtr[index];
      }


      /**
       * Returns the (index)th element of the array by value.
       *
       * @param index Indicates the selected element.
       *
       * @return Value of the (index)th element of the array.
       */
      Type const& operator()(size_t index) const {
        this->checkBounds(index);
        return m_dataPtr[index];
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
      Type const& operator[](size_t index) const {
        return this->operator()(index);
      }


      /**
       * Increments each element of *this by the value of the
       * corresponding element of arg.
       *
       * @param arg Array1D of values to be added to the elements of
       * *this.
       *
       * @exception ValueException on incompatible array sizes
       *
       * @return Reference to *this.
       */
      template <class Type2>
      Array1D<Type>&
      operator+=(const Array1D<Type2>& arg);


      /**
       * Increments each element of *this by a constant.
       *
       * @param arg Value by which array elements will be incremented.
       *
       * @return Reference to *this.
       */
      Array1D<Type>&
      operator+=(const Type arg);


      /**
       * Decrements each element of *this by the value of the
       * corresponding element of arg.
       *
       * @param arg Array1D of values to be subtracted from the elements
       * of *this.
       *
       * @exception ValueException on incompatible array sizes
       *
       * @return Reference to *this.
       */
      template <class Type2>
      Array1D<Type>&
      operator-=(const Array1D<Type2>& arg);


      /**
       * Decrements each element of *this by a constant.
       *
       * @param arg Value by which array elements will be decremented.
       *
       * @return Reference to *this.
       */
      Array1D<Type>&
      operator-=(const Type arg);


      /**
       * Multiplies each element of *this by the value of the
       * corresponding element of arg.
       *
       * @param arg Array1D of values by which the elements of *this
       * should be multiplied.
       *
       * @exception ValueException on incompatible array sizes
       *
       * @return Reference to *this.
       */
      template <class Type2>
      Array1D<Type>&
      operator*=(const Array1D<Type2>& arg);


      /**
       * Multiplies each element of *this by a constant.
       *
       * @param arg Value by which array elements will be multiplied.
       *
       * @return Reference to *this.
       */
      Array1D<Type>&
      operator*=(const Type arg);


      /**
       * Divides each element of *this by the value of the corresponding
       * element of arg.
       *
       * @param arg Array1D of values by which the elements of *this
       * should be divided.
       *
       * @exception ValueException on incompatible array sizes
       *
       * @return Reference to *this.
       */
      template <class Type2>
      Array1D<Type>&
      operator/=(const Array1D<Type2>& arg);


      /**
       * Divides each element of *this by a constant.
       *
       * @param arg Value by which array elements will be divided.
       *
       * @return Reference to *this.
       */
      Array1D<Type>&
      operator/=(const Type arg);


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
      Array1D<Type>&
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
      Array1D<Type>&
      operator>>=(int numberOfBits);


    private:
      /* ======== Private member functions ======== */

      /**
       * Allocate memory for array data and initialize reference count.
       */
      void
      allocate(size_t arraySize);

      /**
       * Optionally (if and only if BRICK_NUMERIC_CHECKBOUNDS is
       * defined) throw an exception if index is beyond the range of
       * this array.
       *
       * @param index The index to check.
       */
      inline void
      checkBounds(size_t index) const;


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
       * Array1D should start.
       */
      inline static const std::string& ioIntro();

      /**
       * Static constant describing how the string representation of an
       * Array1D should end.
       */
      inline static const char& ioOutro();

      /**
       * Static constant describing how the the data portion of the
       * string representation of an Array1D should start.
       */
      inline static const char& ioOpening();

      /**
       * Static constant describing how the the data portion of the
       * string representation of an Array1D should end.
       */
      inline static const char& ioClosing();

      /**
       * Static constant describing how individual elements should be
       * separated in the string representation of Array1D.
       */
      inline static const char& ioSeparator();


      /* ======== Private data members ======== */
      size_t m_size;
      Type* m_dataPtr;
      common::ReferenceCount m_referenceCount;
    };

    /* Non-member functions which should maybe wind up in a different file */

    /**
     * Elementwise addition of Array1D instances.
     *
     * @param array0 First argument for addition.
     *
     * @param array1 Second argument for addition.
     *
     * @exception ValueException on incompatible array sizes
     *
     * @return Array1D instance in which the value of each element is
     * the sum of the values of the corresponding elements of the two
     * Array1D arguments.
     */
    template <class Type>
    Array1D<Type>
    operator+(const Array1D<Type>& array0, const Array1D<Type>& array1);


    /**
     * Elementwise subtraction of Array1D instances.
     *
     * @param array0 First argument for subtraction.
     *
     * @param array1 Second argument for subtraction.
     *
     * @exception ValueException on incompatible array sizes
     *
     * @return Array1D instance in which the value of each element is
     * the difference of the values of the corresponding elements of the two
     * Array1D arguments.
     */
    template <class Type>
    Array1D<Type>
    operator-(const Array1D<Type>& array0, const Array1D<Type>& array1);


    /**
     * Elementwise multiplication of Array1D instances.
     *
     * @param array0 First argument for multiplication.
     *
     * @param array1 Second argument for multiplication.
     *
     * @exception ValueException on incompatible array sizes
     *
     * @return Array1D instance in which the value of each element is
     * the product of the values of the corresponding elements of the two
     * Array1D arguments.
     */
    template <class Type>
    Array1D<Type>
    operator*(const Array1D<Type>& array0, const Array1D<Type>& array1);


    /**
     * Elementwise division of Array1D instances.
     *
     * @param array0 First argument for division.
     *
     * @param array1 Second argument for division.
     *
     * @exception ValueException on incompatible array sizes
     *
     * @return Array1D instance in which the value of each element is
     * the dividend of the values of the corresponding elements of the two
     * Array1D arguments.
     */
    template <class Type>
    Array1D<Type>
    operator/(const Array1D<Type>& array0, const Array1D<Type>& array1);


    /**
     * Addition of Array1D and scalar.
     *
     * @param array Array1D argument of the addition.
     *
     * @param scalar Scalar argument of the addition.
     *
     * @return Array1D instance in which the value of each element is
     * the sum of the corresponding element of the Array1D argument and
     * the scalar argument.
     */
    template <class Type>
    Array1D<Type> operator+(const Array1D<Type>& array, Type scalar);


    /**
     * Subtraction of Array1D and scalar.
     *
     * @param array0 Array1D argument of the subtraction.
     *
     * @param scalar Scalar argument of the subtraction.
     *
     * @return Array1D instance in which the value of each element is
     * the difference of the corresponding element of the Array1D
     * argument and the scalar argument.
     */
    template <class Type>
    Array1D<Type> operator-(const Array1D<Type>& array0, Type scalar);


    /**
     * Multiplication of Array1D and scalar.
     *
     * @param array0 Array1D argument of the multiplication.
     *
     * @param scalar Scalar argument of the multiplication.
     *
     * @return Array1D instance in which the value of each element is
     * the product of the corresponding element of the Array1D argument
     * and the scalar argument.
     */
    template <class Type>
    Array1D<Type> operator*(const Array1D<Type>& array0, Type scalar);


    /**
     * Division of Array1D and scalar.
     *
     * @param array0 Array1D argument of the division.
     *
     * @param scalar Scalar argument of the division.
     *
     * @return Array1D instance in which the value of each element is
     * the difference of the corresponding element of the Array1D
     * argument and the scalar argument.
     */
    template <class Type>
    Array1D<Type> operator/(const Array1D<Type>& array0, Type scalar);


    /**
     * Addition of scalar and Array1D.
     *
     * @param scalar Scalar argument of the addition.
     *
     * @param array0 Array1D argument of the addition.
     *
     * @return Array1D instance in which the value of each element is
     * the sum of the scalar argument and the corresponding element of
     * the Array1D argument.
     */
    template <class Type>
    inline Array1D<Type> operator+(Type scalar, const Array1D<Type>& array0);


    /**
     * Subtraction of scalar and Array1D.
     *
     * @param scalar Scalar argument of the subtraction.
     *
     * @param array0 Array1D argument of the subtraction.
     *
     * @return Array1D instance in which the value of each element is
     * the difference of the scalar argument and the corresponding
     * element of the Array1D argument.
     */
    template <class Type>
    inline Array1D<Type> operator-(Type scalar, const Array1D<Type>& array0);


    /**
     * Multiplication of scalar and Array1D.
     *
     * @param scalar Scalar argument of the multiplication.
     *
     * @param array0 Array1D argument of the multiplication.
     *
     * @return Array1D instance in which the value of each element is
     * the product of the scalar argument and the corresponding element
     * of the Array1D argument.
     */
    template <class Type>
    inline Array1D<Type> operator*(Type scalar, const Array1D<Type>& array0);


    /**
     * Division of scalar and Array1D.
     *
     * @param scalar Scalar argument of the division.
     *
     * @param array0 Array1D argument of the division.
     *
     * @return Array1D instance in which the value of each element is
     * the dividend of the scalar argument and the corresponding element
     * of the Array1D argument.
     */
    template <class Type>
    inline Array1D<Type> operator/(Type scalar, const Array1D<Type>& array0);


    /**
     * Elementwise comparison of an Array1D with a constant.
     *
     * @param array0 An Array1D instance.
     *
     * @param arg Value to which the elements of array0 should be
     * compared.
     *
     * @return An Array1D<bool> in which each element has value "true"
     * if the corresponding element of array0 is equal to arg.
     */
    template <class Type>
    Array1D<bool>
    operator==(const Array1D<Type>& array0, const Type arg);


    /**
     * Elementwise comparison of an Array1D with another array.
     *
     * @param array0 An Array1D instance.
     *
     * @param array1 A second Array1D instance with the same size as
     * array0.
     *
     * @return An Array1D<bool> in which each element has value "true"
     * if the corresponding element of array0 is equal to the
     * corresponding element of array1.
     */
    template <class Type>
    Array1D<bool>
    operator==(const Array1D<Type>& array0, const Array1D<Type>& array1);


    /**
     * Elementwise comparison of an Array1D with a constant.
     *
     * @param array0 An Array1D instance.
     *
     * @param arg Value to which array elements should be compared.
     *
     * @return An Array1D<bool> in which each element has value "true"
     * if the corresponding element of array0 is greater than arg.
     */
    template <class Type>
    Array1D<bool>
    operator>(const Array1D<Type>& array0, const Type arg);


    /**
     * Elementwise comparison of an Array1D with a constant.
     *
     * @param array0 An Array1D instance.
     *
     * @param arg Value to which array elements should be compared.
     *
     * @return An Array1D<bool> in which each element has value "true"
     * if the corresponding element of array0 is greater than or equal
     * to arg.
     */
    template <class Type>
    Array1D<bool>
    operator>=(const Array1D<Type>& array0, const Type arg);


    /**
     * Elementwise comparison of an Array1D with a constant.
     *
     * @param array0 An Array1D instance.
     *
     * @param arg Value to which array elements should be compared.
     *
     * @return An Array1D<bool> in which each element has value "true"
     * if the corresponding element of array0 is less than arg.
     */
    template <class Type>
    Array1D<bool>
    operator<(const Array1D<Type>& array0, const Type arg);


    /**
     * Elementwise comparison of an Array1D with a constant.
     *
     * @param array0 An Array1D instance.
     *
     * @param arg Value to which array elements should be compared.
     *
     * @return An Array1D<bool> in which each element has value "true"
     * if the corresponding element of array0 is less than or equal to
     * arg.
     */
    template <class Type>
    Array1D<bool>
    operator<=(const Array1D<Type>& array0, const Type arg);


    /**
     * Outputs a text representation of an Array1D instance to a
     * std::ostream.  The output format looks like this:
     *
     * Array1D([1, 2, 4, 8, 16])
     *
     * Where the array elements are output using
     * operator<<(std::ostream&, const Type&)
     * and each element is separated from its neighbors by a comma and
     * whitespace.
     *
     * @param stream Reference to the the output stream.
     *
     * @param array0 const Reference to the Array1D to be output.
     *
     * @exception IOException on invalid stream.
     *
     * @return Reference to output stream.
     */
    template <class Type>
    std::ostream& operator<<(std::ostream& stream, const Array1D<Type>& array0);

    /**
     * Sets the value of an Array1D instance from a std::istream.
     * The input format is as described for
     * operator<<(std::ostream&, const Array1D<Type>&) above.
     *
     * @param stream Reference to the the input stream.
     *
     * @param array0 Reference to the Array1D which will take the input.
     *
     * @return Reference to input stream.
     */
    template <class Type>
    std::istream&
    operator>>(std::istream& stream, Array1D<Type>& array0);

  } // namespace numeric

} // namespace brick

// Include file containing definitions of inline and template
// functions.
#include <brick/numeric/array1D_impl.hh>

#endif /* #ifdef BRICK_NUMERIC_ARRAY1D_HH */
