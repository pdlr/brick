/**
***************************************************************************
* @file brick/numeric/array3D.hh
*
* Header file declaring Array3D class template.
*
* Copyright (C) 2001-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_ARRAY3D_HH
#define BRICK_NUMERIC_ARRAY3D_HH

#include <iostream>
#include <brick/numeric/array1D.hh>
#include <brick/numeric/array2D.hh>
#include <brick/common/exception.hh>

namespace brick {

  namespace numeric {
    
    /**
     ** The Array3D class template represents a 3D array of arbitrary type.
     ** This class has internal reference counting.
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
    template <class Type>
    class Array3D {
    public:
      /* ******** Public typedefs ******** */

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

      /* ******** Public member functions ******** */

      /** 
       * Default constructor initializes to zero size.
       */
      Array3D();

      /** 
       * Construct an arrayShape0 x arrayShape1 x arrayShape2 array.
       * The dimensions are arranged from most minor to most major.
       * If you like to think of the 3D array as being composed of 2D
       * slices, you might want to create the array like this:
       *
       *   Array3D<Type> myArray(numSlices, numRows, numColumns);
       *
       * In this way, stepping through the allocated data from beginning
       * to end will traverse the first row of the first slice, followed
       * subsequent rows of the first slice, followed by the first row
       * of the second slice, and so on.
       * 
       * @param arrayShape0 Number of elements in the first dimension of the array.
       * @param arrayShape1 Number of elements in the second dimension of the array.
       * @param arrayShape2 Number of elements in the third dimension of the array.
       */
      Array3D(size_t arrayShape0, size_t arrayShape1, size_t arrayShape2);

      /**
       * Construct from an initialization string.  The idea is to make
       * constructor very flexible about interpreting string syntax.
       * For now, though, the input string must have the format:
       *
       * "[[[#, #, #, ...], [#, #, #, ...], [#, #, #, ...], ...],
       *   [[#, #, #, ...], [#, #, #, ...], [#, #, #, ...], ...],
       *   [[#, #, #, ...], [#, #, #, ...], [#, #, #, ...], ...], ...]"
       *
       * Where "#" indicates text which can be converted to the element
       * type of the array using the stream input operator.
       * 
       * @param inputString The argument specifies the string from which
       * the array will be constructed.
       */
      explicit
      Array3D(const std::string& inputString);

      /**
       * The copy constructor does a shallow copy.  The newly created
       * array points to the same data as copied array.
       * 
       * @param source Array which is to be copied.
       */
      Array3D(const Array3D<Type> &source);

      /**
       * Construct an array around external data.  Arrays constructed
       * in this way will not implement reference counting, and will
       * not delete dataPtr when done.  The elements of the array are 
       * traversed with axis 0 being the most major axis and axis 2 being
       * the most minor axis.  That is, the first element in the C style
       * array is at index (0, 0, 0), the second is at (0, 0, 1), followed
       * by (0, 0, 2), (0, 0, 3), ..., (0, 0, N), (0, 1, 0), (0, 1, 1),
       * and so on.
       *
       * @param arrayShape0 Number of elements in the first dimension of the array.
       * @param arrayShape1 Number of elements in the second dimension of the array.
       * @param arrayShape2 Number of elements in the third dimension of the array.
       * @param dataPtr A C-style array of Type into which the newly
       * constructed Array3D should index.
       */
      Array3D(size_t arrayShape0, size_t arrayShape1, size_t arrayShape2, Type* const dataPtr);

      /**
       ** Destroys the Array3D instance and deletes the internal data
       ** store if no remaining arrays point to it.
       **/
      virtual
      ~Array3D();

      /** 
       * Return begin() iterator for Standard Library algorithms.
       * 
       * @return Iterator pointing to the first element of the Array3D.
       */
      iterator
      begin() {return m_dataPtr;}

      /** 
       * Return begin() const_iterator for Standard Library algorithms.
       * 
       * @return Const iterator pointing to the first element of the
       * array3D.
       */
      const_iterator
      begin() const {return m_dataPtr;}

      /** 
       * Reset the array to zero size, abandoning all contents.  This is
       * equivalent to this->reinit(0, 0, 0);
       */
      void
      clear() {this->reinit(0, 0, 0);}


      /** 
       * Optionally throw an exception if the shape of *this is
       * different than specified.
       *
       * @param arrayShape0 The required size along the first dimension.
       *
       * @param arrayShape1 The required size along the second dimension.
       *
       * @param arrayShape2 The required size along the third dimension.
       *
       */
      inline void 
      checkDimension(size_t arrayShape0, size_t arrayShape1, size_t arrayShape2) const;


      /** 
       * Allocates a new array and deep copies the contents of *this.
       * 
       * @return A new array which is a (deep) copy of *this.
       */
      Array3D<Type>
      copy() const;

      /** 
       * Deep copies the contents of source.  It is an error if source
       * does not have the same size as *this.  Note that source need
       * not have exactly the same number of rows and columns as *this,
       * as long as their product is right.
       *
       * @param source The array to be copied.
       */
      template <class Type2>
      void
      copy(const Array3D<Type2>& source);

      /**
       * Copies elements from dataPtr.  There must be valid data at all
       * addresses from dataPtr to (dataPtr + this->size());
       * 
       * @param dataPtr Pointer to the data to be copied.
       */
      template <class Type2>
      void
      copy(const Type2* dataPtr);

      /**
       * Returns a pointer to the internal data store.  This is ugly
       * but often necessary for interfacing with external libraries.
       * Data is stored contiguously in raster order.  The first element
       * corresponds to indices (0, 0, 0), the second to (0, 0, 1), the
       * (shape2 + 1)th element corresponds to index (0, 1, 0), and so on.
       * 
       * @return Pointer to the internal data store.
       */
      Type*
      data() {return m_dataPtr;}

      /**
       * This version of data(void) is appropriate for const Array2D,
       * and returns a pointer-to-const.
       * 
       * @return Const pointer to the internal data store.
       */
      const Type*
      data() const {return m_dataPtr;}

      /**
       * Just like data(void), which is documented above, but returns a
       * pointer to the (index)th element instead of the first element.
       * Note that "data(index)" is synonymous with "data() + index" .
       *
       * @param index Indicates to which element of the array the return
       * value should point.
       * @return pointer to the (index)th element of the array.
       */
      Type*
      data(size_t index) {
        checkBounds(index);
        return m_dataPtr + index;
      }

      /**
       * This version of data(size_t) is appropriate for const Array1D,
       * and returns a pointer-to-const.
       * 
       * @param index Indicates to which element of the array the return
       * value should point.
       * @return Const pointer to the (index)th element of the array.
       */
      const Type*
      data(size_t index) const {
        checkBounds(index);
        return m_dataPtr+index;
      }

      /**
       * Just like data(void), but returns a pointer to the element indexed by
       * (index0, index1, index2).
       *
       * @param index0 First index of selected element.
       * @param index1 Second index of selected element.
       * @param index2 Third index of selected element.
       * @return Pointer to the selected element of the array.
       */
      Type*
      data(size_t index0, size_t index1, size_t index2) {
        checkBounds(index0, index1, index2);
        return (m_dataPtr + index2 + (index1 * m_shape2)
                + (index0 * m_shape1Times2));
      }

      /**
       * This version of data(size_t, size_t, size_t) is appropriate for
       * const Array3D, and returns a pointer-to-const.
       *
       * @param index0 First index of selected element.
       * @param index1 Second index of selected element.
       * @param index2 Third index of selected element.
       * @return Pointer to the selected element of the array.
       */
      const Type*
      data(size_t index0, size_t index1, size_t index2) const {
        checkBounds(index0, index1, index2);
        return (m_dataPtr + index2 + (index1 * m_shape2)
                + (index0 * m_shape1Times2));
      }


      /** 
       * This member function returns true if the array instance
       * contains no elements.  It has complexity O(1).
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
       * the array3D.
       */
      iterator
      end() {return m_dataPtr + m_size;}

      /** 
       * Return end() const_iterator for Standard Library algorithms.
       * 
       * @return Const iterator pointing just past the last element of
       * the array3D.
       */
      const_iterator
      end() const {return m_dataPtr + m_size;}


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
       * @param index0 This argument, and the next two, specify which
       * element value should be returned.
       * 
       * @param index1 This argument, the previous, and the next,
       * specify which element value should be returned.
       * 
       * @param index2 This argument, and the previous two, specify
       * which element value should be returned.
       * 
       * @return The return value is a copy of the requested element.
       */
      Type
      getElement(size_t index0, size_t index1, size_t index2) const {
        return this->operator()(index0, index1, index2);
      }        


      /** 
       * This member function sets the value of the array from an input
       * stream.  The array is modified only if the read was successful,
       * otherwise the array is not modified, and failbit is set in the
       * stream state.  Because of this nice behavior, readFromStream is
       * quite slow.
       * 
       * @param inputStream This is the stream from which to read the
       * array.
       * @return The return value is a reference to inputStream.
       */
      std::istream&
      readFromStream(std::istream& inputStream);
    
      /** 
       * Changes the shape of the array and reallocates storage.  The
       * current array contents are lost.  After a successful call to
       * reinit(), the array will have arrayShape0 x arrayShape1 x arrayShape2
       * elements.
       * 
       * @param arrayShape0 Requested dimension along axis 0.
       * @param arrayShape1 Requested dimension along axis 1.
       * @param arrayShape2 Requested dimension along axis 2.
       */
      void
      reinit(size_t arrayShape0, size_t arrayShape1, size_t arrayShape2);


      /** 
       * Behaves just like member function reinit() unless *this
       * already has the requested shape, in which case this member
       * function does nothing, or *this has the requested size (but a
       * different shape), in which case this member function behaves
       * just like member function resize().
       * 
       * @param arrayShape0 Requested dimension along axis 0.
       * @param arrayShape1 Requested dimension along axis 1.
       * @param arrayShape2 Requested dimension along axis 2.
       */
      inline void
      reinitIfNecessary(size_t arrayShape0, size_t arrayShape1,
                        size_t arrayShape2);

      
      /**
       * Changes the shape of the array.  It is an error if arrayShape0 *
       * arrayShape1 * arrayShape2 is not equal to this->size().  Any one argument
       * may be specified as -1, but no more than one.  If one of the
       * arguments is -1, its value will be calculated based on the total
       * size of the array and the value of the other arguments.
       * 
       * @param arrayShape0 Requested dimension along axis 0.
       * @param arrayShape1 Requested dimension along axis 1.
       * @param arrayShape2 Requested dimension along axis 2.
       */
      void
      reshape(int arrayShape0, int arrayShape1, int arrayShape2);


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
       * @param index0 This argument, and the next two, specify which
       * element value should be set.
       * 
       * @param index1 This argument, the previous, and the next,
       * specify which element value should be set.
       * 
       * @param index2 This argument, and the previous two, specify
       * which element value should be set.
       * 
       * @param value This argument will be copied into the selected
       * array element.
       * 
       * @return The return value is a reference to the newly created
       * array element.
       */
      Type&
      setElement(size_t index0, size_t index1, size_t index2,
                 const Type& value) {
        return this->operator()(index0, index1, index2) = value;
      }


      /** 
       * Returns a 3 element Array1D containing the dimension of *this
       * along each axis.  That is, the first element of the returned
       * array will be this->shape0(), the second element will be
       * this->shape1(), and the third element will be this->shape2().
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
       * @return Dimension along the selected axis.
       */
      size_t
      shape(size_t axis) const;

      /** 
       * Returns the array dimension along axis 0.  This is synonymous
       * with shape(0), but may execute faster.
       * 
       * @return Number of elements.
       */
      size_t
      shape0() const {return m_shape0;}

      /** 
       * Returns the array dimension along axis 1.  This is synonymous
       * with shape(1), but may execute faster.
       * 
       * @return Number of elements.
       */
      size_t
      shape1() const {return m_shape1;}

      /** 
       * Returns the array dimension along axis 2.  This is synonymous
       * with shape(2), but may execute faster.
       * 
       * @return Number of elements.
       */
      size_t
      shape2() const {return m_shape2;}

      /** 
       * Returns the number of elements in the array.  This is the
       * equal to the product of the elements of this->shape().
       * 
       * @return Number of elements.
       */
      size_t
      size() const {return m_size;}
  
      /**
       * Returns an Array2D<Type> which addresses an entire (rows x
       * columns) slice of *this.  The returned Array2D references the
       * same data as the selected slice of *this.
       * 
       * @param index Specifies which slice to reference.
       * @return An Array2D instance referencing the selected slice.
       */
      Array2D<Type>
      slice(size_t index);

      /**
       * Returns a const Array2D<Type> which addresses an entire (rows x
       * columns) slice of *this.  The returned Array2D references the
       * same data as the selected slice of *this.
       * 
       * @param index Specifies which slice to reference.
       * @return An Array2D instance referencing the selected slice.
       */
      const Array2D<Type>
      slice(size_t index) const;

      /** 
       * Assignment operator shallow copies the contents of source.
       * After the copy, both arrays reference the same data.
       * 
       * @param source The Array3D instance to be copied.
       * @return Reference to *this.
       */
      Array3D<Type>&
      operator=(const Array3D<Type>& source);

      /** 
       * Assign value to every element in the array.
       * 
       * @param value The value to be copied.
       * @return Reference to *this.
       */
      Array3D<Type>&
      operator=(Type value);

      /** 
       * Returns the (index)th element of the array by reference.
       * Synonymous with operator()(size_t).
       * 
       * @param index Indicates the selected element.
       * @return Reference to the (index)th element of the array.
       */
      inline Type&
      operator()(size_t index) {
        checkBounds(index);
        return m_dataPtr[index];
      }

      /** 
       * Returns the (index)th element of the array by value.
       * Synonymous with operator()(size_t) const.
       * 
       * @param index Indicates the selected element.
       * @return Value of the (index)th element of the array.
       */
      inline Type
      operator()(size_t index) const {
        checkBounds(index);
        return m_dataPtr[index];
      }
    
      /** 
       * Returns one element of the array by reference.  Elements are
       * indexed as described in the documentation of 
       * Array3D::Array3D(size_t, size_t, size_t).
       *
       * @param index0 Indicates the position of the selected element
       * along the first axis.
       * @param index1 Indicates the position of the selected element
       * along the third axis.
       * @param index2 Indicates the position of the selected element
       * along the third axis.
       * @return Reference to the selected element of the array.
       */
      inline Type&
      operator()(size_t index0, size_t index1, size_t index2) {
        checkBounds(index0, index1, index2);
        return m_dataPtr[index2 + (index1 * m_shape2)
                         + (index0 * m_shape1Times2)];
      }

      /** 
       * Returns one element of the array by value.  Elements are
       * indexed as described in the documentation of 
       * Array3D::Array3D(size_t, size_t, size_t).
       *
       * @param index0 Indicates the position of the selected element
       * along the first axis.
       * @param index1 Indicates the position of the selected element
       * along the third axis.
       * @param index2 Indicates the position of the selected element
       * along the third axis.
       * @return Value of the selected element of the array.
       */
      inline Type
      operator()(size_t index0, size_t index1, size_t index2) const {
        checkBounds(index0, index1, index2);
        return m_dataPtr[index2 + (index1 * m_shape2)
                         + (index0 * m_shape1Times2)];
      }

      /** 
       * Returns the (index)th element of the array by reference.
       * Synonymous with operator()(size_t).
       * 
       * @param index Indicates the selected element.
       * @return Reference to the (index)th element of the array.
       */
      inline Type&
      operator[](size_t index) {return this->operator()(index);}

      /** 
       * Returns the (index)th element of the array by value.
       * Synonymous with operator()(size_t) const.
       * 
       * @param index Indicates the selected element.
       * @return Value of the (index)th element of the array.
       */
      inline Type
      operator[](size_t index) const {return this->operator()(index);}

      /** 
       * Multiplies each element of *this by the value of the
       * corresponding element of arg.
       * 
       * @param arg Array3D of values to be added to the elements of
       * *this.     
       * @return Reference to *this.
       */
      template <class Type2>
      Array3D<Type>&
      operator*=(const Array3D<Type2>& arg);

      /** 
       * Divides each element of *this by the value of the
       * corresponding element of arg.
       * 
       * @param arg Array3D of values to be added to the elements of
       * *this.     
       * @return Reference to *this.
       */
      template <class Type2>
      Array3D<Type>&
      operator/=(const Array3D<Type2>& arg);

      /** 
       * Increments each element of *this by the value of the
       * corresponding element of arg.
       * 
       * @param arg Array3D of values to be added to the elements of
       * *this.     
       * @return Reference to *this.
       */
      template <class Type2>
      Array3D<Type>&
      operator+=(const Array3D<Type2>& arg);

      /** 
       * Decrements each element of *this by the value of the
       * corresponding element of arg.
       * 
       * @param arg Array3D of values to be subtracted from the elements
       * of *this.
       * @return Reference to *this.
       */
      template <class Type2>
      Array3D<Type>&
      operator-=(const Array3D<Type2>& arg);

      /** 
       * Increments each element of *this by a constant.
       * 
       * @param arg Value by which array elements will be incremented.
       * @return Reference to *this.
       */
      Array3D<Type>&
      operator+=(const Type arg);

      /** 
       * Decrements each element of *this by a constant.
       * 
       * @param arg Value by which array elements will be decremented.
       * @return Reference to *this.
       */
      Array3D<Type>&
      operator-=(const Type arg);

      /** 
       * Multiplies each element of *this by a constant.
       * 
       * @param arg Value by which array elements will be multiplied.
       * @return Reference to *this.
       */
      Array3D<Type>&
      operator*=(const Type arg);

      /** 
       * Divides each element of *this by a constant.
       * 
       * @param arg Value by which array elements will be divided.
       * @return Reference to *this.
       */
      Array3D<Type>&
      operator/=(const Type arg);

    private:
      // Private member functions
      void allocate();

      // Optionally throw an exception if index is beyond the range of
      // this array.
      inline void
      checkBounds(size_t index) const;

      // Optionally throw an exception if any index is beyond the range of
      // this array.
      inline void
      checkBounds(size_t index0, size_t index1, size_t index2) const;

      void deAllocate();

      // Constants to help with formatting.  We use the initialization
      // on first use paradigm for the string constants to avoid
      // headaches.

      /**
       * Static constant describing how the string representation of an
       * Array3D should start.
       */
      static const std::string& ioIntro(); 

      /**
       * Static constant describing how the string representation of an
       * Array3D should end.
       */
      static const std::string& ioOutro();

      /**
       * Static constant describing how the the data portion of the
       * string representation of an Array3D should start.
       */
      static const char ioOpening = '[';

      /**
       * Static constant describing how the the data portion of the
       * string representation of an Array3D should end.
       */
      static const char ioClosing = ']';

      /**
       * Static constant describing how individual elements should be
       * separated in the string representation of Array3D.
       */
      static const char ioSeparator = ',';
    
      // Private member variables
      size_t m_shape0;
      size_t m_shape1;
      size_t m_shape2;
      size_t m_shape1Times2;
      size_t m_size;
      Type* m_dataPtr;
      size_t* m_refCountPtr;
      bool m_isAllocated;
    
    };


    /* =========== Non-member functions related to Array3D =========== */

//   /** 
//    * This function returns the maximum element of the an Array3D
//    * instance.
//    * 
//    * @param array This argument is the Array3D instance in which to
//    * search for the largest element.
//    * 
//    * @return The return value is the value of the largest element.
//    */
//   template <class Type> 
//   Type maximum(const Array3D<Type>& array);


//   /** 
//    * This function returns the minimum element of the an Array3D
//    * instance.
//    * 
//    * @param array This argument is the Array3D instance in which to
//    * search for the smallest element.
//    * 
//    * @return The return value is the value of the smallest element.
//    */
//   template <class Type> 
//   Type minimum(const Array3D<Type>& array);


//   /** 
//    * This function computes the sum of the elements of its argument.
//    * The summation is accumulated into a variable of type
//    * NumericTraits<Type>::SumType, which for now defaults to Type, but
//    * ultimately should be a type which (for fixed point and integral
//    * types) has enough bits to hold the sum of at least 65535
//    * elements.
//    * 
//    * @param array0 This argument is the array to be summed.
//    *
//    * @return The summation of all the elements of array0.
//    */
//   template <class Type> 
//   Type sum(const Array3D<Type>& array);


  
    /** 
     * Elementwise addition of Array3D instances.
     * 
     * @param array0 First argument for addition.
     *
     * @param array1 Second argument for addition.
     *
     * @exception ValueException on incompatible array sizes
     *
     * @return Array3D instance in which the value of each element is
     * the sum of the values of the corresponding elements of the two
     * Array3D arguments.
     */
    template <class Type>
    Array3D<Type>
    operator+(const Array3D<Type>& array0, const Array3D<Type>& array1);


    /** 
     * Elementwise subtraction of Array3D instances.
     * 
     * @param array0 First argument for subtraction.
     *
     * @param array1 Second argument for subtraction.
     *
     * @exception ValueException on incompatible array sizes
     *
     * @return Array3D instance in which the value of each element is
     * the difference of the values of the corresponding elements of the two
     * Array3D arguments.
     */
    template <class Type>
    Array3D<Type>
    operator-(const Array3D<Type>& array0, const Array3D<Type>& array1);


    /** 
     * Elementwise multiplication of Array3D instances.
     * 
     * @param array0 First argument for multiplication.
     *
     * @param array1 Second argument for multiplication.
     *
     * @exception ValueException on incompatible array sizes
     *
     * @return Array3D instance in which the value of each element is
     * the product of the values of the corresponding elements of the two
     * Array3D arguments.
     */
    template <class Type>
    Array3D<Type>
    operator*(const Array3D<Type>& array0, const Array3D<Type>& array1);

  
    /** 
     * Elementwise division of Array3D instances.
     * 
     * @param array0 First argument for division.
     *
     * @param array1 Second argument for division.
     *
     * @exception ValueException on incompatible array sizes
     *
     * @return Array3D instance in which the value of each element is
     * the dividend of the values of the corresponding elements of the two
     * Array3D arguments.
     */
    template <class Type>
    Array3D<Type>
    operator/(const Array3D<Type>& array0, const Array3D<Type>& array1);


    /** 
     * Addition of Array3D and scalar.
     * 
     * @param array0 Array3D argument of the addition.
     *
     * @param scalar Scalar argument of the addition.
     *
     * @return Array3D instance in which the value of each element is
     * the sum of the corresponding element of the Array3D argument and
     * the scalar argument.
     */
    template <class Type>
    Array3D<Type>
    operator+(const Array3D<Type>& array0, Type scalar);

  
    /** 
     * Subtraction of Array3D and scalar.
     * 
     * @param array0 Array3D argument of the subtraction.
     *
     * @param scalar Scalar argument of the subtraction.
     *
     * @return Array3D instance in which the value of each element is
     * the difference of the corresponding element of the Array3D
     * argument and the scalar argument.
     */
    template <class Type>
    Array3D<Type>
    operator-(const Array3D<Type>& array0, Type scalar);

  
    /** 
     * Multiplication of Array3D and scalar.
     * 
     * @param array0 Array3D argument of the multiplication.
     *
     * @param scalar Scalar argument of the multiplication.
     *
     * @return Array3D instance in which the value of each element is
     * the product of the corresponding element of the Array3D argument
     * and the scalar argument.
     */
    template <class Type>
    Array3D<Type>
    operator*(const Array3D<Type>& array0, Type scalar);

  
    /** 
     * Division of Array3D and scalar.
     * 
     * @param array0 Array3D argument of the division.
     *
     * @param scalar Scalar argument of the division.
     *
     * @return Array3D instance in which the value of each element is
     * the difference of the corresponding element of the Array3D
     * argument and the scalar argument.
     */
    template <class Type>
    Array3D<Type>
    operator/(const Array3D<Type>& array0, Type scalar);

  
    /** 
     * Addition of scalar and Array3D.
     * 
     * @param scalar Scalar argument of the addition.
     *
     * @param array0 Array3D argument of the addition.
     *
     * @return Array3D instance in which the value of each element is
     * the sum of the scalar argument and the corresponding element of
     * the Array3D argument.
     */
    template <class Type>
    inline Array3D<Type>
    operator+(Type scalar, const Array3D<Type>& array0);

  
    /** 
     * Multiplication of scalar and Array3D.
     * 
     * @param scalar Scalar argument of the multiplication.
     *
     * @param array0 Array3D argument of the multiplication.
     *
     * @return Array3D instance in which the value of each element is
     * the product of the scalar argument and the corresponding element
     * of the Array3D argument.
     */
    template <class Type>
    inline Array3D<Type>
    operator*(Type scalar, const Array3D<Type>& array0);


    /** 
     * Elementwise comparison of an Array3D with a constant.
     *
     * @param array0 An Array3D instance.
     * 
     * @param arg Value to which the elements of array0 should be
     * compared.
     * 
     * @return An Array3D<bool> in which each element has value "true"
     * if the corresponding element of array0 is equal to arg.
     */
    template <class Type>
    Array3D<bool>
    operator==(const Array3D<Type>& array0, const Type arg);

    
    /** 
     * Elementwise comparison of an Array3D with another array.
     * 
     * @param array0 An Array3D instance.
     *
     * @param array1 A second Array3D instance with the same size as
     * array0.
     * 
     * @return An Array3D<bool> in which each element has value "true"
     * if the corresponding element of array0 is equal to the
     * corresponding element of array1.
     */
    template <class Type>
    Array3D<bool>
    operator==(const Array3D<Type>& array0, const Array3D<Type>& array1);
    

    /** 
     * Elementwise comparison of Array2D with a constant.
     * 
     * @param array0 Array2D instance.
     *
     * @param arg Value to which array0 elements should be compared.
     *
     * @return Array2D<bool> in which each element has value "true" if
     * the corresponding element of array0 is less than arg.
     */
    template <class Type>
    Array3D<bool>
    operator<(const Array3D<Type>& array0, Type arg);

  
    /** 
     * Elementwise comparison of Array3D with a constant.
     * 
     * @param array0 Array3D instance.
     *
     * @param arg Value to which array0 elements should be compared.
     *
     * @return Array3D<bool> in which each element has value "true" if
     * the corresponding element of array0 is less than or equal to arg.
     */
    template <class Type>
    Array3D<bool>
    operator<=(const Array3D<Type>& array0, Type arg);


    /** 
     * Elementwise comparison of Array3D with a constant.
     * 
     * @param array0 Array3D instance.
     *
     * @param arg Value to which array0 elements should be compared.
     *
     * @return Array3D<bool> in which each element has value "true" if
     * the corresponding element of array0 is greater than arg.
     */
    template <class Type>
    Array3D<bool>
    operator>(const Array3D<Type>& array0, Type arg);

  
    /** 
     * Elementwise comparison of Array3D with a constant.
     * 
     * @param array0 Array3D instance.
     *
     * @param arg Value to which array0 elements should be compared.
     *
     * @return Array3D<bool> in which each element has value "true" if
     * the corresponding element of array0 is greater than or equal to arg.
     */
    template <class Type>
    Array3D<bool>
    operator>=(const Array3D<Type>& array0, Type arg);


    /** 
     * This operator outputs a text representation of an Array3D
     * instance to a std::ostream.  The output format looks like this:
     *
     * Array3D([[[1, 2, 4, 8, 16],
     *           [5, 1, 6, 7, 2]],
     *          [[8, 2, 4, 1, 0],
     *           [3, 3, 4, 5, 1]]])
     *
     * Where the array elements are output using
     * operator<<(std::ostream&, const Type&)
     * and each element is separated from its neighbors by a comma and
     * whitespace.
     * 
     * @param stream Reference to the the output stream.
     *
     * @param array0 Const reference to the Array3D to be output.
     *
     * @return Reference to output stream.
     */
    template <class Type>
    std::ostream& operator<<(std::ostream& stream, const Array3D<Type>& array0);

  
    /** 
     * This operator sets the value of an Array3D instance from a
     * std::istream.  The input format is as described for
     * operator<<(std::ostream&, const Array3D<Type>&) above.
     * 
     * @param stream Reference to the the input stream.
     *
     * @param array0 Reference to the Array3D which will take the input.
     *
     * @return Reference to the input stream.
     */
    template <class Type>
    std::istream&
    operator>>(std::istream& stream, Array3D<Type>& array0);
  
  } // namespace numeric

} // namespace brick

#include <brick/numeric/array3D_impl.hh>

#endif /* #ifdef BRICK_NUMERIC_ARRAY3D_HH */
