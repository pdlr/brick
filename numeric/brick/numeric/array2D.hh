/**
***************************************************************************
* @file brick/numeric/array2D.hh
*
* Header file declaring Array2D class.
*
* Copyright (C) 2001-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_ARRAY2D_HH
#define BRICK_ARRAY2D_HH

// Early inclusion of <algorithm> to allow inlining of member templates.
#include <algorithm>

#include <iostream>
#include <brick/common/exception.hh>
#include <brick/numeric/array1D.hh>

namespace brick {

  namespace numeric {
    
    /**
     ** The Array2D class template represents a 2D array of arbitrary type.
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
     ** array1 and array2 have different shapes.
     **/
    template <class Type>
    class Array2D {
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
      Array2D();


      /** 
       * Constructs a "arrayRows x arrayColumns" element array.
       * 
       * @param arrayRows Number of rows in the array after successful
       * construction.
       *
       * @param arrayColumns Number of columns in the array after successful
       * construction.
       */
      Array2D(size_t arrayRows, size_t arrayColumns);


      /**
       * Construct from an initialization string.  The idea is to make
       * constructor very flexible about interpreting string syntax.
       * For now, though, the input string must have the format:
       *
       * "[[#, #, #, ...], [#, #, #, ...], [#, #, #, ...], ...]"
       *
       * Where "#" indicates text that can be converted to the element
       * type of the array using the stream input operator.
       * 
       * @param inputString The argument specifies the string from which
       * the array will be constructed.
       */
      explicit
      Array2D(const std::string& inputString);

      
      /** 
       * The copy constructor does a shallow copy.  The newly created
       * array points to the same data as copied array.
       * 
       * @param source The Array2D instance to be copied.
       */
      Array2D(const Array2D<Type> &source);

      
      /**
       * Construct an array around external data.  Arrays constructed in
       * this way will not implement reference counting, and will not
       * delete dataPtr when done.  The elements of the Array are
       * organized in row-major order.
       * 
       * @param arrayRows Number of rows in the array after successful
       * construction.
       *
       * @param arrayColumns Number of columns in the array after successful
       * construction.
       *
       * @param dataPtr A C-style array of Type into which the newly
       * constructed Array2D should index.
       */
      Array2D(size_t arrayRows, size_t arrayColumns, Type* const dataPtr);


      /** 
       * Construct an array around external data that was allocated by
       * an Array?D instance.  Arrays constructed in this way _do_
       * implement reference counting, and will delete dataPtr when
       * done.  This constructor is provided primarily so that other
       * dimensionality array classes can return Array2D instances that
       * reference their data without being friend classes.  Caveat
       * emptor.
       * 
       * @param arrayRows This argument specifies the number of rows in the
       * array.
       * 
       * @param arrayColumns This argument specifies the number of columns in
       * the array.
       * 
       * @param dataPtr This argument is a C-style array containing the
       * data to which the new Array2D instance should refer.
       * 
       * @param referenceCount ReferenceCount instance indicating
       * the number of Array classes currently using dataPtr.
       */
      Array2D(size_t arrayRows, size_t arrayColumns, Type* const dataPtr,
              common::ReferenceCount const& referenceCount);

      
      /**
       * Destroys the Array2D instance and deletes the internal data
       * store if no remaining arrays point to it.
       */
      virtual
      ~Array2D();

      
      /** 
       * Return begin() iterator for Standard Library algorithms.
       * 
       * @return Iterator pointing to the first element of the Array2D.
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
       * does not have the same size as *this.  Note that source need
       * not have exactly the same number of rows and columns as *this,
       * as long as their product is right.
       *
       * @param source The array to be copied.
       *
       * @exception ValueException thrown when array sizes do not match.
       */
      template <class Type2> void
      copy(const Array2D<Type2>& source);

      
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
       * stored contiguously in raster order.  The first element
       * corresponds to indices (0, 0), the second to (0, 1), the
       * (columns + 1)th element corresponds to index (1, 0), and so
       * on.
       *
       * @return Pointer to the internal data store.
       */
      Type*
      data() {return this->getData();}

      
      /**
       * This is an alias for member function getData().  This version
       * of data(void) is appropriate for const Array2D, and returns a
       * pointer-to-const.
       * 
       * @return Const pointer to the internal data store.
       */
      const Type*
      data() const {return this->getData();}


      /** 
       * This is an alias for member function getData().  Just like
       * data(void), which is documented above, but returns a pointer
       * to the (index)th element instead of the first element.  Note
       * that "data(index)" is synonymous with "data() + index" .
       * 
       * @param index Indicates to which element of the array the return
       * value should point.
       *
       * @return Pointer to the (index)th element of the array.
       */
      Type*
      data(size_t index) {return this->getData(index);}


      /**
       * This is an alias for member function getData().  This version
       * of data(size_t) is appropriate for const Array2D, and returns
       * a pointer-to-const.
       * 
       * @param index Indicates to which element of the array the return
       * value should point.
       *
       * @return Const pointer to the (index)th element of the array.
       */
      const Type*
      data(size_t index) const {return this->getData(index);}


      /** 
       * This is an alias for member function getData().  Just like
       * data(void), which is documented above, but returns a pointer
       * to the element indexed by (row, column).
       * 
       * @param rowIndex The row of the element to which the return value
       * should point.
       *
       * @param columnIndex The column of the element to which the return
       * value should point.
       *
       * @return Pointer to the selected element of the array.
       */
      Type*
      data(size_t rowIndex, size_t columnIndex) {
        return this->getData(rowIndex, columnIndex);
      }


      /**
       * This is an alias for member function getData().  This version
       * of data(size_t, size_t) is appropriate for const Array2D, and
       * returns a pointer-to-const.  @param rowIndex The row of the
       * element to which the return value should point.
       *
       * @param columnIndex The column of the element to which the return
       * value should point.
       *
       * @return Const pointer to the selected element of the array.
       */
      const Type*
      data(size_t rowIndex, size_t columnIndex) const {
        return this->getData(rowIndex, columnIndex);
      }

    
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
       * Returns the array dimension along the "second" axis (the number
       * of columns in the array.)  This is synonymous with shape(1), but
       * may execute faster.
       * 
       * @return Number of columns.
       */
      size_t
      getColumns() const {return m_columns;}


      /**
       * Returns a pointer to the internal data store.  This is ugly
       * but often necessary for interfacing with external libraries.
       * Data is stored contiguously in raster order.  The first
       * element corresponds to indices (0, 0), the second to (0, 1),
       * the (columns + 1)th element corresponds to index (1, 0), and
       * so on.
       *
       * @return Pointer to the internal data store.
       */
      Type*
      getData() {return m_dataPtr;}

      
      /**
       * This version
       * of getData(void) is appropriate for const Array2D, and returns a
       * pointer-to-const.
       * 
       * @return Const pointer to the internal data store.
       */
      const Type*
      getData() const {return m_dataPtr;}

      /** 
       * Just like
       * getData(void), which is documented above, but returns a pointer
       * to the (index)th element instead of the first element.  Note
       * that "getData(index)" is synonymous with "getData() + index" .
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
       * This version
       * of getData(size_t) is appropriate for const Array2D, and returns
       * a pointer-to-const.
       * 
       * @param index Indicates to which element of the array the return
       * value should point.
       *
       * @return Const pointer to the (index)th element of the array.
       */
      const Type*
      getData(size_t index) const {
        this->checkBounds(index);
        return m_dataPtr+index;
      }

      /** 
       * Just like
       * getData(void), which is documented above, but returns a pointer
       * to the element indexed by (row, column).
       * 
       * @param rowIndex The row of the element to which the return value
       * should point.
       *
       * @param columnIndex The column of the element to which the return
       * value should point.
       *
       * @return Pointer to the selected element of the array.
       */
      Type*
      getData(size_t rowIndex, size_t columnIndex) {
        this->checkBounds(rowIndex, columnIndex);
        return m_dataPtr + columnIndex + (rowIndex * m_columns);
      }

      /**
       * This version
       * of getData(size_t, size_t) is appropriate for const Array2D, and
       * returns a pointer-to-const.  @param rowIndex The row of the
       * element to which the return value should point.
       *
       * @param columnIndex The column of the element to which the return
       * value should point.
       *
       * @return Const pointer to the selected element of the array.
       */
      const Type*
      getData(size_t rowIndex, size_t columnIndex) const {
        this->checkBounds(rowIndex, columnIndex);
        return m_dataPtr + columnIndex + (rowIndex * m_columns);
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
      getElement(size_t rowIndex, size_t columnIndex) const {
        return this->operator()(rowIndex, columnIndex);
      }        


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
       * Returns an Array1D<Type> that addresses an entire row of
       * *this.  The returned Array1D references the same data as the
       * selected row of *this, and is valid only as long as *this is
       * valid.
       * 
       * @param index Specifies which row to reference.
       *
       * @return An Array1D instance referencing the selected row.
       */
      Array1D<Type>
      getRow(size_t index);

      
      /**
       * Returns a const Array1D<Type> that addresses an entire row of
       * *this.  The returned Array1D references the same data as the
       * selected row of *this, and is valid only as long as *this is
       * valid.
       * 
       * @param index Specifies which row to reference.
       *
       * @return An Array1D instance referencing the selected row.
       */
      const Array1D<Type>
      getRow(size_t index) const;

      
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
       * Returns an Array1D, with size equal to this->getSize(), which
       * references the same data as *this.  In other words, ravel()
       * returns a flattened version of *this.
       * 
       * @return Array1D referencing the same data as *this.
       */
      Array1D<Type>
      ravel();

      
      /**
       * Returns a const Array1D, with size equal to this->getSize(),
       * which references the same data as *this.  In other words,
       * ravel() returns a flattened version of *this.
       * 
       * @return Array1D referencing the same data as *this.
       */
      const Array1D<Type>
      ravel() const;

      
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
       * Changes the shape of the array.  It is an error if
       * (arrayRows * arrayColumns) is not equal to this->getSize().  Either
       * arrayRows or arrayColumns may be specified as -1, but not both.
       * If one of the arguments is -1, its value will be calculated based
       * on the total size of the array and the value of the other argument.
       * 
       * @param arrayRows Requested row dimension.
       *
       * @param arrayColumns Requested column dimension.
       *
       * @exception ValueException thrown when old and new array sizes
       * do not match.
       */
      void
      reshape(int arrayRows, int arrayColumns);

      
      /**
       * This function is an alias for member function getRow().
       * 
       * @param index Specifies which row to reference.
       *
       * @return An Array1D instance referencing the selected row.
       */
      inline Array1D<Type>
      row(size_t index) {return this->getRow(index);}

      
      /**
       * This function is an alias for member function getRow() const.
       * 
       * @param index Specifies which row to reference.
       *
       * @return An Array1D instance referencing the selected row.
       */
      const Array1D<Type>
      row(size_t index) const {return this->getRow(index);}

      
      /** 
       * Return a Standard Library style iterator that points to the
       * beginning of the specified row.
       *
       * @param rowIndex Specifies which row if the array should be
       * referenced.  Row 0 is the top row.
       *
       * @return Iterator pointing to the first element of the
       * specified row of the Array2D instance.
       */
      iterator
      rowBegin(size_t rowIndex) {return m_dataPtr + rowIndex * m_columns;}

      
      /** 
       * Return a Standard Library style const iterator that points
       * to the beginning of the specified row.
       *
       * @param rowIndex Specifies which row if the array should be
       * referenced.  Row 0 is the top row.
       *
       * @return Const iterator pointing to the first element of the
       * specified row of the Array2D instance.
       */
      const_iterator
      rowBegin(size_t rowIndex) const {
	return m_dataPtr + rowIndex * m_columns;
      }

      
      /** 
       * Return a Standard Library style iterator that points one
       * element past the end of the specified row.
       *
       * @param rowIndex Specifies which row if the array should be
       * referenced.  Row 0 is the top row.
       *
       * @return Iterator pointing one element past the end of the
       * specified row.
       */
      iterator
      rowEnd(size_t rowIndex) {
	return m_dataPtr + (rowIndex + 1) * m_columns;
      }


      /** 
       * Return a Standard Library style const iterator that points
       * one element past the end of the specified row.
       *
       * @param rowIndex Specifies which row if the array should be
       * referenced.  Row 0 is the top row.
       *
       * @return Const iterator pointing one element past the end of
       * the specified row.
       */
      const_iterator
      rowEnd(size_t rowIndex) const {
	return m_dataPtr + (rowIndex + 1) * m_columns;
      }

      
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
       * Compute matrix transpose.  The resulting array does not
       * reference the same memory as *this.
       * 
       * @return Transposed copy of *this.  
       */
      Array2D<Type>
      transpose() const;

      
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
       * Assign value to every element in the array.
       * 
       * @param value The value to be copied.
       *
       * @return Reference to *this.
       */
      Array2D<Type>&
      operator=(Type value);


      /** 
       * Returns the (index)th element of the array by reference.
       * Elements are indexed in row-major order.
       *
       * @param index Indicates the selected element.
       *
       * @return Reference to the (index)th element of the array.
       */
      Type&
      operator()(size_t index) {
        this->checkBounds(index);
        return m_dataPtr[index];
      }


      /** 
       * Returns the (index)th element of the array by value.  Elements
       * are indexed in row-major order.
       * 
       * @return Value of the (index)th element of the array.
       */
      Type operator()(size_t index) const {
        this->checkBounds(index);
        return m_dataPtr[index];
      }
  

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
        return m_dataPtr[columnIndex + rowIndex * m_columns];
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
        return m_dataPtr[columnIndex + rowIndex * m_columns];
      }


      /** 
       * Returns the (index)th element of the array by reference.
       * Synonymous with operator()(size_t).
       * 
       * @param index Indicates the selected element.
       *
       * @return Reference to the (index)th element of the array.
       */
      Type&
      operator[](size_t index) {return this->operator()(index);}
  

      /** 
       * Returns the (index)th element of the array by value.
       * Synonymous with operator()(size_t) const.
       * 
       * @param index Indicates the selected element.
       *
       * @return Value of the (index)th element of the array.
       */
      Type
      operator[](size_t index) const {return this->operator()(index);}


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

    private:
      /* ******** Private member functions ******** */
      /** 
       * Allocate memory for array data and initialize reference count.
       */
      void
      allocate(size_t arrayRows, size_t arrayColumns);
    

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
      Type* m_dataPtr;
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
  
  } // namespace numeric

} // namespace brick


/*******************************************************************
 * Member function definitions follow.  This would be a .cc file
 * if it weren't templated.
 *******************************************************************/

#include <algorithm>
#include <functional>
#include <sstream>
#include <vector>
#include <brick/common/expect.hh>
#include <brick/common/functional.hh>
#include <brick/numeric/functional.hh>
#include <brick/numeric/numericTraits.hh>

namespace brick {

  namespace numeric {
    
    // Static constant describing how the string representation of an
    // Array2D should start.
    template <class Type>
    const std::string&
    Array2D<Type>::
    ioIntro()
    {
      static const std::string intro = "Array2D(";
      return intro;
    }

    // Static constant describing how the string representation of an
    // Array2D should end.
    template <class Type>
    const char&
    Array2D<Type>::
    ioOutro()
    {
      static const char outro = ')';
      return outro;
    }


    // Static constant describing how the the data portion of the
    // string representation of an Array1D should start.
    template <class Type>
    const char&
    Array2D<Type>::
    ioOpening()
    {
      static const char opening = '[';
      return opening;
    }


    // Static constant describing how the the data portion of the
    // string representation of an Array2D should end.
    template <class Type>
    const char&
    Array2D<Type>::
    ioClosing()
    {
      static const char closing = ']';
      return closing;
    }


    // Static constant describing how individual elements should be
    // separated in the string representation of Array2D.
    template <class Type>
    const char&
    Array2D<Type>::
    ioSeparator()
    {
      static const char separator = ',';
      return separator;
    }

    // Non-static member functions below.

    template <class Type>
    Array2D<Type>::
    Array2D()
      : m_rows(0),
        m_columns(0),
        m_size(0),
        m_dataPtr(0),
        m_referenceCount(0)
    {
      // Empty
    }


    template <class Type>
    Array2D<Type>::
    Array2D(size_t arrayRows, size_t arrayColumns)
      : m_rows(arrayRows),
        m_columns(arrayColumns),
        m_size(0),           // This will be set in the call to allocate().
        m_dataPtr(0),        // This will be set in the call to allocate().
        m_referenceCount(0)  // This will be set in the call to allocate().
    {
      this->allocate(arrayRows, arrayColumns);
    }

  
    // Construct from an initialization string.
    template <class Type>
    Array2D<Type>::
    Array2D(const std::string& inputString)
      : m_rows(0),
        m_columns(0),
        m_size(0),
        m_dataPtr(0),
        m_referenceCount(0)
    {
      // We'll use the stream input operator to parse the string.
      std::istringstream inputStream(inputString);

      // Now read the string into an array.
      Array2D<Type> inputArray;
      inputStream >> inputArray;
      if(!inputStream) {
        std::ostringstream message;
        message << "Couldn't parse input string: \"" << inputString << "\".";
        BRICK_THROW(common::ValueException, "Array2D::Array2D(const std::string&)",
                   message.str().c_str());                 
      }

      // If all went well, copy into *this.
      *this = inputArray;
    }


    /* When copying from a Array2D do a shallow copy */
    /* Update reference count if the array we're copying has */
    /* valid data. */
    template <class Type>
    Array2D<Type>::
    Array2D(const Array2D<Type>& source)
      : m_rows(source.m_rows),
        m_columns(source.m_columns),
        m_size(source.m_size),
        m_dataPtr(source.m_dataPtr),
        m_referenceCount(source.m_referenceCount)
    {
      // Empty.
    }

  
    /* Here's a constructor for getting image data into the array */
    /* cheaply. */
    template <class Type>
    Array2D<Type>::
    Array2D(size_t arrayRows, size_t arrayColumns, Type* const dataPtr)
      : m_rows(arrayRows),
        m_columns(arrayColumns),
        m_size(arrayRows * arrayColumns),
        m_dataPtr(dataPtr),
        m_referenceCount(0)
    {
      // Empty
    }

  
    // Construct an array around external data that was allocated by
    // an Array?D instance.
    template <class Type>
    Array2D<Type>::
    Array2D(size_t arrayRows, size_t arrayColumns, Type* const dataPtr,
            common::ReferenceCount const& referenceCount)
      : m_rows(arrayRows),
        m_columns(arrayColumns),
        m_size(arrayRows*arrayColumns),
        m_dataPtr(dataPtr),
        m_referenceCount(referenceCount)
    {
      // Empty.
    }

  
    template <class Type>
    Array2D<Type>::
    ~Array2D()
    {
      deAllocate();
    }

  
    template <class Type>
    inline void Array2D<Type>::
    checkDimension(size_t
#ifdef BRICK_NUMERIC_CHECKBOUNDS
                   arrayRows
#endif /* #ifdef BRICK_NUMERIC_CHECKBOUNDS */
                   , size_t
#ifdef BRICK_NUMERIC_CHECKBOUNDS
                   arrayColumns
#endif /* #ifdef BRICK_NUMERIC_CHECKBOUNDS */
      ) const
    {
#ifdef BRICK_NUMERIC_CHECKBOUNDS
      if(arrayRows != this->rows()
         || arrayColumns != this->columns()) {
        std::ostringstream message;
        message << "Size mismatch: required dimension is ("
                << arrayRows << ", " << arrayColumns << ") "
                << " while *this has dimension "
                << this->rows() << ", " << this->columns() << ").";
        BRICK_THROW(common::IndexException, "Array2D::checkDimension()",
                    message.str().c_str());
      }
#endif /* #ifdef BRICK_NUMERIC_CHECKBOUNDS */
    }


    template <class Type>
    Array2D<Type> Array2D<Type>::
    copy() const
    {
      Array2D<Type> newArray(m_rows, m_columns);
      newArray.copy(*this);
      return newArray;
    }

    
    template <class Type> template <class Type2>
    void Array2D<Type>::
    copy(const Array2D<Type2>& source)
    {
      if(source.size() != m_size) {
        std::ostringstream message;
        message << "Mismatched array sizes. Source array has "
                << source.size() << " elements, while destination array has "
                << m_size << " elements.";
        BRICK_THROW(common::ValueException, "Array2D::copy(const Array2D&)",
                    message.str().c_str());
      }
      if(m_size != 0) {
        this->copy(source.data());
      }
    }


    template <class Type> template <class Type2>
    void Array2D<Type>::
    copy(const Type2* dataPtr)
    {
      if(dataPtr == 0) {
        BRICK_THROW(common::ValueException, "Array2D::copy(const Type2*)",
                    "Argument is a NULL pointer.");
      }
      std::transform(dataPtr, dataPtr + m_size, m_dataPtr,
                     StaticCastFunctor<Type2, Type>());
    }


    template <class Type>
    Array1D<Type> Array2D<Type>::
    getRow(size_t index)
    {
      this->checkBounds(index, 0);
      return Array1D<Type>(this->columns(),
                           m_dataPtr + (index * this->m_columns));
    }

  
    template <class Type>
    const Array1D<Type> Array2D<Type>::
    getRow(size_t index) const
    {
      this->checkBounds(index, 0);
      return Array1D<Type>(this->columns(),
                           m_dataPtr + (index * this->m_columns));
    }


    template <class Type>
    const Array1D<Type> Array2D<Type>::
    ravel() const
    {
      if(this->isReferenceCounted()) {
        return Array1D<Type>(m_size, m_dataPtr, m_referenceCount);
      }
      return Array1D<Type>(m_size, m_dataPtr);
    }

    
    template <class Type>
    Array1D<Type> Array2D<Type>::
    ravel()
    {
      if(this->isReferenceCounted()) {
        return Array1D<Type>(m_size, m_dataPtr, m_referenceCount);
      }
      return Array1D<Type>(m_size, m_dataPtr);
    }


    template <class Type>
    std::istream&
    Array2D<Type>::
    readFromStream(std::istream& inputStream)
    {
      // Most of the time, InputType will be the same as Type.
      // TextOutputType is the type one would use to write out a value
      // of Type.  Not surprisingly, this is the type you need to use
      // when reading those values back in.
      typedef typename NumericTraits<Type>::TextOutputType InputType;

      // If stream is in a bad state, we can't read from it.
      if (!inputStream){
        return inputStream;
      }
    
      // It's a lot easier to use a try block than to be constantly
      // testing whether the IO has succeeded, so we tell inputStream to
      // complain if anything goes wrong.
      std::ios_base::iostate oldExceptionState = inputStream.exceptions();
      inputStream.exceptions(
        std::ios_base::badbit | std::ios_base::failbit | std::ios_base::eofbit);

      // Now on with the show.
      try{
        common::Expect::FormatFlag flags = common::Expect::SkipWhitespace;

        // Skip any preceding whitespace.
        inputStream >> common::Expect("", flags);
              
        // We won't require the input format to start with "Array2D(", but
        // if it does we read it here.
        bool foundIntro = false;
        if(inputStream.peek() == ioIntro()[0]) {
          foundIntro = true;
          inputStream >> common::Expect(ioIntro(), flags);
        }

        // OK.  We've dispensed with the intro.  What's left should be of
        // the format "[row, row, row, ...]".  We require the square
        // brackets to be there.
        inputStream >> common::Expect(&(ioOpening()), 1, flags);

        // Read the data.  We'll use the Array1D<Type> stream operator to
        // read each row.
        Array1D<Type> inputValue;
        std::vector< Array1D<Type> > inputBuffer;
        while(1) {
          // Read the next row.
          inputStream >> inputValue;
          inputBuffer.push_back(inputValue);

          // Read the separator, or else the closing character.
          char inChar = 0;
          inputStream >> inChar;
          if(inChar == ioClosing()) {
            // Found a closing.  Stop here.
            break;
          }
          if(inChar != ioSeparator()) {
            // Missing separator?  Fail here.
            inputStream.clear(std::ios_base::failbit);
          }
        }
    
        // If we found an intro, we expect the corresponding outro.
        if(foundIntro) {
          inputStream >> common::Expect(&(ioOutro()), 1, flags);
        }

        // Now we're done with all of the parsing, verify that all rows
        // have the same length.
        size_t arrayRows = inputBuffer.size();
        size_t arrayColumns = ((inputBuffer.size() != 0)
                               ? inputBuffer[0].size() : 0);
        for(size_t index = 1; index < arrayRows; ++index) {
          if(inputBuffer[index].size() != arrayColumns) {
            // Inconsistent row lengths!  Fail here.
            inputStream.clear(std::ios_base::failbit);
          }
        }

        // And finally, copy the data.
        this->reinit(arrayRows, arrayColumns);
        for(size_t index = 0; index < arrayRows; ++index) {
          std::copy(inputBuffer[index].begin(), inputBuffer[index].end(),
                    this->begin() + (index * arrayColumns));
        }
      } catch(std::ios_base::failure) {
        // Empty
      }
      inputStream.exceptions(oldExceptionState);
      return inputStream;
    }
  

    template <class Type>
    void Array2D<Type>::
    reinit(size_t arrayRows, size_t arrayColumns)
    {
      this->allocate(arrayRows, arrayColumns);
    }


    template <class Type>
    void Array2D<Type>::
    reinitIfNecessary(size_t arrayRows, size_t arrayColumns)
    {
      if(this->size() != arrayRows * arrayColumns) {
        this->reinit(arrayRows, arrayColumns);
      } else {
        if(this->rows() != arrayRows) {
          this->reshape(arrayRows, arrayColumns);
        }
      }
    }


    /* After reshaping, matrix is still row major order */
    template <class Type>
    void Array2D<Type>::
    reshape(int arrayRows, int arrayColumns)
    {
      // If one axis is specified as -1, it will be automatically 
      // chosen to match the number of elements in the array.
      if((arrayRows == -1) && (arrayColumns != 0)) {
        arrayRows = static_cast<int>(this->size()) / arrayColumns;
      } else
        if((arrayColumns == -1) && (arrayRows != 0)) {
          arrayColumns = static_cast<int>(this->size()) / arrayRows;
        }
      if((arrayRows * arrayColumns) != static_cast<int>(this->size())) {
        std::ostringstream message;
        message << "Can't reshape a(n) " << this->size()
                << " element array to have " << arrayRows << " rows and "
                << arrayColumns << " columns.";
        BRICK_THROW(common::ValueException, "Array2D::reshape()",
                    message.str().c_str());
      }
      m_rows = arrayRows;
      m_columns = arrayColumns;
    }

  
    template <class Type>
    Array1D<size_t> Array2D<Type>::
    shape() const
    {
      Array1D<size_t> rc(2);
      rc(0) = this->rows();
      rc(1) = this->columns();
      return rc;
    }


    template <class Type>
    size_t Array2D<Type>::
    shape(size_t axis) const
    {
      size_t result;
      switch(axis) {
      case 0:
        result = this->rows();
        break;
      case 1:
        result = this->columns();
        break;
      default:
        std::ostringstream message;
        message << "Invalid Axis: "<< axis << ".";
        BRICK_THROW(common::ValueException, "Array2D::shape(size_t)",
                    message.str().c_str());
        result = 0;
        break;
      }
      return result;
    }


    template <class Type>
    Array2D<Type>& Array2D<Type>::
    operator=(Type value)
    {
      std::fill(m_dataPtr, m_dataPtr + m_size, value);
      return *this;
    }

  
    template <class Type>
    Array2D<Type>& Array2D<Type>::
    operator=(const Array2D<Type>& source)
    {
      // Check for self-assignment
      if(&source != this) {
        this->deAllocate();
        m_rows = source.m_rows;
        m_columns = source.m_columns;
        m_size = source.m_size;
        m_dataPtr = source.m_dataPtr;
        m_referenceCount = source.m_referenceCount;
      }
      return *this;
    }


    template <class Type> template <class Type2>
    Array2D<Type>&
    Array2D<Type>::
    operator+=(const Array2D<Type2>& arg)
    {
      if(m_size != arg.size()) {
        std::ostringstream message;
        message << "Mismatched array sizes. Argument array is "
                << arg.rows() << " x " << arg.columns()
                << ", while destination array is "
                << m_rows << " x " << m_columns << ".";
        BRICK_THROW(common::ValueException, "Array2D::operator+=()",
                    message.str().c_str());
      }
      std::transform(m_dataPtr, m_dataPtr + m_size, arg.data(), m_dataPtr,
                     std::plus<Type>());
      return *this;
    }


    template <class Type> template <class Type2>
    Array2D<Type>&
    Array2D<Type>::
    operator-=(const Array2D<Type2>& arg)
    {
      if(m_size != arg.size()) {
        std::ostringstream message;
        message << "Mismatched array sizes. Argument array is "
                << arg.rows() << " x " << arg.columns()
                << ", while destination array is "
                << m_rows << " x " << m_columns << ".";
        BRICK_THROW(common::ValueException, "Array2D::operator-=()",
                    message.str().c_str());
      }
      std::transform(m_dataPtr, m_dataPtr + m_size, arg.data(), m_dataPtr,
                     std::minus<Type>());
      return *this;
    }


    template <class Type> template <class Type2>
    Array2D<Type>&
    Array2D<Type>::
    operator*=(const Array2D<Type2>& arg)
    {
      if(m_size != arg.size()) {
        std::ostringstream message;
        message << "Mismatched array sizes. Argument array is "
                << arg.rows() << " x " << arg.columns()
                << ", while destination array is "
                << m_rows << " x " << m_columns << ".";
        BRICK_THROW(common::ValueException, "Array2D::operator*=()",
                    message.str().c_str());
      }
      std::transform(m_dataPtr, m_dataPtr + m_size, arg.data(), m_dataPtr,
                     std::multiplies<Type>());
      return *this;
    }


    template <class Type> template <class Type2>
    Array2D<Type>&
    Array2D<Type>::
    operator/=(const Array2D<Type2>& arg)
    {
      if(m_size != arg.size()) {
        std::ostringstream message;
        message << "Mismatched array sizes. Argument array is "
                << arg.rows() << " x " << arg.columns()
                << ", while destination array is "
                << m_rows << " x " << m_columns << ".";
        BRICK_THROW(common::ValueException, "Array2D::operator/=()",
                    message.str().c_str());
      }
      std::transform(m_dataPtr, m_dataPtr + m_size, arg.data(), m_dataPtr,
                     std::divides<Type>());
      return *this;
    }


    template <class Type>
    Array2D<Type>&
    Array2D<Type>::
    operator*=(Type arg)
    {
      std::transform(m_dataPtr, m_dataPtr + m_size, m_dataPtr,
                     std::bind2nd(std::multiplies<Type>(), arg));
      return *this;
    }


    template <class Type>
    Array2D<Type>&
    Array2D<Type>::
    operator/=(Type arg)
    {
      std::transform(m_dataPtr, m_dataPtr + m_size, m_dataPtr,
                     std::bind2nd(std::divides<Type>(), arg));
      return *this;
    }


    template <class Type>
    Array2D<Type>&
    Array2D<Type>::
    operator+=(Type arg)
    {
      std::transform(m_dataPtr, m_dataPtr + m_size, m_dataPtr,
                     std::bind2nd(std::plus<Type>(), arg));
      return *this;
    }


    template <class Type>
    Array2D<Type>&
    Array2D<Type>::
    operator-=(Type arg)
    {
      std::transform(m_dataPtr, m_dataPtr + m_size, m_dataPtr,
                     std::bind2nd(std::minus<Type>(), arg));
      return *this;
    }


    template <class Type>
    Array2D<Type> Array2D<Type>::
    transpose() const
    {
      Array2D<Type> newMx(m_columns, m_rows);

      // Waiting for row & column iterators
      Type *tPtr0 = newMx.m_dataPtr;
      for(size_t j = 0; j < m_columns; ++j) {
        // const Type *tPtr1 = this->data(0, j);
        const Type *tPtr1 = this->data(j);
        for(size_t i = 0; i < m_rows; ++i) {
          *tPtr0 = *tPtr1;
          ++tPtr0;
          tPtr1 += m_columns;
        }
      }
      return newMx;
    }


    template <class Type>
    void Array2D<Type>::
    allocate(size_t arrayRows, size_t arrayColumns)
    {
      // Make sure to release any shared resources.
      this->deAllocate();
      
      // Check array size.  It doesn't make sense to allocate memory
      // for a zero size array.
      m_rows = arrayRows;
      m_columns = arrayColumns;
      m_size = m_rows * m_columns;
      if(m_size > 0) {
        // Allocate data storage.  new() should throw an exception if
        // we run out of memory.
        m_dataPtr = new(Type[m_size]);

        // Set reference count to show that exactly one Array is pointing
        // to this data.
        m_referenceCount.reset(1);
      }
      return;
    }


    template <class Type>
    inline void Array2D<Type>::
    checkBounds(size_t
#ifdef BRICK_NUMERIC_CHECKBOUNDS
                index
#endif /* #ifdef BRICK_NUMERIC_CHECKBOUNDS */
      ) const
    {
#ifdef BRICK_NUMERIC_CHECKBOUNDS
      if(index >= m_size) {
        std::ostringstream message;
        message << "Index " << index << " is invalid for a(n) " << m_rows
                << " x " << m_columns << " array.";
        BRICK_THROW(common::IndexException, "Array2D::checkBounds(size_t)",
                    message.str().c_str());
      }
#endif /* #ifdef BRICK_NUMERIC_CHECKBOUNDS */
    }


    template <class Type>
    inline void Array2D<Type>::
    checkBounds(size_t
#ifdef BRICK_NUMERIC_CHECKBOUNDS
                rowIndex
#endif /* #ifdef BRICK_NUMERIC_CHECKBOUNDS */
                , size_t
#ifdef BRICK_NUMERIC_CHECKBOUNDS
                columnIndex
#endif /* #ifdef BRICK_NUMERIC_CHECKBOUNDS */
      ) const
    {
#ifdef BRICK_NUMERIC_CHECKBOUNDS
      if(rowIndex >= m_rows) {
        std::ostringstream message;
        message << "Row index " << rowIndex << " is invalid for a(n) "
                << m_rows << " x " << m_columns << " array.";
        BRICK_THROW(common::IndexException,
                    "Array2D::checkBounds(size_t, size_t)",
                    message.str().c_str());
      }
      if(columnIndex >= m_columns) {
        std::ostringstream message;
        message << "Column index " << columnIndex << " is invalid for a(n) "
                << m_rows << " x " << m_columns << " array.";
        BRICK_THROW(common::IndexException,
                    "Array2D::checkBounds(size_t, size_t)",
                    message.str().c_str());
      }
#endif
    }


    template <class Type>
    void Array2D<Type>::
    deAllocate()
    {
      // Are we responsible for deallocating the contents of this array?
      if(m_referenceCount.isCounted()) {
        // If yes, are we currently the only array pointing to this data?
        if(!m_referenceCount.isShared()) {
          // If yes, then delete the data.
          delete[] m_dataPtr;
        }
      }
      // Abandon our pointers to data.  Reference would take care of
      // itself, but it's cleaner conceptually to wipe it here, rather
      // than waiting for a subsequent call to ~ReferenceCount() or
      // ReferenceCount.allocate().
      m_dataPtr = 0;
      m_size = 0;
      m_rows = 0;
      m_columns = 0;
      m_referenceCount.reset(0);
    }

    /* Non-member functions that will ultimately wind up in a different file */

    template <class Type>
    inline Array2D<Type>
    sqrt(const Array2D<Type>& array0)
    {
      return squareRoot(array0);
    }

    // This function returns an Array2D instance of the same shape and
    // element type as its input, in which each element contains the
    // square root of the corresponding element of the input array.
    template <class Type>
    Array2D<Type>
    squareRoot(const Array2D<Type>& array0)
    {
      Array2D<Type> result(array0.rows(), array0.columns());
      std::transform(array0.begin(), array0.end(),
                     result.begin(), SquareRootFunctor<Type>());
      return result;
    }

    template <class Type>
    Array2D<Type> operator+(const Array2D<Type>& array0,
                            const Array2D<Type>& array1)
    {
      if((array0.rows() != array1.rows())
         || (array0.columns() != array1.columns())) {
        std::ostringstream message;
        message << "Array sizes do not match.  Array0 is " << array0.rows()
                << " x " << array0.columns() << ", while array1 is "
                << array1.rows() << " x " << array1.columns() << ".";
        BRICK_THROW(common::ValueException, "Array2D::operator+()",
                    message.str().c_str());
      }
      Array2D<Type> result(array0.rows(), array0.columns());
      std::transform(array0.begin(), array0.end(), array1.begin(),
                     result.begin(), std::plus<Type>());
      return result;
    }

    template <class Type>
    Array2D<Type> operator-(const Array2D<Type>& array0,
                            const Array2D<Type>& array1)
    {
      if((array0.rows() != array1.rows())
         || (array0.columns() != array1.columns())) {
        std::ostringstream message;
        message << "Array sizes do not match.  Array0 is " << array0.rows()
                << " x " << array0.columns() << ", while array1 is "
                << array1.rows() << " x " << array1.columns() << ".";
        BRICK_THROW(common::ValueException, "Array2D::operator-()", message.str().c_str());
      }
      Array2D<Type> result(array0.rows(), array0.columns());
      std::transform(array0.begin(), array0.end(), array1.begin(),
                     result.begin(), std::minus<Type>());
      return result;
    }

    template <class Type>
    Array2D<Type> operator*(const Array2D<Type>& array0,
                            const Array2D<Type>& array1)
    {
      if((array0.rows() != array1.rows())
         || (array0.columns() != array1.columns())) {
        std::ostringstream message;
        message << "Array sizes do not match.  Array0 is " << array0.rows()
                << " x " << array0.columns() << ", while array1 is "
                << array1.rows() << " x " << array1.columns() << ".";
        BRICK_THROW(common::ValueException, "Array2D::operator*()", message.str().c_str());
      }
      Array2D<Type> result(array0.rows(), array0.columns());
      std::transform(array0.begin(), array0.end(), array1.begin(),
                     result.begin(), std::multiplies<Type>());
      return result;
    }

    template <class Type>
    Array2D<Type> operator/(const Array2D<Type>& array0,
                            const Array2D<Type>& array1)
    {
      if((array0.rows() != array1.rows())
         || (array0.columns() != array1.columns())) {
        std::ostringstream message;
        message << "Array sizes do not match.  Array0 is " << array0.rows()
                << " x " << array0.columns() << ", while array1 is "
                << array1.rows() << " x " << array1.columns() << ".";
        BRICK_THROW(common::ValueException, "Array2D::operator/()", message.str().c_str());
      }
      Array2D<Type> result(array0.rows(), array0.columns());
      std::transform(array0.begin(), array0.end(), array1.begin(),
                     result.begin(), std::divides<Type>());
      return result;
    }

    template <class Type>
    Array2D<Type> operator+(const Array2D<Type>& array0, Type scalar)
    {
      Array2D<Type> result(array0.rows(), array0.columns());
      std::transform(array0.begin(), array0.end(), result.begin(),
                     std::bind2nd(std::plus<Type>(), scalar));
      return result;
    }

    template <class Type>
    Array2D<Type> operator-(const Array2D<Type>& array0, Type scalar)
    {
      Array2D<Type> result(array0.rows(), array0.columns());
      std::transform(array0.begin(), array0.end(), result.begin(),
                     std::bind2nd(std::minus<Type>(), scalar));
      return result;
    }

    template <class Type>
    Array2D<Type> operator*(const Array2D<Type>& array0, Type scalar)
    {
      Array2D<Type> result(array0.rows(), array0.columns());
      std::transform(array0.begin(), array0.end(), result.begin(),
                     std::bind2nd(std::multiplies<Type>(), scalar));
      return result;
    }

    template <class Type>
    Array2D<Type> operator/(const Array2D<Type>& array0, Type scalar)
    {
      Array2D<Type> result(array0.rows(), array0.columns());
      std::transform(array0.begin(), array0.end(), result.begin(),
                     std::bind2nd(std::divides<Type>(), scalar));
      return result;
    }

    template <class Type>
    inline Array2D<Type> operator+(Type scalar, const Array2D<Type>& array0)
    {
      return array0 + scalar;
    }

    template <class Type>
    inline Array2D<Type> operator*(Type scalar, const Array2D<Type>& array0)
    {
      return array0 * scalar;
    }


    // Elementwise comparison of an Array2D with a constant.
    template <class Type>
    Array2D<bool>
    operator==(const Array2D<Type>& array0, const Type arg)
    {
      Array2D<bool> result(array0.rows(), array0.columns());
      std::transform(array0.begin(), array0.end(), result.data(),
                     std::bind2nd(std::equal_to<Type>(), arg));
      return result;
    }

    
    // Elementwise comparison of an Array2D with another array.
    template <class Type>
    Array2D<bool>
    operator==(const Array2D<Type>& array0, const Array2D<Type>& array1)
    {
      array0.checkDimension(array1.rows(), array1.columns());
      Array2D<bool> result(array0.rows(), array0.columns());
      std::transform(array0.begin(), array0.end(), array1.begin(),
                     result.begin(), std::equal_to<Type>());
      return result;
    }

  
    template <class Type>
    Array2D<bool> operator>(const Array2D<Type>& array0, Type arg)
    {
      Array2D<bool> result(array0.rows(), array0.columns());
      std::transform(array0.begin(), array0.end(), result.begin(),
                     std::bind2nd(std::greater<Type>(), arg));
      return result;
    }

    template <class Type>
    Array2D<bool> operator<(const Array2D<Type>& array0, Type arg)
    {
      Array2D<bool> result(array0.rows(), array0.columns());
      std::transform(array0.begin(), array0.end(), result.begin(),
                     std::bind2nd(std::less<Type>(), arg));
      return result;
    }

    template <class Type>
    Array2D<bool> operator>=(const Array2D<Type>& array0, Type arg)
    {
      Array2D<bool> result(array0.rows(), array0.columns());
      std::transform(array0.begin(), array0.end(), result.begin(),
                     std::bind2nd(std::greater_equal<Type>(), arg));
      return result;
    }

    template <class Type>
    Array2D<bool> operator<=(const Array2D<Type>& array0, Type arg)
    {
      Array2D<bool> result(array0.rows(), array0.columns());
      std::transform(array0.begin(), array0.end(), result.begin(),
                     std::bind2nd(std::less_equal<Type>(), arg));
      return result;
    }

    template <class Type>
    std::ostream& operator<<(std::ostream& stream, const Array2D<Type>& array0)
    {
      // Most of the time, OutputType will be the same as Type.
      typedef typename NumericTraits<Type>::TextOutputType OutputType;

      stream << "Array2D([[";
      for(size_t row = 0; row < array0.rows(); ++row) {
        if (array0.columns() > 0) {
          for(size_t column = 0; column < array0.columns() - 1; ++column) {
            stream << static_cast<OutputType>(array0(row, column)) << ", ";
          }
          stream << static_cast<OutputType>(array0(row, array0.columns() - 1));
          if(row != array0.rows() - 1) {
            stream << "],\n";
            stream << "         [";
          }
        }
      }	
      stream << "]])";
      stream.flush();
      return stream;
    }

  
    template <class Type>
    std::istream& operator>>(std::istream& inputStream, Array2D<Type>& array0)
    {
      return array0.readFromStream(inputStream);
    }
  
  } // namespace numeric

} // namespace brick

#endif /* #ifdef BRICK_ARRAY2D_HH */
