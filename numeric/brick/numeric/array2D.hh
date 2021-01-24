/**
***************************************************************************
* @file brick/numeric/array2D.hh
*
* Header file declaring Array2D class.
*
* Copyright (C) 2001-2011 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_ARRAY2D_HH
#define BRICK_NUMERIC_ARRAY2D_HH

// Early inclusion of <algorithm> to allow inlining of member templates.
#include <algorithm>

#include <iostream>
#include <brick/common/exception.hh>
#include <brick/numeric/array1D.hh>
#include <brick/numeric/index2D.hh>

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
       * WARNING: The rowStep argument is newly added, and much of the
       * brick source code assumes that rows are contiguous.  If you
       * set it to anything other than the default, be on the lookout
       * for bugs.
       *
       * @param arrayRows Number of rows in the array after successful
       * construction.
       *
       * @param arrayColumns Number of columns in the array after successful
       * construction.
       *
       * @param rowStep Number of elements between the start of
       * consecutive rows in the array.  Setting this argument to zero
       * has the same effect as setting it ot the same value as
       * arrayColumns.  This argument is useful for creating arrays
       * with rows aligned on, for example, 16 byte boundaries.
       */
      Array2D(size_t arrayRows, size_t arrayColumns,
              size_t rowStep = 0);


      /**
       * Construct from an initialization string.  The idea is to make
       * this constructor very flexible about interpreting string syntax.
       * For now, though, the input string must have the format:
       *
       * "[[#, #, #, ...], [#, #, #, ...], [#, #, #, ...], ...]"
       *
       * Where "#" indicates text that can be converted to the element
       * type of the array using the stream input operator.
       *
       * WARNING: The rowStep argument is newly added, and much of the
       * brick source code assumes that rows are contiguous.  If you
       * set it to anything other than the default, be on the lookout
       * for bugs.
       *
       * @param inputString The argument specifies the string from which
       * the array will be constructed.
       *
       * @param rowStep Number of elements between the start of
       * consecutive rows in the array.  Setting this argument to zero
       * has the same effect as setting it ot the same value as
       * arrayColumns.  This argument is useful for creating arrays
       * with rows aligned on, for example, 16 byte boundaries.
       */
      explicit
      Array2D(const std::string& inputString,
              size_t rowStep = 0);


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
       * WARNING: The rowStep argument is newly added, and much of the
       * brick source code assumes that rows are contiguous.  If you
       * set it to anything other than the default, be on the lookout
       * for bugs.
       *
       * @param arrayRows Number of rows in the array after successful
       * construction.
       *
       * @param arrayColumns Number of columns in the array after successful
       * construction.
       *
       * @param dataPtr A C-style array of Type into which the newly
       * constructed Array2D should index.
       *
       * @param rowStep Number of elements between the start of
       * consecutive rows in the array.  Setting this argument to zero
       * has the same effect as setting it ot the same value as
       * arrayColumns.  This argument is useful for creating arrays
       * with rows aligned on, for example, 16 byte boundaries.
       */
      Array2D(size_t arrayRows, size_t arrayColumns, Type* const dataPtr,
              size_t rowStep = 0);


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
       * Construct an array around external data that was allocated by
       * an Array?D instance.  Arrays constructed in this way _do_
       * implement reference counting, and will delete dataPtr when
       * done.  This constructor is provided primarily so that other
       * dimensionality array classes can return Array2D instances that
       * reference their data without being friend classes.  Caveat
       * emptor.
       *
       * WARNING: The rowStep argument is newly added, and much of the
       * brick source code assumes that rows are contiguous.  If you
       * set it to anything other than the default, be on the lookout
       * for bugs.
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
       * @param rowStep Number of elements between the start of
       * consecutive rows in the array.  Setting this argument to zero
       * has the same effect as setting it to the same value as
       * arrayColumns.  This argument is useful for creating arrays
       * with rows aligned on, for example, 16 byte boundaries.
       *
       * @param referenceCount ReferenceCount instance indicating
       * the number of Array classes currently using dataPtr.
       */
      Array2D(size_t arrayRows, size_t arrayColumns, Type* const dataPtr,
              size_t rowStep, common::ReferenceCount const& referenceCount);


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
      Array2D(std::initializer_list< std::initializer_list<Type> > initializer);


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
       * Just like getData(void), which is documented above, but
       * returns a pointer to the element indexed by (row, column).
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
        return m_dataPtr + columnIndex + (rowIndex * m_rowStep);
      }

      /**
       * This version of getData(size_t, size_t) is appropriate for
       * const Array2D, and returns a pointer-to-const.  @param
       * rowIndex The row of the element to which the return value
       * should point.
       *
       * @param columnIndex The column of the element to which the return
       * value should point.
       *
       * @return Const pointer to the selected element of the array.
       */
      const Type*
      getData(size_t rowIndex, size_t columnIndex) const {
        this->checkBounds(rowIndex, columnIndex);
        return m_dataPtr + columnIndex + (rowIndex * m_rowStep);
      }


      /**
       * Just like getData(void), which is documented above, but
       * returns a pointer to the element indexed by (idx.getRow(),
       * idx.getColumn()).
       *
       * @param idx The coordinates of the element to which the return
       * value should point.
       *
       * @return Pointer to the selected element of the array.
       */
      Type*
      getData(Index2D const& idx) {
        return this->getData(idx.getRow(), idx.getColumn());
      }


      /**
       * This version of getData(Index2D const&) is appropriate for
       * const Array2D, and returns a pointer-to-const.
       *
       * @param idx The coordinates of the element to which the return
       * value should point.
       *
       * @return Const pointer to the selected element of the array.
       */
      const Type*
      getData(Index2D const& idx) const {
        return this->getData(idx.getRow(), idx.getColumn());
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
       * Return an Array2D instance referencing a subset of *this.
       * The returned array will reference the same memory, but may
       * have different start, end, and rowstep.  The returned array
       * _will not_ be reference counted, so that its contents are
       * valid only as long as the contents of the original array are
       * valid.  The operation of this function is illustrated in the
       * following example:
       *
       * @code
       *   // Create an array in which every element is set to 1.
       *   Array2D<int>* fullArrayP = new Array2D<int>(10, 20);
       *   *fullArrayP = 1;
       *
       *   // Get an array that refers to a 4x4 subset of fullArray.
       *   Index2D corner0(5, 10);
       *   Index2D corner1(9, 14);
       *   Array2D<int> subRegion = fullArrayP->getRegion(corner0, corner1);
       *
       *   // Set just the elements within that subset to 100.
       *   subRegion = 100;
       *
       *   for(unsigned int rr = 0, rr < fullArray.rows(); ++rr) {
       *     for(unsigned int cc = 0, cc < fullArray.columns(); ++cc) {
       *       if((rr >= corner0.getRow()) &&
       *          (rr <  corner1.getRow()) &&
       *          (cc >= corner0.getColumn()) &&
       *          (cc <  corner1.getColumn())) {
       *         std::assert(100 == (*fullArrayP)(rr, cc));
       *       } else {
       *         std::assert(1   == (*fullArrayP)(rr, cc));
       *       }
       *     }
       *   }
       *
       *   // Contents of *fullArrayP are deleted here.
       *   delete fullArrayP;
       *
       *   // ERROR! The data to which subRegion refers was deleted
       *   // along with fullArrayP.
       *   subRegion(2, 3) = 90;
       * @endCode
       *
       * @param corner0 This argument and the next define the
       * subregion.  Corner0 is included in the region, but corner1 is
       * not.  It is an error if corner0 is not above and to the left
       * of corner1.  It is an error if corner0.getRow() is less than
       * zero, corner0.getColumn() is less than zero, corner0.getRow()
       * is greater than this->getRows(), or corner0.getColumn() is
       * greater than this->getColumns().
       *
       * @param corner1 This argument and the previous define the
       * subregion.  It is an error if corner1.getRow() is less than
       * zero, corner1.getColumn() is less than zero, corner1.getRow()
       * is greater than this->getRows(), or corner1.getColumn() is
       * greater than this->getColumns().
       *
       * @return The return value is a shallow copy of the selected
       * region of the original array.
       */
      Array2D<Type>
      getRegion(Index2D const& corner0, Index2D const& corner1);


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
       * Returns the spacing (in elements, not bytes) between the
       * beginning of subsequent rows of the array.
       *
       * @return Number elements between consecutive rows.
       */
      size_t
      getRowStep() const {return m_rowStep;}


      /**
       * Returns the number of elements in the array.  This is the
       * product of rows() and columns().
       *
       * @return Number of elements.
       */
      size_t
      getSize() const {return m_size;}


      /**
       * Returns the number of elements in the storage area used by
       * the array.  This is the product of rows() and getRowStep().
       *
       * @return Number of elements.
       */
      size_t
      getStorageSize() const {return m_storageSize;}


      /**
       * Indicates whether or not there is unused space between rows
       * of the array.
       *
       * @return The return value is true if rows are contiguous,
       * false if there is padding.
       */
      bool
      isContiguous() const {return m_columns == m_rowStep;}


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
      reinit(size_t arrayRows, size_t arrayColumns,
             size_t rowStep = 0);


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
      reinitIfNecessary(size_t arrayRows, size_t arrayColumns,
                        size_t rowStep = 0);


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
      reshape(int arrayRows, int arrayColumns, int rowStep = 0);


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
      rowBegin(size_t rowIndex) {return m_dataPtr + rowIndex * m_rowStep;}


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
	return m_dataPtr + rowIndex * m_rowStep;
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
	return m_dataPtr + (rowIndex + 1) * m_rowStep;
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
	return m_dataPtr + (rowIndex + 1) * m_rowStep;
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
      Type const& operator()(size_t index) const {
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
      Type const&
      operator()(size_t rowIndex, size_t columnIndex) const {
        this->checkBounds(rowIndex, columnIndex);
        return m_dataPtr[columnIndex + rowIndex * m_rowStep];
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
      Type const&
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
      /* ******** Private member functions ******** */
      /**
       * Allocate memory for array data and initialize reference count.
       */
      void
      allocate(size_t arrayRows, size_t arrayColumns, size_t rowStep = 0);


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


#if 0
      /**
       * Computes the size, in number of elements, of the storage
       * required for a (rows x columns) array in which rows are
       * aligned on rowStep element boundaries.
       *
       * @param rows Number of rows in the target array.
       *
       * @param columns Number of columns in the target array.
       *
       * @param rowStep Number of elements between subsequent rows
       * (including end-of-row padding).  Zero indicates no padding
       * (same as rowStep == columns).
       *
       * @return The return value is the size (number of Type
       * instances, not number of bytes) of the required storage.
       */
      size_t
      computeStorageSize(size_t rows, size_t columns, size_t rowStep);
#endif


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
      size_t m_rowStep;
      size_t m_size;
      size_t m_storageSize;
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

#include <brick/numeric/array2D_impl.hh>

#endif /* #ifdef BRICK_NUMERIC_ARRAY2D_HH */
