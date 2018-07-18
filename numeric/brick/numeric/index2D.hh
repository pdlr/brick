/**
***************************************************************************
* @file brick/numeric/index2D.hh
*
* Header file declaring Index2D class.
*
* Copyright (C) 2001-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_INDEX2D_HH
#define BRICK_NUMERIC_INDEX2D_HH

#include <iostream>
#include <brick/common/exception.hh>

namespace brick {

  namespace numeric {

    /**
     ** The Index2D class represents a 2 dimensional index in (row,
     ** column) format, such as (0, 1), (23, 7), or (-4, 2).  An
     ** alternate way to think of this class is as a point with
     ** integer coordinates in a (U, V) coordinate system, where the
     ** origin lies at row = 0, column = 0, the U axis points along
     ** the first row, and the V axis points along the first column.
     ** Note that by this definition, the value of U is equal to the
     ** column index, and the value of V is equal to the row index.
     ** This leads to some confusing ordering issues.  For example,
     ** the value returned by Index2D(1, 2).getU() is 2, and the value
     ** returned by Index2D(1, 2).getV() is 1.
     **
     ** WARNING: be particularly careful with
     ** Index2D::operator[](size_t), as illustrated by this example,
     ** in which all assertions pass:
     **
     ** @code
     **   Index2D exampleIndex(1, 2);
     **   assert(exampleIndex[0] == exampleIndex.getU());
     **   assert(exampleIndex[0] == 2);
     **   assert(exampleIndex[1] == exampleIndex.getV());
     **   assert(exampleIndex[1] == 1);
     ** @endcode
     **/
    class Index2D {
    public:

      /**
       * The default constructor initializes to (0, 0).
       */
      Index2D()
        : m_column(0), m_row(0) {}


      /**
       * This constructor explicitly sets the indices.
       *
       * @param uIndex The first component of the Index2D.
       *
       * @param uIndex The second component of the Index2D.
       */
      Index2D(int row, int column)
        : m_column(column), m_row(row) {}


      /**
       * The copy constructor deep copies its argument.
       *
       * @param other This argument is the Index2D instance to be
       * copied.
       */
      Index2D(const Index2D& other)
        : m_column(other.m_column), m_row(other.m_row) {}


      /**
       * The destructor destroys the Index2D instance.
       */
      ~Index2D() {}


      /**
       * This member function explicitly sets the sets the indices.
       *
       * @param row The first component of the Index2D.
       *
       * @param column The second component of the Index2D.
       */
      inline void setValue(int row, int column) {
        m_column = column; m_row = row;
      }


      /**
       * This member function returns the first component of the Index2D
       * by value.
       *
       * @return The return value is the first (U) coordinate.
       */
      inline int getColumn() const {return m_column;}


      /**
       * This member function returns the second component of the
       * Index2D by value.
       *
       * @return The return value is the second (V) coordinate.
       */
      inline int getRow() const {return m_row;}


      /**
       * This member function returns the first component of the Index2D
       * by value.  It is a synonym for member function getColumn().
       *
       * @return The return value is the first (U) coordinate.
       */
      inline int getU() const {return this->getColumn();}


      /**
       * This member function returns the second component of the
       * Index2D by value.  It is a synonym for member function getRow().
       *
       * @return The return value is the second (V) coordinate.
       */
      inline int getV() const {return this->getRow();}


      /**
       * The assignment operator deep copies its argument.
       *
       * @param other This argument is the Index2D instance to be
       * copied.
       *
       * @return The return value is a reference to *this after the
       * operation has been performed.
       */
      Index2D& operator=(const Index2D& other) {
        // Self-assignment is OK.
        setValue(other.m_row, other.m_column); return *this;
      }


      /**
       * WARNING: this operator has confusing semantics, as
       * illustrated by this example, in which all assertions pass:
       *
       * @code
       *   Index2D exampleIndex(1, 2);
       *   assert(exampleIndex[0] == exampleIndex.getU());
       *   assert(exampleIndex[0] == 2);
       *   assert(exampleIndex[1] == exampleIndex.getU());
       *   assert(exampleIndex[1] == 1);
       * @endcode
       *
       * The indexing operator returns a reference to the U or V
       * component of *this as if *this were a two element array.
       * Out of bounds indices will return the V component.
       *
       * @param index This argument is the index into *this.  An index
       * of 0 means to return the "U" component, where the U axis is
       * parallel to the first row of the indexed array.  In other
       * words, an index of 0 means to return the column value.
       *
       * @return The return value is the selected component of *this.
       */
      int operator[](size_t index) const {
        switch(index) {
        case 0: return this->getColumn(); break;
        }
        return this->getRow();
      }


      /**
       * This operator multiplies each component of the Index2D instance
       * by a scalar.
       *
       * @param scalar This argument is the scalar by which to multiply.
       *
       * @return The return value is a reference to *this after the
       * operation has been performed.
       */
      Index2D& operator*=(int scalar) {
        m_column *= scalar; m_row *= scalar; return *this;
      }


      /**
       * This operator divides each component of the Index2D instance
       * by a scalar.
       *
       * @param scalar This argument is the scalar by which to divide.
       *
       * @return The return value is a reference to *this after the
       * operation has been performed.
       */
      Index2D& operator/=(int scalar) {
        if (scalar == 0 ) {
          BRICK_THROW(common::ValueException, "Index2D::operator/=()",
                      "Bad scalar: Divide by Zero\n");
        }
        m_column /= scalar; m_row /= scalar; return *this;
      }


      /**
       * This operator adds a scalar to each component of the Index2D
       * instance.
       *
       * @param scalar This argument is the scalar to be added.
       *
       * @return The return value is a reference to *this after the
       * operation has been performed.
       */
      Index2D& operator+=(const Index2D& other) {
        m_column += other.m_column; m_row += other.m_row; return *this;
      }


      /**
       * This operator subtracts a scalar from each component of the
       * Index2D instance.
       *
       * @param scalar This argument is the scalar to be subtracted.
       *
       * @return The return value is a reference to *this after the
       * operation has been performed.
       */
      Index2D& operator-=(const Index2D& other) {
        m_column -= other.m_column; m_row -= other.m_row; return *this;
      }


      /**
       * This operator returns an Index2D equal to *this, but with each
       * element negated.
       *
       * @return The return value is a negated copy of *this.
       */
      Index2D operator-() {
        return Index2D(-m_row, -m_column);
      }

    private:

      // Private data members.
      int m_column;
      int m_row;
    }; // class Index2D


    /* ============== Non-member function declarations ============== */

    /**
     * This operator returns the elementwise sum of two Index2D instances.
     *
     * @param index0 This is the first of the two Index2D instances to
     * be added.
     * @param index1 This is the second of the two Index2D instances to
     * be added.
     * @return An Index2D instance in which the value of each element is
     * equal to the sum of the corresponding elements of the two arguments.
     */
    Index2D
    operator+(const Index2D& index0, const Index2D& index1);

    /**
     * This operator returns the elementwise difference of two Index2D
     * instances.
     *
     * @param index0 This is the first of the two Index2D instances to
     * be subtracted.
     * @param index1 This is the second of the two Index2D instances to
     * be subtracted.
     * @return An Index2D instance in which the value of each element is
     * equal to the difference of the corresponding elements of the two
     * arguments.
     */
    Index2D
    operator-(const Index2D& index0, const Index2D& index1);

    /**
     * This operator returns the elementwise product of two Index2D instances.
     *
     * @param index0 This is the first of the two Index2D instances to
     * be multiplied.
     * @param index1 This is the second of the two Index2D instances to
     * be multiplied.
     * @return An Index2D instance in which the value of each element is
     * equal to the product of the corresponding elements of the two arguments.
     */
    Index2D
    operator*(const Index2D& index0, const Index2D& index1);

    /**
     * This operator returns the elementwise dividend of two Index2D instances.
     *
     * @param index0 This is the Index2D instance whose element values
     * are to be divided.
     * @param index1 This is the Index2D instance by whose elements
     * the first argument is to be divided.
     * @return An Index2D instance in which the value of each element is
     * equal to the corresponding value of the first argument divided by
     * the corresponding value of the second argument.
     */
    Index2D
    operator/(const Index2D& index0, const Index2D& index1);

    /**
     * This operator adds a scalar and an Index2D.
     *
     * @param index0 This is the Index2D instance to which the scalar
     * should be added.
     * @param scalar0 This is amount which should be added to each
     * element of argument index0.
     * @return An Index2D instance in which the value of each element is
     * equal to the corresponding value of the first argument plus the
     * value of the second argument.
     */
    Index2D operator+(const Index2D& index0, int scalar0);

    /**
     * This operator subtracts a scalar from an Index2D.
     *
     * @param index0 This is the Index2D instance from which the scalar
     * should be subtracted.
     * @param scalar0 This is amount which should be subtracted from each
     * element of argument index0.
     * @return An Index2D instance in which the value of each element is
     * equal to the corresponding value of the first argument minus the
     * value of the second argument.
     */
    Index2D operator-(const Index2D& index0, int scalar0);

    /**
     * This operator multiplies an Index2D by scalar.
     *
     * @param index0 This is the Index2D instance which is to be
     * multiplied by the scalar.
     * @param scalar0 This is amount by which should argument index0 is
     * to be multiplied.
     * @return An Index2D instance in which the value of each element is
     * equal to the corresponding value of the first argument multiplied by
     * the value of the second argument.
     */
    Index2D operator*(const Index2D& index0, int scalar0);

    /**
     * This operator divides an Index2D by scalar.
     *
     * @param index0 This is the Index2D instance which is to be
     * divided by the scalar.
     * @param scalar0 This is amount by which should argument index0 is
     * to be divided.
     * @return An Index2D instance in which the value of each element is
     * equal to the corresponding value of the first argument divided by
     * the value of the second argument.
     */
    Index2D operator/(const Index2D& index0, int scalar0);

    /**
     * This operator checks the supplied indexs for equality.
     *
     * @param  index0  First Index2D instance to compare.
     * @param  index1  Second Index2D instance to compare.
     * @return  Result of comparing @p index0 to @p index1 for equality.
     */
    bool operator==(const Index2D& index0, const Index2D& index1);

    /**
     * This operator checks the supplied indexs for inequality.
     *
     * @param  index0  First Index2D instance to compare.
     * @param  index1  Second Index2D instance to compare.
     * @return  Result of comparing @p index0 to @p index1 for
     *          inequality.
     */
    bool operator!=(const Index2D& index0, const Index2D& index1);


    /**
     * This operator adds a scalar value to each element of an Index2D
     * instance.
     *
     * @param scalar0 Scalar argument of the addition.
     *
     * @param index0 Index2D argument of the addition.
     *
     * @return Index2D instance in which the value of each element is
     * the sum of the scalar argument and the corresponding element of
     * the Index2D argument.
     */
    Index2D operator+(int scalar0, const Index2D& index0);


    /**
     * This operator multiplies a scalar value with each element of a
     * Index2D instance.
     *
     * @param scalar0 Scalar argument of the multiplication.
     *
     * @param index0 Index2D argument of the multiplication.
     *
     * @return Index2D instance in which the value of each element is
     * the product of the scalar argument and the corresponding element
     * of the Index2D argument.
     */
    Index2D operator*(int scalar0, const Index2D& index0);


    /**
     * This function outputs a text representation of an Index2D
     * instance to a std::ostream.  The output format looks like this:
     *
     * Index2D(1.0, 2.0)
     *
     * @param stream This argument is a reference to the the output
     * stream.
     *
     * @param index0 This argument is a const reference to the
     * Index2D instance to be output.
     *
     * @return The return value is a reference to the input stream after
     * the write has taken place.
     */
    std::ostream& operator<<(std::ostream& stream, const Index2D& index0);


    /**
     * This function sets the value of an Index2D instance from a
     * std::istream.  The input format is as described for
     * operator<<(std::ostream&, const Index2D&) above.
     *
     * @param stream This argument is a reference to the the input
     * stream from which to read.
     *
     * @param index0 This argument is a reference to the Index2D
     * which will take the input.
     *
     * @return The return value is a reference to the input stream after
     * the read has taken place.
     */
    std::istream& operator>>(std::istream& stream, Index2D& index0);

  } // namespace numeric

} // namespace brick


/* ============ Definitions of inline & template functions ============ */

namespace brick {

  namespace numeric {

    inline Index2D operator+(int scalar0, const Index2D& index0)
    {
      return index0 + scalar0;
    }

    inline Index2D operator*(int scalar0, const Index2D& index0)
    {
      return index0 * scalar0;
    }

  } // namespace numeric

} // namespace brick

#endif // #ifndef BRICK_NUMERIC_INDEX2D_HH
