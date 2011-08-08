/**
***************************************************************************
* @file brick/numeric/index3D.hh
*
* Header file declaring Index3D class.
*
* Copyright (C) 2001-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_INDEX3D_HH
#define BRICK_NUMERIC_INDEX3D_HH

#include <iostream>
#include <brick/common/exception.hh>

namespace brick {

  namespace numeric {
    /**
     ** The Index3D class represents a 3 dimensional index, such as
     ** (0, 1, 4), (23, 7, -11), or (-4, 2, 0).  By convention, we
     ** refer to these coordinates as (U, V, W).  We also refer to
     ** slice, row, and column, where column == U, row == V, and slice
     ** == W.
     **/
    class Index3D {
    public:

      /** 
       * The default constructor initializes to (0, 0, 0).
       */
      Index3D()
        : m_column(0), m_row(0), m_slice(0) {}
    

      /** 
       * This constructor explicitly sets the indices.
       * 
       * @param slice The third (W) component of the Index3D.
       * 
       * @param row The second (U) component of the Index3D.
       *
       * @param column The first (U) component of the Index3D.
       */
      Index3D(int slice, int row, int column)
        : m_column(column), m_row(row), m_slice(slice) {}


      /** 
       * The copy constructor deep copies its argument.
       * 
       * @param other This argument is the Index3D instance to be
       * copied.
       */
      Index3D(const Index3D& other)
        : m_column(other.m_column), m_row(other.m_row), m_slice(other.m_slice)
        {}


      /** 
       * The destructor destroys the Index3D instance.
       */
      ~Index3D() {}


      /** 
       * This member function explicitly sets the sets the indices.
       * 
       * @param slice The third (W) component of the Index3D.
       * 
       * @param row The second (U) component of the Index3D.
       *
       * @param column The first (U) component of the Index3D.
       */
      inline void setValue(int slice, int row, int column) {
        m_column = column; m_row = row; m_slice = slice;
      }


      /** 
       * This member function returns the first component of the Index3D
       * by value.
       * 
       * @return The return value is the first (U) coordinate.
       */
      inline int getColumn() const {return m_column;}


      /** 
       * This member function returns the second component of the
       * Index3D by value.
       * 
       * @return The return value is the second (V) coordinate.
       */
      inline int getRow() const {return m_row;}

    
      /** 
       * This member function returns the third component of the
       * Index3D by value.
       * 
       * @return The return value is the third (W) coordinate.
       */
      inline int getSlice() const {return m_slice;}

    
      /** 
       * This member function returns the first component of the Index3D
       * by value.  It is a synonym for member function getColumn().
       * 
       * @return The return value is the first (U) coordinate.
       */
      inline int getU() const {return this->getColumn();}


      /** 
       * This member function returns the second component of the
       * Index3D by value.  It is a synonym for member function getRow().
       * 
       * @return The return value is the second (V) coordinate.
       */
      inline int getV() const {return this->getRow();}

    
      /** 
       * This member function returns the third component of the
       * Index3D by value.  It is a synonym for member function
       * getSlice().
       * 
       * @return The return value is the third (W) coordinate.
       */
      inline int getW() const {return this->getSlice();}

    
      /** 
       * The assignment operator deep copies its argument.
       * 
       * @param other This argument is the Index3D instance to be
       * copied.
       *
       * @return The return value is a reference to *this after the
       * operation has been performed.
       */
      Index3D& operator=(const Index3D& other) {
        // Self-assignment is OK.
        this->setValue(other.m_slice, other.m_row, other.m_column);
        return *this;
      }


      
      /** 
       * The indexing operator returns a reference to the U, V, or W
       * component of *this as if *this were a three element array.
       * Out of bounds indices will return the W component.
       * 
       * @param index This argument is the index into *this.
       * 
       * @return The return value is the selected component of *this.
       */
      int operator[](size_t index) const {
        switch(index) {
        case 0: return this->getColumn(); break;
        case 1: return this->getRow(); break;
        }
        return this->getSlice();
      }
        

      /** 
       * This operator multiplies each component of the Index3D instance
       * by a scalar.
       * 
       * @param scalar This argument is the scalar by which to multiply.
       * 
       * @return The return value is a reference to *this after the
       * operation has been performed.
       */
      Index3D& operator*=(int scalar) {
        m_column *= scalar; m_row *= scalar; m_slice *= scalar; return *this;
      }


      /** 
       * This operator divides each component of the Index3D instance
       * by a scalar.
       * 
       * @param scalar This argument is the scalar by which to divide.
       * 
       * @return The return value is a reference to *this after the
       * operation has been performed.
       */
      Index3D& operator/=(int scalar) {
        if (scalar == 0 ) {
          BRICK_THROW(common::ValueException, "Index3D::operator/=()",
                      "Bad scalar: Divide by Zero\n");
        }
        m_column /= scalar; m_row /= scalar; m_slice /= scalar; return *this;
      }


      /** 
       * This operator adds a scalar to each component of the Index3D
       * instance.
       * 
       * @param scalar This argument is the scalar to be added.
       * 
       * @return The return value is a reference to *this after the
       * operation has been performed.
       */
      Index3D& operator+=(const Index3D& other) {
        m_column += other.m_column; m_row += other.m_row;
        m_slice += other.m_slice; return *this;
      }


      /** 
       * This operator subtracts a scalar from each component of the
       * Index3D instance.
       * 
       * @param scalar This argument is the scalar to be subtracted.
       * 
       * @return The return value is a reference to *this after the
       * operation has been performed.
       */
      Index3D& operator-=(const Index3D& other) {
        m_column -= other.m_column; m_row -= other.m_row;
        m_slice -= other.m_slice; return *this;
      }

    
      /** 
       * This operator returns an Index3D equal to *this, but with each
       * element negated.
       * 
       * @return The return value is a negated copy of *this.
       */
      Index3D operator-() {
        return Index3D(-m_row, -m_column, -m_slice);
      }
    
    private:

      // Private data members.
      int m_column;
      int m_row;
      int m_slice;
    }; // class Index3D


    /* ============== Non-member function declarations ============== */
  
    /** 
     * This operator returns the elementwise sum of two Index3D instances.
     * 
     * @param index0 This is the first of the two Index3D instances to
     * be added.
     * @param index1 This is the second of the two Index3D instances to
     * be added.
     * @return An Index3D instance in which the value of each element is
     * equal to the sum of the corresponding elements of the two arguments.
     */
    Index3D
    operator+(const Index3D& index0, const Index3D& index1);

    /** 
     * This operator returns the elementwise difference of two Index3D
     * instances.
     * 
     * @param index0 This is the first of the two Index3D instances to
     * be subtracted.
     * @param index1 This is the second of the two Index3D instances to
     * be subtracted.
     * @return An Index3D instance in which the value of each element is
     * equal to the difference of the corresponding elements of the two
     * arguments.
     */
    Index3D
    operator-(const Index3D& index0, const Index3D& index1);
  
    /** 
     * This operator returns the elementwise product of two Index3D instances.
     * 
     * @param index0 This is the first of the two Index3D instances to
     * be multiplied.
     * @param index1 This is the second of the two Index3D instances to
     * be multiplied.
     * @return An Index3D instance in which the value of each element is
     * equal to the product of the corresponding elements of the two arguments.
     */
    Index3D
    operator*(const Index3D& index0, const Index3D& index1);

    /** 
     * This operator returns the elementwise dividend of two Index3D instances.
     * 
     * @param index0 This is the Index3D instance whose element values
     * are to be divided.
     * @param index1 This is the Index3D instance by whose elements
     * the first argument is to be divided.
     * @return An Index3D instance in which the value of each element is
     * equal to the corresponding value of the first argument divided by
     * the corresponding value of the second argument.
     */
    Index3D
    operator/(const Index3D& index0, const Index3D& index1);

    /** 
     * This operator adds a scalar and an Index3D.
     * 
     * @param index0 This is the Index3D instance to which the scalar
     * should be added.
     * @param scalar0 This is amount which should be added to each
     * element of argument index0.
     * @return An Index3D instance in which the value of each element is
     * equal to the corresponding value of the first argument plus the
     * value of the second argument.
     */
    Index3D operator+(const Index3D& index0, int scalar0);

    /** 
     * This operator subtracts a scalar from an Index3D.
     * 
     * @param index0 This is the Index3D instance from which the scalar
     * should be subtracted.
     * @param scalar0 This is amount which should be subtracted from each
     * element of argument index0.
     * @return An Index3D instance in which the value of each element is
     * equal to the corresponding value of the first argument minus the
     * value of the second argument.
     */
    Index3D operator-(const Index3D& index0, int scalar0);

    /** 
     * This operator multiplies an Index3D by scalar.
     * 
     * @param index0 This is the Index3D instance which is to be
     * multiplied by the scalar.
     * @param scalar0 This is amount by which should argument index0 is
     * to be multiplied.
     * @return An Index3D instance in which the value of each element is
     * equal to the corresponding value of the first argument multiplied by
     * the value of the second argument.
     */
    Index3D operator*(const Index3D& index0, int scalar0);

    /** 
     * This operator divides an Index3D by scalar.
     * 
     * @param index0 This is the Index3D instance which is to be
     * divided by the scalar.
     * @param scalar0 This is amount by which should argument index0 is
     * to be divided.
     * @return An Index3D instance in which the value of each element is
     * equal to the corresponding value of the first argument divided by
     * the value of the second argument.
     */
    Index3D operator/(const Index3D& index0, int scalar0);

    /**
     * This operator checks the supplied indexs for equality.
     *
     * @param  index0  First Index3D instance to compare.
     * @param  index1  Second Index3D instance to compare.
     * @return  Result of comparing @p index0 to @p index1 for equality.
     */
    bool operator==(const Index3D& index0, const Index3D& index1);

    /**
     * This operator checks the supplied indexs for inequality.
     *
     * @param  index0  First Index3D instance to compare.
     * @param  index1  Second Index3D instance to compare.
     * @return  Result of comparing @p index0 to @p index1 for
     *          inequality.
     */
    bool operator!=(const Index3D& index0, const Index3D& index1);


    /** 
     * This operator adds a scalar value to each element of an Index3D
     * instance.
     * 
     * @param scalar0 Scalar argument of the addition.
     *
     * @param index0 Index3D argument of the addition.
     *
     * @return Index3D instance in which the value of each element is
     * the sum of the scalar argument and the corresponding element of
     * the Index3D argument.
     */
    Index3D operator+(int scalar0, const Index3D& index0);


    /** 
     * This operator multiplies a scalar value with each element of a
     * Index3D instance.
     * 
     * @param scalar0 Scalar argument of the multiplication.
     *
     * @param index0 Index3D argument of the multiplication.
     *
     * @return Index3D instance in which the value of each element is
     * the product of the scalar argument and the corresponding element
     * of the Index3D argument.
     */
    Index3D operator*(int scalar0, const Index3D& index0);

  
    /** 
     * This function outputs a text representation of an Index3D
     * instance to a std::ostream.  The output format looks like this:
     *
     * Index3D(1.0, 2.0)
     *
     * @param stream This argument is a reference to the the output
     * stream.
     *
     * @param index0 This argument is a const reference to the
     * Index3D instance to be output.
     *
     * @return The return value is a reference to the input stream after
     * the write has taken place.
     */
    std::ostream& operator<<(std::ostream& stream, const Index3D& index0);

  
    /** 
     * This function sets the value of an Index3D instance from a
     * std::istream.  The input format is as described for
     * operator<<(std::ostream&, const Index3D&) above.
     * 
     * @param stream This argument is a reference to the the input
     * stream from which to read.
     *
     * @param index0 This argument is a reference to the Index3D
     * which will take the input.
     *
     * @return The return value is a reference to the input stream after
     * the read has taken place.
     */
    std::istream& operator>>(std::istream& stream, Index3D& index0);

  } // namespace numeric

} // namespace brick


/* ============ Definitions of inline & template functions ============ */

namespace brick {

  namespace numeric {
    
    inline Index3D operator+(int scalar0, const Index3D& index0)
    {
      return index0 + scalar0;
    }
  
    inline Index3D operator*(int scalar0, const Index3D& index0)
    {
      return index0 * scalar0;
    }

  } // namespace numeric

} // namespace brick

#endif // #ifndef BRICK_NUMERIC_INDEX3D_HH
