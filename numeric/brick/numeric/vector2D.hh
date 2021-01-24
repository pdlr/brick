/**
***************************************************************************
* @file brick/numeric/vector2D.hh
*
* Header file declaring Vector2D class template.
*
* Copyright (C) 2001-2011 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_VECTOR2D_HH
#define BRICK_NUMERIC_VECTOR2D_HH

#include <iostream>
#include <brick/common/exception.hh>

namespace brick {

  namespace numeric {

    /**
     ** The Vector2D class template represents a 2D vector in which
     ** the elements are of a user specified type.
     **/
    template <class Type>
    class Vector2D {
    public:

      /**
       * The default constructor initializes to (0, 0).
       */
      inline
      Vector2D();


      /**
       * This constructor explicitly sets the 2D coordinates.
       *
       * @param xCoord The first component of the Vector2D.
       *
       * @param yCoord The second component of the Vector2D.
       */
      inline
      Vector2D(Type const& xCoord, Type const& yCoord);


      /**
       * This constructor explicitly sets 2D homogeneous coordinates.
       * The equivalent 2D coordinates are (xCoord/alpha, yCoord/alpha).
       *
       * @param xCoord The first component of the Vector2D.
       *
       * @param yCoord The second component of the Vector2D.
       *
       * @param alpha The projective scale parameter.
       */
      inline
      Vector2D(Type const& xCoord, Type const& yCoord, Type const& alpha);


      /**
       * The copy constructor deep copies its argument.
       *
       * @param other This argument is the Vector2D instance to be
       * copied.
       */
      inline
      Vector2D(Vector2D const& vec);


      /**
       * The destructor destroys the Vector2D instance.
       */
      inline
      ~Vector2D();


      /**
       * Resets the vector to (0.0, 0.0).
       *
       * @return A reference to *this.
       */
      inline Vector2D&
      clear();


      /**
       * This member function returns the X component of the Vector2D.
       *
       * @return The return value is the X coordinate.
       */
      inline Type const&
      getX() const;


      /**
       * This member function returns the Y component of the Vector2D.
       *
       * @return The return value is the Y coordinate.
       */
      inline Type const&
      getY() const;


      /**
       * This member function explicitly sets the sets coordinates of
       * the Vector2D instance.
       *
       * @param xCoord The first component of the Vector2D.
       *
       * @param yCoord The second component of the Vector2D.
       *
       * @return A reference to the updated Vector2D instance, *this.
       */
      inline Vector2D<Type> const&
      setValue(Type const& xCoord, Type const& yCoord);


      /**
       * This member function sets the X component of the Vector2D.
       *
       * @param newX This parameter is the value to which the X
       * coordinate will be set.
       *
       * @return The return value is the new X coordinate.
       */
      inline Type const&
      setX(Type const& newX);


      /**
       * This member function sets the Y component of the Vector2D.
       *
       * @param newY This parameter is the value to which the Y
       * coordinate will be set.
       *
       * @return The return value is the new Y coordinate.
       */
      inline Type const&
      setY(Type const& newY);


      /**
       * This member function explicitly sets 2D homogeneous
       * coordinates.  The equivalent 2D coordinates are
       * (xCoord/alpha, yCoord/alpha).
       *
       * @param xCoord The first component of the Vector2D.
       *
       * @param yCoord The second component of the Vector2D.
       *
       * @param alpha The projective scale parameter.
       *
       * @return A reference to the updated Vector2D instance, *this.
       */
      inline Vector2D<Type> const&
      setValue(Type const& xCoord, Type const& yCoord, Type const& alpha);

#if 0
      /**
       * This deprecated member function returns the X component of
       * the Vector2D by reference.
       *
       * @return The return value is a reference to the X coordinate.
       */
      inline Type&
      x() {return m_x;}
#endif /* #if 0 */

      /**
       * This member function returns the X component of the Vector2D
       * by value.  It is a synonym for memer function getX().
       *
       * @return The return value the X coordinate.
       */
      inline Type const& x() const {return m_x;}


#if 0
      /**
       * This deprecated member function returns the Y component of
       * the Vector2D by value.
       *
       * @return The return value the Y coordinate.
       */
      inline Type& y() {return m_y;}
#endif /* #if 0 */


      /**
       * This member function returns the Y component of the Vector2D
       * by value.  It is a synonym for memer function getY().
       *
       * @return The return value the Y coordinate.
       */
      inline Type const& y() const {return m_y;}


      /**
       * The assignment operator deep copies its argument.
       *
       * @param other This argument is the Vector2D instance to be
       * copied.
       *
       * @return The return value is a reference to *this after the
       * operation has been performed.
       */
      inline Vector2D&
      operator=(Vector2D const& vec);


      /**
       * The indexing operator returns a reference to the x, or y
       * component of *this as if *this were a two element array.
       * Out-of-bounds indices will return this->y().
       *
       * @param index This argument is the index into *this.
       *
       * @return The return value is the selected component of *this.
       */
      inline Type&
      operator[](size_t index);


      /**
       * The indexing operator returns the value of the x, or y
       * component of *this as if *this were a two element array.
       * Out of bounds indices will return this-y().
       *
       * @param index This argument is the index into *this.
       *
       * @return The return value is the selected component of *this.
       */
      inline Type
      operator[](size_t index) const;


      /**
       * This operator multiplies each component of the Vector2D instance
       * by a scalar.
       *
       * @param scalar This argument is the scalar by which to multiply.
       *
       * @return The return value is a reference to *this after the
       * operation has been performed.
       */
      inline Vector2D&
      operator*=(Type const& scalar);


      /**
       * This operator divides each component of the Vector2D instance
       * by a scalar.
       *
       * @param scalar This argument is the scalar by which to divide.
       *
       * @return The return value is a reference to *this after the
       * operation has been performed.
       */
      inline Vector2D&
      operator/=(Type const& scalar);


      /**
       * This operator adds a scalar to each component of the Vector2D
       * instance.
       *
       * @param scalar This argument is the scalar to be added.
       *
       * @return The return value is a reference to *this after the
       * operation has been performed.
       */
      inline Vector2D&
      operator+=(Vector2D const& vec);


      /**
       * This operator subtracts a scalar from each component of the
       * Vector2D instance.
       *
       * @param scalar This argument is the scalar to be subtracted.
       *
       * @return The return value is a reference to *this after the
       * operation has been performed.
       */
      inline Vector2D&
      operator-=(Vector2D const& vec);


      /**
       * This operator returns a Vector2D equal to *this, but with each
       * element negated.
       *
       * @return The return value is a negated copy of *this.
       */
      inline Vector2D
      operator-();


    private:
      // Private member functions.
      void
      normalize(Type const& alpha);

      // Private data members.
      Type m_x;
      Type m_y;
    }; // class Vector2D


    /* ============== Non-member function declarations ============== */

    /**
     * This operator returns the elementwise sum of two Vector2D instances.
     *
     * @param vector0 This is the first of the two Vector2D instances to
     * be added.
     * @param vector1 This is the second of the two Vector2D instances to
     * be added.
     * @return A Vector2D instance in which the value of each element is
     * equal to the sum of the corresponding elements of the two arguments.
     */
    template <class Type>
    Vector2D<Type>
    operator+(Vector2D<Type> const& vector0, Vector2D<Type> const& vector1);

    /**
     * This operator returns the elementwise difference of two Vector2D
     * instances.
     *
     * @param vector0 This is the first of the two Vector2D instances to
     * be subtracted.
     * @param vector1 This is the second of the two Vector2D instances to
     * be subtracted.
     * @return A Vector2D instance in which the value of each element is
     * equal to the difference of the corresponding elements of the two
     * arguments.
     */
    template <class Type>
    Vector2D<Type>
    operator-(Vector2D<Type> const& vector0, Vector2D<Type> const& vector1);

    /**
     * This operator returns the elementwise product of two Vector2D instances.
     *
     * @param vector0 This is the first of the two Vector2D instances to
     * be multiplied.
     * @param vector1 This is the second of the two Vector2D instances to
     * be multiplied.
     * @return A Vector2D instance in which the value of each element is
     * equal to the product of the corresponding elements of the two arguments.
     */
    template <class Type>
    Vector2D<Type>
    operator*(Vector2D<Type> const& vector0, Vector2D<Type> const& vector1);

    /**
     * This operator returns the elementwise dividend of two Vector2D instances.
     *
     * @param vector0 This is the Vector2D instance whose element values
     * are to be divided.
     * @param vector1 This is the Vector2D instance by whose elements
     * the first argument is to be divided.
     * @return A Vector2D instance in which the value of each element is
     * equal to the corresponding value of the first argument divided by
     * the corresponding value of the second argument.
     */
    template <class Type>
    Vector2D<Type>
    operator/(Vector2D<Type> const& vector0, Vector2D<Type> const& vector1);

    /**
     * This operator adds a scalar and a Vector2D.
     *
     * @param vector0 This is the Vector2D instance to which the scalar
     * should be added.
     * @param scalar0 This is amount which should be added to each
     * element of argument vector0.
     * @return A Vector2D instance in which the value of each element is
     * equal to the corresponding value of the first argument plus the
     * value of the second argument.
     */
    template <class Type>
    Vector2D<Type> operator+(Vector2D<Type> const& vector0,
                             Type const& scalar0);

    /**
     * This operator subtracts a scalar from a Vector2D.
     *
     * @param vector0 This is the Vector2D instance from which the scalar
     * should be subtracted.
     * @param scalar0 This is amount which should be subtracted from each
     * element of argument vector0.
     * @return A Vector2D instance in which the value of each element is
     * equal to the corresponding value of the first argument minus the
     * value of the second argument.
     */
    template <class Type>
    Vector2D<Type> operator-(Vector2D<Type> const& vector0,
                             Type const& scalar0);

    /**
     * This operator multiplies a Vector2D by scalar.
     *
     * @param vector0 This is the Vector2D instance which is to be
     * multiplied by the scalar.
     * @param scalar0 This is amount by which should argument vector0 is
     * to be multiplied.
     * @return A Vector2D instance in which the value of each element is
     * equal to the corresponding value of the first argument multiplied by
     * the value of the second argument.
     */
    template <class Type>
    Vector2D<Type> operator*(Vector2D<Type> const& vector0,
                             Type const& scalar0);

    /**
     * This operator divides a Vector2D by scalar.
     *
     * @param vector0 This is the Vector2D instance which is to be
     * divided by the scalar.
     * @param scalar0 This is amount by which should argument vector0 is
     * to be divided.
     * @return A Vector2D instance in which the value of each element is
     * equal to the corresponding value of the first argument divided by
     * the value of the second argument.
     */
    template <class Type>
    Vector2D<Type> operator/(Vector2D<Type> const& vector0,
                             Type const& scalar0);

    /**
     * This operator checks the supplied vectors for equality.
     *
     * @param  vector0  First vector to compare.
     * @param  vector1  Second vector to compare.
     * @return  Result of comparing @p vector0 to @p vector1 for equality.
     */
    template <class Type>
    bool operator==(Vector2D<Type> const& vector0,
                    Vector2D<Type> const& vector1);

    /**
     * This operator checks the supplied vectors for inequality.
     *
     * @param  vector0  First vector to compare.
     * @param  vector1  Second vector to compare.
     * @return  Result of comparing @p vector0 to @p vector1 for
     *          inequality.
     */
    template <class Type>
    bool operator!=(Vector2D<Type> const& vector0,
                    Vector2D<Type> const& vector1);


    /**
     * This operator adds a scalar value to each element of a Vector2D
     * instance.
     *
     * @param scalar0 Scalar argument of the addition.
     *
     * @param vector0 Vector argument of the addition.
     *
     * @return Vector2D instance in which the value of each element is
     * the sum of the scalar argument and the corresponding element of
     * the Vector2D argument.
     */
    template <class Type>
    Vector2D<Type> operator+(Type const& scalar0,
                             Vector2D<Type> const& vector0);


    /**
     * This operator multiplies a scalar value with each element of a
     * Vector2D instance.
     *
     * @param scalar0 Scalar argument of the multiplication.
     *
     * @param vector0 Vector argument of the multiplication.
     *
     * @return Vector2D instance in which the value of each element is
     * the product of the scalar argument and the corresponding element
     * of the Vector2D argument.
     */
    template <class Type>
    Vector2D<Type> operator*(Type const& scalar0,
                             Vector2D<Type> const& vector0);


    /**
     * This function outputs a text representation of a Vector2D
     * instance to a std::ostream.  The output format looks like this:
     *
     * Vector2D(1.0, 2.0)
     *
     * @param stream This argument is a reference to the the output
     * stream.
     *
     * @param vector0 This argument is a const reference to the
     * Vector2D instance to be output.
     *
     * @return The return value is a reference to the input stream after
     * the write has taken place.
     */
    template <class Type>
    std::ostream& operator<<(std::ostream& stream,
                             Vector2D<Type> const& vector0);

    /**
     * This function sets the value of a Vector2D instance from a
     * std::istream.  The input format is as described for
     * operator<<(std::ostream&, Vector2D const&) above.
     *
     * @param stream This argument is a reference to the the input
     * stream from which to read.
     *
     * @param vector0 This argument is a reference to the Vector2D
     * which will take the input.
     *
     * @return The return value is a reference to the input stream after
     * the read has taken place.
     */
    template <class Type>
    std::istream& operator>>(std::istream& stream, Vector2D<Type>& vector0);

  } // namespace numeric

} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/numeric/vector2D_impl.hh>

#endif /* #ifndef BRICK_NUMERIC_VECTOR2D_HH */
