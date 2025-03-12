/**
***************************************************************************
* @file brick/numeric/transform2D.hh
*
* Header file declaring Transform2D class template.
*
* Copyright (C) 2001-2011 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_TRANSFORM2D_HH
#define BRICK_NUMERIC_TRANSFORM2D_HH

#include <brick/numeric/array2D.hh>
#include <brick/numeric/vector2D.hh>

namespace brick {

  namespace numeric {

    // Forward declaration.
    template <class Type>
    class Transform2DFunctor;


    /**
     ** The Transform2D class template represents a homogeneous
     ** coordinate transform from one 2D coordinate system to another
     ** 2D coordinate system.
     **/
    template <class Type>
    class Transform2D {
    public:

      /**
       * Default constructor.  Initializes to identity.
       */
      inline
      Transform2D();


      /**
       * Build a Transform2D instance by explicitly setting element values
       * as if setting the elements of a 3x3 transformation matrix:
       *
       *    [[a00, a01, a02],
       *     [a10, a11, a12],
       *     [a20, a21, a22]]
       *
       * @param a00 The value of one element of the transformation matrix.
       * @param a01 The value of one element of the transformation matrix.
       * @param a02 The value of one element of the transformation matrix.
       * @param a10 The value of one element of the transformation matrix.
       * @param a11 The value of one element of the transformation matrix.
       * @param a12 The value of one element of the transformation matrix.
       * @param a20 The value of one element of the transformation matrix.
       * @param a21 The value of one element of the transformation matrix.
       * @param a22 The value of one element of the transformation matrix.
       * @param doNormalize If true, the matrix will be rescaled so that its
       * lower right element is 1.0.
       */
      inline
      Transform2D(Type const& a00, Type const& a01, Type const& a02,
                  Type const& a10, Type const& a11, Type const& a12,
                  Type const& a20, Type const& a21, Type const& a22,
                  bool doNormalize = false);


      /**
       * Build a Transform2D instance from a sequence that specifies
       * element values in row major order.
       *
       *   [a0, a1, a2,
       *    a3, a4, a5,
       *    a6, a7, a8]]
       *
       * @param sequence A sequence of 8 or more elements that
       * specifies the elements of the transform.  If the sequence has
       * only 8 elements, then a value of Type(1.0) is used for the
       * 9th element.
       *
       * @param doNormalize If true, the matrix will be rescaled so that its
       * lower right element is 1.0.
       */
      Transform2D(std::initializer_list<Type> sequence,
                  bool doNormalize = false);


      /**
       * Build a Transform2D from a homogeneous 3x3 matrix.
       *
       * @param source A 2D array containing the elements of the desired
       * homogeneous transformation.
       * @param doNormalize If true, the resulting transform will be
       * rescaled so that its lower right element is 1.0.
       */
      Transform2D(Array2D<Type> const& source, bool doNormalize = false);


      /**
       * The copy constructor simply duplicates its argument.
       *
       * @param src This is the Transform2D instance to be copied.
       */
      inline
      Transform2D(Transform2D<Type> const& src);


      /**
       * Destructor.
       */
      ~Transform2D() {}


      /**
       * This member function returns a functor which makes it easier to
       * transform arrays of points using algorithms such as
       * std::transform().  For example:
       *
       * @code
       * std::transform(myPoints.begin(), myPoints.end(), myNewPoints.begin(),
       *                myTransform.getFunctor());
       * @endcode
       *
       * @return The return value is a functor instance which will
       * transform points according to the current value of *this.  The
       * functor will contain a copy of *this, so that subsequent
       * changes to *this will not affect the functor.
       */
      Transform2DFunctor<Type>
      getFunctor() const;


      /**
       * This member function returns one element from the matrix
       * representation of the coordinate transform by value.  For
       * example, calling myTransform2D.getValue<1, 2>() will return
       * the element from the 2nd row, 3rd column of the matrix
       * representation of the coordinate transformation.
       *
       * @return The value of the requested element.
       */
      template <size_t row, size_t column>
      Type const&
      getValue() const;


      /**
       * This member function returns the inverse of *this.  It is an
       * error if *this is not invertible.
       *
       * @return The return value is the inverse of *this.
       */
      Transform2D<Type>
      invert() const;


      /**
       * Change the Transform2D value by explicitly setting element values
       * as if setting the elements of a 3x3 transformation matrix:
       *    [[a00, a01, a02],
       *     [a10, a11, a12],
       *     [a20, a21, a22]]
       *
       * @param a00 The value of one element of the transformation matrix.
       * @param a01 The value of one element of the transformation matrix.
       * @param a02 The value of one element of the transformation matrix.
       * @param a10 The value of one element of the transformation matrix.
       * @param a11 The value of one element of the transformation matrix.
       * @param a12 The value of one element of the transformation matrix.
       * @param a20 The value of one element of the transformation matrix.
       * @param a21 The value of one element of the transformation matrix.
       * @param a22 The value of one element of the transformation matrix.
       * @param doNormalize If true, the matrix will be rescaled so that its
       * lower right element is 1.0.
       */
      void
      setTransform(Type const& a00, Type const& a01, Type const& a02,
                   Type const& a10, Type const& a11, Type const& a12,
                   Type const& a20, Type const& a21, Type const& a22,
                   bool doNormalize = false);


      /**
       * This member function sets one element of the matrix
       * representation of the coordinate transform.  Note that indexing
       * is zero-based.  For example, calling myTransform2D.setValue(1,
       * 2, 5.0) will set the element from the 2nd row, 3rd column of
       * the matrix representation of the coordinate transformation to
       * 5.0.  If blindingly fast execution is important, consider using
       * the setValue<size_t, size_t>(Type) member function instead.
       *
       * @param row This argument specifies the row of the matrix
       * element to be modified.
       *
       * @param column This argument specifies the row of the matrix
       * element to be modified.
       *
       * @param value This argument specifies the value to be copied
       * into the matrix element at the specified row & column.
       */
      void
      setValue(size_t row, size_t column, Type const& value);


      /**
       * This member function sets one element from the matrix
       * representation of the coordinate transform.  Note that indexing
       * is zero-based.  For example, calling myTransform2D.setValue<1,
       * 2>(5.0) will set the element from the 2nd row, 3rd column of
       * the matrix representation of the coordinate transformation to
       * 5.0.
       *
       * @param value This argument specifies the value to be copied
       * into the matrix element at the specified row & column.
       */
      template <size_t row, size_t column>
      void
      setValue(Type const& value);


      /**
       * This member function returns one element from the matrix
       * representation of the coordinate transform by value.
       * For example, calling myTransform2D.value<1, 2>() will return
       * the element from the 2nd row, 3rd column of the matrix
       * representation of the coordinate transformation.
       *
       * @return The value of the requested element.
       */
      template <size_t row, size_t column>
      inline Type const&
      value() const {return this->getValue<row, column>();}


      /**
       * This operator returns one element from the matrix
       * representation of the coordinate transform by value.
       * If blindingly fast execution is important, consider using
       * value<size_t, size_t>() member function.
       *
       * @param row The row of the requested element.
       * @param column The column of the requested element.
       * @return The value of the requested element.
       */
      Type const&
      operator()(size_t row, size_t column) const;

      /**
       * This operator takes a point and applies the coordinate
       * transform, returning the result.
       *
       * @param vector0 The point to be transformed.
       * @return The result of the transformation.
       */
      Vector2D<Type>
      operator*(const Vector2D<Type>& vector0) const;

      /**
       * The assignment operator simply duplicates its argument.
       *
       * @param source This is the Transform2D instance to be copied.
       * @return A reference to *this.
       */
      Transform2D<Type>&
      operator=(Transform2D<Type> const& source);

    private:
      void normalize();

      Type m_00, m_01, m_02;
      Type m_10, m_11, m_12;
      Type m_20, m_21, m_22;

    }; // class Transform2D


    /* ============== Helper classes ============== */


    /**
     ** This helper class works with Transform2D<Type>::getFunctor()
     **/
    template <class Type>
    class Transform2DFunctor
    {

    public:

      /**
       * The constructor deep-copies its argument.
       *
       *  @param transform This is the transform instance to be copied.
       */
      Transform2DFunctor(Transform2D<Type> const& transform)
        : m_transform(transform) {}


      /**
       * The application operator transforms its argument.
       *
       *  @param vec This is the vector to be transformed.
       *
       *  @return The transformed vector.
       */
      inline Vector2D<Type>
      operator()(Vector2D<Type> const& vec) const {return m_transform * vec;}

    private:
      Transform2D<Type> m_transform;
    };


    /* ============== Non-member functions ============== */

    /**
     * This operator composes two Transform2D instances.  The resulting
     * transform satisfies the equation:
     *
     *  (transform0 * transform1) * v0 = transform0 * (transform1 * v0),
     *
     * where v0 is a Vector2D instance.
     *
     * @param transform0 This is the first of the two Transform2D instances to
     * be composed.
     * @param transform1 This is the second of the two Transform2D instances to
     * be composed.
     * @return A Transform2D instance which is equivalent to the composition
     * of the two arguments.
     */
    template <class Type>
    Transform2D<Type>
    operator*(Transform2D<Type> const& transform0,
              Transform2D<Type> const& transform1);


    /**
     * Outputs a text representation of a Transform2D instance to a
     * std::ostream.  The output format looks like this:
     *
     * Transform2D(1.0, 2.0, 3.0,
     *             5.0, 6.0, 7.0,
     *             9.0, 10.0, 1.0)
     *
     * @param stream Reference to the the output stream.
     *
     * @param transform0 Const reference to the Transform2D instance to be
     * output.
     *
     * @return Reference to the output stream.
     */
    template <class Type>
    std::ostream&
    operator<<(std::ostream& stream, Transform2D<Type> const& transform0);


    /**
     * Sets the value of a Transform2D instance from a std::istream.
     * The input format is as described for
     * operator<<(std::ostream&, const Transform2D&), above.
     *
     * @param stream Reference to the the input stream.
     *
     * @param transform0 Reference to the Transform2D that will take
     * the input.
     *
     * @return Reference to the input stream.
     */
    template <class Type>
    std::istream&
    operator>>(std::istream& stream, Transform2D<Type>& transform0);

  } // namespace numeric

} // namespace brick

// Include file containing definitions of inline and template
// functions.
#include <brick/numeric/transform2D_impl.hh>

#endif /* #ifndef BRICK_NUMERIC_TRANSFORM2D_HH */
