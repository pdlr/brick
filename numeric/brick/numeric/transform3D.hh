/**
***************************************************************************
* @file brick/numeric/transform3D.hh
*
* Header file declaring Transform3D class template.
*
* Copyright (C) 2001-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_TRANSFORM3D_HH
#define BRICK_NUMERIC_TRANSFORM3D_HH

#include <brick/numeric/array2D.hh>
#include <brick/numeric/vector3D.hh>

namespace brick {

  namespace numeric {
    
    // Forward declaration.
    template <class Type>
    class Transform3DFunctor;

  
    /**
     ** The Transform3D class represents a homogeneous coordinate
     ** transform from one 3D coordinate system to another 3D
     ** coordinate system.
     **/
    template <class Type>
    class Transform3D {
    public:

      /** 
       * Default constructor.  Initializes to identity.
       */
      inline
      Transform3D();

      
      /** 
       * Build a Transform3D instance by explicitly setting element values
       * as if setting the elements of a 4x4 transformation matrix:
       *    [[a00, a01, a02, a03],
       *     [a10, a11, a12, a13],
       *     [a20, a21, a22, a23],
       *     [a30, a31, a32, a33]]
       * 
       * @param a00 The value of one element of the transformation matrix.
       * @param a01 The value of one element of the transformation matrix.
       * @param a02 The value of one element of the transformation matrix.
       * @param a03 The value of one element of the transformation matrix.
       * @param a10 The value of one element of the transformation matrix.
       * @param a11 The value of one element of the transformation matrix.
       * @param a12 The value of one element of the transformation matrix.
       * @param a13 The value of one element of the transformation matrix.
       * @param a20 The value of one element of the transformation matrix.
       * @param a21 The value of one element of the transformation matrix.
       * @param a22 The value of one element of the transformation matrix.
       * @param a23 The value of one element of the transformation matrix.
       * @param a30 The value of one element of the transformation matrix.
       * @param a31 The value of one element of the transformation matrix.
       * @param a32 The value of one element of the transformation matrix.
       * @param a33 The value of one element of the transformation matrix.
       * @param doNormalize If true, the matrix will be rescaled so that its
       * lower right element is 1.0.
       */
      inline
      Transform3D(
        Type const& a00, Type const& a01, Type const& a02, Type const& a03,
        Type const& a10, Type const& a11, Type const& a12, Type const& a13,
        Type const& a20, Type const& a21, Type const& a22, Type const& a23,
        Type const& a30, Type const& a31, Type const& a32, Type const& a33,
        bool doNormalize = false);


      /** 
       * Build a Transform3D from a homogeneous 4x4 matrix.
       *
       * @param source A 2D array containing the elements of the desired
       * homogeneous transformation.
       */
      Transform3D(const Array2D<Type>& source);


      /** 
       * The copy constructor simply duplicates its argument.
       * 
       * @param src This is the Transform3D instance to be copied.
       */
      inline
      Transform3D(const Transform3D& src);


      /** 
       * Destructor.
       */
      ~Transform3D() {}


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
      Transform3DFunctor<Type>
      getFunctor() const;
    

      /** 
       * This member function returns one element from the matrix
       * representation of the coordinate transform by value.  Note
       * that indexing is zero-based.  For example, calling
       * myTransform3D.value<1, 2>() will return the element from the
       * 2nd row, 3rd column of the matrix representation of the
       * coordinate transformation.  Indexing using this function is
       * faster than operator()(), but the indices must be known at
       * compile time.  This is an alias for member function template
       * value<row, column>().
       * 
       * @return The value of the requested element.
       */
      template <size_t row, size_t column>
      inline
      Type const&
      getValue() const;


      /** 
       * This member function returns the inverse of *this.  It is an
       * error if *this is not invertible.
       * 
       * @return The return value is the inverse of *this.
       */
      Transform3D
      invert() const;

    
      /** 
       * Change the Transform3D value by explicitly setting element values
       * as if setting the elements of a 4x4 transformation matrix:
       *    [[a00, a01, a02, a03],
       *     [a10, a11, a12, a13],
       *     [a20, a21, a22, a23],
       *     [a30, a31, a32, a33]]
       * 
       * @param a00 The value of one element of the transformation matrix.
       * @param a01 The value of one element of the transformation matrix.
       * @param a02 The value of one element of the transformation matrix.
       * @param a03 The value of one element of the transformation matrix.
       * @param a10 The value of one element of the transformation matrix.
       * @param a11 The value of one element of the transformation matrix.
       * @param a12 The value of one element of the transformation matrix.
       * @param a13 The value of one element of the transformation matrix.
       * @param a20 The value of one element of the transformation matrix.
       * @param a21 The value of one element of the transformation matrix.
       * @param a22 The value of one element of the transformation matrix.
       * @param a23 The value of one element of the transformation matrix.
       * @param a30 The value of one element of the transformation matrix.
       * @param a31 The value of one element of the transformation matrix.
       * @param a32 The value of one element of the transformation matrix.
       * @param a33 The value of one element of the transformation matrix.
       * @param doNormalize If true, the matrix will be rescaled so that its
       * lower right element is 1.0.
       */
      void
      setTransform(
        Type const& a00, Type const& a01, Type const& a02, Type const& a03,
        Type const& a10, Type const& a11, Type const& a12, Type const& a13,
        Type const& a20, Type const& a21, Type const& a22, Type const& a23,
        Type const& a30, Type const& a31, Type const& a32, Type const& a33,
        bool doNormalize = false);



      /** 
       * This member function sets one element of the matrix
       * representation of the coordinate transform.  Note that indexing
       * is zero-based.  For example, calling myTransform3D.setValue(1,
       * 2, 5.0) will set the element from the 2nd row, 3rd column of
       * the matrix representation of the coordinate transformation to
       * 5.0.  If blindingly fast execution is important, consider using
       * the setValue<size_t, size_t>(Type) member function instead.
       *
       * It is not permitted to set the value of the lower-right element
       * of the matrix: setValue(3, 3, ...) will throw an
       * IndexException.
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
       * is zero-based.  For example, calling myTransform3D.setValue<1,
       * 2>(5.0) will set the element from the 2nd row, 3rd column of
       * the matrix representation of the coordinate transformation to
       * 5.0.
       *
       * It is not permitted to set the value of the lower-right element
       * of the matrix: setValue<3, 3>(Type) will throw an
       * IndexException.
       * 
       * @param value This argument specifies the value to be copied
       * into the matrix element at the specified row & column.
       */
      template <size_t row, size_t column>
      void
      setValue(Type const& value);


      /** 
       * This member function is an alias for member function setTransform().
       */
      void
      setValue(
        Type const& a00, Type const& a01, Type const& a02, Type const& a03,
        Type const& a10, Type const& a11, Type const& a12, Type const& a13,
        Type const& a20, Type const& a21, Type const& a22, Type const& a23,
        Type const& a30, Type const& a31, Type const& a32, Type const& a33,
        bool doNormalize = false);

      
      /** 
       * This member function returns one element from the matrix
       * representation of the coordinate transform by value.  Note that
       * indexing is zero-based.  For example, calling
       * myTransform3D.value<1, 2>() will return the element from the
       * 2nd row, 3rd column of the matrix representation of the
       * coordinate transformation.
       * 
       * @return The value of the requested element.
       */
      template <size_t row, size_t column>
      Type const&
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
      Vector3D<Type>
      operator*(const Vector3D<Type>& vector0) const;

    
      /** 
       * The assignment operator simply duplicates its argument.
       * 
       * @param source This is the Transform3D instance to be copied.
       * @return A reference to *this.
       */
      Transform3D<Type>&
      operator=(const Transform3D& source);

    
    private:
      void normalize();
      
      Type m_00, m_01, m_02, m_03;
      Type m_10, m_11, m_12, m_13;
      Type m_20, m_21, m_22, m_23;
      Type m_30, m_31, m_32, m_33;

    }; // class Transform3D


    /* ============== Helper classes ============== */


    /**
     ** This helper class works with Transform3D::getFunctor()
     **/
    template <class Type>
    class Transform3DFunctor
      : public std::unary_function< Vector3D<Type>, Vector3D<Type> > {

    public:

      /** 
       * The constructor deep-copies its argument.
       * 
       *  @param transform This is the transform instance to be copied.
       */
      Transform3DFunctor(Transform3D<Type> const& transform)
        : m_transform(transform) {}

    
      /** 
       * The application operator transforms its argument.
       * 
       *  @param vec This is the vector to be transformed.
       * 
       *  @return The transformed vector.
       */
      inline Vector3D<Type>
      operator()(Vector3D<Type> const& vec) const {return m_transform * vec;}
    
    private:
      Transform3D<Type> m_transform;
    };
  
  
    /* ============== Non-member functions ============== */
  
    /** 
     * This operator composes two Transform3D instances.  The resulting
     * transform satisfies the equation:
     *
     *  (transform0 * transform1) * v0 = transform0 * (transform1 * v0),
     *
     * where v0 is a Vector3D instance.
     * 
     * @param transform0 This is the first of the two Transform3D instances to
     * be composed.
     * @param transform1 This is the second of the two Transform3D instances to
     * be composed.
     * @return A Transform3D instance which is equivalent to the composition
     * of the two arguments.
     */
    template <class Type>
    Transform3D<Type>
    operator*(const Transform3D<Type>& transform0, const Transform3D<Type>& transform1);

  
    /** 
     * Outputs a text representation of a Transform3D instance to a
     * std::ostream.  The output format looks like this:
     *
     * Transform3D(1.0, 2.0, 3.0, 4.0,
     *             5.0, 6.0, 7.0, 8.0,
     *             9.0, 10.0, 11.0, 12.0,
     *             13.0, 14.0, 15.0, 1.0)
     *
     * @param stream Reference to the the output stream.
     *
     * @param transform0 Const reference to the Transform3D instance to be
     * output.
     *
     * @return Reference to the output stream.
     */
    template <class Type>
    std::ostream&
    operator<<(std::ostream& stream, const Transform3D<Type>& transform0);

  
    /** 
     * Sets the value of a Transform3D instance from a std::istream.
     * The input format is as described for
     * operator<<(std::ostream&, const Transform3D&), above.
     * 
     * @param stream Reference to the the input stream.
     *
     * @param transform0 Reference to the Transform3D that will take
     * the input.
     *
     * @return Reference to the input stream.
     */
    template <class Type>
    std::istream&
    operator>>(std::istream& stream, Transform3D<Type>& transform0);
  
  } // namespace numeric

} // namespace brick

#include <brick/numeric/transform3D_impl.hh>

#endif /* #ifndef BRICK_NUMERIC_TRANSFORM3D_HH */
