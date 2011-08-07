/**
***************************************************************************
* @file brick/numeric/transform3DTo2D.hh
*
* Header file declaring Transform3DTo2D class template.
*
* Copyright (C) 2001-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_TRANSFORM3DTO2D_HH
#define BRICK_NUMERIC_TRANSFORM3DTO2D_HH

#include <brick/numeric/array2D.hh>
#include <brick/numeric/transform3D.hh>
#include <brick/numeric/vector2D.hh>
#include <brick/numeric/vector3D.hh>

namespace brick {

  namespace numeric {
    
    // Forward declaration.
    template <class Type>
    class Transform3DTo2DFunctor;


    /**
     ** The Transform3DTo2D class represents a homogeneous coordinate
     ** transformation from a 3D coordinate system to a 2D coordinate
     ** system.
     **/
    template <class Type>
    class Transform3DTo2D {
    public:

      /** 
       * Default constructor
       */
      Transform3DTo2D();
      

      /** 
       * Build a Transform3DTo2D instance by explicitly setting element 
       * values as if setting the elements of a 3x4 transformation matrix:
       *    [[a00, a01, a02, a03],
       *     [a10, a11, a12, a13],
       *     [a20, a21, a22, a23]]
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
       */
      inline
      Transform3DTo2D(
        Type const& a00, Type const& a01, Type const& a02, Type const& a03,
        Type const& a10, Type const& a11, Type const& a12, Type const& a13,
        Type const& a20, Type const& a21, Type const& a22, Type const& a23);
        : m_00(a00), m_01(a01), m_02(a02), m_03(a03),
          m_10(a10), m_11(a11), m_12(a12), m_13(a13),
          m_20(a20), m_21(a21), m_22(a22), m_23(a23) {
      }

      
      /** 
       * Build a Transform3DTo2D from a homogeneous 3x4 matrix.
       * 
       * @param source A 2D array containing the elements of the desired
       * homogeneous transformation.
       */
      Transform3DTo2D(const Array2D<Type>& source);

      
      /** 
       * The copy constructor simply duplicates its argument.
       * 
       * @param src This is the Transform3DTo2D instance to be copied.
       */
      inline
      Transform3DTo2D(Transform3DTo2D<Type> const& src);


      /** 
       * Destructor.
       */
      ~Transform3DTo2D() {}


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
      Transform3DTo2DFunctor<Type>
      getFunctor() const;
    
    
      /** 
       * Change the Transform3DTo2D value by explicitly setting element values
       * as if setting the elements of a 3x4 transformation matrix:
       *    [[a00, a01, a02, a03],
       *     [a10, a11, a12, a13],
       *     [a20, a21, a22, a23]]
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
       */
      void
      setTransform(
        Type const& a00, Type const& a01, Type const& a02, Type const& a03,
        Type const& a10, Type const& a11, Type const& a12, Type const& a13,
        Type const& a20, Type const& a21, Type const& a22, Type const& a23);


      /** 
       * This member function returns one element from the matrix
       * representation of the coordinate transform by value.  For
       * example, calling myTransform3DTo2D.value<2, 3>() will return
       * the element from the 2nd row, 3rd column of the matrix
       * representation of the coordinate transformation.
       * 
       * @return The value of the requested element.
       */
      template <size_t row, size_t column>
      Type const&
      value() const;


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
       * Applies the coordinate transformation to a Vector3D point and
       * returns the result. 
       * 
       * @param vector0 The point to be transformed.
       * @return The result of transforming vector0.
       */
      Vector2D
      operator*(const Vector3D& vector0) const;

      /** 
       * The assignment operator simply duplicates its argument.
       * 
       * @param source This is the Transform3D instance to be copied.
       * @return A reference to *this.
       */
      Transform3DTo2D<Type>&
      operator=(const Transform3DTo2D<Type>& source);

    private:

      Type m_00, m_01, m_02, m_03;
      Type m_10, m_11, m_12, m_13;
      Type m_20, m_21, m_22, m_23;

    }; // class Transform3DTo2D


    /* ============== Helper classes ============== */


    /**
     ** This helper class works with Transform3DTo2D::getFunctor()
     **/
    template <class Type>
    class Transform3DTo2DFunctor
      : public std::unary_function<Vector3D, Vector2D> {

    public:

      /** 
       * The constructor deep-copies its argument.
       * 
       *  @param transform This is the transform instance to be copied.
       */
      Transform3DTo2DFunctor(Transform3DTo2D<Type> const& transform)
        : m_transform(transform) {}

    
      /** 
       * The application operator transforms its argument.
       * 
       *  @param vec This is the vector to be transformed.
       * 
       *  @return The transformed vector.
       */
      inline Vector2D
      operator()(const Vector3D& vec) const {return m_transform * vec;}
    
    private:
      Transform3DTo2D<Type> m_transform;
    };
  

    /* ================ Non member functions below ================ */

    /** 
     * This operator composes a Transform3DTo2D instance with a Transform3D
     * instance.  The resulting Transform3DTo2D instance the equation:
     *
     *  (transform0 * transform1) * v0 = transform0 * (transform1 * v0),
     *
     * where v0 is a Vector3D instance.
     * 
     * @param transform0 This is the Transform3DTo2D instance to be
     * composed.
     *
     * @param transform1 This is the Transform3D instance to be composed.
     *
     * @return A Transform3DTo2D instance which is equivalent to the 
     * composition of the two arguments.
     */
    template <class Type>
    Transform3DTo2D<Type>
    operator*(Transform3DTo2D<Type> const& transform0,
              Transform3D<Type> const& transform1);


    /** 
     * Outputs a text representation of a Transform3D instance to a
     * std::ostream.  The output format looks like this:
     *
     * Transform3DTo2D(1.0, 2.0, 3.0, 4.0,
     *                 5.0, 6.0, 7.0, 8.0,
     *                 9.0, 10.0, 11.0, 1.0)
     *
     * @param stream Reference to the the output stream.
     *
     * @param transform0 Const reference to the Transform3DTo2D instance 
     * to be output.
     *
     * @return Reference to the output stream after writing.
     */
    template <class Type>
    std::ostream&
    operator<<(std::ostream& stream, Transform3DTo2D<Type> const& transform0);

  
    /** 
     * Sets the value of a Transform3DTo2D instance from a std::istream.
     * The input format is as described for
     * operator<<(std::ostream&, const Transform3DTo2D&), above.
     * 
     * @param stream Reference to the the input stream.
     *
     * @param transform0 Reference to the Transform3DTo2D that will take
     * the input.
     *
     * @return Reference to the input stream.
     */
    template <class Type>
    std::istream&
    operator>>(std::istream& stream, Transform3DTo2D<Type>& transform0);
  
  } // namespace numeric

}; // namespace brick


#include <transform3DTo2D_impl.hh>

#endif /* #ifndef BRICK_NUMERIC_TRANSFORM3DTO2D_HH */
