/**
***************************************************************************
* @file transform3DTo2D.cpp
*
* Source file defining Transform3DTo2D class template.
*
* Copyright (C) 2001-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_TRANSFORM3DTO2D_IMPL_HH
#define BRICK_NUMERIC_TRANSFORM3DTO2D_IMPL_HH

// This file is included by transform3DTo2D.hh, and should not be directly
// included by user code, so no need to include transform3DTo2D.hh here.
// 
// #include <brick/numeric/transform3DTo2D.hh>

namespace brick {

  namespace numeric {
    
    // Default constructor
    template <class Type>
    Transform3DTo2D<Type>::
    Transform3DTo2D()
      : m_00(1.0), m_01(0.0), m_02(0.0), m_03(0.0),
        m_10(0.0), m_11(1.0), m_12(0.0), m_13(0.0),
        m_20(0.0), m_21(0.0), m_22(1.0), m_23(1.0)
    {
      // Empty.
    }


    // Build a Transform3DTo2D instance by explicitly setting element 
    // values as if setting the elements of a 3x4 transformation matrix.
    template <class Type>
    inline
    Transform3DTo2D<Type>::
    Transform3DTo2D(
      Type const& a00, Type const& a01, Type const& a02, Type const& a03,
      Type const& a10, Type const& a11, Type const& a12, Type const& a13,
      Type const& a20, Type const& a21, Type const& a22, Type const& a23)
      : m_00(a00), m_01(a01), m_02(a02), m_03(a03),
        m_10(a10), m_11(a11), m_12(a12), m_13(a13),
        m_20(a20), m_21(a21), m_22(a22), m_23(a23)
    {
      // Empty.
    }

    
    // Build a Transform3DTo2D from a homogeneous 3x4 matrix.
    template <class Type>
    Transform3DTo2D<Type>::
    Transform3DTo2D(const Array2D<Type>& source)
    {
      if((source.rows() != 3) || (source.columns() != 4)) {
        std::ostringstream message;
        message << "Can't create a Transform3DTo2D from a " << source.rows()
                << " x " << source.columns() << "Array2D<Type> instance.";
        BRICK_THROW(ValueException, "Transform3DTo2D::Transform3DTo2D()",
                  message.str().c_str());
      }
      m_00 = source(0); m_01 = source(1); m_02 = source(2); m_03 = source(3);
      m_10 = source(4); m_11 = source(5); m_12 = source(6); m_13 = source(7);
      m_20 = source(8); m_21 = source(9); m_22 = source(10); m_23 = source(11);
    }


    // The copy constructor simply duplicates its argument.
    template <class Type>
    inline
    Transform3DTo2D<Type>::
    Transform3DTo2D(Transform3DTo2D<Type> const& src)
      : m_00(src.m_00), m_01(src.m_01), m_02(src.m_02), m_03(src.m_03),
        m_10(src.m_10), m_11(src.m_11), m_12(src.m_12), m_13(src.m_13),
        m_20(src.m_20), m_21(src.m_21), m_22(src.m_22), m_23(src.m_23)
    {
      // Empty.
    }


    // This member function returns a functor which makes it easier to
    // transform arrays of points using algorithms such as
    // std::transform().
    template <class Type>
    Transform3DTo2DFunctor<Type>
    Transform3DTo2D<Type>::
    getFunctor() const
    {
      return Transform3DTo2DFunctor<Type>(*this);
    }    
  
  
    // Change the Transform3DTo2D value by explicitly setting element values
    // as if setting the elements of a 4x4 transformation matrix:
    //    [[a00, a01, a02, a03],
    //     [a10, a11, a12, a13],
    //     [a20, a21, a22, a23],
    //     [a30, a31, a32, a33]]
    template <class Type>
    void
    Transform3DTo2D<Type>::
    setTransform(
      Type const& a00, Type const& a01, Type const& a02, Type const& a03,
      Type const& a10, Type const& a11, Type const& a12, Type const& a13,
      Type const& a20, Type const& a21, Type const& a22, Type const& a23)
    {
      m_00 = a00; m_01 = a01; m_02 = a02; m_03 = a03;
      m_10 = a10; m_11 = a11; m_12 = a12; m_13 = a13;
      m_20 = a20; m_21 = a21; m_22 = a22; m_23 = a23;
    }


    // This operator returns one element from the matrix
    // representation of the coordinate transform by value.
    // The case statements should optimize away, since row and column
    // are both known at compile time.
    template <class Type>
    template <size_t ROW, size_t COLUMN>
    Type const&
    Transform3DTo2D<Type>::
    value() const
    {
      // // Avoid ugly duplication of code using ugly const_cast.
      // return const_cast<Transform3D*>(this)->value<ROW, COLUMN>();
      switch(ROW) {
      case 0:
        switch(COLUMN) {
        case 0: return m_00; break;
        case 1: return m_01; break;
        case 2: return m_02; break;
        case 3: return m_03; break;
        default: break;
        }
        break;
      case 1:
        switch(COLUMN) {
        case 0: return m_10; break;
        case 1: return m_11; break;
        case 2: return m_12; break;
        case 3: return m_13; break;
        default: break;
        }
        break;
      case 2:
        switch(COLUMN) {
        case 0: return m_20; break;
        case 1: return m_21; break;
        case 2: return m_22; break;
        case 3: return m_23; break;
        default: break;
        }
        break;
      default:
        break;
      }
      std::ostringstream message;
      message << "Indices (" << ROW << ", " << COLUMN << ") are out of bounds.";
      BRICK_THROW(IndexException, "Transform3DTo2D::value<size_t, size_t>()",
                message.str().c_str());
      return 0.0; // Dummy return to keep the compiler happy.
    }

    
    // This operator returns one element from the matrix
    // representation of the coordinate transform by value.
    template <class Type>
    Type const&
    Transform3DTo2D<Type>::
    operator()(size_t row, size_t column) const
    {
      // // Avoid ugly duplication of code using ugly const_cast.
      // return const_cast<Transform3DTo2D*>(this)->operator()(row, column);
      switch(row) {
      case 0:
        switch(column) {
        case 0: return m_00; break;
        case 1: return m_01; break;
        case 2: return m_02; break;
        case 3: return m_03; break;
        default: break;
        }
        break;
      case 1:
        switch(column) {
        case 0: return m_10; break;
        case 1: return m_11; break;
        case 2: return m_12; break;
        case 3: return m_13; break;
        default: break;
        }
        break;
      case 2:
        switch(column) {
        case 0: return m_20; break;
        case 1: return m_21; break;
        case 2: return m_22; break;
        case 3: return m_23; break;
        default: break;
        }
        break;
      default:
        break;
      }
      std::ostringstream message;
      message << "Indices (" << row << ", " << column << ") are out of bounds.";
      BRICK_THROW(IndexException, "Transform3DTo2D::operator()(size_t, size_t)",
                  message.str().c_str());
      return m_23; // Dummy return to keep the compiler happy.
    }

    
    // This operator takes a point and applies the coordinate
    // transform, returning the result.
    template <class Type>
    Vector2D
    Transform3DTo2D<Type>::
    operator*(const Vector3D& vector0) const
    {
      return Vector2D(
        m_00 * vector0.x() + m_01 * vector0.y() + m_02 * vector0.z() + m_03,
        m_10 * vector0.x() + m_11 * vector0.y() + m_12 * vector0.z() + m_13,
        m_20 * vector0.x() + m_21 * vector0.y() + m_22 * vector0.z() + m_23);
    }

    
    // The assignment operator simply duplicates its argument.
    template <class Type>
    Transform3DTo2D<Type>&
    Transform3DTo2D<Type>::
    operator=(const Transform3DTo2D<Type>& source)
    {
      if(&source != this) {
        m_00 = source.m_00; m_01 = source.m_01;
        m_02 = source.m_02; m_03 = source.m_03;
        m_10 = source.m_10; m_11 = source.m_11;
        m_12 = source.m_12; m_13 = source.m_13;
        m_20 = source.m_20; m_21 = source.m_21;
        m_22 = source.m_22; m_23 = source.m_23;
      }
      return *this;
    }

    
    /* ================ Non member functions below ================ */


    // This operator composes a Transform3DTo2D instance with a Transform3D
    // instance.
    template <class Type>
    Transform3DTo2D<Type>
    operator*(Transform3DTo2D<Type> const& transform0,
              Transform3D const& transform1)
    {
      // We'd rather use the templated getValue() member here, as it's
      // faster than operator()(), but for some reason the compiler is
      // choking on it, so we default to what works.  We repair most
      // of the damage by specializing for double and float below.
#if 0
      Type a00 = (transform0.value<0, 0>() * transform1.value<0, 0>()
                    + transform0.value<0, 1>() * transform1.value<1, 0>()
                    + transform0.value<0, 2>() * transform1.value<2, 0>()
                    + transform0.value<0, 3>() * transform1.value<3, 0>());
      Type a01 = (transform0.value<0, 0>() * transform1.value<0, 1>()
                    + transform0.value<0, 1>() * transform1.value<1, 1>()
                    + transform0.value<0, 2>() * transform1.value<2, 1>()
                    + transform0.value<0, 3>() * transform1.value<3, 1>());
      Type a02 = (transform0.value<0, 0>() * transform1.value<0, 2>()
                    + transform0.value<0, 1>() * transform1.value<1, 2>()
                    + transform0.value<0, 2>() * transform1.value<2, 2>()
                    + transform0.value<0, 3>() * transform1.value<3, 2>());
      Type a03 = (transform0.value<0, 0>() * transform1.value<0, 3>()
                    + transform0.value<0, 1>() * transform1.value<1, 3>()
                    + transform0.value<0, 2>() * transform1.value<2, 3>()
                    + transform0.value<0, 3>() * transform1.value<3, 3>());
      Type a10 = (transform0.value<1, 0>() * transform1.value<0, 0>()
                    + transform0.value<1, 1>() * transform1.value<1, 0>()
                    + transform0.value<1, 2>() * transform1.value<2, 0>()
                    + transform0.value<1, 3>() * transform1.value<3, 0>());
      Type a11 = (transform0.value<1, 0>() * transform1.value<0, 1>()
                    + transform0.value<1, 1>() * transform1.value<1, 1>()
                    + transform0.value<1, 2>() * transform1.value<2, 1>()
                    + transform0.value<1, 3>() * transform1.value<3, 1>());
      Type a12 = (transform0.value<1, 0>() * transform1.value<0, 2>()
                    + transform0.value<1, 1>() * transform1.value<1, 2>()
                    + transform0.value<1, 2>() * transform1.value<2, 2>()
                    + transform0.value<1, 3>() * transform1.value<3, 2>());
      Type a13 = (transform0.value<1, 0>() * transform1.value<0, 3>()
                    + transform0.value<1, 1>() * transform1.value<1, 3>()
                    + transform0.value<1, 2>() * transform1.value<2, 3>()
                    + transform0.value<1, 3>() * transform1.value<3, 3>());
      Type a20 = (transform0.value<2, 0>() * transform1.value<0, 0>()
                    + transform0.value<2, 1>() * transform1.value<1, 0>()
                    + transform0.value<2, 2>() * transform1.value<2, 0>()
                    + transform0.value<2, 3>() * transform1.value<3, 0>());
      Type a21 = (transform0.value<2, 0>() * transform1.value<0, 1>()
                    + transform0.value<2, 1>() * transform1.value<1, 1>()
                    + transform0.value<2, 2>() * transform1.value<2, 1>()
                    + transform0.value<2, 3>() * transform1.value<3, 1>());
      Type a22 = (transform0.value<2, 0>() * transform1.value<0, 2>()
                    + transform0.value<2, 1>() * transform1.value<1, 2>()
                    + transform0.value<2, 2>() * transform1.value<2, 2>()
                    + transform0.value<2, 3>() * transform1.value<3, 2>());
      Type a23 = (transform0.value<2, 0>() * transform1.value<0, 3>()
                    + transform0.value<2, 1>() * transform1.value<1, 3>()
                    + transform0.value<2, 2>() * transform1.value<2, 3>()
                    + transform0.value<2, 3>() * transform1.value<3, 3>());
#else /* #if 0 */
      Type a00 = (transform0(0, 0) * transform1(0, 0)
                  + transform0(0, 1) * transform1(1, 0)
                  + transform0(0, 2) * transform1(2, 0)
                  + transform0(0, 3) * transform1(3, 0));
      Type a01 = (transform0(0, 0) * transform1(0, 1)
                  + transform0(0, 1) * transform1(1, 1)
                  + transform0(0, 2) * transform1(2, 1)
                  + transform0(0, 3) * transform1(3, 1));
      Type a02 = (transform0(0, 0) * transform1(0, 2)
                  + transform0(0, 1) * transform1(1, 2)
                  + transform0(0, 2) * transform1(2, 2)
                  + transform0(0, 3) * transform1(3, 2));
      Type a03 = (transform0(0, 0) * transform1(0, 3)
                  + transform0(0, 1) * transform1(1, 3)
                  + transform0(0, 2) * transform1(2, 3)
                  + transform0(0, 3) * transform1(3, 3));
      Type a10 = (transform0(1, 0) * transform1(0, 0)
                  + transform0(1, 1) * transform1(1, 0)
                  + transform0(1, 2) * transform1(2, 0)
                  + transform0(1, 3) * transform1(3, 0));
      Type a11 = (transform0(1, 0) * transform1(0, 1)
                  + transform0(1, 1) * transform1(1, 1)
                  + transform0(1, 2) * transform1(2, 1)
                  + transform0(1, 3) * transform1(3, 1));
      Type a12 = (transform0(1, 0) * transform1(0, 2)
                  + transform0(1, 1) * transform1(1, 2)
                  + transform0(1, 2) * transform1(2, 2)
                  + transform0(1, 3) * transform1(3, 2));
      Type a13 = (transform0(1, 0) * transform1(0, 3)
                  + transform0(1, 1) * transform1(1, 3)
                  + transform0(1, 2) * transform1(2, 3)
                  + transform0(1, 3) * transform1(3, 3));
      Type a20 = (transform0(2, 0) * transform1(0, 0)
                  + transform0(2, 1) * transform1(1, 0)
                  + transform0(2, 2) * transform1(2, 0)
                  + transform0(2, 3) * transform1(3, 0));
      Type a21 = (transform0(2, 0) * transform1(0, 1)
                  + transform0(2, 1) * transform1(1, 1)
                  + transform0(2, 2) * transform1(2, 1)
                  + transform0(2, 3) * transform1(3, 1));
      Type a22 = (transform0(2, 0) * transform1(0, 2)
                  + transform0(2, 1) * transform1(1, 2)
                  + transform0(2, 2) * transform1(2, 2)
                  + transform0(2, 3) * transform1(3, 2));
      Type a23 = (transform0(2, 0) * transform1(0, 3)
                  + transform0(2, 1) * transform1(1, 3)
                  + transform0(2, 2) * transform1(2, 3)
                  + transform0(2, 3) * transform1(3, 3));
#endif /* #if 0 */

      return Transform3DTo2D<Type>(a00, a01, a02, a03,
                             a10, a11, a12, a13,
                             a20, a21, a22, a23);
    }


    template <>
    Transform3DTo2D<double>
    operator*(Transform3DTo2D<double> const& transform0,
              Transform3D const& transform1)
    {
      double a00 = (transform0.value<0, 0>() * transform1.value<0, 0>()
                    + transform0.value<0, 1>() * transform1.value<1, 0>()
                    + transform0.value<0, 2>() * transform1.value<2, 0>()
                    + transform0.value<0, 3>() * transform1.value<3, 0>());
      double a01 = (transform0.value<0, 0>() * transform1.value<0, 1>()
                    + transform0.value<0, 1>() * transform1.value<1, 1>()
                    + transform0.value<0, 2>() * transform1.value<2, 1>()
                    + transform0.value<0, 3>() * transform1.value<3, 1>());
      double a02 = (transform0.value<0, 0>() * transform1.value<0, 2>()
                    + transform0.value<0, 1>() * transform1.value<1, 2>()
                    + transform0.value<0, 2>() * transform1.value<2, 2>()
                    + transform0.value<0, 3>() * transform1.value<3, 2>());
      double a03 = (transform0.value<0, 0>() * transform1.value<0, 3>()
                    + transform0.value<0, 1>() * transform1.value<1, 3>()
                    + transform0.value<0, 2>() * transform1.value<2, 3>()
                    + transform0.value<0, 3>() * transform1.value<3, 3>());
      double a10 = (transform0.value<1, 0>() * transform1.value<0, 0>()
                    + transform0.value<1, 1>() * transform1.value<1, 0>()
                    + transform0.value<1, 2>() * transform1.value<2, 0>()
                    + transform0.value<1, 3>() * transform1.value<3, 0>());
      double a11 = (transform0.value<1, 0>() * transform1.value<0, 1>()
                    + transform0.value<1, 1>() * transform1.value<1, 1>()
                    + transform0.value<1, 2>() * transform1.value<2, 1>()
                    + transform0.value<1, 3>() * transform1.value<3, 1>());
      double a12 = (transform0.value<1, 0>() * transform1.value<0, 2>()
                    + transform0.value<1, 1>() * transform1.value<1, 2>()
                    + transform0.value<1, 2>() * transform1.value<2, 2>()
                    + transform0.value<1, 3>() * transform1.value<3, 2>());
      double a13 = (transform0.value<1, 0>() * transform1.value<0, 3>()
                    + transform0.value<1, 1>() * transform1.value<1, 3>()
                    + transform0.value<1, 2>() * transform1.value<2, 3>()
                    + transform0.value<1, 3>() * transform1.value<3, 3>());
      double a20 = (transform0.value<2, 0>() * transform1.value<0, 0>()
                    + transform0.value<2, 1>() * transform1.value<1, 0>()
                    + transform0.value<2, 2>() * transform1.value<2, 0>()
                    + transform0.value<2, 3>() * transform1.value<3, 0>());
      double a21 = (transform0.value<2, 0>() * transform1.value<0, 1>()
                    + transform0.value<2, 1>() * transform1.value<1, 1>()
                    + transform0.value<2, 2>() * transform1.value<2, 1>()
                    + transform0.value<2, 3>() * transform1.value<3, 1>());
      double a22 = (transform0.value<2, 0>() * transform1.value<0, 2>()
                    + transform0.value<2, 1>() * transform1.value<1, 2>()
                    + transform0.value<2, 2>() * transform1.value<2, 2>()
                    + transform0.value<2, 3>() * transform1.value<3, 2>());
      double a23 = (transform0.value<2, 0>() * transform1.value<0, 3>()
                    + transform0.value<2, 1>() * transform1.value<1, 3>()
                    + transform0.value<2, 2>() * transform1.value<2, 3>()
                    + transform0.value<2, 3>() * transform1.value<3, 3>());
      return Transform3DTo2D<double>(a00, a01, a02, a03,
                             a10, a11, a12, a13,
                             a20, a21, a22, a23);
    }


    template <>
    Transform3DTo2D<float>
    operator*(Transform3DTo2D<float> const& transform0,
              Transform3D const& transform1)
    {
      float a00 = (transform0.value<0, 0>() * transform1.value<0, 0>()
                    + transform0.value<0, 1>() * transform1.value<1, 0>()
                    + transform0.value<0, 2>() * transform1.value<2, 0>()
                    + transform0.value<0, 3>() * transform1.value<3, 0>());
      float a01 = (transform0.value<0, 0>() * transform1.value<0, 1>()
                    + transform0.value<0, 1>() * transform1.value<1, 1>()
                    + transform0.value<0, 2>() * transform1.value<2, 1>()
                    + transform0.value<0, 3>() * transform1.value<3, 1>());
      float a02 = (transform0.value<0, 0>() * transform1.value<0, 2>()
                    + transform0.value<0, 1>() * transform1.value<1, 2>()
                    + transform0.value<0, 2>() * transform1.value<2, 2>()
                    + transform0.value<0, 3>() * transform1.value<3, 2>());
      float a03 = (transform0.value<0, 0>() * transform1.value<0, 3>()
                    + transform0.value<0, 1>() * transform1.value<1, 3>()
                    + transform0.value<0, 2>() * transform1.value<2, 3>()
                    + transform0.value<0, 3>() * transform1.value<3, 3>());
      float a10 = (transform0.value<1, 0>() * transform1.value<0, 0>()
                    + transform0.value<1, 1>() * transform1.value<1, 0>()
                    + transform0.value<1, 2>() * transform1.value<2, 0>()
                    + transform0.value<1, 3>() * transform1.value<3, 0>());
      float a11 = (transform0.value<1, 0>() * transform1.value<0, 1>()
                    + transform0.value<1, 1>() * transform1.value<1, 1>()
                    + transform0.value<1, 2>() * transform1.value<2, 1>()
                    + transform0.value<1, 3>() * transform1.value<3, 1>());
      float a12 = (transform0.value<1, 0>() * transform1.value<0, 2>()
                    + transform0.value<1, 1>() * transform1.value<1, 2>()
                    + transform0.value<1, 2>() * transform1.value<2, 2>()
                    + transform0.value<1, 3>() * transform1.value<3, 2>());
      float a13 = (transform0.value<1, 0>() * transform1.value<0, 3>()
                    + transform0.value<1, 1>() * transform1.value<1, 3>()
                    + transform0.value<1, 2>() * transform1.value<2, 3>()
                    + transform0.value<1, 3>() * transform1.value<3, 3>());
      float a20 = (transform0.value<2, 0>() * transform1.value<0, 0>()
                    + transform0.value<2, 1>() * transform1.value<1, 0>()
                    + transform0.value<2, 2>() * transform1.value<2, 0>()
                    + transform0.value<2, 3>() * transform1.value<3, 0>());
      float a21 = (transform0.value<2, 0>() * transform1.value<0, 1>()
                    + transform0.value<2, 1>() * transform1.value<1, 1>()
                    + transform0.value<2, 2>() * transform1.value<2, 1>()
                    + transform0.value<2, 3>() * transform1.value<3, 1>());
      float a22 = (transform0.value<2, 0>() * transform1.value<0, 2>()
                    + transform0.value<2, 1>() * transform1.value<1, 2>()
                    + transform0.value<2, 2>() * transform1.value<2, 2>()
                    + transform0.value<2, 3>() * transform1.value<3, 2>());
      float a23 = (transform0.value<2, 0>() * transform1.value<0, 3>()
                    + transform0.value<2, 1>() * transform1.value<1, 3>()
                    + transform0.value<2, 2>() * transform1.value<2, 3>()
                    + transform0.value<2, 3>() * transform1.value<3, 3>());
      return Transform3DTo2D<float>(a00, a01, a02, a03,
                             a10, a11, a12, a13,
                             a20, a21, a22, a23);
    }
    

    template <class Type>
    std::ostream&
    operator<<(std::ostream& stream, const Transform3DTo2D<Type>& transform0)
    {
      stream << "Transform3DTo2D("
             << transform0.value<0, 0>() << ", "
             << transform0.value<0, 1>() << ", "
             << transform0.value<0, 2>() << ", "
             << transform0.value<0, 3>() << ",\n"
             << transform0.value<1, 0>() << ", "
             << transform0.value<1, 1>() << ", "
             << transform0.value<1, 2>() << ", "
             << transform0.value<1, 3>() << ",\n"
             << transform0.value<2, 0>() << ", "
             << transform0.value<2, 1>() << ", "
             << transform0.value<2, 2>() << ", "
             << transform0.value<2, 3>() << ")";
      return stream;
    }

    
    template <class Type>
    std::istream&
    operator>>(std::istream& stream, Transform3DTo2D<Type>& transform0)
    {
      // If stream is in a bad state, we can't read from it.
      if (!stream){
        return stream;
      }
    
      // It's a lot easier to use a try block than to be constantly
      // testing whether the IO has succeeded, so we tell stream to
      // complain if anything goes wrong.
      std::ios_base::iostate oldExceptionState = stream.exceptions();
      stream.exceptions(
        std::ios_base::badbit | std::ios_base::failbit | std::ios_base::eofbit);

      // Now on with the show.
      try{
        // Construct an InputStream instance so we can use our
        // convenience functions.
        InputStream inputStream(stream);

        // Advance to the next relevant character.
        inputStream.skipWhiteSpace();
      
        // Read the header.
        inputStream.expect("Transform3DTo2D(");

        // Read all the data except the last element.
        std::vector<Type> inputValues(12);
        for(size_t index = 0; index < (inputValues.size() - 1); ++index) {
          // Read the value.
          inputStream >> inputValues[index];

          // Read punctuation before the next value.
          inputStream.expect(",");
        }

        // Read the final value.
        inputStream >> inputValues[inputValues.size() - 1];

        // Read the closing parenthesis.
        inputStream.expect(")");

        // And update the transform.
        transform0.setTransform(inputValues[0], inputValues[1],
                                inputValues[2], inputValues[3], 
                                inputValues[4], inputValues[5], 
                                inputValues[6], inputValues[7], 
                                inputValues[8], inputValues[9], 
                                inputValues[10], inputValues[11]);
      } catch(std::ios_base::failure) {
        // Empty
      }
      stream.exceptions(oldExceptionState);
      return stream;
    }

  } // namespace numeric

} // namespace brick

#endif /* #ifndef BRICK_NUMERIC_TRANSFORM3DTO2D_IMPL_HH */