/**
***************************************************************************
* @file brick/numeric/transform3D_impl.hh
*
* Header file providing implementations for inline and template
* functions declared in transform3D.hh.
*
* Copyright (C) 2001-2013 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_TRANSFORM3D_IMPL_HH
#define BRICK_NUMERIC_TRANSFORM3D_IMPL_HH

// This file is included by transform3D.hh, and should not be directly
// included by user code, so no need to include transform3D.hh here.
//
// #include <brick/numeric/transform3D.hh>

namespace brick {

  namespace numeric {

    // Default constructor.  Initializes to identity.
    template <class Type>
    inline
    Transform3D<Type>::
    Transform3D()
      : m_00(1.0), m_01(0.0), m_02(0.0), m_03(0.0),
        m_10(0.0), m_11(1.0), m_12(0.0), m_13(0.0),
        m_20(0.0), m_21(0.0), m_22(1.0), m_23(0.0),
        m_30(0.0), m_31(0.0), m_32(0.0), m_33(1.0)
    {
      // Empty.
    }


    // Build a Transform3D instance by explicitly setting element values
    // as if setting the elements of a 4x4 transformation matrix.
    template <class Type>
    inline
    Transform3D<Type>::
    Transform3D(
      Type const& a00, Type const& a01, Type const& a02, Type const& a03,
      Type const& a10, Type const& a11, Type const& a12, Type const& a13,
      Type const& a20, Type const& a21, Type const& a22, Type const& a23,
      Type const& a30, Type const& a31, Type const& a32, Type const& a33,
      bool doNormalize)
      : m_00(a00), m_01(a01), m_02(a02), m_03(a03),
        m_10(a10), m_11(a11), m_12(a12), m_13(a13),
        m_20(a20), m_21(a21), m_22(a22), m_23(a23),
        m_30(a30), m_31(a31), m_32(a32), m_33(a33)
    {
      if(doNormalize) {
        this->normalize();
      }
    }


    // Build a Transform3D instance from a sequence that specifies
    // element values in row major order.
    template <class Type>
    Transform3D<Type>::
    Transform3D(std::initializer_list<Type> sequence,
                bool doNormalize)
      : m_00(Type(1)), m_01(Type(0)), m_02(Type(0)), m_03(Type(0)),
        m_10(Type(0)), m_11(Type(1)), m_12(Type(0)), m_13(Type(0)),
        m_20(Type(0)), m_21(Type(0)), m_22(Type(1)), m_23(Type(0)),
        m_30(Type(0)), m_31(Type(0)), m_32(Type(0)), m_33(Type(1))
    {
      if(sequence.size() < 15) {
        BRICK_THROW(common::ValueException, "Transform3D::Transform3D()",
                    "Initialization list is not long enough.  Expected at "
                    "least 15 elements.");
      }
      auto iter = sequence.begin();
      m_00 = *iter; ++iter;
      m_01 = *iter; ++iter;
      m_02 = *iter; ++iter;
      m_03 = *iter; ++iter;
      m_10 = *iter; ++iter;
      m_11 = *iter; ++iter;
      m_12 = *iter; ++iter;
      m_13 = *iter; ++iter;
      m_20 = *iter; ++iter;
      m_21 = *iter; ++iter;
      m_22 = *iter; ++iter;
      m_23 = *iter; ++iter;
      m_30 = *iter; ++iter;
      m_31 = *iter; ++iter;
      m_32 = *iter; ++iter;
      if(iter != sequence.end()) {
        m_33 = *iter;
      }
      if(doNormalize) {
        this->normalize();
      }
    }


    // Build a Transform3D from a homogeneous 4x4 matrix.
    template <class Type>
    Transform3D<Type>::
    Transform3D(const Array2D<Type>& source)
      : m_00(1.0), m_01(0.0), m_02(0.0), m_03(0.0),
        m_10(0.0), m_11(1.0), m_12(0.0), m_13(0.0),
        m_20(0.0), m_21(0.0), m_22(1.0), m_23(0.0),
        m_30(0.0), m_31(0.0), m_32(0.0), m_33(1.0)
    {
      if((source.rows() != 4) || (source.columns() != 4)) {
        std::ostringstream message;
        message << "Can't create a Transform3D from a " << source.rows()
                << " x " << source.columns() << "Array2D instance.";
        BRICK_THROW(common::ValueException, "Transform3D::Transform3D()",
                  message.str().c_str());
      }
      m_00 = source(0); m_01 = source(1); m_02 = source(2); m_03 = source(3);
      m_10 = source(4); m_11 = source(5); m_12 = source(6); m_13 = source(7);
      m_20 = source(8); m_21 = source(9); m_22 = source(10); m_23 = source(11);
      m_30 = source(12); m_31 = source(13); m_32 = source(14);
      m_33 = source(15);
    }


    // The copy constructor simply duplicates its argument.
    template <class Type>
    inline
    Transform3D<Type>::
    Transform3D(Transform3D<Type> const& src)
      : m_00(src.m_00), m_01(src.m_01), m_02(src.m_02), m_03(src.m_03),
        m_10(src.m_10), m_11(src.m_11), m_12(src.m_12), m_13(src.m_13),
        m_20(src.m_20), m_21(src.m_21), m_22(src.m_22), m_23(src.m_23),
        m_30(src.m_30), m_31(src.m_31), m_32(src.m_32), m_33(src.m_33)
    {
      // Empty.
    }


    // This member function returns a functor which makes it easier to
    // transform arrays of points using algorithms such as
    // std::transform().
    template <class Type>
    Transform3DFunctor<Type>
    Transform3D<Type>::
    getFunctor() const
    {
      return Transform3DFunctor<Type>(*this);
    }


    // This member function returns one element from the matrix
    // representation of the coordinate transform by value.
    template <class Type>
    template <size_t row, size_t column>
    inline
    Type const&
    Transform3D<Type>::
    getValue() const
    {
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
      case 3:
        switch(column) {
        case 0: return m_30; break;
        case 1: return m_31; break;
        case 2: return m_32; break;
        case 3: return m_33; break;
        default: break;
        }
        break;
      default:
        break;
      }
      std::ostringstream message;
      message << "Indices (" << row << ", " << column << ") are out of bounds.";
      BRICK_THROW(common::IndexException,
                  "Transform3D::value<size_t, size_t>()",
                  message.str().c_str());
      return m_33; // Dummy return to keep the compiler happy.
    }


    // This member function returns the inverse of *this.
    template <class Type>
    Transform3D<Type>
    Transform3D<Type>::
    invert() const
    {
      // We use the cofactor method for now, since it's easier to code
      // than Gauss-Jordan elimination.  We suspect that it's less
      // efficient, however.

      // Notation for determinant values is detRRRCCC, where the
      // Rs indicate the involved rows, from top to bottom, and the Cs
      // indicate the involved columns, from left to right.

      Type det1201 = m_10 * m_21 - m_11 * m_20;
      Type det1202 = m_10 * m_22 - m_12 * m_20;
      Type det1203 = m_10 * m_23 - m_13 * m_20;
      Type det1212 = m_11 * m_22 - m_12 * m_21;
      Type det1213 = m_11 * m_23 - m_13 * m_21;
      Type det1223 = m_12 * m_23 - m_13 * m_22;

      Type det1301 = m_10 * m_31 - m_11 * m_30;
      Type det1302 = m_10 * m_32 - m_12 * m_30;
      Type det1303 = m_10 * m_33 - m_13 * m_30;
      Type det1312 = m_11 * m_32 - m_12 * m_31;
      Type det1313 = m_11 * m_33 - m_13 * m_31;
      Type det1323 = m_12 * m_33 - m_13 * m_32;

      Type det2301 = m_20 * m_31 - m_21 * m_30;
      Type det2302 = m_20 * m_32 - m_22 * m_30;
      Type det2303 = m_20 * m_33 - m_23 * m_30;
      Type det2312 = m_21 * m_32 - m_22 * m_31;
      Type det2313 = m_21 * m_33 - m_23 * m_31;
      Type det2323 = m_22 * m_33 - m_23 * m_32;

      Type det012012 = (m_00 * det1212 - m_01 * det1202 + m_02 * det1201);
      Type det012013 = (m_00 * det1213 - m_01 * det1203 + m_03 * det1201);
      Type det012023 = (m_00 * det1223 - m_02 * det1203 + m_03 * det1202);
      Type det012123 = (m_01 * det1223 - m_02 * det1213 + m_03 * det1212);

      Type det013012 = (m_00 * det1312 - m_01 * det1302 + m_02 * det1301);
      Type det013013 = (m_00 * det1313 - m_01 * det1303 + m_03 * det1301);
      Type det013023 = (m_00 * det1323 - m_02 * det1303 + m_03 * det1302);
      Type det013123 = (m_01 * det1323 - m_02 * det1313 + m_03 * det1312);

      Type det023012 = (m_00 * det2312 - m_01 * det2302 + m_02 * det2301);
      Type det023013 = (m_00 * det2313 - m_01 * det2303 + m_03 * det2301);
      Type det023023 = (m_00 * det2323 - m_02 * det2303 + m_03 * det2302);
      Type det023123 = (m_01 * det2323 - m_02 * det2313 + m_03 * det2312);

      Type det123012 = (m_10 * det2312 - m_11 * det2302 + m_12 * det2301);
      Type det123013 = (m_10 * det2313 - m_11 * det2303 + m_13 * det2301);
      Type det123023 = (m_10 * det2323 - m_12 * det2303 + m_13 * det2302);
      Type det123123 = (m_11 * det2323 - m_12 * det2313 + m_13 * det2312);

      // Type det01230123 = (
      //   m_00 * det123123 - m_01 * det123023
      //   + m_02 * det123013 - m_03 * det123012
      //   - m_10 * det023123 + m_11 * det023023
      //   - m_12 * det023013 + m_13 * det023012
      //   + m_20 * det013123 - m_21 * det013023
      //   + m_22 * det013013 - m_23 * det013012
      //   - m_30 * det012123 + m_31 * det012023
      //   - m_32 * det012013 + det012012);
      Type det01230123 = (
        m_00 * det123123 - m_01 * det123023
        + m_02 * det123013 - m_03 * det123012);

      // Note that in general, roundoff error will make us pass these
      // tests, even for singular matrices.
      if(det01230123 == Type(0.0)) {
        BRICK_THROW(common::ValueException, "Transform3D::invert()",
                  "Transform is not invertible.");
      }

      return Transform3D<Type>(
        det123123 / det01230123, -det023123 / det01230123,
        det013123 / det01230123, -det012123 / det01230123,
        -det123023 / det01230123, det023023 / det01230123,
        -det013023 / det01230123, det012023 / det01230123,
        det123013 / det01230123, -det023013 / det01230123,
        det013013 / det01230123, -det012013 / det01230123,
        -det123012 / det01230123, det023012 / det01230123,
        -det013012 / det01230123, det012012 / det01230123);
    }


    // This operator takes a point and applies only the rotational
    // part of the coordinate transform, returning the result.
    template <class Type>
    Vector3D<Type>
    Transform3D<Type>::
    rotate(const Vector3D<Type>& vector0) const
    {
      return Vector3D<Type>(
        m_00 * vector0.x() + m_01 * vector0.y() + m_02 * vector0.z(),
        m_10 * vector0.x() + m_11 * vector0.y() + m_12 * vector0.z(),
        m_20 * vector0.x() + m_21 * vector0.y() + m_22 * vector0.z());
    }


    // Change the Transform3D value by explicitly setting element values
    // as if setting the elements of a 4x4 transformation matrix:
    //    [[a00, a01, a02, a03],
    //     [a10, a11, a12, a13],
    //     [a20, a21, a22, a23],
    //     [a30, a31, a32, a33]]
    template <class Type>
    void
    Transform3D<Type>::
    setTransform(
      Type const& a00, Type const& a01, Type const& a02, Type const& a03,
      Type const& a10, Type const& a11, Type const& a12, Type const& a13,
      Type const& a20, Type const& a21, Type const& a22, Type const& a23,
      Type const& a30, Type const& a31, Type const& a32, Type const& a33,
      bool doNormalize)
    {
      m_00 = a00; m_01 = a01; m_02 = a02; m_03 = a03;
      m_10 = a10; m_11 = a11; m_12 = a12; m_13 = a13;
      m_20 = a20; m_21 = a21; m_22 = a22; m_23 = a23;
      m_30 = a30; m_31 = a31; m_32 = a32; m_33 = a33;
      if(doNormalize) {
        this->normalize();
      }
    }


    // This member function sets one element from the matrix
    // representation of the coordinate transform to the specified
    // value.
    template <class Type>
    void
    Transform3D<Type>::
    setValue(size_t row, size_t column, Type const& inValue)
    {
      switch(row) {
      case 0:
        switch(column) {
        case 0: m_00 = inValue; return; break;
        case 1: m_01 = inValue; return; break;
        case 2: m_02 = inValue; return; break;
        case 3: m_03 = inValue; return; break;
        default: break;
        }
        break;
      case 1:
        switch(column) {
        case 0: m_10 = inValue; return; break;
        case 1: m_11 = inValue; return; break;
        case 2: m_12 = inValue; return; break;
        case 3: m_13 = inValue; return; break;
        default: break;
        }
        break;
      case 2:
        switch(column) {
        case 0: m_20 = inValue; return; break;
        case 1: m_21 = inValue; return; break;
        case 2: m_22 = inValue; return; break;
        case 3: m_23 = inValue; return; break;
        default: break;
        }
        break;
      case 3:
        switch(column) {
        case 0: m_30 = inValue; return; break;
        case 1: m_31 = inValue; return; break;
        case 2: m_32 = inValue; return; break;
        case 3: m_33 = inValue; return; break;
        default: break;
        }
        break;
      default:
        break;
      }
      std::ostringstream message;
      message << "Indices (" << row << ", " << column << ") are out of bounds.";
      BRICK_THROW(common::IndexException,
                  "Transform3D::operator()(size_t, size_t)",
                  message.str().c_str());
    }


    // This operator sets one element of the matrix representation of
    // the coordinate transform.  The case statements should optimize
    // away, since row and column are both known at compile time.
    template <class Type>
    template <size_t row, size_t column>
    void
    Transform3D<Type>::
    setValue(Type const& inValue)
    {
      switch(row) {
      case 0:
        switch(column) {
        case 0: m_00 = inValue; return; break;
        case 1: m_01 = inValue; return; break;
        case 2: m_02 = inValue; return; break;
        case 3: m_03 = inValue; return; break;
        default: break;
        }
        break;
      case 1:
        switch(column) {
        case 0: m_10 = inValue; return; break;
        case 1: m_11 = inValue; return; break;
        case 2: m_12 = inValue; return; break;
        case 3: m_13 = inValue; return; break;
        default: break;
        }
        break;
      case 2:
        switch(column) {
        case 0: m_20 = inValue; return; break;
        case 1: m_21 = inValue; return; break;
        case 2: m_22 = inValue; return; break;
        case 3: m_23 = inValue; return; break;
        default: break;
        }
        break;
      case 3:
        switch(column) {
        case 0: m_30 = inValue; return; break;
        case 1: m_31 = inValue; return; break;
        case 2: m_32 = inValue; return; break;
        case 3: m_33 = inValue; return; break;
        default: break;
        }
        break;
      default:
        break;
      }
      std::ostringstream message;
      message << "Indices (" << row << ", " << column << ") are out of bounds.";
      BRICK_THROW(common::IndexException,
                  "Transform3D::setValue<size_t, size_t>(Type const&)",
                  message.str().c_str());
    }


    // This member function is an alias for member function setTransform().
    template <class Type>
    void
    Transform3D<Type>::
    setValue(
        Type const& a00, Type const& a01, Type const& a02, Type const& a03,
        Type const& a10, Type const& a11, Type const& a12, Type const& a13,
        Type const& a20, Type const& a21, Type const& a22, Type const& a23,
        Type const& a30, Type const& a31, Type const& a32, Type const& a33,
        bool doNormalize)
    {
      this->setTransform(a00, a01, a02, a03, a10, a11, a12, a13,
                         a20, a21, a22, a23, a30, a31, a32, a33, doNormalize);
    }


    // This operator returns one element from the matrix
    // representation of the coordinate transform by value.
    template <class Type>
    Type const&
    Transform3D<Type>::
    operator()(size_t row, size_t column) const
    {
      // // Avoid ugly duplication of code using ugly const_cast.
      // return const_cast<Transform3D*>(this)->operator()(row, column);
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
      case 3:
        switch(column) {
        case 0: return m_30; break;
        case 1: return m_31; break;
        case 2: return m_32; break;
        case 3: return m_33; break;
        default: break;
        }
        break;
      default:
        break;
      }
      std::ostringstream message;
      message << "Indices (" << row << ", " << column << ") are out of bounds.";
      BRICK_THROW(common::IndexException,
                  "Transform3D::operator()(size_t, size_t)",
                  message.str().c_str());
      return m_33; // Dummy return to keep the compiler happy.
    }


    // This operator takes a point and applies the coordinate
    // transform, returning the result.
    template <class Type>
    Vector3D<Type>
    Transform3D<Type>::
    operator*(Vector3D<Type> const& vector0) const
    {
      return Vector3D<Type>(
        m_00 * vector0.x() + m_01 * vector0.y() + m_02 * vector0.z() + m_03,
        m_10 * vector0.x() + m_11 * vector0.y() + m_12 * vector0.z() + m_13,
        m_20 * vector0.x() + m_21 * vector0.y() + m_22 * vector0.z() + m_23,
        m_30 * vector0.x() + m_31 * vector0.y() + m_32 * vector0.z() + m_33);
    }


    // The assignment operator simply duplicates its argument.
    template <class Type>
    Transform3D<Type>&
    Transform3D<Type>::
    operator=(const Transform3D<Type>& source)
    {
      if(&source != this) {
        m_00 = source.m_00; m_01 = source.m_01;
        m_02 = source.m_02; m_03 = source.m_03;
        m_10 = source.m_10; m_11 = source.m_11;
        m_12 = source.m_12; m_13 = source.m_13;
        m_20 = source.m_20; m_21 = source.m_21;
        m_22 = source.m_22; m_23 = source.m_23;
        m_30 = source.m_30; m_31 = source.m_31;
        m_32 = source.m_32; m_33 = source.m_33;
      }
      return *this;
    }


    template <class Type>
    void
    Transform3D<Type>::
    normalize()
    {
      if(m_33 == Type(0.0)) {
        BRICK_THROW(common::ValueException, "Transform3D::normalize()",
                  "Invalid normalization constant. "
                  "The bottom right element of a homogeneous transformation "
                  "cannot be equal to 0.0.");
      }
      if(m_33 != Type(1.0)) {
        m_00 /= m_33;
        m_01 /= m_33;
        m_02 /= m_33;
        m_03 /= m_33;
        m_10 /= m_33;
        m_11 /= m_33;
        m_12 /= m_33;
        m_13 /= m_33;
        m_20 /= m_33;
        m_21 /= m_33;
        m_22 /= m_33;
        m_23 /= m_33;
        m_30 /= m_33;
        m_31 /= m_33;
        m_32 /= m_33;
        m_33 /= m_33;
      }
    }


    /* ============== Non-member functions ============== */


    // This operator composes two Transform3D instances.  The resulting
    // transform satisfies the equation:
    //   (transform0 * transform1) * v0 = transform0 * (transform1 * v0),
    // where v0 is a Vector3D instance.
    //
    //
    template <class Type>
    Transform3D<Type>
    operator*(const Transform3D<Type>& transform0, const Transform3D<Type>& transform1)
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
      Type a30 = (transform0.value<3, 0>() * transform1.value<0, 0>()
                    + transform0.value<3, 1>() * transform1.value<1, 0>()
                    + transform0.value<3, 2>() * transform1.value<2, 0>()
                    + transform0.value<3, 3>() * transform1.value<3, 0>());
      Type a31 = (transform0.value<3, 0>() * transform1.value<0, 1>()
                    + transform0.value<3, 1>() * transform1.value<1, 1>()
                    + transform0.value<3, 2>() * transform1.value<2, 1>()
                    + transform0.value<3, 3>() * transform1.value<3, 1>());
      Type a32 = (transform0.value<3, 0>() * transform1.value<0, 2>()
                    + transform0.value<3, 1>() * transform1.value<1, 2>()
                    + transform0.value<3, 2>() * transform1.value<2, 2>()
                    + transform0.value<3, 3>() * transform1.value<3, 2>());
      Type a33 = (transform0.value<3, 0>() * transform1.value<0, 3>()
                    + transform0.value<3, 1>() * transform1.value<1, 3>()
                    + transform0.value<3, 2>() * transform1.value<2, 3>()
                    + transform0.value<3, 3>() * transform1.value<3, 3>());
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
      Type a30 = (transform0(3, 0) * transform1(0, 0)
                  + transform0(3, 1) * transform1(1, 0)
                  + transform0(3, 2) * transform1(2, 0)
                  + transform0(3, 3) * transform1(3, 0));
      Type a31 = (transform0(3, 0) * transform1(0, 1)
                  + transform0(3, 1) * transform1(1, 1)
                  + transform0(3, 2) * transform1(2, 1)
                  + transform0(3, 3) * transform1(3, 1));
      Type a32 = (transform0(3, 0) * transform1(0, 2)
                  + transform0(3, 1) * transform1(1, 2)
                  + transform0(3, 2) * transform1(2, 2)
                  + transform0(3, 3) * transform1(3, 2));
      Type a33 = (transform0(3, 0) * transform1(0, 3)
                  + transform0(3, 1) * transform1(1, 3)
                  + transform0(3, 2) * transform1(2, 3)
                  + transform0(3, 3) * transform1(3, 3));
#endif /* #if 0 */
      return Transform3D<Type>(a00, a01, a02, a03,
                         a10, a11, a12, a13,
                         a20, a21, a22, a23,
                         a30, a31, a32, a33);
    }


    template <class Type>
    std::ostream&
    operator<<(std::ostream& stream, const Transform3D<Type>& transform0)
    {
      // We'd rather use the templated getValue() member here, as it's
      // faster than operator()(), but for some reason the compiler is
      // choking on it, so we default to what works.  We repair most
      // of the damage by specializing for double and float below.
#if 0
      stream << "Transform3D("
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
             << transform0.value<2, 3>() << ",\n"
             << transform0.value<3, 0>() << ", "
             << transform0.value<3, 1>() << ", "
             << transform0.value<3, 2>() << ", "
             << transform0.value<3, 3>() << ")";
#else /* #if 0 */
      stream << "Transform3D("
             << transform0(0, 0) << ", "
             << transform0(0, 1) << ", "
             << transform0(0, 2) << ", "
             << transform0(0, 3) << ",\n"
             << transform0(1, 0) << ", "
             << transform0(1, 1) << ", "
             << transform0(1, 2) << ", "
             << transform0(1, 3) << ",\n"
             << transform0(2, 0) << ", "
             << transform0(2, 1) << ", "
             << transform0(2, 2) << ", "
             << transform0(2, 3) << ",\n"
             << transform0(3, 0) << ", "
             << transform0(3, 1) << ", "
             << transform0(3, 2) << ", "
             << transform0(3, 3) << ")";
#endif /* #if 0 */
      return stream;
    }


    template <class Type>
    std::istream&
    operator>>(std::istream& stream, Transform3D<Type>& transform0)
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
        common::Expect::FormatFlag flags = common::Expect::SkipWhitespace();

        // Skip any preceding whitespace.
        stream >> common::Expect("", flags);

        // Read the "Transform3D(" part.
        stream >> common::Expect("Transform3D(", flags);

        // Read all the data except the last element.
        std::vector<Type> inputValues(16);
        for(size_t index = 0; index < (inputValues.size() - 1); ++index) {
          // Read the value.
          stream >> inputValues[index];

          // Read punctuation before the next value.
          stream >> common::Expect(",", flags);
        }

        // Read the final value.
        stream >> inputValues[inputValues.size() - 1];

        // Read the closing parenthesis.
        stream >> common::Expect(")", flags);

        // And update the transform.
        transform0.setTransform(inputValues[0], inputValues[1],
                                inputValues[2], inputValues[3],
                                inputValues[4], inputValues[5],
                                inputValues[6], inputValues[7],
                                inputValues[8], inputValues[9],
                                inputValues[10], inputValues[11],
                                inputValues[12], inputValues[13],
                                inputValues[14], inputValues[15]);
      } catch(std::ios_base::failure) {
        // Empty
      }
      stream.exceptions(oldExceptionState);
      return stream;
    }

  } // namespace numeric

} // namespace brick

#endif /* #ifndef BRICK_NUMERIC_TRANSFORM3D_IMPL_HH */
