/**
***************************************************************************
* @file brick/numeric/transform3D.cc
*
* Source file providing implementations for symbols declared in
* transform3D.hh.
*
* Copyright (C) 2001-2013 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/numeric/transform3D.hh>

namespace brick {

  namespace numeric {

    template <>
    inline Transform3D<double>
    operator*(const Transform3D<double>& transform0,
              const Transform3D<double>& transform1)
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
      double a30 = (transform0.value<3, 0>() * transform1.value<0, 0>()
                    + transform0.value<3, 1>() * transform1.value<1, 0>()
                    + transform0.value<3, 2>() * transform1.value<2, 0>()
                    + transform0.value<3, 3>() * transform1.value<3, 0>());
      double a31 = (transform0.value<3, 0>() * transform1.value<0, 1>()
                    + transform0.value<3, 1>() * transform1.value<1, 1>()
                    + transform0.value<3, 2>() * transform1.value<2, 1>()
                    + transform0.value<3, 3>() * transform1.value<3, 1>());
      double a32 = (transform0.value<3, 0>() * transform1.value<0, 2>()
                    + transform0.value<3, 1>() * transform1.value<1, 2>()
                    + transform0.value<3, 2>() * transform1.value<2, 2>()
                    + transform0.value<3, 3>() * transform1.value<3, 2>());
      double a33 = (transform0.value<3, 0>() * transform1.value<0, 3>()
                    + transform0.value<3, 1>() * transform1.value<1, 3>()
                    + transform0.value<3, 2>() * transform1.value<2, 3>()
                    + transform0.value<3, 3>() * transform1.value<3, 3>());
      return Transform3D<double>(a00, a01, a02, a03,
                                 a10, a11, a12, a13,
                                 a20, a21, a22, a23,
                                 a30, a31, a32, a33);
    }


    template <>
    inline Transform3D<float>
    operator*(const Transform3D<float>& transform0,
              const Transform3D<float>& transform1)
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
      float a30 = (transform0.value<3, 0>() * transform1.value<0, 0>()
                   + transform0.value<3, 1>() * transform1.value<1, 0>()
                   + transform0.value<3, 2>() * transform1.value<2, 0>()
                   + transform0.value<3, 3>() * transform1.value<3, 0>());
      float a31 = (transform0.value<3, 0>() * transform1.value<0, 1>()
                   + transform0.value<3, 1>() * transform1.value<1, 1>()
                   + transform0.value<3, 2>() * transform1.value<2, 1>()
                   + transform0.value<3, 3>() * transform1.value<3, 1>());
      float a32 = (transform0.value<3, 0>() * transform1.value<0, 2>()
                   + transform0.value<3, 1>() * transform1.value<1, 2>()
                   + transform0.value<3, 2>() * transform1.value<2, 2>()
                   + transform0.value<3, 3>() * transform1.value<3, 2>());
      float a33 = (transform0.value<3, 0>() * transform1.value<0, 3>()
                   + transform0.value<3, 1>() * transform1.value<1, 3>()
                   + transform0.value<3, 2>() * transform1.value<2, 3>()
                   + transform0.value<3, 3>() * transform1.value<3, 3>());
      return Transform3D<float>(a00, a01, a02, a03,
                                a10, a11, a12, a13,
                                a20, a21, a22, a23,
                                a30, a31, a32, a33);
    }


#if 0
    template <>
    std::ostream&
    operator<<(std::ostream& stream, const Transform3D<double>& transform0)
    {
      stream << "Transform3D<double>("
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
      return stream;
    }


    template <>
    std::ostream&
    operator<<(std::ostream& stream, const Transform3D<float>& transform0)
    {
      stream << "Transform3D<float>("
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
      return stream;
    }
#endif /* #if 0 */


  } // namespace numeric

} // namespace brick
