/**
***************************************************************************
* @file brick/numeric/quaternion.hh
* Header file declaring Quaternion class template.
*
* (C) Copyright 1996-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_QUATERNION_HH
#define BRICK_NUMERIC_QUATERNION_HH

namespace brick {

  namespace numeric {
    
    /**
     ** This class implements a quaternion.  Among other things,
     ** quaternions are useful for expressing rigid rotations in 3D
     ** space.  To do this, people generally Set the elements of the
     ** quaternion equal to [c, s*x, s*y, s*z], where c = cos(theta/2),
     ** s = sin(theta/2), theta is the angle of rotation, and [x, y, z]
     ** is a unit vector along the axis of rotation.  Note that in this
     ** representation of rotation, the quaternion must have unit
     ** magnitude.  This class does not enforce unit magnitude.  If you
     ** want the quaternion to have unit magnitude, you must either set
     ** the values appropriately, or call the normalize() method.
     **/
    template <class Type>
    class Quaternion {
    public:

      /** 
       * The default constructor initializes a quaternion with values
       * [1.0, 0.0, 0.0, 0.0].
       */
      Quaternion() 
        : m_s(1.0), m_i(0.0), m_j(0.0), m_k(0.0), m_isNormalized(true) {}

      /** 
       * This constructor explicitly sets the values of the quaternion.
       * 
       * @param ss This parameter specifies the value of s, the real
       * component of the quaternion.
       *
       * @param ii This parameter specifies the value of i, the first
       * imaginary component of the quaternion.
       *
       * @param jj This parameter specifies the value of j, the second
       * imaginary component of the quaternion.
       *
       * @param kk This parameter specifies the value of k, the third
       * imaginary component of the quaternion.
       *
       * @param isAlreadyNormalized This parameter is used to indicate whether
       * the provided values are normalized.  If you know that the
       * resulting quaternion will have magnitude equal to one, you can
       * set this argument to true.
       */
      Quaternion(double ss, double ii, double jj, double kk,
                 bool isAlreadyNormalized=false) 
        : m_s(ss), m_i(ii), m_j(jj), m_k(kk),
          m_isNormalized(isAlreadyNormalized) {}

      /** 
       * The copy constructor deep copies its argument.
       * 
       * @param source This parameter specifies the Quaternion instance
       * to be copied.
       */
      Quaternion(const Quaternion<Type> &source) :
        m_s(source.m_s), m_i(source.m_i), m_j(source.m_j), m_k(source.m_k),
        m_isNormalized(source.m_isNormalized) {}

      /** 
       * The destructor destroys the Quaternion instance.
       */
      virtual
      inline
      ~Quaternion() {}

      /** 
       * This member function returns the real component of the
       * Quaternion.
       * 
       * @return The return value is the real component of the
       * Quaternion.
       */
      inline
      double
      s() const { return m_s; };

      /** 
       * This member function returns the first imaginary component of
       * the Quaternion.
       * 
       * @return The return value is the first imaginary component of
       * the Quaternion.
       */
      inline
      double
      i() const { return m_i; };

      /** 
       * This member function returns the second imaginary component of
       * the Quaternion.
       * 
       * @return The return value is the second imaginary component of
       * the Quaternion.
       */
      inline
      double
      j() const { return m_j; };

      /** 
       * This member function returns the third imaginary component of
       * the Quaternion.
       * 
       * @return The return value is the third imaginary component of
       * the Quaternion.
       */
      inline
      double
      k() const { return m_k; };

      /** 
       * This member function normalizes the Quaternion, first computing
       * the magnitude of the Quaternion, and then dividing each element
       * by that value.  Strictly speaking, the magnitude of a
       * Quaternion q is equal to dot(q, conjugate(q)).  This is the
       * same as the magnitude of a 4 element vector [s, i, j, k].
       */
      bool
      isNormalized() {return m_isNormalized;}


      /** 
       * This member function computes and returns the magnitude of
       * the quaternion.
       * 
       * @return The return value is the computed magnitude.
       */
      Type
      getMagnitude();

      
      /** 
       * This member function normalizes the Quaternion, first computing
       * the magnitude of the Quaternion, and then dividing each element
       * by that value.  Strictly speaking, the magnitude of a
       * Quaternion q is equal to dot(q, conjugate(q)).  This is the
       * same as the magnitude of a 4 element vector [s, i, j, k].
       */
      void
      normalize();

      /** 
       * This member function sets the values of the Quaternion
       * components explicitly.
       * 
       * @param ss This parameter specifies the real component of the
       * Quaternion.
       *
       * @param ii This parameter specifies the first imaginary component
       * of the Quaternion.
       *
       * @param jj This parameter specifies the second imaginary
       * component of the Quaternion.
       *
       * @param kk This parameter specifies the third imaginary component
       * of the Quaternion.
       *
       * @param isAlreadyNormalized This parameter is used to indicate
       * whether the provided values are normalized.  If you know that
       * the resulting quaternion will have magnitude equal to one,
       * you can set this argument to true.
       */
      inline
      void
      setValue(double ss, double ii, double jj, double kk,
               bool isAlreadyNormalized = false) {
        m_s = ss; m_i = ii; m_j = jj; m_k = kk;
        m_isNormalized = isAlreadyNormalized;
      }

      /** 
       * The assignment operator simply deep copies its argument.
       * 
       * @param source This argument is the Quaternion instance to be copied.
       *
       * @return The return value is a reference to *this.
       */
      inline
      Quaternion<Type>&
      operator=(const Quaternion<Type> &source) {
        setValue(source.m_s, source.m_i, source.m_j, source.m_k);
        m_isNormalized = source.m_isNormalized;
        return *this;
      }
    
    private:
      // Components of the quaternion are listed in order.
      double m_s;
      double m_i;
      double m_j;
      double m_k;

      // A flag to avoid repeatedly re-normalizing.
      bool m_isNormalized;
    };

    /* ======================= Non-member functions ===================== */

    /** 
     * This function returns the conjugate of a Quaternion, in which the
     * sign of each imaginary component has been reversed.  Note that
     * for quaternions which represent rotation, the rotation
     * represented by the conjugate is also the inverse rotation.  That
     * is, if you rotate by a quaternion, q, and then rotate by
     * conjugate(q), you wind up exactly where you started.
     * 
     * @param source This argument specifies the Quaternion whose
     * conjugate is to be computed.
     *
     * @return The return value is a Quaternion which is conjugate to
     * the argument.
     */
    template <class Type>
    inline
    Quaternion<Type>
    conjugate(const Quaternion<Type>& source) {
      return Quaternion<Type>(source.s(), -(source.i()), -(source.j()),
                              -(source.k()));
    }
  
  } // namespace numeric

} // namespace brick

// Implementation lives in a different include file.
#include <brick/numeric/quaternion_impl.hh>

#endif /* #ifndef BRICK_NUMERIC_QUATERNION_HH */
