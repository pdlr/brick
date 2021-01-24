/**
***************************************************************************
* @file brick/numeric/quaternion.cpp
*
* Source file defining Quaternion class.
*
* (C) Copyright 1996-2011 David LaRose, dlr@davidlarose.com
* See accompanying LICENSE file for details.
*
***************************************************************************
**/


// This file is included by quaternion.hh, and should not be directly
// included by user code, so no need to include quaternion.hh here.
//
// #include <brick/numeric/quaternion.hh>

#include <cmath>

#include <brick/numeric/mathFunctions.hh>

namespace brick {

  namespace numeric {


    // This member function computes and returns the magnitude of the
    // quaternion.
    template <class Type>
    Type
    Quaternion<Type>::
    getMagnitude()
    {
      return brick::numeric::squareRoot(m_s*m_s + m_i*m_i + m_j*m_j + m_k*m_k);
    }


    // This member function normalizes the Quaternion, first computing
    // the magnitude of the Quaternion, and then dividing each element
    // by that value.
    template <class Type>
    void
    Quaternion<Type>::
    normalize()
    {
      // Skip all the effort if the Quaternion already has unit
      // magnitude.
      if(!m_isNormalized) {
	// First compute the magnitude of the Quaternion.
	Type norm = this->getMagnitude();

	// Next some error checking.
	if(norm == Type(0.0)) {
	  // BRICK_THROW(ValueException, "Quaternion::normalize()",
	  //           "Can't normalize a zero magnitude Quaternion.");
	  m_s = Type(1.0);
	  m_i = Type(0.0);
	  m_j = Type(0.0);
	  m_k = Type(0.0);
	  m_isNormalized = true;
	} else {
	  // Finally, normalize.
	  m_s /= norm;
	  m_i /= norm;
	  m_j /= norm;
	  m_k /= norm;
	  m_isNormalized = true;
	}
      }
    }

  } // namespace numeric

} // namespace brick
