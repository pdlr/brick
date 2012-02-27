/**
***************************************************************************
* @file brick/numeric/amanatidesWoo2DIterator.hh
*
* Header file defining inline and template functions declared in
* AmanatidesWoo2DIterator.hh.
*
* Copyright (C) 2004-2007,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_AMANATIDESWOO2DITERATOR_IMPL_HH
#define BRICK_NUMERIC_AMANATIDESWOO2DITERATOR_IMPL_HH

// This file is included by amanatidesWoo2DIterator.hh, and should not
// be directly included by user code, so no need to include
// amanatidesWoo2DIterator.hh here.
// 
// #include <brick/numeric/amanatidesWoo2DIterator.hh>

namespace brick {

  namespace numeric {
    
    // The class constructor is initialized with all of the internal
    // variables of the voxel traversal algorithm, plus the starting
    // value of the ray parameter.
    template <class ARRAY2D, class FLOAT_TYPE>
    AmanatidesWoo2DIterator<ARRAY2D, FLOAT_TYPE>::
    AmanatidesWoo2DIterator(ARRAY2D& data,
                            int startU, int startV,
                            int stepU, int stepV,
                            FLOAT_TYPE tMaxU, FLOAT_TYPE tMaxV,
                            FLOAT_TYPE tDeltaU, FLOAT_TYPE tDeltaV,
                            FLOAT_TYPE tStart)
      : m_data(data),
        m_inBounds(true),
        m_stepU(stepU),
        m_stepV(stepV),
        m_tDeltaU(tDeltaU),
        m_tDeltaV(tDeltaV),
        m_tEntry(tStart),
        m_tMaxU(tMaxU),
        m_tMaxV(tMaxV),
        m_U(startU),
        m_uLimit(stepU > 0 ? static_cast<int>(data.columns()) : -1),
        m_V(startV),
        m_vLimit(stepV > 0 ? static_cast<int>(data.rows()) : -1)
    {
      if((m_U < 0)
         || (m_V < 0)
         || (m_U >= static_cast<int>(m_data.columns()))
         || (m_V >= static_cast<int>(m_data.rows()))) {
        m_inBounds = false;      
      }
    }

    // Copy constructor.
    template<class ARRAY2D, class FLOAT_TYPE>
    AmanatidesWoo2DIterator<ARRAY2D, FLOAT_TYPE>::
    AmanatidesWoo2DIterator(const AmanatidesWoo2DIterator& source)
      : m_data(source.m_data),
        m_inBounds(source.m_inBounds),
        m_stepU(source.m_stepU),
        m_stepV(source.m_stepV),
        m_tDeltaU(source.m_tDeltaU),
        m_tDeltaV(source.m_tDeltaV),
        m_tEntry(source.m_tEntry),
        m_tMaxU(source.m_tMaxU),
        m_tMaxV(source.m_tMaxV),
        m_U(source.m_U),
        m_uLimit(source.m_uLimit),
        m_V(source.m_V),
        m_vLimit(source.m_vLimit)
    {
      // Empty
    }

    // This operator returns a reference to the Array2D element at the
    // current pixel.
    template<class ARRAY2D, class FLOAT_TYPE>
    inline typename ARRAY2D::value_type& // element_type?
    AmanatidesWoo2DIterator<ARRAY2D, FLOAT_TYPE>::
    operator*()
    {
      return m_data(m_V, m_U);
    }

    // This operator returns a pointer to the Array2D element at the
    // current pixel.
    template<class ARRAY2D, class FLOAT_TYPE>
    inline typename ARRAY2D::value_type* // element_type?
    AmanatidesWoo2DIterator<ARRAY2D, FLOAT_TYPE>::
    operator->()
    {
      return &(this->operator*());
    }

    // The pre-increment operator increments the iterator so that it
    // points to the next pixel along the path.
    template<class ARRAY2D, class FLOAT_TYPE>
    inline AmanatidesWoo2DIterator<ARRAY2D, FLOAT_TYPE>&
    AmanatidesWoo2DIterator<ARRAY2D, FLOAT_TYPE>::
    operator++()
    {
      if(m_tMaxU < m_tMaxV) {
        m_tEntry = m_tMaxU;
        m_U += m_stepU;
        m_tMaxU += m_tDeltaU;
        if(m_U == m_uLimit) {
          m_inBounds = false;
        }
      } else {
        m_tEntry = m_tMaxV;
        m_V += m_stepV;
        m_tMaxV += m_tDeltaV;
        if(m_V == m_vLimit) {
          m_inBounds = false;
        }
      }
      return *this;
    }

    // The post-increment operator increments the iterator so that it
    // points to the next pixel along the path.  It differs from the
    // pre-increment operator in its return value.
    template<class ARRAY2D, class FLOAT_TYPE>
    inline AmanatidesWoo2DIterator<ARRAY2D, FLOAT_TYPE>
    AmanatidesWoo2DIterator<ARRAY2D, FLOAT_TYPE>::
    operator++(int)
    {
      AmanatidesWoo2DIterator<ARRAY2D, FLOAT_TYPE> thisCopy(*this);
      ++this;
      return thisCopy;
    }

    // This is the assignment operator.  It copies the value of its
    // argument into *this.
    template<class ARRAY2D, class FLOAT_TYPE>
    AmanatidesWoo2DIterator<ARRAY2D, FLOAT_TYPE>&
    AmanatidesWoo2DIterator<ARRAY2D, FLOAT_TYPE>::
    operator=(const AmanatidesWoo2DIterator& source)
    {
      m_data = source.m_data;
      m_inBounds = source.m_inBounds;
      m_stepU = source.m_stepU;
      m_stepV = source.m_stepV;
      m_tDeltaU = source.m_tDeltaU;
      m_tDeltaV = source.m_tDeltaV;
      m_tEntry = source.m_tEntry;
      m_tMaxU = source.m_tMaxU;
      m_tMaxV = source.m_tMaxV;
      m_U = source.m_U;
      m_uLimit = source.m_uLimit;
      m_V = source.m_V;
      m_vLimit = source.m_vLimit;
      return *this;
    }

    // The equality operator returns true if *this currently
    // references the same pixel as the argument.
    template<class ARRAY2D, class FLOAT_TYPE>
    inline bool
    AmanatidesWoo2DIterator<ARRAY2D, FLOAT_TYPE>::
    operator==(const AmanatidesWoo2DIterator& other)
    {
      // Return true if both refer to valid pixels or if both refer to
      // invalid pixels.
      return (m_inBounds == other.m_inBounds);
    }

    // The inequality operator returns true if *this currently
    // references the a different pixel than the argument.
    template<class ARRAY2D, class FLOAT_TYPE>
    inline bool
    AmanatidesWoo2DIterator<ARRAY2D, FLOAT_TYPE>::
    operator!=(const AmanatidesWoo2DIterator& other)
    {
      return !(this->operator==(other));
    }
  
  } // namespace numeric

} // namespace brick

#endif /* #ifndef BRICK_NUMERIC_AMANATIDESWOO2DITERATOR_IMPL_HH */
