/**
***************************************************************************
* @file brick/numeric/amanatidesWoo3DIterator_impl.hh
*
* Header file defining inline and template functions declared in
* AmanatidesWoo3DIterator.hh.
*
* Copyright (C) 2004-2007,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_AMANATIDESWOO3DITERATOR_IMPL_HH
#define BRICK_NUMERIC_AMANATIDESWOO3DITERATOR_IMPL_HH

#include <iostream>
#include <brick/common/exception.hh>

// This file is included by amanatidesWoo3DIterator.hh, and should not
// be directly included by user code, so no need to include
// amanatidesWoo3DIterator.hh here.
// 
// #include <brick/numeric/amanatidesWoo3DIterator.hh>

namespace brick {

  namespace numeric {
    
    // The class constructor is initialized with all of the internal
    // variables of the voxel traversal algorithm, plus the starting
    // value of the ray parameter.
    template<class ARRAY3D>
    AmanatidesWoo3DIterator<ARRAY3D>::
    AmanatidesWoo3DIterator(ARRAY3D& data,
                            int startU, int startV, int startW,
                            int stepU, int stepV, int stepW,
                            double tMaxU, double tMaxV, double tMaxW,
                            double tDeltaU, double tDeltaV, double tDeltaW,
                            double tStart)
      : m_data(data),
        m_inBounds(true),
        m_stepU(stepU),
        m_stepV(stepV),
        m_stepW(stepW),
        m_tDeltaU(tDeltaU),
        m_tDeltaV(tDeltaV),
        m_tDeltaW(tDeltaW),
        m_tEntry(tStart),
        m_tMaxU(tMaxU),
        m_tMaxV(tMaxV),
        m_tMaxW(tMaxW),
        m_U(startU),
        m_uLimit(stepU > 0 ? static_cast<int>(data.shape()[2]) : -1),
        m_V(startV),
        m_vLimit(stepV > 0 ? static_cast<int>(data.shape()[1]) : -1),
        m_W(startW),
        m_wLimit(stepW > 0 ? static_cast<int>(data.shape()[0]) : -1)
    {
      if((m_U < 0)
         || (m_V < 0)
         || (m_W < 0)
         || (m_U >= static_cast<int>(data.shape()[2]))
         || (m_V >= static_cast<int>(data.shape()[1]))
         || (m_W >= static_cast<int>(data.shape()[0]))) {
        m_inBounds = false;      
      }
    }

    // Copy constructor.
    template<class ARRAY3D>
    AmanatidesWoo3DIterator<ARRAY3D>::
    AmanatidesWoo3DIterator(const AmanatidesWoo3DIterator& source)
      : m_data(source.m_data),
        m_inBounds(source.m_inBounds),
        m_stepU(source.m_stepU),
        m_stepV(source.m_stepV),
        m_stepW(source.m_stepW),
        m_tDeltaU(source.m_tDeltaU),
        m_tDeltaV(source.m_tDeltaV),
        m_tDeltaW(source.m_tDeltaW),
        m_tEntry(source.m_tEntry),
        m_tMaxU(source.m_tMaxU),
        m_tMaxV(source.m_tMaxV),
        m_tMaxW(source.m_tMaxW),
        m_U(source.m_U),
        m_uLimit(source.m_uLimit),
        m_V(source.m_V),
        m_vLimit(source.m_vLimit),
        m_W(source.m_W),
        m_wLimit(source.m_wLimit)
    {
      // Empty
    }

    // This operator returns a reference to the Array3D element at the
    // current voxel.
    template<class ARRAY3D>
    inline typename ARRAY3D::value_type& // element_type?
    AmanatidesWoo3DIterator<ARRAY3D>::
    operator*()
    {
      return m_data(m_W, m_V, m_U);
    }

    // This operator returns a pointer to the Array3D element at the
    // current voxel.
    template<class ARRAY3D>
    inline typename ARRAY3D::value_type* // element_type?
    AmanatidesWoo3DIterator<ARRAY3D>::
    operator->()
    {
      return &(this->operator*());
    }

    // The pre-increment operator increments the iterator so that it
    // points to the next voxel along the path.
    template<class ARRAY3D>
    inline AmanatidesWoo3DIterator<ARRAY3D>&
    AmanatidesWoo3DIterator<ARRAY3D>::
    operator++()
    {
      if(m_tMaxU < m_tMaxV) {
        if(m_tMaxW < m_tMaxU) {
          m_tEntry = m_tMaxW;
          m_W += m_stepW;
          m_tMaxW += m_tDeltaW;
          if(m_W == m_wLimit) {
            m_inBounds = false;
          }
        } else {
          m_tEntry = m_tMaxU;
          m_U += m_stepU;
          m_tMaxU += m_tDeltaU;
          if(m_U == m_uLimit) {
            m_inBounds = false;
          }
        }
      } else {
        if(m_tMaxW < m_tMaxV) {
          m_tEntry = m_tMaxW;
          m_W += m_stepW;
          m_tMaxW += m_tDeltaW;
          if(m_W == m_wLimit) {
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
      }
      return *this;
    }

    // The post-increment operator increments the iterator so that it
    // points to the next voxel along the path.  It differs from the
    // pre-increment operator in its return value.
    template<class ARRAY3D>
    inline AmanatidesWoo3DIterator<ARRAY3D>
    AmanatidesWoo3DIterator<ARRAY3D>::
    operator++(int)
    {
      AmanatidesWoo3DIterator<ARRAY3D> thisCopy(*this);
      ++this;
      return thisCopy;
    }

    // This is the assignment operator.  It copies the value of its
    // argument into *this.
    template<class ARRAY3D>
    AmanatidesWoo3DIterator<ARRAY3D>&
    AmanatidesWoo3DIterator<ARRAY3D>::
    operator=(const AmanatidesWoo3DIterator& source)
    {
      m_data = source.m_data;
      m_inBounds = source.m_inBounds;
      m_stepU = source.m_stepU;
      m_stepV = source.m_stepV;
      m_stepW = source.m_stepW;
      m_tDeltaU = source.m_tDeltaU;
      m_tDeltaV = source.m_tDeltaV;
      m_tDeltaW = source.m_tDeltaW;
      m_tEntry = source.m_tEntry;
      m_tMaxU = source.m_tMaxU;
      m_tMaxV = source.m_tMaxV;
      m_tMaxW = source.m_tMaxW;
      m_U = source.m_U;
      m_uLimit = source.m_uLimit;
      m_V = source.m_V;
      m_vLimit = source.m_vLimit;
      m_W = source.m_W;
      m_wLimit = source.m_wLimit;
      return *this;
    }

    // The equality operator returns true if *this currently
    // references the same voxel as the argument.
    template<class ARRAY3D>
    inline bool
    AmanatidesWoo3DIterator<ARRAY3D>::
    operator==(const AmanatidesWoo3DIterator& other)
    {
      // Return true if both refer to valid voxels or if both refer to
      // invalid voxels.
      return (m_inBounds == other.m_inBounds);
    }

    // The inequality operator returns true if *this currently
    // references the a different voxel than the argument.
    template<class ARRAY3D>
    inline bool
    AmanatidesWoo3DIterator<ARRAY3D>::
    operator!=(const AmanatidesWoo3DIterator& other)
    {
      return !(this->operator==(other));
    }
  
  } // namespace numeric

} // namespace brick

#endif /* #ifndef BRICK_NUMERIC_AMANATIDESWOO3DITERATOR_IMPL_HH */
