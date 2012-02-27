/**
***************************************************************************
* @file brick/numeric/amanatidesWoo2D_impl.hh
*
* Header file defining inline and template functions declared in
* AmanatidesWoo2D.hh.
*
* Copyright (C) 2004-2007,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_AMANATIDESWOO2D_IMPL_HH
#define BRICK_NUMERIC_AMANATIDESWOO2D_IMPL_HH

// This file is included by amanatidesWoo2D.hh, and should not be
// directly included by user code, so no need to include
// amanatidesWoo2D.hh here.
// 
// #include <brick/numeric/amanatidesWoo2D.hh>

namespace brick {

  namespace numeric {
    
    template <class ARRAY2D, class FLOAT_TYPE, class INT_TYPE>
    AmanatidesWoo2D<ARRAY2D, FLOAT_TYPE, INT_TYPE>::
    AmanatidesWoo2D(ARRAY2D& data,
                    const Transform2D<FLOAT_TYPE>& pixelTworld,
                    const Vector2D<FLOAT_TYPE>& rayOrigin,
                    const Vector2D<FLOAT_TYPE>& rayDirection,
                    bool downstreamOnly)
      : m_data(data),
        m_initialU(-1),  // Initialize to illegal values.  These will be 
        m_initialV(-1),  // replaced with legal ones if there's a valid
        // intersection between the ray and the pixel
        // array.
        m_stepU(),
        m_stepV(),
        m_tDeltaU(),
        m_tDeltaV(),
        m_tMaxU(),
        m_tMaxV(),
        m_tStart(),
        m_validIntersection(true)    
    {
      // First convert everything into pixel coordinates.
      Vector2D<FLOAT_TYPE> rayOriginPixel = pixelTworld * rayOrigin;
      Vector2D<FLOAT_TYPE> rayDirectionPixel =
        (pixelTworld * (rayOrigin + rayDirection)) - rayOriginPixel;

      // Now find points entry and exit from the CT volume.  These are
      // expressed in as parameter values tEntry and tExit such that
      // (rayOriginPixel + tEntry * rayDirectionPixel) is the entry
      // point and (rayOriginPixel + tExit * rayDirectionPixel) is the
      // exit point.
      std::pair<FLOAT_TYPE, FLOAT_TYPE> tEntry_tExit =
        this->findEntryAndExitPoints(rayOriginPixel, rayDirectionPixel, m_data);

      // Sometimes we want to disallow any pixels which are "behind" the
      // ray origin.  That is, sometimes we want to only trace those
      // parts of the ray which correspond to positive values of the
      // parameter t.
      if(downstreamOnly && (tEntry_tExit.first < 0.0)) {
        tEntry_tExit.first = 0.0;
      }

      // Sometimes the there's no intersection with the pixel array.  In
      // this case, tExit will be less than tEntry, and we don't bother
      // doing any more calculation.
      if(tEntry_tExit.second <= tEntry_tExit.first) {
        m_validIntersection = false;
        return;
      }

      // Now make the point of entry explicit.
      Vector2D<FLOAT_TYPE> entryPoint =
        rayOriginPixel + tEntry_tExit.first * rayDirectionPixel;

      // Correct for rounding error which sometimes makes entryPoint
      // have values like -1.0e-14.
      if(entryPoint.x() < 0.0) {
        entryPoint.setX(0.0);
      }
      if(entryPoint.y() < 0.0) {
        entryPoint.setY(0.0);
      }
    
      // Finally, assign the variables described in the Amanatides' and
      // Woo's paper.

      // Since we've already converted to pixel coords, finding the
      // first pixel on the path is trivial.
      m_initialU = static_cast<INT_TYPE>(entryPoint.x());
      m_initialV = static_cast<INT_TYPE>(entryPoint.y());

      // There's an similar case to the one handled by the rounding
      // error test above.  Member function findEntryAndExitPoints()
      // uses m_data.columns() and m_data.rows() to define the maximum U
      // and V coordinates, respectively.  This means that it's very
      // possible for m_initialU to be equal to m_data.columns(), which
      // will cause an out-of-bounds index. Correct for this now.
      if(m_initialU == static_cast<INT_TYPE>(m_data.columns())) {
        --m_initialU;
      }
      if(m_initialV == static_cast<INT_TYPE>(m_data.rows())) {
        --m_initialV;
      }

      // Sanity check.  This sometimes fails if the ray is parallel to
      // one axis of the data array.
      if((m_initialU >= static_cast<INT_TYPE>(m_data.columns()))
         || (m_initialV >= static_cast<INT_TYPE>(m_data.rows()))) {
        m_validIntersection = false;
        return;
      }
    
      // m_tStart is just the same as tEntry.
      m_tStart = tEntry_tExit.first;
    
      // The remaining variables depend on whether U & V will be
      // increasing or decreasing as we travel along the ray, so we need
      // if clauses.  Please see the declaration for documentation on
      // what each of these member variables means.
      if(rayDirectionPixel.x() > 0.0) {
        m_stepU = 1;
        m_tDeltaU = 1.0 / rayDirectionPixel.x();
        m_tMaxU = m_tStart + (((m_initialU + 1) - entryPoint.x())
                              / rayDirectionPixel.x());
      } else if(rayDirectionPixel.x() < 0.0) {
        m_stepU = -1;
        m_tDeltaU = -(1.0 / rayDirectionPixel.x());
        m_tMaxU = m_tStart + ((m_initialU - entryPoint.x())
                              / rayDirectionPixel.x());
      } else { // rayDirectionPixel.x() == 0.0;
        m_stepU = 0;
        m_tDeltaU = std::numeric_limits<FLOAT_TYPE>::max();
        m_tMaxU = std::numeric_limits<FLOAT_TYPE>::max();
      }
      if(rayDirectionPixel.y() > 0.0) {
        m_stepV = 1;
        m_tDeltaV = 1.0 / rayDirectionPixel.y();
        m_tMaxV = m_tStart + (((m_initialV + 1) - entryPoint.y())
                              / rayDirectionPixel.y());
      } else if(rayDirectionPixel.y() < 0.0) {
        m_stepV = -1;
        m_tDeltaV = -(1.0 / rayDirectionPixel.y());
        m_tMaxV = m_tStart + ((m_initialV - entryPoint.y())
                              / rayDirectionPixel.y());
      } else { // rayDirectionPixel.y() == 0.0;
        m_stepV = 0;
        m_tDeltaV = std::numeric_limits<FLOAT_TYPE>::max();
        m_tMaxV = std::numeric_limits<FLOAT_TYPE>::max();
      }
    }


    // The copy constructor deep copies its argument.
    template <class ARRAY2D, class FLOAT_TYPE, class INT_TYPE>
    AmanatidesWoo2D<ARRAY2D, FLOAT_TYPE, INT_TYPE>::
    AmanatidesWoo2D(const AmanatidesWoo2D& source)
      : m_data(source.m_data),
        m_initialU(source.m_initialU),
        m_initialV(source.m_initialV),
        m_stepU(source.m_stepU),
        m_stepV(source.m_stepV),
        m_tDeltaU(source.m_tDeltaU),
        m_tDeltaV(source.m_tDeltaV),
        m_tMaxU(source.m_tMaxU),
        m_tMaxV(source.m_tMaxV),
        m_tStart(source.m_tStart),
        m_validIntersection(source.m_validIntersection)
    {
      // Empty
    }

    template <class ARRAY2D, class FLOAT_TYPE, class INT_TYPE>
    AmanatidesWoo2D<ARRAY2D, FLOAT_TYPE, INT_TYPE>::
    ~AmanatidesWoo2D() {};

    template <class ARRAY2D, class FLOAT_TYPE, class INT_TYPE>
    typename AmanatidesWoo2D<ARRAY2D, FLOAT_TYPE, INT_TYPE>::iterator
    AmanatidesWoo2D<ARRAY2D, FLOAT_TYPE, INT_TYPE>::
    begin()
    {
      return AmanatidesWoo2DIterator<ARRAY2D, FLOAT_TYPE, INT_TYPE>(
        m_data, m_initialU, m_initialV, m_stepU, m_stepV, m_tMaxU, m_tMaxV,
        m_tDeltaU, m_tDeltaV, m_tStart);
    }

    template <class ARRAY2D, class FLOAT_TYPE, class INT_TYPE>
    typename AmanatidesWoo2D<ARRAY2D, FLOAT_TYPE, INT_TYPE>::iterator
    AmanatidesWoo2D<ARRAY2D, FLOAT_TYPE, INT_TYPE>::
    end()
    {
      // AmanatidesWoo2DIterator<ARRAY2D, FLOAT_TYPE, INT_TYPE>::operator==(...) considers all
      // iterators with illegal pixel coordinates to be equal, so all we
      // have to do here is return an iterator which references an
      // illegal pixel.
      return AmanatidesWoo2DIterator<ARRAY2D, FLOAT_TYPE, INT_TYPE>(
        m_data, -1, -1, m_stepU, m_stepV, m_tMaxU, m_tMaxV,
        m_tDeltaU, m_tDeltaV, m_tStart);
    }

    template <class ARRAY2D, class FLOAT_TYPE, class INT_TYPE>
    inline bool
    AmanatidesWoo2D<ARRAY2D, FLOAT_TYPE, INT_TYPE>::
    validIntersection()
    {
      return m_validIntersection;
    }

    template <class ARRAY2D, class FLOAT_TYPE, class INT_TYPE>
    AmanatidesWoo2D<ARRAY2D, FLOAT_TYPE, INT_TYPE>&
    AmanatidesWoo2D<ARRAY2D, FLOAT_TYPE, INT_TYPE>::
    operator=(const AmanatidesWoo2D& source)
    {
      m_data = source.m_data;
      m_initialU = source.m_initialU;
      m_initialV = source.m_initialV;
      m_stepU = source.m_stepU;
      m_stepV = source.m_stepV;
      m_tDeltaU = source.m_tDeltaU;
      m_tDeltaV = source.m_tDeltaV;
      m_tMaxU = source.m_tMaxU;
      m_tMaxV = source.m_tMaxV;
      m_tStart = source.m_tStart;
      m_validIntersection = source.m_validIntersection;
      return *this;
    }
    
    template <class ARRAY2D, class FLOAT_TYPE, class INT_TYPE>
    std::pair<FLOAT_TYPE, FLOAT_TYPE>
    AmanatidesWoo2D<ARRAY2D, FLOAT_TYPE, INT_TYPE>::
    findEntryAndExitPoints(const Vector2D<FLOAT_TYPE>& rayOriginPixel,
                           const Vector2D<FLOAT_TYPE>& rayDirectionPixel,
                           const ARRAY2D& data)
    {
      // First find intersection with each boundary line.

      // ... Find the intersection with the line U = 0, or else a really
      // small number if rayDirection is parallel to the rows of the
      // array.
      FLOAT_TYPE tIntersectU0 = findIntersection(
        rayOriginPixel, rayDirectionPixel, Vector2D<FLOAT_TYPE>(1.0, 0.0), 0.0,
        -std::numeric_limits<FLOAT_TYPE>::max());

      // ... Find the intersection with the line U = data.columns(),
      // or else a really big number if rayDirection is parallel to the
      // rows of the array.
      FLOAT_TYPE tIntersectU1 = findIntersection(
        rayOriginPixel, rayDirectionPixel, Vector2D<FLOAT_TYPE>(1.0, 0.0), data.columns(),
        std::numeric_limits<FLOAT_TYPE>::max());

      // ... Find the intersection with the line V = 0, or else a really
      // small number if rayDirection is parallel to the columns of the
      // array.
      FLOAT_TYPE tIntersectV0 = findIntersection(
        rayOriginPixel, rayDirectionPixel, Vector2D<FLOAT_TYPE>(0.0, 1.0), 0.0,
        -std::numeric_limits<FLOAT_TYPE>::max());
    
      // ... Find the intersection with the line V = data.rows(), or
      // else a really big number if rayDirection is parallel to the
      // columns of the array.
      FLOAT_TYPE tIntersectV1 = findIntersection(
        rayOriginPixel, rayDirectionPixel, Vector2D<FLOAT_TYPE>(0.0, 1.0), data.rows(),
        std::numeric_limits<FLOAT_TYPE>::max());

      // Now find the closer and farther of each pair of intersections.
      FLOAT_TYPE tMinU = std::min(tIntersectU0, tIntersectU1);
      FLOAT_TYPE tMinV = std::min(tIntersectV0, tIntersectV1);
      FLOAT_TYPE tMaxU = std::max(tIntersectU0, tIntersectU1);
      FLOAT_TYPE tMaxV = std::max(tIntersectV0, tIntersectV1);

      // Compute closest point which could possibly intersect with the volume.
      FLOAT_TYPE tEntry = std::max(tMinU, tMinV);
      // Compute farthest point which could possibly intersect with the volume.
      FLOAT_TYPE tExit = std::min(tMaxU, tMaxV);

      // Return our findings.
      return std::make_pair(tEntry, tExit);
    }

    template <class ARRAY2D, class FLOAT_TYPE, class INT_TYPE>
    FLOAT_TYPE
    AmanatidesWoo2D<ARRAY2D, FLOAT_TYPE, INT_TYPE>::
    findIntersection(const Vector2D<FLOAT_TYPE>& rayOrigin,
                     const Vector2D<FLOAT_TYPE>& rayDirection,
                     const Vector2D<FLOAT_TYPE>& bVector,
                     FLOAT_TYPE cConstant,
                     FLOAT_TYPE defaultValue)
    {
      // bVector and cConstant describe the desired plane:
      //   dot(x, bVector) = cConstant.
      // Now, x is constrained to lie on the specified line:
      //   x = rayOrigin + t * rayDirection.
      // Substituting, we have:
      //   dot(rayOrigin + t * rayDirection, bVector) = cConstant.
      // Which we rewrite:
      //   dot(rayOrigin, bVector) + t * dot(rayDirection, bVector)
      //     = cConstant.
      // Solving for t, we have:
      //   t = (cConstant - dot(rayOrigin, bVector))
      //        / dot(rayDirection, bVector).
      FLOAT_TYPE denominator = dot<FLOAT_TYPE>(rayDirection, bVector);
      if(denominator == 0.0) {
        return defaultValue;
      }
      // else
      return (cConstant - dot<FLOAT_TYPE>(rayOrigin, bVector)) / denominator;
    }      

  } // namespace numeric

} // namespace brick

#endif /* #ifndef BRICK_NUMERIC_AMANATIDESWOO2D_IMPL_HH */
