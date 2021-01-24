/**
***************************************************************************
* @file brick/numeric/amanatidesWoo3D_impl.hh
*
* Header file defining inline and template functions declared in
* AmanatidesWoo3D.hh.
*
* Copyright (C) 2004-2007,2012 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_AMANATIDESWOO3D_IMPL_HH
#define BRICK_NUMERIC_AMANATIDESWOO3D_IMPL_HH

// This file is included by amanatidesWoo3D.hh, and should not be
// directly included by user code, so no need to include
// amanatidesWoo3D.hh here.
//
// #include <brick/numeric/amanatidesWoo3D.hh>

namespace brick {

  namespace numeric {

    template <class ARRAY3D, class FLOAT_TYPE, class INT_TYPE>
    AmanatidesWoo3D<ARRAY3D, FLOAT_TYPE, INT_TYPE>::
    AmanatidesWoo3D(ARRAY3D& data,
                    const Transform3D<FLOAT_TYPE>& voxelTworld,
                    const Vector3D<FLOAT_TYPE>& rayOrigin,
                    const Vector3D<FLOAT_TYPE>& rayDirection,
                    bool downstreamOnly)
      : m_data(data),
        m_initialU(-1),  // Initialize to illegal values.  These will be
        m_initialV(-1),  // replaced with legal ones if there's a valid
        m_initialW(-1),  // intersection between the ray and the voxel
        // array.
        m_stepU(),
        m_stepV(),
        m_stepW(),
        m_tDeltaU(),
        m_tDeltaV(),
        m_tDeltaW(),
        m_tMaxU(),
        m_tMaxV(),
        m_tMaxW(),
        m_tStart(),
        m_validIntersection(true)
    {
      // First convert everything into voxel coordinates.
      Vector3D<FLOAT_TYPE> rayOriginVoxel = voxelTworld * rayOrigin;
      Vector3D<FLOAT_TYPE> rayDirectionVoxel =
        (voxelTworld * (rayOrigin + rayDirection)) - rayOriginVoxel;

      // Now find points entry and exit from the CT volume.  These are
      // expressed in as parameter values tEntry and tExit such that
      // (rayOriginVoxel + tEntry * rayDirectionVoxel) is the entry
      // point and (rayOriginVoxel + tExit * rayDirectionVoxel) is the
      // exit point.
      std::pair<FLOAT_TYPE, FLOAT_TYPE> tEntry_tExit =
        this->findEntryAndExitPoints(rayOriginVoxel, rayDirectionVoxel, m_data);

      // Sometimes we want to disallow any voxels which are "behind" the
      // ray origin.  That is, sometimes we want to only trace those
      // parts of the ray which correspond to positive values of the
      // parameter t.
      if(downstreamOnly && (tEntry_tExit.first < 0.0)) {
        tEntry_tExit.first = 0.0;
      }

      // Sometimes the there's no intersection with the voxel array.  In
      // this case, tExit will be less than tEntry, and we don't bother
      // doing any more calculation.
      if(tEntry_tExit.second <= tEntry_tExit.first) {
        m_validIntersection = false;
        return;
      }

      // Now make the point of entry explicit.
      Vector3D<FLOAT_TYPE> entryPoint =
        rayOriginVoxel + tEntry_tExit.first * rayDirectionVoxel;

      // Correct for rounding error which sometimes makes entryPoint
      // have values like -1.0e-14.
      if(entryPoint.x() < 0.0) {
        entryPoint.setX(0.0);
      }
      if(entryPoint.y() < 0.0) {
        entryPoint.setY(0.0);
      }
      if(entryPoint.z() < 0.0) {
        entryPoint.setZ(0.0);
      }

      // Finally, assign the variables described in the Amanatides' and
      // Woo's paper.

      // Since we've already converted to voxel coords, finding the
      // first voxel on the path is trivial.
      m_initialU = static_cast<INT_TYPE>(entryPoint.x());
      m_initialV = static_cast<INT_TYPE>(entryPoint.y());
      m_initialW = static_cast<INT_TYPE>(entryPoint.z());

      // There's an similar case to the one handled by the rounding
      // error test above.  Member function findEntryAndExitPoints()
      // uses m_data.shape to define the maximum U, V, and W
      // coordinates, respectively.  This means that it's very possible
      // for m_initialU to be equal to m_data.columns(), which will
      // cause an out-of-bounds index. Correct for this now.
      if(m_initialU == static_cast<INT_TYPE>(m_data.shape()[2])) {
        --m_initialU;
      }
      if(m_initialV == static_cast<INT_TYPE>(m_data.shape()[1])) {
        --m_initialV;
      }
      if(m_initialW == static_cast<INT_TYPE>(m_data.shape()[0])) {
        --m_initialW;
      }

      // Sanity check.  This sometimes fails if the ray is parallel to
      // one axis of the data array.
      if((m_initialU >= static_cast<INT_TYPE>(m_data.shape()[2]))
         || (m_initialV >= static_cast<INT_TYPE>(m_data.shape()[1]))
         || (m_initialW >= static_cast<INT_TYPE>(m_data.shape()[0]))) {
        m_validIntersection = false;
        return;
      }

      // m_tStart is just the same as tEntry.
      m_tStart = tEntry_tExit.first;

      // The remaining variables depend on whether U & V will be
      // increasing or decreasing as we travel along the ray, so we need
      // if clauses.  Please see the declaration for documentation on
      // what each of these member variables means.
      if(rayDirectionVoxel.x() > 0.0) {
        m_stepU = 1;
        m_tDeltaU = 1.0 / rayDirectionVoxel.x();
        m_tMaxU = m_tStart + (((m_initialU + 1) - entryPoint.x())
                              / rayDirectionVoxel.x());
      } else if(rayDirectionVoxel.x() < 0.0) {
        m_stepU = -1;
        m_tDeltaU = -(1.0 / rayDirectionVoxel.x());
        m_tMaxU = m_tStart + ((m_initialU - entryPoint.x())
                              / rayDirectionVoxel.x());
      } else { // rayDirectionVoxel.x() == 0.0;
        m_stepU = 0;
        m_tDeltaU = std::numeric_limits<FLOAT_TYPE>::max();
        m_tMaxU = std::numeric_limits<FLOAT_TYPE>::max();
      }
      if(rayDirectionVoxel.y() > 0.0) {
        m_stepV = 1;
        m_tDeltaV = 1.0 / rayDirectionVoxel.y();
        m_tMaxV = m_tStart + (((m_initialV + 1) - entryPoint.y())
                              / rayDirectionVoxel.y());
      } else if(rayDirectionVoxel.y() < 0.0) {
        m_stepV = -1;
        m_tDeltaV = -(1.0 / rayDirectionVoxel.y());
        m_tMaxV = m_tStart + ((m_initialV - entryPoint.y())
                              / rayDirectionVoxel.y());
      } else { // rayDirectionVoxel.y() == 0.0;
        m_stepV = 0;
        m_tDeltaV = std::numeric_limits<FLOAT_TYPE>::max();
        m_tMaxV = std::numeric_limits<FLOAT_TYPE>::max();
      }
      if(rayDirectionVoxel.z() > 0.0) {
        m_stepW = 1;
        m_tDeltaW = 1.0 / rayDirectionVoxel.z();
        m_tMaxW = m_tStart + (((m_initialW + 1) - entryPoint.z())
                              / rayDirectionVoxel.z());
      } else if(rayDirectionVoxel.z() < 0.0) {
        m_stepW = -1;
        m_tDeltaW = -(1.0 / rayDirectionVoxel.z());
        m_tMaxW = m_tStart + ((m_initialW - entryPoint.z())
                              / rayDirectionVoxel.z());
      } else { // rayDirectionVoxel.z() == 0.0;
        m_stepW = 0;
        m_tDeltaW = std::numeric_limits<FLOAT_TYPE>::max();
        m_tMaxW = std::numeric_limits<FLOAT_TYPE>::max();
      }
    }

    template <class ARRAY3D, class FLOAT_TYPE, class INT_TYPE>
    AmanatidesWoo3D<ARRAY3D, FLOAT_TYPE, INT_TYPE>::
    AmanatidesWoo3D(const AmanatidesWoo3D& source)
      : m_data(source.m_data),
        m_initialU(source.m_initialU),
        m_initialV(source.m_initialV),
        m_initialW(source.m_initialW),
        m_stepU(source.m_stepU),
        m_stepV(source.m_stepV),
        m_stepW(source.m_stepW),
        m_tDeltaU(source.m_tDeltaU),
        m_tDeltaV(source.m_tDeltaV),
        m_tDeltaW(source.m_tDeltaW),
        m_tMaxU(source.m_tMaxU),
        m_tMaxV(source.m_tMaxV),
        m_tMaxW(source.m_tMaxW),
        m_tStart(source.m_tStart),
        m_validIntersection(source.m_validIntersection)
    {
      // Empty
    }

    template <class ARRAY3D, class FLOAT_TYPE, class INT_TYPE>
    AmanatidesWoo3D<ARRAY3D, FLOAT_TYPE, INT_TYPE>::
    ~AmanatidesWoo3D() {}

    template <class ARRAY3D, class FLOAT_TYPE, class INT_TYPE>
    typename AmanatidesWoo3D<ARRAY3D, FLOAT_TYPE, INT_TYPE>::iterator
    AmanatidesWoo3D<ARRAY3D, FLOAT_TYPE, INT_TYPE>::
    begin()
    {
      return AmanatidesWoo3DIterator<ARRAY3D, FLOAT_TYPE, INT_TYPE>(
        m_data, m_initialU, m_initialV, m_initialW, m_stepU, m_stepV, m_stepW,
        m_tMaxU, m_tMaxV, m_tMaxW, m_tDeltaU, m_tDeltaV, m_tDeltaW, m_tStart);
    }

    template <class ARRAY3D, class FLOAT_TYPE, class INT_TYPE>
    typename AmanatidesWoo3D<ARRAY3D, FLOAT_TYPE, INT_TYPE>::iterator
    AmanatidesWoo3D<ARRAY3D, FLOAT_TYPE, INT_TYPE>::
    end()
    {
      // AmanatidesWoo3DIterator<ARRAY3D, FLOAT_TYPE, INT_TYPE>::operator==(...) considers all
      // iterators with illegal voxel coordinates to be equal, so all we
      // have to do here is return an iterator which references an
      // illegal voxel.
      return AmanatidesWoo3DIterator<ARRAY3D, FLOAT_TYPE, INT_TYPE>(
        m_data, -1, -1, -1, m_stepU, m_stepV, m_stepW,
        m_tMaxU, m_tMaxV, m_tMaxW, m_tDeltaU, m_tDeltaV, m_tDeltaW, m_tStart);
    }

    template <class ARRAY3D, class FLOAT_TYPE, class INT_TYPE>
    inline bool
    AmanatidesWoo3D<ARRAY3D, FLOAT_TYPE, INT_TYPE>::
    validIntersection()
    {
      return m_validIntersection;
    }

    template <class ARRAY3D, class FLOAT_TYPE, class INT_TYPE>
    AmanatidesWoo3D<ARRAY3D, FLOAT_TYPE, INT_TYPE>&
    AmanatidesWoo3D<ARRAY3D, FLOAT_TYPE, INT_TYPE>::
    operator=(const AmanatidesWoo3D& source)
    {
      m_data = source.m_data;
      m_initialU = source.m_initialU;
      m_initialV = source.m_initialV;
      m_initialW = source.m_initialW;
      m_stepU = source.m_stepU;
      m_stepV = source.m_stepV;
      m_stepW = source.m_stepW;
      m_tDeltaU = source.m_tDeltaU;
      m_tDeltaV = source.m_tDeltaV;
      m_tDeltaW = source.m_tDeltaW;
      m_tMaxU = source.m_tMaxU;
      m_tMaxV = source.m_tMaxV;
      m_tMaxW = source.m_tMaxW;
      m_tStart = source.m_tStart;
      m_validIntersection = source.m_validIntersection;
      return *this;
    }

    template <class ARRAY3D, class FLOAT_TYPE, class INT_TYPE>
    std::pair<FLOAT_TYPE, FLOAT_TYPE>
    AmanatidesWoo3D<ARRAY3D, FLOAT_TYPE, INT_TYPE>::
    findEntryAndExitPoints(const Vector3D<FLOAT_TYPE>& rayOriginVoxel,
                           const Vector3D<FLOAT_TYPE>& rayDirectionVoxel,
                           const ARRAY3D& data)
    {
      // First find intersection with each boundary line.

      // ... Find the intersection with the plane U = 0, or else a
      // really small number if rayDirection is parallel to the U axis.
      FLOAT_TYPE tIntersectU0 = findIntersection(
        rayOriginVoxel, rayDirectionVoxel, Vector3D<FLOAT_TYPE>(1.0, 0.0, 0.0),
        0.0, -std::numeric_limits<FLOAT_TYPE>::max());

      // ... Find the intersection with the plane U = data.shape()[2],
      // or else a really big number if rayDirection is parallel to the
      // U axis.
      FLOAT_TYPE tIntersectU1 = findIntersection(
        rayOriginVoxel, rayDirectionVoxel, Vector3D<FLOAT_TYPE>(1.0, 0.0, 0.0),
        data.shape()[2], std::numeric_limits<FLOAT_TYPE>::max());

      // ... Find the intersection with the plane V = 0, or else a
      // really small number if rayDirection is parallel to the V axis.
      FLOAT_TYPE tIntersectV0 = findIntersection(
        rayOriginVoxel, rayDirectionVoxel, Vector3D<FLOAT_TYPE>(0.0, 1.0, 0.0),
        0.0, -std::numeric_limits<FLOAT_TYPE>::max());

      // ... Find the intersection with the plane V = data.shape()[1],
      // or else a really big number if rayDirection is parallel to the
      // V axis.
      FLOAT_TYPE tIntersectV1 = findIntersection(
        rayOriginVoxel, rayDirectionVoxel, Vector3D<FLOAT_TYPE>(0.0, 1.0, 0.0),
        data.shape()[1], std::numeric_limits<FLOAT_TYPE>::max());

      // ... Find the intersection with the plane W = 0, or else a
      // really small number if rayDirection is parallel to the W axis.
      FLOAT_TYPE tIntersectW0 = findIntersection(
        rayOriginVoxel, rayDirectionVoxel, Vector3D<FLOAT_TYPE>(0.0, 0.0, 1.0),
        0.0, -std::numeric_limits<FLOAT_TYPE>::max());

      // ... Find the intersection with the plane W = data.shape()[0],
      // or else a really big number if rayDirection is parallel to the
      // W axis.
      FLOAT_TYPE tIntersectW1 = findIntersection(
        rayOriginVoxel, rayDirectionVoxel, Vector3D<FLOAT_TYPE>(0.0, 0.0, 1.0),
        data.shape()[0], std::numeric_limits<FLOAT_TYPE>::max());

      // Now find the closer and farther of each pair of intersections.
      FLOAT_TYPE tMinU = std::min(tIntersectU0, tIntersectU1);
      FLOAT_TYPE tMinV = std::min(tIntersectV0, tIntersectV1);
      FLOAT_TYPE tMinW = std::min(tIntersectW0, tIntersectW1);
      FLOAT_TYPE tMaxU = std::max(tIntersectU0, tIntersectU1);
      FLOAT_TYPE tMaxV = std::max(tIntersectV0, tIntersectV1);
      FLOAT_TYPE tMaxW = std::max(tIntersectW0, tIntersectW1);

      // Compute closest point which could possibly intersect with the volume.
      FLOAT_TYPE tEntry = std::max(tMinU, std::max(tMinV, tMinW));
      // Compute farthest point which could possibly intersect with the volume.
      FLOAT_TYPE tExit = std::min(tMaxU, std::min(tMaxV, tMaxW));

      // Return our findings.
      return std::make_pair(tEntry, tExit);
    }

    template <class ARRAY3D, class FLOAT_TYPE, class INT_TYPE>
    FLOAT_TYPE
    AmanatidesWoo3D<ARRAY3D, FLOAT_TYPE, INT_TYPE>::
    findIntersection(const Vector3D<FLOAT_TYPE>& rayOrigin,
                     const Vector3D<FLOAT_TYPE>& rayDirection,
                     const Vector3D<FLOAT_TYPE>& bVector,
                     FLOAT_TYPE cConstant,
                     FLOAT_TYPE defaultValue)
    {
      // bVector and cConstant describe the desired plane:
      //   dot(x, bVector) = cConstant.
      // Now, x is constrained to lie on the specified line:
      //   x = rayOrigin + alpha * rayDirection.
      // Substituting, we have:
      //   dot(rayOrigin + alpha * rayDirection, bVector) = cConstant.
      // Which we rewrite:
      //   dot(rayOrigin, bVector) + alpha * dot(rayDirection, bVector)
      //     = cConstant.
      // Solving for alpha, we have:
      //   alpha = (cConstant - dot(rayOrigin, bVector))
      //            / dot(rayDirection, bVector).
      FLOAT_TYPE denominator = dot<FLOAT_TYPE>(rayDirection, bVector);
      if(denominator == 0.0) {
        return defaultValue;
      }
      // else
      return (cConstant - dot<FLOAT_TYPE>(rayOrigin, bVector)) / denominator;
    }

  } // namespace numeric

} // namespace brick

#endif /* #ifndef BRICK_NUMERIC_AMANATIDESWOO3D_IMPL_HH */
