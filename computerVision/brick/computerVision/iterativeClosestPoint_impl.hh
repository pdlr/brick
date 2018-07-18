/**
***************************************************************************
* @file brick/computerVision/iterativeClosestPoint_impl.hh
*
* Header file defining a class template implementing a derivative of
* Besl's and McKay's Iterative Closest Point algorithm.
*
* Copyright (C) 2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_ITERATIVECLOSESTPOINT_IMPL_HH
#define BRICK_COMPUTERVISION_ITERATIVECLOSESTPOINT_IMPL_HH

#include <brick/common/mathFunctions.hh>
#include <brick/computerVision/registerPoints3D.hh>
#include <brick/numeric/utilities.hh>

// This file is included by iterativeClosestPoint.hh, and should not
// be directly included by user code, so no need to include
// iterativeClosestPoint.hh here.
//
// #include <brick/computerVision/iterativeClosestPoint.hh>

namespace brick {

  namespace computerVision {

    // The default constructor.
    template <unsigned int Dimension, class Type, class FloatType>
    IterativeClosestPoint<Dimension, Type, FloatType>::
    IterativeClosestPoint()
      : m_convergenceThreshold(0.1), // TBD: set this and add better term crit.
        m_distanceThreshold(1.0),    // TBD: set this and add better criterion.
        m_iterationCount(0),
        m_modelTree()
    {
      // Empty.
    }


    // The destructor cleans up any system resources.
    template <unsigned int Dimension, class Type, class FloatType>
    IterativeClosestPoint<Dimension, Type, FloatType>::
    ~IterativeClosestPoint()
    {
      // Empty.
    }


    template <unsigned int Dimension, class Type, class FloatType>
    brick::numeric::Transform3D<FloatType>
    IterativeClosestPoint<Dimension, Type, FloatType>::
    getTransform()
    {
      brick::numeric::Transform3D<FloatType> modelFromQueryEstimate =
        brick::numeric::identity<FloatType>(Dimension, Dimension);
      return this->getTransform(modelFromQueryEstimate);
    }


    template <unsigned int Dimension, class Type, class FloatType>
    template <class Iter>
    brick::numeric::Transform3D<FloatType>
    IterativeClosestPoint<Dimension, Type, FloatType>::
    registerPoints(
      Iter beginIter, Iter endIter,
      brick::numeric::Transform3D<FloatType> const& modelFromQueryEstimate)
    {
      std::vector<Type>        queryPoints(beginIter, endIter);
      std::vector<FloatType>   weights;
      std::vector<Type>        selectedQueryPoints;
      std::vector<Type const*> matchingModelPointAddresses;

      weights.reserve(            queryPoints.size());
      selectedQueryPoints.reserve(queryPoints.size());
      matchingModelPointAddresses.reserve(queryPoints.size());

      brick::numeric::Transform3D<FloatType> modelFromQuery =
        modelFromQueryEstimate;

      m_iterationCount = 0;
      unsigned int pointCount = 0;
      FloatType rmsError(0);
      bool isConverged = false;
      while(!isConverged) {

        // Pick the starting point for our next ICP iteration.
        brick::numeric::Transform3D<FloatType> modelFromQueryHypothesis =
          modelFromQuery;

        // Set up the next ICP registration by selecting query points
        // and corresponding model points.
        bool isQuerySetChanged = this->selectQueryPoints(
          selectedQueryPoints, queryPoints);
        bool isMatchingSetChanged = this->findMatches(
          matchingModelPointAddresses, weights, pointCount, rmsError,
          selectedQueryPoints, modelFromQueryHypothesis);

        // Check for convergence.
        if(!isQuerySetChanged && !isMatchingSetChanged) {
          isConverged = true;
          continue;
        }

        // Re-estimate coordinate transformation.
        modelFromQuery = this->estimateTransformModelFromQuery(
          selectedQueryPoints, matchingModelPointAddresses, weights);

        ++m_iterationCount;
      }

      return modelFromQuery;
    }



    template <unsigned int Dimension, class Type, class FloatType>
      template <class Iter>
      void
    IterativeClosestPoint<Dimension, Type, FloatType>::
      setModelPoints(Iter beginIter, Iter endIter)
    {
      this->m_modelTree.clear();
      this->m_modelTree.addSamples(beginIter, endIter);
    }


    /* ================ Protected ================= */


    template <unsigned int Dimension, class Type, class FloatType>
    brick::numeric::Transform3D<FloatType>
    IterativeClosestPoint<Dimension, Type, FloatType>::
    estimateTransformModelFromQuery(
      std::vector<Type> const& selectedQueryPoints,
      std::vector<Type const*> const& matchingModelPointAddresses,
      std::vector<FloatType> const& weights)
    {
      std::vector<Type> matchingModelPoints(matchingModelPointAddresses.size());
      for(unsigned int ii = 0; ii < matchingModelPointAddresses.size(); ++ii) {
        matchingModelPoints[ii] = *(matchingModelPointAddresses[ii]);
      }

      return registerPoints3D<FloatType>(selectedQueryPoints.begin(),
                                         selectedQueryPoints.end(),
                                         matchingModelPoints.begin(),
                                         weights.begin(),
                                         true);
    }


    template <unsigned int Dimension, class Type, class FloatType>
    bool
    IterativeClosestPoint<Dimension, Type, FloatType>::
    findMatches(std::vector<Type const*>& matchingModelPointAddresses,
                std::vector<FloatType>& weights,
		unsigned int& count,
		FloatType& rmsError,
                std::vector<Type> const& queryPoints,
                brick::numeric::Transform3D<FloatType> const& modelFromQuery)
    {
      count = 0;
      rmsError = FloatType(0);
      bool isMatchingSetChanged = false;

      // Make sure output arrays are initialized.
      if(matchingModelPointAddresses.size() != queryPoints.size()) {
        matchingModelPointAddresses.resize(queryPoints.size());
        isMatchingSetChanged = true;
      }
      if(weights.size() != queryPoints.size()) {
        weights.resize(queryPoints.size());
        isMatchingSetChanged = true;
      }

      // Iterate over all query points.
      for(unsigned int ii = 0; ii < queryPoints.size(); ++ii) {

        // Transform the query point using our best-so-far estimate of
        // the final coordinate transformation.
        brick::numeric::Vector3D<FloatType> transformedQueryPoint =
          modelFromQuery * queryPoints[ii];

        // Find the model point that is nearest to the current
        // transformed query point.
        FloatType distance;
        Type const* matchingPointPtr = &(this->m_modelTree.findNearest(
                                           transformedQueryPoint, distance));

        // Notice if this particular point match has changed since the
        // last ICP iteration, and remember the match.
        isMatchingSetChanged |= (matchingPointPtr
                                 != matchingModelPointAddresses[ii]);
        matchingModelPointAddresses[ii] = matchingPointPtr;

        // Checks of normals, etc., go here.
        if(distance < this->m_distanceThreshold) {
          weights[ii] = 1.0;
          rmsError += distance * distance;
          ++count;
        } else {
          weights[ii] = 0.0;
        }
      }

      // std::sort(distances.begin(), distance.end());
      typedef unsigned int UInt;
      count = std::max(count, UInt(1));
      rmsError = brick::common::squareRoot(rmsError / count);

      // Tell the calling context whether this ICP iteration has any
      // new point-to-point correspondences.
      return isMatchingSetChanged;
    }


    template <unsigned int Dimension, class Type, class FloatType>
    bool
    IterativeClosestPoint<Dimension, Type, FloatType>::
    selectQueryPoints(std::vector<Type>& selectedQueryPoints,
                      std::vector<Type> const& allQueryPoints)
    {
      selectedQueryPoints = allQueryPoints;
      return false;
    }

  } // namespace computerVision

} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_ITERATIVECLOSESTPOINT_IMPL_HH */
