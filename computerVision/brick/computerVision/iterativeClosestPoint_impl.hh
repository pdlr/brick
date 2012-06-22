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
      std::vector<Type>      queryPoints(beginIter, endIter);
      std::vector<FloatType> weights;
      std::vector<Type>      selectedQueryPoints;
      std::vector<Type>      matchingModelPoints;

      weights.reserve(            queryPoints.size());
      selectedQueryPoints.reserve(queryPoints.size());
      matchingModelPoints.reserve(queryPoints.size());

      brick::numeric::Transform3D<FloatType> modelFromQuery =
        modelFromQueryEstimate;

      bool isConverged = false;
      while(!isConverged) {
        this->selectQueryPoints(selectedQueryPoints, queryPoints);
        isConverged |= this->findMatches(
          matchingModelPoints, weights, selectedQueryPoints,
          modelFromQuery);
        modelFromQuery = this->estimateTransformModelFromQuery(
          selectedQueryPoints, matchingModelPoints, weights);
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
      std::vector<Type> const& matchingModelPoints,
      std::vector<FloatType> const& weights)
    {
      return registerPoints3D<FloatType>(selectedQueryPoints.begin(),
					 selectedQueryPoints.end(),
					 matchingModelPoints.begin(),
					 weights.begin(),
					 true);
    }


    template <unsigned int Dimension, class Type, class FloatType>
    bool
    IterativeClosestPoint<Dimension, Type, FloatType>::
    findMatches(std::vector<Type>& matchingModelPoints,
                std::vector<FloatType>& weights,
                std::vector<Type> const& queryPoints,
                brick::numeric::Transform3D<FloatType> const& modelFromQuery)
    {
      unsigned int count = 0;
      double rmsError = 0.0;

      if(matchingModelPoints.size() != queryPoints.size()) {
        matchingModelPoints.resize(queryPoints.size());
      }
      if(weights.size() != queryPoints.size()) {
        weights.resize(queryPoints.size());
      }

      for(unsigned int ii = 0; ii < queryPoints.size(); ++ii) {
        FloatType distance;
	brick::numeric::Vector3D<FloatType> transformedQueryPoint =
	  modelFromQuery * queryPoints[ii];
        matchingModelPoints[ii] = this->m_modelTree.findNearest(
          transformedQueryPoint, distance);

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
      if(rmsError < m_convergenceThreshold) {
        return true;
      }
      return false;
    }
    

    template <unsigned int Dimension, class Type, class FloatType>
    void
    IterativeClosestPoint<Dimension, Type, FloatType>::
    selectQueryPoints(std::vector<Type>& selectedQueryPoints,
                      std::vector<Type> const& allQueryPoints)
    {
      selectedQueryPoints = allQueryPoints;
    }


    
    
  } // namespace computerVision
  
} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_ITERATIVECLOSESTPOINT_IMPL_HH */
