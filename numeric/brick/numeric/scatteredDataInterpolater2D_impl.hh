/**
***************************************************************************
* @file brick/numeric/scatteredDataInterpolater2D_impl.hh
*
* Header file defining inline and template functions from ScatteredDataInterpolater2D.hh.
*
* Copyright (C) 2014 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_NUMERIC_SCATTEREDDATAINTERPOLATER2D_IMPL_HH
#define BRICK_NUMERIC_SCATTEREDDATAINTERPOLATER2D_IMPL_HH

// This file is included by scatteredDataInterpolater2D.hh, and should
// not be directly included by user code, so no need to include
// scatteredDataInterpolater2D.hh here.
// 
// #include <brick/numeric/scatteredDataInterpolater2D.hh>

#include <cmath>
#include <algorithm>
#include <brick/numeric/functional.hh>
#include <brick/numeric/utilities.hh>

namespace brick {

  namespace numeric {

    // This constructor builds a ScatteredDataInterpolater2D instance
    // of unspecified length and width.
    template <class Type, class FloatType, class TestType>
    ScatteredDataInterpolater2D<Type, FloatType, TestType>::
    ScatteredDataInterpolater2D(size_t numberOfLevels, 
                                bool isMeanCentered,
                                bool isIsotropic)
      : m_bSpline2D(isIsotropic),
        m_isMeanCentered(isMeanCentered),
        m_meanValue(static_cast<Type>(0.0)),
        m_numberOfLevels(numberOfLevels),
        m_testFunctor()
    {
      if(0 == numberOfLevels) {
        BRICK_THROW(brick::common::ValueException, 
                    "ScatteredDataInterpolater::ScatteredDataInterpolater()",
                    "Argument numberOfLevels must be greater than 0.");
      }
    }

    
    // The copy constructor does a deep copy.
    template <class Type, class FloatType, class TestType>
    ScatteredDataInterpolater2D<Type, FloatType, TestType>::
    ScatteredDataInterpolater2D(
      ScatteredDataInterpolater2D<Type, FloatType, TestType> const& other)
      : m_bSpline2D(other.m_bSpline2D),
        m_isMeanCentered(other.m_isMeanCentered),
        m_meanValue(other.m_meanValue),
        m_numberOfLevels(other.m_numberOfLevels),
        m_testFunctor(other.m_testFunctor)
    {
      // Empty.
    }

    
    // This function allows the spline parameters to be automatically
    // set in order to approximate an irregularly sampled function.
    template <class Type, class FloatType, class TestType>
    template <class CoordIter, class ObsIter>
    void
    ScatteredDataInterpolater2D<Type, FloatType, TestType>::
    approximate(CoordIter sBegin, CoordIter sEnd,
                CoordIter tBegin,
                ObsIter observationsBegin,
                FloatType buffer)
    {
      CoordIter tEnd = tBegin + (sEnd - sBegin);
    
      // Establish bounds for reconstruction
      Vector2D<FloatType> corner0(*std::min_element(sBegin, sEnd) - buffer,
                                  *std::min_element(tBegin, tEnd) - buffer);
      Vector2D<FloatType> corner1(*std::max_element(sBegin, sEnd) + buffer,
                                  *std::max_element(tBegin, tEnd) + buffer);

      // Do the interpolation.
      this->approximate(sBegin, sEnd, tBegin, observationsBegin,
                        corner0, corner1);
    }


    // This function allows the spline parameters to be automatically
    // set in order to approximate an irregularly sampled function.
    template <class Type, class FloatType, class TestType>
    template <class CoordIter, class ObsIter>
    void
    ScatteredDataInterpolater2D<Type, FloatType, TestType>::
    approximate(CoordIter sBegin, CoordIter sEnd,
                CoordIter tBegin,
                ObsIter observationsBegin,
                Vector2D<FloatType> const& corner0,
                Vector2D<FloatType> const& corner1)
    {
      // Make a local -- and mutable -- copy of the input values.
      size_t numberOfObservations = sEnd - sBegin;
      Array1D<Type> shiftedObservations(numberOfObservations);
      std::copy(observationsBegin, observationsBegin + numberOfObservations,
                shiftedObservations.begin());

      // Start by computing the mean value, if appropriate, and
      // subtracting it out of the input data.
      if(this->m_isMeanCentered) {
        this->m_meanValue = mean<Type>(shiftedObservations);
        shiftedObservations -= this->m_meanValue;
      } else {
        this->m_meanValue = static_cast<Type>(0.0);
      }

      // Create and initial, very course, B-spline approximation to
      // the possibly-shifted input data.
      this->m_bSpline2D.setNumberOfNodes(4, 4);
      this->m_bSpline2D.approximateScatteredData(
        sBegin, sEnd, tBegin, shiftedObservations.begin(), corner0, corner1);
      
      // Iterate until there are no more refinements to be done.
      Array1D<Type> residuals(shiftedObservations.size());
      for(size_t levelNumber = 1; levelNumber < this->m_numberOfLevels; 
          ++levelNumber) {

        // Subtract the best-so-far interpolation from the observations
        // to get a residual, which will be interpolated below.
        CoordIter sIter = sBegin;
        CoordIter tIter = tBegin;
        typename Array1D<Type>::const_iterator obsIter =
          shiftedObservations.begin();
        typename Array1D<Type>::iterator residualsIter =
          residuals.begin();

        bool isTerminationOk = true;
        while(sIter != sEnd) {
          *residualsIter = *obsIter - this->m_bSpline2D(*sIter, *tIter);
          isTerminationOk = (isTerminationOk
                             && this->m_testFunctor(*residualsIter));
          ++residualsIter;
          ++obsIter;
          ++sIter;
          ++tIter;
        }

        // If all residuals are acceptable, then there's no sense in
        // continuing to iterate.
        if(isTerminationOk) {
          break;
        }
        
        // Resample the B-spline to the next higher resolution.
        this->m_bSpline2D.promote();
        
        // Recover the dimensions of the current spline control grid.
        size_t numberOfNodesS;
        size_t numberOfNodesT;
        this->m_bSpline2D.getNumberOfNodes(numberOfNodesS,
                                           numberOfNodesT);

        // Interpolate the residual at the current resolution.
        BSpline2D<Type, FloatType> residualInterpolator(
          this->m_bSpline2D.getIsIsotropic());
        residualInterpolator.setNumberOfNodes(numberOfNodesS,
                                              numberOfNodesT);
        residualInterpolator.approximateScatteredData(
          sBegin, sEnd, tBegin, residuals.begin(), corner0, corner1);

        // Add the residual spline to the overall interpolation
        // function, making the interpolation more accurate.
        this->m_bSpline2D += residualInterpolator;
      } // for(size_t levelNumber...)
    }


    // This member function returns the maximum values for the spline
    // parameters S and T.
    template <class Type, class FloatType, class TestType>
    void
    ScatteredDataInterpolater2D<Type, FloatType, TestType>::
    getMaximumSAndTValues(FloatType& maximumS, FloatType& maximumT) const
    {
      this->m_bSpline2D.getMaximumSAndTValues(maximumS, maximumT);
    }


    // This member function returns the minimum values for the spline
    // parameters S and T.
    template <class Type, class FloatType, class TestType>
    void
    ScatteredDataInterpolater2D<Type, FloatType, TestType>::
    getMinimumSAndTValues(FloatType& minimumS, FloatType& minimumT) const
    {
      this->m_bSpline2D.getMinimumSAndTValues(minimumS, minimumT);
    }


    // The assigment operator does a deep copy.
    template <class Type, class FloatType, class TestType>
    ScatteredDataInterpolater2D<Type, FloatType, TestType>&
    ScatteredDataInterpolater2D<Type, FloatType, TestType>::
    operator=(ScatteredDataInterpolater2D<Type, FloatType, TestType> const& other)
    {
      if(&other != this) {
        this->m_bSpline2D = other.m_bSpline2D;
        this->m_isMeanCentered = other.m_isMeanCentered;
        this->m_meanValue = other.m_meanValue;
        this->m_numerOfLevels = other.m_numberOfLevels;
        this->m_testFunctor = other.m_testFunctor;
      }
    }

    
    // This operator evaluates the spline at the specified values of
    // spline parameters s and t.
    template <class Type, class FloatType, class TestType>
    Type
    ScatteredDataInterpolater2D<Type, FloatType, TestType>::
    operator()(FloatType sValue, FloatType tValue) const
    {
      return this->m_bSpline2D(sValue, tValue) + this->m_meanValue;
    }


  } // namespace numeric

} // namespace brick

#endif /* #ifndef BRICK_NUMERIC_SCATTEREDDATAINTERPOLATER2D_IMPL_HH */
