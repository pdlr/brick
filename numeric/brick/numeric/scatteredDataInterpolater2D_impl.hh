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
    template <class Type, class FloatType>
    ScatteredDataInterpolater2D<Type, FloatType>::
    ScatteredDataInterpolater2D(size_t numberOfLevels, 
                                bool isMeanCentered,
                                bool isIsotropic)
      : m_bSpline2D(isIsotropic),
        m_isMeanCentered(isMeanCentered),
        m_meanValue(static_cast<FloatType>(0.0)),
        m_numberOfLevels(numberOfLevels)
    {
      if(0 == numberOfLevels) {
        BRICK_THROW(brick::common::ValueException, 
                    "ScatteredDataInterpolater::ScatteredDataInterpolater()",
                    "Argument numberOfLevels must be greater than 0.");
      }
    }

    
    // The copy constructor does a deep copy.
    template <class Type, class FloatType>
    ScatteredDataInterpolater2D<Type, FloatType>::
    ScatteredDataInterpolater2D(
      ScatteredDataInterpolater2D<Type, FloatType> const& other)
      : m_bSpline2D(other.m_bSpline2D),
        m_isMeanCentered(other.m_isMeanCentered),
        m_meanValue(other.m_meanValue),
        m_numberOfLevels(other.m_numberOfLevels)
    {
      // Empty.
    }

    
    // This function allows the spline parameters to be automatically
    // set in order to approximate an irregularly sampled function.
    template <class Type, class FloatType>
    template <class CoordIter, class ObsIter>
    void
    ScatteredDataInterpolater2D<Type, FloatType>::
    approximateScatteredData(CoordIter sBegin,
                             CoordIter sEnd,
                             CoordIter tBegin,
                             ObsIter observationsBegin,
                             FloatType buffer)
    {
      // Make a local -- and mutable -- copy of the input values.
      size_t numberOfObservations = sEnd - sBegin;
      Array1D<Type> shiftedObservations(numberOfObservations);
      std::copy(observationsBegin, observationsBegin + numberOfObservations,
                shiftedObservations.begin());

      // Start by computing the mean value, if appropriate, and
      // subtracting it out of the input data.
      if(this->m_isMeanCentered) {
        this->m_meanValue = mean(observationsBegin, observationsEnd);
        shiftedObservations -= this->m_meanValue;
      } else {
        this->m_meanValue = static_cast<FloatType>(0.0);
      }

      // Create and initial, very course, B-spline approximation to
      // the possibly-shifted input data.
      this->bSpline2D.setNumberOfControlPoints(4, 4);
      this->bSpline2D.approximateScatteredData(
        sBegin, sEnd, tBegin, shiftedObservations.begin(), buffer);
      
      // Iterate until there are no more refinements to be done.
      Array1D<Type> residuals(shiftedObservations.size());
      for(size_t levelNumber = 0; levelNumber < this->m_numberOfLevels; 
          ++levelNumber) {

        // Resample the B-spline to the next higher resolution.
        this->m_bSpline2D.promote();
        
        // Subtract the best-so-far interpolation from the observations
        // to get a residual, which will be interpolated below.
        CoordIter sIter = sBegin;
        CoordIter tIter = tBegin;
        Array1D<Type>::const_iter obsIter = shiftedObservations.begin();
        Array1D<Type>::iter residualsIter = residuals.begin();
        while(sIter != sEnd) {
          *residualsIter = *obsIter - this->bSpline2D(*sIter, *tIter);
          ++residualsIter;
          ++obsIter;
          ++sIter;
          ++tIter;
        }

        // Recover the dimensions of the current spline control grid.
        size_t numberOfControlPointsS;
        size_t numberOfControlPointsT;
        this->bSpline2D.getNumberOfControlPoints(numberOfControlPointsS,
                                                 numberOfControlPointsT);

        // Interpolate the residual at the current resolution.
        BSpline2D<Type, FloatType> residualInterpolator(this->isIsotropic);
        residualInterpolator.setNumberOfControlPoints(numberOfControlPointsS,
                                                      numberOfControlPointsT);
        residualInterpolator.approximateScatteredData(
          sBegin, sEnd, tBegin, residuals.begin(), buffer);

        // Add the residual spline to the overall interpolation
        // function, making the interpolation more accurate.
        this->m_bSpline2D += residualInterpolator;
      } // for(size_t levelNumber...)
    }


    // This member function returns the maximum values for the spline
    // parameters S and T.
    template <class Type, class FloatType>
    void
    ScatteredDataInterpolater2D<Type, FloatType>::
    getMaximumSAndTValues(FloatType& maximumS, FloatType& maximumT) const
    {
      this->m_bSpline2D.getMaximumSAndTValues(maximumS, maximumT);
    }


    // This member function returns the minimum values for the spline
    // parameters S and T.
    template <class Type, class FloatType>
    void
    ScatteredDataInterpolater2D<Type, FloatType>::
    getMinimumSAndTValues(FloatType& minimumS, FloatType& minimumT) const
    {
      this->m_bSpline2D.getMinimumSAndTValues(minimumS, minimumT);
    }


    // The assigment operator does a deep copy.
    template <class Type, class FloatType>
    ScatteredDataInterpolater2D<Type, FloatType>&
    ScatteredDataInterpolater2D<Type, FloatType>::
    operator=(ScatteredDataInterpolater2D<Type, FloatType> const& other)
    {
      if(&other != this) {
        this->m_bSpline2D = other.m_bSpline2D;
        this->m_isMeanCentered = other.m_isMeanCentered;
        this->m_meanValue = other.m_meanValue;
        this->m_numerOfLevels = other.m_numberOfLevels;
      }
    }

    
    // This operator evaluates the spline at the specified values of
    // spline parameters s and t.
    template <class Type, class FloatType>
    Type
    ScatteredDataInterpolater2D<Type, FloatType>::
    operator()(FloatType sValue, FloatType tValue) const
    {
      return this->m_bSpline2D(sValue, tValue) + this->m_meanValue;
    }


  } // namespace numeric

} // namespace brick

#endif /* #ifndef BRICK_NUMERIC_SCATTEREDDATAINTERPOLATER2D_IMPL_HH */
