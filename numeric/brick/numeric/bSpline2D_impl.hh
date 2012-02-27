/**
***************************************************************************
* @file brick/numeric/bSpline2D_impl.hh
*
* Header file defining inline and template functions from BSpline2D.hh.
*
* Copyright (C) 2006-2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_NUMERIC_BSPLINE2D_IMPL_HH
#define BRICK_NUMERIC_BSPLINE2D_IMPL_HH

// This file is included by bSpline2D.hh, and should not be directly included
// by user code, so no need to include bSpline2D.hh here.
// 
// #include <brick/numeric/bSpline2D.hh>

#include <cmath>
#include <algorithm>
#include <brick/numeric/functional.hh>
#include <brick/numeric/utilities.hh>

namespace brick {

  namespace numeric {

    // This constructor builds a BSpline2D instance of unspecified length.
    template <class Type>
    BSpline2D<Type>::
    BSpline2D()
      : m_basisArray(4),
        m_controlGrid(),
        m_minimumXY(0.0, 0.0),
        m_maximumXY(1.0, 1.0),
        m_numberOfNodesS(0),
        m_numberOfNodesT(0),
        m_xyCellOrigin(),
        m_xyCellSize()
    {
      // Temporary storage for polynomial coefficients.
      Array1D<double> basisCoefficients(4);

      // // The following basis functions are copied from Lee,
      // // Wolberg, and Shin, "Scattered Data Interpolation with Multilevel B-Splines, IEEE Transactions on Visualization and Computer Graphics, Vol 3, 228-244, 1997.  We present them here without further
      
      // // justification, although we suspect you could them recursively
      // // using a 2D version of the procedure used in
      // // BSpline::computeBasisFunction() from file bSpline.h.
      
      // First cubic spline basis component is
      // B(s) = (1 - s)**3 / 6.
      // This expands to
      // B(s) = -(1/6)s**3 + (1/2)s**2 - (1/2)s + 1/6.
      basisCoefficients(0) = 1.0 / 6.0;
      basisCoefficients(1) = -0.5;
      basisCoefficients(2) = 0.5;
      basisCoefficients(3) = -1.0 / 6.0;
      m_basisArray(0) = basisCoefficients.copy();

      // Second cubic spline basis component is
      // B(s) = (1/2)t**3 - t**2 + 2/3.
      basisCoefficients(0) = 2.0 / 3.0;
      basisCoefficients(1) = 0.0;
      basisCoefficients(2) = -1.0;
      basisCoefficients(3) = 0.5;
      m_basisArray(1) = basisCoefficients.copy();

      // Third cubic spline basis component is
      // B(s) = -(1/2)t**3 + (1/2)t**2 + (1/2)t + 1/6.
      basisCoefficients(0) = 1.0 / 6.0;
      basisCoefficients(1) = 0.5;
      basisCoefficients(2) = 0.5;
      basisCoefficients(3) = -0.5;
      m_basisArray(2) = basisCoefficients.copy();

      // Fourth cubic spline basis component is
      // B(s) = (1/6)t**3.
      basisCoefficients(0) = 0.0;
      basisCoefficients(1) = 0.0;
      basisCoefficients(2) = 0.0;
      basisCoefficients(3) = 1.0 / 6.0;
      m_basisArray(3) = basisCoefficients.copy();

      // // What we would do if we had data...
      // m_xyCellSize.setValue((m_maximumXY.x() - m_minimumXY.x())
      //                       / m_numberOfNodesS,
      //                       (m_maximumXY.y() - m_minimumXY.y())
      //                       / m_numberOfNodesT);
      // m_xyCellOrigin = m_minimumXY - m_xyCellSize;
      // m_controlGrid.reinit(
      //   std::ceil((m_maximumXY.x() - m_minimumXY.x()) / m_xyCellSize.x()) + 3,
      //   std::ceil((m_maximumXY.y() - m_minimumXY.y()) / m_xyCellSize.y()) + 3);
      // m_controlGrid = 0.0;

      // What we would do instead, because we don't have data yet.
      m_xyCellSize.setValue(1.0, 1.0);
      m_xyCellOrigin.setValue(0.0, 0.0);
      m_controlGrid.reinit(3, 3);
      m_controlGrid = 0.0;
    }

    
    // The copy constructor does a deep copy.
    template <class Type>
    BSpline2D<Type>::
    BSpline2D(const BSpline2D& other)
      : m_basisArray(4),
        m_controlGrid(other.m_controlGrid.copy()),
        m_minimumXY(other.m_minimumXY),
        m_maximumXY(other.m_maximumXY),
        m_numberOfNodesS(other.m_numberOfNodesS),
        m_numberOfNodesT(other.m_numberOfNodesT),
        m_xyCellOrigin(other.m_xyCellOrigin),
        m_xyCellSize(other.m_xyCellSize)
    {
      // Deep copy basis coefficients.
      for(size_t index0 = 0; index0 < m_basisArray.size(); ++index0) {
        m_basisArray[index0] = (other.m_basisArray[index0]).copy();
      }
    }

    
    // This function allows the spline parameters to be automatically
    // set in order to approximate an irregularly sampled function.
    template <class Type>
    template <class CoordIter, class ObsIter>
    void
    BSpline2D<Type>::
    approximateScatteredData(CoordIter sBegin,
                             CoordIter sEnd,
                             CoordIter tBegin,
                             ObsIter observationsBegin,
                             double buffer)
    {
      CoordIter tEnd = tBegin + (sEnd - sBegin);
      
      // Set bounds for reconstruction
      m_minimumXY.setValue(*std::min_element(sBegin, sEnd) - buffer,
                           *std::min_element(tBegin, tEnd) - buffer);
      m_maximumXY.setValue(*std::max_element(sBegin, sEnd) + buffer,
                           *std::max_element(tBegin, tEnd) + buffer);
      
      // Chose the origin and spacing of the control grid.
      m_xyCellSize.setValue(
        (m_maximumXY.x() - m_minimumXY.x()) / (m_numberOfNodesS - 3),
        (m_maximumXY.y() - m_minimumXY.y()) / (m_numberOfNodesT - 3));

      // Note(xxx): Here is where we diverge from the reference
      // implementation.
      // double cellSize = std::min(
      //   (m_maximumXY.x() - m_minimumXY.x()) / (m_numberOfNodesS - 3),
      //   (m_maximumXY.y() - m_minimumXY.y()) / (m_numberOfNodesT - 3));
      // m_xyCellSize.setValue(cellSize, cellSize);
      m_xyCellOrigin = m_minimumXY - m_xyCellSize;
      
      // Allocate some space for intermediate grids, as described in
      // Lee, Wolberg, and Shin.
      Array2D<Type> deltaGrid(m_controlGrid.rows(), m_controlGrid.columns());
      Array2D<double> omegaGrid(m_controlGrid.rows(), m_controlGrid.columns());
      deltaGrid = static_cast<Type>(0.0);
      omegaGrid = 0.0;
    
      // This code implements the algorithm on page 231 of the paper.
      size_t iIndex;
      size_t jIndex;
      Array2D<double> wArray(4, 4);
      double powersOfS[4];
      double powersOfT[4];
      while(sBegin != sEnd) {
        double w2Sum = 0.0;
        // This call sets the value of powersOfS and powersOfT,
        // and returns by reference the indices of the control grid
        // cell into which (s, t) falls.
        this->decomposeSamplePoint(*sBegin, *tBegin, iIndex, jIndex,
                                   powersOfS, powersOfT);

        // Now on with Lee, Wolberg, and Shin's algorithm.
        for(size_t kIndex = 0; kIndex < 4; ++kIndex) {
          double B_k = std::inner_product(
            powersOfS, powersOfS + 4, (this->m_basisArray)[kIndex].data(),
            static_cast<double>(0));
          for(size_t lIndex = 0; lIndex < 4; ++lIndex) {
            double B_l = std::inner_product(
              powersOfT, powersOfT + 4, (this->m_basisArray)[lIndex].data(),
              static_cast<double>(0));
            double wValue = B_k * B_l;
            wArray(kIndex, lIndex) = wValue;
            w2Sum += wValue * wValue;            
          }
        }
        for(size_t kIndex = 0; kIndex < 4; ++kIndex) {
          for(size_t lIndex = 0; lIndex < 4; ++lIndex) {
            size_t index0 = iIndex + kIndex - 1;
            size_t index1 = jIndex + lIndex - 1;
            Type phi = (wArray(kIndex, lIndex) / w2Sum) * (*observationsBegin);
            double wValue = wArray(kIndex, lIndex);

            // xxx
            if(index0 >= deltaGrid.rows() || index1 >= deltaGrid.columns()) {
              int* foo = 0;
              *foo = 0;
            }
            std::cout << "(" << index0 << ", " << index1 << ") vs ("
                      << deltaGrid.rows() << ", " << deltaGrid.columns()
                      << ")" << std::endl;
            deltaGrid(index0, index1) += phi * wValue * wValue;
            omegaGrid(index0, index1) += wValue * wValue;
          }
        }
        ++sBegin;
        ++tBegin;
        ++observationsBegin;
      }

      // Final averaging step in case neighboring input points want
      // different control grid values.
      for(size_t index0 = 0; index0 < m_controlGrid.size(); ++index0) {
        if(omegaGrid[index0] == 0.0) {
          m_controlGrid[index0] = static_cast<Type>(0.0);
        } else {
          m_controlGrid[index0] = deltaGrid[index0] / omegaGrid[index0];
        }
      }
    }
    

    // This member function returns the maximum values for the spline
    // parameters S and T.
    template <class Type>
    void
    BSpline2D<Type>::
    getMaximumSAndTValues(double& maximumS, double& maximumT)
    {
      maximumS = m_maximumXY.x();
      maximumT = m_maximumXY.y();
    }


    // This member function returns the minimum values for the spline
    // parameters S and T.
    template <class Type>
    void
    BSpline2D<Type>::
    getMinimumSAndTValues(double& minimumS, double& minimumT)
    {
      minimumS = m_minimumXY.x();
      minimumT = m_minimumXY.y();
    }


    // This member function sets the values of the control points of
    // the spline.  If the spline is periodic, then the value of the
    template <class Type>
    void
    BSpline2D<Type>::
    setControlPoints(const Array2D<Type>& controlPoints)
    {
      BRICK_THROW(brick::common::NotImplementedException,
                  "BSpline2D::setControlPoints()",
                  "Not yet written...");
      if(controlPoints.rows() <= 3 || controlPoints.columns() <= 3) {
        BRICK_THROW(brick::common::ValueException,
                    "BSpline2D::setControlPoints()",
                    "For a bicubic spline, control points array must be at "
                    "least 4x4.");
      }
      m_controlGrid = controlPoints.copy();
      // Error(xxx): m_xyCellOrigin = ??;
      m_minimumXY.setValue(0.0, 0.0);
      m_maximumXY.setValue(static_cast<double>(controlPoints.columns() - 3),
                           static_cast<double>(controlPoints.rows() - 3));
    }

    
    // This member function both specifies the number of nodes in
    // the spline and sets the node positions so that the spline is
    // "uniform".
    template <class Type>
    void
    BSpline2D<Type>::
    setNumberOfNodes(size_t numberOfNodesS,
                     size_t numberOfNodesT)
    {
      if(numberOfNodesS <= 3 || numberOfNodesT <= 3) {
        BRICK_THROW(brick::common::ValueException,
                    "BSpline2D::setNumberOfNodes()",
                    "Both arguments must be greater than 3.");
      }
      m_numberOfNodesS = numberOfNodesS;
      m_numberOfNodesT = numberOfNodesT;
      m_controlGrid.reinit(m_numberOfNodesS, m_numberOfNodesT);
      m_controlGrid = 0.0;
    }

    
    // The assigment operator does a deep copy.
    template <class Type>
    BSpline2D<Type>&
    BSpline2D<Type>::
    operator=(const BSpline2D<Type>& other)
    {
      if(&other != this) {
        m_controlGrid = other.m_controlGrid.copy();
        m_minimumXY = other.m_minimumXY;
        m_maximumXY = other.m_maximumXY;
        m_xyCellOrigin = other.m_xyCellOrigin;
        m_xyCellSize = other.m_xyCellSize;
      }
    }

    
    // This operator evaluates the spline at the specified value of
    // spline parameter s.
    template <class Type>
    Type
    BSpline2D<Type>::
    operator()(double sValue, double tValue)
    {
      // This call sets the value of powersOfS and powersOfT,
      // and returns by reference the indices of the control grid
      // cell into which (s, t) falls.
      size_t iIndex;
      size_t jIndex;
      double powersOfS[4];
      double powersOfT[4];
      this->decomposeSamplePoint(sValue, tValue, iIndex, jIndex,
                                 powersOfS, powersOfT);
      
      // Interpolate by adding spline basis functions from the
      // surrounding control points.
      int index0 = iIndex - 1;
      int index1 = jIndex - 1;
      double functionValue = 0.0;
      for(size_t kIndex = 0; kIndex < 4; ++kIndex) {
        size_t i0PlusK = index0 + kIndex;
        // double B_k = dot<double>((this->m_basisArray)[kIndex], powersOfS);
        double B_k = std::inner_product(
          powersOfS, powersOfS + 4, (this->m_basisArray)[kIndex].data(),
          static_cast<double>(0));
        for(size_t lIndex = 0; lIndex < 4; ++lIndex) {
          size_t i1PlusL = index1 + lIndex;
          // double B_l = dot<double>((this->m_basisArray)[lIndex], powersOfT);
          double B_l = std::inner_product(
            powersOfT, powersOfT + 4, (this->m_basisArray)[lIndex].data(),
            static_cast<double>(0));
          functionValue += (B_k * B_l * m_controlGrid(i0PlusK, i1PlusL));
        }
      }
      return functionValue;
    }


    template <class Type>
    void
    BSpline2D<Type>::
    decomposeSamplePoint(double sValue, double tValue,
                         size_t& iIndex, size_t& jIndex,
                         double* powersOfS, double* powersOfT)
    {
      // Note(xxx): consider using std::modf() here.
      
      // Find the integer coords of the control grid cell in which the
      // input point lies.
      double iTmp = (sValue - m_xyCellOrigin.x()) / m_xyCellSize.x();
      double jTmp = (tValue - m_xyCellOrigin.y()) / m_xyCellSize.y();
      int iCoord = static_cast<int>(std::floor(iTmp));
      int jCoord = static_cast<int>(std::floor(jTmp));

      // Find real valued coords within the cell, along with all of
      // the powers of those coords we'll be wanting to plug into
      // spline basis functions.
      powersOfS[0] = 1.0;
      powersOfS[1] = iTmp - iCoord;
      powersOfS[2] = powersOfS[1] * powersOfS[1];
      powersOfS[3] = powersOfS[2] * powersOfS[1];
      powersOfT[0] = 1.0;
      powersOfT[1] = jTmp - jCoord;
      powersOfT[2] = powersOfT[1] * powersOfT[1];
      powersOfT[3] = powersOfT[2] * powersOfT[1];

      iIndex = static_cast<size_t>(iCoord);
      jIndex = static_cast<size_t>(jCoord);
    }


  } // namespace numeric

} // namespace brick

#endif /* #ifndef BRICK_NUMERIC_BSPLINE2D_IMPL_HH */
