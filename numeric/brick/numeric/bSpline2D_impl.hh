/**
***************************************************************************
* @file brick/numeric/bSpline2D_impl.hh
*
* Header file defining inline and template functions from BSpline2D.hh.
*
* Copyright (C) 2006-2014 David LaRose, dlr@cs.cmu.edu
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

    // This constructor builds a BSpline2D instance of unspecified
    // length and width.
    template <class Type, class FloatType>
    BSpline2D<Type, FloatType>::
    BSpline2D(bool isIsotropic)
      : m_basisArray(4),
        m_controlGrid(),
        m_isIsotropic(isIsotropic),
        m_minimumXY(0.0, 0.0),
        m_maximumXY(0.0, 0.0),
        m_xyCellOrigin(0.0, 0.0),
        m_xyCellSize(0.0, 0.0)
    {
      // Temporary storage for polynomial coefficients.
      Array1D<FloatType> basisCoefficients(4);

      // // The following basis functions are copied from Lee,
      // // Wolberg, and Shin, "Scattered Data Interpolation with
      // // Multilevel B-Splines, IEEE Transactions on Visualization and
      // // Computer Graphics, Vol 3, 228-244, 1997.  We present them
      // // here without further justification, although we suspect
      // // you could derive them recursively using a 2D version of
      // // the procedure used in BSpline::computeBasisFunction() from
      // // file bSpline.hh.
      
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

      // Control grid, etc., are left in an invalid state until we
      // have some data.
    }

    
    // The copy constructor does a deep copy.
    template <class Type, class FloatType>
    BSpline2D<Type, FloatType>::
    BSpline2D(BSpline2D<Type, FloatType> const& other)
      : m_basisArray(other.m_basisArray.size()),
        m_controlGrid(other.m_controlGrid.copy()),
        m_isIsotropic(other.m_isIsotropic),
        m_minimumXY(other.m_minimumXY),
        m_maximumXY(other.m_maximumXY),
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
    template <class Type, class FloatType>
    template <class CoordIter, class ObsIter>
    void
    BSpline2D<Type, FloatType>::
    approximateScatteredData(CoordIter sBegin,
                             CoordIter sEnd,
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
      this->approximateScatteredData(sBegin, sEnd, tBegin, observationsBegin,
                                     corner0, corner1);
    }


    // This function allows the spline parameters to be automatically
    // set in order to approximate an irregularly sampled function.
    template <class Type, class FloatType>
    template <class CoordIter, class ObsIter>
    void
    BSpline2D<Type, FloatType>::
    approximateScatteredData(CoordIter sBegin,
                               CoordIter sEnd,
                               CoordIter tBegin,
                               ObsIter observationsBegin,
                               Vector2D<FloatType> corner0,
                               Vector2D<FloatType> corner1)
    {
      // Sanity check input corners (cleanupCorners() is declared in
      // brick/numeric/utilities.hh).
      cleanupCorners(corner0, corner1);
      
      // Remember the specified bounds for reconstruction
      m_minimumXY = corner0;
      m_maximumXY = corner1;

      // Recover control grid size, and sanity check it.      
      size_t numberOfNodesS = m_controlGrid.columns();
      size_t numberOfNodesT = m_controlGrid.rows();
      if ((numberOfNodesS <= 3) || (numberOfNodesT <= 3)){
        BRICK_THROW(brick::common::StateException,
                    "BSpline2D::approximateScatteredData()",
                    "Control grid has not been initialized.  Please call "
                    "member function setNumberOfNodes() before calling "
                    "approximateScatteredData().");
      }

      // Chose the origin and spacing of the control grid.
      this->m_xyCellSize.setValue(
        ((m_maximumXY.x() - m_minimumXY.x())
         / static_cast<FloatType>(numberOfNodesS - 3)),
        ((m_maximumXY.y() - m_minimumXY.y())
         / static_cast<FloatType>(numberOfNodesT - 3)));

      if(this->m_isIsotropic) {
        FloatType cellSize = std::max(this->m_xyCellSize.x(),
                                      this->m_xyCellSize.y());
        this->m_xyCellSize.setValue(cellSize, cellSize);
      }

      this->m_xyCellOrigin = this->m_minimumXY - this->m_xyCellSize;
      
      // Allocate some space for intermediate grids, as described in
      // Lee, Wolberg, and Shin.
      Array2D<Type> deltaGrid(this->m_controlGrid.rows(),
                              this->m_controlGrid.columns());
      Array2D<FloatType> omegaGrid(this->m_controlGrid.rows(),
                                   this->m_controlGrid.columns());
      deltaGrid = static_cast<Type>(0.0);
      omegaGrid = static_cast<FloatType>(0.0);
    
      // This code implements the algorithm on page 231 of the paper.
      size_t iIndex;
      size_t jIndex;
      Array2D<FloatType> weightArray(4, 4);
      FloatType powersOfS[4];
      FloatType powersOfT[4];

      // Iterate over each observation (each scattered data point).
      while(sBegin != sEnd) {
        FloatType weightSquaredSum = 0.0;

        // Sanity check input data.
        if(*sBegin < this->m_minimumXY.x()
           || *sBegin >= this->m_maximumXY.x()
           || *tBegin < this->m_minimumXY.y()
           || *tBegin >= this->m_maximumXY.y()) {

          BRICK_THROW(brick::common::ValueException,
                      "BSpline2D::approximateScatteredData()",
                      "Input datum is out of bounds.");
        }
        
        // This call sets the value of powersOfS and powersOfT,
        // and returns by reference the indices of the control grid
        // cell into which (s, t) falls.
        this->decomposeSamplePoint(*sBegin, *tBegin, iIndex, jIndex,
                                   powersOfS, powersOfT);

        // Now on with Lee, Wolberg, and Shin's algorithm.  The four
        // basis polynomials define the weights with which each
        // control point in the neighborhood will affect the BSpline
        // value at the position of the currently selected
        // observation.  Here we compute the values of the basis
        // polynomials, and multiply them to get those weights.  Later
        // we will use the weights to solve for the most appropriate
        // control point values.
        for(size_t kIndex = 0; kIndex < 4; ++kIndex) {

          // Compute one basis function value.
          FloatType B_k = std::inner_product(
            powersOfS, powersOfS + 4, (this->m_basisArray)[kIndex].data(),
            static_cast<FloatType>(0));
  
          for(size_t lIndex = 0; lIndex < 4; ++lIndex) {

            // Compute the second basis function value.
            FloatType B_l = std::inner_product(
              powersOfT, powersOfT + 4, (this->m_basisArray)[lIndex].data(),
              static_cast<FloatType>(0));

            // Multiply to get the relevant weight.  Indexing into
            // weightArray is (row, column), not (k, l).
            FloatType weight = B_k * B_l;
            weightArray(lIndex, kIndex) = weight;
            weightSquaredSum += weight * weight;            
          }
        }

        // Here we solve for control point values, assuming that each
        // control point is within range of at most one of the
        // scattered data points.  We also keep some statistics that
        // will be used later to resolve control points for which this
        // assumption is not true.
        for(size_t kIndex = 0; kIndex < 4; ++kIndex) {
          for(size_t lIndex = 0; lIndex < 4; ++lIndex) {
            size_t index0 = iIndex + kIndex - 1;
            size_t index1 = jIndex + lIndex - 1;

            // Solve directly for the control point value (phi)
            // following Equation 3 of [1].  Indexing into
            // weightArray is (row, column), not (k, l).
            FloatType weight = weightArray(lIndex, kIndex);
            Type phi = (weight / weightSquaredSum) * (*observationsBegin);

            // Sanity check.  This test should not pass because the
            // extent of the spline was set based on calls to
            // std::min_element() and std::max_element() at the
            // beginning of this function.
            if(index0 >= deltaGrid.columns() || index1 >= deltaGrid.rows()) {
              BRICK_THROW(brick::common::LogicException,
                          "BSpline2D::approximateScatteredData()",
                          "Spline bounds appear to have been set incorrectly.");
            }

            // Delta and omega are the numerator and denominator of
            // Equation 5 of [1].  Essentially, these keep a weighted
            // average of the phi values from each scattered data
            // observation that affects each control point.  Next, we
            // will use these to compute a final value of phi.
            FloatType weightSquared = weight * weight;
            deltaGrid(index1, index0) += phi * weightSquared;
            omegaGrid(index1, index0) += weightSquared;
          }
        }
        ++sBegin;
        ++tBegin;
        ++observationsBegin;
      }

      // Final averaging step in case neighboring input points want
      // different control grid values.
      for(size_t index0 = 0; index0 < m_controlGrid.size(); ++index0) {
        if(omegaGrid[index0] == static_cast<FloatType>(0.0)) {
          this->m_controlGrid[index0] = static_cast<Type>(0.0);
        } else {
          this->m_controlGrid[index0] = deltaGrid[index0] / omegaGrid[index0];
        }
      }
    }
    

    // This member queries the size of the underlying B-Spline
    // control grid.
    template <class Type, class FloatType>
    void
    BSpline2D<Type, FloatType>::
    getNumberOfNodes(size_t& numberOfNodesS,
                     size_t& numberOfNodesT) const
    {
      numberOfNodesS = this->m_controlGrid.columns();
      numberOfNodesT = this->m_controlGrid.rows();
    }
      
      
    // This member function returns the maximum values for the spline
    // parameters S and T.
    template <class Type, class FloatType>
    void
    BSpline2D<Type, FloatType>::
    getMaximumSAndTValues(FloatType& maximumS, FloatType& maximumT) const
    {
      maximumS = m_maximumXY.x();
      maximumT = m_maximumXY.y();
    }


    // This member function returns the minimum values for the spline
    // parameters S and T.
    template <class Type, class FloatType>
    void
    BSpline2D<Type, FloatType>::
    getMinimumSAndTValues(FloatType& minimumS, FloatType& minimumT) const
    {
      minimumS = m_minimumXY.x();
      minimumT = m_minimumXY.y();
    }


    // This member function doubles the resolution of the B-spline
    // control grid without changing the shape of the interpolated
    // function or altering its range.
    template <class Type, class FloatType>
    void 
    BSpline2D<Type, FloatType>::
    promote()
    {
      // Check that *this is a valid B-spline.
      size_t numberOfNodesS = this->m_controlGrid.columns();
      size_t numberOfNodesT = this->m_controlGrid.rows();
      if ((numberOfNodesS <= 3) || (numberOfNodesT <= 3)){
        BRICK_THROW(brick::common::StateException,
                    "BSpline2D::promote()",
                    "Control grid has not been initialized.  Please call "
                    "member function setNumberOfNodes() before calling "
                    "promote().");
      }

      // Member m_basisArray is unchanged.
      
      // Member m_controlGrid nearly doubles in height and width.  To
      // see why it doesn't actually double, consider that the valid
      // interpolation range is interior to the rectangle formed by
      // the second row of the control grid, the sencond column, the
      // next-to-last row, and the next-to-last column.  The very
      // outside border of control points is only there to support
      // interpolation within this valid area.  When the control grid
      // resolution doubles, the required supporting border halves in
      // width.  There is no need to represent the entire span of the
      // original control grid.  Put another way, a spline of order N
      // requires N + 3 control points.  We're doubling the order of
      // the spline, not the number of control points.
      size_t newNumberOfNodesS = 2 * (numberOfNodesS - 3) + 3;
      size_t newNumberOfNodesT = 2 * (numberOfNodesT - 3) + 3;
      Array2D<Type> newControlGrid(newNumberOfNodesT, newNumberOfNodesS);

      // To help with clean indexing below.
      size_t const inRowsMinus1 = this->m_controlGrid.rows() - 1;
      size_t const inColumnsMinus1 = this->m_controlGrid.columns() - 1;

      // Assign values to the higher resolution control grid.  This
      // step is briefly described on page 234 of [1], which
      // references more detail in [2].  One could also derive these
      // numbers by solving a system of linear equations requiring the
      // output of the resampled spline to match that of the original
      // spline.
      //
      // 2. T. Lyche and K. Morken, "Making the Oslo Algorithm More
      // Efficient," SIAM Journel of Numerical Analysis, Vol. 23, No
      // 3, June 1986.
      for(size_t ii = 0; ii < inColumnsMinus1; ++ii) {
        for(size_t jj = 0; jj < inRowsMinus1; ++jj) {

          // To help with clean indexing below.
          size_t const twoTimesIi = (ii << 1);
          size_t const twoTimesJj = (jj << 1);

          if(0 != ii && 0 != jj) {
            newControlGrid(twoTimesJj - 1, twoTimesIi - 1) =
              (this->m_controlGrid(jj - 1, ii - 1) 
               + this->m_controlGrid(jj + 1, ii - 1)
               + this->m_controlGrid(jj - 1, ii + 1)
               + this->m_controlGrid(jj + 1, ii + 1)
               + 6.0 * (this->m_controlGrid(jj, ii - 1)
                        + this->m_controlGrid(jj - 1, ii)
                        + this->m_controlGrid(jj + 1, ii)
                        + this->m_controlGrid(jj, ii + 1))
               + 36.0 * this->m_controlGrid(jj, ii)) / 64.0;
          }
          if(0 != ii) {
            newControlGrid(twoTimesJj, twoTimesIi - 1) =
              (this->m_controlGrid(jj, ii - 1) 
               + this->m_controlGrid(jj + 1, ii - 1)
               + this->m_controlGrid(jj, ii + 1) 
               + this->m_controlGrid(jj + 1, ii + 1)
               + 6.0 * (this->m_controlGrid(jj, ii)
                        + this->m_controlGrid(jj + 1, ii))) / 16.0;
          }
          if(0 != jj) {
            newControlGrid(twoTimesJj - 1, twoTimesIi) =
              (this->m_controlGrid(jj - 1, ii) 
               + this->m_controlGrid(jj + 1, ii)
               + this->m_controlGrid(jj - 1, ii + 1) 
               + this->m_controlGrid(jj + 1, ii + 1)
               + 6.0 * (this->m_controlGrid(jj, ii)
                        + this->m_controlGrid(jj, ii + 1))) / 16.0;
          }
          newControlGrid(twoTimesJj, twoTimesIi) =
            (this->m_controlGrid(jj, ii) 
             + this->m_controlGrid(jj + 1, ii)
             + this->m_controlGrid(jj, ii + 1) 
             + this->m_controlGrid(jj + 1, ii + 1)) / 4.0;
        }
      }
      this->m_controlGrid = newControlGrid;

      // Member m_minimumXY is unchanged.
      // Member m_maximumXY is unchanged.

      // New cell size is dictated by the new resolution, and should
      // be exactly half of the previous cell size.
      m_xyCellSize.setValue(
        ((m_maximumXY.x() - m_minimumXY.x())
         / static_cast<FloatType>(newNumberOfNodesS - 3)),
        ((m_maximumXY.y() - m_minimumXY.y())
         / static_cast<FloatType>(newNumberOfNodesT - 3)));

      if(this->m_isIsotropic) {
        FloatType cellSize = std::max(this->m_xyCellSize.x(),
                                      this->m_xyCellSize.y());
        this->m_xyCellSize.setValue(cellSize, cellSize);
      }

      m_xyCellOrigin = m_minimumXY - m_xyCellSize;
    }


    // This member function sets the values of the control points of
    // the spline.
    template <class Type, class FloatType>
    void
    BSpline2D<Type, FloatType>::
    setControlPoints(Array2D<Type> const& controlPoints)
    {
      if(controlPoints.rows() <= 3 || controlPoints.columns() <= 3) {
        BRICK_THROW(brick::common::ValueException,
                    "BSpline2D::setControlPoints()",
                    "For a bicubic spline, control points array must be at "
                    "least 4x4.");
      }
      m_controlGrid = controlPoints.copy();
      m_minimumXY.setValue(0.0, 0.0);
      m_maximumXY.setValue(static_cast<FloatType>(controlPoints.columns() - 3),
                           static_cast<FloatType>(controlPoints.rows() - 3));
      m_xyCellOrigin.setValue(-1.0, -1.0);
      m_xyCellSize.setValue(1.0, 1.0);
    }

    
    // This member function both specifies the number of nodes in
    // the spline and sets the node positions so that the spline is
    // "uniform".
    template <class Type, class FloatType>
    void
    BSpline2D<Type, FloatType>::
    setNumberOfNodes(size_t numberOfNodesS,
                     size_t numberOfNodesT)
    {
      if(numberOfNodesS <= 3 || numberOfNodesT <= 3) {
        BRICK_THROW(brick::common::ValueException,
                    "BSpline2D::setNumberOfNodes()",
                    "Both arguments must be greater than 3.");
      }
      m_controlGrid.reinit(numberOfNodesT, numberOfNodesS);
      m_controlGrid = 0.0;
      m_minimumXY.setValue(0.0, 0.0);
      m_maximumXY.setValue(static_cast<FloatType>(numberOfNodesS - 3),
                           static_cast<FloatType>(numberOfNodesT - 3));
      m_xyCellOrigin.setValue(-1.0, -1.0);
      m_xyCellSize.setValue(1.0, 1.0);
    }

    
    // The assigment operator does a deep copy.
    template <class Type, class FloatType>
    BSpline2D<Type, FloatType>&
    BSpline2D<Type, FloatType>::
    operator=(BSpline2D<Type, FloatType> const& other)
    {
      if(&other != this) {
        this->m_basisArray.reinit(other.m_basisArray.size());
        this->m_controlGrid = other.m_controlGrid.copy();
        this->m_isIsotropic = other.m_isIsotropic;
        this->m_minimumXY = other.m_minimumXY;
        this->m_maximumXY = other.m_maximumXY;
        this->m_xyCellOrigin = other.m_xyCellOrigin;
        this->m_xyCellSize = other.m_xyCellSize;

        // Deep copy basis coefficients.
        for(size_t index0 = 0; index0 < this->m_basisArray.size(); ++index0) {
          this->m_basisArray[index0] = (other.m_basisArray[index0]).copy();
        }
      }
      return *this;
    }


    // This operator updates a spline so that its output is exactly
    // the sum of its original output and the output of another
    // spline.
    template <class Type, class FloatType>
    BSpline2D<Type, FloatType>&
    BSpline2D<Type, FloatType>::
    operator+=(BSpline2D<Type, FloatType> const& other)
    {
      if(&other != this) {
        // Member m_basisArray is assumed to be the same between the
        // two splines.

        this->m_controlGrid += other.m_controlGrid;

        // Several members are assumed to be unchanged.
        // this->m_isIsotropic
        // this->m_minimumXY
        // this->m_maximumXY
        // this->m_xyCellOrigin
        // this->m_xyCellSize
      }
      return *this;
    }
    
    
    // This operator evaluates the spline at the specified values of
    // spline parameters s and t.
    template <class Type, class FloatType>
    Type
    BSpline2D<Type, FloatType>::
    operator()(FloatType sValue, FloatType tValue) const
    {
      // This call sets the value of powersOfS and powersOfT,
      // and returns by reference the indices of the control grid
      // cell into which (s, t) falls.
      size_t iIndex;
      size_t jIndex;
      FloatType powersOfS[4];
      FloatType powersOfT[4];
      this->decomposeSamplePoint(sValue, tValue, iIndex, jIndex,
                                 powersOfS, powersOfT);
      
      // Interpolate by adding spline basis functions from the
      // surrounding control points.
      int index0 = iIndex - 1;
      int index1 = jIndex - 1;
      FloatType functionValue = 0.0;
      for(size_t kIndex = 0; kIndex < 4; ++kIndex) {
 
        size_t i0PlusK = index0 + kIndex;

        //  FloatType B_k = dot<FloatType>((this->m_basisArray)[kIndex], powersOfS);
        FloatType B_k = std::inner_product(
          powersOfS, powersOfS + 4, (this->m_basisArray)[kIndex].data(),
          static_cast<FloatType>(0));

        for(size_t lIndex = 0; lIndex < 4; ++lIndex) {

          size_t i1PlusL = index1 + lIndex;

          // FloatType B_l = dot<FloatType>((this->m_basisArray)[lIndex], powersOfT);
          FloatType B_l = std::inner_product(
            powersOfT, powersOfT + 4, (this->m_basisArray)[lIndex].data(),
            static_cast<FloatType>(0));

          // Indexing into control grid is (row, column), not (k, l).
          functionValue += (B_k * B_l * m_controlGrid(i1PlusL, i0PlusK));
        }
      }
      return functionValue;
    }


    template <class Type, class FloatType>
    void
    BSpline2D<Type, FloatType>::
    decomposeSamplePoint(FloatType sValue, FloatType tValue,
                         size_t& iIndex, size_t& jIndex,
                         FloatType* powersOfS, FloatType* powersOfT) const
    {
      // Note(xxx): consider using std::modf() here.
      
      // Find the integer coords of the control grid cell in which the
      // input point lies.
      FloatType iTmp = (sValue - m_xyCellOrigin.x()) / m_xyCellSize.x();
      FloatType jTmp = (tValue - m_xyCellOrigin.y()) / m_xyCellSize.y();

      int iCoord = static_cast<int>(roundToFloor(iTmp));
      int jCoord = static_cast<int>(roundToFloor(jTmp));

      // Find real valued coords within the cell, along with all of
      // the powers of those coords we'll be wanting to plug into
      // spline basis functions.
      powersOfS[0] = FloatType(1.0);
      powersOfS[1] = iTmp - FloatType(iCoord);
      powersOfS[2] = powersOfS[1] * powersOfS[1];
      powersOfS[3] = powersOfS[2] * powersOfS[1];
      powersOfT[0] = FloatType(1.0);
      powersOfT[1] = jTmp - FloatType(jCoord);
      powersOfT[2] = powersOfT[1] * powersOfT[1];
      powersOfT[3] = powersOfT[2] * powersOfT[1];

      iIndex = static_cast<size_t>(iCoord);
      jIndex = static_cast<size_t>(jCoord);
    }


  } // namespace numeric

} // namespace brick

#endif /* #ifndef BRICK_NUMERIC_BSPLINE2D_IMPL_HH */
