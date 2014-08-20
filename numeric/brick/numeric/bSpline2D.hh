/**
***************************************************************************
* @file brick/numeric/bSpline2D.hh
*
* Header file declaring the BSpline2D class.
*
* Copyright (C) 2006-2014 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_NUMERIC_BSPLINE2D_HH
#define BRICK_NUMERIC_BSPLINE2D_HH

#include <vector>
#include <brick/numeric/array1D.hh>
#include <brick/numeric/array2D.hh>
#include <brick/numeric/index2D.hh>
#include <brick/numeric/vector2D.hh>
#include <brick/numeric/polynomial.hh>

namespace brick {

  namespace numeric {

    /**
     ** Warning: This class is very new, and its test suite is
     ** incomplete.  It almost certainly contain bugs, and its
     ** interface may change.
     **
     ** This class template implements a bicubic 2D B-spline.  Other
     ** orders (e.g., biquadratic) are currently not supported.  This
     ** class is templated on control point type so that you can, for
     ** example, create a spline with 2D or 3D values by specifying
     ** control points of type Vector2D or Vector3D.  Note that the
     ** splines represented by this class are all 2D in the sense that
     ** you can represent them using a parametric function of two
     ** parameters.  If you want a spline that maps to a 1D parametric
     ** function (such as a curve in 2D or 3D space) you need to use
     ** the similar class, BSpline.  BSpline2D currently only supports
     ** non-periodic splines, and currently only supports uniform node
     ** spacing.  For now knot multiplicities are fixed at 1, meaning
     ** that folds and discontinuities are not supported.
     **
     ** The second template paramter, FloatType, controls the
     ** precision with which internal calculations are carried out.
     **
     ** This class implements the 2D spline used by
     ** ScatteredDataInterpolater, which is an implementation of the
     ** scattered data interpolation algorithm of Lee, Wolberg, and
     ** Shin [1].
     **
     ** 1. Seungyong Lee, George Wolberg, and Sung Yong Shin,
     **    "Scattered Data Interpolation with Multilevel B-Splines," IEEE
     **    Transactions on Visualization and Computer Graphics, Vol. 3,
     **    No. 3, July–September 1997.
     **/
    template <class Type, class FloatType = double>
    class BSpline2D {
    public:

      /** 
       * This constructor builds a BSpline2D instance of unspecified
       * length and width.  Currently only bicubic splines are
       * supported; we do not provide any way to construct splines of
       * order different that 3.
       * 
       * @param isIsotropic If this argument is true, the spatial
       * resolution of the spline control grid will be the same in
       * each of the two interpolated axes.  If this argument is set
       * to false, the spacing of spline control points will differ
       * between the two axes so that after a call to
       * approximateScatteredData(), the resolution will be higher in
       * the axis-aligned direction along which the data are most
       * tightly grouped.
       */
      explicit BSpline2D(bool isIsotropic = false);


      /** 
       * The copy constructor does a deep copy.
       * 
       * @param other This argument is the BSpline2D instance to be copied.
       */
      BSpline2D(BSpline2D<Type, FloatType> const& other);


      /** 
       * This function allows the spline parameters to be
       * automatically set in order to approximate an irregularly
       * sampled function.  If you have some randomly distributed
       * observations of a function value, and you want to approximate
       * them with a spline, this the right member function to call.
       * First you must create a spline, decide how many control
       * points it should have along each axis (set this using member
       * function setNumberOfNodes()), and then call
       * approximateScatteredData().
       * 
       * @param sBegin This iterator specifies the beginning of a
       * sequence of (possibly non-uniformly distributed) S
       * coordinates of points at which observations of the
       * to-be-approximated function were made.
       *
       * @param sEnd This iterator specifies the end of the sequence
       * of S coordinates.
       *
       * @param tBegin This iterator specifies the beginning of a
       * sequence of (possibly non-uniformly distributed) T
       * coordinates corresponding to the S coordinate sequence
       * described above.  The sequence of T coordinates must have at
       * least as many elements as the sequence of S coordinates.
       *
       * @param observationsBegin This iterator specifies the
       * beginning of a sequence of observations corresponding to the
       * sequences of S and T coordinates described above.  The spline
       * control points will be set so that, for 0 <= N <
       * (observedSPositionsEnd - observedSPositionsBegin), the value
       * of the spline at (S, T) = (*(observedSPositionsBegin + N),
       * *(observedTPositionsBegin + N)) approximates
       * *(observationsBegin + N) as closely as possible.
       *
       * @param buffer This argument specifies an amount by which the
       * valid range for the approximated function will extend past
       * the minimum and maximum coordinates of the observed
       * positions.  This is useful because sometimes it is necessary
       * (although not terribly accurate) to extrapolate beyond the
       * range of the input observations.
       */
      template <class CoordIter, class ObsIter>
      void
      approximateScatteredData(CoordIter sBegin,
                               CoordIter sEnd,
                               CoordIter tBegin,
                               ObsIter observationsBegin,
                               FloatType buffer = 1.0E-5);

      
      /** 
       * This function allows the spline parameters to be
       * automatically set in order to approximate an irregularly
       * sampled function.  If you have some randomly distributed
       * observations of a function value, and you want to approximate
       * them with a spline, this the right member function to call.
       * First you must create a spline, decide how many control
       * points it should have along each axis (set this using member
       * function setNumberOfNodes()), and then call
       * approximateScatteredData().
       * 
       * @param sBegin This iterator specifies the beginning of a
       * sequence of (possibly non-uniformly distributed) S
       * coordinates of points at which observations of the
       * to-be-approximated function were made.
       *
       * @param sEnd This iterator specifies the end of the sequence
       * of S coordinates.
       *
       * @param tBegin This iterator specifies the beginning of a
       * sequence of (possibly non-uniformly distributed) T
       * coordinates corresponding to the S coordinate sequence
       * described above.  The sequence of T coordinates must have at
       * least as many elements as the sequence of S coordinates.
       *
       * @param observationsBegin This iterator specifies the
       * beginning of a sequence of observations corresponding to the
       * sequences of S and T coordinates described above.  The spline
       * control points will be set so that, for 0 <= N <
       * (observedSPositionsEnd - observedSPositionsBegin), the value
       * of the spline at (S, T) = (*(observedSPositionsBegin + N),
       * *(observedTPositionsBegin + N)) approximates
       * *(observationsBegin + N) as closely as possible.
       *
       * @param corner0 This argument and the next define a
       * rectangular region in parameter space over which the
       * interpolated function will be valid.  It is an error if any
       * of the input (S, T) coordinates lie outside of this region.
       *
       * @param corner1 This argument and the previous define a
       * rectangular region in parameter space over which the
       * interpolated function will be valid.  It is an error if any
       * of the input (S, T) coordinates lie outside of this region.
       */
      template <class CoordIter, class ObsIter>
      void
      approximateScatteredData(CoordIter sBegin,
                               CoordIter sEnd,
                               CoordIter tBegin,
                               ObsIter observationsBegin,
                               Vector2D<FloatType> corner0,
                               Vector2D<FloatType> corner1);


      /**
       * Indicates whether the spacing of the spline control grid is
       * the same in both S and T directions.
       *
       * @return The return value is true if the spatial resolution of
       * the spline control grid is constrained to be the same along
       * each of the two interpolated axes.  If the return value is
       * false, the spacing of spline control points can differ
       * between the two axes.
       */
      bool getIsIsotropic() {return this->m_isIsotropic;}

      
      /** 
       * This member queries the size of the underlying B-Spline
       * control grid.
       * 
       * @param numberOfNodesS This argument returns how many nodes
       * the spline has along the S axis.
       * 
       * @param numberOfNodesT This argument returns how many nodes
       * the spline has along the T axis.
       */
      void
      getNumberOfNodes(size_t& numberOfNodesS,
                       size_t& numberOfNodesT) const;
      
      
      /** 
       * This member function returns the maximum value for the spline
       * parameters S and T.  Calling operator()(FloatType, FloatType) with
       * arguments greater than or equal to those reported by
       * getMaximumSAndTValues() is an error.
       * 
       * @param maximumS This argument is used to return the maximum
       * value of spline parameter S by reference.
       * 
       * @param maximumT This argument is used to return the maximum
       * value of spline parameter T by reference.
       */
      void
      getMaximumSAndTValues(FloatType& maximumS, FloatType& maximumT) const;


      /** 
       * This member function returns the minimum value for the spline
       * parameters S and T.  Calling operator()(FloatType, FloatType) with
       * arguments less than those reported by getMinimumSAndTValues()
       * is an error.
       * 
       * @param minimumS This argument is used to return the minimum
       * value of spline parameter S by reference.
       * 
       * @param minimumT This argument is used to return the minimum
       * value of spline parameter T by reference.
       */
      void
      getMinimumSAndTValues(FloatType& minimumS, FloatType& minimumT) const;


      /**
       * This member function doubles the resolution of the B-spline
       * control grid without changing the shape of the interpolated
       * function or altering its range.  Member function promote() is
       * primarily useful for adding the function represented by *this
       * to the function represented by a B-spline with twice the
       * resolution.
       */
      void 
      promote();


      /** 
       * This member function sets the values of the control points of
       * the spline.
       * 
       * @param controlPoints This argument specifies the control
       * point values for the spline.  It must have shape
       * (numberOfNodesS, numberOfNodesT), where numberOfNodesS and
       * numberOfNodesT are the same values passed to member function
       * setNumberOfNodes.
       */
      void
      setControlPoints(Array2D<Type> const& controlPoints);


      /** 
       * This member function both specifies the number of nodes in
       * the spline and sets the node positions so that the spline is
       * "uniform".  The node positions will be set so that the first
       * node lies at spline parameter (s, t) = (-1.0, -1.0) and
       * subsequent nodes lie at (0.0, 0.0), (0.0, 1.0), (1.0, 0.0),
       * etc.  Note that the actual number of nodes in the spline is
       * equal to numberOfNodesS * numberOfNodesT, because the nodes
       * form a 2D array.  Further note that the valid range for
       * interpolated s and t values starts at 0.0, and goes to
       * (numberOfNodesS - 3.0) and (numberOfNodesT - 3.0),
       * respectively.
       * 
       * @param numberOfNodesS This argument specifies how many nodes
       * the spline should have along the S axis.
       * 
       * @param numberOfNodesT This argument specifies how many nodes
       * the spline should have along the T axis.
       */
      void
      setNumberOfNodes(size_t numberOfNodesS,
                       size_t numberOfNodesT);
      

      /** 
       * The assigment operator does a deep copy.
       * 
       * @param other This argument is the BSpline2D instance to be copied.
       */
      BSpline2D<Type, FloatType>&
      operator=(BSpline2D<Type, FloatType> const& other);

      
      /** 
       * This operator updates a spline so that its output is exactly
       * the sum of its original output and the output of another
       * spline.  The argument of this function must have the same
       * minimumXY, maximumXY, xyCellSize, and controlGrid size as
       * *this.
       * 
       * @param other This argument is the BSpline2D instance to be
       * added to *this.
       */
      BSpline2D<Type, FloatType>&
      operator+=(BSpline2D<Type, FloatType> const& other);

      
      /** 
       * This operator evaluates the spline at the specified values of
       * spline parameters S and T.
       * 
       * @return The return value is the calculated spline value.
       */
      Type
      operator()(FloatType sValue, FloatType tValue) const;
      

    protected:

      /** 
       * This protected member function returns the integer part of
       * sValue and tValue, while -- as a side effect -- setting
       * member variables m_powersOfS and m_powersOfT so that the i^th
       * (counting from zero) element of m_powersOfS is equal to the
       * i^th power of the fractional part of sValue, and the i^th
       * element of m_powersOfT is equal to the i^th power of the
       * fractional part of tValue.
       * 
       * @param sValue This argument is the value of spline parameter
       * S to be decomposed.
       * 
       * @param tValue This argument  is the value of spline parameter
       * T to be decomposed.
       * 
       * @param iIndex This argument is used to return the integer
       * part of sValue.  For uniform spline (which *this always is),
       * this is equivalent to the index of the span into which sValue
       * falls.
       * 
       * @param jIndex This argument is used to return the integer
       * part of tValue.  For uniform spline (which *this always is),
       * this is equivalent to the index of the span into which tValue
       * falls.
       * 
       * @param powersOfS This argument must point to a four-element
       * array, into which the first four powers of the scalar offset,
       * s, will be written.
       *
       * @param powersOfT This argument must point to a four-element
       * array, into which the first four powers of the scalar offset,
       * t, will be written.
       */
      void
      decomposeSamplePoint(FloatType sValue, FloatType tValue,
                           size_t& iIndex, size_t& jIndex,
                           FloatType* powersOfS, FloatType* powersOfT) const;


      Array1D< Array1D<FloatType> > m_basisArray;
      Array2D<Type> m_controlGrid;
      bool m_isIsotropic;
      Vector2D<FloatType> m_minimumXY;
      Vector2D<FloatType> m_maximumXY;
      Vector2D<FloatType> m_xyCellOrigin;
      Vector2D<FloatType> m_xyCellSize;
    };
    
  } // namespace numeric
  
} // namespace brick

// Include file containing definitions of inline and template
// functions.
#include <brick/numeric/bSpline2D_impl.hh>

#endif /* #ifndef BRICK_NUMERIC_BSPLINE2D_HH */
