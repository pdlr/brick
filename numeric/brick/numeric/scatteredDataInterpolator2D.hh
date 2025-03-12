/**
***************************************************************************
* @file brick/numeric/scatteredDataInterpolator2D.hh
*
* Header file declaring the ScatteredDataInterpolator2D class.
*
* Copyright (C) 2014 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_NUMERIC_SCATTEREDDATAINTERPOLATOR2D_HH
#define BRICK_NUMERIC_SCATTEREDDATAINTERPOLATOR2D_HH

#include <vector>
#include <brick/numeric/bSpline2D.hh>

namespace brick {

  namespace numeric {

    /**
     ** This functor always returns false.  It is used as a default
     ** argument in ScatteredDataInterpolator2D::approximate().  Users
     ** can replace it with a functor that returns true to indicate
     ** when approximation residuals are sufficiently small.
     **/
    template <class ResidualType>
    struct FailFunctor {
      bool operator()(ResidualType const&) {return false;}
    };


    /**
     ** Warning: This class is very new, and its test suite is
     ** incomplete.  It almost certainly contain bugs, and its
     ** interface may change.
     **
     ** This class template implements the scattered data
     ** interpolation algorithm of Lee, Wolberg, and Shin [1].  This
     ** algorithm constructs bicubic 2D B-splines of progressively
     ** finer resolution to approximate input data.
     **
     ** Template argument Type describes the quantity to be
     ** interpolated.  This is the type of the scattered data.
     **
     ** Template argument FloatType describes the type that will be
     ** used for internal calculations, and is typically set to float
     ** or double.
     **
     ** Use this class as follows:
     **
     ** @code
     **   TBD
     ** @endcode
     **
     ** 1. Seungyong Lee, George Wolberg, and Sung Yong Shin,
     **    "Scattered Data Interpolation with Multilevel B-Splines," IEEE
     **    Transactions on Visualization and Computer Graphics, Vol. 3,
     **    No. 3, July–September 1997.
     **/
    template < class Type, class FloatType = double,
               class TestType = FailFunctor<Type> >
    class ScatteredDataInterpolator2D {
    public:

      /**
       * This constructor builds a ScatteredDataInterpolator2D instance
       * of unspecified length and width.
       *
       * @param numberOfLevels This argument specifies how many levels
       * of spline interpolation are to be performed.  The ultimate
       * size of the grid of spline control points is approximately
       * 2^numberOfLevels by 2^numberOfLevels.  Keep in mind that each
       * interpolated point draws from a 4x4 support region in the
       * grid of spline control points.
       *
       * @param isMeanCentered If this argument is true, the
       * interpolating function will include a DC shift so that its
       * "resting" value is the mean value of the data points to be
       * interpolated.  That is, if this argument is true, the
       * interpolating function will relax to the mean of the input
       * data in regions where there are no samples to interpolate.
       * If this argument is false, then the interpolating function
       * will relax to zero in areas where there are no data.
       *
       * @param isIsotropic If this argument is true, the spatial
       * resolution of the interpolating function will be the same in
       * each of the two interpolated axes.  If this argument is set
       * to false, the spacing of spline control points will differ
       * between the two axes so that resolution will be higher in the
       * axis-aligned direction along which the data are most tightly
       * grouped.
       */
      ScatteredDataInterpolator2D(size_t numberOfLevels = 8,
                                  bool isMeanCentered = true,
                                  bool isIsotropic = true);


      /**
       * The copy constructor does a deep copy.
       *
       * @param other This argument is the ScatteredDataInterpolator2D
       * instance to be copied.
       */
      ScatteredDataInterpolator2D(
        ScatteredDataInterpolator2D<Type, FloatType, TestType> const& other);


      /**
       * This function specifies the data to be interpolated.  With
       * each call to this function, any previous approximation is
       * discarded, and the interpolating function is re-estimated.
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
       * sequences of S and T coordinates described above.  The
       * interpolating function will be set so that, for 0 <= N <
       * (observedSPositionsEnd - observedSPositionsBegin), the value
       * of the interpolating function at (S, T) =
       * (*(observedSPositionsBegin + N), *(observedTPositionsBegin +
       * N)) approximates *(observationsBegin + N) as closely as
       * possible.
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
      approximate(CoordIter sBegin, CoordIter sEnd,
                  CoordIter tBegin,
                  ObsIter observationsBegin,
                  FloatType buffer = 1.0E-5);


      /**
       * This function specifies the data to be interpolated.  With
       * each call to this function, any previous approximation is
       * discarded, and the interpolating function is re-estimated.
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
       * sequences of S and T coordinates described above.  The
       * interpolating function will be set so that, for 0 <= N <
       * (observedSPositionsEnd - observedSPositionsBegin), the value
       * of the interpolating function at (S, T) =
       * (*(observedSPositionsBegin + N), *(observedTPositionsBegin +
       * N)) approximates *(observationsBegin + N) as closely as
       * possible.
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
      approximate(CoordIter sBegin, CoordIter sEnd,
                  CoordIter tBegin,
                  ObsIter observationsBegin,
                  Vector2D<FloatType> const& corner0,
                  Vector2D<FloatType> const& corner1);


      /**
       * This member function returns the maximum value for the
       * interpolating function parameters S and T.  Calling
       * operator()(FloatType, FloatType) with arguments greater than
       * or equal to those reported by getMaximumSAndTValues() is an
       * error.
       *
       * @param maximumS This argument is used to return the maximum
       * value of parameter S by reference.
       *
       * @param maximumT This argument is used to return the maximum
       * value of parameter T by reference.
       */
      void
      getMaximumSAndTValues(FloatType& maximumS, FloatType& maximumT) const;


      /**
       * This member function returns the minimum value for the
       * interpolating function parameters S and T.  Calling
       * operator()(FloatType, FloatType) with arguments less than
       * those reported by getMinimumSAndTValues() is an error.
       *
       * @param minimumS This argument is used to return the minimum
       * value of parameter S by reference.
       *
       * @param minimumT This argument is used to return the minimum
       * value of parameter T by reference.
       */
      void
      getMinimumSAndTValues(FloatType& minimumS, FloatType& minimumT) const;


      /**
       * This member function specifies a functor that is used to
       * test the quality of the approximation.  At each iteration (up
       * to the number of iterations specified by constructor argument
       * numberOfLevels), a residual is computed for each observation
       * by subtracting the interpolated value from the observation.
       * Each residual is evaluated by passing it to this functor, and
       * if the return value of this functor is true for each
       * residual, then iteration will terminate.  This allows the
       * caller to avoid additional computation once the interpolated
       * function is sufficiently close to the input data.
       *
       * For example, if the interpolated type is double, one might
       * pass std::bind2nd(std::less<double>(), 0.1) as the value of
       * this argument.  Or in C++11, one might pass [](double
       * x){return x < 0.1;} as the value of this argument.  The
       * default value of this argument always returns false, and will
       * never terminate the iteration early.
       *
       * @param testFunctor This argument is the functor to be applied
       * to each residual.
       */
      void
      setTestFunctor(TestType const& testFunctor) {
        this->m_testFunctor = testFunctor;
      }


      /**
       * The assigment operator does a deep copy.
       *
       * @param other This argument is the ScatteredDataInterpolator2D
       * instance to be copied.
       */
      ScatteredDataInterpolator2D<Type, FloatType, TestType>&
      operator=(
        ScatteredDataInterpolator2D<Type, FloatType, TestType> const& other);


      /**
       * This operator evaluates the interpolating function at the
       * specified values of spline parameters S and T.
       *
       * @return The return value is the calculated spline value.
       */
      Type
      operator()(FloatType sValue, FloatType tValue) const;


    protected:

      BSpline2D<Type, FloatType> m_bSpline2D;
      bool m_isMeanCentered;
      Type m_meanValue;
      size_t m_numberOfLevels;
      TestType m_testFunctor;
    };

  } // namespace numeric

} // namespace brick

// Include file containing definitions of inline and template
// functions.
#include <brick/numeric/scatteredDataInterpolator2D_impl.hh>

#endif /* #ifndef BRICK_NUMERIC_SCATTEREDDATAINTERPOLATOR2D_HH */
