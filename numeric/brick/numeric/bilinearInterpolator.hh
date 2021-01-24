/**
***************************************************************************
* @file brick/numeric/bilinearInterpolator.hh
*
* Header file declaring BilinearInterpolator class.
*
* Copyright (C) 1999-2007,2011,2015 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_BILINEARINTERPOLATOR_HH
#define BRICK_NUMERIC_BILINEARINTERPOLATOR_HH

#include <brick/common/types.hh>
#include <brick/numeric/array2D.hh>

namespace brick {

  namespace numeric {

    /**
     ** This class provides an easy way to bilinearly interpolate
     ** among elements of an Array2D instance.  Use it like this:
     **
     ** @code
     **   BilinearInterpolator interpolator(myArray);
     **   double myInterpolatedValue = interpolator(3.9, 2.1);
     ** @endcode
     **
     ** Template argument TypeIn describes the contents of the input array.
     **
     ** Template argument TypeOut describes the type of the
     ** interpolation result.
     **
     ** Template argument FloatType describes the type of to be used
     ** to for internal calculations, as well as the type that will be
     ** used to "index" into the interpolated array.
     **/
    template <class TypeIn,
              class TypeOut = brick::common::Float64,
              class FloatType = brick::common::Float64>
    class BilinearInterpolator
    {
    public:

      /**
       * Constructor sets the array to be interpolated.
       *
       * @param image0 This argument can be an Array2D instance of any
       * type that supports basic math functions.
       */
      BilinearInterpolator(const Array2D<TypeIn>& image0)
        : m_array(image0) {}


      /**
       * Destructor.
       */
      ~BilinearInterpolator() {}


      /**
       * The application operator performs a bilinear interpolation
       * among the elements of the array specified at construction,
       *
       * @return The interpolated value.
       */
      inline TypeOut
      operator()(FloatType row, FloatType column) const;

    private:

      inline void
      checkBounds(FloatType row, FloatType column) const;

      Array2D<TypeIn> m_array;
    };

  } // namespace numeric

} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/numeric/bilinearInterpolator_impl.hh>

#endif /* #ifndef BRICK_NUMERIC_BILINEARINTERPOLATOR_HH */
