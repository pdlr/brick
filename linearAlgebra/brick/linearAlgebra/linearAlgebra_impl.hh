/**
***************************************************************************
* @file brick/linearAlgebra/linearAlgebra_impl.hh
*
* Header file defining inline and template functions declared in
* brick/linearAlgebra/linearAlgebra.hh
*
* Copyright (C) 2001-2017 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
***************************************************************************
**/

#ifndef BRICK_LINEARALGEBRA_LINEARALGEBRA_IMPL_HH
#define BRICK_LINEARALGEBRA_LINEARALGEBRA_IMPL_HH

// This file is included by linearAlgebra.hh, and should not be
// directly included by user code, so no need to include
// linearAlgebra.hh here.
//
// #include <brick/linearAlgebra/linearAlgebra.hh>

#include <brick/common/exception.hh>
#include <brick/numeric/utilities.hh>
#include <brick/numeric/numericTraits.hh>

namespace brick {

  namespace linearAlgebra {

    // This function computes the best linear fit between the two input
    // arrays.
    template <class FloatType>
    std::pair<FloatType, FloatType>
    linearFit(brick::numeric::Array1D<FloatType> const& array0,
              brick::numeric::Array1D<FloatType> const& array1)
    {
      // We're looking for constants a and b which most nearly (in the
      // least squares sense) satisfy the equation
      //
      //   a * array0 + b = array1
      //
      // Which can be rewritten
      //
      //   [array0[0], 1]   [a] = [array1[0]]
      //   [array0[1], 1] * [b]   [array1[1]]
      //   [array0[2], 1]         [array1[2]]
      //   ...                    ...
      //
      // Solving this using the Moore-Penrose pseudoinverse gives
      //
      //                                             -1
      //   [a]  =  [dot(array0, array0), sum(array0)]  * [dot(array0, array1)]
      //   [b]     [sum(array0),         N          ]    [sum(array1)        ]
      //
      // For compactness, define A and x, so we can write
      //
      //   [a] = inverse(A) * x
      //   [b]

      // First some argument checking.
      if(array0.size() != array1.size()) {
        std::ostringstream message;
        message << "Arguments array0 and array1 must have the same size, "
                << "but are of size " << array0.size()
                << " and " << array1.size() << " respectively." << std::endl;
        BRICK_THROW(brick::common::ValueException,
                    "linearFit(Array1D const&, Array1D const&)",
                    message.str().c_str());
      }
      if(array0.size() <= 1) {
        BRICK_THROW(brick::common::ValueException,
                    "linearFit(Array1D const&, Array1D const&)",
                    "Arguments cannot have size zero or one.");
      }

      // Compute the right side of the equation above.  A and x are so
      // small that we simply represent their elements individually.
      // Note that A(1, 0) == A(0, 1).

      FloatType a00 = brick::numeric::dot<FloatType>(array0, array0);
      FloatType a01 = brick::numeric::sum<FloatType>(array0);
      FloatType a11 = array0.size();

      FloatType x0 = brick::numeric::dot<FloatType>(array0, array1);
      FloatType x1 = brick::numeric::sum<FloatType>(array1);

      // Inverse of a 2x2 matrix is easy to compute via the cofactor
      // method.  Rather than do so explicitly, though, we simply
      // multiply the vector x by the matrix of cofactors [a11, -a10;
      // -a01, a00], and then divide by the determinant. This saves us
      // two divisions.
      FloatType determinant = a00 * a11 - a01 * a01;
      if(brick::common::absoluteValue(determinant)
         < brick::numeric::NumericTraits<FloatType>::epsilon()) {
        BRICK_THROW(RankException,
                    "linearFit(Array1D const&, Array1D const&)",
                    "Problem solving for slope and offset. Perhaps line is "
                    "vertical, or not enough points have been supplied.");
      }

      return std::make_pair((a11 * x0 - a01 * x1) / determinant,
                            (a00 * x1 - a01 * x0) / determinant);
    }

  } // namespace linearAlgebra

} // namespace brick

#endif // #ifndef BRICK_LINEARALGEBRA_LINEARALGEBRA_IMPL_HH
