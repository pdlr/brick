/**
***************************************************************************
* @file brick/computerVision/fitPolynomial_impl.hh
*
* Header file defining inline and template functions declared in
* fitPolynomial.hh.
*
* Copyright (C) 2019 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
***************************************************************************
**/

#ifndef BRICK_COMPUTERVISION_FITPOLYNOMIAL_IMPL_HH
#define BRICK_COMPUTERVISION_FITPOLYNOMIAL_IMPL_HH

// This file is included by fitPolynomial.hh, and should not be
// directly included by user code, so no need to include
// fitPolynomial.hh here.
//
// #include <brick/computerVision/fitPolynomial.hh>

#include <brick/linearAlgebra/linearAlgebra.hh>
#include <brick/numeric/array1D.hh>
#include <brick/numeric/array2D.hh>

namespace brick {

  namespace computerVision {

    // This function computes the best polynomial fit between the two input
    // sequences.
    template <class FloatType, class Iter0Type, class Iter1Type>
    brick::numeric::Polynomial<FloatType>
    fitPolynomial(Iter0Type xBegin, Iter0Type xEnd, Iter1Type yBegin,
                  std::size_t order)
    {
      using brick::linearAlgebra::linearLeastSquares;
      using brick::numeric::Array1D;
      using brick::numeric::Array2D;

      // Sadly, we need the size of the input sequences up front so we
      // can allocate storage.
      std::size_t numberOfElements = xEnd - xBegin;

      // Set up a system of linear equations to solve for the polynomial
      // coefficients.
      Array2D<FloatType> AMatrix(numberOfElements, order + 1);
      Array1D<FloatType> bVector(numberOfElements);

      FloatType constexpr one{1.0};
      for(std::size_t ii = 0; ii < numberOfElements; ++ii) {

        // Each row of the "A" matrix has powers of one of the X values.
        // Multiplying these powers of X by the polynomial coefficients,
        // and summing the results, should give us the least squares
        // approximation of the corresponding Y value.
        FloatType const xx{*xBegin};
        FloatType const yy{*yBegin};
        AMatrix(ii, 0) = one;
        for(std::size_t jj = 1; jj <= order; ++jj) {
          AMatrix(ii, jj) = AMatrix(ii, jj - 1) * xx;
        }
        bVector(ii) = yy;
        ++xBegin;
        ++yBegin;
      }

      // All set.  Now we can solve for the coefficients that most nearly
      // match the input data.
      Array1D<FloatType> coefficients = linearLeastSquares(AMatrix, bVector);
      return brick::numeric::Polynomial<FloatType>(coefficients);
    }

  } // namespace computerVision

} // namespace brick

// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/fitPolynomial_impl.hh>

#endif // #ifndef BRICK_COMPUTERVISION_FITPOLYNOMIAL_IMPL_HH
