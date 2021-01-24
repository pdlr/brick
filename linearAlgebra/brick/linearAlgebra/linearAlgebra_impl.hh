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
#include <brick/numeric/maxRecorder.hh>

namespace brick {

  namespace linearAlgebra {

    // Namespace for local symbols.  We don't make this anonymous
    // because we want to avoid reproducing these symbols in every
    // compilation that includes this file.
    //
    // Definitions are at the bottom of this file.
    namespace privateCode {

      template<class FloatType>
      bool
      hasLargerMagnitude(FloatType const& arg0, FloatType const& arg1) {
        FloatType abs0 = arg0 < FloatType(0) ? -arg0 : arg0;
        FloatType abs1 = arg1 < FloatType(1) ? -arg1 : arg1;
        return abs0 > abs1;
      }

      template<class FloatType>
      bool
      hasLargerMagnitude(std::complex<FloatType> const& arg0,
                         std::complex<FloatType> const& arg1) {
        FloatType abs0 = std::abs(arg0);
        FloatType abs1 = std::abs(arg1);
        return abs0 > abs1;
      }

      template<class FloatType>
      std::size_t
      selectPivotRow(brick::numeric::Array2D<FloatType> const& AA,
                     std::size_t startRow);

      template<class FloatType>
      void
      swapRows(brick::numeric::Array2D<FloatType>& AA,
               std::size_t row0,
               std::size_t row1);

    } // namespace privateCode


    // This function accepts a square Array2D instance and returns an
    // Array2D instance such that the matrix product of the two is
    // equal to the identity matrix.
    template <class FloatType>
    brick::numeric::Array2D<FloatType>
    inverse(brick::numeric::Array2D<FloatType> const& AA)
    {
      // First argument checking.
      if(AA.columns() != AA.rows()) {
        BRICK_THROW(brick::common::ValueException,
                    "inverse(Array2D<> const&)",
                    "Input array is not square.");
      }

      // Now set up some linear equations to solve.
      brick::numeric::Array2D<FloatType> AInverse =
        brick::numeric::identity<FloatType>(AA.rows(), AA.rows());
      brick::numeric::Array2D<FloatType> ACopy = AA.copy();

      // And solve for the inverse matrix.
      linearSolveInPlace(ACopy, AInverse); //Modifies AInverse.
      return AInverse;
    }


    // This function computes the best linear fit between the two input
    // arrays.
    template <class FloatType>
    std::pair<FloatType, FloatType>
    linearFit(brick::numeric::Array1D<FloatType> const& array0,
              brick::numeric::Array1D<FloatType> const& array1)
    {
      // We're looking for constants a and b that most nearly (in the
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


    // This function solves the system of equations A*x = b, where A and
    // b are known Array2D<double> instances.
    template<class FloatType>
    brick::numeric::Array1D<FloatType>
    linearLeastSquares(brick::numeric::Array2D<FloatType> const& AA,
                       brick::numeric::Array1D<FloatType> const& bb)
    {
      // First some argument checking.
      if(AA.size() == 0) {
        BRICK_THROW(brick::common::ValueException,
                    "linearLeastSquares(Array2D const&, Array1D const&)",
                    "Input array AA must have nonzero size.");
      }
      if(AA.columns() > AA.rows()) {
        BRICK_THROW(brick::common::ValueException,
                    "linearLeastSquares(Array2D const&, Array1D const&)",
                    "Input array AA must be square or have more rows "
                    "than columns.");
      }
      if(AA.rows() != bb.size()) {
        BRICK_THROW(brick::common::
                    ValueException,
                    "linearLeastSquares(Array2D const&, Array1D const&)",
                    "The number of rows in input array AA must be "
                    "the same as the number of elements in bb.");
      }

      // This two-line implementation works with non-built-in types, but
      // runs much slower than the specializations below.
      brick::numeric::Array2D<FloatType> APInv = pseudoinverse(AA);
      return brick::numeric::matrixMultiply<FloatType>(APInv, bb);
    }


    // Specializations are implemented in linearAlgebra.cc.
    template <>
    brick::numeric::Array1D<brick::common::Float32>
    linearLeastSquares(
      brick::numeric::Array2D<brick::common::Float32> const& AA,
      brick::numeric::Array1D<brick::common::Float32> const& bb);

    template <>
    brick::numeric::Array1D<brick::common::Float64>
    linearLeastSquares(
      brick::numeric::Array2D<brick::common::Float64> const& AA,
      brick::numeric::Array1D<brick::common::Float64> const& bb);


    // WARNING:  linearSolveInPlace() destructively modifies
    // both arguments!
    //
    // This function solves the system of equations A*x = b, where A is
    // a known matrix, and b is a known vector.
    template <class FloatType>
    void
    linearSolveInPlace(brick::numeric::Array2D<FloatType>& AA,
                       brick::numeric::Array1D<FloatType>& bb)
    {
      brick::numeric::Array2D<FloatType> bMatrix(bb.size(), 1, bb.data());
      linearSolveInPlace(AA, bMatrix);
    }


    // This function is identical to linearSolveInPlace(Array2D&, Array1D&),
    // except that b (and therefore x) is not constrained to be a vector.
    template <class FloatType>
    void
    linearSolveInPlace(brick::numeric::Array2D<FloatType>& AA,
                       brick::numeric::Array2D<FloatType>& bb)
    {
      // We use Gauss-Jordan elimination with partial pivoting to
      // handle the general case.  This implementation is slow, but
      // dispatch to LAPACK speeds us up in template specializations
      // for common types.  Because we use partial pivoting, rather
      // than full pivoting, we are more vulnerable to numerical
      // imprecision, but again the LAPACK specializations save us in
      // most cases.

      // First some argument checking.
      if(AA.rows() == 0) {
        BRICK_THROW(brick::common::ValueException,
                    "linearSolveInPlace(Array2D<Float32>&, Array2D<Float32>&)",
                    "Input array AA must not be empty.");
      }
      if(AA.rows() != bb.rows()) {
        BRICK_THROW(brick::common::ValueException,
                    "linearSolveInPlace(Array2D<Float32>&, Array2D<Float32>&)",
                    "Input arrays AA and bb must have the same number of "
                    "rows.");
      }
      if(AA.rows() != AA.columns()) {
        BRICK_THROW(brick::common::ValueException,
                    "linearSolveInPlace(Array2D<Float32>&, Array2D<Float32>&)",
                    "Input array AA must be square.");
      }

      // Make the A matrix upper triangular, and apply the same
      // transformations to bb.
      for(std::size_t rr = 0; rr < AA.rows() - 1; ++rr) {
        std::size_t pivotRow = privateCode::selectPivotRow(AA, rr);
        privateCode::swapRows(AA, rr, pivotRow);
        privateCode::swapRows(bb, rr, pivotRow);

        // Now subtract the current row from the rows below it to
        // triangularize.
        for(std::size_t subRow = rr + 1; subRow < AA.rows(); ++subRow) {
          FloatType const lambda = AA(subRow, rr) / AA(rr, rr);
          if(lambda != FloatType(0)) {

            // We can restrict operation to the nonzero region of AA,
            // and save a few flops.
            AA(subRow, rr) = FloatType(0);
            for(size_t cc = rr + 1; cc < AA.columns(); ++cc) {
              AA(subRow, cc) -= lambda * AA(rr, cc);
            }

            // But we have to operate on the entire row of bb, since
            // bb is not being triangularized.  There aren't any
            // regions that we know to already be zero.
            bb.getRow(subRow) -= lambda * bb.getRow(rr);
          }
        }
      }

      // Now that AA is upper triangular, work back up through the
      // rows, making it identity.  Apply the same transformations to
      // BB.  Note that we count on unsigned integer rollover to
      // terminate this for loop.
      for(std::size_t rr = AA.rows() - 1; rr < AA.rows(); --rr) {
        // Diagonal elements should not be zero because if they were
        // we would have thrown an exception in selectPivotRow() while
        // making upper triangular.

        // Here we rescale the diagonal element to 1.0, and apply the
        // rescaling to bb, too.  We have to rescale the entire row of
        // bb because our transformations haven't made it triangular /
        // partially diagonal.
        FloatType const lambda = AA(rr, rr);
        AA(rr, rr) = FloatType(1);
        bb.getRow(rr) /= lambda;

        // Now subtract the current row from the rows above it to
        // diagonalize AA.  Note that we count on unsigned integer
        // rollover to terminate this for loop.
        for(std::size_t subRow = rr - 1; subRow < AA.rows(); --subRow) {
          FloatType const beta = AA(subRow, rr);
          AA(subRow, rr) = FloatType(0);
          bb.getRow(subRow) -= beta * bb.getRow(rr);
        }
      }

      // All done.  AA is now identity, and bb has been transformed
      // identically, so our new equation is: I * X = bb', where bb'
      // is the transformed version of bb.  From this, it's clear that
      // bb' is the solution for X.
    }

    // Specializations are implemented in linearAlgebra.cc.
    template <>
    void
    linearSolveInPlace(brick::numeric::Array2D<brick::common::Float32>& AA,
                       brick::numeric::Array2D<brick::common::Float32>& bb);

    template <>
    void
    linearSolveInPlace(brick::numeric::Array2D<brick::common::Float64>& AA,
                       brick::numeric::Array2D<brick::common::Float64>& bb);



    // This function accepts an Array2D<Float64> instance having at least
    // as many rows as columns, and returns the Moore-Penrose
    // pseudoinverse.
    template <class FloatType>
    brick::numeric::Array2D<FloatType>
    pseudoinverse(brick::numeric::Array2D<FloatType> const& AA)
    {
      using brick::numeric::Array2D;
      using brick::numeric::matrixMultiply;

      Array2D<FloatType> ATranspose = AA.transpose();
      Array2D<FloatType> ATA = matrixMultiply<FloatType>(ATranspose, AA);
      return matrixMultiply<FloatType>(inverse(ATA), ATranspose);
    }



    // Namespace for local symbols.
    namespace privateCode {

      template<class FloatType>
      std::size_t
      selectPivotRow(brick::numeric::Array2D<FloatType> const& AA,
                     std::size_t startRow)
      {
        FloatType largestSoFar(0);
        std::size_t pivotRow = AA.rows();
        for(std::size_t rr = startRow; rr < AA.rows(); ++rr) {
          FloatType const& testValue = AA(rr, startRow);
          if(privateCode::hasLargerMagnitude(testValue, largestSoFar)) {
            largestSoFar = testValue;
            pivotRow = rr;
          }
        }

        if(pivotRow == AA.rows()) {
          BRICK_THROW(brick::common::ValueException,
                      "brick::linearAlgebra::selectPivotRow()",
                      "Matrix appears to be singular.");
        }
        return pivotRow;
      }


      template<class FloatType>
      void
      swapRows(brick::numeric::Array2D<FloatType>& AA,
               std::size_t row0,
               std::size_t row1)
      {
        for(std::size_t cc = 0; cc < AA.columns(); ++cc) {
          std::swap(AA(row0, cc), AA(row1, cc));
        }
      }

    } // namespace privateCode

  } // namespace linearAlgebra

} // namespace brick

#endif // #ifndef BRICK_LINEARALGEBRA_LINEARALGEBRA_IMPL_HH
