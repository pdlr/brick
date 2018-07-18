/**
***************************************************************************
* @file brick/computerVision/eightPointAlgorithm_impl.hh
*
* Header file defining inline and template functions from
* eightPointAlgorithm.hh.  function template.
*
* Copyright (C) 2008,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_EIGHTPOINTALGORITHM_IMPL_HH
#define BRICK_COMPUTERVISION_EIGHTPOINTALGORITHM_IMPL_HH

// This file is included by eightPointAlgorithm.hh, and should not be
// directly included by user code, so no need to include
// eightPointAlgorithm.hh here.
//
// #include <brick/computerVision/eightPointAlgorithm.hh>

#include <cmath>
#include <brick/linearAlgebra/linearAlgebra.hh>
#include <brick/numeric/utilities.hh>

namespace brick {

  namespace computerVision {

    template<class FloatType, class Iterator>
    brick::numeric::Array2D<FloatType>
    eightPointAlgorithm(Iterator sequence0Begin, Iterator sequence0End,
                        Iterator sequence1Begin)
    {
      brick::numeric::Array1D<FloatType> eigenvalues;
      return eightPointAlgorithm(
        sequence0Begin, sequence0End, sequence1Begin, eigenvalues);
    }


    template<class FloatType, class Iterator>
    brick::numeric::Array2D<FloatType>
    eightPointAlgorithm(Iterator sequence0Begin, Iterator sequence0End,
                        Iterator sequence1Begin,
                        brick::numeric::Array1D<FloatType>& eigenvalues)
    {
      // Find out how many points we have.  Even if this subtraction
      // is O(N), it will be dominated by the matrix multiplication
      // below.
      size_t numberOfCorrespondences = sequence0End - sequence0Begin;


      // Following Hartley, precondition the data by translating and
      // scaling (independently for each image) so that the points
      // roughly form a unit circle.  This greatly improves the
      // stability of the math below.
      brick::numeric::Array2D<FloatType> inputArray(numberOfCorrespondences, 3);
      brick::numeric::Array2D<FloatType> inputPrimeArray(numberOfCorrespondences, 3);
      Iterator begin0 = sequence0Begin;
      Iterator begin1 = sequence1Begin;
      size_t pointIndex = 0;
      while(begin0 != sequence0End) {
        brick::numeric::Array1D<FloatType> tmpRow = inputArray.getRow(pointIndex);
        tmpRow[0] = begin0->x();
        tmpRow[1] = begin0->y();
        tmpRow[2] = 1.0;

        tmpRow = inputPrimeArray.getRow(pointIndex);
        tmpRow[0] = begin1->x();
        tmpRow[1] = begin1->y();
        tmpRow[2] = 1.0;

        ++pointIndex;
        ++begin0;
        ++begin1;
      }

      brick::numeric::Array2D<FloatType> KKInv;
      brick::numeric::Array2D<FloatType> normalizedPoints;
      normalizePointSequence(inputArray, normalizedPoints, KKInv);

      brick::numeric::Array2D<FloatType> KPrimeInv;
      brick::numeric::Array2D<FloatType> normalizedPrimePoints;
      normalizePointSequence(inputPrimeArray, normalizedPrimePoints, KPrimeInv);

      // For each pair of points u, u', the fundamental matrix, F,
      // satisfies the equation:
      //
      //   transpose(u') * F * u = 0
      //
      // where the points u are drawn from sequence0, and the points
      // u' are drawn from sequence1.
      //
      // We rearrange this equation to get
      //
      //   ||f_00*u_x*u'_x + f_01*u_y*u'_x + f_02*u'_x
      //     + f_10*u_x*u'_y + f_11*u_y*u'_y + f_12*u'_y
      //     + f_20*u_x + f_21*u_y + f_22|| = 0
      //
      // where f_ij is the element of F in the i^th row and the j^th
      // column, u_x & u_y are the x & y components of u,
      // respectively, and u'_x & u'_y are the x & y components of u',
      // respectively.
      //
      // or,
      //
      //   ||A * vec(F)|| = 0
      //
      // With the matrix A as specified in the code below.

      brick::numeric::Array2D<FloatType> AMatrix(numberOfCorrespondences, 9);
      for(size_t rowIndex = 0; rowIndex < numberOfCorrespondences; ++rowIndex) {
        brick::numeric::Array1D<FloatType> currentRow =
          AMatrix.getRow(rowIndex);
        const brick::numeric::Array1D<FloatType>& uu =
          normalizedPoints.getRow(rowIndex);
        const brick::numeric::Array1D<FloatType>& uPrime =
          normalizedPrimePoints.getRow(rowIndex);
        currentRow[0] = uu[0] * uPrime[0];
        currentRow[1] = uu[1] * uPrime[0];
        currentRow[2] = uu[2] * uPrime[0];
        currentRow[3] = uu[0] * uPrime[1];
        currentRow[4] = uu[1] * uPrime[1];
        currentRow[5] = uu[2] * uPrime[1];
        currentRow[6] = uu[0] * uPrime[2];
        currentRow[7] = uu[1] * uPrime[2];
        currentRow[8] = uu[2] * uPrime[2];
      }

      // Solve for the F that minimizes the residual in the least
      // squares sense.
      brick::numeric::Array2D<FloatType> ATA =
        brick::numeric::matrixMultiply<FloatType>(
          AMatrix.transpose(), AMatrix);
      brick::numeric::Array2D<FloatType> eigenvectors;
      brick::linearAlgebra::eigenvectorsSymmetric(
        ATA, eigenvalues, eigenvectors);
      brick::numeric::Array2D<FloatType> FMatrix(3, 3);
      FMatrix[0] = eigenvectors(0, 8);
      FMatrix[1] = eigenvectors(1, 8);
      FMatrix[2] = eigenvectors(2, 8);
      FMatrix[3] = eigenvectors(3, 8);
      FMatrix[4] = eigenvectors(4, 8);
      FMatrix[5] = eigenvectors(5, 8);
      FMatrix[6] = eigenvectors(6, 8);
      FMatrix[7] = eigenvectors(7, 8);
      FMatrix[8] = eigenvectors(8, 8);

      // Good.  Now we have an estimate for F.  Here we enforce that F
      // must not be full rank.
      brick::numeric::Array2D<FloatType> uArray;
      brick::numeric::Array1D<FloatType> sigmaArray;
      brick::numeric::Array2D<FloatType> vTransposeArray;
      brick::linearAlgebra::singularValueDecomposition(
        FMatrix, uArray, sigmaArray, vTransposeArray);
      uArray[2] = 0.0;
      for(size_t ii = 0; ii < sigmaArray.size(); ++ii) {
        vTransposeArray.getRow(ii) *= sigmaArray(ii);
      }
      FMatrix = brick::numeric::matrixMultiply<FloatType>(
        uArray, vTransposeArray);

      // Transform back to unnormalized coordinates.
      FMatrix =
        brick::numeric::matrixMultiply<FloatType>(
          brick::numeric::matrixMultiply<FloatType>(
            KPrimeInv.transpose(), FMatrix), KKInv);

      return FMatrix;
    }


    // This function is used internally by eightPointAlgorithm() to
    // translate and scale input points so that their mean lies at the
    // origin and they have isotropic unit variance.
    template <class FloatType>
    void
    normalizePointSequence(
      brick::numeric::Array2D<FloatType> const& inputPoints,
      brick::numeric::Array2D<FloatType>& outputPoints,
      brick::numeric::Array2D<FloatType>& transform)
    {
      // Following Hartley (and quoting almost directly from the
      // paper), let u_i = [u_i,x , u_i,y , 1.0]^T, and form the
      // matrix
      //
      //   E = sum(u_i * u_i^T).
      //
      // Since this matrix is symmetric and positive definite, we may
      // take its Choleski factorization
      //
      //   E = K * K^T
      //
      // where K is upper triangular. It follows that
      //
      //   sum(inv(K / sqrt(N)) * u_i * u_i^T * inv(K / sqrt(N))^T) = NI,
      //
      // where I is the identity matrix, [and N is the number of
      // elements in the sum]. [...] Note that K is upper triangular,
      // and so it represents an affine transformation.
      //
      // In other words transforming by K/sqrt(N) puts the centroid of
      // our points at the origin, and normalizes their variance.
      // Just what we need.
      brick::numeric::Array2D<FloatType> EE =
        brick::numeric::matrixMultiply<FloatType>(
          inputPoints.transpose(), inputPoints);
      FloatType sqrtN = std::sqrt(inputPoints.rows());

      brick::numeric::Array2D<FloatType> KK;
      brick::linearAlgebra::choleskyFactorization(EE, KK, false);
      KK /= sqrtN;
      brick::numeric::Array2D<FloatType> KKInv =
        brick::linearAlgebra::inverse(KK);

      // OK, now we've found our normalizing transform.  Apply it
      // to the input points.
      outputPoints = brick::numeric::matrixMultiply<FloatType>(
        inputPoints, KKInv.transpose());
      transform = KKInv;
    }

  } // namespace computerVision

} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_EIGHTPOINTALGORITHM_IMPL_HH */
