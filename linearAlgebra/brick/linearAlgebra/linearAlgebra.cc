/**
***************************************************************************
* @file brick/linearAlgebra/linearAlgebra.cpp
*
* Source file defining linear algebra functions.  Many of these depend
* on the LAPACK and BLAS libraries.
*
* Copyright (C) 2001-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
***************************************************************************
**/

#include <brick/common/exception.hh>
#include <brick/linearAlgebra/linearAlgebra.hh>
#include <brick/linearAlgebra/clapack.hh>
#include <brick/numeric/utilities.hh>

// Using directives for this source file only.
using namespace brick::common;
using namespace brick::numeric;

namespace brick {

  /**
   ** This namespace contains functions which implement common linear
   ** algebra tasks, such as eigenvector computations, SVD, etc.
   **/
  namespace linearAlgebra {

    void
    choleskyFactorization(Array2D<Float64> const& inputArray,
                          Array2D<Float64>& kArray,
                          bool isUpperTriangular)
    {
      // Argument checking.
      if(inputArray.size() == 0) {
        BRICK_THROW(brick::common::ValueException,
                    "choleskyFactorization()",
                    "Argument inputArray cannot have zero size.");
      }
      if(inputArray.rows() != inputArray.columns()) {
        BRICK_THROW(brick::common::ValueException,
                    "choleskyFactorization()",
                    "Argument inputArray must be square.");
      }
    
      // Transpose A to match LAPACK's convention.  Since inputArray
      // is symmetric, we don't really need to transpose!  Also, we
      // only need to copy the upper/lower triangular part, depending
      // on the value of isUpperTriangular.  For convenience, we copy
      // the whole matrix, though.
      size_t dimension = inputArray.rows();
      Array2D<Float64> aColumnMajor(dimension, dimension);
      aColumnMajor = 0.0;
      {
        Array2D<Float64>::const_iterator inPtr = inputArray.begin();
        Array2D<Float64>::iterator outPtr = aColumnMajor.begin();
        if(isUpperTriangular) {
          for(size_t index0 = 0; index0 < dimension; ++index0) {
            inPtr += index0;
            outPtr += index0;
            for(size_t index1 = index0; index1 < dimension; ++index1) {
              *(outPtr++) = *(inPtr++);
            }
          }
        } else {
          for(size_t index0 = 0; index0 < dimension; ++index0) {
            for(size_t index1 = 0; index1 <= index0; ++index1) {
              *(outPtr++) = *(inPtr++);
            }
            inPtr += (dimension - index0 - 1);
            outPtr += (dimension - index0 - 1);
          }
        }        
      }


      // Set up arguments for the LAPACK call.

      // Set UPLO opposite of what you might expect because LAPACK is
      // column major.
      char UPLO = isUpperTriangular ? 'L' : 'U';
      Int32 N = static_cast<Int32>(dimension);
      Int32 LDA = static_cast<Int32>(dimension);
      Int32 INFO;

      // Dispatch to the LAPACK routine that actually does the
      // factorization.
      dpotrf_(&UPLO, &N, aColumnMajor.data(), &LDA, &INFO);

      // Check for errors.
      if(INFO < 0L) {
        std::ostringstream message;
        message << "Call to dpotrf_ returns " << INFO
                << ".  Something is wrong.";
        BRICK_THROW(brick::common::ValueException, "choleskyFactorization()",
                    message.str().c_str());
      } else if(INFO > 0L) {
        std::ostringstream message;
        BRICK_THROW(brick::common::ValueException, "choleskyFactorization()",
                    "Input matrix is not positive definite.");
      }

      // Recover the result.
      kArray = aColumnMajor;
    }

    
    Float64
    determinant(Array2D<Float64> const& A)
    {
      // First argument checking.
      if(A.columns() != A.rows()) {
        BRICK_THROW(brick::common::ValueException,
                    "determinant(Array2D<Float64> const&)",
                    "Input array is not square.");
      }

      // In this routine, we take advantage of the fact that the
      // determinant of a matrix is related to the product of the
      // diagonal elements of its LU factorization.
    
      // Start by computing the LU factorization of A.
      Array2D<Float64> AColumnMajor = A.transpose();
      Int32 M = static_cast<Int32>(A.rows());
      Int32 N = static_cast<Int32>(A.columns());
      Int32 LDA = static_cast<Int32>(A.rows());
      Array1D<Int32> IPIV(M); // Really should be min(M, N), but A is square.
      Int32 info;
      dgetrf_(&M, &N, AColumnMajor.data(), &LDA, IPIV.data(), &info);

      // Check for errors.
      if(info < 0L) {
        std::ostringstream message;
        message << "Call to dgetrf_ returns " << info
                << ".  Something is wrong.";
        BRICK_THROW(brick::common::ValueException,
                    "determinant(Array2D<Float64> const&)",
                    message.str().c_str());
      }

      // Compute the product of the diagonal elements, and find out if
      // the determinant is equal to + or - the product of the diagonal
      // element. Do this second thing by counting how many row-swaps
      // were conducteed.
      Float64 determinant = 1.0;
      size_t swapCount = 0;
      for(size_t index = 0; index < AColumnMajor.rows(); ++index) {
        determinant *= AColumnMajor(index, index);
        if(IPIV[index] != static_cast<Int32>(index + 1)) {
          // Note(xxx): check that this is right.
          ++swapCount;
        }
      }

      // Finally, correct the sign
      if((swapCount % 2) == 1) {
        determinant *= -1.0;
      }
      return determinant;
    }

  
    Array1D<Float64>
    eigenvaluesSymmetric(Array2D<Float64> const& inputArray)
    {
      // Argument checking.
      if(inputArray.size() == 0) {
        BRICK_THROW(brick::common::ValueException,
                    "eigenvaluesSymmetric()",
                    "Argument inputArray cannot have zero size.");
      }
      if(inputArray.rows() != inputArray.columns()) {
        BRICK_THROW(brick::common::ValueException,
                    "eigenvaluesSymmetric()",
                    "Argument inputArray must be square.");
      }
    
      // Transpose A to match LAPACK's convention.  Since inputArray is
      // symmetric, we don't really need to transpose!  Also, we only
      // need to copy the upper triangular part.
      size_t dimension = inputArray.rows();
      Array2D<Float64> aColumnMajor(dimension, dimension);
      {
        Array2D<Float64>::const_iterator inPtr = inputArray.begin();
        Array2D<Float64>::iterator outPtr = aColumnMajor.begin();
        for(size_t index0 = 0; index0 < dimension; ++index0) {
          inPtr += index0;
          outPtr += index0;
          for(size_t index1 = index0; index1 < dimension; ++index1) {
            *(outPtr++) = *(inPtr++);
          }
        }
      }

      // Set up storage for return values.
      Array1D<Float64> eigenvalues(dimension);

      // Set up arguments for the LAPACK call.
      char JOBZ = 'N';  // Compute eigenvalues only.
      char UPLO = 'L';  // Get input from upper triangle of A.
      Int32 N = static_cast<Int32>(inputArray.rows());
      Int32 LDA = N;
      Float64 WORK;
      Int32 LWORK = -1;
      Int32 INFO;

      // Call once to request optimal workspace size.
      dsyev_(&JOBZ,  &UPLO, &N, aColumnMajor.data(), &LDA,
             eigenvalues.data(), &WORK, &LWORK, &INFO);

      // Check for errors.
      if(INFO != 0L) {
        std::ostringstream message;
        message << "First call to dsyev_ returns " << INFO
                << ".  Something is wrong.";
        BRICK_THROW(brick::common::ValueException,
                    "eigenvaluesSymmetric()",
                    message.str().c_str());
      }
    
      // Resize workspace.
      LWORK = static_cast<Int32>(WORK);
      Array1D<Float64> doubleWorkSpace(static_cast<size_t>(LWORK));

      // Call again to really compute the eigenvectors.
      dsyev_(&JOBZ,  &UPLO, &N, aColumnMajor.data(), &LDA,
             eigenvalues.data(), doubleWorkSpace.data(), &LWORK, &INFO);

      // Check for errors.
      if(INFO != 0L) {
        std::ostringstream message;
        message << "Second call to dsyev_ returns " << INFO
                << ".  Something is wrong.";
        BRICK_THROW(brick::common::ValueException,
                    "eigenvaluesSymmetric()",
                    message.str().c_str());
      }

      // Our convention differs from LAPACK's about the order of eigenvalues.
      std::reverse(eigenvalues.begin(), eigenvalues.end());
      return eigenvalues;
    }


    void
    eigenvectors(Array2D<Float64> const& inputArray,
                 Array1D< std::complex<Float64> >& eigenvalues,
                 Array2D< std::complex<Float64> >& eigenvectors,
                 bool isSortRequired)
    {
      // Argument checking.
      if(inputArray.size() == 0) {
        BRICK_THROW(brick::common::ValueException,
                    "eigenvectors()",
                    "Argument inputArray cannot have zero size.");
      }
      if(inputArray.rows() != inputArray.columns()) {
        BRICK_THROW(brick::common::ValueException,
                    "eigenvectors()",
                    "Argument inputArray must be square.");
      }
    
      // Transpose A to match LAPACK's convention.
      size_t dimension = inputArray.rows();
      Array2D<Float64> aTransposeColumnMajor = inputArray.copy();

      // Set up storage for return values.
      Array1D<Float64> eigenvaluesReal(dimension);
      Array1D<Float64> eigenvaluesImag(dimension);
      Array2D<Float64> leftEigenvectorsTranspose(dimension, dimension);
      Array2D<Float64> rightEigenvectorsTranspose(dimension, dimension);

      // Set up arguments for the LAPACK call.
      char JOBVL = 'V'; // Don compute left eigenvectors.
      char JOBVR = 'N'; // Don't compute right eigenvectors.
      Int32 N = static_cast<Int32>(dimension);
      Int32 LDA = static_cast<Int32>(dimension);
      Int32 LDVL = static_cast<Int32>(dimension);
      Int32 LDVR = static_cast<Int32>(dimension);
      Float64 WORK;
      Int32 LWORK = -1;
      Int32 INFO;

      // Call once to request optimal workspace size.
      dgeev_(&JOBVL, &JOBVR, &N, aTransposeColumnMajor.data(), &LDA,
             eigenvaluesReal.data(), eigenvaluesImag.data(),
             leftEigenvectorsTranspose.data(), &LDVL,
             rightEigenvectorsTranspose.data(), &LDVR,
             &WORK, &LWORK, &INFO);

      // Check for errors.
      if(INFO != 0L) {
        std::ostringstream message;
        message << "First call to dgeev_ returns " << INFO
                << ".  Something is wrong.";
        BRICK_THROW(brick::common::ValueException, "eigenvectors()",
                    message.str().c_str());
      }
    
      // Resize workspace.
      LWORK = static_cast<Int32>(WORK);
      Array1D<Float64> doubleWorkSpace(static_cast<size_t>(LWORK));

      // Call again to really compute the eigenvectors.
      dgeev_(&JOBVL, &JOBVR, &N, aTransposeColumnMajor.data(), &LDA,
             eigenvaluesReal.data(), eigenvaluesImag.data(),
             leftEigenvectorsTranspose.data(), &LDVL,
             rightEigenvectorsTranspose.data(), &LDVR,
             doubleWorkSpace.data(), &LWORK, &INFO);

      // Check for errors.
      if(INFO != 0L) {
        std::ostringstream message;
        message << "Second call to dgeev_ returns " << INFO
                << ".  Something is wrong.";
        BRICK_THROW(brick::common::ValueException, "eigenvectors()",
                    message.str().c_str());
      }

      // Check the sizes of the return references.
      if(eigenvalues.size() != dimension) {
        eigenvalues.reinit(dimension);
      }
      if(eigenvectors.rows() != dimension
         || eigenvectors.columns() != dimension) {
        eigenvectors.reinit(dimension, dimension);
      }

      // Copy eigenvalues and eigenvectors.
      for(size_t ii = 0; ii < dimension; ++ii) {
        eigenvalues[ii].real() = eigenvaluesReal[ii];
        eigenvalues[ii].imag() = eigenvaluesImag[ii];
        if(eigenvaluesImag[ii] == 0.0) {
          for(size_t jj = 0; jj < dimension; ++jj) {
            eigenvectors(jj, ii).real() = leftEigenvectorsTranspose(ii, jj);
            eigenvectors(jj, ii).imag() = 0.0;
          }
        } else {
          eigenvalues[ii + 1].real() = eigenvaluesReal[ii + 1];
          eigenvalues[ii + 1].imag() = eigenvaluesImag[ii + 1];
          for(size_t jj = 0; jj < dimension; ++jj) {
            eigenvectors(jj, ii).real() =
              leftEigenvectorsTranspose(ii, jj);
            eigenvectors(jj, ii).imag() =
              leftEigenvectorsTranspose(ii + 1, jj);
            eigenvectors(jj, ii + 1).real() =
              leftEigenvectorsTranspose(ii, jj);
            eigenvectors(jj, ii + 1).imag() =
              -leftEigenvectorsTranspose(ii + 1, jj);
          }
          ++ii;
        }
      }

      // Do eigenvalues have to be in descending order of magnitude.
      if(isSortRequired) {
        Array1D<double> magnitudes2(eigenvalues.size());
        for(size_t ii = 0; ii < eigenvalues.size(); ++ii) {
          magnitudes2[ii] = (eigenvalues[ii].real() * eigenvalues[ii].real()
                             + eigenvalues[ii].imag() * eigenvalues[ii].imag());
        }
        Array1D<size_t> outputOrder = argsort(magnitudes2);
        std::reverse(outputOrder.begin(), outputOrder.end());

        Array1D< std::complex<Float64> > sortedEigenvalues(eigenvalues.size());
        Array2D< std::complex<Float64> > sortedEigenvectors(
          eigenvectors.rows(), eigenvectors.columns());
        for(size_t ii = 0; ii < dimension; ++ii) {
          sortedEigenvalues[ii] = eigenvalues[outputOrder[ii]];
          for(size_t jj = 0; jj < dimension; ++jj) {
            sortedEigenvectors(jj, ii) = eigenvectors(jj, outputOrder[ii]);
          }
        }
        eigenvalues = sortedEigenvalues;
        eigenvectors = sortedEigenvectors;
      }
    }

    
    void
    eigenvectorsSymmetric(Array2D<Float64> const& inputArray,
                          Array1D<Float64>& eigenvalues,
                          Array2D<Float64>& eigenvectors)
    {
      // Argument checking.
      if(inputArray.size() == 0) {
        BRICK_THROW(brick::common::ValueException,
                    "eigenvectorsSymmetric()",
                    "Argument inputArray cannot have zero size.");
      }
      if(inputArray.rows() != inputArray.columns()) {
        BRICK_THROW(brick::common::ValueException,
                    "eigenvectorsSymmetric()",
                    "Argument inputArray must be square.");
      }
    
      // Transpose A to match LAPACK's convention.  Since inputArray is
      // symmetric, we don't really need to transpose!  Also, we only
      // need to copy the upper triangular part.
      size_t dimension = inputArray.rows();
      Array2D<Float64> aColumnMajor(dimension, dimension);
      {
        Array2D<Float64>::const_iterator inPtr = inputArray.begin();
        Array2D<Float64>::iterator outPtr = aColumnMajor.begin();
        for(size_t index0 = 0; index0 < dimension; ++index0) {
          inPtr += index0;
          outPtr += index0;
          for(size_t index1 = index0; index1 < dimension; ++index1) {
            *(outPtr++) = *(inPtr++);
          }
        }
      }

      // Set up storage for return values.
      Array1D<Float64> eigenvaluesTmp(dimension);

      // Set up arguments for the LAPACK call.
      char JOBZ = 'V';  // Compute eigenvalues and eigenvectors.
      char UPLO = 'L';  // Get input from lower triangle of A.
      Int32 N = static_cast<Int32>(dimension);
      Int32 LDA = static_cast<Int32>(dimension);
      Float64 WORK;
      Int32 LWORK = -1;
      Int32 INFO;

      // Call once to request optimal workspace size.
      dsyev_(&JOBZ,  &UPLO, &N, aColumnMajor.data(), &LDA,
             eigenvaluesTmp.data(), &WORK, &LWORK, &INFO);

      // Check for errors.
      if(INFO != 0L) {
        std::ostringstream message;
        message << "First call to dsyev_ returns " << INFO
                << ".  Something is wrong.";
        BRICK_THROW(brick::common::ValueException,
                    "eigenvectorsSymmetric()",
                    message.str().c_str());
      }
    
      // Resize workspace.
      LWORK = static_cast<Int32>(WORK);
      Array1D<Float64> doubleWorkSpace(static_cast<size_t>(LWORK));

      // Call again to really compute the eigenvectors.
      dsyev_(&JOBZ,  &UPLO, &N, aColumnMajor.data(), &LDA,
             eigenvaluesTmp.data(), doubleWorkSpace.data(), &LWORK, &INFO);

      // Check for errors.
      if(INFO != 0L) {
        std::ostringstream message;
        message << "Second call to dsyev_ returns " << INFO
                << ".  Something is wrong.";
        BRICK_THROW(brick::common::ValueException,
                    "eigenvectorsSymmetric()",
                    message.str().c_str());
      }

      // Recover the result.  Eigenvectors are left in aColumnMajor, but
      // unfortunately they're transposed and in the reverse order from
      // our convention.  How awkward.

      // Check the sizes of the return references.
      if(eigenvalues.size() != dimension) {
        eigenvalues.reinit(dimension);
      }
      if(eigenvectors.rows() != dimension
         || eigenvectors.columns() != dimension) {
        eigenvectors.reinit(dimension, dimension);
      }

      // Copy eigenvalues in reverse order.
      {
        Array1D<Float64>::iterator outPtr = eigenvalues.end() - 1;
        Array1D<Float64>::iterator inPtr0 = eigenvaluesTmp.begin();
        Array1D<Float64>::iterator inPtr1 = eigenvaluesTmp.end();
        while(inPtr0 != inPtr1) {
          *outPtr = *inPtr0;
          ++inPtr0;
          --outPtr;
        }
      }

      // Transpose and reverse and copy, oh my!
      for(size_t index0 = 0; index0 < dimension; ++index0) {
        Float64* outPtr = eigenvectors.data(index0);
        Float64* inPtr = aColumnMajor.data((dimension - index0 - 1) * dimension);
        for(size_t index1 = 0; index1 < dimension; ++index1) {
          *outPtr = *inPtr;
          ++inPtr;
          outPtr += dimension;
        }
      }
    }


    Array2D<Float64>
    inverse(Array2D<Float64> const& A)
    {
      // First argument checking.
      if(A.columns() != A.rows()) {
        BRICK_THROW(brick::common::ValueException,
                    "inverse(Array2D<Float64> const&)",
                    "Input array is not square.");
      }
      // Now set up some linear equations to solve.
      Array2D<Float64> AInverse = identity<Float64>(A.rows(), A.rows());
      Array2D<Float64> AA = A.transpose();
      // And solve for the inverse matrix.
      linearSolveInPlace(AA, AInverse); //Modifies AInverse.
      return AInverse;
    }


    // This function computes the best linear fit between the two input
    // arrays.
    std::pair<Float64, Float64>
    linearFit(Array1D<Float64> const& array0,
              Array1D<Float64> const& array1)
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

      // First some argument checking.
      if(array0.size() != array1.size()) {
        std::ostringstream message;
        message << "Arguments array0 and array1 must have the same size, "
                << "but are of size " << array0.size()
                << " and " << array1.size() << " respectively." << std::endl;
        BRICK_THROW(brick::common::ValueException,
                    "linearFit(Array1D<Float64> const&, Array1D<Float64> const&)",
                    message.str().c_str());                 
      }
      if(array0.size() == 0) {
        BRICK_THROW(brick::common::ValueException,
                    "linearFit(Array1D<Float64> const&, Array1D<Float64> const&)",
                    "Arguments cannot have zero size.");
      }

      // Do linear regression.
      Array2D<Float64> AMatrix(2, 2);
      AMatrix(0, 0) = dot<Float64>(array0, array0);
      AMatrix(0, 1) = sum<Float64>(array0);
      AMatrix(1, 0) = sum<Float64>(array0);
      AMatrix(1, 1) = array0.size();

      Array1D<Float64> bVector(2);
      bVector(0) = dot<Float64>(array0, array1);
      bVector(1) = sum<Float64>(array1);
    
      Array2D<Float64> AInverse = inverse(AMatrix);

      Array1D<Float64> result = matrixMultiply<Float64>(AInverse, bVector);
      return std::make_pair(result(0), result(1));
    }

  
    // This function solves the system of equations A*x = b, where A and
    // b are known Array2D<double> instances.
    Array1D<Float64>
    linearLeastSquares(Array2D<Float64> const& A, Array1D<Float64> const& b)
    {
      // First some argument checking.
      if(A.size() == 0) {
        BRICK_THROW(brick::common::
                    ValueException,
                    "linearLeastSquares(Array2D<Float64> const&, Array1D<Float64> const&)",
                    "Input array A must have nonzero size.");
      }
      if(A.rows() != b.size()) {
        BRICK_THROW(brick::common::
                    ValueException,
                    "linearLeastSquares(Array2D<Float64> const&, Array1D<Float64> const&)",
                    "The number of rows in input array A must be the same as the number "
                    "of elements in b.");
      }

      // This two-line implementation is slow.
      // Array2D<Float64> APInv = pseudoinverse(A);
      // return matrixMultiply<Float64>(APInv, b);

      // Set up scalar arguments for the LAPACK routine.
      char trans = 'N';
      Int32 rows = static_cast<Int32>(A.rows());
      Int32 columns = static_cast<Int32>(A.columns());
      Int32 nrhs = static_cast<Int32>(1);
      Int32 ldb = static_cast<Int32>(std::max(rows, columns));
      Float64 temporaryWorkspace;
      Int32 lwork = -1;  // Request info on optimal workspace size.
      Int32 info;

      // Set up array arguments for the LAPACK routine.
      Array2D<Float64> AColumnMajor = A.transpose();
      Array1D<Float64> bCopy(ldb);
      std::copy(b.begin(), b.end(), bCopy.begin());
    
      // Now invoke the LAPACK routine to find the optimal workspace
      // size.
      dgels_(&trans, &rows, &columns, &nrhs, AColumnMajor.data(), &rows,
             bCopy.data(), &ldb, &temporaryWorkspace, &lwork, &info);

      // Check for errors.
      if(info != 0L) {
        std::ostringstream message;
        message << "First call to dgels_ returns " << info
                << ".  Something is wrong.";
        BRICK_THROW(brick::common::ValueException,
                    "linearLeastSquares()",
                    message.str().c_str());
      }
    
      // Resize workspace.
      lwork = static_cast<Int32>(temporaryWorkspace);
      Array1D<Float64> doubleWorkspace(static_cast<size_t>(lwork));

      // Call again to solve the system of equations.
      dgels_(&trans, &rows, &columns, &nrhs, AColumnMajor.data(), &rows,
             bCopy.data(), &ldb, doubleWorkspace.data(), &lwork, &info);
    
      // Check for errors.
      if(info != 0L) {
        std::ostringstream message;
        message << "Second call to dgels_ returns " << info
                << ".  Something is wrong.";
        BRICK_THROW(brick::common::ValueException,
                    "linearLeastSquares()",
                    message.str().c_str());
      }

      if(bCopy.size() == A.columns()) {
        return bCopy;
      } else {
        Array1D<Float64> returnValue(A.columns());
        std::copy(bCopy.begin(), bCopy.begin() + A.columns(),
                  returnValue.begin());
        return returnValue;
      }
    }


    // WARNING:  linearSolveInPlace() destructively modifies 
    // both arguments!
    //
    // This function solves the system of equations A*x = b, where A is
    // a known matrix, and b is a known vector.
    void
    linearSolveInPlace(Array2D<Float64>& A, Array1D<Float64>& b)
    {
      Array2D<Float64> bMatrix(b.size(), 1, b.data());
      linearSolveInPlace(A, bMatrix);
    }


    // This function is identical to linearSolveInPlace(Array2D<Float64>,
    // Array1D<Float64>&), except that b (and therefore x) is not
    // constrained to be a vector.
    void
    linearSolveInPlace(Array2D<Float64> &A, Array2D<Float64> &b)
    {
      // First some argument checking.
      if(A.rows() != b.rows()) {
        BRICK_THROW(brick::common::ValueException,
                    "linearSolveInPlace(Array2D<Float64>&, Array2D<Float64>&)",
                    "Input arrays A and b must have the same number of rows.");
      }
      if(A.rows() != A.columns()) {
        BRICK_THROW(brick::common::ValueException,
                    "linearSolveInPlace(Array2D<Float64>&, Array2D<Float64>&)",
                    "Input array A must be square.");
      }

      // Grr.  Have to transpose to match lapack.
      Array2D<Float64> AColumnMajor = A.transpose();
    
      // Now invoke the LAPACK routine.
      Int32 rows = static_cast<Int32>(A.rows());
      Int32 xColumns = static_cast<Int32>(b.columns());
      Int32 *iPiv = new Int32[rows];
      Int32 info;
      dgesv_(&rows, &xColumns, AColumnMajor.data(), &rows, iPiv, b.data(), &rows,
             &info );
      // Clean up
      delete[] iPiv;
      // And do some error checking.
      if(info != 0L) {
        std::ostringstream message;
        message << "Call to dgesv_ returns " << info << ".  Something is wrong."
                << "  Perhaps the the input equations are poorly conditioned "
                << "or perhaps there is no solution.";
        BRICK_THROW(brick::common::ValueException,
                    "linearSolveInPlace(Array2D<Float64>&, Array2D<Float64>&)",
                    message.str().c_str());
      }
    }


    // This function solves the system of equations A*x = b, where A is
    // a known tridiagonal matrix and b is a known vector.
    Array1D<Float64>
    linearSolveTridiagonal(Array1D<Float64> const& subDiagonal,
                           Array1D<Float64> const& centerDiagonal,
                           Array1D<Float64> const& superDiagonal,
                           Array1D<Float64> const& bVector)
    {
      // First some argument checking.
      if(subDiagonal.size() != superDiagonal.size()) {
        BRICK_THROW(brick::common::ValueException,
                    "linearSolveTridiagonal()",
                    "Input arguments subDiagonal and superDiagonal "
                    "must have the same size.");
      }
      if(centerDiagonal.size() != bVector.size()) {
        BRICK_THROW(brick::common::ValueException,
                    "linearSolveTridiagonal()",
                    "Input arguments centerDiagonal and bVector "
                    "must have the same size.");
      }
      if(centerDiagonal.size() != subDiagonal.size() + 1) {
        BRICK_THROW(brick::common::ValueException,
                    "linearSolveTridiagonal()",
                    "Input argument centerDiagonal must have one more element "
                    "than input argument subDiagonal.");
      }

    
      // Now invoke the LAPACK routine.
      Int32 N = static_cast<Int32>(centerDiagonal.size());
      Int32 NRHS = 1;
      Array1D<Float64> subDiagonalCopy = subDiagonal.copy();
      Array1D<Float64> centerDiagonalCopy = centerDiagonal.copy();
      Array1D<Float64> superDiagonalCopy = superDiagonal.copy();
      Array1D<Float64> xVector = bVector.copy();
      Int32 INFO;

      dgtsv_(&N, &NRHS, subDiagonalCopy.data(), centerDiagonalCopy.data(),
             superDiagonalCopy.data(), xVector.data(), &N, &INFO);

      // And do some error checking.
      if(INFO != 0L) {
        std::ostringstream message;
        message << "Call to dgtsv_ returns " << INFO << ".  Something is wrong."
                << "  Perhaps the the input equations are poorly conditioned "
                << "or perhaps there is no solution.";
        BRICK_THROW(brick::common::ValueException,
                    "linearSolveTridiagonal()",
                    message.str().c_str());
      }

      return xVector;
    }
    

    // This function computes the QR factorization of a general
    // matrix.
    void
    qrFactorization(Array2D<Float64> const& inputArray,
                    Array2D<Float64>& qArray,
                    Array2D<Float64>& rArray)
    {
      // Argument checking.
      if(inputArray.size() == 0) {
        qArray.clear();
        rArray.clear();
        return;
      }

      // Set up scalar arguments for the LAPACK routine.
      Int32 rows = static_cast<Int32>(inputArray.rows());
      Int32 columns = static_cast<Int32>(inputArray.columns());
      Int32 lda = rows;
      Float64 temporaryWorkspace;
      Int32 lwork = -1;  // Request info on optimal workspace size.
      Int32 info;

      // Set up array arguments for the LAPACK routine.
      Array2D<Float64> AColumnMajor = inputArray.transpose();
      Array1D<Float64> tauArray = Array1D<Float64>(std::min(rows, columns));
    
      // Now invoke the LAPACK routine to find the optimal workspace
      // size.
      dgeqrf_(&rows, &columns, AColumnMajor.data(), &lda, tauArray.data(),
              &temporaryWorkspace, &lwork, &info);

      // Check for errors.
      if(info != 0L) {
        std::ostringstream message;
        message << "First call to dgeqrf_ returns " << info
                << ".  Something is wrong.";
        BRICK_THROW(brick::common::ValueException, "qrFactorization()", message.str().c_str());
      }
    
      // Resize workspace.
      lwork = static_cast<Int32>(temporaryWorkspace);
      Array1D<Float64> doubleWorkspace(static_cast<size_t>(lwork));

      // Call again to solve the system of equations.
      dgeqrf_(&rows, &columns, AColumnMajor.data(), &lda, tauArray.data(),
              doubleWorkspace.data(), &lwork, &info);
    
      // Check for errors.
      if(info != 0L) {
        std::ostringstream message;
        message << "Second call to dgeqrf_ returns " << info
                << ".  Something is wrong.";
        BRICK_THROW(brick::common::ValueException, "qrFactorization()", message.str().c_str());
      }

      // Prepare reference arguments to receive the result.
      if((rArray.rows() != inputArray.rows())
         || (rArray.columns() != inputArray.columns())) {
        rArray.reinit(inputArray.rows(), inputArray.columns());
      }
      rArray = 0.0;
      if((qArray.rows() != inputArray.rows())
         || (qArray.columns() != inputArray.rows())) {
        qArray.reinit(inputArray.rows(), inputArray.rows());
      }
      qArray = 0.0;

      // Recover the upper trapezoidal matrix R.
      size_t numNonzeroRows = std::min(inputArray.rows(), inputArray.columns());
      for(size_t row = 0; row < numNonzeroRows; ++row) {
        for(size_t column = row; column < rArray.columns(); ++column) {
          rArray(row, column) = AColumnMajor(column, row);
        }
      }

      // Recover the orthogonal matrix Q.  Here's a quote from the lapack code:
      //
      // *  The matrix Q is represented as a product of elementary reflectors
      // *
      // *     Q = H(1) H(2) . . . H(k), where k = min(m,n).
      // *
      // *  Each H(i) has the form
      // *
      // *     H(i) = I - tau * v * v'
      // *
      // *  where tau is a real scalar, and v is a real vector with
      // *  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
      // *  and tau in TAU(i).

      // Start with identity.
      for(size_t row = 0; row < qArray.rows(); ++row) {
        qArray(row, row) = static_cast<Float64>(1.0);
      }
      
      // Now subtract out the tau * v * v^T terms.
      Array1D<Float64> vArray(qArray.rows());
      Array2D<Float64> reflector(qArray.rows(), qArray.rows());
      for(size_t ii = 0; ii < numNonzeroRows; ++ii) {
        // Get ready to compute the next factor.  Could do this more
        // efficiently.
        reflector = 0.0;
        for(size_t kk = 0; kk < reflector.rows(); ++kk) {
          reflector(kk, kk) = 1.0;
        }
        
        // Compute tau and v.
        Float64 tau = tauArray[ii];
        // We don't need a loop here because the zeros from last
        // iteration stick around.
        // 
        // for(size_t jj = 0; jj < ii; ++jj) {
        //   vArray[jj] = 0.0;
        // }
        if(ii > 0) {
          vArray[ii - 1] = 0.0;
        }
        vArray[ii] = 1.0;
        for(size_t jj = ii + 1; jj < qArray.rows(); ++jj) {
          vArray[jj] = AColumnMajor(ii, jj);
        }

        // Subtract tau * v * v^T.
        for(size_t row = ii; row < qArray.rows(); ++row) {
          for(size_t column = ii; column < qArray.columns(); ++column) {
            reflector(row, column) -= tau * vArray[row] * vArray[column];
          }
        }

        qArray = matrixMultiply<Float64>(qArray, reflector);
      }

      // Make sure diagonal elements of rArray are non-negative, as
      // promised.
      for(size_t ii = 0; ii < inputArray.rows(); ++ii) {
        if(rArray(ii, ii) == 0.0) {
          continue;
        }
        if(rArray(ii, ii) < 0.0) {
          for(size_t jj = 0; jj < inputArray.columns(); ++jj) {
            rArray(ii, jj) *= -1.0;
          }
          for(size_t jj = 0; jj < qArray.rows(); ++jj) {
            qArray(jj, ii) *= -1.0;
          }
        }
      }
    }

    
    // This function accepts an Array2D<Float64> instance having at least
    // as many rows as columns, and returns the Moore-Penrose
    // pseudoinverse.
    Array2D<Float64>
    pseudoinverse(Array2D<Float64> const& A)
    {
      Array2D<Float64> ATranspose = A.transpose();
      Array2D<Float64> ATA = matrixMultiply<Float64>(ATranspose, A);
      return matrixMultiply<Float64>(inverse(ATA), ATranspose);
    }


    // This function computes the singular value decomposition of a
    // matrix.
    void
    singularValueDecomposition(Array2D<Float64> const& inputArray,
                               Array2D<Float64>& uArray,
                               Array1D<Float64>& sigmaArray,
                               Array2D<Float64>& vTransposeArray,
                               bool isNullSpaceRequired,
                               bool isFullRangeRequired)
    {
      // Argument checking.
      if(inputArray.size() == 0) {
        BRICK_THROW(brick::common::ValueException,
                    "singularValueDecomposition()",
                    "Argument inputArray cannot have zero size.");
      }

      // Note(xxx): Fix things so that isFullRangeRequired is not ignored.
      isFullRangeRequired = isNullSpaceRequired;
      
      // 
      // We do not transpose A. Instead, we let LAPACK operate on the
      // untransposed matrix as if it were column major and then swap U
      // and VT at the end.  We make a copy here since the contents of
      // whatever array we pass to LAPACK will be destroyed.
      Array2D<Float64> aColumnMajor = inputArray.copy();

      // Set up storage for return values.
      size_t numberOfSingularValues =
        std::min(inputArray.rows(), inputArray.columns());

      size_t numberOfVRows;
      if(isNullSpaceRequired) {
        numberOfVRows = inputArray.columns();
      } else {
        numberOfVRows = numberOfSingularValues;
      }
      if((vTransposeArray.rows() != numberOfVRows)
         || (vTransposeArray.columns() != inputArray.columns())) {
        vTransposeArray.reinit(numberOfVRows, inputArray.columns());
      }
      // Shallow copy.
      Array2D<Float64> uColumnMajor = vTransposeArray;

      size_t numberOfUColumns;
      if(isFullRangeRequired) {
        numberOfUColumns = inputArray.rows();
      } else {
        numberOfUColumns = numberOfSingularValues;
      }
      if((uArray.rows() != inputArray.rows())
         || (uArray.columns() != numberOfUColumns)) {
        uArray.reinit(inputArray.rows(), numberOfUColumns);
      }
      // Shallow copy.
      Array2D<Float64> vTColumnMajor = uArray;

      if(sigmaArray.size() != numberOfSingularValues) {
        sigmaArray.reinit(numberOfSingularValues);
      }
      Array1D<Int32> integerWorkSpace(8 * numberOfSingularValues);
    
      // Set up arguments for the LAPACK call.
      char JOBZ;
      if(isNullSpaceRequired) {
        JOBZ = 'A';
      } else {
        JOBZ = 'S';  // Compute only min(M, N) columns of U and rows of VT.
      }
      Int32 M = static_cast<Int32>(inputArray.columns());
      Int32 N = static_cast<Int32>(inputArray.rows());
      Int32 LDA = M;
      Int32 LDU = M;
      Int32 LDVT = static_cast<Int32>(numberOfSingularValues);
      Float64 WORK;
      Int32 LWORK = -1;  // Request info on optimal workspace size.
      Int32 INFO;

      // Call once to request optimal workspace size.
      dgesdd_(&JOBZ, &M, &N, aColumnMajor.data(), &LDA,
              sigmaArray.data(), uColumnMajor.data(), &LDU,
              vTColumnMajor.data(), &LDVT,
              &WORK, &LWORK,
              integerWorkSpace.data(), &INFO);

      // Check for errors.
      if(INFO != 0L) {
        std::ostringstream message;
        message << "First call to dgesdd_ returns " << INFO
                << ".  Something is wrong.";
        BRICK_THROW(brick::common::ValueException,
                    "singularValueDecomposition()",
                    message.str().c_str());
      }
    
      // Resize workspace.
      LWORK = static_cast<Int32>(WORK);
      Array1D<Float64> doubleWorkSpace(static_cast<size_t>(LWORK));

      // Call again to do the SVD.
      dgesdd_(&JOBZ, &M, &N, aColumnMajor.data(), &LDA,
              sigmaArray.data(), uColumnMajor.data(), &LDU,
              vTColumnMajor.data(), &LDVT,
              doubleWorkSpace.data(), &LWORK,
              integerWorkSpace.data(), &INFO);

      // Check for errors.
      if(INFO != 0L) {
        std::ostringstream message;
        message << "Second call to dgesdd_ returns " << INFO
                << ".  Something is wrong.";
        BRICK_THROW(brick::common::ValueException,
                    "singularValueDecomposition()",
                    message.str().c_str());
      }
    } 



    void
    singularValueDecompositionSimple(Array2D<Float64> const& inputArray,
                                     Array2D<Float64>& uArray,
                                     Array1D<Float64>& sigmaArray,
                                     Array2D<Float64>& vTransposeArray)
    {
      // Argument checking.
      if(inputArray.size() == 0) {
        BRICK_THROW(brick::common::ValueException,
                    "singularValueDecomposition()",
                    "Argument inputArray cannot have zero size.");
      }
    
      // We do not transpose A. Instead, we let LAPACK operate on the
      // untransposed matrix as if it were column major and then swap U
      // and VT at the end.  We make a copy here since the contents of
      // whatever array we pass to LAPACK will be destroyed.
      Array2D<Float64> aColumnMajor = inputArray.copy();

      // Set up storage for return values.
      size_t numberOfSingularValues =
        std::min(inputArray.rows(), inputArray.columns());

      // Make sure V is appropriately sized.
      if((vTransposeArray.rows() != numberOfSingularValues)
         || (vTransposeArray.columns() != inputArray.columns())) {
        vTransposeArray.reinit(numberOfSingularValues, inputArray.columns());
      }

      // Shallow copy.
      Array2D<Float64> uColumnMajor = vTransposeArray;

      // Make sure U is appropriately sized.
      if((uArray.rows() != inputArray.rows())
         || (uArray.columns() != numberOfSingularValues)) {
        uArray.reinit(inputArray.rows(), numberOfSingularValues);
      }
    
      // Shallow copy.
      Array2D<Float64> vTColumnMajor = uArray;

      if(sigmaArray.size() != numberOfSingularValues) {
        sigmaArray.reinit(numberOfSingularValues);
      }

      // Integer workspace.
      Array1D<Int32> integerWorkspace(8 * numberOfSingularValues);

      // Set up arguments for the LAPACK call.
      char JOBZ = 'S';  // Compute only min(M, N) columns of U.
      Int32 M = static_cast<Int32>(inputArray.columns());
      Int32 N = static_cast<Int32>(inputArray.rows());
      Int32 LDA = M;
      Int32 LDU = M;
      Int32 LDVT = static_cast<Int32>(numberOfSingularValues);
      Float64 WORK;
      Int32 LWORK = -1;  // Request info on optimal workspace size.
      Int32 INFO;

      // Call once to request optimal workspace size.
      dgesdd_(&JOBZ, &M, &N, aColumnMajor.data(), &LDA,
              sigmaArray.data(), uColumnMajor.data(), &LDU,
              vTColumnMajor.data(), &LDVT,
              &WORK, &LWORK, integerWorkspace.data(), &INFO);

      // Check for errors.
      if(INFO != 0L) {
        std::ostringstream message;
        message << "First call to dgesdd_ returns " << INFO
                << ".  Something is wrong.";
        BRICK_THROW(brick::common::ValueException,
                    "singularValueDecomposition()",
                    message.str().c_str());
      }
    
      // Resize workspace.
      LWORK = static_cast<Int32>(WORK);
      Array1D<Float64> doubleWorkSpace(static_cast<size_t>(LWORK));

      // Call again to do the SVD.
      dgesdd_(&JOBZ, &M, &N, aColumnMajor.data(), &LDA,
              sigmaArray.data(), uColumnMajor.data(), &LDU,
              vTColumnMajor.data(), &LDVT,
              doubleWorkSpace.data(), &LWORK, integerWorkspace.data(), &INFO);

      // Check for errors.
      if(INFO != 0L) {
        std::ostringstream message;
        message << "Second call to dgesdd_ returns " << INFO
                << ".  Something is wrong.";
        BRICK_THROW(brick::common::ValueException,
                    "singularValueDecomposition()",
                    message.str().c_str());
      }
    } 

  
//   // This function computes the singular value decomposition of a
//   // matrix, but is less efficient than the routine above: It
//   // explicitly transposes the matrix before passing it to LAPACK.
//   void
//   singularValueDecomposition(
//     Array2D<double> const& inputArray,
//     Array2D<double>& uArray,
//     Array1D<double>& sigmaArray,
//     Array2D<double>& vTransposeArray)
//   {
//     // Argument checking.
//     if(inputArray.size() == 0) {
//       BRICK_THROW(brick::common::ValueException,
//                  "singularValueDecomposition()",
//                  "Argument inputArray cannot have zero size.");
//     }
    
//     // Transpose A to match LAPACK's convention.
//     Array2D<double> aColumnMajor = inputArray.transpose();

//     // Set up storage for return values.
//     size_t numberOfSingularValues =
//       std::min(inputArray.rows(), inputArray.columns());
//     Array2D<double> uColumnMajor(numberOfSingularValues, inputArray.rows());
//     Array2D<double> vTColumnMajor(
//       inputArray.columns(), numberOfSingularValues);
//     if(sigmaArray.size() != numberOfSingularValues) {
//       sigmaArray.reinit(numberOfSingularValues);
//     }
//     Array1D<Int32> integerWorkSpace(8 * numberOfSingularValues);
    
//     // Set up arguments for the LAPACK call.
//     char JOBZ = 'S';  // Compute only min(M, N) columns of U and rows of VT.
//     Int32 M = static_cast<Int32>(inputArray.rows());
//     Int32 N = static_cast<Int32>(inputArray.columns());
//     Int32 LDA = M;
//     Int32 LDU = M;
//     Int32 LDVT = static_cast<Int32>(numberOfSingularValues);
//     double WORK;
//     Int32 LWORK = -1;  // Request info on optimal workspace size.
//     Int32 INFO;

//     // Call once to request optimal workspace size.
//     dgesdd_(&JOBZ, &M, &N, aColumnMajor.data(), &LDA,
//             sigmaArray.data(), uColumnMajor.data(), &LDU,
//             vTColumnMajor.data(), &LDVT,
//             &WORK, &LWORK,
//             integerWorkSpace.data(), &INFO);

//     // Check for errors.
//     if(INFO != 0L) {
//       std::ostringstream message;
//       message << "First call to dgesdd_ returns " << INFO
//               << ".  Something is wrong.";
//       BRICK_THROW(brick::common::ValueException,
//                  "singularValueDecomposition()",
//                  message.str().c_str());
//     }
    
//     // Resize workspace.
//     LWORK = static_cast<Int32>(WORK);
//     Array1D<double> doubleWorkSpace(static_cast<size_t>(LWORK));

//     // Call again to do the SVD.
//     dgesdd_(&JOBZ, &M, &N, aColumnMajor.data(), &LDA,
//             sigmaArray.data(), uColumnMajor.data(), &LDU,
//             vTColumnMajor.data(), &LDVT,
//             doubleWorkSpace.data(), &LWORK,
//             integerWorkSpace.data(), &INFO);

//     // Check for errors.
//     if(INFO != 0L) {
//       std::ostringstream message;
//       message << "Second call to dgesdd_ returns " << INFO
//               << ".  Something is wrong.";
//       BRICK_THROW(brick::common::ValueException,
//                  "singularValueDecomposition()",
//                  message.str().c_str());
//     }

//     // Recover the result.
//     uArray = uColumnMajor.transpose();
//     vTransposeArray = vTColumnMajor.transpose();
//   } 
  

    // This function computes the singular values a matrix without
    // computing the associated U and V matrices.
    Array1D<Float64>
    singularValues(Array2D<Float64> const& inputArray)
    {
      // Argument checking.
      if(inputArray.size() == 0) {
        BRICK_THROW(brick::common::ValueException,
                    "singularValues()",
                    "Argument inputArray cannot have zero size.");
      }
    
      // Transpose A to match LAPACK's convention.
      //
      // Actually, it's slightly quicker not to transpose, and singular
      // values remain the same, so we just copy.
      // 
      // Array2D<Float64> aColumnMajor = inputArray.transpose();
      Array2D<Float64> aColumnMajor = inputArray.copy();

      // Set up storage for return values.
      size_t numberOfSingularValues =
        std::min(inputArray.rows(), inputArray.columns());
      Array1D<Float64> sigmaArray(numberOfSingularValues);
      Array1D<Int32> integerWorkSpace(8 * numberOfSingularValues);
    
      // Set up arguments for the LAPACK call.
      char JOBZ = 'N';  // Compute no columns of U and no rows of VT.
      // Since we didn't transpose, we have to reverse rows & columns.
      // Int32 M = static_cast<Int32>(inputArray.rows());
      // Int32 N = static_cast<Int32>(inputArray.columns());
      Int32 M = static_cast<Int32>(inputArray.columns());
      Int32 N = static_cast<Int32>(inputArray.rows());
      Int32 LDA = M;
      Float64 U;         // U is not referenced by the LAPACK call.
      Int32 LDU = 1;
      Float64 VT;        // VT is not referenced by the LAPACK call.
      Int32 LDVT = 1;
      Float64 WORK;
      Int32 LWORK = -1;  // Request info on optimal workspace size.
      Int32 INFO;

      // Call once to request optimal workspace size.
      dgesdd_(&JOBZ, &M, &N, aColumnMajor.data(), &LDA,
              sigmaArray.data(), &U, &LDU, &VT, &LDVT,
              &WORK, &LWORK,
              integerWorkSpace.data(), &INFO);

      // Check for errors.
      if(INFO != 0L) {
        std::ostringstream message;
        message << "First call to dgesdd_ returns " << INFO
                << ".  Something is wrong.";
        BRICK_THROW(brick::common::ValueException,
                    "singularValues()",
                    message.str().c_str());
      }
    
      // Resize workspace.
      LWORK = static_cast<Int32>(WORK);

      // Sanity check to correct a bug in LAPACK3.
      Int32 minimumLWORK =
        3 * std::min(M,N) + std::max(std::max(M, N), 6 * std::min(M, N));
      if(LWORK < minimumLWORK) {
        // std::cerr << "Warning: changing LWORK from " << LWORK << " to "
        //           << minimumLWORK << "." << std::endl;
        LWORK = minimumLWORK;
      }
      Array1D<Float64> doubleWorkSpace(static_cast<size_t>(LWORK));
      
      // Call again to do the SVD.
      dgesdd_(&JOBZ, &M, &N, aColumnMajor.data(), &LDA,
              sigmaArray.data(), &U, &LDU, &VT, &LDVT,
              doubleWorkSpace.data(), &LWORK,
              integerWorkSpace.data(), &INFO);

      // Check for errors.
      if(INFO != 0L) {
        std::ostringstream message;
        message << "Second call to dgesdd_ returns " << INFO
                << ".  Something is wrong.";
        BRICK_THROW(brick::common::ValueException,
                    "singularValues()",
                    message.str().c_str());
      }

      // Return the result.
      return sigmaArray;
    }
  
  } // namespace linearAlgebra
  
} // namespace brick
