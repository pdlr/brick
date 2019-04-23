/**
***************************************************************************
* @file brick/linearAlgebra/linearAlgebra.hh
*
* Header file declaring linear algebra functions.  Many of these depend
* on the LAPACK and BLAS libraries.
*
* Copyright (C) 2001-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
***************************************************************************
**/

#ifndef BRICK_LINEARALGEBRA_LINEARALGEBRA_HH
#define BRICK_LINEARALGEBRA_LINEARALGEBRA_HH

#include <complex>
#include <brick/common/types.hh>
#include <brick/numeric/array1D.hh>
#include <brick/numeric/array2D.hh>

namespace brick {

  namespace linearAlgebra {

    /**
     ** This exception class is thrown when the rank of a matrix is
     ** not as expected.
     **/
    class RankException;  // Forward declaration to help Doxygen.
    BRICK_DECLARE_EXCEPTION_TYPE(RankException, brick::common::ValueException);


    /**
     * This function computes the Cholesky factorization of a
     * symmetric, positive definite matrix.  That is, for a symmetric,
     * positive definite matrix E, it computes the upper triangular
     * matrix K such that E == K^T * K (if argument isUpperTriangular
     * is true), or the lower triangular matrix such that E = K * K^T
     * (if argument isUpperTriangular is false).
     *
     * @param inputArray This argument is the matrix to be factored.
     *
     * @param kArray This argument will be filled in with the upper
     * triangular matrix K.
     *
     * @param isUpperTriangular This argument specifies whether the
     * result should be returned as an upper triangular or lower
     * triangular matrix.
     */
    void
    choleskyFactorization(
      brick::numeric::Array2D<brick::common::Float64> const& inputArray,
      brick::numeric::Array2D<brick::common::Float64>& kArray,
      bool isUpperTriangular = true);


    /**
     * This function computes the determinant of a square
     * Array2D<Float64> instance.
     *
     * @param A This argument represents the matrix from which to
     * compute the determinant.
     * @return The return value is a Float64 representing the determinant
     * of the input argument.
     */
    brick::common::Float64
    determinant(brick::numeric::Array2D<brick::common::Float64> const& A);


    /**
     * This function computes the eigenvalues of a symmetric real matrix.
     *
     * @param inputArray This argument is an Array2D<Float64> instance
     * representing a symmetric matrix.  It must be square and have
     * non-zero size.  Only the data in the upper triangular portion of
     * the array (including the diagonal) will be examined.  The
     * remaining elements are assumed to be symmetric.
     *
     * @return The return value is an Array1D<Float64> instance
     * containing the eigenvalues, sorted into descending order.
     */
    brick::numeric::Array1D<brick::common::Float64>
    eigenvaluesSymmetric(
      brick::numeric::Array2D<brick::common::Float64> const& inputArray);


    /**
     * This function computes the eigenvalues and eigenvectors of a
     * real matrix.
     *
     * @param inputArray This argument is an Array2D<Float64> instance
     * representing input matrix.  It must be square and have non-zero
     * size.
     *
     * @param eigenvalues This argument is used to return an array of
     * (possibly complex) eigenvalues.
     *
     * @param eigenvectors This argument is used to return the
     * eigenvectors.  On return, the first column contains the
     * eigenvector corresponding to the first eigenvalue, the second
     * column contains the eigenvector corresponding to the second
     * eigenvalue, and so on.
     *
     * @param isSortRequired Setting this argument to true will cause
     * the eigenvalues to sorted in decreasing order of magnitude.
     * The eigenvector array will also be manipulated so the i^th
     * column of eigenvectors still corresponds to the i^th element of
     * eigenvalues.
     */
    void
    eigenvectors(
      brick::numeric::Array2D<brick::common::Float64> const& inputArray,
      brick::numeric::Array1D< std::complex<brick::common::Float64> >& eigenvalues,
      brick::numeric::Array2D< std::complex<brick::common::Float64> >& eigenvectors,
      bool isSortRequired = false);


    /**
     * This function computes the eigenvalues and eigenvectors of a
     * symmetric real matrix.
     *
     * @param inputArray This argument is an Array2D<Float64> instance
     * representing a symmetric matrix.  It must be square and have
     * non-zero size.  Only the data in the upper triangular portion of
     * the array (including the diagonal) will be examined.  The
     * remaining elements are assumed to be symmetric.
     *
     * @param eigenvalues This argument is used to return an
     * Array1D<Float64> instance containing the eigenvalues, sorted into
     * descending order.
     *
     * @param eigenvectors This argument is used to return the
     * eigenvectors.  On return, the first column contains the
     * eigenvector corresponding to the first eigenvalue, the second
     * column contains the eigenvector corresponding to the second
     * eigenvalue, and so on.
     */
    void
    eigenvectorsSymmetric(
      brick::numeric::Array2D<brick::common::Float64> const& inputArray,
      brick::numeric::Array1D<brick::common::Float64>& eigenvalues,
      brick::numeric::Array2D<brick::common::Float64>& eigenvectors);


    /**
     * This function accepts a square Array2D<Float64> instance and
     * returns an Array2D<Float64> instance such that the matrix product
     * of the two is equal to the identity matrix.  It is an error if
     * the argument is not invertible.
     *
     * @param A This argument is the matrix to be inverted.
     * @return The return value is the inverse of argument A.
     */
    brick::numeric::Array2D<brick::common::Float64>
    inverse(brick::numeric::Array2D<brick::common::Float64> const& A);


    /**
     * This function computes the best linear fit between the two input
     * arrays.  That is, it solves for scalars a and b that minimize
     * the quantity
     *
     *   e = |(a * array0 + b) - array1|^2,
     *
     * where array0 and array1 are the two input arrays.  Put another
     * way, this means we're solving (in the least squares sense) for
     * the variables a and b that most nearly satisfy the following
     * equation:
     *
     *   a * array0 + b = array1
     *
     * @param array0 This argument is the first input array.
     * @param array1 This argument is the second input array.
     * @return The return value is a pair of Float64s in which the first
     * element is the variable a and the second is the variable b.
     */
    template <class FloatType>
    std::pair<FloatType, FloatType>
    linearFit(brick::numeric::Array1D<FloatType> const& array0,
              brick::numeric::Array1D<FloatType> const& array1);


    /**
     * This function solves the system of equations A*x = b, where A and
     * b are known Array2D<Float64> instances.  The contents of the
     * arguments are not modified.  If the solution fails, a
     * ValueException will be generated.
     *
     * @param A This argument specifies the A matrix in the system "Ax =
     * b," and must either be square or have more rows than columns.
     *
     * @param b This argument specifies the b vector in the system "Ax =
     * b."  It must have the same number of elements as argument A has
     * rows.
     *
     * @return The return value is the vector x that most nearly
     * satisfies the equation.  "Nearly" is defined in the least-squares
     * sense.
     */
    brick::numeric::Array1D<brick::common::Float64>
    linearLeastSquares(brick::numeric::Array2D<brick::common::Float64> const& A,
                       brick::numeric::Array1D<brick::common::Float64> const& b);


    /**
     * This function solves the system of equations A*x = b, where A is
     * a known matrix, and b is a known vector.  The contents of both
     * arguments are modified as part of the process.  If the solution
     * fails, a ValueException will be generated.
     *
     * @param A This argument specifies the A matrix in the system "Ax =
     * b," and must be square.
     *
     * @param b This argument specifies the b vector in the system "Ax =
     * b."  It must have the same size as x, and will be replaced with
     * the recovered value of x.
     */
    void
    linearSolveInPlace(brick::numeric::Array2D<brick::common::Float64>& A,
                       brick::numeric::Array1D<brick::common::Float64>& b);


    /**
     * This function is identical to linearSolveInPlace(Array2D<Float64>,
     * Array1D<Float64>&), except that b (and therefore x) is not
     * constrained to be a vector.
     *
     * @param A This argument specifies the A matrix in the system "Ax =
     * b," and must be square.  In the current implementation, it will
     * be replaced the LU decomposition of A, however this behavior may
     * change in the future.
     *
     * @param b This argument specifies the b matrix in the system "Ax =
     * b."  It must have the same size as x, and will be replaced with
     * the recovered value of x.
     */
    void
    linearSolveInPlace(brick::numeric::Array2D<brick::common::Float64>& A,
                       brick::numeric::Array2D<brick::common::Float64>& b);


    /**
     * This function solves the system of equations A*x = b, where A is
     * a known tridiagonal matrix and b is a known vector.  The contents
     * of the arguments are not modified.  If the solution fails, a
     * ValueException will be generated.
     *
     * @param subDiagonal This argument specifies the lower diagonal
     * of the A matrix in the system "A*x = b."
     *
     * @param centerDiagonal This argument specifies the center diagonal
     * of the A matrix in the system "A*x = b."  It's length must be 1
     * greater than the length of argument subDiagonal.
     *
     * @param superDiagonal This argument specifies the upper diagonal
     * of the A matrix in the system "A*x = b."  It's length must be the
     * same as the length of argument subDiagonal.
     *
     * @param b This argument specifies the b vector in the system "Ax =
     * b."  It must have the same length as argument centerDiagonal.
     *
     * @return The return value is the vector x that most nearly
     * satisfies the equation.  "Nearly" is defined in the least-squares
     * sense.
     */
    brick::numeric::Array1D<brick::common::Float64>
    linearSolveTridiagonal(
      brick::numeric::Array1D<brick::common::Float64> const& subDiagonal,
      brick::numeric::Array1D<brick::common::Float64> const& centerDiagonal,
      brick::numeric::Array1D<brick::common::Float64> const& superDiagonal,
      brick::numeric::Array1D<brick::common::Float64> const& bVector);


    /**
     * This function accepts an Array2D<Float64> instance having at least
     * as many rows as columns, and returns the Moore-Penrose
     * pseudoinverse.  That is, for an input matrix A, it returns
     * inverse(matrixMultiply(transpose(A), A)).  It is an error if the
     * rank of A is less than A.columns().
     *
     * @param A This argument is the matrix to be psuedoinverted.
     *
     * @return The return value is the pseudoinverse of argument A.
     */
    brick::numeric::Array2D<brick::common::Float64>
    pseudoinverse(brick::numeric::Array2D<brick::common::Float64> const& A);


    /**
     * This function computes the QR factorization of a general
     * matrix.  That is, for an MxN matrix A (M rows and N columns),
     * it computes the MxM orthogonal matrix Q (that is, Q^T * Q == I)
     * and the MxN upper trapezoidal (upper triangular M >= N) matrix
     * R such that A = Q*R, and the diagonal elements of R are
     * non-negative.
     *
     * @param inputArray This argument is the matrix to be factored.
     *
     * @param qArray This argument will be filled in with the
     * orthogonal matrix Q.  If qArray is already MxM, its
     * contents will be overwritten.  If qArray is not MxM, it will be
     * resized (i.e., new memory will be allocated).
     *
     * @param rArray This argument will be filled in with the
     * upper trapezoidal matrix R.  If rArray is already MxN, its
     * contents will be overwritten.  If rArray is not MxN, it will be
     * resized (i.e., new memory will be allocated).
     */
    void
    qrFactorization(brick::numeric::Array2D<brick::common::Float64> const& inputArray,
                    brick::numeric::Array2D<brick::common::Float64>& qArray,
                    brick::numeric::Array2D<brick::common::Float64>& rArray);


    /**
     * This function computes the singular value decomposition of a
     * matrix.  After a successful call,
     * matrixMultiply(matrixMultiply(uArray, sigmaArray), vTransposeArray)
     *  == inputArray.
     *
     * @param inputArray This argument is the matrix to be decomposed.
     *
     * @param uArray This argument will be filled in with an
     * orthonormal basis spanning the range of the input matrix, one
     * basis vector per column.  If inputArray has M rows and N
     * columns, and argument isFullRangeRequired (actually
     * isFullRangeRequired is currently ignored, and argument
     * isNullSpaceRequired is used instead) is set to false, then this
     * array will have M rows and min(M, N) columns when the function
     * returns.  If isFullRangeRequired (again, actually
     * isNullSpaceRequired) is set to true, then this array will have
     * M rows and M columns when the function returns.
     *
     * @param sigmaArray This argument will be filled in with the
     * singular values of the matrix, in descending order.  If
     * inputArray has M rows and N columns, then this array will have
     * min(M, N) elements when the function returns.
     *
     * @param vTransposeArray This argument will be filled in with an
     * orthonormal basis spanning the domain of the input matrix, one
     * basis vector per row.  If inputArray has M rows and N columns
     * and argument isNullSpaceRequired is set to false, then this
     * array will have min(M, N) rows and N columns when the function
     * returns.  If isNullSpaceRequired is set to true, then this
     * array will have N rows and N columns when the function returns.
     *
     * @param isNullSpaceRequired This argument specifies how many
     * rows of vTransposeArray will be returned, and temporarily also
     * controls the behavior that would normally be controlled by
     * argument isFullRangeRequired.  Please see the documentation of
     * arguments vTransposeArray and isFullRangeRequired for more
     * information.
     *
     * @param isFullRangeRequired Ostensibly, this argument is
     * specifies how many columns of uArray will be returned, however
     * it is currently ignored, and the value of isNullSpaceRequired
     * is used instead.  Please see the documentation of argument
     * uArray for more information.
     */
    void
    singularValueDecomposition(
      brick::numeric::Array2D<brick::common::Float64> const& inputArray,
      brick::numeric::Array2D<brick::common::Float64>& uArray,
      brick::numeric::Array1D<brick::common::Float64>& sigmaArray,
      brick::numeric::Array2D<brick::common::Float64>& vTransposeArray,
      bool isNullSpaceRequired=false,
      bool isFullRangeRequired=false);


    /**
     * This function computes the singular values a matrix without
     * computing the associated U and V matrices.
     *
     * @param inputArray This argument is the matrix from which to compute
     * the singular values.
     *
     * @return The return value is an Array1D of singular values in
     * descending order.
     */
    brick::numeric::Array1D<brick::common::Float64>
    singularValues(brick::numeric::Array2D<brick::common::Float64> const& inputArray);

  } // namespace linearAlgebra

} // namespace brick

// Include file containing definitions of inline and template
// functions.
#include <brick/linearAlgebra/linearAlgebra_impl.hh>

#endif // #ifndef BRICK_LINEARALGEBRA_LINEARALGEBRA_HH
