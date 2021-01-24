/**
***************************************************************************
* @file brick/linearAlgebra/clapack.h
* Klugey header file to declare the LAPACK routines we need.
*
* Copyright (C) 2001-2011 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/


#ifndef BRICK_LINEARALGEBRA_CLAPACK_HH
#define BRICK_LINEARALGEBRA_CLAPACK_HH

#include <brick/common/types.hh>

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * This is a declaration for the LAPACK routine dgeev(), which
   * computes the eigenvalues and eigenvectors of a general square
   * matrix.
   */
  void dgeev_(char* JOBVL, char* JOBVR, brick::common::Int32* N,
              brick::common::Float64* A, brick::common::Int32* LDA,
              brick::common::Float64* WR, brick::common::Float64* WI,
              brick::common::Float64* VL, brick::common::Int32* LDVL,
              brick::common::Float64* VR, brick::common::Int32* LDVR,
              brick::common::Float64* WORK, brick::common::Int32* LWORK,
              brick::common::Int32* INFO);


  /**
   * This is a declaration for the LAPACK routine dgels(), which
   * computes the solution of a general system of linear equations.
   */
  void dgels_(char* TRANS, brick::common::Int32* M, brick::common::Int32* N,
              brick::common::Int32* NRHS,
              brick::common::Float64* A, brick::common::Int32* LDA,
              brick::common::Float64* B, brick::common::Int32* LDB,
              brick::common::Float64* WORK, brick::common::Int32* LWORK,
              brick::common::Int32* INFO);


  /**
   * This is a declaration for the LAPACK routine dgesdd(), which
   * computes the singular value decomposition of a matrix using a
   * fast divide & conquer algorithm.  There's another SVD routine in
   * LAPACK that's slower, but uses less space.
   */
  void dgesdd_(char* JOBZ, brick::common::Int32* M, brick::common::Int32* N,
               brick::common::Float64* A, brick::common::Int32* LDA,
               brick::common::Float64* S,
               brick::common::Float64* U, brick::common::Int32* LDU,
               brick::common::Float64* VT, brick::common::Int32* LDVT,
               brick::common::Float64* WORK, brick::common::Int32* LWORK,
               brick::common::Int32* IWORK, brick::common::Int32* INFO);


  /**
   * This is a declaration for the LAPACK routine dgesv(), which
   * computes the solution to a systems of linear equations.
   */
  void dgesv_(brick::common::Int32* N, brick::common::Int32* NRHS,
              brick::common::Float64 *A, brick::common::Int32* LDA,
              brick::common::Int32* IPIV,
              brick::common::Float64* B, brick::common::Int32* LDB,
              brick::common::Int32* INFO);


  /**
   * This is a declaration for the LAPACK routine dgeqrf(), which
   * computes the QR factorization of a general MxN matrix.
   */
  void dgeqrf_(brick::common::Int32* M, brick::common::Int32* N,
               brick::common::Float64* A, brick::common::Int32* LDA,
               brick::common::Float64* TAU, brick::common::Float64* WORK,
               brick::common::Int32* LWORK, brick::common::Int32* INFO);


  /**
   * This is a declaration for the LAPACK routine dgetrf(), which
   * computes LU decomposition of a general MxN matrix.
   */
  void dgetrf_(brick::common::Int32* M, brick::common::Int32* N,
               brick::common::Float64* A, brick::common::Int32* LDA,
               brick::common::Int32* IPIV, brick::common::Int32* INFO);


  /**
   * This is a declaration for the LAPACK routine dlarnv(), which
   * computes a vector of random real numbers from a uniform
   * distribution.
   */
  void dlarnv_(brick::common::Int32* IDIST, brick::common::Int32* ISEED,
               brick::common::Int32* N, brick::common::Float64* X);


  /**
   * This is a declaration for the LAPACK routine dgtsv(), which
   * computes the solution of a general tridiagonal system of linear
   * equations.
   */
  void dgtsv_(brick::common::Int32* N, brick::common::Int32* NRHS,
              brick::common::Float64* DL, brick::common::Float64* D,
              brick::common::Float64* DU,
              brick::common::Float64* B, brick::common::Int32* LDB,
              brick::common::Int32* INFO);


  /**
   * This is a declaration for the LAPACK routine dpotrf(), which
   * computes the Cholesky factorization of a symmetric positive
   * definite matrix.
   */
  void dpotrf_(char* UPLO, brick::common::Int32* N, brick::common::Float64* A,
               brick::common::Int32* LDA, brick::common::Int32* INFO);


  /**
   * This is a declaration for the LAPACK routine dsyev(), which
   * computes the eigenvalues and eigenvectors of a real symmetric
   * matrix.
   */
  void dsyev_(char* JOBZ,  char* UPLO, brick::common::Int32* N,
              brick::common::Float64* A, brick::common::Int32* LDA,
              brick::common::Float64* W, brick::common::Float64* WORK,
              brick::common::Int32* LWORK, brick::common::Int32* INFO);


  /**
   * This is a declaration for the LAPACK routine dgels(), which
   * computes the solution of a general system of linear equations.
   */
  void sgels_(char* TRANS, brick::common::Int32* M, brick::common::Int32* N,
              brick::common::Int32* NRHS,
              brick::common::Float32* A, brick::common::Int32* LDA,
              brick::common::Float32* B, brick::common::Int32* LDB,
              brick::common::Float32* WORK, brick::common::Int32* LWORK,
              brick::common::Int32* INFO);


  /**
   * This is a declaration for the LAPACK routine sgesv(), which
   * computes the solution to a systems of linear equations.
   */
  void sgesv_(brick::common::Int32* N, brick::common::Int32* NRHS,
              brick::common::Float32 *A, brick::common::Int32* LDA,
              brick::common::Int32* IPIV,
              brick::common::Float32* B, brick::common::Int32* LDB,
              brick::common::Int32* INFO);

#ifdef __cplusplus
}
#endif
#endif /* BRICK_LINEARALGEBRA_CLAPACK_HH */
