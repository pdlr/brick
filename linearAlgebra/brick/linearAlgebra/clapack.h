/**
***************************************************************************
* @file dlrLinearAlgebra/clapack.h
* Klugey header file to declare the LAPACK routines we need.
*
* Copyright (C) 2001-2004 David LaRose, dlr@cs.cmu.edu
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
  void dgeev_(char* JOBVL, char* JOBVR, brick::Int32* N,
              brick::Float64* A, dlr::Int32* LDA,
              brick::Float64* WR, dlr::Float64* WI, 
              brick::Float64* VL, dlr::Int32* LDVL,
              brick::Float64* VR, dlr::Int32* LDVR,
              brick::Float64* WORK, dlr::Int32* LWORK,
              brick::Int32* INFO);

  
  /**
   * This is a declaration for the LAPACK routine dgels(), which
   * computes the solution of a general system of linear equations.
   */
  void dgels_(char* TRANS, brick::Int32* M, dlr::Int32* N, dlr::Int32* NRHS,
              brick::Float64* A, dlr::Int32* LDA,
              brick::Float64* B, dlr::Int32* LDB,
              brick::Float64* WORK, dlr::Int32* LWORK, dlr::Int32* INFO);

  
  /**
   * This is a declaration for the LAPACK routine dgesdd(), which
   * computes the singular value decomposition of a matrix using a
   * fast divide & conquer algorithm.  There's another SVD routine in
   * LAPACK that's slower, but uses less space.
   */
  void dgesdd_(char* JOBZ, brick::Int32* M, dlr::Int32* N,
               brick::Float64* A, dlr::Int32* LDA, dlr::Float64* S,
               brick::Float64* U, dlr::Int32* LDU,
               brick::Float64* VT, dlr::Int32* LDVT,
               brick::Float64* WORK, dlr::Int32* LWORK,
               brick::Int32* IWORK, dlr::Int32* INFO);

  
  /**
   * This is a declaration for the LAPACK routine dgesv(), which
   * computes the solution to a systems of linear equations.
   */
  void dgesv_(brick::Int32* N, dlr::Int32* NRHS,
              brick::Float64 *A, dlr::Int32* LDA,
              brick::Int32* IPIV,
              brick::Float64* B, dlr::Int32* LDB,
              brick::Int32* INFO);

  
  /**
   * This is a declaration for the LAPACK routine dgeqrf(), which
   * computes the QR factorization of a general MxN matrix.
   */
  void dgeqrf_(brick::Int32* M, dlr::Int32* N,
               brick::Float64* A, dlr::Int32* LDA,
               brick::Float64* TAU, dlr::Float64* WORK,
               brick::Int32* LWORK, dlr::Int32* INFO);


  /**
   * This is a declaration for the LAPACK routine dgetrf(), which
   * computes LU decomposition of a general MxN matrix.
   */
  void dgetrf_(brick::Int32* M, dlr::Int32* N,
               brick::Float64* A, dlr::Int32* LDA,
               brick::Int32* IPIV, dlr::Int32* INFO);


  /**
   * This is a declaration for the LAPACK routine dlarnv(), which
   * computes a vector of random real numbers from a uniform
   * distribution.
   */
  void dlarnv_(brick::Int32* IDIST, dlr::Int32* ISEED, dlr::Int32* N,
               brick::Float64* X);


  /**
   * This is a declaration for the LAPACK routine dgtsv(), which
   * computes the solution of a general tridiagonal system of linear
   * equations.
   */
  void dgtsv_(brick::Int32* N, dlr::Int32* NRHS,
              brick::Float64* DL, dlr::Float64* D, dlr::Float64* DU,
              brick::Float64* B, dlr::Int32* LDB, dlr::Int32* INFO);
  

  /**
   * This is a declaration for the LAPACK routine dpotrf(), which
   * computes the Cholesky factorization of a symmetric positive
   * definite matrix.
   */
  void dpotrf_(char* UPLO, brick::Int32* N, dlr::Float64* A,
               brick::Int32* LDA, dlr::Int32* INFO);

  
  /**
   * This is a declaration for the LAPACK routine dsyev(), which
   * computes the eigenvalues and eigenvectors of a real symmetric
   * matrix.
   */
  void dsyev_(char* JOBZ,  char* UPLO, brick::Int32* N,
              brick::Float64* A, dlr::Int32* LDA, dlr::Float64* W,
              brick::Float64* WORK, dlr::Int32* LWORK, dlr::Int32* INFO);

#ifdef __cplusplus
}
#endif
#endif /* BRICK_LINEARALGEBRA_CLAPACK_HH */
