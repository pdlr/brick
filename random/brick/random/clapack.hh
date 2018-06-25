/**
***************************************************************************
* @file brick/random/clapack.hh
* Klugey header file to declare the LAPACK routines we need for
* brick::random.
*
* Copyright (C) 2001-2004 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/


#ifndef BRICK_RANDOM_CLAPACK_HH
#define BRICK_RANDOM_CLAPACK_HH

#include <brick/common/types.hh>

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * This is a declaration for the LAPACK routine dlarnv(), which
   * computes a vector of random real numbers from a uniform
   * distribution.
   */
  void dlarnv_(brick::common::Int32* IDIST, brick::common::Int32* ISEED,
               brick::common::Int32* N, brick::common::Float64* X);


#ifdef __cplusplus
}
#endif
#endif /* BRICK_RANDOM_CLAPACK_HH */
