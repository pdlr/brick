/**
***************************************************************************
* @file brick/random/pseudoRandom.hh
*
* Header file declaring PseudoRandom class.
*
* Copyright (C) 1997-2011 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_RANDOM_PSEUDORANDOM_HH
#define BRICK_RANDOM_PSEUDORANDOM_HH

#include <brick/common/types.hh>

namespace brick {

  /**
   ** This namespace contains classes for generating pseudo-random
   ** numbers (and, eventually, random numbers).
   **/
  namespace random {

    /**
     ** The PseudoRandom class generates psuedo-random numbers.
     **/
    class PseudoRandom {
    public:

      /**
       * The default constructor initializes the random number generator
       * with a seed derived from the system clock.  If you create two
       * PseudoRandom instances with a millisecond or so of each other,
       * you may want to check that they have different seeds by calling
       * the getCurrentSeed() member function.
       */
      PseudoRandom();


      /**
       * This constructor sets the seed of the random number generator.
       * This is useful if you need a repeatable sequence of
       * pseudo-random numbers.
       *
       * @param seed This argument is a long long integer.  Note that
       * currently only the low-order 47 bits are used.  The remaining
       * bits are ignored.
       */
      PseudoRandom(brick::common::Int64 seed);


      /**
       * The destructor cleans up any allocated resources.
       */
      ~PseudoRandom() {}


      /**
       * This member function returns a Float64 drawn from a Gaussian
       * distribution with the specified mean and standard deviation.
       *
       * @param mu This argument specifies the mean of the desired
       * distribution.
       *
       * @param sigma This argument specifies the standard deviation of
       * the desired distribution.
       *
       * @return The return value is a sample drawn from a Gaussian
       * distribution.
       */
      brick::common::Float64
      gaussian(brick::common::Float64 mu, brick::common::Float64 sigma);


      /**
       * This member function returns the current state of the random
       * number generator.  Note that the result of this call will
       * change each time a random number is generated.  You can pick up
       * at any given place in the sequence of random numbers by
       * recording the seed at that point and setting it later.
       *
       * @return The return value is a long long int.  Currently only
       * the low-order 47 bits are meaningful.  The remaining bits will
       * be set to zero.
       */
      brick::common::Int64
      getCurrentSeed();


      /**
       * This member function returns a Float64 drawn from a Gaussian
       * distribution with mean equal to 0.0 and and standard deviation
       * equal to 1.0.
       *
       * @return The return value is a sample drawn from a normal
       * distribution.
       */
      brick::common::Float64
      normal();


      /**
       * This member function sets the seed for the random number
       * generator.  This is useful if you need a repeatable sequence of
       * pseudoRandom numbers.
       *
       * @param seed This argument is a long long integer.  Note that
       * currently only the low-order 47 bits are used.  The remaining
       * bits are ignored.
       */
      void
      setCurrentSeed(brick::common::Int64 seed);


      /**
       * This member function returns a Float64, x, drawn from a uniform
       * distribution with bounds such that:
       *
       *   lowerBound <= x < upperBound
       *
       * @param lowerBound This argument specifies the lower bound of
       * the uniform distribution.
       *
       * @param upperBound This argument specifies the upper bound of
       * the uniform distribution.
       *
       * @return The return value is a sample drawn from a uniform
       * distribution.
       */
      brick::common::Float64
      uniform(brick::common::Float64 lowerBound,
              brick::common::Float64 upperBound);


      /**
       * This member function returns a integer, x, drawn from a uniform
       * distribution with bounds such that:
       *
       *   lowerBound <= x < upperBound
       *
       * @param lowerBound This argument specifies the lower bound of
       * the uniform distribution.
       *
       * @param upperBound This argument specifies the upper bound of
       * the uniform distribution.
       *
       * @return The return value is a sample drawn from a uniform
       * distribution.
       */
      int
      uniformInt(int lowerBound, int upperBound);

    private:

      /**
       * This private member function sets the the internal state which
       * will be used to seed the LAPACK function call
       *
       * @param seed0 This argument is a seed integer, and must be in
       * the range [0,4096).
       *
       * @param seed1 This argument is an arbitrary integer, and must be
       * in the range [0,4096).
       *
       * @param seed2 This argument is an arbitrary integer, and must be
       * in the range [0,4096).
       *
       * @param seed3 This argument is an arbitrary integer, must be
       * odd, and must be in the range [0,4096).
       */
      void
      setLapackSeed(brick::common::Int32 seed0, brick::common::Int32 seed1,
                    brick::common::Int32 seed2, brick::common::Int32 seed3);


      /**
       * The four seed values used by the pseudo-random number generator.
       */
      brick::common::Int32 m_seed[4];

    }; // class PseudoRandom

  } // namespace random

} // namespace brick

#endif /* #ifndef BRICK_RANDOM_PSEUDORANDOM_HH */
