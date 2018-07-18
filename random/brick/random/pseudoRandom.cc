/**
***************************************************************************
* @file brick/random/pseudoRandom.cc
*
* Source file defining PseudoRandom class.
*
* Copyright (C) 1997-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/common/exception.hh>
#include <brick/portability/timeUtilities.hh>
#include <brick/random/clapack.hh>
#include <brick/random/pseudoRandom.hh>

namespace brick {

  namespace random {

    // The default constructor initializes the random number generator
    // with a seed derived from the system clock.
    PseudoRandom::
    PseudoRandom()
    {
      // Note that time resolution is not critical here.  The only
      // impact of low resolution is that two PseudorRandom instances
      // created at nearly the same time may wind up on the same
      // sequence.  If you're worried about that, you can use the
      // getCurrentSeed() member function to make sure it hasn't
      // happened.
      brick::common::Float64 currentTime = brick::portability::getCurrentTime();
      brick::common::Int32 seconds =
        static_cast<brick::common::Int32>(currentTime);
      brick::common::Int32 uSeconds =
        static_cast<brick::common::Int32>(1000000.0 * (currentTime - seconds));

      // All seeds must be in the range [0, 4096).
      brick::common::Int32 seed0 =
        static_cast<brick::common::Int32>(seconds & 0x00000fff);
      brick::common::Int32 seed1 =
        static_cast<brick::common::Int32>((seconds & 0x00fff000) >> 12);
      brick::common::Int32 seed2 =
        static_cast<brick::common::Int32>((uSeconds & 0x00000fff));
      brick::common::Int32 seed3 =
        static_cast<brick::common::Int32>((uSeconds & 0x00fff000) >> 12);

      // seed3 must be odd.
      if(seed3 % 2 == 0) {
        seed3 -= 1;
      }

      // All seeds must be >= zero.
      if(seed3 < 0) {
        seed3 = -seed3;
      }

      // try {
      this->setLapackSeed(seed0, seed1, seed2, seed3);
      // } catch(const ValueException&) {
      // This won't happen because we've chosen the seed values wisely.
      // }
    }


    // This constructor sets the seed of the random number generator.
    PseudoRandom::
    PseudoRandom(brick::common::Int64 seed)
    {
      this->setCurrentSeed(seed);
    }


    // This member function returns a Float64 drawn from a Gaussian
    // distribution with the specified mean and standard deviation.
    brick::common::Float64
    PseudoRandom::
    gaussian(brick::common::Float64 mu, brick::common::Float64 sigma)
    {
      brick::common::Float64 returnValue = this->normal();
      return (returnValue * sigma) + mu;
    }


    // This member function returns the current state of the random
    // number generator.
    brick::common::Int64
    PseudoRandom::
    getCurrentSeed()
    {
      brick::common::Int64 seed = 0;
      seed += static_cast<brick::common::Int64>(m_seed[0]) << 35;
      seed += static_cast<brick::common::Int64>(m_seed[1]) << 23;
      seed += static_cast<brick::common::Int64>(m_seed[2]) << 11;
      seed += static_cast<brick::common::Int64>((m_seed[3] & 0x0ffe) >> 1);
      return seed;
    }


    // This member function returns a Float64 drawn from a Gaussian
    // distribution with equal to 0.0 and and standard deviation equal to
    // 1.0.
    brick::common::Float64
    PseudoRandom::
    normal()
    {
      // Tells lapack we want a normal distribution.
      brick::common::Int32 idist = 3;
      // The nuber we'll return.
      brick::common::Float64 returnValue;
      // Only one value needed.
      brick::common::Int32 size = 1;

      // Call the lapack routine
      dlarnv_(&idist, m_seed, &size, &returnValue);

      return returnValue;
    }


    //This member function sets the seed for the random number
    // generator.
    void
    PseudoRandom::
    setCurrentSeed(brick::common::Int64 seed)
    {
      // All seeds must be in the range [0,4096).
      brick::common::Int32 seed0 =
        static_cast<brick::common::Int32>((seed & 0x00007ff800000000LL) >> 35);
      brick::common::Int32 seed1 =
        static_cast<brick::common::Int32>((seed & 0x00000007ff800000LL) >> 23);
      brick::common::Int32 seed2 =
        static_cast<brick::common::Int32>((seed & 0x00000000007ff800LL) >> 11);

      // Seed3 must be odd,
      brick::common::Int32 seed3 =
        static_cast<brick::common::Int32>(((seed & 0x00000000000007ffLL) << 1)
                                          + 1);

      // try {
      this->setLapackSeed(seed0, seed1, seed2, seed3);
      // } catch(const brick::common::ValueException&) {
      // This won't happen because we've chosen the seed values wisely.
      // }
    }


    // This member function returns a Float64, x, drawn from a uniform
    // distribution.
    brick::common::Float64
    PseudoRandom::
    uniform(brick::common::Float64 lowerBound,
            brick::common::Float64 upperBound)
    {
      brick::common::Int32 idist = 1;   // Tells lapack we want uniform [0,1)
      brick::common::Float64 returnValue;   // The number we'll return
      brick::common::Int32 size = 1;    // Only one value needed.

      // Call the lapack routine
      dlarnv_(&idist, m_seed, &size, &returnValue);

      brick::common::Float64 range = upperBound - lowerBound;
      returnValue = returnValue * range + lowerBound;
      return returnValue;
    }


    // This member function returns a integer, x, drawn from a uniform
    // distribution.
    int
    PseudoRandom::
    uniformInt(int lowerBound, int upperBound)
    {
      return static_cast<int>(
        this->uniform(static_cast<brick::common::Float64>(lowerBound),
                      static_cast<brick::common::Float64>(upperBound)));
    }

    void
    PseudoRandom::
    setLapackSeed(brick::common::Int32 seed0, brick::common::Int32 seed1,
                  brick::common::Int32 seed2, brick::common::Int32 seed3)
    {
      if(seed3 % 2 == 0) {
        BRICK_THROW(brick::common::ValueException,
                    "PseudoRandom::setLapackSeed(Int32, Int32, Int32, Int32)",
                    "Seed 3 must be odd.");
      }
      if(seed0 >= 4096 || seed1 >= 4096 || seed2 >= 4096 || seed3 >= 4096
         || seed0 < 0 || seed1 < 0 || seed2 < 0 || seed3 < 0 ) {
        BRICK_THROW(brick::common::ValueException,
                    "PseudoRandom::PseudoRandom(Int32, Int32, Int32, Int32)",
                    "All seeds must be >= 0 and < 4096.");
      }
      m_seed[0] = seed0;
      m_seed[1] = seed1;
      m_seed[2] = seed2;
      m_seed[3] = seed3;
    }


  } // namespace random

} // namespace brick
