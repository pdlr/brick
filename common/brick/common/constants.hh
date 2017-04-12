/**
***************************************************************************
* @file brick/common/constants.hh
*
* Header file declaring useful physical constants.
*
* Copyright (C) 2009-2010 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_COMMON_CONSTANTS_HH
#define BRICK_COMMON_CONSTANTS_HH

#include <brick/common/types.hh>

#if __cplusplus <= 199711L
#define BRICK_CONSTEXPR const
#else
#define BRICK_CONSTEXPR constexpr
#endif

namespace brick {

  namespace common {
    
    namespace constants {

      // Basic physical constants.
      BRICK_CONSTEXPR brick::common::Float64 e =
        2.71828182845904523536028747135266249775724709369995;
      BRICK_CONSTEXPR brick::common::Float64 pi =
        3.14159265358979323846264338327950288419716939937510;
      BRICK_CONSTEXPR brick::common::Float64 rootTwo =
        1.41421356237309504880168872420969807856967187537694;

      // Empirical constants.
      BRICK_CONSTEXPR brick::common::Float64 avogadro = 6.02214179E23;

      // Constants that derive from the basic ones.

      BRICK_CONSTEXPR brick::common::Float64 degreesPerRadian = 180.0 / pi;
      BRICK_CONSTEXPR brick::common::Float64 piOverFour = pi / 4.0;
      BRICK_CONSTEXPR brick::common::Float64 piOverTwo = pi / 2.0;
      BRICK_CONSTEXPR brick::common::Float64 radiansPerDegree = pi / 180.0;
      BRICK_CONSTEXPR brick::common::Float64 rootTwoOverTwo = rootTwo / 2.0;
      BRICK_CONSTEXPR brick::common::Float64 rootOverTwo = rootTwoOverTwo;
      BRICK_CONSTEXPR brick::common::Float64 threePiOverFour = piOverTwo + piOverFour;
      BRICK_CONSTEXPR brick::common::Float64 threePiOverTwo = pi + piOverTwo;
      BRICK_CONSTEXPR brick::common::Float64 twoPi = 2.0 * pi;
      
    } // namespace contants

  } // namespace common

} // namespace brick

#endif /* #ifndef BRICK_COMMON_CONSTANTS_HH */
