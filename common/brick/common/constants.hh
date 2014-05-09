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

namespace brick {

  namespace common {
    
    namespace constants {

      // Basic physical constants.
      const brick::common::Float64 e =
        2.71828182845904523536028747135266249775724709369995;
      const brick::common::Float64 pi =
        3.14159265358979323846264338327950288419716939937510;
      const brick::common::Float64 rootTwo =
        1.41421356237309504880168872420969807856967187537694;

      // Empirical constants.
      const brick::common::Float64 avogadro = 6.02214179E23;

      // Constants that derive from the basic ones.
      const brick::common::Float64 degreesPerRadian = 180.0 / pi;
      const brick::common::Float64 piOverTwo = pi / 2.0;
      const brick::common::Float64 radiansPerDegree = pi / 180.0;
      const brick::common::Float64 rootTwoOverTwo = rootTwo / 2.0;
      const brick::common::Float64 rootOverTwo = rootTwoOverTwo;
      const brick::common::Float64 twoPi = 2.0 * pi;

    } // namespace contants

  } // namespace common

} // namespace brick

#endif /* #ifndef BRICK_COMMON_CONSTANTS_HH */
