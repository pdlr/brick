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
      const brick::common::Float64 e = 2.7182818284590451;
      const brick::common::Float64 pi = 3.141592653589793238;

      // Empirical constants.
      const brick::common::Float64 avogadro = 6.02214179E23;

      // Constants that derive from the basic ones.
      const brick::common::Float64 degreesPerRadian = 180.0 / pi;
      const brick::common::Float64 piOverTwo = 1.5707963267948966;
      const brick::common::Float64 radiansPerDegree = pi / 180.0;
      const brick::common::Float64 rootTwo = 1.4142135623730951;
      const brick::common::Float64 rootOverTwo = 0.70710678118654757;
      const brick::common::Float64 twoPi = 2.0 * pi;

    } // namespace contants

  } // namespace common

} // namespace brick

#endif /* #ifndef BRICK_COMMON_CONSTANTS_HH */
