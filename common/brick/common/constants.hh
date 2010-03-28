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

namespace brick {

  namespace common {
    
    namespace constants {

      // Basic physical constants.
      const double e = 2.7182818284590451;
      const double pi = 3.141592653589793238;

      // Empirical constants.
      const double avogadro = 6.02214179E23;

      // Constants that derive from the basic ones.
      const double degreesPerRadian = 180.0 / pi;
      const double piOverTwo = 1.5707963267948966;
      const double radiansPerDegree = pi / 180.0;
      const double rootTwo = 1.4142135623730951;
      const double rootOverTwo = 0.70710678118654757;
      
    } // nameapace contants

  } // namespace common

} // namespace brick

#endif /* #ifndef BRICK_COMMON_CONSTANTS_HH */
