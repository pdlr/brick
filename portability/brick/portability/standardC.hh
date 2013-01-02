/**
***************************************************************************
* @file brick/portability/standardC.hh
*
* Header file declaring routines that conform to ANSI C, but are
* missing from various compilers.
*
* (C) Copyright 2004-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying LICENSE file for details.
*
***************************************************************************
**/

#ifndef BRICK_PORTABILITY_STANDARDC_HH
#define BRICK_PORTABILITY_STANDARDC_HH

#include <cstdio>

#ifdef _WIN32

#ifndef HAVE_SNPRINTF
#define HAVE_SNPRINTF 0
#endif

#else /* #ifdef _WIN32 */

#ifndef HAVE_SNPRINTF
#define HAVE_SNPRINTF 1
#endif

// Some compilers declare close() in unistd.h, rather than fcntl.h.
// No big deal, but unistd.h is not portable, so we include it here.
// Now client code can include standardC.hh and not have to include
// unistd.h.
#include <unistd.h>

#endif /* #ifdef _WIN32 ... #else */


namespace brick {

  /**
   ** This namespace contains as much of the platform-specific code as
   ** possible, making the other libraries more readable.
   **/
  namespace portability {
    
#if HAVE_SNPRINTF
  using std::snprintf;
#else

  /** 
   * This function is a crippled replacement for std::snprintf.  It is
   * useful on architectures which do not supply a fully featured
   * version.
   * 
   * @param targetString This argument is the buffer into which the
   * formatted string should be written.
   * 
   * @param targetSize This argument Is the size of the target buffer
   * in bytes.
   * 
   * @param format This argument is the format string in normal printf
   * style.
   * 
   * @return The return value is an integer indicating how many bytes
   * were written.
   */
  int snprintf(char* targetString, size_t targetSize, const char* format,
               ...);
  
#endif

  } // namespace portability
  
} // namespace brick

#endif /* #ifndef BRICK_PORTABILITY_STANDARDC_HH */
