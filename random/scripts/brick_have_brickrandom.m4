dnl
dnl BRICK_HAVE_BRICKRANDOM()
dnl
dnl Call this autoconf macro to check whether the brick::random library
dnl development files are installed.

AC_DEFUN([BRICK_HAVE_BRICKRANDOM],
    [BRICK_SAVED_CPPFLAGS="${CPPFLAGS}"
     CPPFLAGS="${CPPFLAGS} -I${includedir} -I${prefix}/include"
     AC_CHECK_HEADER(brick/random/pseudoRandom.h, 
         [],
         [AC_MSG_ERROR([Unable to find brick/random/pseudoRandom.h.  Please
              check that the brick::random utility library is installed.
              If it is, and this test still fails, try modifying the CPPFLAGS
              environment variable so that the compiler can more easily find
              this include file.])])
     AC_MSG_CHECKING([for -lbrickRandom])
     BRICK_SAVED_LDFLAGS="${LDFLAGS}"
     BRICK_SAVED_LIBS="${LIBS}"
     LDFLAGS="${LDFLAGS} -L${libdir} -L${prefix}/lib"
     LIBS="-lbrickRandom -lbrickPortability -lbrickCommon"
     AC_TRY_LINK([#include <brick/random/pseudoRandom.h>],
         [brick::random::PseudoRandom dummy();],
         [AC_MSG_RESULT([yes])],
         [AC_MSG_ERROR([Unable to build a test program against 
          libbrickRandom.  Please check that the brick::random utility 
          library is installed.  If it is, and this test still fails, try 
          modifying the LDFLAGS environment variable so that the compiler 
          can more easily find this library.  Please also check that any 
          dependencies are installed, as described in ./README.TXT.  See 
	  ./config.log for details of the failure.])])
     LIBS="${BRICK_SAVED_LIBS}"
     LDFLAGS="$BRICK_SAVED_LDFLAGS"
     CPPFLAGS="$BRICK_SAVED_CPPFLAGS"])
