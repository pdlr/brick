dnl
dnl BRICK_HAVE_BRICKUTILITIES()
dnl
dnl Call this autoconf macro to check whether the brick::utilities library
dnl development files are installed on your system.

AC_DEFUN([BRICK_HAVE_BRICKUTILITIES],
    [BRICK_SAVED_CPPFLAGS="${CPPFLAGS}"
     CPPFLAGS="${CPPFLAGS} -I${includedir} -I${prefix}/include"
     AC_CHECK_HEADER(brick/utilities/timeUtilities.hh, 
         [],
         [AC_MSG_ERROR([Unable to find brick/utilities/timeUtilities.hh. 
              Please check that the brick::utilities utility library is 
              installed.  If it is, and this test still fails, try modifying 
              the CPPFLAGS environment variable so that the compiler can 
              more easily find the missing include file.])])
     AC_MSG_CHECKING([for -lbrickUtilities])
     BRICK_SAVED_LDFLAGS="${LDFLAGS}"
     BRICK_SAVED_LIBS="${LIBS}"
     LDFLAGS="${LDFLAGS} -L${libdir} -L${prefix}/lib"
     LIBS="-lbrickUtilities -lbrickPortability -lbrickCommon"
     AC_TRY_LINK([#include <brick/utilities/timeUtilities.hh>],
         [brick::utilities::getCurrentTime();],
         [AC_MSG_RESULT([yes])],
         [AC_MSG_ERROR([Unable to build a test program against 
          libbrickUtilities.  Please check that the brickUtilities 
          library is installed.  If it is, and this test still fails, try 
          modifying the LDFLAGS environment variable so that the compiler 
          can more easily find the missing library.  Please also check that any 
          dependencies are installed, as described in ./README.TXT.  Please see
          ./config.log for details of the failure.])])
     LIBS="${BRICK_SAVED_LIBS}"
     LDFLAGS="$BRICK_SAVED_LDFLAGS"
     CPPFLAGS="$BRICK_SAVED_CPPFLAGS"])
