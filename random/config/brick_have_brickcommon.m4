dnl
dnl BRICK_HAVE_BRICKCOMMON()
dnl
dnl Call this autoconf macro to check whether the brick::common library
dnl development files are installed on your system.

AC_DEFUN([BRICK_HAVE_BRICKCOMMON],
    [BRICK_SAVED_CPPFLAGS="${CPPFLAGS}"
     CPPFLAGS="${CPPFLAGS} -I${includedir} -I${prefix}/include"
     AC_CHECK_HEADER(brick/common/types.hh, 
         [],
         [AC_MSG_ERROR([Unable to find brick/common/types.hh.  Please check 
              that the brickCommon utility library is installed.  If it is, 
              and this test still fails, try modifying the CPPFLAGS environment
              variable so that the compiler can more easily find the missing
              include file.
              Please see ./config.log for details of the failure.])])
     AC_MSG_CHECKING([for -lbrickCommon])
     BRICK_SAVED_LDFLAGS="${LDFLAGS}"
     BRICK_SAVED_LIBS="${LIBS}"
     LDFLAGS="${LDFLAGS} -L${libdir} -L${prefix}/lib"
     LIBS="-lbrickCommon"
     AC_TRY_LINK([#include <brick/common/types.hh>],
         [brick::common::Float64 f = 3.0; f = 2.0*f;],
         [AC_MSG_RESULT([yes])],
         [AC_MSG_ERROR([Unable to build a test program against 
          libbrickCommon.  Please check that the brickCommon utility 
          library is installed.  If it is, and this test still fails, try 
          modifying the LDFLAGS environment variable so that the compiler 
          can more easily find the missing library.  
          Please see ./config.log for details of the failure.])])
     LIBS="${BRICK_SAVED_LIBS}"
     LDFLAGS="$BRICK_SAVED_LDFLAGS"
     CPPFLAGS="$BRICK_SAVED_CPPFLAGS"])
