dnl
dnl BRICK_HAVE_BRICKNUMERIC()
dnl
dnl Call this autoconf macro to check whether the brick::numeric library
dnl development files are installed on your system.

AC_DEFUN([BRICK_HAVE_BRICKNUMERIC],
    [BRICK_SAVED_CPPFLAGS="${CPPFLAGS}"
     CPPFLAGS="${CPPFLAGS} -I${includedir} -I${prefix}/include"
     AC_CHECK_HEADER(brick/numeric/index2D.hh, 
         [],
         [AC_MSG_ERROR([Unable to find brick/numeric/index2D.hh.  Please
              check that the brick::numeric utility library is installed.
              If it is, and this test still fails, try modifying the CPPFLAGS
              environment variable so that the compiler can more easily find
              the missing include file.])])
     AC_MSG_CHECKING([for -lbrickNumeric])
     BRICK_SAVED_LDFLAGS="${LDFLAGS}"
     BRICK_SAVED_LIBS="${LIBS}"
     LDFLAGS="${LDFLAGS} -L${libdir} -L${prefix}/lib"
     LIBS="-lbrickNumeric -lbrickPortability -lbrickCommon"
     AC_TRY_LINK([#include <brick/numeric/index2D.hh>],
         [brick::numeric::Index2D dummy; dummy.setValue(0, 0);],
         [AC_MSG_RESULT([yes])],
         [AC_MSG_ERROR([Unable to build a test program against libbrickNumeric.
          Please check that the brick::numeric utility library is installed. 
          If it is, and this test still fails, try modifying the LDFLAGS 
          environment variable so that the compiler can more easily find 
          the missing library.  Please also check that any dependencies are 
          installed, as described in ./README.TXT.  Please see ./config.log for 
          details of the failure.])])
     LIBS="${BRICK_SAVED_LIBS}"
     LDFLAGS="$BRICK_SAVED_LDFLAGS"
     CPPFLAGS="$BRICK_SAVED_CPPFLAGS"])
