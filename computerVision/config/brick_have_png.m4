dnl
dnl BRICK_HAVE_PNG()
dnl
dnl Call this autoconf macro to check whether the brick::common library
dnl development files are installed on your system.

AC_DEFUN([BRICK_HAVE_PNG],
    if test "x$PNG_LIBS" == x; then
        PNG_LIBS="-lpng"
    fi

    [BRICK_SAVED_CPPFLAGS="${CPPFLAGS}"
     CPPFLAGS="${CPPFLAGS} -I${includedir} -I${prefix}/include"
     AC_CHECK_HEADER(png.h, 
         [],
         [AC_MSG_ERROR([Unable to find png.h.  Please check 
              that the libpng12 development files are installed.  If they are, 
              and this test still fails, try modifying the CPPFLAGS environment
              variable so that the compiler can more easily find the missing
              include file.
              Please see ./config.log for details of the failure.])])
     AC_MSG_CHECKING([for -lpng])
     BRICK_SAVED_LDFLAGS="${LDFLAGS}"
     BRICK_SAVED_LIBS="${LIBS}"
     LDFLAGS="${LDFLAGS} -L${libdir} -L${prefix}/lib"
     LIBS="${PNG_LIBS}"
     AC_TRY_LINK([#include <png.h>],
         [png_infop* infoP = 0; infoP += 1;],
         [AC_MSG_RESULT([yes])],
         [AC_MSG_ERROR([Unable to build a test program against 
          libpng.  Please check that the libpng12 development files
          are installed.  If they are, and this test still fails, try 
          modifying the LDFLAGS environment variable so that the compiler 
          can more easily find the missing library.  
          Please see ./config.log for details of the failure.])])
     LIBS="${BRICK_SAVED_LIBS}"
     LDFLAGS="$BRICK_SAVED_LDFLAGS"
     CPPFLAGS="$BRICK_SAVED_CPPFLAGS"])
