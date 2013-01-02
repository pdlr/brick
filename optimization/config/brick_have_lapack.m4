dnl
dnl BRICK_HAVE_LAPACK()
dnl
dnl Call this autoconf macro to check whether the lapack library
dnl development files are installed on your system.

AC_DEFUN([BRICK_HAVE_LAPACK],
    if test "x$LAPACK_LIBS" == x; then
        LAPACK_LIBS="-llapack -lblas"
    fi

    [BRICK_SAVED_CPPFLAGS="${CPPFLAGS}"
     CPPFLAGS="${CPPFLAGS} -I${includedir} -I${prefix}/include"
     AC_MSG_CHECKING([for -llapack])
     BRICK_SAVED_LDFLAGS="${LDFLAGS}"
     BRICK_SAVED_LIBS="${LIBS}"
     LDFLAGS="${LDFLAGS} -L${libdir} -L${prefix}/lib"
     LIBS="${LAPACK_LIBS}"
     AC_TRY_LINK(
         [
           #include <stdint.h>
           extern "C" {void dpotrf_(char*, int32_t*, int64_t*, int32_t*, int32_t*);}
         ],
         [
           char* arg0 = 0;
           int32_t* arg1 = 0;
           int64_t* arg2 = 0;
           int32_t* arg3 = 0;
           int32_t* arg4 = 0;
           dpotrf_(arg0, arg1, arg2, arg3, arg4);
         ],
         [AC_MSG_RESULT([yes])],
         [AC_MSG_ERROR([Unable to build a test program against 
          liblapack.  Please check that the lapack and blas development 
          libraries are installed and accessible.  If they are, and this
          test still fails, try modifying the LDFLAGS environment variable
          so that the compiler can more easily find the missing libraries.  
          Please see ./config.log for details of the failure.])])
     LIBS="${BRICK_SAVED_LIBS}"
     LDFLAGS="$BRICK_SAVED_LDFLAGS"
     CPPFLAGS="$BRICK_SAVED_CPPFLAGS"])
