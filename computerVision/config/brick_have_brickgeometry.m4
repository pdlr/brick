dnl
dnl BRICK_HAVE_BRICKGEOMETRY()
dnl
dnl Call this autoconf macro to check whether the brickGeometry library
dnl development files are installed on your system.

AC_DEFUN([BRICK_HAVE_BRICKGEOMETRY],
    [BRICK_SAVED_CPPFLAGS="${CPPFLAGS}"
     CPPFLAGS="${CPPFLAGS} -I${includedir} -I${prefix}/include"
     AC_CHECK_HEADER(brick/geometry/ray3D.h, 
         [],
         [AC_MSG_ERROR([Unable to find brick/geometry/ray3D.h.  
              Please check that the dlrGeometry utility library is 
              installed.  If it is, and this test still fails, try 
              modifying the CPPFLAGS environment variable so that the 
              compiler can more easily find this include file.  Please see 
              ./config.log for details of the failure.])])
     AC_MSG_CHECKING([for -lbrickGeometry])
     BRICK_SAVED_LDFLAGS="${LDFLAGS}"
     BRICK_SAVED_LIBS="${LIBS}"
     LDFLAGS="${LDFLAGS} -L${libdir} -L${prefix}/lib"
     LIBS="-lbrickGeometry -lbrickLinearAlgebra -lbrickNumeric -lbrickPortability -lbrickCommon ${LAPACK_LIBS} ${BLAS_LIBS} ${BRICK_SAVED_LIBS} ${FLIBS}"
     AC_TRY_LINK([#include <brick/geometry/ray3D.h>],
         [dlr::geometry::Ray3D<double> dummy;],
         [AC_MSG_RESULT([yes])],
         [AC_MSG_ERROR([Unable to build a test program against 
          libdlrGeometry.  Please check that the dlrGeometry 
          utility library is installed.  If it is, and this test still fails,
          try modifying the LDFLAGS environment variable so that the compiler 
          can more easily find this library.  Please see ./config.log for 
          details of the failure.])])
     LIBS="${BRICK_SAVED_LIBS}"
     LDFLAGS="$BRICK_SAVED_LDFLAGS"
     CPPFLAGS="$BRICK_SAVED_CPPFLAGS"])
