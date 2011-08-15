dnl
dnl DLR_HAVE_DLRLINEARALGEBRA()
dnl
dnl Call this autoconf macro to check whether the dlrCommon library
dnl development files are installed on your system.

AC_DEFUN([BRICK_HAVE_BRICK_LINEARALGEBRA],
    [BRICK_SAVED_CPPFLAGS="${CPPFLAGS}"
     CPPFLAGS="${CPPFLAGS} -I${includedir} -I${prefix}/include"
     AC_CHECK_HEADER(brick/linearAlgebra/linearAlgebra.hh, 
         [],
         [AC_MSG_ERROR([Unable to find brick/linearAlgebra/linearAlgebra.h.  
              Please check that the brickLinearAlgebra utility library is 
              installed.  If it is, and this test still fails, try 
              modifying the CPPFLAGS environment variable so that the 
              compiler can more easily find this include file.  Please see 
              ./config.log for details of the failure.])])
     AC_MSG_CHECKING([for -lbrickLinearAlgebra])
     BRICK_SAVED_LDFLAGS="${LDFLAGS}"
     BRICK_SAVED_LIBS="${LIBS}"
     LDFLAGS="${LDFLAGS} -L${libdir} -L${prefix}/lib"
     LIBS="-lbrickLinearAlgebra -lbrickNumeric -lbrickPortability -lbrickCommon ${LAPACK_LIBS} ${BLAS_LIBS} ${BRICK_SAVED_LIBS} ${FLIBS}"
     AC_TRY_LINK([#include <brick/linearAlgebra/linearAlgebra.hh>],
         [brick::linearAlgebra::determinant(brick::numeric::Array2D<brick::common::Float64>());],
         [AC_MSG_RESULT([yes])],
         [AC_MSG_ERROR([Unable to build a test program against 
          brickLinearAlgebra.  Please check that the brickLinearAlgebra utility 
          library is installed.  If it is, and this test still fails, try 
          modifying the LDFLAGS environment variable so that the compiler 
          can more easily find this library.  Please see ./config.log for 
          details of the failure.])])
     LIBS="${BRICK_SAVED_LIBS}"
     LDFLAGS="${BRICK_SAVED_LDFLAGS}"
     CPPFLAGS="${BRICK_SAVED_CPPFLAGS}"])
