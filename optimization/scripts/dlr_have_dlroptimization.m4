dnl
dnl DLR_HAVE_DLROPTIMIZATION()
dnl
dnl Call this autoconf macro to check whether the dlrCommon library
dnl development files are installed on your system.

AC_DEFUN([DLR_HAVE_DLROPTIMIZATION],
    [DLR_SAVED_CPPFLAGS="${CPPFLAGS}"
     CPPFLAGS="${CPPFLAGS} -I${includedir} -I${prefix}/include"
     AC_CHECK_HEADER(dlrOptimization/optimizerBFGS.h, 
         [],
         [AC_MSG_ERROR([Unable to find dlrOptimization/optimizerBFGS.h.  
              Please check that the dlrOptimization utility library is 
              installed.  If it is, and this test still fails, try 
              modifying the CPPFLAGS environment variable so that the 
              compiler can more easily find this include file.  Please see 
              ./config.log for details of the failure.])])])
