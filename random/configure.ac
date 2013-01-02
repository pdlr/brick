dnl
dnl configure.in
dnl   Process this file with autoconf to produce a configure script.
dnl
dnl $Revision: 940 $
dnl $Date: 2007-07-09 16:23:56 -0400 (Mon, 09 Jul 2007) $
dnl


dnl Initialization.
AC_INIT(brickrandom, 1.0)
AC_CONFIG_SRCDIR(brick/random/Makefile.am)
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS(config.h)
# AM_CONFIG_HEADER(config.h:config-h.in)
AM_INIT_AUTOMAKE(1.7)
AM_MAINTAINER_MODE


dnl Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_LIBTOOL


dnl Use C++ compiler to build programs when running ./configure.
AC_LANG([C++])


dnl Checks for libraries.
BRICK_HAVE_LAPACK
BRICK_HAVE_BRICKCOMMON
BRICK_HAVE_BRICKPORTABILITY

dnl Checks for header files.
AC_HEADER_STDC


dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST


dnl Other stuff.
AC_SUBST(ac_aux_dir)


dnl Output.
AC_OUTPUT(Makefile \
          doc/Makefile \
          brick/random/Makefile)

echo \
"------------------------------------------------------------------------
Configuration:

  Source code location:       ${srcdir}
  C compiler:                 ${CC}
  C compiler flags:           ${CFLAGS}
  C++ compiler:               ${CXX}
  C++ compiler flags:         ${CXXFLAGS}
  Linker flags:               ${LDFLAGS}
  Host System Type:           ${host}
  Install path:               ${prefix}

  See config.h for further configuration information.
------------------------------------------------------------------------"
