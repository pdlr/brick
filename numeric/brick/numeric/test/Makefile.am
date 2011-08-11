#
# Makefile.am
#   Process this file with automake to produce Makefile.in.
#
# $Revision: 1277 $
# $Date: 2010-01-15 00:35:26 -0500 (Fri, 15 Jan 2010) $
#


# Set VPATH so we can build in a different directory.

srcdir			= @srcdir@
VPATH			= @srcdir@


# Specify which tests should be run.

TESTS			= \
	array1DTest \
	array2DTest \
	array3DTest \
	index3DTest \
	solveQuadraticTest \
	transform2DTest \
	transform3DTest \
	utilitiesTest


# Specify which tests should be build from source.

check_PROGRAMS          = \
	array1DTest \
	array2DTest \
	array3DTest \
	index3DTest \
	solveQuadraticTest \
	transform2DTest \
	transform3DTest \
	utilitiesTest


# For each test which should be built from source, specify which 
# source files.

array1DTest_SOURCES		= array1DTest.cc arrayTestCommon.hh
array2DTest_SOURCES		= array2DTest.cc arrayTestCommon.hh
array3DTest_SOURCES		= array3DTest.cc arrayTestCommon.hh
index3DTest_SOURCES		= index3DTest.cc
solveQuadraticTest_SOURCES	= solveQuadraticTest.cc
transform2DTest_SOURCES		= transform2DTest.cc arrayTestCommon.hh
transform3DTest_SOURCES		= transform3DTest.cc arrayTestCommon.hh
utilitiesTest_SOURCES		= utilitiesTest.cc


# Set the include path.

AM_CPPFLAGS		= -I$(top_srcdir)


# Explicitly set library locations to avoid linking with previously
# installed versions of the libraries.

AM_LDFLAGS		= -L$(top_builddir)/brick/numeric


# Specify which libraries are needed to build these tests.

LDADD			= -lbrickNumeric \
			  -lbrickPortability \
			  -lbrickCommon \
			  -lbrickTest \
			  -lbrickTestAutoMain


# List the files that need to be cleaned.

MAINTAINERCLEANFILES	= Makefile.in
