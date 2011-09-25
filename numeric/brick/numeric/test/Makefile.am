#
# Makefile.am
#   Process this file with automake to produce Makefile.in.
#

# Set VPATH so we can build in a different directory.

srcdir			= @srcdir@
VPATH			= @srcdir@


# Specify which tests should be run.

TESTS			= \
	array1DTest \
	array2DTest \
	array3DTest \
	derivativeRiddersTest \
	index3DTest \
	rotationsTest \
	solveQuadraticTest \
	transform2DTest \
	transform3DTest \
	utilitiesTest


# Specify which tests should be build from source.

check_PROGRAMS          = \
	array1DTest \
	array2DTest \
	array3DTest \
	derivativeRiddersTest \
	index3DTest \
	rotationsTest \
	solveQuadraticTest \
	transform2DTest \
	transform3DTest \
	utilitiesTest


# For each test which should be built from source, specify which 
# source files.

array1DTest_SOURCES		= array1DTest.cc arrayTestCommon.hh
array2DTest_SOURCES		= array2DTest.cc arrayTestCommon.hh
array3DTest_SOURCES		= array3DTest.cc arrayTestCommon.hh
derivativeRiddersTest_SOURCES	= derivativeRiddersTest.cc
index3DTest_SOURCES		= index3DTest.cc
rotationsTest_SOURCES		= rotationsTest.cc
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
