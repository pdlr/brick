#
# Makefile.am
#   Process this file with automake to produce Makefile.in.
#


# Set VPATH so we can build in a different directory.

srcdir			= @srcdir@
VPATH			= @srcdir@


# Specify which tests should be run.

TESTS = \
	plane3DTest \
	ray2DTest \
	triangle3DTest \
	utilities2DTest \
	utilities3DTest


# Specify which tests should be build from source.

check_PROGRAMS = \
	plane3DTest \
	ray2DTest \
	triangle3DTest \
	utilities2DTest \
	utilities3DTest


# For each test which should be built from source, specify which 
# source files.

plane3DTest_SOURCES		= plane3DTest.cc
ray2DTest_SOURCES		= ray2DTest.cc
triangle3DTest_SOURCES		= triangle3DTest.cc
utilities2DTest_SOURCES		= utilities2DTest.cc
utilities3DTest_SOURCES		= utilities3DTest.cc


# Set the include path.

AM_CPPFLAGS		= -I$(top_srcdir) -DDLR_TEST_USE_AUTOMATIC_MAIN


# Explicitly set library locations to avoid linking with previously
# installed versions of the libraries.

AM_LDFLAGS		= -L$(top_builddir)/dlrGeometry


# Specify which libraries are needed to build these tests.

LDADD = \
	-lbrickGeometry \
	-lbrickLinearAlgebra \
	-lbrickNumeric \
	-lbrickPortability \
	-lbrickTest \
	-lbrickTestAutoMain \
	-lbrickCommon \
	-llapack -lblas


# List the files that need to be cleaned.

MAINTAINERCLEANFILES	= Makefile.in bar.pgm bar.ppm foo.pgm foo.ppm
