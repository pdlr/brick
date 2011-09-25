#
# Makefile.am
#   Process this file with automake to produce Makefile.in.
#

# Set VPATH so we can build in a different directory.

srcdir			= @srcdir@
VPATH			= @srcdir@


# Specify which tests should be run.

TESTS			= \
	linearAlgebraTest


# Specify which tests should be build from source.

check_PROGRAMS          = \
	linearAlgebraTest


# For each test which should be built from source, specify which 
# source files.

linearAlgebraTest_SOURCES		= linearAlgebraTest.cc


# Set the include path.

AM_CPPFLAGS		= -I$(top_srcdir)


# Explicitly set library locations to avoid linking with previously
# installed versions of the libraries.

AM_LDFLAGS		= -L$(top_builddir)/brick/linearAlgebra


# Specify which libraries are needed to build these tests.

LDADD			= \
	-lbrickLinearAlgebra \
	-lbrickNumeric \
	-lbrickPortability \
	-lbrickCommon \
	-lbrickTest \
	-lbrickTestAutoMain \
	-llapack \
	-lblas


# List the files that need to be cleaned.

MAINTAINERCLEANFILES	= Makefile.in
