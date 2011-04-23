#
# Makefile.am
#   Process this file with automake to produce Makefile.in.
#


# Set VPATH so we can build in a different directory.

srcdir			= @srcdir@
VPATH			= @srcdir@


# Specify which tests should be run.

TESTS = \
	byteOrderTest \
	expectTest \
	traceableTest


# Specify which tests should be build from source.

check_PROGRAMS = \
	byteOrderTest \
	expectTest \
	traceableTest


# For each test which should be built from source, specify which 
# source files.

byteOrderTest_SOURCES = \
	byteOrderTest.cc

expectTest_SOURCES = \
	expectTest.cc

traceableTest_SOURCES = \
	traceableTest.cc


# Set the include path.

AM_CPPFLAGS		= -I$(top_srcdir)


# Explicitly set library locations to avoid linking with previously
# installed versions of the libraries.

AM_LDFLAGS		= -L$(top_builddir)/brick/common


# Specify which libraries are needed to build these tests.

LDADD			= -lbrickCommon


# List the files that need to be cleaned.

MAINTAINERCLEANFILES	= Makefile.in
