#
# Makefile.am
#   Process this file with automake to produce Makefile.in.
#
# $Revision: 977 $
# $Date: 2007-12-30 02:50:58 -0500 (Sun, 30 Dec 2007) $
#


# Set VPATH so we can build in a different directory.

srcdir			= @srcdir@
VPATH			= @srcdir@


# Set the include path.

AM_CPPFLAGS		= -I$(top_srcdir)


# Define the list of header files to install and the location to
# install them.

brickTestdir				= $(includedir)/brick/test
brickTest_HEADERS			= \
	autoregister.hh \
	functors.hh \
	runnableObject.hh \
	testException.hh \
	testFixture.hh \
	testMacros.hh


# Define the libraries to build.

lib_LTLIBRARIES				= \
	libbrickTest.la \
	libbrickTestAutoMain.la

libbrickTest_la_SOURCES			= \
	autoregister.cc autoregister.hh

libbrickTestAutoMain_la_SOURCES		= \
	automaticMain.cc autoregister.hh


# Define the files to build to perform a check.


# List the files that need to be cleaned.

MAINTAINERCLEANFILES    		= Makefile.in

