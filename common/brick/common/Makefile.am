#
# Makefile.am
#   Process this file with automake to produce Makefile.in.
#
# $Revision: 1324 $
# $Date: 2010-03-07 20:26:18 -0500 (Sun, 07 Mar 2010) $
#


# Set VPATH so we can build in a different directory.

srcdir			= @srcdir@
VPATH			= @srcdir@


# We need to recurse into the following directories.

SUBDIRS			= test


# Set the include path.

AM_CPPFLAGS		= -I$(top_srcdir)


# Define the list of header files to install and the location to
# install them.

brickCommondir            = $(includedir)/brick/common
brickCommon_HEADERS       = \
	argumentDescription.hh \
	byteOrder.hh \
	compileTimestamp.hh \
	constants.hh \
	exception.hh \
	expect.hh \
	mathFunctions.hh \
	referenceCount.hh \
	stridedPointer.hh \
	traceable.hh \
	triple.hh \
	types.hh


# Define the library to build.

lib_LTLIBRARIES         = libbrickCommon.la
libbrickCommon_la_SOURCES = \
	compileTimestamp.cc compileTimestamp.hh \
	exception.cc exception.hh \
	expect.cc expect.hh \
	traceable.cc traceable.hh


# List the files that need to be cleaned.

MAINTAINERCLEANFILES    = Makefile.in
