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


# We need to recurse into the following directories.

SUBDIRS = test


# Set the include path.

AM_CPPFLAGS		= -I$(top_srcdir)


# Define the list of header files to install and the location to
# install them.

brickLinearAlgebradir		= $(includedir)/brick/numeric
brickNumeric_HEADERS	= \
	clapack.h \
	linearAlgebra.hh

# Define the library to build.

lib_LTLIBRARIES		= libbricklinearAlgebra.la
libbrickLinearAlgebra_la_SOURCES = \
	linearAlgebra.cc linearAlgebra.hh


# List the files that need to be cleaned.

MAINTAINERCLEANFILES    = Makefile.in
