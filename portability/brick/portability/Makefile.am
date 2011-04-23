#
# Makefile.am
#   Process this file with automake to produce Makefile.in.
#
# $Revision: 975 $
# $Date: 2007-12-30 01:57:17 -0500 (Sun, 30 Dec 2007) $
#


# Set the include path.

AM_CPPFLAGS = -I$(top_builddir) -I$(top_srcdir) -I..


# Define the list of header files to install and the location to
# install them.

brickPortabilitydir		= $(includedir)/brick/portability
brickPortability_HEADERS 	= \
	filesystem.hh \
	standardC.hh \
	timeUtilities.hh


# Define the library to build.

lib_LTLIBRARIES			= libbrickPortability.la
libbrickPortability_la_SOURCES	= \
	filesystem.cc filesystem.hh \
	standardC.cc standardC.hh \
	timeUtilities.cc timeUtilities.hh


# List the files that need to be cleaned.

MAINTAINERCLEANFILES    = Makefile.in

