#
# Makefile.am
#   Process this file with automake to produce Makefile.in.
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

