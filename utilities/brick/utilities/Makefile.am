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


# We need to recurse into the following directories.

SUBDIRS			= test

# Set the include path.

AM_CPPFLAGS		= -I$(top_srcdir)


# Define the list of header files to install and the location to
# install them.

brickUtilitiesdir		= $(includedir)/brick/utilities
brickUtilities_HEADERS		= \
	date.hh \
	exception.hh \
	imageIO.hh \
	frequencyGoverner.hh \
	lockFile.hh \
	optionDescription.hh \
	optionParser.hh \
	path.hh \
	pythonIO.hh \
	stringManipulation.hh \
	timeUtilities.hh


# Define the library to build.

lib_LTLIBRARIES		= libbrickUtilities.la
libbrickUtilities_la_SOURCES = \
	date.cc date.hh \
	lockFile.cc lockFile.hh \
	optionDescription.cc optionDescription.hh \
	optionParser.cc optionParser.hh \
	path.cc path.hh \
	stringManipulation.cc stringManipulation.hh

# List the files that need to be cleaned.

MAINTAINERCLEANFILES    = Makefile.in
