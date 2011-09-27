#
# Makefile.am
#   Process this file with automake to produce Makefile.in.
#

# Set VPATH so we can build in a different directory.

srcdir			= @srcdir@
VPATH			= @srcdir@


# Set the include path.

AM_CPPFLAGS		= -I$(top_srcdir)


# Define the list of header files to install and the location to
# install them.

brickRandomdir		= $(includedir)/brick/random
brickRandom_HEADERS	= clapack.hh pseudoRandom.hh


# Define the library to build.

lib_LTLIBRARIES		= \
	libbrickRandom.la
libbrickRandom_la_SOURCES = \
	pseudoRandom.cc pseudoRandom.hh


# Define the files to build to perform a check.

TESTS			= 
check_PROGRAMS		= 


# List the files that need to be cleaned.

MAINTAINERCLEANFILES	= Makefile.in
